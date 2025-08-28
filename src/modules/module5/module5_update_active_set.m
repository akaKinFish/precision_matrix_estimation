function [A_mask_new, update_stats] = module5_update_active_set(A_mask_current, Gamma_current, G_smooth, params)
% MODULE5_UPDATE_ACTIVE_SET - Semi-dynamic active set updates using KKT conditions
%
% Syntax:
%   [A_mask_new, update_stats] = module5_update_active_set(A_mask_current, Gamma_current, G_smooth, params)
%
% Description:
%   Updates active set masks based on KKT optimality conditions with 
%   hysteresis to prevent oscillations. Uses持续性要求 for both adding 
%   and removing edges from the active set.
%   
%   KKT conditions for L1 penalty:
%   - For active edges (i,j) ∈ A: no constraint on |G_ij|
%   - For inactive edges (i,j) ∉ A: must have |G_ij| ≤ λ₂
%   
%   Adding rule: |G_ij| > λ₂ + η for t_hold consecutive calls
%   Removing rule: |Γ_ij| < ε_prune AND |G_ij| < λ₂ - η for t_hold consecutive calls
%
% Input Arguments:
%   A_mask_current - (cell array, Fx1) Current active set masks (logical pxp)
%   Gamma_current  - (cell array, Fx1) Current precision matrices  
%   G_smooth       - (cell array, Fx1) Smooth gradients (no L1 component)
%   params         - (struct) Parameters:
%     .lambda2        - (double) L1 penalty parameter
%     .t_hold         - (integer) Consecutive iterations for changes (default: 2)
%     .eta_margin     - (double) KKT tolerance margin (default: 0.05 * lambda2)
%     .eps_prune      - (double) Magnitude threshold for removal (default: 1e-6)
%     .protect_diagonal - (logical) Never remove diagonal elements (default: true)
%     .max_additions  - (integer) Maximum edges to add per call (default: Inf)
%     .max_removals   - (integer) Maximum edges to remove per call (default: Inf)
%     .verbose        - (logical) Display update statistics (default: false)
%
% Output Arguments:
%   A_mask_new - (cell array, Fx1) Updated active set masks
%   update_stats - (struct) Contains:
%     .edges_added          - (integer) Total edges added across frequencies
%     .edges_removed        - (integer) Total edges removed across frequencies
%     .frequencies_changed  - (integer) Number of frequencies with changes
%     .kkt_violations       - (integer) Number of KKT violations found
%     .add_candidates       - (integer) Edges eligible for adding
%     .remove_candidates    - (integer) Edges eligible for removal
%     .total_active_edges   - (integer) Final total active edges
%     .sparsity_ratio       - (double) Final overall sparsity ratio
%
% Examples:
%   % Basic active set update
%   params.lambda2 = 0.01;
%   [A_new, stats] = module5_update_active_set(A_current, Gamma, G_smooth, params);
%   
%   % Monitor active set evolution
%   if stats.edges_added > 0 || stats.edges_removed > 0
%       fprintf('Active set changed: +%d edges, -%d edges\n', ...
%               stats.edges_added, stats.edges_removed);
%   end
%   
%   % Check for excessive changes (may indicate instability)
%   total_changes = stats.edges_added + stats.edges_removed;
%   if total_changes > 0.1 * stats.total_active_edges
%       warning('Large active set changes detected');
%   end
%
% Mathematical Background:
%   For L1-penalized problems, the KKT conditions state that at optimality:
%   - Active variables: ∂f/∂x_i + λ₂ sign(x_i) = 0
%   - Inactive variables: |∂f/∂x_i| ≤ λ₂
%   
%   The hysteresis mechanism prevents chattering between active/inactive
%   states by requiring sustained violations before making changes.
%
% See also: MODULE5_PROXIMAL_MAIN, MODULE3_ACTIVE_SET_MAIN
%
% Author: [Your Name]
% Date: [Current Date]
% Version: 1.0

% ==================== Input Validation ====================
if nargin < 4
    error('module5_update_active_set:insufficient_input', ...
          'All 4 arguments are required');
end

% Set default parameters
default_params = struct();
default_params.t_hold = 2;
default_params.eta_margin = [];  % Will set to 0.05 * lambda2
default_params.eps_prune = 1e-6;
default_params.protect_diagonal = true;
default_params.max_additions = Inf;
default_params.max_removals = Inf;
default_params.verbose = false;

% Merge parameters
param_names = fieldnames(default_params);
for i = 1:length(param_names)
    if ~isfield(params, param_names{i})
        params.(param_names{i}) = default_params.(param_names{i});
    end
end

% Set eta_margin default based on lambda2
if isempty(params.eta_margin)
    params.eta_margin = 0.05 * params.lambda2;
end

% Validate inputs
F = length(A_mask_current);
if length(Gamma_current) ~= F || length(G_smooth) ~= F
    error('module5_update_active_set:length_mismatch', ...
          'All cell arrays must have length F=%d', F);
end

p = size(A_mask_current{1}, 1);

% ==================== Initialize Persistent State for Hysteresis ====================
% Use persistent variables to track consecutive violation counts
persistent add_violation_count remove_violation_count p_last F_last

if isempty(add_violation_count) || isempty(p_last) || p_last ~= p || F_last ~= F
    add_violation_count = zeros(p, p, F);
    remove_violation_count = zeros(p, p, F);
    p_last = p;
    F_last = F;
end

% （可选入口）允许显式重置
if isfield(params,'reset_counters') && params.reset_counters
    add_violation_count(:) = 0;
    remove_violation_count(:) = 0;
end

% ==================== Initialize Output and Statistics ====================
A_mask_new = A_mask_current;  % Start with current masks

update_stats = struct();
update_stats.edges_added = 0;
update_stats.edges_removed = 0;
update_stats.frequencies_changed = 0;
update_stats.kkt_violations = 0;
update_stats.add_candidates = 0;
update_stats.remove_candidates = 0;

% ==================== Process Each Frequency ====================
for f = 1:F
    A_f = A_mask_current{f};
    G_f = G_smooth{f};
    Gamma_f = Gamma_current{f};
    
    freq_changed = false;
    
    % ==================== Check for Edges to Add ====================
    % Find inactive edges with strong KKT violations: |G_ij| > λ₂ + η
    inactive_mask = ~A_f;
    if params.protect_diagonal
        inactive_mask = inactive_mask & ~eye(p);  % Don't consider diagonal
    end
    
    kkt_add_violations = (abs(G_f) > params.lambda2 + params.eta_margin) & inactive_mask;
    
    % Update violation counters for adding
    add_violation_count(:, :, f) = add_violation_count(:, :, f) + kkt_add_violations;
    add_violation_count(:, :, f) = add_violation_count(:, :, f) .* inactive_mask;  % Reset for active edges
    
    % Find edges with sustained violations
    sustained_add = (add_violation_count(:, :, f) >= params.t_hold);
    
    % Apply adding with limits
    [add_i, add_j] = find(triu(sustained_add, 1));  % Upper triangle only
    n_add_candidates = length(add_i);
    n_to_add = min(n_add_candidates, params.max_additions);
    
    if n_to_add > 0
        % Add edges (with Hermitian pairing)
        for k = 1:n_to_add
            i = add_i(k);
            j = add_j(k);
            A_f(i, j) = true;
            A_f(j, i) = true;  % Hermitian pairing
        end
        
        update_stats.edges_added = update_stats.edges_added + n_to_add;
        freq_changed = true;
        
        % Reset counters for added edges
        for k = 1:n_to_add
            i = add_i(k);
            j = add_j(k);
            add_violation_count(i, j, f) = 0;
            add_violation_count(j, i, f) = 0;
        end
    end
    
    update_stats.add_candidates = update_stats.add_candidates + n_add_candidates;
    
    % ==================== Check for Edges to Remove ====================
    % Find active edges that are weak AND satisfy KKT: |Γ_ij| < ε_prune AND |G_ij| < λ₂ - η
    active_mask = A_f;
    if params.protect_diagonal
        active_mask = active_mask & ~eye(p);  % Don't consider diagonal for removal
    end
    
    magnitude_small = (abs(Gamma_f) < params.eps_prune);
    kkt_remove_ok = (abs(G_f) < params.lambda2 - params.eta_margin);
    kkt_remove_violations = magnitude_small & kkt_remove_ok & active_mask;
    
    % Update violation counters for removing
    remove_violation_count(:, :, f) = remove_violation_count(:, :, f) + kkt_remove_violations;
    remove_violation_count(:, :, f) = remove_violation_count(:, :, f) .* active_mask;  % Reset for inactive edges
    
    % Find edges with sustained removal conditions
    sustained_remove = (remove_violation_count(:, :, f) >= params.t_hold);
    
    % Apply removal with limits
    [remove_i, remove_j] = find(triu(sustained_remove, 1));  % Upper triangle only
    n_remove_candidates = length(remove_i);
    n_to_remove = min(n_remove_candidates, params.max_removals);
    
    if n_to_remove > 0
        % Remove edges (with Hermitian pairing)
        for k = 1:n_to_remove
            i = remove_i(k);
            j = remove_j(k);
            A_f(i, j) = false;
            A_f(j, i) = false;  % Hermitian pairing
        end
        
        update_stats.edges_removed = update_stats.edges_removed + n_to_remove;
        freq_changed = true;
        
        % Reset counters for removed edges
        for k = 1:n_to_remove
            i = remove_i(k);
            j = remove_j(k);
            remove_violation_count(i, j, f) = 0;
            remove_violation_count(j, i, f) = 0;
        end
    end
    
    update_stats.remove_candidates = update_stats.remove_candidates + n_remove_candidates;
    
    % ==================== Update Statistics ====================
    if freq_changed
        update_stats.frequencies_changed = update_stats.frequencies_changed + 1;
    end
    
    % Count total KKT violations (for monitoring)
    all_kkt_add = sum(sum(kkt_add_violations));
    all_kkt_remove = sum(sum(kkt_remove_violations));
    update_stats.kkt_violations = update_stats.kkt_violations + all_kkt_add + all_kkt_remove;
    
    % Store updated mask
    A_mask_new{f} = A_f;
end

% ==================== Compute Final Statistics ====================
% Count total active edges
total_active = 0;
total_possible = 0;

for f = 1:F
    A_f = A_mask_new{f};
    active_f = sum(sum(triu(A_f, 1)));  % Count upper triangle only
    total_active = total_active + active_f;
    total_possible = total_possible + p * (p - 1) / 2;  % Upper triangle per frequency
end

update_stats.total_active_edges = total_active;
update_stats.sparsity_ratio = 1 - (total_active / total_possible);

% ==================== Verbose Output ====================
if params.verbose && (update_stats.edges_added > 0 || update_stats.edges_removed > 0)
    fprintf('  Active set update:\n');
    fprintf('    Added: %d edges\n', update_stats.edges_added);
    fprintf('    Removed: %d edges\n', update_stats.edges_removed);
    fprintf('    Frequencies changed: %d/%d\n', update_stats.frequencies_changed, F);
    fprintf('    Total active edges: %d\n', update_stats.total_active_edges);
    fprintf('    Sparsity ratio: %.3f\n', update_stats.sparsity_ratio);
    fprintf('    KKT violations detected: %d\n', update_stats.kkt_violations);
end

end