function [combined_masks, active_set_stats] = module3_combined_active_set(edge_proxies, threshold_value, options)
% MODULE3_COMBINED_ACTIVE_SET - Construct combined active set from edge and node activity
%
% Syntax:
%   combined_masks = module3_combined_active_set(edge_proxies, threshold_value)
%   [combined_masks, active_set_stats] = module3_combined_active_set(edge_proxies, threshold_value, options)
%
% Description:
%   Constructs the final active set by combining edge-level and node-level
%   activity decisions. The algorithm:
%   1) Identifies active edges: A_edge = {(i,j,f) : c_ij(f) >= τ}
%   2) Computes node activity: r_i(f) = max_{j≠i} c_ij(f)
%   3) Identifies active nodes: A_node = {(i,f) : r_i(f) >= τ}
%   4) Forms combined set: A = {(i,j,f) ∈ A_edge : (i,f),(j,f) ∈ A_node}
%
% Input Arguments:
%   edge_proxies - (cell array, Fx1) Edge proxy matrices from proxy computation
%   threshold_value - (double) Threshold τ for activity determination
%
% Name-Value Arguments:
%   force_diagonal_active - (logical) Always include diagonal in active set
%                          Default: true
%   symmetrize_masks     - (logical) Ensure output masks are symmetric
%                          Default: true
%   min_active_per_node  - (integer) Minimum active edges per active node
%                          Default: 1
%   verbose              - (logical) Display detailed progress. Default: false
%
% Output Arguments:
%   combined_masks - (struct) Contains active set masks:
%     .edge_masks    - (logical array, pxpxF) Edge-level activity masks
%     .node_masks    - (logical array, pxF) Node-level activity masks
%     .combined_masks - (logical array, pxpxF) Final combined activity masks
%   
%   active_set_stats - (struct) Detailed statistics:
%     .active_edges_per_frequency    - (array, Fx1) Number of active edges per frequency
%     .active_nodes_per_frequency    - (array, Fx1) Number of active nodes per frequency
%     .sparsity_ratios_per_frequency - (array, Fx1) Sparsity ratios per frequency
%     .total_active_edges            - Total active edges across all frequencies
%     .average_edges_per_frequency   - Average number of active edges
%     .total_active_nodes            - Total active nodes across all frequencies
%     .average_nodes_per_frequency   - Average number of active nodes
%     .overall_sparsity_ratio        - Overall sparsity ratio (proportion of edges kept)
%     .edge_reduction_ratio          - Edge reduction ratio (1 - overall_sparsity_ratio)
%
% Examples:
%   % Basic usage
%   masks = module3_combined_active_set(edge_proxies, threshold_value);
%   
%   % With detailed statistics and custom options
%   options.force_diagonal_active = false;
%   options.verbose = true;
%   [masks, stats] = module3_combined_active_set(edge_proxies, tau, options);
%
% Mathematical Background:
%   Edge activity: A_edge = {(i,j,f) : c_ij(f) >= τ}
%   Node activity: r_i(f) = max_{j≠i} c_ij(f), A_node = {(i,f) : r_i(f) >= τ}
%   Combined activity: A = A_edge ∩ {(i,j,f) : (i,f),(j,f) ∈ A_node}
%
%   This ensures that active edges connect active nodes, preventing
%   spurious connections to weakly connected nodes.
%
% See also: MODULE3_EDGE_PROXY_COMPUTATION, MODULE3_THRESHOLD_DETERMINATION
%
% Author: [Your Name]
% Date: [Current Date]
% Version: 1.0

% ==================== Input Validation ====================
if nargin < 2
    error('module3_combined_active_set:insufficient_input', ...
          'edge_proxies and threshold_value are required');
end

if nargin < 3
    options = struct();
end

% Validate edge_proxies
if ~iscell(edge_proxies)
    error('module3_combined_active_set:invalid_input', ...
          'edge_proxies must be a cell array');
end

F = numel(edge_proxies);
if F == 0
    error('module3_combined_active_set:empty_input', ...
          'edge_proxies cannot be empty');
end

% Get dimensions and validate matrices
p = size(edge_proxies{1}, 1);
for f = 1:F
    if ~isnumeric(edge_proxies{f}) || ~isequal(size(edge_proxies{f}), [p p])
        error('module3_combined_active_set:invalid_matrix', ...
              'All proxy matrices must be %dx%d numeric, matrix %d is %dx%d', ...
              p, p, f, size(edge_proxies{f}, 1), size(edge_proxies{f}, 2));
    end
    
    if any(~isfinite(edge_proxies{f}(:)))
        error('module3_combined_active_set:non_finite_values', ...
              'Proxy matrix %d contains non-finite values', f);
    end
end

% Validate threshold_value
if ~isnumeric(threshold_value) || ~isscalar(threshold_value) || ~isfinite(threshold_value)
    error('module3_combined_active_set:invalid_threshold', ...
          'threshold_value must be a finite numeric scalar');
end

% ==================== Parameter Setup ====================
defaults = struct();
defaults.force_diagonal_active = true;
defaults.symmetrize_masks = true;
defaults.min_active_per_node = 1;
defaults.verbose = false;

field_names = fieldnames(defaults);
for i = 1:numel(field_names)
    fname = field_names{i};
    if ~isfield(options, fname)
        options.(fname) = defaults.(fname);
    end
end

% ==================== Initialize Output Structures ====================
combined_masks = struct();
combined_masks.edge_masks = false(p, p, F);
combined_masks.node_masks = false(p, F);
combined_masks.combined_masks = false(p, p, F);

if nargout > 1 || options.verbose
    active_set_stats = struct();
    active_set_stats.active_edges_per_frequency = zeros(F, 1);
    active_set_stats.active_nodes_per_frequency = zeros(F, 1);
    active_set_stats.sparsity_ratios_per_frequency = zeros(F, 1);
end

if options.verbose
    fprintf('Constructing combined active set for %d frequencies (%dx%d matrices)\n', F, p, p);
    fprintf('Threshold: τ = %.6f\n', threshold_value);
end

% ==================== Build Active Sets ====================
total_edge_active = 0;
total_node_active = 0;

for f = 1:F
    proxy_matrix = edge_proxies{f};
    
    % Step 1: Edge active set
    edge_mask = (proxy_matrix >= threshold_value);
    if options.symmetrize_masks
        edge_mask = edge_mask | edge_mask';
    end
    if options.force_diagonal_active
        edge_mask(1:p+1:end) = true;
    else
        edge_mask(1:p+1:end) = false;
    end
    combined_masks.edge_masks(:, :, f) = edge_mask;
    
    % Step 2: Node active set
    node_activity = zeros(p, 1);
    for i = 1:p
        other_indices = [1:i-1, i+1:p];
        if ~isempty(other_indices)
            node_activity(i) = max(proxy_matrix(i, other_indices));
        end
    end
    node_mask = (node_activity >= threshold_value);
    
    % Apply minimum active edges constraint
    if options.min_active_per_node > 0
        for i = 1:p
            if node_mask(i)
                edge_count = sum(edge_mask(i, :)) - 1; % Subtract diagonal
                if edge_count < options.min_active_per_node
                    node_mask(i) = false;
                end
            end
        end
    end
    combined_masks.node_masks(:, f) = node_mask;
    
    % Step 3: Combined active set
    combined_mask = false(p, p);
    for i = 1:p
        for j = 1:p
            if edge_mask(i, j) && node_mask(i) && node_mask(j)
                combined_mask(i, j) = true;
            end
        end
    end
    
    if options.symmetrize_masks
        combined_mask = combined_mask | combined_mask';
    end
    combined_masks.combined_masks(:, :, f) = combined_mask;
    
    % Count statistics for this frequency
    off_diagonal_combined = combined_mask & ~eye(p);
    edge_count = sum(off_diagonal_combined(:)) / 2;
    node_count = sum(node_mask);
    
    total_edge_active = total_edge_active + edge_count;
    total_node_active = total_node_active + node_count;
    
    if nargout > 1 || options.verbose
        active_set_stats.active_edges_per_frequency(f) = edge_count;
        active_set_stats.active_nodes_per_frequency(f) = node_count;
        
        total_possible_edges_per_freq = p * (p - 1) / 2;
        active_set_stats.sparsity_ratios_per_frequency(f) = edge_count / total_possible_edges_per_freq;
    end
end

% ==================== Final Statistics Computation ====================
if nargout > 1 || options.verbose
    % CRITICAL FIX: Add all required fields that tests expect
    active_set_stats.total_active_edges = total_edge_active;
    active_set_stats.average_edges_per_frequency = total_edge_active / F;
    active_set_stats.total_active_nodes = total_node_active;
    active_set_stats.average_nodes_per_frequency = total_node_active / F;
    
    % CRITICAL: overall_sparsity_ratio = proportion of edges kept (what tests expect)
    max_possible_edges = F * p * (p - 1) / 2;
    active_set_stats.overall_sparsity_ratio = total_edge_active / max_possible_edges;
    active_set_stats.edge_reduction_ratio = 1 - active_set_stats.overall_sparsity_ratio;
    
    % Additional statistics for completeness
    if F > 0
        active_set_stats.average_sparsity_per_frequency = mean(active_set_stats.sparsity_ratios_per_frequency);
        active_set_stats.std_sparsity_per_frequency = std(active_set_stats.sparsity_ratios_per_frequency);
    else
        active_set_stats.average_sparsity_per_frequency = 0;
        active_set_stats.std_sparsity_per_frequency = 0;
    end
    
    if options.verbose
        fprintf('\n=== Combined Active Set Statistics ===\n');
        fprintf('Total active edges: %d / %d (%.3f%%)\n', ...
                total_edge_active, max_possible_edges, ...
                100 * active_set_stats.overall_sparsity_ratio);
        fprintf('Average edges per frequency: %.1f\n', active_set_stats.average_edges_per_frequency);
        fprintf('Average nodes per frequency: %.1f\n', active_set_stats.average_nodes_per_frequency);
        fprintf('Edge reduction ratio: %.3f\n', active_set_stats.edge_reduction_ratio);
    end
end

% ==================== Final Validation ====================
for f = 1:F
    edge_mask = combined_masks.edge_masks(:, :, f);
    node_mask = combined_masks.node_masks(:, f);
    combined_mask = combined_masks.combined_masks(:, :, f);
    
    % Validate dimensions and types
    if ~isequal(size(edge_mask), [p p]) || ~isequal(size(node_mask), [p 1]) || ...
       ~isequal(size(combined_mask), [p p])
        error('module3_combined_active_set:dimension_error', ...
              'Mask dimensions incorrect at frequency %d', f);
    end
    
    if ~islogical(edge_mask) || ~islogical(node_mask) || ~islogical(combined_mask)
        error('module3_combined_active_set:type_error', ...
              'All masks must be logical at frequency %d', f);
    end
    
    % Check symmetry if requested
    if options.symmetrize_masks
        edge_asymmetry = nnz(edge_mask ~= edge_mask');
        combined_asymmetry = nnz(combined_mask ~= combined_mask');
        
        if edge_asymmetry > 0
            warning('module3_combined_active_set:edge_asymmetry', ...
                    'Edge mask not symmetric at f=%d (%d asymmetric elements)', f, edge_asymmetry);
        end
        
        if combined_asymmetry > 0
            warning('module3_combined_active_set:combined_asymmetry', ...
                    'Combined mask not symmetric at f=%d (%d asymmetric elements)', f, combined_asymmetry);
        end
    end
    
    % Validate that combined mask is subset of edge mask
    if any(combined_mask(:) & ~edge_mask(:))
        error('module3_combined_active_set:logic_error', ...
              'Combined mask contains edges not in edge mask at frequency %d', f);
    end
end

if options.verbose
    fprintf('\n✓ All validation checks passed\n');
end

end