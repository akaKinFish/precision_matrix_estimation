function [l1_subgradients, computation_stats] = module4_l1_subgradient(precision_matrices, params)
% MODULE4_L1_SUBGRADIENT - Compute subgradients of L1 penalty terms
%
% Syntax:
%   [l1_subgradients, computation_stats] = module4_l1_subgradient(precision_matrices, params)
%
% Description:
%   Computes subgradients of the L1 penalty ||Γ_ω||_L1 for each precision matrix.
%   Handles complex-valued matrices with proper Hermitian symmetry constraints.
%   
%   This function is primarily for analysis - in practice, L1 penalties are
%   handled via proximal operators (soft thresholding) in the optimization loop.
%   
%   Subgradient formula:
%   For z ≠ 0: ∂|z| = z/|z|
%   For z = 0: ∂|z| ∈ {u : |u| ≤ 1}
%   
%   For Hermitian constraint: ∂|Γ_ij| = conj(∂|Γ_ji|)
%
% Input Arguments:
%   precision_matrices - (cell array, Fx1) Precision matrices {Γ_ω}
%   params - (struct) Parameters:
%     .lambda2                  - (double) L1 penalty weight (default: 0.01)
%     .penalize_diagonal        - (logical) Include diagonal entries (default: false)
%     .zero_threshold           - (double) Threshold for zero detection (default: 1e-12)
%     .force_hermitian          - (logical) Enforce Hermitian subgradients (default: true)
%     .subgradient_selection    - (string) Method for zero subgradients:
%                                'zero'|'random'|'minimal' (default: 'zero')
%     .verbose                  - (logical) Display computation info (default: false)
%
% Output Arguments:
%   l1_subgradients - (cell array, Fx1) Subgradient matrices
%   computation_stats - (struct) Contains:
%     .computation_time         - (double) Total computation time
%     .zero_entries_count       - (integer array, Fx1) Count of zero entries per freq
%     .nonzero_entries_count    - (integer array, Fx1) Count of nonzero entries per freq
%     .hermitian_violations     - (double array, Fx1) Symmetry constraint violations
%     .subgradient_norms        - (double array, Fx1) Frobenius norms of subgradients
%
% Examples:
%   % Basic subgradient computation
%   params = struct('lambda2', 0.02, 'penalize_diagonal', false);
%   [subgrads, stats] = module4_l1_subgradient(gammas, params);
%   
%   % Check subgradient properties
%   for f = 1:length(subgrads)
%       fprintf('Freq %d: %d zeros, %d nonzeros\n', f, ...
%               stats.zero_entries_count(f), stats.nonzero_entries_count(f));
%   end
%   
%   % Verify subgradient optimality conditions
%   % At optimum: smooth_gradient + l1_subgradient = 0
%   combined_grad = smooth_grads{1} + subgrads{1};
%   optimality_residual = norm(combined_grad, 'fro');
%
% Mathematical Background:
%   The L1 penalty is non-smooth, so we compute subgradients instead of gradients.
%   For complex matrices with Hermitian constraints, the subgradient must also
%   satisfy Hermitian symmetry: S_ji = conj(S_ij).
%   
%   This ensures that the resulting update maintains Hermitian structure.
%
% See also: MODULE4_OBJECTIVE_GRADIENT_MAIN, MODULE4_OBJECTIVE_EVALUATION
%
% Author: [Your Name]
% Date: [Current Date]
% Version: 1.0

% ==================== Input Validation ====================
if nargin < 1
    error('module4_l1_subgradient:insufficient_input', ...
          'precision_matrices is required');
end

if nargin < 2
    params = struct();
end

if ~iscell(precision_matrices) || isempty(precision_matrices)
    error('module4_l1_subgradient:invalid_input', ...
          'precision_matrices must be a non-empty cell array');
end

% Set default parameters
defaults = struct();
defaults.lambda2 = 0.01;
defaults.penalize_diagonal = false;
defaults.zero_threshold = 1e-12;
defaults.force_hermitian = true;
defaults.subgradient_selection = 'zero';  % 'zero', 'random', 'minimal'
defaults.verbose = false;

field_names = fieldnames(defaults);
for i = 1:numel(field_names)
    fname = field_names{i};
    if ~isfield(params, fname)
        params.(fname) = defaults.(fname);
    end
end

% Validate parameters
if params.lambda2 < 0
    error('module4_l1_subgradient:invalid_lambda2', ...
          'lambda2 must be non-negative, got %.6f', params.lambda2);
end

valid_selections = {'zero', 'random', 'minimal'};
if ~ismember(params.subgradient_selection, valid_selections)
    error('module4_l1_subgradient:invalid_selection', ...
          'subgradient_selection must be one of: %s', strjoin(valid_selections, ', '));
end

% ==================== Initialize ====================
F = length(precision_matrices);
p = size(precision_matrices{1}, 1);

l1_subgradients = cell(F, 1);
computation_stats = struct();
computation_stats.computation_time = 0;
computation_stats.zero_entries_count = zeros(F, 1);
computation_stats.nonzero_entries_count = zeros(F, 1);
computation_stats.hermitian_violations = zeros(F, 1);
computation_stats.subgradient_norms = zeros(F, 1);

main_tic = tic;

% Early return if lambda2 is zero
if params.lambda2 == 0
    for f = 1:F
        l1_subgradients{f} = zeros(size(precision_matrices{f}));
    end
    computation_stats.computation_time = toc(main_tic);
    if params.verbose
        fprintf('  L1 subgradient: λ₂=0, returning zero subgradients\n');
    end
    return;
end

if params.verbose
    fprintf('Computing L1 subgradients (λ₂=%.6f)...\n', params.lambda2);
end

% ==================== Process Each Frequency ====================
for f = 1:F
    Gamma_f = precision_matrices{f};
    
    % Validate matrix
    if ~isnumeric(Gamma_f) || ~isequal(size(Gamma_f), [p, p])
        error('module4_l1_subgradient:invalid_matrix', ...
              'precision_matrices{%d} must be %dx%d numeric matrix', f, p, p);
    end
    
    % Initialize subgradient matrix
    subgrad_f = zeros(p, p);
    zero_count = 0;
    nonzero_count = 0;
    
    % Compute subgradient for each entry
    for i = 1:p
        for j = 1:p
            % Check if this entry should be penalized
            should_penalize = params.penalize_diagonal || (i ~= j);
            
            if should_penalize
                gamma_ij = Gamma_f(i, j);
                abs_gamma_ij = abs(gamma_ij);
                
                if abs_gamma_ij > params.zero_threshold
                    % Non-zero entry: ∂|z| = z/|z|
                    subgrad_f(i, j) = params.lambda2 * (gamma_ij / abs_gamma_ij);
                    nonzero_count = nonzero_count + 1;
                else
                    % Zero (or near-zero) entry: ∂|z| ∈ {u : |u| ≤ 1}
                    zero_count = zero_count + 1;
                    
                    switch params.subgradient_selection
                        case 'zero'
                            subgrad_f(i, j) = 0;
                        case 'random'
                            % Random subgradient with |u| ≤ λ₂
                            random_phase = 2 * pi * rand();
                            random_magnitude = params.lambda2 * rand();
                            subgrad_f(i, j) = random_magnitude * exp(1i * random_phase);
                        case 'minimal'
                            % Minimal norm subgradient (usually zero)
                            subgrad_f(i, j) = 0;
                    end
                end
            else
                % Not penalized (typically diagonal entries)
                subgrad_f(i, j) = 0;
            end
        end
    end
    
    % ==================== Enforce Hermitian Constraint ====================
    if params.force_hermitian
        % For Hermitian matrices: S_ji = conj(S_ij)
        subgrad_f_hermitian = zeros(p, p);
        
        for i = 1:p
            for j = 1:p
                if i <= j
                    % Upper triangular: keep as computed
                    subgrad_f_hermitian(i, j) = subgrad_f(i, j);
                    
                    if i ~= j
                        % Lower triangular: conjugate of upper
                        subgrad_f_hermitian(j, i) = conj(subgrad_f(i, j));
                    end
                else
                    % Already handled in upper triangular case
                end
            end
        end
        
        % Measure Hermitian violation
        hermitian_error = norm(subgrad_f - subgrad_f_hermitian, 'fro');
        computation_stats.hermitian_violations(f) = hermitian_error;
        
        subgrad_f = subgrad_f_hermitian;
    end
    
    l1_subgradients{f} = subgrad_f;
    computation_stats.zero_entries_count(f) = zero_count;
    computation_stats.nonzero_entries_count(f) = nonzero_count;
    computation_stats.subgradient_norms(f) = norm(subgrad_f, 'fro');
    
    if params.verbose && f <= 3  % Show details for first few frequencies
        fprintf('  Freq %d: %d zeros, %d nonzeros, ||∂||_F = %.6f\n', ...
                f, zero_count, nonzero_count, computation_stats.subgradient_norms(f));
    end
end

% ==================== Final Statistics ====================
computation_stats.computation_time = toc(main_tic);
computation_stats.total_zero_entries = sum(computation_stats.zero_entries_count);
computation_stats.total_nonzero_entries = sum(computation_stats.nonzero_entries_count);
computation_stats.max_hermitian_violation = max(computation_stats.hermitian_violations);

if params.verbose
    fprintf('  L1 subgradient summary:\n');
    fprintf('    Total computation time: %.3fs\n', computation_stats.computation_time);
    fprintf('    Total zero entries: %d\n', computation_stats.total_zero_entries);
    fprintf('    Total nonzero entries: %d\n', computation_stats.total_nonzero_entries);
    fprintf('    Max Hermitian violation: %.2e\n', computation_stats.max_hermitian_violation);
end

end