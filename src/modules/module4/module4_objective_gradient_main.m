function gradient_results = module4_objective_gradient_main(input_data, gradient_params)
% MODULE4_OBJECTIVE_GRADIENT_MAIN - Compute smooth part gradient for sparse precision estimation
%
% Syntax:
%   gradient_results = module4_objective_gradient_main(input_data, gradient_params)
%
% Description:
%   Computes the gradient of the smooth part of the objective function for 
%   the sparse precision matrix estimation problem. Handles log-determinant,
%   trace, and cross-frequency smoothing terms with proper Hermitian symmetry.
%   
%   The smooth objective is:
%   f = Σ_ω[-log det(Γ_ω) + tr(Σ_ω Γ_ω)] + λ₁Σ_{ω,ω'} k_{ω,ω'} ||Γ_ω - Γ_{ω'}||²_{W^Γ}
%   
%   Gradient:
%   G_ω = -Γ_ω⁻¹ + Σ_ω + 2λ₁Σ_{ω'} k_{ω,ω'} W^Γ(Γ_ω - Γ_{ω'})
%
% Input Arguments:
%   input_data - (struct) Required fields:
%     .precision_matrices         - (cell array, Fx1) Current precision estimates Γ_ω
%     .whitened_covariances      - (cell array, Fx1) Whitened covariances Σ_ω
%     .smoothing_kernel          - (double array, FxF) Kernel matrix k_{ω,ω'}
%     .weight_matrix             - (double array, pxp) Weight matrix W^Γ
%
%   gradient_params - (struct) Optional parameters:
%     .lambda1                   - (double) Smoothing parameter (default: 0.01)
%     .penalize_diagonal         - (logical) Include diagonal in L1 penalty (default: false)
%     .use_graph_laplacian       - (logical) Use Laplacian method for smoothing (default: true)
%     .chol_tolerance            - (double) Tolerance for Cholesky decomposition (default: 1e-12)
%     .symmetrization_tolerance  - (double) Tolerance for Hermitian check (default: 1e-10)
%     .force_hermitian           - (logical) Force Hermitian symmetry (default: true)
%     .verbose                   - (logical) Display computation progress (default: false)
%
% Output Arguments:
%   gradient_results - (struct) Contains:
%     .smooth_gradients          - (cell array, Fx1) Gradient matrices G_ω
%     .gradient_components       - (struct) Individual gradient components
%       .logdet_gradients        - (cell array, Fx1) Log-determinant gradients
%       .trace_gradients         - (cell array, Fx1) Trace gradients  
%       .smoothing_gradients     - (cell array, Fx1) Smoothing gradients
%     .computation_stats         - (struct) Computation statistics
%     .hermitian_violations      - (double array, Fx1) Hermitian symmetry errors
%     .success                   - (logical) Overall success indicator
%
% Examples:
%   % Basic usage
%   input_data.precision_matrices = gamma_estimates;
%   input_data.whitened_covariances = sigma_whitened;  
%   input_data.smoothing_kernel = kernel_matrix;
%   input_data.weight_matrix = weight_matrix;
%   results = module4_objective_gradient_main(input_data, struct());
%   
%   % Advanced usage with custom parameters
%   params.lambda1 = 0.05;
%   params.use_graph_laplacian = false;
%   params.verbose = true;
%   results = module4_objective_gradient_main(input_data, params);
%
% Mathematical Background:
%   This function computes only the smooth (differentiable) part of the objective.
%   The non-smooth L1 penalty is handled separately via proximal operators.
%   
%   For Hermitian matrices in complex domain, gradients are computed with respect
%   to the Frobenius inner product ⟨A,B⟩ = Re(tr(A^H B)).
%
% See also: MODULE4_OBJECTIVE_EVALUATION, MODULE4_LOG_DET_GRADIENT,
%           MODULE4_TRACE_GRADIENT, MODULE4_SMOOTHING_GRADIENT
%
% Author: [Your Name]
% Date: [Current Date]
% Version: 1.0

% ==================== Input Validation ====================
if nargin < 1
    error('module4_objective_gradient_main:insufficient_input', ...
          'At least input_data is required');
end

if nargin < 2
    gradient_params = struct();
end

% Validate input_data structure
required_fields = {'precision_matrices', 'whitened_covariances', 'smoothing_kernel', 'weight_matrix'};
for i = 1:length(required_fields)
    if ~isfield(input_data, required_fields{i})
        error('module4_objective_gradient_main:missing_field', ...
              'Required field "%s" not found in input_data', required_fields{i});
    end
end

% Extract and validate precision matrices
Gammas = input_data.precision_matrices;
if ~iscell(Gammas) || isempty(Gammas)
    error('module4_objective_gradient_main:invalid_precision', ...
          'precision_matrices must be a non-empty cell array');
end

F = length(Gammas);  % Number of frequencies
p = size(Gammas{1}, 1);  % Number of nodes

% Validate matrix dimensions and properties
for f = 1:F
    if ~isnumeric(Gammas{f}) || ~ismatrix(Gammas{f})
        error('module4_objective_gradient_main:invalid_matrix_type', ...
              'precision_matrices{%d} must be a numeric matrix', f);
    end
    
    if ~isequal(size(Gammas{f}), [p, p])
        error('module4_objective_gradient_main:dimension_mismatch', ...
              'precision_matrices{%d} must be %dx%d, got %dx%d', ...
              f, p, p, size(Gammas{f}, 1), size(Gammas{f}, 2));
    end
    
    % Check positive definiteness via Cholesky
    [~, chol_flag] = chol(Gammas{f});
    if chol_flag ~= 0
        error('module4_objective_gradient_main:not_positive_definite', ...
              'precision_matrices{%d} is not positive definite', f);
    end
end

% Validate whitened covariances
Sigmas = input_data.whitened_covariances;
if ~iscell(Sigmas) || length(Sigmas) ~= F
    error('module4_objective_gradient_main:invalid_covariances', ...
          'whitened_covariances must be a cell array of length %d', F);
end

for f = 1:F
    if ~isnumeric(Sigmas{f}) || ~isequal(size(Sigmas{f}), [p, p])
        error('module4_objective_gradient_main:covariance_dimension_mismatch', ...
              'whitened_covariances{%d} must be %dx%d', f, p, p);
    end
end

% Validate smoothing kernel
K = input_data.smoothing_kernel;
if ~isnumeric(K) || ~isequal(size(K), [F, F])
    error('module4_objective_gradient_main:invalid_kernel', ...
          'smoothing_kernel must be a %dx%d numeric matrix', F, F);
end

% Check kernel symmetry
if norm(K - K', 'fro') > 1e-12
    warning('module4_objective_gradient_main:kernel_not_symmetric', ...
            'smoothing_kernel is not symmetric (error: %.2e)', norm(K - K', 'fro'));
    K = (K + K') / 2;  % Symmetrize
end

% Validate weight matrix
Wg = input_data.weight_matrix;
if ~isnumeric(Wg) || ~isequal(size(Wg), [p, p])
    error('module4_objective_gradient_main:invalid_weight_matrix', ...
          'weight_matrix must be a %dx%d numeric matrix', p, p);
end

% Check weight matrix is Hermitian PSD
if norm(Wg - Wg', 'fro') > 1e-12
    warning('module4_objective_gradient_main:weight_not_hermitian', ...
            'weight_matrix is not Hermitian (error: %.2e)', norm(Wg - Wg', 'fro'));
    Wg = (Wg + Wg') / 2;  % Force Hermitian
end

min_eig_W = min(real(eig(Wg)));
if min_eig_W < -1e-12
    error('module4_objective_gradient_main:weight_not_psd', ...
          'weight_matrix is not positive semi-definite (min eigenvalue: %.2e)', min_eig_W);
end

% ==================== Parameter Setup ====================
defaults = struct();
defaults.lambda1 = 0.01;
defaults.penalize_diagonal = false;
defaults.use_graph_laplacian = true;
defaults.chol_tolerance = 1e-12;
defaults.symmetrization_tolerance = 1e-10;
defaults.force_hermitian = true;
defaults.verbose = false;

field_names = fieldnames(defaults);
for i = 1:numel(field_names)
    fname = field_names{i};
    if ~isfield(gradient_params, fname)
        gradient_params.(fname) = defaults.(fname);
    end
end

% Validate parameters
if gradient_params.lambda1 < 0
    error('module4_objective_gradient_main:invalid_lambda1', ...
          'lambda1 must be non-negative, got %.6f', gradient_params.lambda1);
end

if gradient_params.chol_tolerance <= 0
    error('module4_objective_gradient_main:invalid_chol_tolerance', ...
          'chol_tolerance must be positive, got %.2e', gradient_params.chol_tolerance);
end

% ==================== Initialize Results Structure ====================
gradient_results = struct();
gradient_results.smooth_gradients = cell(F, 1);
gradient_results.gradient_components = struct();
gradient_results.gradient_components.logdet_gradients = cell(F, 1);
gradient_results.gradient_components.trace_gradients = cell(F, 1);
gradient_results.gradient_components.smoothing_gradients = cell(F, 1);
gradient_results.computation_stats = struct();
gradient_results.hermitian_violations = zeros(F, 1);
gradient_results.success = false;

% Initialize computation statistics
stats = struct();
stats.computation_times = zeros(1, 4);  % [logdet, trace, smoothing, total]
stats.method_used = gradient_params.use_graph_laplacian;
stats.lambda1_used = gradient_params.lambda1;
stats.hermitian_enforcement_count = 0;
stats.cholesky_failures = 0;

total_tic = tic;

if gradient_params.verbose
    fprintf('=== Module 4: Objective Gradient Computation ===\n');
    fprintf('Processing %d frequencies, %d nodes\n', F, p);
    if gradient_params.use_graph_laplacian
        method_str = 'Laplacian';
    else
        method_str = 'Direct';
    end
    fprintf('Lambda1: %.6f | Method: %s\n', gradient_params.lambda1, method_str);
    fprintf('--------------------------------------------\n');
end

try
    % ==================== Compute Log-Determinant Gradients ====================
    if gradient_params.verbose, fprintf('Computing log-determinant gradients... '); end
    logdet_tic = tic;
    
    [logdet_grads, logdet_stats] = module4_log_det_gradient(Gammas, gradient_params);
    gradient_results.gradient_components.logdet_gradients = logdet_grads;
    
    stats.computation_times(1) = toc(logdet_tic);
    stats.cholesky_failures = stats.cholesky_failures + logdet_stats.cholesky_failures;
    
    if gradient_params.verbose
        fprintf('completed (%.3fs)\n', stats.computation_times(1));
    end
    
    % ==================== Compute Trace Gradients ====================
    if gradient_params.verbose, fprintf('Computing trace gradients... '); end
    trace_tic = tic;
    
    [trace_grads, trace_stats] = module4_trace_gradient(Sigmas, gradient_params);
    gradient_results.gradient_components.trace_gradients = trace_grads;
    
    stats.computation_times(2) = toc(trace_tic);
    stats.hermitian_enforcement_count = stats.hermitian_enforcement_count + trace_stats.hermitian_enforcement_count;
    
    if gradient_params.verbose
        fprintf('completed (%.3fs)\n', stats.computation_times(2));
    end
    
    % ==================== Compute Smoothing Gradients ====================
    if gradient_params.verbose, fprintf('Computing smoothing gradients... '); end
    smoothing_tic = tic;
    
    [smoothing_grads, smoothing_stats] = module4_smoothing_gradient(Gammas, K, Wg, gradient_params);
    gradient_results.gradient_components.smoothing_gradients = smoothing_grads;
    
    stats.computation_times(3) = toc(smoothing_tic);
    
    if gradient_params.verbose
        fprintf('completed (%.3fs)\n', stats.computation_times(3));
    end
    
    % ==================== Combine Gradients ====================
    if gradient_params.verbose, fprintf('Combining gradient components... '); end
    combine_tic = tic;
    
    for f = 1:F
        % Combine all gradient components
        G_f = logdet_grads{f} + trace_grads{f} + smoothing_grads{f};
        
        % Force Hermitian symmetry if requested
        if gradient_params.force_hermitian
            G_f_symmetric = (G_f + G_f') / 2;
            hermitian_error = norm(G_f - G_f_symmetric, 'fro');
            gradient_results.hermitian_violations(f) = hermitian_error;
            
            if hermitian_error > gradient_params.symmetrization_tolerance
                stats.hermitian_enforcement_count = stats.hermitian_enforcement_count + 1;
                if gradient_params.verbose && hermitian_error > 1e-8
                    fprintf('\n  Warning: Large Hermitian violation at freq %d: %.2e', f, hermitian_error);
                end
            end
            
            G_f = G_f_symmetric;
        end
        
        gradient_results.smooth_gradients{f} = G_f;
    end
    
    stats.computation_times(4) = toc(combine_tic);
    
    if gradient_params.verbose
        fprintf('completed (%.3fs)\n', stats.computation_times(4));
    end
    
    % ==================== Final Statistics ====================
    stats.total_computation_time = toc(total_tic);
    stats.max_hermitian_violation = max(gradient_results.hermitian_violations);
    stats.mean_hermitian_violation = mean(gradient_results.hermitian_violations);
    
    gradient_results.computation_stats = stats;
    gradient_results.success = true;
    
    if gradient_params.verbose
        fprintf('--------------------------------------------\n');
        fprintf('Module 4 Summary:\n');
        fprintf('- Total time: %.3fs\n', stats.total_computation_time);
        fprintf('- Max Hermitian violation: %.2e\n', stats.max_hermitian_violation);
        fprintf('- Hermitian enforcements: %d\n', stats.hermitian_enforcement_count);
        fprintf('- Cholesky failures: %d\n', stats.cholesky_failures);
        fprintf('- Success: YES\n');
        fprintf('============================================\n');
    end
    
catch ME
    stats.total_computation_time = toc(total_tic);
    gradient_results.computation_stats = stats;
    gradient_results.success = false;
    
    if gradient_params.verbose
        fprintf('\n✗ Module 4 failed: %s\n', ME.message);
    end
    
    rethrow(ME);
end

end