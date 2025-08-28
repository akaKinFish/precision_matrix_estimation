function [alpha_auto, lambda1_auto, selection_stats] = module5_step_size_selection(Gamma_init, K_smooth, W_matrix, params)
% MODULE5_STEP_SIZE_SELECTION - Automatic step size and lambda1 selection via Gershgorin bounds
%
% Syntax:
%   [alpha_auto, lambda1_auto] = module5_step_size_selection(Gamma_init, K_smooth, W_matrix, params)
%   [alpha_auto, lambda1_auto, selection_stats] = module5_step_size_selection(Gamma_init, K_smooth, W_matrix, params)
%
% Description:
%   Computes data-driven step size α and smoothing parameter λ₁ using 
%   Gershgorin circle theorem to bound the Lipschitz constant of the 
%   smooth objective function.
%   
%   Theory:
%   L ≤ L_logdet + L_smooth
%   L_logdet = max_ω ||Γ_ω^(-1)||₂²
%   L_smooth ≤ 2λ₁ K_max ρ(W^Γ)
%   
%   Data-driven setting:
%   λ₁ = δ/(2 K_max R_max)
%   α = 1/(L_logdet + δ)
%
% Input Arguments:
%   Gamma_init - (cell array, Fx1) Initial precision matrices
%   K_smooth   - (double array, FxF) Smoothing kernel k_{ω,ω'}
%   W_matrix   - (double array, pxp) Weight matrix W^Γ
%   params     - (struct) Parameters:
%     .delta              - (double) Safety margin (default: 0.9)
%     .power_iter_tol     - (double) Power iteration tolerance (default: 1e-6)
%     .power_iter_max     - (integer) Max power iterations (default: 50)
%     .use_exact_spectral - (logical) Use exact eigenvalue computation (default: false)
%     .verbose           - (logical) Display computation details (default: false)
%
% Output Arguments:
%   alpha_auto    - (double) Computed step size α
%   lambda1_auto  - (double) Computed smoothing parameter λ₁  
%   selection_stats - (struct) Contains:
%     .L_logdet           - (double) Log-determinant Lipschitz constant
%     .K_max             - (double) Maximum kernel row sum
%     .R_max             - (double) Gershgorin radius bound for W^Γ
%     .lambda1_theory    - (double) Theoretical λ₁ bound
%     .alpha_theory      - (double) Theoretical α bound
%     .condition_numbers - (double array, Fx1) Condition numbers of Γ_ω
%     .spectral_norms    - (double array, Fx1) ||Γ_ω^(-1)||₂ estimates
%
% Examples:
%   % Basic automatic selection
%   [alpha, lambda1] = module5_step_size_selection(Gamma_cells, kernel, W, struct());
%   
%   % With detailed statistics
%   params.verbose = true;
%   params.delta = 0.85;
%   [alpha, lambda1, stats] = module5_step_size_selection(Gamma_cells, kernel, W, params);
%   fprintf('Computed: α=%.4e, λ₁=%.4e\n', alpha, lambda1);
%   fprintf('Theory bounds: L_logdet=%.2e, K_max=%.2e, R_max=%.2e\n', ...
%           stats.L_logdet, stats.K_max, stats.R_max);
%
% Mathematical Background:
%   Uses Gershgorin's circle theorem: for Hermitian matrix A,
%   ρ(A) ≤ max_i Σ_j |A_ij| = R_max
%   
%   This avoids expensive eigenvalue computation while providing rigorous
%   upper bounds on the Lipschitz constant for guaranteed convergence.
%
% See also: MODULE5_PROXIMAL_MAIN, MODULE6_HYPERPARAMETER_CONFIG
%
% Author: [Your Name]
% Date: [Current Date]
% Version: 1.0

% ==================== Input Validation ====================
if nargin < 3
    error('module5_step_size_selection:insufficient_input', ...
          'Gamma_init, K_smooth, and W_matrix are required');
end

if nargin < 4
    params = struct();
end

% Set default parameters
default_params = struct();
default_params.delta = 0.9;
default_params.power_iter_tol = 1e-6;
default_params.power_iter_max = 50;
default_params.use_exact_spectral = false;
default_params.verbose = false;

% Merge parameters
param_names = fieldnames(default_params);
for i = 1:length(param_names)
    if ~isfield(params, param_names{i})
        params.(param_names{i}) = default_params.(param_names{i});
    end
end

% Validate inputs
F = length(Gamma_init);
p = size(Gamma_init{1}, 1);

if ~isequal(size(K_smooth), [F, F])
    error('module5_step_size_selection:kernel_size_mismatch', ...
          'K_smooth must be %dx%d', F, F);
end

if ~isequal(size(W_matrix), [p, p])
    error('module5_step_size_selection:weight_size_mismatch', ...
          'W_matrix must be %dx%d', p, p);
end

% ==================== Initialize Statistics ====================
selection_stats = struct();
selection_stats.condition_numbers = zeros(F, 1);
selection_stats.spectral_norms = zeros(F, 1);

if params.verbose
    fprintf('\n--- Automatic Parameter Selection ---\n');
    fprintf('Computing Lipschitz bounds via Gershgorin theorem...\n');
end

% ==================== Compute L_logdet = max_ω ||Γ_ω^(-1)||₂² ====================
compute_tic = tic;
L_logdet = 0;

for f = 1:F
    Gamma_f = Gamma_init{f};
    
    % Validate positive definiteness
    [~, chol_flag] = chol(Gamma_f);
    if chol_flag ~= 0
        error('module5_step_size_selection:not_psd', ...
              'Gamma_init{%d} is not positive definite', f);
    end
    
    % Compute condition number
    try
        cond_f = cond(Gamma_f);
        selection_stats.condition_numbers(f) = cond_f;
        
        if cond_f > 1e12
            warning('module5_step_size_selection:ill_conditioned', ...
                    'Gamma_init{%d} is ill-conditioned (cond=%.2e)', f, cond_f);
        end
    catch
        selection_stats.condition_numbers(f) = Inf;
    end
    
    % Compute spectral norm of inverse: ||Γ_f^(-1)||₂
    if params.use_exact_spectral
        % Exact computation via eigenvalues (expensive)
        eigenvals = eig(Gamma_f);
        min_eigenval = min(real(eigenvals));
        spectral_norm_inv = 1 / min_eigenval;  % ||Γ^(-1)||₂ = 1/λ_min(Γ)
    else
        % Power iteration estimate (faster)
        spectral_norm_inv = estimate_spectral_norm_inverse(Gamma_f, params);
    end
    
    selection_stats.spectral_norms(f) = spectral_norm_inv;
    L_logdet = max(L_logdet, spectral_norm_inv^2);
end

selection_stats.L_logdet = L_logdet;
compute_time = toc(compute_tic);

if params.verbose
    fprintf('L_logdet computation: %.6e (%.3fs)\n', L_logdet, compute_time);
end

% ==================== Compute K_max = max_ω Σ_{ω'} k_{ω,ω'} ====================
K_row_sums = sum(abs(K_smooth), 2);  % Row sums of |k_{ω,ω'}|
K_max = max(K_row_sums);
selection_stats.K_max = K_max;

if params.verbose
    fprintf('K_max (max kernel row sum): %.6e\n', K_max);
end

% ==================== Compute R_max via Gershgorin ====================
% For Hermitian matrix W^Γ, Gershgorin radius bound
R_max = max( sum(abs(W_matrix), 2) );
selection_stats.R_max = R_max;

if params.verbose
    fprintf('R_max (Gershgorin radius bound): %.6e\n', R_max);
end

% ==================== Data-Driven Parameter Selection ====================
% Compute λ₁ = δ/(2 K_max R_max)
if K_max > 0 && R_max > 0
    lambda1_auto = params.delta / (2 * K_max * R_max);
else
    lambda1_auto = params.delta * 0.01;  % Fallback for degenerate cases
    warning('module5_step_size_selection:degenerate_bounds', ...
            'K_max or R_max is zero, using fallback λ₁');
end

% Compute α = 1/(L_logdet + δ)
alpha_auto = 1 / (L_logdet + params.delta);

% Store theoretical bounds
selection_stats.lambda1_theory = lambda1_auto;
selection_stats.alpha_theory = alpha_auto;

% ==================== Validation and Safety Checks ====================
% Check parameter reasonableness
if lambda1_auto > 1.0
    warning('module5_step_size_selection:large_lambda1', ...
            'Computed λ₁=%.2e is large, may over-smooth', lambda1_auto);
end

if alpha_auto > 1.0
    warning('module5_step_size_selection:large_step_size', ...
            'Computed α=%.2e is large, may be unstable', alpha_auto);
    alpha_auto = min(alpha_auto, 1.0);  % Cap at 1.0
end

if alpha_auto < 1e-8
    warning('module5_step_size_selection:small_step_size', ...
            'Computed α=%.2e is very small, convergence may be slow', alpha_auto);
end

% ==================== Final Summary ====================
if params.verbose
    fprintf('\n--- Parameter Selection Results ---\n');
    fprintf('Safety margin δ: %.3f\n', params.delta);
    fprintf('Theoretical bounds:\n');
    fprintf('  L_logdet: %.6e\n', L_logdet);
    fprintf('  K_max:    %.6e\n', K_max);
    fprintf('  R_max:    %.6e\n', R_max);
    fprintf('Selected parameters:\n');
    fprintf('  λ₁ = %.6e\n', lambda1_auto);
    fprintf('  α  = %.6e\n', alpha_auto);
    fprintf('Average condition number: %.2e\n', mean(selection_stats.condition_numbers));
    fprintf('Max spectral norm ||Γ^(-1)||₂: %.2e\n', max(selection_stats.spectral_norms));
end

end

% ==================== Helper Function ====================
function spectral_norm_inv = estimate_spectral_norm_inverse(Gamma, params)
% Estimate ||Γ^(-1)||₂ via power iteration without explicit inverse computation
%
% Uses: ||A^(-1)||₂ = 1/λ_min(A) where λ_min is smallest eigenvalue

% Power iteration for smallest eigenvalue (inverse iteration)
p = size(Gamma, 1);
v = randn(p, 1) + 1i * randn(p, 1);  % Random complex starting vector
v = v / norm(v);

prev_lambda = Inf;
for iter = 1:params.power_iter_max
    try
        % Solve Γ * w = v (equivalent to w = Γ^(-1) * v)
        w = Gamma \ v;
        w = w / norm(w);
        
        % Rayleigh quotient: λ ≈ v^H * Γ * v
        lambda_est = real(w' * Gamma * w);
        
        % Check convergence
        if abs(lambda_est - prev_lambda) < params.power_iter_tol * abs(lambda_est)
            break;
        end
        
        prev_lambda = lambda_est;
        v = w;
        
    catch ME
        % Linear solve failed - matrix might be singular
        warning('module5_step_size_selection:power_iteration_failed', ...
                'Power iteration failed at step %d: %s', iter, ME.message);
        lambda_est = 1e-12;  % Conservative fallback
        break;
    end
end

% Spectral norm of inverse: ||Γ^(-1)||₂ = 1/λ_min(Γ)
if lambda_est > 1e-15
    spectral_norm_inv = 1 / lambda_est;
else
    spectral_norm_inv = 1e15;  % Very large for nearly singular matrix
    warning('module5_step_size_selection:near_singular', ...
            'Matrix appears nearly singular (λ_min ≈ %.2e)', lambda_est);
end

end