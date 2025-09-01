function output = module6_hyperparameters(input, varargin)
% MODULE6_HYPERPARAMETERS - Automatic hyperparameter configuration using Gershgorin bounds
%
% Syntax:
%   output = module6_hyperparameters(input)
%   output = module6_hyperparameters(input, Name, Value)
%
% Description:
%   Automatically derives smoothing weight λ₁ and step-size α using safe upper bounds.
%   This module computes:
%     1) L_logdet: upper bound for the log-det Hessian block (max spectral norm of Γ^{-1})
%     2) K_max: maximum kernel row sum for frequency smoothing
%     3) R_max: Gershgorin-based upper bound for spectral radius of W^Γ (off-diagonal row sums)
%   Then derives λ₁ and α with safety margins to ensure convergence for proximal-gradient.
%
% Important:
%   - 'use_gershgorin' ONLY affects how L_logdet is estimated:
%       * true  : Gershgorin upper bound on ||Γ^{-1}||₂ (conservative)
%       * false : exact value via 1 / λ_min(Σ) (no explicit matrix inverse)
%   - R_max always uses Gershgorin off-diagonal row sums (diag(W^Γ) is assumed zero or ignored).
%
% Input Arguments:
%   input - (struct) with fields:
%           .whitened_covariances - {F×1 cell} preprocessed covariance matrices Σ̃_ω (Hermitian)
%           .kernel_matrix        - (F×F double) smoothing kernel K
%           .weight_matrix        - (p×p double) edge weight matrix W^Γ (diagonal ignored)
%           .active_set_mask      - {F×1 cell} optional binary masks per frequency (p×p)
%
% Name-Value Arguments:
%   safety_margin        - (double) δ ∈ (0,1). Default: 0.9
%   lambda2_scale        - (double) positive scale for empirical λ₂ suggestion. Default: 0.01
%   initialization_method- (string) initial Γ method when Gershgorin is used. Default: 'inverse'
%   verbose              - (logical) print diagnostics. Default: true
%   use_gershgorin       - (logical) for L_logdet estimation. Default: true
%
% Output Arguments:
%   output - (struct) with fields:
%       .lambda1            - (double) smoothing parameter λ₁
%       .alpha              - (double) step-size α
%       .lambda2_suggested  - (double) suggested λ₂ (empirical)
%       .diagnostics        - (struct) detail diagnostics:
%           .L_logdet
%           .K_max
%           .R_max
%           .L_logdet_per_frequency
%           .condition_numbers
%           .gershgorin_bounds
%           .kernel_row_sums
%           .weight_row_sums
%           .lambda2_policy
%           .active_density
%           .theoretical_convergence_rate   % equals alpha * L_logdet
%       .computation_info   - (struct) timing, methods, sizes
%       .quality_assessment - (struct) simple sanity checks
%
% See also: MODULE4_OBJECTIVE_GRADIENT, MODULE5_PROXIMAL_MAIN, MODULE1_PREPROCESSING
%
% Author: Precision Matrix Project Team
% Date: 2025-08-28
% Version: 1.1

%% Input validation and parameter parsing
start_time = tic;

if ~isstruct(input)
    error('module6_hyperparameters:invalid_input', 'Input must be a structure');
end

required_fields = {'whitened_covariances', 'kernel_matrix', 'weight_matrix'};
for i = 1:numel(required_fields)
    if ~isfield(input, required_fields{i})
        error('module6_hyperparameters:missing_field', ...
              'Required field "%s" not found in input', required_fields{i});
    end
end

% Parse optional parameters
p = inputParser;
addParameter(p, 'safety_margin', 0.9, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
addParameter(p, 'lambda2_scale', 0.01, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'initialization_method', 'inverse', @(x) ischar(x) || isstring(x));
addParameter(p, 'verbose', true, @islogical);
addParameter(p, 'use_gershgorin', true, @islogical);
parse(p, varargin{:});

delta          = p.Results.safety_margin;
lambda2_scale  = p.Results.lambda2_scale;
init_method    = char(p.Results.initialization_method);
verbose        = p.Results.verbose;
use_gershgorin = p.Results.use_gershgorin;

% Extract input data
Sigma_whitened = input.whitened_covariances;
K       = input.kernel_matrix;
W_gamma = input.weight_matrix;

% Basic sizes
F = numel(Sigma_whitened);
if F == 0
    error('module6_hyperparameters:empty_covariances', ...
          'Input.whitened_covariances must be a non-empty cell array');
end
[n_nodes, m_nodes] = size(Sigma_whitened{1});
if n_nodes ~= m_nodes
    error('module6_hyperparameters:dimension_mismatch', ...
          'First covariance must be square, got %dx%d', n_nodes, m_nodes);
end

% Check all frequency matrices: size, finite, Hermitian up to tolerance
tolH = 1e-12;
for f = 1:F
    A = Sigma_whitened{f};
    if ~isnumeric(A) || ndims(A) ~= 2
        error('module6_hyperparameters:invalid_cov_matrix', ...
              'Covariance at freq %d is not a 2D numeric matrix', f);
    end
    [r, c] = size(A);
    if r ~= n_nodes || c ~= n_nodes
        error('module6_hyperparameters:dimension_mismatch', ...
              'Covariance %d size mismatch: expected %dx%d, got %dx%d', f, n_nodes, n_nodes, r, c);
    end
    if ~all(isfinite(A(:)))
        error('module6_hyperparameters:nonfinite', ...
              'Covariance %d contains NaN/Inf', f);
    end
    % Hermitian check with relative tolerance
    relErr = norm(A - A', 'fro') / max(1, norm(A, 'fro'));
    if relErr > tolH
        warning('module6_hyperparameters:not_hermitian', ...
                'Covariance %d deviates from Hermitian (rel.err=%.2e)', f, relErr);
        % symmetrize to be safe
        Sigma_whitened{f} = (A + A')/2;
    end
end

% Validate kernel and weight matrix sizes
if size(K,1) ~= F || size(K,2) ~= F
    error('module6_hyperparameters:dimension_mismatch', ...
          'Kernel matrix must be F×F, got %dx%d for F=%d', size(K,1), size(K,2), F);
end
if size(W_gamma,1) ~= n_nodes || size(W_gamma,2) ~= n_nodes
    error('module6_hyperparameters:dimension_mismatch', ...
          'Weight matrix must be p×p, got %dx%d for p=%d', size(W_gamma,1), size(W_gamma,2), n_nodes);
end

if verbose
    fprintf('Module 6: Hyperparameter Configuration\n');
    fprintf('=====================================\n');
    fprintf('Problem dimensions: %d nodes, %d frequencies\n', n_nodes, F);
    fprintf('Safety margin δ: %.2f\n', delta);
    fprintf('L_logdet method: %s\n', tern(use_gershgorin, ...
        'Gershgorin upper bound on ||Γ^{-1}||₂', ...
        'exact via 1 / λ_min(Σ)'));
end

%% Step 1: Compute L_logdet (maximum spectral norm of Γ^{-1})
if verbose, fprintf('\nStep 1: Computing L_logdet...\n'); end

% Diagnostics buffers
condition_numbers = zeros(F,1);
gershgorin_bounds = zeros(F,1);
L_logdet_values   = zeros(F,1);

% When using Gershgorin for L_logdet we need an initial Γ^{-1} estimate
Gamma_init = [];
if use_gershgorin
    Gamma_init = cell(F,1);
    for f = 1:F
        A = Sigma_whitened{f};
        % Small diagonal regularization for numerical robustness
        reg_param = 1e-10 * trace(real(A)) / n_nodes;
        A_reg = A + reg_param * eye(n_nodes);
        switch lower(init_method)
            case 'inverse'
                % NOTE: only used here to form a matrix for Gershgorin bound
                Gamma_init{f} = inv(A_reg);
            case 'pseudoinverse'
                Gamma_init{f} = pinv(A_reg);
            case 'diagonal'
                d = diag(A_reg);
                d = real(d);
                d = max(d, eps);
                Gamma_init{f} = diag(1 ./ d);
            otherwise
                error('module6_hyperparameters:invalid_init_method', ...
                      'Unknown initialization method: %s', init_method);
        end
    end
end

for f = 1:F
    A = Sigma_whitened{f};
    % Condition number for diagnostics (on Σ)
    condition_numbers(f) = cond(A);

    if use_gershgorin
        % Gershgorin bound on ||Γ^{-1}||₂ using current Γ^{-1} estimate
        G = Gamma_init{f};
        diagG = diag(G);
        offsum = sum(abs(G),2) - abs(diagG);
        max_bound = max(abs(diagG) + offsum);
        gershgorin_bounds(f) = max_bound;
        L_logdet_values(f)   = max_bound^2;
    else
        % Exact: ||Γ^{-1}||₂^2 = 1 / λ_min(Σ)^2  (no explicit inverse)
        reg_param = 1e-10 * trace(real(A)) / n_nodes;
        A_reg = A + reg_param * eye(n_nodes);
        % Use eig to get smallest eigenvalue (Σ is Hermitian PSD after symmetrization)
        lam_min = min(eig((A_reg + A_reg')/2, 'vector'));
        lam_min = max(lam_min, eps); % protect against numerical non-positivity
        L_logdet_values(f) = (1 / lam_min)^2;
    end
end

L_logdet = max(L_logdet_values);

if verbose
    fprintf('  L_logdet = %.4e (max over %d frequencies)\n', L_logdet, F);
    fprintf('  Mean condition number of Σ: %.2e\n', mean(condition_numbers));
    if use_gershgorin
        fprintf('  Gershgorin ||Γ^{-1}||₂ bound range: [%.2e, %.2e]\n', ...
            min(gershgorin_bounds), max(gershgorin_bounds));
    end
end

%% Step 2: Compute K_max (maximum kernel row sum)
if verbose, fprintf('\nStep 2: Computing K_max...\n'); end

if ~issymmetric(K)
    if verbose
        warning('module6_hyperparameters:kernel_not_symmetric', ...
            'Kernel matrix is not symmetric, symmetrizing');
    end
    K = (K + K')/2;
end
kernel_row_sums = sum(abs(K), 2);
K_max = max(kernel_row_sums);

if verbose
    fprintf('  K_max = %.4f (max row sum)\n', K_max);
    fprintf('  Kernel row sum range: [%.4f, %.4f]\n', min(kernel_row_sums), max(kernel_row_sums));
end

%% Step 3: Compute R_max (maximum off-diagonal row sum of W^Γ)
if verbose, fprintf('\nStep 3: Computing R_max...\n'); end

if ~issymmetric(W_gamma)
    if verbose
        warning('module6_hyperparameters:weight_not_symmetric', ...
            'Weight matrix is not symmetric, symmetrizing');
    end
    W_gamma = (W_gamma + W_gamma')/2;
end

% Gershgorin off-diagonal row-sum bound (diag ignored)
off_diag_row_sums = sum(abs(W_gamma), 2) - abs(diag(W_gamma));
R_max = max(off_diag_row_sums);

if verbose
    fprintf('  R_max = %.4f (max off-diagonal row sum)\n', R_max);
    fprintf('  Off-diagonal row sum range: [%.4f, %.4f]\n', ...
        min(off_diag_row_sums), max(off_diag_row_sums));
end

%% Step 4: Set λ₁ and α with safety margins
if verbose, fprintf('\nStep 4: Computing hyperparameters...\n'); end

% λ₁: if K_max*R_max≈0, fall back to a small positive number
if K_max * R_max < eps
    warning('module6_hyperparameters:zero_product', ...
        'K_max * R_max is near zero, using fallback λ₁');
    lambda1 = 1e-6;  % conservative fallback
else
    lambda1 = delta / (2 * K_max * R_max);
end

% α: unified safe formula
alpha = 1 / (max(L_logdet, eps) + delta);

% λ₂ suggestion
lambda2_policy  = '';
active_density  = NaN;
if isfield(input, 'active_set_mask') && ~isempty(input.active_set_mask) ...
        && iscell(input.active_set_mask) && numel(input.active_set_mask) == F
    total_active_edges = 0;
    total_possible_edges = 0;
    valid_masks = true;
    for f = 1:F
        mask = input.active_set_mask{f};
        if ~isnumeric(mask) && ~islogical(mask)
            valid_masks = false; break;
        end
        [r,c] = size(mask);
        if r ~= n_nodes || c ~= n_nodes
            valid_masks = false; break;
        end
        mask = logical(mask);
        total_active_edges = total_active_edges + sum(mask(triu(true(n_nodes),1)));
        total_possible_edges = total_possible_edges + n_nodes*(n_nodes-1)/2;
    end
    if valid_masks && total_possible_edges > 0
        active_density = total_active_edges / total_possible_edges;
        lambda2_suggested = lambda2_scale * active_density;
        lambda2_policy = 'density_scaled';
    else
        lambda2_suggested = lambda2_scale * sqrt(log(n_nodes) / n_nodes);
        lambda2_policy = 'size_scaled';
    end
else
    lambda2_suggested = lambda2_scale * sqrt(log(n_nodes) / n_nodes);
    lambda2_policy = 'size_scaled';
end

if verbose
    fprintf('  λ₁ = %.6e (smoothing)\n', lambda1);
    fprintf('  α  = %.6e (step size)\n', alpha);
    fprintf('  λ₂ = %.6e (suggested)\n', lambda2_suggested);
    fprintf('  Rationale: λ₁ controls smoothness Lipschitz, α enforces descent, λ₂ is empirical.\n');
end

%% Step 5: Diagnostics and safety checks
if verbose, fprintf('\nStep 5: Validation and diagnostics...\n'); end

theoretical_rate = alpha * L_logdet;  % should be < 1 for safety

% Kernel spectrum summary (optional)
kernel_eigenvals = eig((K+K')/2);
kernel_effective_rank = sum(kernel_eigenvals > 1e-10 * max(kernel_eigenvals));

% Simple connectivity statistics for W
weight_connectivity = (nnz(W_gamma) - n_nodes) / max(1, (n_nodes^2 - n_nodes)); % offdiag density approx.

stability_ok = alpha > 0 && lambda1 > 0 && isfinite(alpha) && isfinite(lambda1);
bounds_reasonable = L_logdet < 1e12 && K_max < 1e8 && R_max < 1e8;

if ~stability_ok
    error('module6_hyperparameters:stability_issue', ...
          'Computed hyperparameters are not stable: α=%.2e, λ₁=%.2e', alpha, lambda1);
end
if ~bounds_reasonable
    warning('module6_hyperparameters:large_bounds', ...
            'Very large bounds detected: L=%.2e, K=%.2e, R=%.2e', L_logdet, K_max, R_max);
end

if verbose
    fprintf('  alpha * L_logdet = %.4f\n', theoretical_rate);
    fprintf('  Kernel effective rank: %d/%d\n', kernel_effective_rank, F);
    fprintf('  Weight connectivity (incl. diag): %.2f%%\n', 100*weight_connectivity);
    fprintf('  Stability check: %s\n', tern(stability_ok, 'PASS', 'FAIL'));
end

%% Assemble output
total_time = toc(start_time);

output = struct();
output.lambda1 = lambda1;
output.alpha = alpha;
output.lambda2_suggested = lambda2_suggested;

% diagnostics
output.diagnostics = struct();
output.diagnostics.L_logdet = L_logdet;
output.diagnostics.K_max = K_max;
output.diagnostics.R_max = R_max;
output.diagnostics.L_logdet_per_frequency = L_logdet_values;
output.diagnostics.condition_numbers = condition_numbers;
output.diagnostics.gershgorin_bounds = gershgorin_bounds;
output.diagnostics.kernel_row_sums = kernel_row_sums;
output.diagnostics.weight_row_sums = off_diag_row_sums;
output.diagnostics.theoretical_convergence_rate = theoretical_rate;
output.diagnostics.kernel_effective_rank = kernel_effective_rank;
output.diagnostics.weight_connectivity = weight_connectivity;
output.diagnostics.lambda2_policy = lambda2_policy;
output.diagnostics.active_density = active_density;
output.diagnostics.stability_check = stability_ok;

% metadata
output.computation_info = struct();
output.computation_info.total_time_seconds = total_time;
output.computation_info.initialization_method = init_method;
output.computation_info.eigenvalue_method = tern(use_gershgorin, 'gershgorin', 'exact-min-eig');
output.computation_info.safety_margin_used = delta;
output.computation_info.n_nodes = n_nodes;
output.computation_info.n_frequencies = F;

% quality gates
output.quality_assessment = struct();
output.quality_assessment.lambda1_reasonable = (lambda1 > 1e-8) && (lambda1 < 1);
output.quality_assessment.alpha_reasonable   = (alpha > 1e-6) && (alpha < 1);
output.quality_assessment.condition_acceptable = all(condition_numbers < 1e12);
output.quality_assessment.overall_quality = ...
    output.quality_assessment.lambda1_reasonable && ...
    output.quality_assessment.alpha_reasonable   && ...
    output.quality_assessment.condition_acceptable;

if verbose
    fprintf('\nModule 6 completed in %.3f seconds\n', total_time);
    fprintf('Configuration quality: %s\n', ...
        tern(output.quality_assessment.overall_quality, 'GOOD', 'NEEDS ATTENTION'));
    if ~output.quality_assessment.overall_quality
        if ~output.quality_assessment.lambda1_reasonable
            fprintf('  - λ₁ = %.2e may be too extreme\n', lambda1);
        end
        if ~output.quality_assessment.alpha_reasonable
            fprintf('  - α = %.2e may be too extreme\n', alpha);
        end
        if ~output.quality_assessment.condition_acceptable
            fprintf('  - Some covariance matrices are ill-conditioned\n');
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s = tern(cond, a, b)
% TERN - Simple ternary helper
if cond, s = a; else, s = b; end
end
