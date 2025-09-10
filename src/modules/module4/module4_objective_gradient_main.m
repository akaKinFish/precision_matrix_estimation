function gradient_results = module4_objective_gradient_main(input_data, gradient_params)
% MODULE4_OBJECTIVE_GRADIENT_MAIN
% Compute the gradient of the smooth part of the objective for {Gamma_f}:
%   -logdet(Gamma_f) + trace(Sigma_f * Gamma_f)
% plus cross-frequency smoothing and (optional) single-frequency spatial smoothing.
%
% Inputs:
%   input_data (struct) with required fields:
%       .precision_matrices    (cell{F,1})  Gamma_f
%       .whitened_covariances  (cell{F,1})  Sigma_f
%       .smoothing_kernel      (F x F)      K
%       .weight_matrix         (p x p)      W^Gamma
%
%   gradient_params (struct) optional fields:
%       % existing:
%       .lambda1                     double (>=0), cross-frequency smooth weight
%       .weight_mode                 'matrix'|'hadamard' (default: 'matrix')
%       .use_graph_laplacian         logical
%       .chol_tolerance              double (default: 1e-12)
%       .symmetrization_tolerance    double (default: 1e-12)
%       .force_hermitian             logical (default: true)
%       .verbose                     logical (default: false)
%
%       % NEW (spatial single-frequency smoothing):
%       .lambda3                     double (>=0, default: 0 -> disabled)
%       .spatial_graph_matrix        (p x p) adjacency/Laplacian for spatial term
%       .spatial_graph_is_laplacian  logical (default: false)
%       .spatial_weight_mode         'node'|'hadamard' (default: 'node')
%
% Outputs:
%   gradient_results (struct)
%       .smooth_gradients    (cell{F,1}) gradients d/dGamma_f of smooth objective
%       .gradient_components (struct)    diagnostics:
%           .logdet_gradients
%           .trace_gradients
%           .smoothing_gradients     (cross-frequency)
%           .spatial_gradients       (NEW, present only if lambda3>0)
%       .computation_stats  (struct)     timing, validation flags
%       .success            (logical)
%
% Notes:
%   - Non-smooth L1 (lambda2) is handled via proximal operator in Module 5.
%   - Gradients are symmetrized to be Hermitian-consistent.

% ==================== Defaults & guards ====================
if nargin < 1
    error('module4_objective_gradient_main:insufficient_input', ...
          'At least input_data is required');
end
if nargin < 2, gradient_params = struct(); end

% Required input fields
req_fields = {'precision_matrices','whitened_covariances','smoothing_kernel','weight_matrix'};
for i = 1:numel(req_fields)
    if ~isfield(input_data, req_fields{i})
        error('module4_objective_gradient_main:missing_field', ...
              'Required field "%s" not found in input_data', req_fields{i});
    end
end

% Extract
Gammas = input_data.precision_matrices;
if ~iscell(Gammas) || isempty(Gammas)
    error('module4_objective_gradient_main:invalid_precision', ...
          'precision_matrices must be a non-empty cell array');
end
F = numel(Gammas);
p = size(Gammas{1},1);

% Validate Gamma
for f = 1:F
    if ~isnumeric(Gammas{f}) || ~ismatrix(Gammas{f})
        error('module4_objective_gradient_main:invalid_matrix_type', ...
              'precision_matrices{%d} must be a numeric matrix', f);
    end
    if ~isequal(size(Gammas{f}), [p,p])
        error('module4_objective_gradient_main:dimension_mismatch', ...
              'precision_matrices{%d} must be %dx%d', f, p, p);
    end
    [~, cf] = chol((Gammas{f}+Gammas{f}')/2);
    if cf ~= 0
        error('module4_objective_gradient_main:not_positive_definite', ...
              'precision_matrices{%d} is not positive definite', f);
    end
end

% Covariances
Sigmas = input_data.whitened_covariances;
if ~iscell(Sigmas) || numel(Sigmas) ~= F
    error('module4_objective_gradient_main:invalid_covariances', ...
          'whitened_covariances must be a cell array of length %d', F);
end
for f = 1:F
    if ~isnumeric(Sigmas{f}) || ~isequal(size(Sigmas{f}), [p,p])
        error('module4_objective_gradient_main:covariance_dimension_mismatch', ...
              'whitened_covariances{%d} must be %dx%d', f, p, p);
    end
end

% Smoothing kernel K (F x F)
K = input_data.smoothing_kernel;
if ~isnumeric(K) || ~isequal(size(K), [F,F])
    error('module4_objective_gradient_main:invalid_kernel', ...
          'smoothing_kernel must be a %dx%d numeric matrix', F, F);
end
if norm(K - K','fro') > 1e-12
    warning('module4_objective_gradient_main:kernel_not_symmetric', ...
            'smoothing_kernel not symmetric; symmetrizing');
    K = (K + K')/2;
end

% Weight matrix for cross-frequency smoothing (p x p)
Wg = input_data.weight_matrix;

% ------------ Read legacy/general params ------------
lambda1      = getf(gradient_params,'lambda1',0);         % cross-frequency
weight_mode  = lower(getf(gradient_params,'weight_mode','matrix')); % 'matrix'|'hadamard'
use_L        = getf(gradient_params,'use_graph_laplacian',false);
chol_tol     = getf(gradient_params,'chol_tolerance',1e-12); %#ok<NASGU>
sym_tol      = getf(gradient_params,'symmetrization_tolerance',1e-12);
force_H      = getf(gradient_params,'force_hermitian',true);
verbose      = getf(gradient_params,'verbose',false); %#ok<NASGU>

% Validate Wg under selected weight_mode
if ~isnumeric(Wg) || ~isequal(size(Wg), [p,p])
    error('module4_objective_gradient_main:invalid_weight_matrix', ...
          'weight_matrix must be a %dx%d numeric matrix', p, p);
end
if norm(Wg - Wg','fro') > 1e-12
    warning('module4_objective_gradient_main:weight_not_hermitian', ...
            'weight_matrix not Hermitian; symmetrizing');
    Wg = (Wg + Wg')/2;
end
switch weight_mode
    case 'matrix'
        minEigW = min(real(eig((Wg+Wg')/2)));
        if minEigW < -1e-12
            error('module4_objective_gradient_main:weight_not_psd', ...
                  'weight_matrix not PSD (min eig %.2e)', minEigW);
        end
    case 'hadamard'
        if any(~isfinite(Wg(:))) || ~isreal(Wg)
            error('module4_objective_gradient_main:weight_mask_invalid', ...
                  'Hadamard weight/mask must be real and finite.');
        end
    otherwise
        error('module4_objective_gradient_main:invalid_weight_mode', ...
              'weight_mode must be matrix|hadamard');
end

% ------------ NEW spatial smoothing params ------------
lambda3      = getf(gradient_params,'lambda3',0);
Gsp          = getf(gradient_params,'spatial_graph_matrix',[]);
isLap        = getf(gradient_params,'spatial_graph_is_laplacian',false);
sp_mode      = lower(getf(gradient_params,'spatial_weight_mode','node')); % 'node'|'hadamard'
use_spatial  = ~isempty(lambda3) && lambda3 > 0;

if use_spatial
    if isempty(Gsp) || ~isnumeric(Gsp) || ~isequal(size(Gsp), [p,p])
        error('module4_objective_gradient_main:invalid_spatial_graph', ...
              'spatial_graph_matrix must be a %dx%d numeric matrix when lambda3>0', p, p);
    end
    if norm(Gsp - Gsp','fro') > 1e-12
        warning('module4_objective_gradient_main:spatial_graph_not_symmetric', ...
                'spatial_graph_matrix not symmetric; symmetrizing');
        Gsp = (Gsp + Gsp')/2;
    end
end

% ==================== Allocate ====================
logdet_grads    = cell(F,1);
trace_grads     = cell(F,1);
smooth_grads    = cell(F,1);   % cross-frequency
spatial_grads   = [];          % optional
total_tic = tic; %#ok<NASGU>

% ==================== Per-frequency basic terms ====================
t1 = tic;
for f = 1:F
    Gf = Gammas{f};
    Gf = (Gf+Gf')/2;  %#ok<NASGU>
    logdet_grads{f} = -inv(Gf);   % -∂logdet/∂G = -G^{-1}
end
t_logdet = toc(t1);

t2 = tic;
for f = 1:F
    trace_grads{f} = (Sigmas{f}+Sigmas{f}')/2;  % ∂tr(SG)/∂G = S
end
t_trace = toc(t2);

% ==================== Cross-frequency smoothing ====================
t3 = tic;
if lambda1 > 0
    try
% function [smoothing_gradients, computation_stats] = module4_smoothing_gradient(precision_matrices, kernel_matrix, weight_matrix, params)
par = struct('lambda1', lambda1,...
    'mode', weight_mode);
        smooth_grads = module4_smoothing_gradient( ...
            Gammas, K, Wg, par);
    catch ME
        error('module4_objective_gradient_main:smoothing_failed', ...
              'Cross-frequency smoothing gradient failed: %s', ME.message);
    end
else
    for f = 1:F, smooth_grads{f} = zeros(p,p); end
end
t_smooth = toc(t3);

% ==================== NEW: Spatial single-frequency smoothing ====================
t_spatial = 0;
if use_spatial
    t4 = tic;
    try
        sp_opts = struct( ...
            'lambda3',                  lambda3, ...
            'spatial_graph_matrix',     Gsp, ...
            'spatial_graph_is_laplacian', isLap, ...
            'spatial_weight_mode',      sp_mode, ...
            'return_gradient',          true, ...
            'validate_inputs',          true, ...
            'enforce_hermitian_grad',   force_H);
        sp_out = module5_spatial_smoothing_singlefreq(Gammas, sp_opts);
    catch ME
        error('module4_objective_gradient_main:spatial_failed', ...
              'Spatial smoothing gradient failed: %s', ME.message);
    end
    t_spatial = toc(t4);

    if isfield(sp_out,'grad')
        spatial_grads = sp_out.grad;
    elseif isfield(sp_out,'gradients')
        spatial_grads = sp_out.gradients;
    else
        error('module4_objective_gradient_main:spatial_missing_grad', ...
              'module5_spatial_smoothing_singlefreq did not return gradients.');
    end

    if ~iscell(spatial_grads) || numel(spatial_grads) ~= F
        error('module4_objective_gradient_main:spatial_bad_shape', ...
              'Spatial gradients must be a cell array of length F.');
    end
end

% ==================== Combine & finalize ====================
t5 = tic;
total_grads = cell(F,1);
for f = 1:F
    Gsum = logdet_grads{f} + trace_grads{f} + smooth_grads{f};
    if use_spatial, Gsum = Gsum + spatial_grads{f}; end
    if force_H
        Gsum = (Gsum + Gsum')/2;
        if sym_tol > 0
            Gsum = (abs(Gsum) > sym_tol).*Gsum;
        end
    end
    total_grads{f} = Gsum;
end
t_combine = toc(t5);

% ==================== Pack results ====================
gradient_results = struct();
gradient_results.smooth_gradients = total_grads;

components = struct();
components.logdet_gradients    = logdet_grads;
components.trace_gradients     = trace_grads;
components.smoothing_gradients = smooth_grads;
if use_spatial, components.spatial_gradients = spatial_grads; end
gradient_results.gradient_components = components;

stats = struct();
stats.computation_times = [t_logdet, t_trace, t_smooth, t_combine];
if use_spatial, stats.spatial_grad_time = t_spatial; end
stats.symmetrization_tolerance = sym_tol;
stats.force_hermitian          = force_H;
gradient_results.computation_stats = stats;

gradient_results.success = true;

% ==================== helpers ====================
function v = getf(s, name, defv)
    if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = defv; end
end

end
