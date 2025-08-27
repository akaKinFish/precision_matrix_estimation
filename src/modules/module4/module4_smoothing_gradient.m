function [smoothing_gradients, computation_stats] = module4_smoothing_gradient(precision_matrices, kernel_matrix, weight_matrix, params)
% MODULE4_SMOOTHING_GRADIENT - Cross-frequency smoothing gradient (revised)
%
% S = (λ1/2) * Σ_{ω,ω'} k_{ω,ω'} ||Γ_ω - Γ_{ω'}||^2_{W}
% ∇_{Γ_ω}S = 2λ1 Σ_{ω'} k_{ω,ω'} W(Γ_ω - Γ_{ω'})
% Laplacian form (equivalent):
%   ∇_{Γ_ω}S = 2λ1 W Σ_{ω'} L_{ω,ω'} Γ_{ω'}, where L = D - K.
%
% New option:
%   .kernel_zero_tol (double) default 1e-12   % drop tiny kernel entries
%
% Other params (same defaults as before):
%   .lambda1 (>=0), .use_graph_laplacian (true), .force_hermitian (true),
%   .symmetrization_tolerance (1e-10), .verbose (false)

% -------- Input checks --------
if nargin < 3
    error('module4_smoothing_gradient:insufficient_input', ...
          'precision_matrices, kernel_matrix, and weight_matrix are required');
end
if nargin < 4, params = struct(); end
if ~iscell(precision_matrices) || isempty(precision_matrices)
    error('module4_smoothing_gradient:invalid_precision_matrices', ...
          'precision_matrices must be a non-empty cell array');
end

F = numel(precision_matrices);
p = size(precision_matrices{1}, 1);

for f = 1:F
    if ~isnumeric(precision_matrices{f}) || ~isequal(size(precision_matrices{f}), [p, p])
        error('module4_smoothing_gradient:precision_dimension_mismatch', ...
              'precision_matrices{%d} must be %dx%d', f, p, p);
    end
end
if ~isnumeric(kernel_matrix) || ~isequal(size(kernel_matrix), [F, F])
    error('module4_smoothing_gradient:invalid_kernel_matrix', ...
          'kernel_matrix must be %dx%d numeric matrix', F, F);
end
if ~isnumeric(weight_matrix) || ~isequal(size(weight_matrix), [p, p])
    error('module4_smoothing_gradient:invalid_weight_matrix', ...
          'weight_matrix must be %dx%d numeric matrix', p, p);
end

% -------- Defaults --------
defaults = struct();
defaults.lambda1 = 0.01;
defaults.use_graph_laplacian = true;
defaults.force_hermitian = true;
defaults.symmetrization_tolerance = 1e-10;
defaults.verbose = false;
defaults.kernel_zero_tol = 1e-12;  % NEW

fn = fieldnames(defaults);
for i = 1:numel(fn)
    f = fn{i};
    if ~isfield(params, f), params.(f) = defaults.(f); end
end
if params.lambda1 < 0
    error('module4_smoothing_gradient:invalid_lambda1', ...
          'lambda1 must be non-negative, got %.6f', params.lambda1);
end

% -------- Init --------
smoothing_gradients = cell(F, 1);
computation_stats = struct();
computation_stats.computation_time = 0;
computation_stats.method_used = ternary(params.use_graph_laplacian,'laplacian','direct');
computation_stats.hermitian_violations = zeros(F,1); % relative error
computation_stats.sparsity_exploited = false;
computation_stats.laplacian_eigenvalues = [];

t_total = tic;

% λ1 == 0 -> zeros
if params.lambda1 == 0
    for f = 1:F, smoothing_gradients{f} = zeros(p,p); end
    computation_stats.computation_time = toc(t_total);
    if params.verbose, fprintf('  Smoothing: λ₁=0, returning zero gradients\n'); end
    return;
end

% -------- Preprocess K and W --------
K = kernel_matrix;
K = (K + K')/2;
K(1:F+1:end) = 0;                                  % zero diagonal
K(abs(K) < params.kernel_zero_tol) = 0;            % drop tiny entries

sparsity_ratio = nnz(K) / numel(K);
computation_stats.sparsity_exploited = sparsity_ratio < 0.5;

W = weight_matrix;
W = (W + W')/2;

% -------- Branch --------
if params.use_graph_laplacian
    [smoothing_gradients, lap_stats] = compute_smoothing_laplacian_method(precision_matrices, K, W, params);
    computation_stats.laplacian_eigenvalues = lap_stats.laplacian_eigenvalues;
else
    [smoothing_gradients] = compute_smoothing_direct_method(precision_matrices, K, W, params);
end

% -------- Hermitian projection & relative error --------
if params.force_hermitian
    for f = 1:F
        G = smoothing_gradients{f};
        Gherm = (G + G')/2;
        rel_err = norm(G - Gherm,'fro') / max(1, norm(G,'fro'));
        computation_stats.hermitian_violations(f) = rel_err;
        smoothing_gradients{f} = Gherm;
    end
end

computation_stats.computation_time = toc(t_total);
end

% -------- Helpers --------
function [gradients] = compute_smoothing_direct_method(Gammas, K, W, params)
F = numel(Gammas);
p = size(Gammas{1},1);
gradients = cell(F,1);

if params.verbose, fprintf('using direct method...\n'); end

for w = 1:F
    Gw = zeros(p,p);
    for wp = 1:F
        kwp = K(w,wp);
        if kwp ~= 0
            D = Gammas{w} - Gammas{wp};
            Gw = Gw + kwp * (W * D);
        end
    end
    gradients{w} = 2 * params.lambda1 * Gw;
end
end

function [gradients, stats] = compute_smoothing_laplacian_method(Gammas, K, W, params)
F = numel(Gammas);
p = size(Gammas{1},1);
gradients = cell(F,1);
stats = struct();

if params.verbose, fprintf('using Laplacian method...\n'); end

% L = D - K
d = sum(K,2);
L = diag(d) - K;

% spectrum (optional; cheap at small F)
try
    stats.laplacian_eigenvalues = eig(L);
catch
    stats.laplacian_eigenvalues = [];
end

% Precompute WG = W*Γ
WG = cell(F,1);
for w = 1:F
    WG{w} = W * Gammas{w};
end

for w = 1:F
    combo = zeros(p,p);
    for wp = 1:F
        Lval = L(w,wp);
        if Lval ~= 0
            combo = combo + Lval * Gammas{wp};
        end
    end
    gradients{w} = 2 * params.lambda1 * (W * combo);
end
end

function y = ternary(cond, a, b)
if cond, y = a; else, y = b; end
end
