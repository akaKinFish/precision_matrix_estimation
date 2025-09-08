function [smoothing_gradients, computation_stats] = module4_smoothing_gradient(precision_matrices, kernel_matrix, weight_matrix, params)
% MODULE4_SMOOTHING_GRADIENT - Cross-frequency smoothing gradient
%
% Supports two weight semantics:
%   weight_mode = 'matrix'  :  ∇_Γw S = 2λ1 * W * ( Σ_{w'} L_{w,w'} Γ_{w'} )         (Laplacian)
%                              or 2λ1 Σ k_{w,w'} W (Γw - Γw')                         (pairwise)
%   weight_mode = 'hadamard':  ∇_Γw S = 2λ1 * (W.^2) .* ( Σ_{w'} L_{w,w'} Γ_{w'} )    (Laplacian)
%                              or 2λ1 Σ k_{w,w'} (W.^2) .* (Γw - Γw')                 (pairwise)
%
% Inputs
%   precision_matrices : cell{F,1}, each p×p
%   kernel_matrix      : F×F, frequency coupling K
%   weight_matrix      : p×p, W (no PSD requirement in hadamard mode)
%   params (struct)    :
%       .lambda1 (>=0)                 default 0.01
%       .use_graph_laplacian (bool)    default true
%       .weight_mode ('matrix'|'hadamard') default 'matrix'
%       .force_hermitian (bool)        default true
%       .symmetrization_tolerance      default 1e-10
%       .verbose (bool)                default false
%       .kernel_zero_tol               default 1e-12  (drop tiny K entries)
%
% Outputs
%   smoothing_gradients : cell{F,1}, each p×p
%   computation_stats   : struct with timings, method flags, etc.

% ---------- Input checks ----------
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
p = size(precision_matrices{1},1);

if ~isnumeric(kernel_matrix) || ~ismatrix(kernel_matrix) || any(size(kernel_matrix)~=[F F])
    error('module4_smoothing_gradient:invalid_kernel_matrix', ...
          'kernel_matrix must be %dx%d', F, F);
end
if ~isnumeric(weight_matrix) || ~ismatrix(weight_matrix) || any(size(weight_matrix)~=[p p])
    error('module4_smoothing_gradient:invalid_weight_matrix', ...
          'weight_matrix must be %dx%d', p, p);
end

% ---------- Defaults ----------
defaults.lambda1                 = 0.01;
defaults.use_graph_laplacian     = true;
defaults.weight_mode             = 'matrix';   % 'matrix' | 'hadamard'
defaults.force_hermitian         = true;
defaults.symmetrization_tolerance= 1e-10;
defaults.verbose                 = false;
defaults.kernel_zero_tol         = 1e-12;

fn = fieldnames(defaults);
for i = 1:numel(fn)
    f = fn{i};
    if ~isfield(params,f) || isempty(params.(f)), params.(f) = defaults.(f); end
end

% ---------- Init stats ----------
t_total = tic;
computation_stats = struct();
computation_stats.method_used = ternary(params.use_graph_laplacian,'laplacian','pairwise');
computation_stats.weight_mode = lower(params.weight_mode);
computation_stats.hermitian_violations = zeros(F,1);
computation_stats.laplacian_eigenvalues = [];

% ---------- Preprocess K and W ----------
K = (kernel_matrix + kernel_matrix')/2;
K(1:F+1:end) = 0;
K(abs(K) < params.kernel_zero_tol) = 0;

W = (weight_matrix + weight_matrix')/2;   % ensure Hermitian numerically
W2 = [];                                  % only used in hadamard

% ---------- Early exit if λ1=0 ----------
smoothing_gradients = cell(F,1);
if params.lambda1 == 0
    for w=1:F, smoothing_gradients{w} = zeros(p,p); end
    computation_stats.computation_time = toc(t_total);
    return;
end

% ---------- Compute gradient ----------
if params.use_graph_laplacian
    % Laplacian path
    d = sum(K,2);
    L = diag(d) - K;

    % Optional spectrum (cheap for small F)
    try
        computation_stats.laplacian_eigenvalues = eig(L);
    catch
        computation_stats.laplacian_eigenvalues = [];
    end

    switch lower(params.weight_mode)
        case 'matrix'
            for w = 1:F
                combo = zeros(p,p);
                for wp = 1:F
                    Lval = L(w,wp);
                    if Lval ~= 0
                        combo = combo + Lval * precision_matrices{wp};
                    end
                end
                smoothing_gradients{w} = 2 * params.lambda1 * (W * combo);
            end

        case 'hadamard'
            if isempty(W2), W2 = W.^2; end
            for w = 1:F
                combo = zeros(p,p);
                for wp = 1:F
                    Lval = L(w,wp);
                    if Lval ~= 0
                        combo = combo + Lval * precision_matrices{wp};
                    end
                end
                smoothing_gradients{w} = 2 * params.lambda1 * (W2 .* combo);
            end

        otherwise
            error('module4_smoothing_gradient:invalid_weight_mode','%s', params.weight_mode);
    end

else
    % Pairwise path (kept for compatibility)
    switch lower(params.weight_mode)
        case 'matrix'
            for w = 1:F
                Acc = zeros(p,p);
                for wp = 1:F
                    kwp = K(w,wp);
                    if kwp ~= 0
                        Acc = Acc + kwp * (precision_matrices{w} - precision_matrices{wp});
                    end
                end
                smoothing_gradients{w} = 2 * params.lambda1 * (W * Acc);
            end

        case 'hadamard'
            if isempty(W2), W2 = W.^2; end
            for w = 1:F
                Acc = zeros(p,p);
                for wp = 1:F
                    kwp = K(w,wp);
                    if kwp ~= 0
                        Acc = Acc + kwp * (precision_matrices{w} - precision_matrices{wp});
                    end
                end
                smoothing_gradients{w} = 2 * params.lambda1 * (W2 .* Acc);
            end

        otherwise
            error('module4_smoothing_gradient:invalid_weight_mode','%s', params.weight_mode);
    end
end

% ---------- Hermitian projection (for numerical symmetry) ----------
if params.force_hermitian
    for w = 1:F
        G = smoothing_gradients{w};
        Gherm = (G + G')/2;
        computation_stats.hermitian_violations(w) = norm(G - Gherm,'fro') / max(1, norm(G,'fro'));
        smoothing_gradients{w} = Gherm;
    end
end

computation_stats.computation_time = toc(t_total);
end

% ---- small ternary helper ----
function y = ternary(cond, a, b)
if cond, y = a; else, y = b; end
end
