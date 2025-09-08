function [alpha_auto, lambda1_auto, selection_stats] = module5_step_size_selection(Gamma_init, K_smooth, W_matrix, params)
% MODULE5_STEP_SIZE_SELECTION - Auto step size α and smoothing λ1 via bounds
%
% Supports:
%   params.weight_mode        : 'matrix' (legacy W*X) | 'hadamard' (W.*X)
%   params.use_graph_laplacian: true|false (use ||L||2 with L=diag(sum K)-K)
%   params.delta              : safety margin (default 0.9)
%   params.use_exact_spectral : exact min-eig (default false; else inverse-iteration)
%
% Bounds:
%   L_total <= L_logdet + L_smooth
%   L_logdet = max_ω ||Γ_ω^{-1}||_2^2
%   If use_graph_laplacian:
%       L_smooth <= 2 * λ1 * ||L||_2 * c_W
%   Else:
%       L_smooth <= 2 * λ1 * K_max * c_W
%   c_W =
%     - matrix   : ||W||_2                (spectral norm)
%     - hadamard : max(abs(W(:)))         (since ∇ uses (W.^2)∘X)
%
% We set:
%   λ1 = δ / (2 * freq_factor * c_W)
%   α  = 1 / (L_logdet + δ)

% -------- Input checks --------
if nargin < 3
    error('module5_step_size_selection:insufficient_input', ...
          'Gamma_init, K_smooth, and W_matrix are required');
end
if nargin < 4, params = struct(); end

% Defaults
def.delta               = 0.9;
def.power_iter_tol      = 1e-6;
def.power_iter_max      = 50;
def.use_exact_spectral  = false;
def.verbose             = false;
def.weight_mode         = 'matrix';     % 'matrix' | 'hadamard'
def.use_graph_laplacian = true;
fn = fieldnames(def);
for i=1:numel(fn)
    if ~isfield(params,fn{i}) || isempty(params.(fn{i})), params.(fn{i}) = def.(fn{i}); end
end

% Dimensions
F = numel(Gamma_init);
p = size(Gamma_init{1},1);

% --- K_smooth ---
if ~isnumeric(K_smooth) || ~ismatrix(K_smooth) || any(size(K_smooth)~=[F F])
    error('module5_step_size_selection:kernel_size_mismatch', ...
          'K_smooth must be %dx%d, got %dx%d.', F, F, size(K_smooth,1), size(K_smooth,2));
end
% --- W_matrix ---
if ~isnumeric(W_matrix) || ~ismatrix(W_matrix) || any(size(W_matrix)~=[p p])
    error('module5_step_size_selection:weight_size_mismatch', ...
          'W_matrix must be %dx%d, got %dx%d.', p, p, size(W_matrix,1), size(W_matrix,2));
end

selection_stats = struct();

% -------- L_logdet = max ||Γ^{-1}||_2^2 --------
L_logdet = 0;
spec_inv = zeros(F,1);
for f=1:F
    G = (Gamma_init{f}+Gamma_init{f}')/2;
    [~,flag] = chol(G);
    if flag~=0
        error('module5_step_size_selection:not_psd','Gamma_init{%d} is not PD', f);
    end
    if params.use_exact_spectral
        lam_min = min(real(eig(G)));
        sni = 1 / max(lam_min, 1e-15);
    else
        sni = estimate_spec_norm_inv(G, params.power_iter_tol, params.power_iter_max);
    end
    spec_inv(f) = sni;
    L_logdet = max(L_logdet, sni^2);
end
selection_stats.L_logdet        = L_logdet;
selection_stats.spectral_norms   = spec_inv;
selection_stats.max_spec_inv     = max(spec_inv);

% -------- Frequency factor: ||L||_2 (if Laplacian) else K_max --------
K = (K_smooth + K_smooth')/2;
if params.use_graph_laplacian
    L = diag(sum(K,2)) - K;
    spec_L = norm(L,2);
    selection_stats.spec_L    = spec_L;
    selection_stats.K_max     = max(sum(abs(K),2)); % keep for compatibility
    freq_factor               = spec_L;
else
    K_row_sums = sum(abs(K),2);
    K_max      = max(K_row_sums);
    selection_stats.K_max     = K_max;
    selection_stats.spec_L    = [];                 % not used
    freq_factor               = K_max;
end
selection_stats.freq_factor = freq_factor;

% -------- Weight factor & "R_max" (compatibility) --------
switch lower(params.weight_mode)
    case 'matrix'
        % true factor used in bound
        cW = norm((W_matrix+W_matrix')/2, 2);
        % compatibility: old code sometimes inspects Gershgorin row-sum as R_max
        R_max = max(sum(abs(W_matrix),2));
    case 'hadamard'
        % gradient contains (W.^2)∘X -> factor is max|W|
        cW = max(abs(W_matrix(:)));
        % compatibility: store R_max = max((W.^2)(:)) to mirror old naming
        R_max = max((abs(W_matrix(:))).^2);
    otherwise
        error('module5_step_size_selection:invalid_weight_mode','%s', params.weight_mode);
end
selection_stats.weight_mode    = params.weight_mode;
selection_stats.weight_factor  = cW;
selection_stats.R_max          = R_max;    % for backward compatibility

% -------- λ1 and α --------
if freq_factor > 0 && cW > 0
    lambda1_auto = params.delta / (2 * freq_factor * cW);
else
    lambda1_auto = params.delta * 0.01; % conservative fallback
    warning('module5_step_size_selection:degenerate_bounds','freq_factor or weight_factor is zero');
end
alpha_auto = 1 / (L_logdet + params.delta);

selection_stats.lambda1_theory = lambda1_auto;
selection_stats.alpha_theory   = alpha_auto;

if params.verbose
    fprintf('\n--- Auto selection ---\n');
    fprintf('weight_mode=%s | use_laplacian=%d\n', params.weight_mode, params.use_graph_laplacian);
    fprintf('L_logdet=%.3e | freq_factor=%.3e | weight_factor=%.3e\n', L_logdet, freq_factor, cW);
    fprintf('=> lambda1=%.3e | alpha=%.3e\n', lambda1_auto, alpha_auto);
end
end

% -------- Helper: estimate ||G^{-1}||_2 via inverse iteration --------
function sni = estimate_spec_norm_inv(G, tol, itmax)
p = size(G,1);
v = randn(p,1) + 1i*randn(p,1);
v = v / norm(v);
prev = Inf; lam = real(v' * G * v);
for t=1:itmax
    w = G \ v;  w = w / norm(w);
    lam_new = real(w' * G * w);
    if abs(lam_new - prev) <= tol * max(1,abs(lam_new)), break; end
    prev = lam_new; v = w; lam = lam_new;
end
lam = max(real(lam), 1e-15);
sni = 1 / lam;
end
