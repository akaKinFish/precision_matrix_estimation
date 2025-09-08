function objective_value = module5_objective(Gamma_cells, Sigma_tilde, aux_data, params)
% MODULE5_OBJECTIVE - Evaluate complete objective function for monitoring
%
% Supports two weight semantics:
%   - 'matrix'  : W * X  (legacy, W must be PSD for theory)
%   - 'hadamard': W .* X (element-wise mask; W can be non-PSD)
%
% Smoothing term is implemented in Laplacian form by default:
%   F_smooth = λ1 * sum_{f,f'} < W∘Γ_f , L_ff' (W∘Γ_{f'}) >
% where ∘ is either matrix-multiply or Hadamard depending on weight_mode,
% and L = diag(sum(K,2)) - K is the frequency Laplacian.
%
% Parameters (params):
%   .mode               : 'simplified'|'joint' (kept for compatibility)
%   .penalize_diagonal  : include diagonals in L1 (default: false)
%   .safe_logdet        : use stable Cholesky-based logdet (default: true)
%   .weight_mode        : 'matrix'|'hadamard' (default: 'matrix')
%   .use_graph_laplacian: true|false (default: true)
%
% Aux data (aux_data):
%   .smoothing_kernel (F×F), .weight_matrix (p×p),
%   .lambda1 (>=0), .lambda2 (>=0)
%

% ---------------- Input checks ----------------
if nargin < 3
    error('module5_objective:insufficient_input', ...
        'Gamma_cells, Sigma_tilde, and aux_data are required');
end
if nargin < 4, params = struct(); end

if ~isfield(params, 'mode'), params.mode = 'simplified'; end
if ~isfield(params, 'penalize_diagonal'), params.penalize_diagonal = false; end
if ~isfield(params, 'safe_logdet'), params.safe_logdet = true; end
if ~isfield(params, 'weight_mode') || isempty(params.weight_mode)
    params.weight_mode = 'matrix';
end
if ~isfield(params, 'use_graph_laplacian') || isempty(params.use_graph_laplacian)
    params.use_graph_laplacian = true;
end

F = numel(Gamma_cells);
if numel(Sigma_tilde) ~= F
    error('module5_objective:length_mismatch', ...
        'Gamma_cells and Sigma_tilde must have same length');
end
p = size(Gamma_cells{1}, 1);

K_smooth = aux_data.smoothing_kernel;
W_matrix = aux_data.weight_matrix;
lambda1  = aux_data.lambda1;
lambda2  = aux_data.lambda2;

% ---------------- Term 1: -logdet + trace(ΣΓ) ----------------
loglik_term = 0;
for f = 1:F
    Gf = Gamma_cells{f};
    Sf = Sigma_tilde{f};

    % PSD check (for safe log-det); if fails, objective is Inf
    [~, chol_flag] = chol((Gf+Gf')/2, 'lower');
    if chol_flag ~= 0
        warning('module5_objective:not_psd','Gamma{%d} not PD', f);
        objective_value = Inf; return;
    end

    if params.safe_logdet
        L = chol((Gf+Gf')/2, 'lower');
        log_det = 2 * sum(log(real(diag(L))));
    else
        log_det = log(det(Gf));
    end
    if ~isfinite(log_det)
        warning('module5_objective:infinite_logdet','Logdet not finite at f=%d', f);
        objective_value = Inf; return;
    end

    loglik_term = loglik_term + (-log_det + real(trace(Sf * Gf)));
end

% ---------------- Term 2: Cross-frequency smoothing ----------------
smoothing_term = 0;
if lambda1 > 0
    K = (K_smooth + K_smooth')/2;     % symmetrize for safety
    if params.use_graph_laplacian
        d = sum(K,2);
        L = diag(d) - K;              % Laplacian
        switch lower(params.weight_mode)
            case 'matrix'
                WG = cellfun(@(G) W_matrix * G, Gamma_cells, 'uni', 0);
                for f = 1:F
                    for fp = 1:F
                        if L(f,fp) ~= 0
                            smoothing_term = smoothing_term + ...
                                lambda1 * real(trace( WG{f}' * (L(f,fp)*WG{fp}) ));
                        end
                    end
                end
            case 'hadamard'
                WG = cellfun(@(G) W_matrix .* G, Gamma_cells, 'uni', 0);
                for f = 1:F
                    for fp = 1:F
                        if L(f,fp) ~= 0
                            smoothing_term = smoothing_term + ...
                                lambda1 * real( trace( WG{f}' * (L(f,fp) * WG{fp}) ) );
                        end
                    end
                end

            otherwise
                error('module5_objective:invalid_weight_mode','%s', params.weight_mode);
        end
    else
        % Pairwise form (keep compatibility). Use 0.5 to avoid double-counting.
        for f1 = 1:F
            for f2 = 1:F
                if K(f1,f2) ~= 0
                    D = Gamma_cells{f1} - Gamma_cells{f2};
                    switch lower(params.weight_mode)
                        case 'matrix'
                            wnorm2 = real(trace(D' * W_matrix * D));
                        case 'hadamard'
                            wnorm2 = sum(sum( abs(W_matrix .* D).^2 ));
                        otherwise
                            error('module5_objective:invalid_weight_mode','%s', params.weight_mode);
                    end
                    smoothing_term = smoothing_term + 0.5 * K(f1,f2) * wnorm2;
                end
            end
        end
        smoothing_term = lambda1 * smoothing_term;
    end
end

% ---------------- Term 3: L1 sparsity ----------------
l1_term = 0;
if lambda2 > 0
    for f = 1:F
        Gf = Gamma_cells{f};
        if params.penalize_diagonal
            l1_term = l1_term + 2*sum(sum(abs(triu(Gf,1)))) + sum(abs(real(diag(Gf))));
        else
            l1_term = l1_term + 2*sum(sum(abs(triu(Gf,1))));
        end
    end
    l1_term = lambda2 * l1_term;
end

% ---------------- Combine ----------------
objective_value = loglik_term + smoothing_term + l1_term;

if ~isfinite(objective_value)
    warning('module5_objective:non_finite_objective', ...
        'Objective not finite (loglik=%.2e, smooth=%.2e, l1=%.2e)', ...
        loglik_term, smoothing_term, l1_term);
end
end
