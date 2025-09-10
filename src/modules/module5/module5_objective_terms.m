function [loglik_term, smooth_term, l1_term, aux] = module5_objective_terms(Gamma_cells, Sigma_cells, aux, params)
% MODULE5_OBJECTIVE_TERMS - Decompose objective into loglik/smooth/l1
%
% Smoothing term semantics (IMPORTANT / unified):
%   We implement cross-frequency smoothing via the **pairwise equivalent** form:
%       0.5 * sum_{f1,f2} K(f1,f2) * Phi( Gamma_{f1} - Gamma_{f2} )
%   where Phi depends on weight_mode:
%       - 'matrix'   : Phi(D) = tr( D' * W * D ), with W PSD (matrix quadratic)
%       - 'hadamard' : Phi(D) = || W ∘ D ||_F^2, with element-wise mask W
%   The factor 0.5 avoids double-counting. This pairwise form is **equivalent**
%   to a graph-Laplacian smoothing when K is interpreted as an adjacency/kernel.
%   Therefore, even when params.use_graph_laplacian = true, the code below uses
%   this pairwise implementation for clarity and consistency with gradients.
%
% Note:
%   This function also optionally computes a single-frequency spatial smoothing
%   *diagnostic term* (spatial_term) without changing the return signature. If
%   params.lambda3 > 0 and a spatial graph is provided, we call
%   module5_spatial_smoothing_singlefreq(...,'return_gradient',false) at the end.

if nargin<4, params = struct(); end
F = numel(Gamma_cells); p = size(Gamma_cells{1},1); %#ok<NASGU>

% ---------------- defaults (existing) ----------------
d.lambda1 = getfield_with_default(aux,'lambda1',0.01);
d.lambda2 = getfield_with_default(aux,'lambda2',1e-3);
d.weight_mode = getfield_with_default(params,'weight_mode','matrix');
d.use_graph_laplacian = getfield_with_default(params,'use_graph_laplacian',true);

% ---------------- NEW: spatial-smoothing related optional params ----------
d.lambda3 = getfield_with_default(params, 'lambda3', 0);                % strength for single-frequency spatial smoothing
d.spatial_graph_matrix = getfield_with_default(params, 'spatial_graph_matrix', []);  % n×n A or L
d.spatial_graph_is_laplacian = getfield_with_default(params, 'spatial_graph_is_laplacian', true);
d.spatial_weight_mode = getfield_with_default(params, 'spatial_weight_mode', 'node'); % 'node' | 'hadamard'

% ---------------- NEW: L1 diagonal control --------------------------------
% If true, include diagonal into L1 penalty; default false keeps legacy behavior.
d.penalize_diagonal = getfield_with_default(params, 'penalize_diagonal', false);

params = set_defaults(params, d);

% Read K and W from aux, and ensure symmetry; **do not** force zero diagonal
K = aux.smoothing_kernel;
W = aux.weight_matrix;
K = (K+K')/2;             % symmetric F×F
W = (W+W')/2;             % Hermitian n×n (PSD if 'matrix' mode)

% ---- 1) loglik ----
loglik_term = 0;
for f=1:F
    G = Gamma_cells{f};
    S = Sigma_cells{f};
    [L,flag] = chol((G+G')/2,'lower');
    if flag~=0, error('objective_terms', 'objective_terms:Gamma_not_PD'); end
    log_det = 2*sum(log(real(diag(L))));
    loglik_term = loglik_term + (-log_det + real(trace(S*G)));
end

% ---- 2) smoothing (cross-frequency; pairwise unified form) ----
smooth_term = 0;
if params.lambda1 > 0
    switch lower(params.weight_mode)
        case 'matrix'
            % pairwise equivalent (0.5 to avoid double counting)
            for f1=1:F
                for f2=1:F
                    if K(f1,f2)~=0
                        D = Gamma_cells{f1} - Gamma_cells{f2};
                        v = real(trace((D') * W * D));
                        smooth_term = smooth_term + 0.5 * K(f1,f2) * v;
                    end
                end
            end
        case 'hadamard'
            W2 = W.^2;
            for f1=1:F
                for f2=1:F
                    if K(f1,f2)~=0
                        D = Gamma_cells{f1} - Gamma_cells{f2};
                        v = sum(sum( W2 .* (real(D).^2 + imag(D).^2) )); % ||W∘D||_F^2
                        smooth_term = smooth_term + 0.5 * K(f1,f2) * v;
                    end
                end
            end
        otherwise
            error('objective_terms','objective_terms:invalid_weight_mode');
    end
    smooth_term = params.lambda1 * smooth_term;
end

% ---- 3) l1 ----
l1_term = 0;
if d.lambda2 > 0
    for f=1:F
        G = Gamma_cells{f};
        % off-diagonal (legacy)
        l1_off = 2*sum(sum(abs(triu(G,1))));
        % (NEW) diagonal optional
        if params.penalize_diagonal
            l1_diag = sum(abs(diag(G)));
        else
            l1_diag = 0;
        end
        l1_term = l1_term + (l1_off + l1_diag);
    end
    l1_term = d.lambda2 * l1_term;
end

% ---- 4) (NEW) single-frequency spatial smoothing *diagnostic* term -------
% Not returned; computed only for logging or external diagnostics.
spatial_term = 0; %#ok<NASGU>
if params.lambda3 > 0 && ~isempty(params.spatial_graph_matrix)
    try
        sp_out = module5_spatial_smoothing_singlefreq(Gamma_cells, struct( ...
            'lambda3', params.lambda3, ...
            'spatial_graph_matrix', params.spatial_graph_matrix, ...
            'spatial_graph_is_laplacian', params.spatial_graph_is_laplacian, ...
            'spatial_weight_mode', params.spatial_weight_mode, ...
            'return_gradient', false, ...
            'validate_inputs', true, ...
            'enforce_hermitian_grad', false ...
        ));
        spatial_term = sp_out.term; %#ok<NASGU>
        % assignin('caller','spatial_term_from_objective_terms',spatial_term); % optional
    catch ME
        warning('module5_objective_terms:spatial_diag_failed', ...
            'Spatial smoothing diagnostic failed: %s', ME.message);
    end
end
end

% --- helpers ---
function v = getfield_with_default(S, name, dv)
    if isfield(S,name) && ~isempty(S.(name)), v = S.(name); else, v = dv; end
end
function S = set_defaults(S, D)
    ff = fieldnames(D);
    for ii=1:numel(ff)
        if ~isfield(S,ff{ii}) || isempty(S.(ff{ii})), S.(ff{ii}) = D.(ff{ii}); end
    end
end
