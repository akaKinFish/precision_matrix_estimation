function out = module5_spatial_smoothing_singlefreq(Gamma_cells, opts)
% MODULE5_SPATIAL_SMOOTHING_SINGLEFREQ
% Spatial (within-frequency) smoothing penalty for precision matrices Γ_f.
%
% This module implements *single-frequency spatial smoothing* and its gradient,
% designed to complement the existing *cross-frequency* smoothing in Module 5.
% It does NOT change any existing interfaces; it is an auxiliary function that
% can be called by Module 4/5 to add a smooth term to the gradient/objective.
%
% ------------------------------- PURPOSE ----------------------------------
% Encourage spatial smoothness of Γ_f across the source graph:
%   - 'node' mode (recommended): row/column graph-Laplacian smoothing
%       R_node(Γ_f) = tr(Γ_f' L Γ_f) + tr(Γ_f L Γ_f')
%     where L is the node Laplacian built from the source adjacency.
%     This penalizes variations of *rows and columns* of Γ_f along the graph.
%
%   - 'hadamard' mode: element-wise weighted smoothing of the *L-filtered* Γ_f
%       S_f = L Γ_f + Γ_f L
%       R_hadamard(Γ_f) = || M ∘ S_f ||_F^2
%     where M is an element-wise weight mask (defaults to all-ones, optionally
%     excluding diagonal). This keeps compatibility with the "Hadamard semantics".
%
% ------------------------------- SIGNATURE --------------------------------
%   out = module5_spatial_smoothing_singlefreq(Gamma_cells, opts)
%
% INPUTS:
%   Gamma_cells  : {F×1} cell, each n×n complex/symmetric (whitened-domain Γ_f)
%   opts         : struct with fields (all optional unless stated):
%       .lambda3                      (double >= 0, default 0)
%       .spatial_graph_matrix         (n×n double)  % either adjacency A or Laplacian L
%       .spatial_graph_is_laplacian   (logical, default: true)
%       .spatial_weight_mode          ('node'|'hadamard', default: 'node')
%       .element_weight_mask          (n×n double, hadamard only; default ones(n))
%       .penalize_diagonal_elements   (logical, hadamard only; default: true)
%       .return_gradient              (logical, default: true)
%       .validate_inputs              (logical, default: true)
%       .enforce_hermitian_grad       (logical, default: true)
%
% OUTPUT:
%   out.term   : scalar, total spatial penalty  λ3 * Σ_f R(Γ_f)
%   out.grad   : {F×1} gradient cells (if return_gradient==true)
%   out.stats  : struct with diagnostic stats (PSD projections, symmetrization, etc.)
%
% USAGE (from Module 4/5):
%   sp = module5_spatial_smoothing_singlefreq(Gamma_cells, struct( ...
%         'lambda3', 0.2, ...
%         'spatial_graph_matrix', G.L, ...
%         'spatial_graph_is_laplacian', true, ...
%         'spatial_weight_mode', 'node', ...
%         'return_gradient', true));
%   total_term = sp.term;
%   total_grad = sp.grad;   % to be added to smooth gradients
%
% NOTES:
%   - If lambda3 == 0 or spatial_graph_matrix is empty, out.term = 0 and grad = zeros.
%   - For numerical robustness, we symmetrize L and (optionally) project to PSD.
%   - Complex matrices are handled in a real-valued objective via real(...) and
%     conjugate-transpose. We assume the spatial graph is real and symmetric.
%
% Author: your team
% Date: 2025-09-09
% Version: 1.0

% --------------------------- Parse & Validate ------------------------------
if nargin < 1
    error('m5_spatial:insufficient_input', ...
        'Gamma_cells is required and must be a {F×1} cell of n×n matrices.');
end
if ~iscell(Gamma_cells) || isempty(Gamma_cells)
    error('m5_spatial:invalid_gamma', 'Gamma_cells must be a non-empty cell.');
end
F = numel(Gamma_cells);
n = size(Gamma_cells{1},1);
for f = 1:F
    Gf = Gamma_cells{f};
    if ~ismatrix(Gf) || size(Gf,1)~=n || size(Gf,2)~=n
        error('m5_spatial:gamma_shape', ...
            'Gamma_cells{%d} must be %dx%d, got %s.', f, n, n, mat2str(size(Gf)));
    end
end

if nargin < 2, opts = struct(); end

% defaults
d.lambda3                    = 0;
d.spatial_graph_matrix       = [];
d.spatial_graph_is_laplacian = true;
d.spatial_weight_mode        = 'node';   % 'node' | 'hadamard'
d.element_weight_mask        = [];       % hadamard only; default ones(n)
d.penalize_diagonal_elements = true;     % hadamard only
d.return_gradient            = true;
d.validate_inputs            = true;
d.enforce_hermitian_grad     = true;

opts = set_defaults_(opts, d);

out = struct('term', 0, 'grad', [], 'stats', struct());
out.stats = struct('sym_L_enforced', false, ...
                   'psd_projection', false, ...
                   'min_eig_L_before', NaN, ...
                   'min_eig_L_after',  NaN);

if opts.lambda3 <= 0 || isempty(opts.spatial_graph_matrix)
    % disabled or no graph provided
    if opts.return_gradient
        out.grad = cell(F,1); for f=1:F, out.grad{f} = zeros(n,n,'like',Gamma_cells{f}); end
    end
    return;
end

% --------------------------- Prepare Graph L -------------------------------
L = opts.spatial_graph_matrix;
if ~ismatrix(L) || size(L,1)~=n || size(L,2)~=n
    error('m5_spatial:graph_shape', ...
        'spatial_graph_matrix must be %dx%d, got %s.', n, n, mat2str(size(L)));
end

% If adjacency supplied, convert to Laplacian
if ~opts.spatial_graph_is_laplacian
    % sanitize adjacency: real & symmetric & zero diagonal
    L = real(L);
    L = 0.5*(L+L');
    L(1:n+1:end) = 0;
    deg = sum(L,2);
    L = diag(deg) - L;  % unnormalized Laplacian
else
    % sanitize Laplacian
    L = real(L);
    L = 0.5*(L+L');
end

% Optional PSD projection (clip tiny negative eigenvalues due to numeric noise)
if opts.validate_inputs
    [V,D] = eig(L);
    evals = real(diag(D));
    out.stats.min_eig_L_before = min(evals);
    if min(evals) < 0
        D = diag(max(evals, 0));
        L = V*D*V';
        out.stats.psd_projection = true;
    end
end
out.stats.sym_L_enforced = true;
out.stats.min_eig_L_after = min(eig(0.5*(L+L')));

% --------------------------- Compute Term/Grad -----------------------------
lambda3 = opts.lambda3;
mode = lower(string(opts.spatial_weight_mode));

switch mode
    case "node"
        % R(Γ_f) = tr(Γ_f' L Γ_f) + tr(Γ_f L Γ_f')
        % For complex Γ_f, use real(...) & conjugate transpose.
        term_total = 0;
        grads = cell(F,1);
        for f=1:F
            G = Gamma_cells{f};
            % penalty value
            term_f = real(trace(G' * L * G)) + real(trace(G * L * G'));
            term_total = term_total + term_f;

            if opts.return_gradient
                % grad = d/dG [ tr(G' L G) + tr(G L G') ] = (L+L')G + G(L+L')
                % with symmetry, L=L', so grad = 2 L G + 2 G L
                Lt = 0.5*(L+L');  % ensure symmetric
                grad_f = 2*(Lt*G + G*Lt);
                if opts.enforce_hermitian_grad
                    grad_f = 0.5*(grad_f + grad_f');  % keep gradient Hermitian
                end
                grads{f} = lambda3 * grad_f;
            end
        end
        out.term = lambda3 * term_total;
        if opts.return_gradient, out.grad = grads; end

    case "hadamard"
        % S_f = L G + G L
        % R(Γ_f) = || M ∘ S_f ||_F^2, with M element-wise mask.
        if isempty(opts.element_weight_mask)
            M = ones(n,n);
        else
            M = opts.element_weight_mask;
            if ~isequal(size(M),[n n])
                error('m5_spatial:mask_shape','element_weight_mask must be %dx%d.', n,n);
            end
            M = real(M);
        end
        if ~opts.penalize_diagonal_elements
            M(1:n+1:end) = 0;
        end
        M2 = M.^2;

        term_total = 0;
        grads = cell(F,1);
        for f=1:F
            G = Gamma_cells{f};
            S = L*G + G*L;                    % "graph-filtered" Γ_f
            S = 0.5*(S + S');                 % symmetrize for stability
            Sf = M .* S;                      % element-wise mask
            % penalty value
            term_f = sum(sum( real(Sf).^2 + imag(Sf).^2 ));
            term_total = term_total + term_f;

            if opts.return_gradient
                % R = || M ∘ (L G + G L) ||_F^2
                % dR/dG = 2 [ L' ∘ I ]( M∘S ) + 2 [ I ∘ L' ]( M∘S )
                % Which simplifies (with symmetric L) to:
                % grad = 2 * ( L * (M2 ∘ S) + (M2 ∘ S) * L )
                T = M2 .* S;
                grad_f = 2*( L*T + T*L );
                if opts.enforce_hermitian_grad
                    grad_f = 0.5*(grad_f + grad_f');
                end
                grads{f} = lambda3 * grad_f;
            end
        end
        out.term = lambda3 * term_total;
        if opts.return_gradient, out.grad = grads; end

    otherwise
        error('m5_spatial:invalid_mode', ...
            'Unsupported spatial_weight_mode: %s. Use ''node'' or ''hadamard''.', opts.spatial_weight_mode);
end

end

% ------------------------------- helpers ----------------------------------
function S = set_defaults_(S, D)
    f = fieldnames(D);
    for i=1:numel(f)
        k = f{i};
        if ~isfield(S,k) || isempty(S.(k)), S.(k) = D.(k); end
    end
end
