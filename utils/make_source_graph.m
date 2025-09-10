function G = make_source_graph(src_positions, varargin)
% MAKE_SOURCE_GRAPH
% Build a source-space graph (adjacency A and Laplacian L) from 3D positions.
%
% ------------------------------- PURPOSE ----------------------------------
% Construct a *node graph* over sources (n nodes) using either:
%   - kNN: connect each node to its k nearest neighbors
%   - radius: connect nodes within a Euclidean distance threshold r
% Optional Gaussian weighting on edges: w_ij = exp(-||xi-xj||^2 / (2*sigma^2))
%
% The output provides:
%   G.A   : n×n symmetric adjacency (zero diagonal)
%   G.L   : n×n symmetric PSD graph Laplacian (L = D - A)
%   G.info: metadata (method, parameters, eigen-stats, symmetry errors, etc.)
%
% ------------------------------- SIGNATURE --------------------------------
%   G = make_source_graph(src_positions, 'Name', Value, ...)
%
% REQUIRED:
%   src_positions : n×d numeric (usually d=3). Each row is the position of a source.
%
% NAME-VALUE PAIRS (all optional):
%   'method'     : 'knn' | 'radius'               (default: 'knn')
%   'k'          : positive integer (for knn)     (default: max(1, round(log(n))))
%   'radius'     : positive double (for radius)   (default: [])
%   'weighted'   : logical                        (default: true)
%   'sigma'      : positive double (Gaussian width for weights; default: median nn-dist)
%   'make_symmetric' : logical                    (default: true)
%   'normalize'  : 'none'|'row'|'sym' (adjacency normalization; default: 'none')
%   'laplacian'  : 'unnormalized'|'sym'|'rw'      (default: 'unnormalized')
%   'verbose'    : logical                        (default: true)
%
% OUTPUT:
%   G.A, G.L, G.info
%
% EXAMPLE:
%   pos = randn(64,3);
%   G = make_source_graph(pos, 'method','knn','k',6,'weighted',true,'normalize','none');
%   % Use G.L in spatial smoothing (node mode):
%   % out = module5_spatial_smoothing_singlefreq(Gamma_cells, struct('lambda3',0.3,'spatial_graph_matrix',G.L,'spatial_graph_is_laplacian',true,'spatial_weight_mode','node'));
%
% Author: your team
% Date: 2025-09-09
% Version: 1.0

p = inputParser;
p.addRequired('src_positions', @(X) isnumeric(X) && ismatrix(X) && size(X,1)>=2);
p.addParameter('method', 'knn', @(s) any(strcmpi(s,{'knn','radius'})));
p.addParameter('k', [], @(x) isempty(x) || (isscalar(x) && x>0 && isfinite(x)));
p.addParameter('radius', [], @(x) isempty(x) || (isscalar(x) && x>0 && isfinite(x)));
p.addParameter('weighted', true, @islogical);
p.addParameter('sigma', [], @(x) isempty(x) || (isscalar(x) && x>0));
p.addParameter('make_symmetric', true, @islogical);
p.addParameter('normalize', 'none', @(s) any(strcmpi(s,{'none','row','sym'})));
p.addParameter('laplacian', 'unnormalized', @(s) any(strcmpi(s,{'unnormalized','sym','rw'})));
p.addParameter('verbose', true, @islogical);
p.parse(src_positions, varargin{:});
opt = p.Results;

X = double(src_positions);
[n, d] = size(X); %#ok<ASGLU>
if isempty(opt.k), opt.k = max(1, round(log(n))); end

% -------- pairwise distances (squared) --------
% robust & vectorized (avoid pdist2 if not available)
D2 = squareform_fast_(X);

% -------- build adjacency A --------
A = zeros(n);
switch lower(opt.method)
    case 'knn'
        % For each node, connect to k nearest (excluding self)
        for i=1:n
            di = D2(i,:);
            di(i) = inf;
            [~, idx] = sort(di, 'ascend');
            neigh = idx(1:min(opt.k, n-1));
            A(i, neigh) = 1;
        end
        if opt.make_symmetric
            A = max(A, A');  % mutualize
        end

    case 'radius'
        if isempty(opt.radius)
            error('make_source_graph:radius_required', ...
                'For method=''radius'', you must specify a positive ''radius''.');
        end
        A = (D2 <= (opt.radius^2)) - eye(n);
        A = double(A>0);
        if opt.make_symmetric
            A = max(A, A');
        end

    otherwise
        error('make_source_graph:invalid_method', 'Unknown method: %s', opt.method);
end

% -------- weights (Gaussian) --------
if opt.weighted
    % choose sigma if not provided: median of nonzero nn distances
    if isempty(opt.sigma)
        nn = zeros(n,1);
        for i=1:n
            di = D2(i,:);
            di(i) = []; di = sort(di,'ascend');
            nn(i) = sqrt(di(1)); %#ok<AGROW>
        end
        opt.sigma = median(nn(nn>0));
        if ~isfinite(opt.sigma) || opt.sigma<=0
            opt.sigma = 1.0;
        end
    end
    % Gaussian weights on existing edges
    W = exp(-D2/(2*opt.sigma^2));
    A = A .* W;
end

% -------- normalize adjacency (optional) --------
switch lower(opt.normalize)
    case 'row'
        rs = sum(A,2) + eps;
        A = A ./ rs;
    case 'sym'
        rs = sum(A,2) + eps;
        Dn = diag(1./sqrt(rs));
        A = Dn * A * Dn;
    case 'none'
        % do nothing
    otherwise
        error('make_source_graph:invalid_normalize','Unknown normalization: %s', opt.normalize);
end

% -------- sanitize: zero diagonal, symmetrize, real --------
A(1:n+1:end) = 0;
A = real(0.5*(A + A'));

% -------- Laplacian --------
deg = sum(A,2);
L = diag(deg) - A;

switch lower(opt.laplacian)
    case 'sym'
        % L_sym = I - D^{-1/2} A D^{-1/2}
        Dn = diag(1./sqrt(max(deg,eps)));
        L = eye(n) - Dn * A * Dn;
    case 'rw'
        % L_rw = I - D^{-1} A
        Dinv = diag(1./max(deg,eps));
        L = eye(n) - Dinv * A;
    case 'unnormalized'
        % already set
    otherwise
        error('make_source_graph:invalid_laplacian','Unknown laplacian: %s', opt.laplacian);
end

% -------- diagnostics --------
% symmetry error, min/max eigenvalues, spectral radius
symerrA = norm(A - A','fro') / max(norm(A,'fro'), eps);
symerrL = norm(L - L','fro') / max(norm(L,'fro'), eps);
eigL = eig(0.5*(L+L'));
mineig = min(real(eigL));
maxeig = max(real(eigL));

G = struct();
G.A = A;
G.L = 0.5*(L + L');   % enforce symmetry
G.info = struct( ...
    'method', lower(opt.method), ...
    'k', opt.k, ...
    'radius', opt.radius, ...
    'weighted', opt.weighted, ...
    'sigma', opt.sigma, ...
    'normalize', lower(opt.normalize), ...
    'laplacian', lower(opt.laplacian), ...
    'symmetry_error_A', symerrA, ...
    'symmetry_error_L', symerrL, ...
    'min_eig_L', mineig, ...
    'max_eig_L', maxeig, ...
    'n_nodes', n, ...
    'dim', size(X,2) ...
);

if opt.verbose
    fprintf('[make_source_graph] n=%d | method=%s | k=%s | radius=%s | weighted=%d | sigma=%.3g\n', ...
        n, lower(opt.method), num2str(opt.k), num2str(opt.radius), opt.weighted, opt.sigma);
    fprintf('  normalize=%s | laplacian=%s | min_eig(L)=%.2e | max_eig(L)=%.2e | symerr(L)=%.2e\n', ...
        lower(opt.normalize), lower(opt.laplacian), mineig, maxeig, symerrL);
end

end

% ------------------------------- helpers ----------------------------------
function D2 = squareform_fast_(X)
% Pairwise squared Euclidean distances (symmetric, zero diagonal).
n = size(X,1);
XX = sum(X.^2,2);
D2 = max(XX + XX' - 2*(X*X'), 0);
D2(1:n+1:end) = 0;
D2 = real(0.5*(D2 + D2'));  % enforce symmetry
end
