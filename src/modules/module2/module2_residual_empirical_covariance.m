function S_xi_xi = module2_residual_empirical_covariance(T_xi_v, S_vv, options)
% S_xi_xi = T_xi_v * S_vv * T_xi_v^H
% Keeps complex Hermitian if S_vv is complex; only collapses to real when both inputs are real.

if nargin < 3 || isempty(options), options = struct(); end
def = struct('enforce_hermitian', true, ...
             'regularization_factor', 0, ...
             'min_eigenvalue_threshold', 1e-12, ...
             'psd_only', true, ...
             'numerical_tolerance', 1e-12, ...
             'verbose', false);
fn = fieldnames(def);
for i=1:numel(fn), if ~isfield(options,fn{i}), options.(fn{i}) = def.(fn{i}); end, end

[p1,p2] = size(T_xi_v);
if p1 ~= p2, error('T_xi_v must be square'); end
p = p1;
if ~isequal(size(S_vv), [p p]), error('S_vv must be %dx%d', p, p); end

% Optional ridge on S_vv
if options.regularization_factor > 0
    lam = options.regularization_factor * trace((S_vv+S_vv')/2) / max(1,p);
    S_vv = S_vv + lam * eye(p);
end

% Core triple product
tmp     = S_vv * T_xi_v';
S_xi_xi = T_xi_v * tmp;

% Hermitianize (does not force real; keeps complex off-diagonals)
if options.enforce_hermitian
    S_xi_xi = (S_xi_xi + S_xi_xi')/2;
end

% PSD/PD projection
try
    [V,D] = eig(S_xi_xi);
    d = real(diag(D));
    neg_tol = 10 * options.numerical_tolerance;
    d(d < 0 & d > -neg_tol) = 0;
    if options.psd_only
        d(d < 0) = 0;
    else
        floor_val = options.min_eigenvalue_threshold;
        if floor_val > 0, d(d < floor_val) = floor_val; end
    end
    S_xi_xi = V * diag(d) * V';
    if options.enforce_hermitian, S_xi_xi = (S_xi_xi + S_xi_xi')/2; end
catch
    % eigen issues: keep the unprojected Hermitian result
end

% IMPORTANT: cast to real ONLY when BOTH inputs are real and imag leakage is tiny
if isreal(T_xi_v) && isreal(S_vv)
    if max(abs(imag(S_xi_xi(:)))) < 1e-13
        S_xi_xi = real(S_xi_xi);
    end
end
end
