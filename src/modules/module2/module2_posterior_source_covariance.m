function Sigma_jj_post = module2_posterior_source_covariance(Sigma_jj_omega, L, Sigma_xi_xi, options)
% POSTERIOR: Sigma_post = (Sigma_jj^{-1} + L' * Sigma_xi^{-1} * L)^{-1}

if nargin < 4 || isempty(options), options = struct(); end
def = struct('regularization_factor',1e-8, ...
             'jitter_max_tries', 6, ...
             'ensure_positive_definite', true, ...
             'min_eigenvalue_ratio', 1e-12, ...
             'verbose', false);
fn = fieldnames(def);
for i=1:numel(fn), if ~isfield(options,fn{i}), options.(fn{i}) = def.(fn{i}); end, end

[n1,n2] = size(Sigma_jj_omega);
if n1~=n2, error('Sigma_jj_omega must be square'); end
n = n1;
[p, nL] = size(L);
if nL ~= n, error('Leadfield columns must match Sigma_jj size'); end
if ~isequal(size(Sigma_xi_xi), [p p]), error('Sigma_xi_xi must be p×p'); end

% Symmetrize inputs for numerical hygiene
Sigma_jj_omega = (Sigma_jj_omega + Sigma_jj_omega')/2;
Sigma_xi_xi    = (Sigma_xi_xi    + Sigma_xi_xi')/2;

% Information matrix: G = Σ_jj^{-1} + L' Σ_xi^{-1} L
B = Sigma_xi_xi \ L;          % p×n   (Σ_xi^{-1} L)
G = (Sigma_jj_omega \ eye(n)) + (L' * B);
G = (G + G')/2;

% Cholesky with jitter
[R, last_jitter, ok] = chol_with_jitter(G, options.regularization_factor, options.jitter_max_tries, options.verbose);
if ~ok
    error('module2_posterior_source_covariance:spd_factorization_failed', ...
          'Failed to factorize G even with jitter (last jitter=%.3e)', last_jitter);
end

% Inverse from Cholesky: G^{-1} = R \ (R' \ I)
Y = R' \ eye(n);
Sigma_jj_post = R \ Y;
Sigma_jj_post = (Sigma_jj_post + Sigma_jj_post')/2;

% Optional eigenvalue floor (usually unnecessary with info-form)
if options.ensure_positive_definite
    [V,D] = eig(Sigma_jj_post);
    d = real(diag(D));
    dmax = max(d);
    if isfinite(dmax) && ~isempty(dmax)
        floor_val = options.min_eigenvalue_ratio * max(dmax, eps);
        idx = d < floor_val;
        if any(idx)
            d(idx) = floor_val;
            Sigma_jj_post = V * diag(d) * V';
            Sigma_jj_post = (Sigma_jj_post + Sigma_jj_post')/2;
        end
    end
end

% Drop tiny imaginary leakage ONLY when all inputs are real
if isreal(L) && isreal(Sigma_jj_omega) && isreal(Sigma_xi_xi)
    if max(abs(imag(Sigma_jj_post(:)))) < 1e-13
        Sigma_jj_post = real(Sigma_jj_post);
    end
end
end

% ---- helper ----
function [R, last_jitter, ok] = chol_with_jitter(A, base_reg, max_tries, verbose)
% Cholesky with geometric jitter escalation.

A = (A + A')/2;
ok = false; last_jitter = 0;
n = size(A,1);
scale = trace(A)/max(1,n); if ~isfinite(scale) || scale <= 0, scale = 1; end
j0 = base_reg * scale;

for k = 0:max_tries
    j = j0 * (10^k);
    At = A + j * eye(n);
    [R, flag] = chol(At);
    if flag == 0
        ok = true; last_jitter = j;
        if verbose, fprintf('  [chol] success with jitter = %.3e (tries=%d)\n', j, k+1); end
        return;
    else
        if verbose, fprintf('  [chol] fail @ jitter = %.3e\n', j); end
    end
end
R = []; last_jitter = j0 * (10^max_tries);
end
