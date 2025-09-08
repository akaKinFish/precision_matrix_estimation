function T_jv = module2_dstf_computation(L, Sigma_jj_omega, Sigma_xi_xi, options)
% Kept for compatibility (NOT used by main path now).

if nargin < 4 || isempty(options), options = struct(); end
def = struct('regularization_factor',1e-8,'jitter_max_tries',6, ...
             'method','information','use_pseudoinverse',false,'verbose',false);
fn = fieldnames(def);
for i=1:numel(fn), if ~isfield(options,fn{i}), options.(fn{i}) = def.(fn{i}); end, end

[p, n] = size(L);
Sigma_jj_omega = (Sigma_jj_omega + Sigma_jj_omega')/2;
Sigma_xi_xi    = (Sigma_xi_xi    + Sigma_xi_xi')/2;

% Information form: G * X = L' Σ_xi^{-1}
B = Sigma_xi_xi \ L;          % p×n
G = (Sigma_jj_omega \ eye(n)) + (L' * B);
G = (G + G')/2;
[X, ok] = solve_spd_with_jitter(G, B', options.regularization_factor, options.jitter_max_tries, options.verbose);
if ~ok
    if options.use_pseudoinverse, X = pinv(G) * B';
    else, error('module2_dstf_computation:spd_solve_failed','G factorization failed'); end
end
T_jv = X; % n×p

% Drop tiny imaginary leakage ONLY when all inputs are real
if isreal(L) && isreal(Sigma_jj_omega) && isreal(Sigma_xi_xi)
    if max(abs(imag(T_jv(:)))) < 1e-13, T_jv = real(T_jv); end
end
end

% ---- local helper ----
function [X, ok] = solve_spd_with_jitter(A, B, base_reg, max_tries, verbose)
A = (A + A')/2;
ok = false; X = [];
n = size(A,1);
scale = trace(A)/max(1,n); if ~isfinite(scale) || scale<=0, scale=1; end
j0 = base_reg * scale;
for k = 0:max_tries
    j = j0 * (10^k);
    At = A + j * eye(n);
    [R, flag] = chol(At);
    if flag == 0
        Y = R' \ B; X = R \ Y; ok = true;
        if verbose, fprintf('  [chol] success with jitter = %.3e (tries=%d)\n', j, k+1); end
        return;
    else
        if verbose, fprintf('  [chol] fail @ jitter = %.3e\n', j); end
    end
end
end
