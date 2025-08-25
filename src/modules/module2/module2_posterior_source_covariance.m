function Sigma_jj_post = module2_posterior_source_covariance(Sigma_jj_omega, L, Sigma_xi_xi, options)
% MODULE2_POSTERIOR_SOURCE_COVARIANCE - Compute Posterior Source Covariance (SPC)
%
% Syntax:
%   Sigma_jj_post = module2_posterior_source_covariance(Sigma_jj_omega, L, Sigma_xi_xi)
%   Sigma_jj_post = module2_posterior_source_covariance(Sigma_jj_omega, L, Sigma_xi_xi, options)
%
% Description:
%   Computes the posterior source covariance for the E-step using the
%   information-form update:
%       Sigma_post = (Sigma_jj^{-1} + L^H * Sigma_xi^{-1} * L)^{-1}
%   This is numerically more stable than the subtractive form
%       Sigma_post = Sigma_jj - Sigma_jj L^H (L Sigma_jj L^H + Sigma_xi)^{-1} L Sigma_jj
%   as it avoids potential loss of positive semidefiniteness due to
%   subtraction/cancellation.
%
% Name-Value Arguments (options):
%   regularization_factor   - (double >=0) base jitter scaling, default 1e-8
%   jitter_max_tries        - (integer >=0) maximum jitter escalations, default 6
%   ensure_positive_definite- (logical) clip eigenvalues if needed, default true
%   min_eigenvalue_ratio    - (double in (0,1)) min eig floor as ratio of max eig, default 1e-12
%   verbose                 - (logical) print diagnostics, default false
%
% Output:
%   Sigma_jj_post - (complex, n×n) posterior source covariance (Hermitian PSD/PD)
%
% Notes:
%   - No explicit inv/cond; solve linear systems via jittered Cholesky.
%   - Inputs are numerically symmetrized to reduce floating asymmetry.
%   - If inputs are real and the result has tiny imaginary leakage, cast to real.

    %% -------- Input validation --------
    if ~isnumeric(Sigma_jj_omega) || ~ismatrix(Sigma_jj_omega)
        error('module2_posterior_source_covariance:invalid_source_covariance', ...
              'Sigma_jj_omega must be a numeric matrix');
    end
    [n1,n2] = size(Sigma_jj_omega);
    if n1 ~= n2
        error('module2_posterior_source_covariance:not_square_source', ...
              'Sigma_jj_omega must be square, got %d×%d', n1, n2);
    end
    n = n1;

    if ~isnumeric(L) || ndims(L) ~= 2
        error('module2_posterior_source_covariance:invalid_leadfield', ...
              'L must be a 2D numeric array');
    end
    [p, n_sources] = size(L);
    if n_sources ~= n
        error('module2_posterior_source_covariance:dimension_mismatch_leadfield', ...
              'Leadfield columns (%d) must match Sigma_jj_omega size (%d)', n_sources, n);
    end

    if ~isnumeric(Sigma_xi_xi) || ~ismatrix(Sigma_xi_xi)
        error('module2_posterior_source_covariance:invalid_noise_covariance', ...
              'Sigma_xi_xi must be a numeric matrix');
    end
    if size(Sigma_xi_xi,1) ~= p || size(Sigma_xi_xi,2) ~= p
        error('module2_posterior_source_covariance:dimension_mismatch_noise', ...
              'Sigma_xi_xi must be p×p where p = size(L,1)');
    end

    if any(isnan(Sigma_jj_omega(:))) || any(isinf(Sigma_jj_omega(:)))
        error('module2_posterior_source_covariance:invalid_source_values', ...
              'Sigma_jj_omega contains NaN or Inf');
    end
    if any(isnan(L(:))) || any(isinf(L(:)))
        error('module2_posterior_source_covariance:invalid_leadfield_values', ...
              'Leadfield L contains NaN or Inf');
    end
    if any(isnan(Sigma_xi_xi(:))) || any(isinf(Sigma_xi_xi(:)))
        error('module2_posterior_source_covariance:invalid_noise_values', ...
              'Sigma_xi_xi contains NaN or Inf');
    end

    %% -------- Parse options --------
    if nargin < 4 || isempty(options); options = struct(); end
    defaults = struct( ...
        'regularization_factor', 1e-8, ...
        'jitter_max_tries', 6, ...
        'ensure_positive_definite', true, ...
        'min_eigenvalue_ratio', 1e-12, ...
        'verbose', false);
    fns = fieldnames(defaults);
    for i = 1:numel(fns)
        if ~isfield(options, fns{i})
            options.(fns{i}) = defaults.(fns{i});
        end
    end
    if ~isscalar(options.regularization_factor) || options.regularization_factor < 0
        error('module2_posterior_source_covariance:invalid_regularization', ...
              'regularization_factor must be a non-negative scalar');
    end
    if ~isscalar(options.jitter_max_tries) || options.jitter_max_tries < 0 || fix(options.jitter_max_tries) ~= options.jitter_max_tries
        error('module2_posterior_source_covariance:invalid_jitter_max_tries', ...
              'jitter_max_tries must be a non-negative integer');
    end
    if ~isscalar(options.min_eigenvalue_ratio) || options.min_eigenvalue_ratio <= 0 || options.min_eigenvalue_ratio >= 1
        error('module2_posterior_source_covariance:invalid_eigenvalue_ratio', ...
              'min_eigenvalue_ratio must be in (0,1)');
    end
    if ~islogical(options.ensure_positive_definite) || ~isscalar(options.ensure_positive_definite)
        error('module2_posterior_source_covariance:invalid_epd_flag', ...
              'ensure_positive_definite must be a logical scalar');
    end
    if ~islogical(options.verbose) || ~isscalar(options.verbose)
        error('module2_posterior_source_covariance:invalid_verbose_flag', ...
              'verbose must be a logical scalar');
    end

    if options.verbose
        fprintf('module2_posterior_source_covariance: n=%d sources, p=%d sensors\n', n, p);
    end

    %% -------- Numerical symmetrization (stability hygiene) --------
    % Reduce floating-point non-Hermitian noise in covariance inputs.
    Sigma_jj_omega = (Sigma_jj_omega + Sigma_jj_omega')/2;
    Sigma_xi_xi    = (Sigma_xi_xi    + Sigma_xi_xi')/2;

    %% -------- Information-form assembly --------
    % Compute Sigma_xi_xi^{-1} * L via linear solves (no explicit inverse)
    B = Sigma_xi_xi \ L;                 % p×n
    % Build the information matrix:
    %   G = Sigma_jj^{-1} + L^H * Sigma_xi^{-1} * L  (n×n, Hermitian)
    I_n = eye(n, class(Sigma_jj_omega));
    G = (Sigma_jj_omega \ I_n) + (L' * B);
    G = (G + G')/2;                      % enforce symmetry numerically

    %% -------- Compute Sigma_post by solving G * X = I --------
    % Solve for the inverse via jittered Cholesky:
    [R, jitter, ok] = chol_with_jitter(G, options.regularization_factor, options.jitter_max_tries, options.verbose);
    if ~ok
        error('module2_posterior_source_covariance:spd_factorization_failed', ...
              'Failed to factorize G even with jitter (last jitter=%.3e)', jitter);
    end
    % Inverse from Cholesky: G = R'R  =>  G^{-1} = R \ (R' \ I)
    Y = R' \ I_n;
    Sigma_jj_post = R \ Y;

    % Ensure Hermitian numerically
    Sigma_jj_post = (Sigma_jj_post + Sigma_jj_post')/2;

    %% -------- Optional eigenvalue clipping (rarely needed with info-form) --------
    if options.ensure_positive_definite
        % In practice, Sigma_post should be PD if Sigma_jj and Sigma_xi are PD.
        % For PSD inputs or severe round-off, clip eigenvalues to a relative floor.
        [V, D] = eig(Sigma_jj_post);
        d = real(diag(D));
        dmax = max(d);
        if ~isempty(dmax) && isfinite(dmax)
            floor_val = options.min_eigenvalue_ratio * max(dmax, eps(class(dmax)));
            idx = d < floor_val;
            if any(idx)
                if options.verbose
                    fprintf('Eigenvalue clipping applied: min %.3e -> floor %.3e (count=%d)\n', ...
                            min(d), floor_val, nnz(idx));
                end
                d(idx) = floor_val;
                Sigma_jj_post = V * diag(d) * V';
                Sigma_jj_post = (Sigma_jj_post + Sigma_jj_post')/2;
            end
        end
    end

    %% -------- Output checks --------
    if any(isnan(Sigma_jj_post(:))) || any(isinf(Sigma_jj_post(:)))
        error('module2_posterior_source_covariance:invalid_output', ...
              'Output contains NaN or Inf');
    end
    if size(Sigma_jj_post,1) ~= n || size(Sigma_jj_post,2) ~= n
        error('module2_posterior_source_covariance:output_dimension_error', ...
              'Output must be n×n');
    end

    % If inputs are real and imaginary leakage is tiny, drop it
    if isreal(L) && isreal(Sigma_jj_omega) && isreal(Sigma_xi_xi)
        imax = max(abs(imag(Sigma_jj_post(:))));
        if imax < 1e-13
            Sigma_jj_post = real(Sigma_jj_post);
        end
    end

    % Optional diagnostic: uncertainty reduction (Sigma_prior - Sigma_post) >= 0
    if options.verbose
        try
            Dred = eig(Sigma_jj_omega - Sigma_jj_post);
            if any(real(Dred) < -1e-10)
                warning('module2_posterior_source_covariance:uncertainty_increase', ...
                        'Uncertainty reduction check failed (numerical tolerance exceeded).');
            end
            fprintf('||Sigma_post||_F = %.6e, trace ratio = %.6f\n', ...
                     norm(Sigma_jj_post,'fro'), trace(Sigma_jj_post)/max(1,trace(Sigma_jj_omega)));
        catch
            fprintf('Skipped uncertainty reduction eigen-check due to numerical issues.\n');
        end
    end
end

%% ===== Helper: Cholesky with jitter for (semi)SPD matrices =====
function [R, last_jitter, ok] = chol_with_jitter(A, base_reg, max_tries, verbose)
% Attempt a Cholesky factorization of A. If it fails, add jitter*I with
% geometrically increasing jitter until success or attempts exhausted.
    ok = false;
    A = (A + A')/2; % enforce symmetry
    scale = trace(A)/max(1,size(A,1));
    if ~isfinite(scale) || scale <= 0, scale = 1; end
    jitter0 = base_reg * scale;
    last_jitter = 0;

    for k = 0:max_tries
        jitter = jitter0 * (10^k);
        At = A + jitter * eye(size(A), class(A));
        [R, flag] = chol(At);
        if flag == 0
            ok = true;
            last_jitter = jitter;
            if verbose
                fprintf('  [chol] success with jitter = %.3e (tries=%d)\n', jitter, k+1);
            end
            return;
        else
            if verbose
                fprintf('  [chol] fail @ jitter = %.3e\n', jitter);
            end
        end
    end
    R = [];
    last_jitter = jitter0 * (10^max_tries);
end
