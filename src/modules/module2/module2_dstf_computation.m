function T_jv = module2_dstf_computation(L, Sigma_jj_omega, Sigma_xi_xi, options)
% MODULE2_DSTF_COMPUTATION - Compute Data-to-Source Transfer Function (DSTF)
%
% Syntax:
%   T_jv = module2_dstf_computation(L, Sigma_jj_omega, Sigma_xi_xi)
%   T_jv = module2_dstf_computation(L, Sigma_jj_omega, Sigma_xi_xi, options)
%
% Description:
%   Computes the DSTF for the E-step:
%       T_jv^(ω) = Σ_jj^(ω) L^H (L Σ_jj^(ω) L^H + Σ_ξξ)^(-1)
%   This implementation prefers the information-form:
%       T = (Σ_jj^{-1} + L^H Σ_ξξ^{-1} L)^{-1} L^H Σ_ξξ^{-1}
%   It avoids explicit matrix inverses by solving linear systems with
%   jittered Cholesky factorization for robustness.
%
% Name-Value Arguments (options):
%   regularization_factor  - (double, >=0) base jitter scaling, default 1e-8
%   jitter_max_tries       - (integer, >=0) max jitter escalations, default 6
%   method                 - 'auto'|'information'|'sensor' (default 'auto')
%   use_pseudoinverse      - (logical) fallback to pinv only if factorization fails
%   verbose                - (logical) print diagnostic info
%
% Output:
%   T_jv - (complex, n×p) Data-to-Source Transfer Function
%
% Notes:
%   - No explicit inv(A) is used. We rely on mldivide and Cholesky.
%   - Inputs are numerically symmetrized to reduce floating asymmetry.
%   - If inputs are real and the result has tiny imaginary parts, cast to real.

    %% -------- Input validation --------
    if ~isnumeric(L) || ndims(L) ~= 2
        error('module2_dstf_computation:invalid_leadfield', ...
              'Leadfield matrix L must be a 2D numeric array');
    end
    [p, n] = size(L);

    if ~isnumeric(Sigma_jj_omega) || ~ismatrix(Sigma_jj_omega)
        error('module2_dstf_computation:invalid_source_covariance', ...
              'Source covariance Sigma_jj_omega must be a numeric matrix');
    end
    if size(Sigma_jj_omega,1) ~= n || size(Sigma_jj_omega,2) ~= n
        error('module2_dstf_computation:dimension_mismatch_source', ...
              'Source covariance must be n×n where n = size(L,2)');
    end

    if ~isnumeric(Sigma_xi_xi) || ~ismatrix(Sigma_xi_xi)
        error('module2_dstf_computation:invalid_noise_covariance', ...
              'Noise covariance Sigma_xi_xi must be a numeric matrix');
    end
    if size(Sigma_xi_xi,1) ~= p || size(Sigma_xi_xi,2) ~= p
        error('module2_dstf_computation:dimension_mismatch_noise', ...
              'Noise covariance must be p×p where p = size(L,1)');
    end

    if any(isnan(L(:))) || any(isinf(L(:)))
        error('module2_dstf_computation:invalid_leadfield_values', ...
              'Leadfield contains NaN or Inf');
    end
    if any(isnan(Sigma_jj_omega(:))) || any(isinf(Sigma_jj_omega(:)))
        error('module2_dstf_computation:invalid_source_values', ...
              'Source covariance contains NaN or Inf');
    end
    if any(isnan(Sigma_xi_xi(:))) || any(isinf(Sigma_xi_xi(:)))
        error('module2_dstf_computation:invalid_noise_values', ...
              'Noise covariance contains NaN or Inf');
    end

    %% -------- Parse options --------
    if nargin < 4 || isempty(options); options = struct(); end
    defaults = struct( ...
        'regularization_factor', 1e-8, ...
        'jitter_max_tries', 6, ...
        'method', 'auto', ...
        'use_pseudoinverse', false, ...
        'verbose', false);
    fns = fieldnames(defaults);
    for i = 1:numel(fns)
        if ~isfield(options, fns{i})
            options.(fns{i}) = defaults.(fns{i});
        end
    end

    if ~isscalar(options.regularization_factor) || options.regularization_factor < 0
        error('module2_dstf_computation:invalid_regularization', ...
              'regularization_factor must be a non-negative scalar');
    end
    if ~isscalar(options.jitter_max_tries) || options.jitter_max_tries < 0 || fix(options.jitter_max_tries) ~= options.jitter_max_tries
        error('module2_dstf_computation:invalid_jitter_max_tries', ...
              'jitter_max_tries must be a non-negative integer');
    end
    if ~ischar(options.method) || ~ismember(lower(options.method), {'auto','information','sensor'})
        error('module2_dstf_computation:invalid_method', ...
              'method must be one of {''auto'',''information'',''sensor''}');
    end
    if ~islogical(options.use_pseudoinverse) || ~isscalar(options.use_pseudoinverse)
        error('module2_dstf_computation:invalid_pseudoinverse_flag', ...
              'use_pseudoinverse must be a logical scalar');
    end
    if ~islogical(options.verbose) || ~isscalar(options.verbose)
        error('module2_dstf_computation:invalid_verbose_flag', ...
              'verbose must be a logical scalar');
    end

    if options.verbose
        fprintf('module2_dstf_computation: p=%d sensors, n=%d sources\n', p, n);
    end

    %% -------- Numerical symmetrization (stability hygiene) --------
    % Reduce floating-point non-Hermitian noise in covariance inputs.
    Sigma_jj_omega = (Sigma_jj_omega + Sigma_jj_omega')/2;
    Sigma_xi_xi    = (Sigma_xi_xi    + Sigma_xi_xi')/2;

    %% -------- Choose computation form --------
    method = lower(options.method);
    if strcmp(method,'auto')
        % Information-form is typically more stable and cheaper if p >> n.
        method = 'information';
    end

    %% -------- Core computation (no explicit inverse) --------
    if strcmp(method,'information')
        % Information form:
        %   G = Σ_jj^{-1} + L^H Σ_ξξ^{-1} L
        %   RHS = L^H Σ_ξξ^{-1}
        %   Solve G * X = RHS, then T = X
        I_n = eye(n, class(Sigma_jj_omega));

        % Compute Sigma_xi_xi^{-1} * L by solving linear systems
        B = Sigma_xi_xi \ L;      % p×n

        % Build G without explicit inversion of Σ_jj
        G = (Sigma_jj_omega \ I_n) + (L' * B);  % n×n, theoretically Hermitian
        G = (G + G')/2;                         % enforce symmetry numerically

        RHS = B';                               % n×p, since L^H Σ_ξξ^{-1} = (Σ_ξξ^{-1} L)^H

        [X, ok, info] = solve_spd_with_jitter(G, RHS, options.regularization_factor, options.jitter_max_tries, options.verbose);
        if ~ok
            if options.use_pseudoinverse
                warning('module2_dstf_computation:using_pinv_fallback', ...
                        'G factorization failed; falling back to pinv (not recommended).');
                X = pinv(G) * RHS;
            else
                error('module2_dstf_computation:spd_solve_failed', ...
                      'Failed to factorize G even with jitter (last jitter=%.3e).', info.last_jitter);
            end
        end

        T_jv = X;  % n×p

    else
        % Sensor form:
        %   A = L Σ_jj L^H + Σ_ξξ
        %   Solve A * Y = L Σ_jj, then T = Y^H
        L_Sigma = L * Sigma_jj_omega;       % p×n
        A = L_Sigma * L' + Sigma_xi_xi;     % p×p
        A = (A + A')/2;

        [Y, ok, info] = solve_spd_with_jitter(A, L_Sigma, options.regularization_factor, options.jitter_max_tries, options.verbose);
        if ~ok
            if options.use_pseudoinverse
                warning('module2_dstf_computation:using_pinv_fallback', ...
                        'A factorization failed; falling back to pinv (not recommended).');
                Y = pinv(A) * L_Sigma;
            else
                error('module2_dstf_computation:spd_solve_failed', ...
                      'Failed to factorize A even with jitter (last jitter=%.3e).', info.last_jitter);
            end
        end

        % A^{-1} (L Σ) is Y; T = (Σ L^H A^{-1}) = (A^{-1} L Σ)^H = Y'
        T_jv = Y';
    end

    %% -------- Output checks and cleanup --------
    if any(isnan(T_jv(:))) || any(isinf(T_jv(:)))
        error('module2_dstf_computation:invalid_output', ...
              'Output T_jv contains NaN or Inf');
    end
    if size(T_jv,1) ~= n || size(T_jv,2) ~= p
        error('module2_dstf_computation:output_dimension_error', ...
              'Output must be n×p');
    end

    % If inputs are real and imaginary leakage is tiny, drop it.
    if isreal(L) && isreal(Sigma_jj_omega) && isreal(Sigma_xi_xi)
        imax = max(abs(imag(T_jv(:))));
        if imax < 1e-13
            T_jv = real(T_jv);
        end
    end

    if options.verbose
        fprintf('DSTF computation completed using method: %s\n', method);
        fprintf('||T_jv||_F = %.6e\n', norm(T_jv,'fro'));
    end
end

%% ===== Helper: SPD linear solve with jittered Cholesky =====
function [X, ok, info] = solve_spd_with_jitter(A, B, base_reg, max_tries, verbose)
% Solve A * X = B where A is expected Hermitian (semi)PD.
% Adds jitter*I progressively until Cholesky succeeds.
    ok = false; X = [];
    info = struct('last_jitter', 0, 'tries', 0);

    if ~ishermitian(A)
        if verbose
            he = norm(A - A', 'fro')/max(1,norm(A,'fro'));
            fprintf('  [sym] A not Hermitian (rel diff=%.2e). Symmetrizing...\n', he);
        end
        A = (A + A')/2;
    end

    scale = trace(A)/max(1,size(A,1));
    if ~isfinite(scale) || scale <= 0, scale = 1; end
    jitter0 = base_reg * scale;

    for k = 0:max_tries
        jitter = jitter0 * (10^k);
        A_try = A + jitter * eye(size(A), class(A));
        [R, flag] = chol(A_try);
        if flag == 0
            Y = R' \ B;
            X = R  \ Y;
            ok = true;
            info.last_jitter = jitter;
            info.tries = k+1;
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
    info.last_jitter = jitter0 * (10^max_tries);
    info.tries = max_tries + 1;
end
