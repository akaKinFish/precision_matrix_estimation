function S_xi_xi = module2_residual_empirical_covariance(T_xi_v, S_vv, options)
% MODULE2_RESIDUAL_EMPIRICAL_COVARIANCE - Compute Residual Empirical Covariance (REC)
%
% Syntax:
%   S_xi_xi = module2_residual_empirical_covariance(T_xi_v, S_vv)
%   S_xi_xi = module2_residual_empirical_covariance(T_xi_v, S_vv, options)
%
% Description:
%   Computes the residual empirical covariance:
%       S_xi_xi = T_xi_v * S_vv * T_xi_v^H
%   where (') denotes conjugate transpose. This captures the empirical noise
%   covariance in the residual space after removing estimated source activity.
%
% Name-Value Arguments (options):
%   enforce_hermitian          - (logical) Symmetrize inputs/outputs numerically. Default: true
%   regularization_factor      - (double >= 0) Ridge on S_vv: S_vv + λ I, λ = factor * tr(S_vv)/p. Default: 0
%   min_eigenvalue_threshold   - (double >= 0) Eigenvalue floor for output (see psd_only). Default: 1e-12
%   psd_only                   - (logical) If true, project to PSD cone (floor negatives to 0, not PD).
%                                If false, floor all eigenvalues below threshold up to threshold (PD). Default: true
%   numerical_tolerance        - (double > 0) Tolerance for checks/projections. Default: 1e-12
%   verbose                    - (logical) Print diagnostics. Default: false
%
% Output:
%   S_xi_xi - (complex, p×p) residual empirical covariance (Hermitian PSD/PD)
%
% Notes:
%   - We always avoid explicit matrix inverses.
%   - If inputs are real and the result has tiny imaginary leakage, we cast to real.

    %% -------- Input validation --------
    if ~isnumeric(T_xi_v) || ~ismatrix(T_xi_v)
        error('module2_residual_empirical_covariance:invalid_transfer_function', ...
              'T_xi_v must be a numeric matrix');
    end
    [p1, p2] = size(T_xi_v);
    if p1 ~= p2
        error('module2_residual_empirical_covariance:not_square_transfer', ...
              'T_xi_v must be square, got %d×%d', p1, p2);
    end
    p = p1;

    if ~isnumeric(S_vv) || ~ismatrix(S_vv)
        error('module2_residual_empirical_covariance:invalid_empirical_covariance', ...
              'S_vv must be a numeric matrix');
    end
    if ~isequal(size(S_vv), [p p])
        error('module2_residual_empirical_covariance:dimension_mismatch', ...
              'S_vv must be %d×%d to match T_xi_v', p, p);
    end

    if any(isnan(T_xi_v(:))) || any(isinf(T_xi_v(:)))
        error('module2_residual_empirical_covariance:invalid_transfer_values', ...
              'T_xi_v contains NaN or Inf values');
    end
    if any(isnan(S_vv(:))) || any(isinf(S_vv(:)))
        error('module2_residual_empirical_covariance:invalid_covariance_values', ...
              'S_vv contains NaN or Inf values');
    end

    %% -------- Parse options --------
    if nargin < 3 || isempty(options), options = struct(); end
    defaults = struct( ...
        'enforce_hermitian', true, ...
        'regularization_factor', 0, ...
        'min_eigenvalue_threshold', 1e-12, ...
        'psd_only', true, ...
        'numerical_tolerance', 1e-12, ...
        'verbose', false);
    fns = fieldnames(defaults);
    for i = 1:numel(fns)
        if ~isfield(options, fns{i}), options.(fns{i}) = defaults.(fns{i}); end
    end

    if ~islogical(options.enforce_hermitian) || ~isscalar(options.enforce_hermitian)
        error('module2_residual_empirical_covariance:invalid_hermitian_flag', ...
              'enforce_hermitian must be a logical scalar');
    end
    if ~isscalar(options.regularization_factor) || options.regularization_factor < 0
        error('module2_residual_empirical_covariance:invalid_regularization', ...
              'regularization_factor must be a non-negative scalar');
    end
    if ~isscalar(options.min_eigenvalue_threshold) || options.min_eigenvalue_threshold < 0
        error('module2_residual_empirical_covariance:invalid_eigenvalue_threshold', ...
              'min_eigenvalue_threshold must be a non-negative scalar');
    end
    if ~isscalar(options.numerical_tolerance) || options.numerical_tolerance <= 0
        error('module2_residual_empirical_covariance:invalid_tolerance', ...
              'numerical_tolerance must be a positive scalar');
    end
    if ~islogical(options.psd_only) || ~isscalar(options.psd_only)
        error('module2_residual_empirical_covariance:invalid_psd_only', ...
              'psd_only must be a logical scalar');
    end
    if ~islogical(options.verbose) || ~isscalar(options.verbose)
        error('module2_residual_empirical_covariance:invalid_verbose_flag', ...
              'verbose must be a logical scalar');
    end

    %% -------- Input hygiene: symmetrize S_vv if needed --------
    if options.enforce_hermitian
        rel_asym_S = norm(S_vv - S_vv', 'fro') / max(1, norm(S_vv, 'fro'));
        if rel_asym_S > 1e-10 && options.verbose
            fprintf('S_vv not Hermitian (rel=%.2e); symmetrizing.\n', rel_asym_S);
        end
        S_vv = (S_vv + S_vv')/2;
    else
        % Keep original behavior: hard error if significantly non-Hermitian
        rel_asym_S = norm(S_vv - S_vv', 'fro') / max(1, norm(S_vv, 'fro'));
        if rel_asym_S > 1e-10
            error('module2_residual_empirical_covariance:not_hermitian_covariance', ...
                  'S_vv is not Hermitian (rel=%.2e). Consider enforce_hermitian=true.', rel_asym_S);
        end
    end

    %% -------- Optional ridge regularization on S_vv --------
    if options.regularization_factor > 0
        lam = options.regularization_factor * trace(S_vv) / max(1, p);
        if options.verbose
            fprintf('Applying ridge on S_vv: lambda = %.3e\n', lam);
        end
        S_vv = S_vv + lam * eye(p, class(S_vv));
    end

    %% -------- Core triple product --------
    % Compute S_xi_xi = T_xi_v * S_vv * T_xi_v^H
    % Order the multiplications to reduce temporary peak memory.
    tmp = S_vv * T_xi_v';     % (p×p)
    S_xi_xi = T_xi_v * tmp;   % (p×p)

    % Enforce Hermitian symmetry for numerical stability
    if options.enforce_hermitian
        S_xi_xi = (S_xi_xi + S_xi_xi')/2;
    end

    %% -------- PSD / PD projection (optional but recommended) --------
    % Project onto PSD cone, or floor to PD if psd_only=false.
    try
        [V, D] = eig(S_xi_xi);
        d = real(diag(D));
        % First, kill tiny negative numerical noise
        neg_tol = 10 * options.numerical_tolerance;
        d(d < 0 & d > -neg_tol) = 0;

        if options.psd_only
            % Pure PSD projection: set negatives to 0 (do not force PD)
            num_neg = sum(d < 0);
            if num_neg > 0 && options.verbose
                fprintf('Projecting to PSD: %d negative eigenvalues set to 0 (min=%.3e)\n', ...
                        num_neg, min(d));
            end
            d(d < 0) = 0;
        else
            % PD floor: raise everything below threshold up to the threshold
            if options.min_eigenvalue_threshold > 0
                floor_val = options.min_eigenvalue_threshold;
                num_raised = sum(d < floor_val);
                if num_raised > 0 && options.verbose
                    fprintf('Raising eigenvalues to floor %.3e (count=%d, min before=%.3e)\n', ...
                            floor_val, num_raised, min(d));
                end
                d(d < floor_val) = floor_val;
            end
        end

        S_xi_xi = V * diag(d) * V';
        if options.enforce_hermitian
            S_xi_xi = (S_xi_xi + S_xi_xi')/2;
        end
    catch
        if options.verbose
            fprintf('Eigen-projection skipped due to numerical issues.\n');
        end
    end

    %% -------- Output hygiene --------
    if any(isnan(S_xi_xi(:))) || any(isinf(S_xi_xi(:)))
        error('module2_residual_empirical_covariance:invalid_output', ...
              'Output S_xi_xi contains NaN or Inf values');
    end
    if ~isequal(size(S_xi_xi), [p p])
        error('module2_residual_empirical_covariance:output_dimension_error', ...
              'Output must be %d×%d, got %d×%d', p, p, size(S_xi_xi,1), size(S_xi_xi,2));
    end

    % Drop tiny imaginary leakage if inputs were real
    if isreal(T_xi_v) && isreal(S_vv)
        imax = max(abs(imag(S_xi_xi(:))));
        if imax < 1e-13
            S_xi_xi = real(S_xi_xi);
        end
    end

    %% -------- Diagnostics (verbose) --------
    if options.verbose
        froN = norm(S_xi_xi, 'fro');
        trS  = real(trace(S_xi_xi));
        try
            ev = real(eig(S_xi_xi));
            minev = min(ev); maxev = max(ev);
            rS   = rank(S_vv, options.numerical_tolerance);
            rX   = rank(S_xi_xi, options.numerical_tolerance);
            fprintf('REC diagnostics:\n');
            fprintf('  - ||S_xi_xi||_F = %.6e, trace = %.6e\n', froN, trS);
            fprintf('  - eig range     = [%.3e, %.3e]\n', minev, maxev);
            fprintf('  - rank(S_vv) = %d, rank(S_xi_xi) = %d (ratio=%.3f)\n', rS, rX, rX/max(1,rS));
            if minev < -options.numerical_tolerance
                warning('module2_residual_empirical_covariance:not_psd', ...
                        'Output not PSD (min eigenvalue = %.3e). Consider psd_only=true or larger threshold.', minev);
            end
        catch
            fprintf('  - Eigen diagnostics skipped due to numerical issues.\n');
        end
        fprintf('  - Is Hermitian: %s\n', mat2str(ishermitian(S_xi_xi)));
        fprintf('  - Is complex:   %s\n', mat2str(~isreal(S_xi_xi)));
    end
end
