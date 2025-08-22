function Sigma_jj_post = module2_posterior_source_covariance(Sigma_jj_omega, L, Sigma_xi_xi, options)
    % MODULE2_POSTERIOR_SOURCE_COVARIANCE - Compute Posterior Source Covariance (SPC)
    %
    % Syntax:
    %   Sigma_jj_post = module2_posterior_source_covariance(Sigma_jj_omega, L, Sigma_xi_xi)
    %   Sigma_jj_post = module2_posterior_source_covariance(Sigma_jj_omega, L, Sigma_xi_xi, options)
    %
    % Description:
    %   Computes the posterior source covariance matrix for the E-step:
    %   Σ_jj,post^(ω) = Σ_jj^(ω) - Σ_jj^(ω) L^H (L Σ_jj^(ω) L^H + Σ_ξξ)^(-1) L Σ_jj^(ω)
    %   
    %   This represents the uncertainty reduction in source estimates after
    %   incorporating sensor measurements, following the Kalman filter update equations.
    %
    % Input Arguments:
    %   Sigma_jj_omega - (complex, n×n) Source prior covariance matrix at frequency ω
    %   L - (double, p×n) Leadfield matrix mapping sources to sensors
    %   Sigma_xi_xi - (complex, p×p) Sensor noise covariance matrix
    %
    % Name-Value Arguments:
    %   regularization_factor - (double) Regularization parameter for numerical stability. Default: 1e-8
    %   condition_threshold - (double) Condition number threshold for stability warning. Default: 1e12
    %   ensure_positive_definite - (logical) Force output to be positive definite. Default: true
    %   min_eigenvalue_ratio - (double) Minimum eigenvalue as ratio of maximum. Default: 1e-12
    %   verbose - (logical) Enable verbose output for debugging. Default: false
    %
    % Output Arguments:
    %   Sigma_jj_post - (complex, n×n) Posterior source covariance matrix
    %
    % Examples:
    %   % Basic usage
    %   L = randn(64, 100);  % 64 sensors, 100 sources
    %   Sigma_jj = eye(100) * 0.5;  % Source prior covariance
    %   Sigma_xi = eye(64) * 0.1;   % Sensor noise covariance
    %   Sigma_post = module2_posterior_source_covariance(Sigma_jj, L, Sigma_xi);
    %   
    %   % Usage with numerical stability options
    %   Sigma_post = module2_posterior_source_covariance(Sigma_jj, L, Sigma_xi, ...
    %                                                     'ensure_positive_definite', true, ...
    %                                                     'min_eigenvalue_ratio', 1e-10);
    %
    % Mathematical Background:
    %   The posterior covariance represents the Bayesian update of source
    %   uncertainty after incorporating sensor measurements. It satisfies:
    %   - Σ_jj,post ≼ Σ_jj (uncertainty reduction)
    %   - Σ_jj,post → Σ_jj as sensor noise increases
    %   - Σ_jj,post → 0 as sensor noise decreases (perfect observations)
    %
    % Numerical Considerations:
    %   The computation involves matrix subtraction which can lead to loss of
    %   positive definiteness due to numerical errors. The function includes
    %   options to enforce matrix properties through eigenvalue clipping.
    %
    % See also: MODULE2_DSTF_COMPUTATION, MODULE2_RESIDUAL_TRANSFER_FUNCTION
    %
    % Author: [Author Name]
    % Date: [Current Date]
    % Version: 1.0
    
    %% Input validation
    if ~isnumeric(Sigma_jj_omega) || ~ismatrix(Sigma_jj_omega)
        error('module2_posterior_source_covariance:invalid_source_covariance', ...
              'Source covariance Sigma_jj_omega must be a numeric matrix');
    end
    
    [n1, n2] = size(Sigma_jj_omega);
    if n1 ~= n2
        error('module2_posterior_source_covariance:not_square_source', ...
              'Source covariance must be square, got %d×%d', n1, n2);
    end
    n = n1;
    
    if ~isnumeric(L) || ndims(L) ~= 2
        error('module2_posterior_source_covariance:invalid_leadfield', ...
              'Leadfield matrix L must be a 2D numeric array');
    end
    
    [p, n_sources] = size(L);
    if n_sources ~= n
        error('module2_posterior_source_covariance:dimension_mismatch_leadfield', ...
              'Leadfield sources (%d) must match source covariance dimension (%d)', ...
              n_sources, n);
    end
    
    if ~isnumeric(Sigma_xi_xi) || ~ismatrix(Sigma_xi_xi)
        error('module2_posterior_source_covariance:invalid_noise_covariance', ...
              'Noise covariance Sigma_xi_xi must be a numeric matrix');
    end
    
    if size(Sigma_xi_xi, 1) ~= p || size(Sigma_xi_xi, 2) ~= p
        error('module2_posterior_source_covariance:dimension_mismatch_noise', ...
              'Noise covariance dimensions (%d×%d) must match leadfield sensors (%d)', ...
              size(Sigma_xi_xi, 1), size(Sigma_xi_xi, 2), p);
    end
    
    % Check for valid covariance matrices (Hermitian)
    if ~ishermitian(Sigma_jj_omega)
        hermitian_error = norm(Sigma_jj_omega - Sigma_jj_omega', 'fro') / norm(Sigma_jj_omega, 'fro');
        if hermitian_error > 1e-10
            error('module2_posterior_source_covariance:not_hermitian_source', ...
                  'Source covariance is not Hermitian (relative error: %.2e)', hermitian_error);
        end
    end
    
    if ~ishermitian(Sigma_xi_xi)
        hermitian_error = norm(Sigma_xi_xi - Sigma_xi_xi', 'fro') / norm(Sigma_xi_xi, 'fro');
        if hermitian_error > 1e-10
            error('module2_posterior_source_covariance:not_hermitian_noise', ...
                  'Noise covariance is not Hermitian (relative error: %.2e)', hermitian_error);
        end
    end
    
    % Check for NaN or Inf values
    if any(isnan(Sigma_jj_omega(:))) || any(isinf(Sigma_jj_omega(:)))
        error('module2_posterior_source_covariance:invalid_source_values', ...
              'Source covariance contains NaN or Inf values');
    end
    
    if any(isnan(L(:))) || any(isinf(L(:)))
        error('module2_posterior_source_covariance:invalid_leadfield_values', ...
              'Leadfield matrix contains NaN or Inf values');
    end
    
    if any(isnan(Sigma_xi_xi(:))) || any(isinf(Sigma_xi_xi(:)))
        error('module2_posterior_source_covariance:invalid_noise_values', ...
              'Noise covariance contains NaN or Inf values');
    end
    
    %% Parse optional arguments
    if nargin < 4
        options = struct();
    end
    
    % Set default options
    default_options = struct(...
        'regularization_factor', 1e-8, ...
        'condition_threshold', 1e12, ...
        'ensure_positive_definite', true, ...
        'min_eigenvalue_ratio', 1e-12, ...
        'verbose', false);
    
    % Merge user options with defaults
    option_names = fieldnames(default_options);
    for i = 1:length(option_names)
        field_name = option_names{i};
        if ~isfield(options, field_name)
            options.(field_name) = default_options.(field_name);
        end
    end
    
    % Validate option types
    if ~isscalar(options.regularization_factor) || ~isnumeric(options.regularization_factor) || options.regularization_factor < 0
        error('module2_posterior_source_covariance:invalid_regularization', ...
              'regularization_factor must be a non-negative scalar');
    end
    
    if ~isscalar(options.min_eigenvalue_ratio) || ~isnumeric(options.min_eigenvalue_ratio) || ...
            options.min_eigenvalue_ratio <= 0 || options.min_eigenvalue_ratio >= 1
        error('module2_posterior_source_covariance:invalid_eigenvalue_ratio', ...
              'min_eigenvalue_ratio must be a scalar in (0, 1)');
    end
    
    %% Main computation
    try
        if options.verbose
            fprintf('Computing posterior source covariance with dimensions: sources=%d, sensors=%d\n', n, p);
        end
        
        % Step 1: Compute intermediate matrices
        L_Sigma_jj = L * Sigma_jj_omega;
        A = L_Sigma_jj * L' + Sigma_xi_xi;  % Sensor space covariance
        
        % Step 2: Ensure A is Hermitian and well-conditioned
        if ~ishermitian(A)
            if options.verbose
                hermitian_error = norm(A - A', 'fro') / norm(A, 'fro');
                fprintf('Matrix A is not Hermitian (relative error: %.2e), symmetrizing...\n', hermitian_error);
            end
            A = (A + A') / 2;
        end
        
        % Step 3: Check condition number and apply regularization if needed
        cond_A = cond(A);
        if options.verbose
            fprintf('Condition number of sensor covariance: %.2e\n', cond_A);
        end
        
        if cond_A > options.condition_threshold
            if options.verbose
                warning('module2_posterior_source_covariance:high_condition_number', ...
                        'High condition number (%.2e) detected, applying regularization', cond_A);
            end
            
            reg_strength = options.regularization_factor * trace(A) / p;
            A = A + reg_strength * eye(p);
            
            if options.verbose
                fprintf('Applied regularization: strength = %.2e\n', reg_strength);
                fprintf('New condition number: %.2e\n', cond(A));
            end
        end
        
        % Step 4: Compute the correction term: Σ_jj L^H A^(-1) L Σ_jj
        A_inv = inv(A);
        correction_term = Sigma_jj_omega * L' * A_inv * L * Sigma_jj_omega;
        
        % Step 5: Compute posterior covariance
        Sigma_jj_post = Sigma_jj_omega - correction_term;
        
        % Step 6: Ensure Hermitian symmetry (critical for numerical stability)
        Sigma_jj_post = (Sigma_jj_post + Sigma_jj_post') / 2;
        
        % Step 7: Enforce positive definiteness if requested
        if options.ensure_positive_definite
            [V, D] = eig(Sigma_jj_post);
            eigenvals = real(diag(D));
            
            % Check if any eigenvalues are negative or too small
            min_allowed_eigenval = options.min_eigenvalue_ratio * max(eigenvals);
            negative_eigenvals = sum(eigenvals < 0);
            small_eigenvals = sum(eigenvals < min_allowed_eigenval & eigenvals >= 0);
            
            if negative_eigenvals > 0 || small_eigenvals > 0
                if options.verbose
                    fprintf('Eigenvalue correction: %d negative, %d too small\n', ...
                            negative_eigenvals, small_eigenvals);
                end
                
                % Clip eigenvalues to ensure positive definiteness
                eigenvals(eigenvals < min_allowed_eigenval) = min_allowed_eigenval;
                
                % Reconstruct matrix
                Sigma_jj_post = V * diag(eigenvals) * V';
                
                % Ensure Hermitian after reconstruction
                Sigma_jj_post = (Sigma_jj_post + Sigma_jj_post') / 2;
                
                if options.verbose
                    fprintf('Applied eigenvalue clipping, min eigenvalue: %.2e\n', min(eigenvals));
                end
            end
        end
        
        if options.verbose
            fprintf('Posterior source covariance computation completed successfully\n');
        end
        
    catch ME
        % Enhanced error handling with context
        switch ME.identifier
            case 'MATLAB:matrix:singular'
                error('module2_posterior_source_covariance:singular_matrix', ...
                      'Sensor covariance matrix is singular. Consider increasing regularization_factor');
            case 'MATLAB:matrix:posdef'
                error('module2_posterior_source_covariance:not_positive_definite', ...
                      'Input covariance matrices are not positive definite');
            otherwise
                rethrow(ME);
        end
    end
    
    %% Output validation
    % Check for NaN or Inf in output
    if any(isnan(Sigma_jj_post(:))) || any(isinf(Sigma_jj_post(:)))
        error('module2_posterior_source_covariance:invalid_output', ...
              'Output contains NaN or Inf values, indicating numerical instability');
    end
    
    % Verify output dimensions
    if size(Sigma_jj_post, 1) ~= n || size(Sigma_jj_post, 2) ~= n
        error('module2_posterior_source_covariance:output_dimension_error', ...
              'Output dimensions (%d×%d) do not match expected (%d×%d)', ...
              size(Sigma_jj_post, 1), size(Sigma_jj_post, 2), n, n);
    end
    
    % Verify uncertainty reduction property: Σ_post ≼ Σ_prior
    try
        eigenvals_diff = eig(Sigma_jj_omega - Sigma_jj_post);
        if any(real(eigenvals_diff) < -1e-10)  % Allow small numerical errors
            warning('module2_posterior_source_covariance:uncertainty_increase', ...
                    'Posterior covariance may not satisfy uncertainty reduction property');
        end
    catch
        % Skip this check if eigenvalue computation fails
        if options.verbose
            fprintf('Skipped uncertainty reduction check due to numerical issues\n');
        end
    end
    
    if options.verbose
        % Report on matrix properties
        frobenius_norm = norm(Sigma_jj_post, 'fro');
        trace_ratio = trace(Sigma_jj_post) / trace(Sigma_jj_omega);
        
        try
            min_eigenval = min(real(eig(Sigma_jj_post)));
            condition_number = cond(Sigma_jj_post);
        catch
            min_eigenval = NaN;
            condition_number = NaN;
        end
        
        fprintf('Output validation:\n');
        fprintf('  - Frobenius norm: %.6f\n', frobenius_norm);
        fprintf('  - Trace reduction ratio: %.6f\n', trace_ratio);
        fprintf('  - Minimum eigenvalue: %.2e\n', min_eigenval);
        fprintf('  - Condition number: %.2e\n', condition_number);
        fprintf('  - Is Hermitian: %s\n', mat2str(ishermitian(Sigma_jj_post)));
    end
end