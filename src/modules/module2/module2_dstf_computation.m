function T_jv = module2_dstf_computation(L, Sigma_jj_omega, Sigma_xi_xi, options)
    % MODULE2_DSTF_COMPUTATION - Compute Data-to-Source Transfer Function (DSTF)
    %
    % Syntax:
    %   T_jv = module2_dstf_computation(L, Sigma_jj_omega, Sigma_xi_xi)
    %   T_jv = module2_dstf_computation(L, Sigma_jj_omega, Sigma_xi_xi, options)
    %
    % Description:
    %   Computes the Data-to-Source Transfer Function (DSTF) for the E-step:
    %   T_jv^(ω) = Σ_jj^(ω) L^H (L Σ_jj^(ω) L^H + Σ_ξξ)^(-1)
    %   
    %   This function implements the core computation for mapping from sensor
    %   space to source space, accounting for both source prior covariance
    %   and sensor noise characteristics.
    %
    % Input Arguments:
    %   L - (double, p×n) Leadfield matrix mapping sources to sensors
    %   Sigma_jj_omega - (complex, n×n) Source prior covariance matrix at frequency ω
    %   Sigma_xi_xi - (complex, p×p) Sensor noise covariance matrix (typically diagonal)
    %
    % Name-Value Arguments:
    %   regularization_factor - (double) Regularization parameter for numerical stability. Default: 1e-8
    %   condition_threshold - (double) Condition number threshold for stability warning. Default: 1e12
    %   use_pseudoinverse - (logical) Use pseudoinverse for near-singular matrices. Default: false
    %   verbose - (logical) Enable verbose output for debugging. Default: false
    %
    % Output Arguments:
    %   T_jv - (complex, n×p) Data-to-Source Transfer Function matrix
    %
    % Examples:
    %   % Basic usage with simple matrices
    %   L = randn(64, 100);  % 64 sensors, 100 sources
    %   Sigma_jj = eye(100) * 0.5;  % Source prior covariance
    %   Sigma_xi = eye(64) * 0.1;   % Sensor noise covariance
    %   T_jv = module2_dstf_computation(L, Sigma_jj, Sigma_xi);
    %   
    %   % Advanced usage with regularization
    %   T_jv = module2_dstf_computation(L, Sigma_jj, Sigma_xi, ...
    %                                   'regularization_factor', 1e-6, ...
    %                                   'verbose', true);
    %
    % Mathematical Background:
    %   The DSTF represents the optimal linear estimator for source activity
    %   given sensor measurements, derived from Bayesian inference under
    %   Gaussian assumptions. The matrix inversion can be numerically
    %   challenging when the matrix L*Sigma_jj*L^H + Sigma_xi_xi is
    %   ill-conditioned, requiring careful regularization strategies.
    %
    % Edge Cases:
    %   - Large source variance (Σ_jj >> Σ_ξξ): T_jv approaches pseudoinverse behavior
    %   - High noise (Σ_ξξ >> L*Σ_jj*L^H): T_jv approaches zero
    %   - Identity leadfield (L=I), diagonal covariances: T_jv = α/(α+σ²)*I
    %
    % See also: MODULE2_POSTERIOR_SOURCE_COVARIANCE, MODULE2_RESIDUAL_TRANSFER_FUNCTION
    %
    % Author: [Author Name]
    % Date: [Current Date]
    % Version: 1.0
    
    %% Input validation
    if ~isnumeric(L) || ndims(L) ~= 2
        error('module2_dstf_computation:invalid_leadfield', ...
              'Leadfield matrix L must be a 2D numeric array');
    end
    
    [p, n] = size(L);
    
    if ~isnumeric(Sigma_jj_omega) || ~ismatrix(Sigma_jj_omega)
        error('module2_dstf_computation:invalid_source_covariance', ...
              'Source covariance Sigma_jj_omega must be a numeric matrix');
    end
    
    if size(Sigma_jj_omega, 1) ~= n || size(Sigma_jj_omega, 2) ~= n
        error('module2_dstf_computation:dimension_mismatch_source', ...
              'Source covariance dimensions (%d×%d) must match leadfield sources (%d)', ...
              size(Sigma_jj_omega, 1), size(Sigma_jj_omega, 2), n);
    end
    
    if ~isnumeric(Sigma_xi_xi) || ~ismatrix(Sigma_xi_xi)
        error('module2_dstf_computation:invalid_noise_covariance', ...
              'Noise covariance Sigma_xi_xi must be a numeric matrix');
    end
    
    if size(Sigma_xi_xi, 1) ~= p || size(Sigma_xi_xi, 2) ~= p
        error('module2_dstf_computation:dimension_mismatch_noise', ...
              'Noise covariance dimensions (%d×%d) must match leadfield sensors (%d)', ...
              size(Sigma_xi_xi, 1), size(Sigma_xi_xi, 2), p);
    end
    
    % Check for NaN or Inf values
    if any(isnan(L(:))) || any(isinf(L(:)))
        error('module2_dstf_computation:invalid_leadfield_values', ...
              'Leadfield matrix contains NaN or Inf values');
    end
    
    if any(isnan(Sigma_jj_omega(:))) || any(isinf(Sigma_jj_omega(:)))
        error('module2_dstf_computation:invalid_source_values', ...
              'Source covariance contains NaN or Inf values');
    end
    
    if any(isnan(Sigma_xi_xi(:))) || any(isinf(Sigma_xi_xi(:)))
        error('module2_dstf_computation:invalid_noise_values', ...
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
        'use_pseudoinverse', false, ...
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
        error('module2_dstf_computation:invalid_regularization', ...
              'regularization_factor must be a non-negative scalar');
    end
    
    if ~isscalar(options.condition_threshold) || ~isnumeric(options.condition_threshold) || options.condition_threshold <= 1
        error('module2_dstf_computation:invalid_condition_threshold', ...
              'condition_threshold must be a scalar greater than 1');
    end
    
    if ~islogical(options.use_pseudoinverse) || ~isscalar(options.use_pseudoinverse)
        error('module2_dstf_computation:invalid_pseudoinverse_flag', ...
              'use_pseudoinverse must be a logical scalar');
    end
    
    if ~islogical(options.verbose) || ~isscalar(options.verbose)
        error('module2_dstf_computation:invalid_verbose_flag', ...
              'verbose must be a logical scalar');
    end
    
    %% Main computation
    try
        if options.verbose
            fprintf('Computing DSTF with dimensions: sensors=%d, sources=%d\n', p, n);
        end
        
        % Step 1: Compute intermediate matrix A = L * Sigma_jj_omega * L^H + Sigma_xi_xi
        L_Sigma_jj = L * Sigma_jj_omega;
        A = L_Sigma_jj * L' + Sigma_xi_xi;
        
        % Step 2: Check numerical properties
        if ~ishermitian(A)
            if options.verbose
                hermitian_error = norm(A - A', 'fro') / norm(A, 'fro');
                fprintf('Matrix A is not Hermitian (relative error: %.2e), symmetrizing...\n', hermitian_error);
            end
            A = (A + A') / 2;  % Force Hermitian symmetry
        end
        
        % Step 3: Condition number analysis and regularization
        cond_A = cond(A);
        if options.verbose
            fprintf('Condition number of A: %.2e\n', cond_A);
        end
        
        if cond_A > options.condition_threshold
            if options.verbose
                warning('module2_dstf_computation:high_condition_number', ...
                        'High condition number (%.2e) detected, applying regularization', cond_A);
            end
            
            % Apply regularization: A_reg = A + ε * trace(A)/p * I
            reg_strength = options.regularization_factor * trace(A) / p;
            A = A + reg_strength * eye(p);
            
            if options.verbose
                fprintf('Applied regularization: strength = %.2e\n', reg_strength);
                fprintf('New condition number: %.2e\n', cond(A));
            end
        end
        
        % Step 4: Matrix inversion with stability checks
        if options.use_pseudoinverse
            % Use pseudoinverse for robust handling of near-singular matrices
            A_inv = pinv(A);
            if options.verbose
                fprintf('Used pseudoinverse for matrix inversion\n');
            end
        else
            % Standard matrix inversion
            A_inv = inv(A);
        end
        
        % Step 5: Compute final DSTF
        T_jv = Sigma_jj_omega * L' * A_inv;
        
        % Ensure numerical consistency for complex matrices
        if ~isreal(T_jv) && isreal(L) && isreal(Sigma_jj_omega)
            % If inputs are real but output is complex due to numerical errors
            if max(abs(imag(T_jv(:)))) < 1e-14
                T_jv = real(T_jv);
            end
        end
        
        if options.verbose
            fprintf('DSTF computation completed successfully\n');
            fprintf('Output matrix dimensions: %d×%d\n', size(T_jv, 1), size(T_jv, 2));
        end
        
    catch ME
        % Enhanced error handling with context
        switch ME.identifier
            case 'MATLAB:matrix:singular'
                error('module2_dstf_computation:singular_matrix', ...
                      'Matrix A = L*Sigma_jj*L^H + Sigma_xi_xi is singular. Consider increasing regularization_factor or using pseudoinverse');
            case 'MATLAB:matrix:posdef'
                error('module2_dstf_computation:not_positive_definite', ...
                      'Matrix A is not positive definite. Check input covariance matrices for validity');
            otherwise
                rethrow(ME);
        end
    end
    
    %% Output validation
    % Check for NaN or Inf in output
    if any(isnan(T_jv(:))) || any(isinf(T_jv(:)))
        error('module2_dstf_computation:invalid_output', ...
              'Output T_jv contains NaN or Inf values, indicating numerical instability');
    end
    
    % Verify output dimensions
    if size(T_jv, 1) ~= n || size(T_jv, 2) ~= p
        error('module2_dstf_computation:output_dimension_error', ...
              'Output dimensions (%d×%d) do not match expected (%d×%d)', ...
              size(T_jv, 1), size(T_jv, 2), n, p);
    end
    
    if options.verbose
        % Report on matrix properties
        max_element = max(abs(T_jv(:)));
        frobenius_norm = norm(T_jv, 'fro');
        
        fprintf('Output validation:\n');
        fprintf('  - Maximum absolute element: %.6f\n', max_element);
        fprintf('  - Frobenius norm: %.6f\n', frobenius_norm);
        fprintf('  - Output is complex: %s\n', mat2str(~isreal(T_jv)));
    end
end