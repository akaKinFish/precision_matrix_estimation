function T_xi_v = module2_residual_transfer_function(T_jv, L, options)
    % MODULE2_RESIDUAL_TRANSFER_FUNCTION - Compute Residual Transfer Function (RTF)
    %
    % Syntax:
    %   T_xi_v = module2_residual_transfer_function(T_jv, L)
    %   T_xi_v = module2_residual_transfer_function(T_jv, L, options)
    %
    % Description:
    %   Computes the Residual Transfer Function (RTF) for the E-step:
    %   T_ξv^(ω) = I_p - L * T_jv^(ω)
    %   
    %   This function maps sensor measurements to the residual space after
    %   removing the estimated source contributions. The RTF is crucial for
    %   computing the residual empirical covariance in the EM algorithm.
    %
    % Input Arguments:
    %   T_jv - (complex, n×p) Data-to-Source Transfer Function matrix
    %   L - (double, p×n) Leadfield matrix mapping sources to sensors
    %
    % Name-Value Arguments:
    %   enforce_identity_diagonal - (logical) Ensure diagonal elements are exactly 1-L*T_jv_diagonal. Default: true
    %   numerical_tolerance - (double) Tolerance for numerical checks. Default: 1e-12
    %   verbose - (logical) Enable verbose output for debugging. Default: false
    %
    % Output Arguments:
    %   T_xi_v - (complex, p×p) Residual Transfer Function matrix
    %
    % Examples:
    %   % Basic usage after computing DSTF
    %   L = randn(64, 100);  % 64 sensors, 100 sources
    %   Sigma_jj = eye(100) * 0.5;  % Source prior covariance
    %   Sigma_xi = eye(64) * 0.1;   % Sensor noise covariance
    %   T_jv = module2_dstf_computation(L, Sigma_jj, Sigma_xi);
    %   T_xi_v = module2_residual_transfer_function(T_jv, L);
    %   
    %   % Usage with numerical stability options
    %   T_xi_v = module2_residual_transfer_function(T_jv, L, ...
    %                                                'enforce_identity_diagonal', true, ...
    %                                                'verbose', true);
    %
    % Mathematical Background:
    %   The RTF represents the linear transformation that maps sensor
    %   measurements to the residual space after optimal source estimation.
    %   Key properties:
    %   - T_ξv is idempotent: T_ξv * T_ξv = T_ξv (projection operator)
    %   - Range(T_ξv) = orthogonal complement of Range(L)
    %   - T_ξv + L * T_jv = I_p (complementary projections)
    %
    % Edge Cases:
    %   - Perfect source estimation (no noise): T_ξv → 0
    %   - No source estimation (high noise): T_ξv → I
    %   - Orthogonal leadfield columns: T_ξv has block structure
    %
    % See also: MODULE2_DSTF_COMPUTATION, MODULE2_RESIDUAL_EMPIRICAL_COVARIANCE
    %
    % Author: [Author Name]
    % Date: [Current Date]
    % Version: 1.0
    
    %% Input validation
    if ~isnumeric(T_jv) || ndims(T_jv) ~= 2
        error('module2_residual_transfer_function:invalid_dstf', ...
              'T_jv must be a 2D numeric array');
    end
    
    [n, p_dstf] = size(T_jv);
    
    if ~isnumeric(L) || ndims(L) ~= 2
        error('module2_residual_transfer_function:invalid_leadfield', ...
              'Leadfield matrix L must be a 2D numeric array');
    end
    
    [p, n_leadfield] = size(L);
    
    % Check dimension consistency
    if n ~= n_leadfield
        error('module2_residual_transfer_function:dimension_mismatch_sources', ...
              'T_jv sources (%d) must match leadfield sources (%d)', n, n_leadfield);
    end
    
    if p_dstf ~= p
        error('module2_residual_transfer_function:dimension_mismatch_sensors', ...
              'T_jv sensors (%d) must match leadfield sensors (%d)', p_dstf, p);
    end
    
    % Check for NaN or Inf values
    if any(isnan(T_jv(:))) || any(isinf(T_jv(:)))
        error('module2_residual_transfer_function:invalid_dstf_values', ...
              'T_jv contains NaN or Inf values');
    end
    
    if any(isnan(L(:))) || any(isinf(L(:)))
        error('module2_residual_transfer_function:invalid_leadfield_values', ...
              'Leadfield matrix contains NaN or Inf values');
    end
    
    %% Parse optional arguments
    if nargin < 3
        options = struct();
    end
    
    % Set default options
    default_options = struct(...
        'enforce_identity_diagonal', true, ...
        'numerical_tolerance', 1e-12, ...
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
    if ~islogical(options.enforce_identity_diagonal) || ~isscalar(options.enforce_identity_diagonal)
        error('module2_residual_transfer_function:invalid_diagonal_flag', ...
              'enforce_identity_diagonal must be a logical scalar');
    end
    
    if ~isscalar(options.numerical_tolerance) || ~isnumeric(options.numerical_tolerance) || options.numerical_tolerance <= 0
        error('module2_residual_transfer_function:invalid_tolerance', ...
              'numerical_tolerance must be a positive scalar');
    end
    
    if ~islogical(options.verbose) || ~isscalar(options.verbose)
        error('module2_residual_transfer_function:invalid_verbose_flag', ...
              'verbose must be a logical scalar');
    end
    
    %% Main computation
    try
        if options.verbose
            fprintf('Computing Residual Transfer Function with dimensions: sensors=%d, sources=%d\n', p, n);
        end
        
        % Step 1: Compute the product L * T_jv
        L_T_jv = L * T_jv;
        
        if options.verbose
            fprintf('Computed L * T_jv with dimensions: %d×%d\n', size(L_T_jv, 1), size(L_T_jv, 2));
            
            % Check numerical properties of L * T_jv
            max_element = max(abs(L_T_jv(:)));
            frobenius_norm = norm(L_T_jv, 'fro');
            fprintf('L * T_jv properties: max_element=%.6f, frobenius_norm=%.6f\n', ...
                    max_element, frobenius_norm);
        end
        
        % Step 2: Compute the identity matrix
        I_p = eye(p);
        
        % Step 3: Compute the residual transfer function
        T_xi_v = I_p - L_T_jv;
        
        % Step 4: Optional diagonal enforcement for numerical stability
        if options.enforce_identity_diagonal
            % Ensure diagonal elements are exactly as expected
            diagonal_vals = diag(T_xi_v);
            expected_diagonal = 1 - diag(L_T_jv);
            
            diagonal_error = max(abs(diagonal_vals - expected_diagonal));
            if diagonal_error > options.numerical_tolerance
                if options.verbose
                    fprintf('Diagonal correction applied: max error = %.2e\n', diagonal_error);
                end
                
                % Correct diagonal elements
                T_xi_v = T_xi_v - diag(diag(T_xi_v)) + diag(expected_diagonal);
            end
        end
        
        if options.verbose
            fprintf('Residual Transfer Function computation completed successfully\n');
            fprintf('Output matrix dimensions: %d×%d\n', size(T_xi_v, 1), size(T_xi_v, 2));
        end
        
    catch ME
        % Enhanced error handling with context
        switch ME.identifier
            case 'MATLAB:matrix:singular'
                error('module2_residual_transfer_function:computation_error', ...
                      'Numerical error in computing L * T_jv');
            otherwise
                rethrow(ME);
        end
    end
    
    %% Output validation and property checks
    % Check for NaN or Inf in output
    if any(isnan(T_xi_v(:))) || any(isinf(T_xi_v(:)))
        error('module2_residual_transfer_function:invalid_output', ...
              'Output T_xi_v contains NaN or Inf values, indicating numerical instability');
    end
    
    % Verify output dimensions
    if size(T_xi_v, 1) ~= p || size(T_xi_v, 2) ~= p
        error('module2_residual_transfer_function:output_dimension_error', ...
              'Output dimensions (%d×%d) do not match expected (%d×%d)', ...
              size(T_xi_v, 1), size(T_xi_v, 2), p, p);
    end
    
    % Mathematical property validation
    if options.verbose
        % Check idempotent property: T_ξv * T_ξv = T_ξv
        T_xi_v_squared = T_xi_v * T_xi_v;
        idempotent_error = norm(T_xi_v_squared - T_xi_v, 'fro') / norm(T_xi_v, 'fro');
        
        % Check complementary projection property: T_ξv + L * T_jv = I
        complementary_sum = T_xi_v + L_T_jv;
        identity_error = norm(complementary_sum - I_p, 'fro') / norm(I_p, 'fro');
        
        % Compute matrix properties
        frobenius_norm = norm(T_xi_v, 'fro');
        max_element = max(abs(T_xi_v(:)));
        
        try
            % Eigenvalue analysis (for real matrices or when feasible)
            if isreal(T_xi_v)
                eigenvals = eig(T_xi_v);
                num_zero_eigenvals = sum(abs(eigenvals) < options.numerical_tolerance);
                num_one_eigenvals = sum(abs(eigenvals - 1) < options.numerical_tolerance);
            else
                eigenvals = [];
                num_zero_eigenvals = NaN;
                num_one_eigenvals = NaN;
            end
        catch
            eigenvals = [];
            num_zero_eigenvals = NaN;
            num_one_eigenvals = NaN;
        end
        
        fprintf('Mathematical property validation:\n');
        fprintf('  - Idempotent error: %.2e (should be ~0)\n', idempotent_error);
        fprintf('  - Complementary projection error: %.2e (should be ~0)\n', identity_error);
        fprintf('  - Frobenius norm: %.6f\n', frobenius_norm);
        fprintf('  - Maximum absolute element: %.6f\n', max_element);
        fprintf('  - Matrix is real: %s\n', mat2str(isreal(T_xi_v)));
        
        if ~isempty(eigenvals)
            fprintf('  - Eigenvalues ≈ 0: %d\n', num_zero_eigenvals);
            fprintf('  - Eigenvalues ≈ 1: %d\n', num_one_eigenvals);
            fprintf('  - Rank (approx): %d\n', p - num_zero_eigenvals);
        end
        
        % Warn about potential numerical issues
        if idempotent_error > 1e-10
            warning('module2_residual_transfer_function:idempotent_violation', ...
                    'Large idempotent property violation (%.2e), check numerical stability', idempotent_error);
        end
        
        if identity_error > 1e-10
            warning('module2_residual_transfer_function:complementary_violation', ...
                    'Large complementary projection violation (%.2e), check input consistency', identity_error);
        end
    end
end