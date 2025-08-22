function S_xi_xi = module2_residual_empirical_covariance(T_xi_v, S_vv, options)
    % MODULE2_RESIDUAL_EMPIRICAL_COVARIANCE - Compute Residual Empirical Covariance (REC)
    %
    % Syntax:
    %   S_xi_xi = module2_residual_empirical_covariance(T_xi_v, S_vv)
    %   S_xi_xi = module2_residual_empirical_covariance(T_xi_v, S_vv, options)
    %
    % Description:
    %   Computes the Residual Empirical Covariance (REC) for the E-step:
    %   S_ξξ^(ω) = T_ξv^(ω) * S_vv^(ω) * (T_ξv^(ω))^H
    %   
    %   This represents the empirical covariance in the residual space after
    %   removing estimated source contributions. The REC is a key component
    %   for estimating the precision matrix in the EM algorithm.
    %
    % Input Arguments:
    %   T_xi_v - (complex, p×p) Residual Transfer Function matrix
    %   S_vv - (complex, p×p) Empirical sensor covariance matrix
    %
    % Name-Value Arguments:
    %   enforce_hermitian - (logical) Force output to be Hermitian. Default: true
    %   regularization_factor - (double) Ridge regularization for numerical stability. Default: 0
    %   min_eigenvalue_threshold - (double) Minimum eigenvalue for positive definiteness. Default: 1e-12
    %   numerical_tolerance - (double) Tolerance for numerical checks. Default: 1e-12
    %   verbose - (logical) Enable verbose output for debugging. Default: false
    %
    % Output Arguments:
    %   S_xi_xi - (complex, p×p) Residual empirical covariance matrix
    %
    % Examples:
    %   % Basic usage in E-step pipeline
    %   L = randn(64, 100);  % 64 sensors, 100 sources
    %   S_vv = randn(64) + 1i*randn(64); S_vv = S_vv * S_vv';  % Empirical covariance
    %   Sigma_jj = eye(100) * 0.5;  % Source prior covariance
    %   Sigma_xi = eye(64) * 0.1;   % Sensor noise covariance
    %   
    %   T_jv = module2_dstf_computation(L, Sigma_jj, Sigma_xi);
    %   T_xi_v = module2_residual_transfer_function(T_jv, L);
    %   S_xi_xi = module2_residual_empirical_covariance(T_xi_v, S_vv);
    %   
    %   % Usage with regularization for numerical stability
    %   S_xi_xi = module2_residual_empirical_covariance(T_xi_v, S_vv, ...
    %                                                    'regularization_factor', 1e-6, ...
    %                                                    'verbose', true);
    %
    % Mathematical Background:
    %   The REC represents the empirical covariance of the residual process
    %   after optimal source estimation. Key properties:
    %   - S_ξξ is Hermitian positive semi-definite
    %   - rank(S_ξξ) ≤ rank(S_vv) (dimensionality reduction possible)
    %   - S_ξξ → S_vv as source estimation quality decreases
    %   - S_ξξ → 0 as source estimation becomes perfect
    %
    % Numerical Considerations:
    %   The triple matrix product can be numerically unstable, especially
    %   when T_ξv has eigenvalues near zero. Regularization options help
    %   maintain numerical stability and positive definiteness.
    %
    % See also: MODULE2_RESIDUAL_TRANSFER_FUNCTION, MODULE2_PRECISION_MATRIX_COMPUTATION
    %
    % Author: [Author Name]
    % Date: [Current Date]
    % Version: 1.0
    
    %% Input validation
    if ~isnumeric(T_xi_v) || ~ismatrix(T_xi_v)
        error('module2_residual_empirical_covariance:invalid_transfer_function', ...
              'T_xi_v must be a numeric matrix');
    end
    
    [p1, p2] = size(T_xi_v);
    if p1 ~= p2
        error('module2_residual_empirical_covariance:not_square_transfer', ...
              'Residual transfer function must be square, got %d×%d', p1, p2);
    end
    p = p1;
    
    if ~isnumeric(S_vv) || ~ismatrix(S_vv)
        error('module2_residual_empirical_covariance:invalid_empirical_covariance', ...
              'S_vv must be a numeric matrix');
    end
    
    if size(S_vv, 1) ~= p || size(S_vv, 2) ~= p
        error('module2_residual_empirical_covariance:dimension_mismatch', ...
              'S_vv dimensions (%d×%d) must match T_xi_v dimensions (%d×%d)', ...
              size(S_vv, 1), size(S_vv, 2), p, p);
    end
    
    % Check for valid covariance matrix (Hermitian)
    if ~ishermitian(S_vv)
        hermitian_error = norm(S_vv - S_vv', 'fro') / norm(S_vv, 'fro');
        if hermitian_error > 1e-10
            error('module2_residual_empirical_covariance:not_hermitian_covariance', ...
                  'Empirical covariance is not Hermitian (relative error: %.2e)', hermitian_error);
        end
    end
    
    % Check for NaN or Inf values
    if any(isnan(T_xi_v(:))) || any(isinf(T_xi_v(:)))
        error('module2_residual_empirical_covariance:invalid_transfer_values', ...
              'T_xi_v contains NaN or Inf values');
    end
    
    if any(isnan(S_vv(:))) || any(isinf(S_vv(:)))
        error('module2_residual_empirical_covariance:invalid_covariance_values', ...
              'S_vv contains NaN or Inf values');
    end
    
    %% Parse optional arguments
    if nargin < 3
        options = struct();
    end
    
    % Set default options
    default_options = struct(...
        'enforce_hermitian', true, ...
        'regularization_factor', 0, ...
        'min_eigenvalue_threshold', 1e-12, ...
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
    if ~islogical(options.enforce_hermitian) || ~isscalar(options.enforce_hermitian)
        error('module2_residual_empirical_covariance:invalid_hermitian_flag', ...
              'enforce_hermitian must be a logical scalar');
    end
    
    if ~isscalar(options.regularization_factor) || ~isnumeric(options.regularization_factor) || options.regularization_factor < 0
        error('module2_residual_empirical_covariance:invalid_regularization', ...
              'regularization_factor must be a non-negative scalar');
    end
    
    if ~isscalar(options.min_eigenvalue_threshold) || ~isnumeric(options.min_eigenvalue_threshold) || options.min_eigenvalue_threshold < 0
        error('module2_residual_empirical_covariance:invalid_eigenvalue_threshold', ...
              'min_eigenvalue_threshold must be a non-negative scalar');
    end
    
    if ~isscalar(options.numerical_tolerance) || ~isnumeric(options.numerical_tolerance) || options.numerical_tolerance <= 0
        error('module2_residual_empirical_covariance:invalid_tolerance', ...
              'numerical_tolerance must be a positive scalar');
    end
    
    if ~islogical(options.verbose) || ~isscalar(options.verbose)
        error('module2_residual_empirical_covariance:invalid_verbose_flag', ...
              'verbose must be a logical scalar');
    end
    
    %% Main computation
    try
        if options.verbose
            fprintf('Computing Residual Empirical Covariance with dimensions: %d×%d\n', p, p);
            
            % Analyze input properties
            transfer_norm = norm(T_xi_v, 'fro');
            covariance_norm = norm(S_vv, 'fro');
            fprintf('Input properties: ||T_xi_v||_F = %.6f, ||S_vv||_F = %.6f\n', ...
                    transfer_norm, covariance_norm);
        end
        
        % Step 1: Apply regularization to empirical covariance if requested
        if options.regularization_factor > 0
            trace_S_vv = trace(S_vv);
            reg_strength = options.regularization_factor * trace_S_vv / p;
            S_vv_reg = S_vv + reg_strength * eye(p);
            
            if options.verbose
                fprintf('Applied ridge regularization: strength = %.2e\n', reg_strength);
            end
        else
            S_vv_reg = S_vv;
        end
        
        % Step 2: Compute the triple matrix product
        % S_ξξ = T_ξv * S_vv * T_ξv^H
        
        % Efficient computation: first compute T_ξv * S_vv, then multiply by T_ξv^H
        intermediate = T_xi_v * S_vv_reg;
        S_xi_xi = intermediate * T_xi_v';
        
        % Always enforce Hermitian symmetry for numerical stability
        S_xi_xi = (S_xi_xi + S_xi_xi') / 2;
        
        % Step 4: Handle potential negative eigenvalues for numerical stability
        if options.min_eigenvalue_threshold > 0
            try
                [V, D] = eig(S_xi_xi);
                eigenvals = real(diag(D));
                
                negative_eigenvals = sum(eigenvals < 0);
                small_eigenvals = sum(eigenvals < options.min_eigenvalue_threshold & eigenvals >= 0);
                
                if negative_eigenvals > 0 || small_eigenvals > 0
                    if options.verbose
                        fprintf('Eigenvalue correction: %d negative, %d below threshold\n', ...
                                negative_eigenvals, small_eigenvals);
                    end
                    
                    % Clip eigenvalues to maintain positive semi-definiteness
                    eigenvals(eigenvals < options.min_eigenvalue_threshold) = options.min_eigenvalue_threshold;
                    
                    % Reconstruct matrix
                    S_xi_xi = V * diag(eigenvals) * V';
                    
                    % Ensure Hermitian after reconstruction
                    if options.enforce_hermitian
                        S_xi_xi = (S_xi_xi + S_xi_xi') / 2;
                    end
                    
                    if options.verbose
                        fprintf('Applied eigenvalue clipping, min eigenvalue: %.2e\n', min(eigenvals));
                    end
                end
            catch
                if options.verbose
                    fprintf('Skipped eigenvalue correction due to numerical issues\n');
                end
            end
        end
        
        if options.verbose
            fprintf('Residual Empirical Covariance computation completed successfully\n');
        end
        
    catch ME
        % Enhanced error handling with context
        switch ME.identifier
            case 'MATLAB:nomem'
                error('module2_residual_empirical_covariance:memory_error', ...
                      'Insufficient memory for triple matrix product. Consider reducing matrix dimensions');
            case 'MATLAB:matrix:singular'
                error('module2_residual_empirical_covariance:singular_matrix', ...
                      'Numerical singularity encountered in computation');
            otherwise
                rethrow(ME);
        end
    end
    
    %% Output validation
    % Check for NaN or Inf in output
    if any(isnan(S_xi_xi(:))) || any(isinf(S_xi_xi(:)))
        error('module2_residual_empirical_covariance:invalid_output', ...
              'Output S_xi_xi contains NaN or Inf values, indicating numerical instability');
    end
    
    % Verify output dimensions
    if size(S_xi_xi, 1) ~= p || size(S_xi_xi, 2) ~= p
        error('module2_residual_empirical_covariance:output_dimension_error', ...
              'Output dimensions (%d×%d) do not match expected (%d×%d)', ...
              size(S_xi_xi, 1), size(S_xi_xi, 2), p, p);
    end
    
    % Mathematical property validation
    if options.verbose
        % Compute matrix properties
        frobenius_norm = norm(S_xi_xi, 'fro');
        max_element = max(abs(S_xi_xi(:)));
        trace_value = trace(S_xi_xi);
        
        try
            min_eigenval = min(real(eig(S_xi_xi)));
            condition_number = cond(S_xi_xi);
            matrix_rank = rank(S_xi_xi, options.numerical_tolerance);
        catch
            min_eigenval = NaN;
            condition_number = NaN;
            matrix_rank = NaN;
        end
        
        % Compute dimensionality reduction ratio
        try
            rank_S_vv = rank(S_vv, options.numerical_tolerance);
            rank_ratio = matrix_rank / rank_S_vv;
        catch
            rank_ratio = NaN;
        end
        
        fprintf('Output validation:\n');
        fprintf('  - Frobenius norm: %.6f\n', frobenius_norm);
        fprintf('  - Maximum absolute element: %.6f\n', max_element);
        fprintf('  - Trace: %.6f\n', trace_value);
        fprintf('  - Minimum eigenvalue: %.2e\n', min_eigenval);
        fprintf('  - Condition number: %.2e\n', condition_number);
        fprintf('  - Matrix rank: %d\n', matrix_rank);
        fprintf('  - Rank ratio (S_ξξ/S_vv): %.3f\n', rank_ratio);
        fprintf('  - Is Hermitian: %s\n', mat2str(ishermitian(S_xi_xi)));
        fprintf('  - Is complex: %s\n', mat2str(~isreal(S_xi_xi)));
        
        % Validate positive semi-definiteness
        if min_eigenval < -options.numerical_tolerance
            warning('module2_residual_empirical_covariance:not_positive_semidefinite', ...
                    'Output matrix is not positive semi-definite (min eigenvalue: %.2e)', min_eigenval);
        end
        
        % Check for significant dimensionality reduction
        if ~isnan(rank_ratio) && rank_ratio < 0.9
            fprintf('  - Significant dimensionality reduction detected (%.1f%% rank retention)\n', ...
                    rank_ratio * 100);
        end
    end
end