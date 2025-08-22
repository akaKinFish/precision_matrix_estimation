function estep_results = module2_estep_main(input_data, estep_params)
    % MODULE2_ESTEP_MAIN - Complete E-step computation for EM algorithm (Full Source Model)
    %
    % Syntax:
    %   estep_results = module2_estep_main(input_data, estep_params)
    %
    % Description:
    %   Performs the complete E-step computation for the EM algorithm, including:
    %   1. Data-to-Source Transfer Function (DSTF) computation
    %   2. Posterior Source Covariance (SPC) computation
    %   3. Residual Transfer Function (RTF) computation
    %   4. Residual Empirical Covariance (REC) computation
    %   5. Initial precision matrix estimation with regularization
    %
    % Input Arguments:
    %   input_data - (struct) Input data structure containing:
    %     .leadfield_matrix - (double, p×n) Leadfield matrix L
    %     .empirical_covariances - (cell, F×1) Cell array of empirical covariance matrices
    %     .source_prior_covariances - (cell, F×1) Source prior covariance matrices
    %     .noise_covariance - (complex, p×p) Sensor noise covariance matrix
    %     .frequencies - (double, F×1) Frequency vector
    %   estep_params - (struct) E-step parameters containing:
    %     .regularization_factor - (double) Regularization strength. Default: 1e-8
    %     .condition_threshold - (double) Condition number threshold. Default: 1e12
    %     .min_eigenvalue_ratio - (double) Minimum eigenvalue ratio. Default: 1e-12
    %     .verbose - (logical) Enable verbose output. Default: false
    %
    % Output Arguments:
    %   estep_results - (struct) E-step computation results containing:
    %     .transfer_functions - (cell, F×1) Data-to-Source Transfer Functions
    %     .posterior_source_covariances - (cell, F×1) Posterior source covariances
    %     .residual_transfer_functions - (cell, F×1) Residual Transfer Functions
    %     .residual_covariances - (cell, F×1) Residual empirical covariances
    %     .initial_precision_matrices - (cell, F×1) Initial precision estimates
    %     .computation_stats - (struct) Computation statistics and diagnostics
    %     .success - (logical) Overall computation success flag
    %
    % Examples:
    %   % Using Module7 simulation data
    %   [true_prec, true_cov, emp_cov, sim_params] = module7_simulation_improved_complex();
    %   input_data.leadfield_matrix = randn(64, 100);
    %   input_data.empirical_covariances = emp_cov;
    %   input_data.source_prior_covariances = repmat({eye(100)*0.5}, length(emp_cov), 1);
    %   input_data.noise_covariance = eye(64) * 0.1;
    %   input_data.frequencies = 8:2:12;
    %   
    %   estep_params = struct('verbose', true, 'regularization_factor', 1e-6);
    %   results = module2_estep_main(input_data, estep_params);
    %
    % Mathematical Background:
    %   The E-step computes the expected sufficient statistics for the EM
    %   algorithm under the current parameter estimates. The computation
    %   follows the classical Kalman filter formulation for linear Gaussian
    %   state-space models, with extensions for frequency-domain analysis.
    %
    % See also: MODULE2_DSTF_COMPUTATION, MODULE2_POSTERIOR_SOURCE_COVARIANCE,
    %          MODULE2_RESIDUAL_TRANSFER_FUNCTION, MODULE2_RESIDUAL_EMPIRICAL_COVARIANCE
    %
    % Author: [Author Name]
    % Date: [Current Date]
    % Version: 2.0 - Full Source Model Only
    
    %% Input validation
    if nargin < 1
        error('module2_estep_main:insufficient_input', ...
              'At least input_data is required');
    end
    
    if nargin < 2
        estep_params = struct();
    end
    
    % Validate input_data structure
    required_fields = {'leadfield_matrix', 'empirical_covariances', 'source_prior_covariances', ...
                       'noise_covariance', 'frequencies'};
    for i = 1:length(required_fields)
        if ~isfield(input_data, required_fields{i})
            error('module2_estep_main:missing_field', ...
                  'Required field "%s" not found in input_data', required_fields{i});
        end
    end
    
    % Extract and validate main inputs
    L = input_data.leadfield_matrix;
    Sigma_emp = input_data.empirical_covariances;
    source_prior_covariances = input_data.source_prior_covariances;
    noise_covariance = input_data.noise_covariance;
    frequencies = input_data.frequencies;
    
    if ~isnumeric(L) || ndims(L) ~= 2
        error('module2_estep_main:invalid_leadfield', ...
              'leadfield_matrix must be a 2D numeric array');
    end
    
    [p, n] = size(L);
    
    if ~iscell(Sigma_emp)
        error('module2_estep_main:invalid_empirical_covariances', ...
              'empirical_covariances must be a cell array');
    end
    
    if ~iscell(source_prior_covariances)
        error('module2_estep_main:invalid_source_covariances', ...
              'source_prior_covariances must be a cell array');
    end
    
    F = length(Sigma_emp);
    
    if length(source_prior_covariances) ~= F
        error('module2_estep_main:source_covariance_mismatch', ...
              'Number of source covariances (%d) must match empirical covariances (%d)', ...
              length(source_prior_covariances), F);
    end
    
    if length(frequencies) ~= F
        error('module2_estep_main:frequency_mismatch', ...
              'Number of frequencies (%d) must match empirical covariances (%d)', ...
              length(frequencies), F);
    end
    
    if ~isnumeric(noise_covariance) || size(noise_covariance, 1) ~= p || size(noise_covariance, 2) ~= p
        error('module2_estep_main:invalid_noise_covariance', ...
              'noise_covariance must be %d×%d numeric matrix', p, p);
    end
    
    %% Parse and validate parameters
    default_params = struct(...
        'regularization_factor', 1e-8, ...
        'condition_threshold', 1e12, ...
        'min_eigenvalue_ratio', 1e-12, ...
        'verbose', false);
    
    % Merge user params with defaults
    param_names = fieldnames(default_params);
    for i = 1:length(param_names)
        param_name = param_names{i};
        if ~isfield(estep_params, param_name)
            estep_params.(param_name) = default_params.(param_name);
        end
    end
    
    % Validate parameter types
    if ~isscalar(estep_params.regularization_factor) || estep_params.regularization_factor < 0
        error('module2_estep_main:invalid_regularization', ...
              'regularization_factor must be a non-negative scalar');
    end
    
    %% Initialize results structure
    estep_results = struct();
    estep_results.transfer_functions = cell(F, 1);
    estep_results.posterior_source_covariances = cell(F, 1);
    estep_results.residual_transfer_functions = cell(F, 1);
    estep_results.residual_covariances = cell(F, 1);
    estep_results.initial_precision_matrices = cell(F, 1);
    estep_results.computation_stats = struct();
    estep_results.success = false;
    
    % Initialize computation statistics
    stats = struct();
    stats.processing_times = zeros(F, 5);  % 5 steps per frequency
    stats.condition_numbers = zeros(F, 3);  % 3 condition numbers per frequency
    stats.regularization_applied = false(F, 1);
    stats.eigenvalue_corrections = zeros(F, 2);  % posterior and residual corrections
    stats.total_computation_time = 0;
    stats.successful_frequencies = 0;
    stats.failed_frequencies = 0;
    stats.error_messages = {};
    
    %% Main computation loop
    overall_start_time = tic;
    
    if estep_params.verbose
        fprintf('========================================\n');
        fprintf('Module 2: E-Step Computation (Full Source Model)\n');
        fprintf('========================================\n');
        fprintf('Processing %d frequencies with %d sensors, %d sources\n', F, p, n);
        fprintf('----------------------------------------\n\n');
    end
    
    try
        for f = 1:F
            if estep_params.verbose
                fprintf('Processing frequency %d/%d (%.2f Hz)\n', f, F, frequencies(f));
            end
            
            freq_start_time = tic;
            
            try
                % Step 1: Validate empirical covariance for this frequency
                S_vv_f = Sigma_emp{f};
                Sigma_jj_f = source_prior_covariances{f};
                
                if ~isnumeric(S_vv_f) || size(S_vv_f, 1) ~= p || size(S_vv_f, 2) ~= p
                    error('Invalid empirical covariance at frequency %d', f);
                end
                
                if ~isnumeric(Sigma_jj_f) || size(Sigma_jj_f, 1) ~= n || size(Sigma_jj_f, 2) ~= n
                    error('Invalid source covariance at frequency %d', f);
                end
                
                % Step 1: Data-to-Source Transfer Function
                step1_time = tic;
                options_dstf = struct(...
                    'regularization_factor', estep_params.regularization_factor, ...
                    'condition_threshold', estep_params.condition_threshold, ...
                    'verbose', false);
                
                T_jv_f = module2_dstf_computation(L, Sigma_jj_f, noise_covariance, options_dstf);
                estep_results.transfer_functions{f} = T_jv_f;
                stats.processing_times(f, 1) = toc(step1_time);
                
                % Step 2: Posterior Source Covariance
                step2_time = tic;
                options_psc = struct(...
                    'regularization_factor', estep_params.regularization_factor, ...
                    'condition_threshold', estep_params.condition_threshold, ...
                    'min_eigenvalue_ratio', estep_params.min_eigenvalue_ratio, ...
                    'verbose', false);
                
                Sigma_jj_post_f = module2_posterior_source_covariance(Sigma_jj_f, L, noise_covariance, options_psc);
                estep_results.posterior_source_covariances{f} = Sigma_jj_post_f;
                stats.processing_times(f, 2) = toc(step2_time);
                
                % Step 3: Residual Transfer Function
                step3_time = tic;
                options_rtf = struct('verbose', false);
                
                T_xi_v_f = module2_residual_transfer_function(T_jv_f, L, options_rtf);
                estep_results.residual_transfer_functions{f} = T_xi_v_f;
                stats.processing_times(f, 3) = toc(step3_time);
                
                % Step 4: Residual Empirical Covariance
                step4_time = tic;
                options_rec = struct(...
                    'regularization_factor', estep_params.regularization_factor, ...
                    'min_eigenvalue_threshold', estep_params.min_eigenvalue_ratio * max(real(eig(S_vv_f))), ...
                    'verbose', false);
                
                S_xi_xi_f = module2_residual_empirical_covariance(T_xi_v_f, S_vv_f, options_rec);
                estep_results.residual_covariances{f} = S_xi_xi_f;
                stats.processing_times(f, 4) = toc(step4_time);
                
                % Step 5: Initial Precision Matrix with Regularization
                step5_time = tic;
                [V, D_eig] = eig(S_xi_xi_f);
                eigenvals = real(diag(D_eig));
                
                % Apply regularization as described in the derivation
                epsilon = estep_params.regularization_factor;
                d_max = max(eigenvals);
                eigenvals_reg = (eigenvals + epsilon * d_max) / (1 + epsilon);
                
                % Count eigenvalue corrections
                negative_count = sum(eigenvals < 0);
                small_count = sum(eigenvals < estep_params.min_eigenvalue_ratio * d_max & eigenvals >= 0);
                stats.eigenvalue_corrections(f, :) = [negative_count, small_count];
                
                % Compute stabilized precision matrix
                Omega_f = V * diag(1 ./ eigenvals_reg) * V';
                
                % Ensure Hermitian symmetry
                Omega_f = (Omega_f + Omega_f') / 2;
                
                estep_results.initial_precision_matrices{f} = Omega_f;
                stats.processing_times(f, 5) = toc(step5_time);
                
                % Record condition numbers
                try
                    stats.condition_numbers(f, 1) = cond(L * Sigma_jj_f * L' + noise_covariance);
                    stats.condition_numbers(f, 2) = cond(Sigma_jj_post_f);
                    stats.condition_numbers(f, 3) = cond(S_xi_xi_f);
                catch
                    stats.condition_numbers(f, :) = [NaN, NaN, NaN];
                end
                
                stats.regularization_applied(f) = epsilon > 0 || negative_count > 0 || small_count > 0;
                stats.successful_frequencies = stats.successful_frequencies + 1;
                
                if estep_params.verbose
                    total_freq_time = toc(freq_start_time);
                    fprintf('  ✓ Completed in %.3f seconds\n', total_freq_time);
                    
                    if stats.regularization_applied(f)
                        fprintf('  ℹ Regularization applied\n');
                    end
                end
                
            catch ME
                stats.failed_frequencies = stats.failed_frequencies + 1;
                stats.error_messages{end+1} = sprintf('Frequency %d: %s', f, ME.message);
                
                if estep_params.verbose
                    fprintf('  ✗ Failed: %s\n', ME.message);
                end
                
                % Fill with empty results for failed frequency
                estep_results.transfer_functions{f} = [];
                estep_results.posterior_source_covariances{f} = [];
                estep_results.residual_transfer_functions{f} = [];
                estep_results.residual_covariances{f} = [];
                estep_results.initial_precision_matrices{f} = [];
            end
        end
        
        stats.total_computation_time = toc(overall_start_time);
        
        % Determine overall success
        success_rate = stats.successful_frequencies / F;
        estep_results.success = success_rate >= 0.8;  % At least 80% success rate
        
        if estep_params.verbose
            fprintf('\n========================================\n');
            fprintf('E-Step Computation Summary\n');
            fprintf('========================================\n');
            fprintf('Total processing time: %.3f seconds\n', stats.total_computation_time);
            fprintf('Successful frequencies: %d/%d (%.1f%%)\n', ...
                    stats.successful_frequencies, F, success_rate * 100);
            
            if stats.failed_frequencies > 0
                fprintf('Failed frequencies: %d\n', stats.failed_frequencies);
                for i = 1:length(stats.error_messages)
                    fprintf('  - %s\n', stats.error_messages{i});
                end
            end
            
            avg_condition = nanmean(stats.condition_numbers(:));
            fprintf('Average condition number: %.2e\n', avg_condition);
            fprintf('Regularization applied: %d/%d frequencies\n', ...
                    sum(stats.regularization_applied), F);
            
            if estep_results.success
                fprintf('✓ E-step computation completed successfully\n');
            else
                fprintf('⚠ E-step computation completed with issues\n');
            end
            fprintf('========================================\n');
        end
        
    catch ME
        stats.total_computation_time = toc(overall_start_time);
        estep_results.success = false;
        
        if estep_params.verbose
            fprintf('\n✗ E-step computation failed: %s\n', ME.message);
        end
        
        rethrow(ME);
    end
    
    %% Store computation statistics
    estep_results.computation_stats = stats;
    
    %% Final validation
    if estep_results.success
        % Verify that all successful frequencies have valid precision matrices
        valid_precisions = 0;
        for f = 1:F
            if ~isempty(estep_results.initial_precision_matrices{f})
                Omega_f = estep_results.initial_precision_matrices{f};
                if ishermitian(Omega_f) && all(real(eig(Omega_f)) > -1e-10)
                    valid_precisions = valid_precisions + 1;
                end
            end
        end
        
        if valid_precisions < stats.successful_frequencies
            warning('module2_estep_main:invalid_precision_matrices', ...
                    'Some precision matrices are not positive definite');
        end
    end
end