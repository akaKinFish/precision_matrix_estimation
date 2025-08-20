function preprocessing_results = module1_preprocessing_main(input_data, varargin)
% MODULE1_PREPROCESSING_MAIN - Complete data preparation and preprocessing pipeline
%
% This function implements the complete Module 1 pipeline for data preparation
% and preprocessing, including data acquisition, diagonal smoothing, whitening
% matrix construction, and covariance whitening.
%
% UPDATED: Now uses CovarianceWhitening class for enhanced functionality
% FIXED: Corrected diagonal_smoothing function call
%
% Usage:
%   results = module1_preprocessing_main(input_data)
%   results = module1_preprocessing_main(input_data, 'param', value, ...)
%
% Inputs:
%   input_data - Structure containing either:
%     For simulation mode:
%       .mode = 'simulation'
%       .sim_results = results from Module 7
%     For EEG data mode:
%       .mode = 'eeg_data'
%       .data_path = path to EEG data file
%       .params = EEG processing parameters
%
%   Optional parameters (name-value pairs):
%     'smoothing_method'    - Diagonal smoothing method (default: 'moving_average')
%     'window_size'         - Smoothing window size (default: 5)
%     'diagonal_loading'    - Enable diagonal loading (default: true)
%     'loading_factor'      - Diagonal loading factor (default: 0.01)
%     'min_power'          - Minimum power threshold (default: 1e-10)
%     'target_diagonal'    - Target diagonal value after whitening (default: 1.0)
%     'diagonal_tolerance' - Tolerance for diagonal errors (default: 0.1)
%     'force_hermitian'    - Force Hermitian symmetry (default: true)
%     'check_psd'          - Check positive semi-definiteness (default: true)
%     'save_intermediate'  - Save intermediate results (default: false)
%     'output_dir'         - Output directory for saving (default: './results/')
%     'verbose'            - Display detailed progress (default: true)
%
% Outputs:
%   preprocessing_results - Structure containing:
%     .Sigma_emp         - Original empirical covariance matrices {F x 1}
%     .Sigma_emp_loaded  - Diagonally loaded covariances (if enabled) {F x 1}
%     .g_smooth          - Smoothed diagonal powers {F x 1}
%     .D                 - Whitening matrices {F x 1}
%     .Sigma_tilde       - Whitened covariance matrices {F x 1}
%     .processing_stats  - Detailed processing statistics
%     .parameters        - Parameters used for processing
%     .timing           - Timing information for each step
%
% File location: src/modules/module1/module1_preprocessing_main.m

    % Start timing
    total_start_time = tic;
    
    % Parse input arguments
    p = inputParser;
    addRequired(p, 'input_data', @isstruct);
    addParameter(p, 'smoothing_method', 'moving_average', @ischar);
    addParameter(p, 'window_size', 5, @(x) isscalar(x) && x > 0);
    addParameter(p, 'diagonal_loading', true, @islogical);
    addParameter(p, 'loading_factor', 0.01, @(x) isscalar(x) && x >= 0);
    addParameter(p, 'min_power', 1e-10, @(x) isscalar(x) && x > 0);
    addParameter(p, 'target_diagonal', 1.0, @(x) isscalar(x) && x > 0);
    addParameter(p, 'diagonal_tolerance', 0.1, @(x) isscalar(x) && x > 0);
    addParameter(p, 'force_hermitian', true, @islogical);
    addParameter(p, 'check_psd', true, @islogical);
    addParameter(p, 'save_intermediate', false, @islogical);
    addParameter(p, 'output_dir', './results/', @ischar);
    addParameter(p, 'verbose', true, @islogical);
    
    parse(p, input_data, varargin{:});
    params = p.Results;
    
    % Initialize results structure
    preprocessing_results = struct();
    preprocessing_results.parameters = params;
    preprocessing_results.processing_stats = struct();
    preprocessing_results.timing = struct();
    
    if params.verbose
        fprintf('========================================\n');
        fprintf('Module 1 Preprocessing Pipeline\n');
        fprintf('========================================\n\n');
    end
    
    %% Step 1: Data Acquisition
    if params.verbose
        fprintf('Step 1/4: Data Acquisition\n');
        fprintf('--------------------------\n');
    end
    
    step1_start = tic;
    try
        switch input_data.mode
            case 'simulation'
                if ~isfield(input_data, 'sim_results')
                    error('simulation mode requires sim_results field');
                end
                
                sim_results = input_data.sim_results;
                
                % Validate required fields
                required_fields = {'Sigma_emp', 'F', 'n', 'T'};
                for i = 1:length(required_fields)
                    if ~isfield(sim_results, required_fields{i})
                        error('Missing required field: %s', required_fields{i});
                    end
                end
                
                Sigma_emp = sim_results.Sigma_emp;
                
                if ~iscell(Sigma_emp)
                    error('Sigma_emp must be a cell array');
                end
                
                if params.verbose
                    fprintf('Loaded simulation data: %d nodes, %d frequencies, %d samples per frequency\n', ...
                            sim_results.n, sim_results.F, sim_results.T);
                    
                    % Check for complex data
                    has_complex = false;
                    for f = 1:length(Sigma_emp)
                        if any(abs(imag(Sigma_emp{f}(:))) > 1e-12)
                            has_complex = true;
                            break;
                        end
                    end
                    
                    if has_complex
                        fprintf('Complex data detected - using enhanced processing\n');
                    end
                end
                
                % Validate matrices
                for f = 1:length(Sigma_emp)
                    matrix = Sigma_emp{f};
                    
                    if size(matrix, 1) ~= size(matrix, 2)
                        error('Matrix %d is not square', f);
                    end
                    
                    if size(matrix, 1) ~= sim_results.n
                        error('Matrix %d has wrong dimension', f);
                    end
                    
                    % Check Hermitian property for complex matrices
                    if any(abs(imag(matrix(:))) > 1e-12)
                        hermitian_error = max(abs(matrix - matrix'));
                        if hermitian_error > 1e-8
                            if params.verbose
                                fprintf('Warning: Matrix %d may not be Hermitian (error: %.2e)\n', f, hermitian_error);
                            end
                        end
                    end
                end
                
                if params.verbose
                    fprintf('Validation completed: %d matrices of size [%d x %d]\n', ...
                            length(Sigma_emp), size(Sigma_emp{1}, 1), size(Sigma_emp{1}, 2));
                    fprintf('Data acquisition completed. Loaded %d frequency points.\n', length(Sigma_emp));
                end
                
            case 'eeg_data'
                error('EEG data mode not implemented yet');
                
            otherwise
                error('module1_preprocessing_main:invalid_mode', ...
                      'Invalid data mode: %s', input_data.mode);
        end
        
        preprocessing_results.timing.data_acquisition = toc(step1_start);
        
        if params.verbose
            fprintf('Data acquisition completed in %.2f seconds\n', ...
                    preprocessing_results.timing.data_acquisition);
            fprintf('Acquired %d frequency matrices of size [%d x %d]\n\n', ...
                    length(Sigma_emp), size(Sigma_emp{1}, 1), size(Sigma_emp{1}, 2));
        end
        
    catch ME
        error('module1_preprocessing_main:data_acquisition_failed', ...
              'Data acquisition failed: %s', ME.message);
    end
    
    % Store results
    preprocessing_results.Sigma_emp = Sigma_emp;
    
    % Save intermediate result if requested
    if params.save_intermediate
        if ~exist(params.output_dir, 'dir')
            mkdir(params.output_dir);
        end
        save(fullfile(params.output_dir, 'step1_empirical_covariances.mat'), 'Sigma_emp');
    end
    
    %% Step 2: Diagonal Smoothing (FIXED)
    if params.verbose
        fprintf('Step 2/4: Diagonal Smoothing\n');
        fprintf('----------------------------\n');
    end
    
    step2_start = tic;
    try
        % FIXED: diagonal_smoothing now only returns 2 output arguments and corrected parameters
        [g_smooth, Sigma_emp_loaded] = diagonal_smoothing(Sigma_emp, ...
            'smoothing_method', params.smoothing_method, ...
            'window_size', params.window_size, ...
            'diagonal_loading', params.diagonal_loading, ...
            'loading_factor', params.loading_factor);
        
        % Create compatible smoothing_stats structure
        smoothing_stats = struct();
        smoothing_stats.method = params.smoothing_method;
        smoothing_stats.window_size = params.window_size;
        smoothing_stats.diagonal_loading = params.diagonal_loading;
        smoothing_stats.loading_factor = params.loading_factor;
        smoothing_stats.success = true;
        smoothing_stats.timestamp = datestr(now);
        
        preprocessing_results.timing.diagonal_smoothing = toc(step2_start);
        
        if params.verbose
            fprintf('Diagonal smoothing completed in %.2f seconds\n\n', ...
                    preprocessing_results.timing.diagonal_smoothing);
        end
        
    catch ME
        error('module1_preprocessing_main:diagonal_smoothing_failed', ...
              'Diagonal smoothing failed: %s', ME.message);
    end
    
    % Store results
    preprocessing_results.g_smooth = g_smooth;
    preprocessing_results.Sigma_emp_loaded = Sigma_emp_loaded;
    preprocessing_results.processing_stats.diagonal_smoothing = smoothing_stats;
    
    % Save intermediate result if requested
    if params.save_intermediate
        save(fullfile(params.output_dir, 'step2_smoothed_powers.mat'), ...
             'g_smooth', 'Sigma_emp_loaded', 'smoothing_stats');
    end
    
    %% Step 3: Whitening Matrix Construction
    if params.verbose
        fprintf('Step 3/4: Whitening Matrix Construction\n');
        fprintf('---------------------------------------\n');
    end
    
    step3_start = tic;
    try
        % Check if whitening_matrix_construction function exists
        if exist('whitening_matrix_construction', 'file')
            [D, whitening_stats] = whitening_matrix_construction(g_smooth, ...
                'min_power', params.min_power, ...
                'verbose', params.verbose);
        else
            % Fallback implementation
            [D, whitening_stats] = construct_whitening_matrices_fallback(g_smooth, params);
        end
        
        preprocessing_results.timing.whitening_construction = toc(step3_start);
        
        if params.verbose
            fprintf('Whitening matrix construction completed in %.2f seconds\n\n', ...
                    preprocessing_results.timing.whitening_construction);
        end
        
    catch ME
        error('module1_preprocessing_main:whitening_construction_failed', ...
              'Whitening matrix construction failed: %s', ME.message);
    end
    
    % Store results
    preprocessing_results.D = D;
    preprocessing_results.processing_stats.whitening_construction = whitening_stats;
    
    % Save intermediate result if requested
    if params.save_intermediate
        save(fullfile(params.output_dir, 'step3_whitening_matrices.mat'), ...
             'D', 'whitening_stats');
    end
    
    %% Step 4: Covariance Whitening (ENHANCED for complex data)
    if params.verbose
        fprintf('Step 4/4: Covariance Whitening\n');
        fprintf('------------------------------\n');
    end
    
    step4_start = tic;
    try
        % Check if CovarianceWhitening class exists
        if exist('CovarianceWhitening', 'file')
            % Use CovarianceWhitening class
            [Sigma_tilde, whitening_quality] = CovarianceWhitening.whiten(Sigma_emp, D, ...
                'target_diagonal', params.target_diagonal, ...
                'diagonal_tolerance', params.diagonal_tolerance, ...
                'force_hermitian', params.force_hermitian, ...
                'check_psd', params.check_psd, ...
                'verbose', params.verbose);
        else
            % Fallback implementation with enhanced complex support
            [Sigma_tilde, whitening_quality] = whiten_covariances_fallback(Sigma_emp, D, params);
        end
        
        preprocessing_results.timing.covariance_whitening = toc(step4_start);
        
        if params.verbose
            fprintf('Covariance whitening completed in %.2f seconds\n\n', ...
                    preprocessing_results.timing.covariance_whitening);
        end
        
    catch ME
        error('module1_preprocessing_main:covariance_whitening_failed', ...
              'Covariance whitening failed: %s', ME.message);
    end
    
    % Store results
    preprocessing_results.Sigma_tilde = Sigma_tilde;
    preprocessing_results.processing_stats.whitening_quality = whitening_quality;
    
    % Save intermediate result if requested
    if params.save_intermediate
        save(fullfile(params.output_dir, 'step4_whitened_covariances.mat'), ...
             'Sigma_tilde', 'whitening_quality');
    end
    
    %% Finalize results
    preprocessing_results.timing.total = toc(total_start_time);
    
    if params.verbose
        fprintf('========================================\n');
        fprintf('Module 1 Preprocessing Complete\n');
        fprintf('========================================\n');
        fprintf('Total processing time: %.2f seconds\n', preprocessing_results.timing.total);
        fprintf('\nTiming breakdown:\n');
        fprintf('  Data acquisition:         %.2f s\n', preprocessing_results.timing.data_acquisition);
        fprintf('  Diagonal smoothing:       %.2f s\n', preprocessing_results.timing.diagonal_smoothing);
        fprintf('  Whitening construction:   %.2f s\n', preprocessing_results.timing.whitening_construction);
        fprintf('  Covariance whitening:     %.2f s\n', preprocessing_results.timing.covariance_whitening);
        fprintf('\n');
    end
    
    % Save final results if requested
    if params.save_intermediate
        save(fullfile(params.output_dir, 'module1_preprocessing_results.mat'), ...
             'preprocessing_results');
        
        if params.verbose
            fprintf('Results saved to: %s\n', params.output_dir);
        end
    end

end

%% Fallback implementations

function [D, whitening_stats] = construct_whitening_matrices_fallback(g_smooth, params)
% Fallback implementation for whitening matrix construction
    
    F = length(g_smooth);
    D = cell(F, 1);
    
    whitening_stats = struct();
    whitening_stats.method = 'diagonal_fallback';
    whitening_stats.min_power = params.min_power;
    
    for f = 1:F
        powers = g_smooth{f};
        
        % Ensure minimum power
        powers = max(real(powers), params.min_power);
        
        % Create diagonal whitening matrix
        D{f} = diag(1 ./ sqrt(powers));
    end
    
    whitening_stats.success = true;
    whitening_stats.condition_numbers = cellfun(@cond, D);
    
    if params.verbose
        fprintf('Using fallback whitening matrix construction\n');
        fprintf('Condition number range: [%.2e, %.2e]\n', ...
                min(whitening_stats.condition_numbers), ...
                max(whitening_stats.condition_numbers));
    end
end

function [Sigma_tilde, whitening_quality] = whiten_covariances_fallback(Sigma_emp, D, params)
% Fallback implementation for covariance whitening with complex support
    
    F = length(Sigma_emp);
    Sigma_tilde = cell(F, 1);
    
    whitening_quality = struct();
    whitening_quality.method = 'enhanced_fallback';
    whitening_quality.diagonal_errors = zeros(F, 1);
    whitening_quality.hermitian_errors = zeros(F, 1);
    whitening_quality.condition_numbers = zeros(F, 1);
    
    for f = 1:F
        % Apply whitening transformation
        Sigma_tilde{f} = D{f} * Sigma_emp{f} * D{f}';
        
        % Force Hermitian if requested (important for complex matrices)
        if params.force_hermitian
            Sigma_tilde{f} = (Sigma_tilde{f} + Sigma_tilde{f}') / 2;
        end
        
        % Compute quality metrics
        diag_vals = diag(Sigma_tilde{f});
        whitening_quality.diagonal_errors(f) = max(abs(real(diag_vals) - params.target_diagonal));
        whitening_quality.hermitian_errors(f) = max(abs(Sigma_tilde{f} - Sigma_tilde{f}'));
        whitening_quality.condition_numbers(f) = cond(Sigma_tilde{f});
    end
    
    whitening_quality.success = true;
    whitening_quality.max_diagonal_error = max(whitening_quality.diagonal_errors);
    whitening_quality.max_hermitian_error = max(whitening_quality.hermitian_errors);
    whitening_quality.avg_condition_number = mean(whitening_quality.condition_numbers);
    
    if params.verbose
        fprintf('Using enhanced fallback covariance whitening\n');
        fprintf('Max diagonal error: %.4f\n', whitening_quality.max_diagonal_error);
        fprintf('Max Hermitian error: %.2e\n', whitening_quality.max_hermitian_error);
        fprintf('Avg condition number: %.2e\n', whitening_quality.avg_condition_number);
    end
end