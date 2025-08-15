function preprocessing_results = module1_preprocessing_main(input_data, varargin)
% MODULE1_PREPROCESSING_MAIN - Complete data preparation and preprocessing pipeline
%
% This function implements the complete Module 1 pipeline for data preparation
% and preprocessing, including data acquisition, diagonal smoothing, whitening
% matrix construction, and covariance whitening.
%
% UPDATED: Now uses CovarianceWhitening class for enhanced functionality
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
    
    % Data acquisition parameters
    addParameter(p, 'smoothing_method', 'moving_average', @ischar);
    addParameter(p, 'window_size', 5, @(x) isscalar(x) && x > 0);
    addParameter(p, 'diagonal_loading', true, @islogical);
    addParameter(p, 'loading_factor', 0.01, @(x) isscalar(x) && x > 0);
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
    preprocessing_results.timing = struct();
    preprocessing_results.processing_stats = struct();
    
    % Create output directory if needed
    if params.save_intermediate && ~exist(params.output_dir, 'dir')
        mkdir(params.output_dir);
    end
    
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
        if strcmp(input_data.mode, 'simulation')
            Sigma_emp = data_acquisition('simulation', input_data.sim_results);
        elseif strcmp(input_data.mode, 'eeg_data')
            Sigma_emp = data_acquisition('eeg_data', input_data.data_path, input_data.params);
        else
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
        save(fullfile(params.output_dir, 'step1_empirical_covariances.mat'), 'Sigma_emp');
    end
    
    %% Step 2: Diagonal Smoothing
    if params.verbose
        fprintf('Step 2/4: Diagonal Smoothing\n');
        fprintf('----------------------------\n');
    end
    
    step2_start = tic;
    try
        [g_smooth, Sigma_emp_loaded, smoothing_stats] = diagonal_smoothing(Sigma_emp, ...
            'method', params.smoothing_method, ...
            'window_size', params.window_size, ...
            'diagonal_loading', params.diagonal_loading, ...
            'loading_factor', params.loading_factor, ...
            'verbose', params.verbose);
        
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
        [D, whitening_stats] = whitening_matrix_construction(g_smooth, ...
            'min_power', params.min_power, ...
            'verbose', params.verbose);
        
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
    
    %% Step 4: Covariance Whitening (UPDATED TO USE CLASS)
    if params.verbose
        fprintf('Step 4/4: Covariance Whitening\n');
        fprintf('------------------------------\n');
    end
    
    step4_start = tic;
    try
        % UPDATED: Use CovarianceWhitening class instead of function
        [Sigma_tilde, whitening_quality] = CovarianceWhitening.whiten(Sigma_emp, D, ...
            'target_diagonal', params.target_diagonal, ...
            'diagonal_tolerance', params.diagonal_tolerance, ...
            'force_hermitian', params.force_hermitian, ...
            'check_psd', params.check_psd, ...
            'verbose', params.verbose);
        
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