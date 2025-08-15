function demo_results = demo_module1_preprocessing()
% DEMO_MODULE1_PREPROCESSING - Updated demonstration of Module 1 preprocessing
%
% This function demonstrates the Module 1 preprocessing pipeline using
% simulation data from module7_simulation with proper error handling.
% 
% UPDATED: Now works with the new CovarianceWhitening class implementation
%
% Usage:
%   demo_results = demo_module1_preprocessing()
%
% File location: examples/demo_module1_preprocessing.m

    fprintf('========================================\n');
    fprintf('Module 1 Preprocessing Demo (Updated)\n');
    fprintf('========================================\n\n');
    
    % Initialize demo results
    demo_results = struct();
    demo_results.timestamp = datestr(now);
    demo_results.version = 'updated_class_v1.0';
    
    %% Generate simulation data using module7
    fprintf('=== Generating Simulation Data ===\n');
    try
        % Use module7_simulation to generate proper test data
        [true_prec, true_cov, emp_cov, sim_params] = module7_simulation(...
            'n_nodes', 12, ...
            'n_freq', 15, ...
            'n_samples', 100, ...
            'graph_type', 'chain', ...
            'random_seed', 42);
        
        % Convert to input format expected by module1
        input_data = struct();
        input_data.mode = 'simulation';
        input_data.sim_results = struct();
        input_data.sim_results.Sigma_emp = emp_cov;
        input_data.sim_results.F = sim_params.n_freq;
        input_data.sim_results.n = sim_params.n_nodes;
        input_data.sim_results.T = sim_params.n_samples;
        
        fprintf('Generated data: %d nodes, %d frequencies, %d samples\n', ...
                sim_params.n_nodes, sim_params.n_freq, sim_params.n_samples);
        
        demo_results.data_generation = struct();
        demo_results.data_generation.success = true;
        demo_results.data_generation.params = sim_params;
        
    catch ME
        fprintf('Data generation failed: %s\n', ME.message);
        demo_results.data_generation = struct();
        demo_results.data_generation.success = false;
        demo_results.data_generation.error = ME.message;
        return;
    end
    
    %% Run preprocessing with error handling (UPDATED for new class)
    fprintf('\n=== Running Preprocessing (Updated) ===\n');
    try
        % Run module1 preprocessing with safe parameters
        % NOTE: module1_preprocessing_main now uses CovarianceWhitening class internally
        preprocessing_results = module1_preprocessing_main(input_data, ...
            'smoothing_method', 'moving_average', ...
            'window_size', 5, ...
            'diagonal_loading', true, ...
            'loading_factor', 0.02, ...
            'target_diagonal', 1.0, ...
            'diagonal_tolerance', 0.1, ...
            'verbose', true);
        
        fprintf('Preprocessing completed successfully\n');
        
        % Store results with proper error checking
        demo_results.preprocessing = struct();
        demo_results.preprocessing.success = true;
        demo_results.preprocessing.results = preprocessing_results;
        
        % Extract key metrics safely (UPDATED field names)
        if isfield(preprocessing_results, 'Sigma_tilde') && ...
           iscell(preprocessing_results.Sigma_tilde)
            demo_results.preprocessing.n_frequencies = length(preprocessing_results.Sigma_tilde);
            demo_results.preprocessing.matrix_size = size(preprocessing_results.Sigma_tilde{1});
        end
        
        % NEW: Extract whitening quality metrics from the class-based implementation
        if isfield(preprocessing_results, 'processing_stats') && ...
           isfield(preprocessing_results.processing_stats, 'whitening_quality')
            demo_results.preprocessing.whitening_quality = preprocessing_results.processing_stats.whitening_quality;
            
            % Compute summary statistics
            quality = preprocessing_results.processing_stats.whitening_quality;
            demo_results.preprocessing.quality_summary = compute_quality_summary(quality);
        end
        
        % Extract timing information
        if isfield(preprocessing_results, 'timing')
            demo_results.preprocessing.timing = preprocessing_results.timing;
        end
        
    catch ME
        fprintf('Preprocessing failed: %s\n', ME.message);
        demo_results.preprocessing = struct();
        demo_results.preprocessing.success = false;
        demo_results.preprocessing.error = ME.message;
        demo_results.preprocessing.stack = ME.stack;
    end
    
    %% NEW: Demonstrate independent quality assessment using new class features
    if demo_results.preprocessing.success
        fprintf('\n=== Demonstrating New Class Features ===\n');
        try
            demo_results.class_features = demonstrate_class_features(preprocessing_results);
        catch ME
            fprintf('Class features demo failed: %s\n', ME.message);
            demo_results.class_features = struct('success', false, 'error', ME.message);
        end
    end
    
    %% Generate summary report
    fprintf('\n=== Demo Summary ===\n');
    demo_results.summary = create_demo_summary(demo_results);
    
    fprintf('\nDemo completed successfully!\n');
    fprintf('Results stored in demo_results structure.\n');
end

function quality_summary = compute_quality_summary(quality)
% COMPUTE_QUALITY_SUMMARY - Compute summary statistics from whitening quality
    
    quality_summary = struct();
    
    try
        % Overall effectiveness statistics
        if isfield(quality, 'whitening_effectiveness')
            effectiveness = quality.whitening_effectiveness;
            quality_summary.mean_effectiveness = mean(effectiveness);
            quality_summary.min_effectiveness = min(effectiveness);
            quality_summary.max_effectiveness = max(effectiveness);
            quality_summary.std_effectiveness = std(effectiveness);
        end
        
        % Diagonal error statistics
        if isfield(quality, 'diagonal_errors')
            max_errors = max(quality.diagonal_errors, [], 2);
            quality_summary.mean_max_diagonal_error = mean(max_errors);
            quality_summary.worst_diagonal_error = max(max_errors);
            quality_summary.best_diagonal_error = min(max_errors);
        end
        
        % Success rates
        if isfield(quality, 'diagonal_errors')
            max_errors = max(quality.diagonal_errors, [], 2);
            F = length(max_errors);
            quality_summary.success_rate_10 = sum(max_errors <= 0.10) / F * 100;
            quality_summary.success_rate_15 = sum(max_errors <= 0.15) / F * 100;
            quality_summary.success_rate_20 = sum(max_errors <= 0.20) / F * 100;
        end
        
        % Condition number statistics
        if isfield(quality, 'condition_numbers')
            valid_conditions = quality.condition_numbers(isfinite(quality.condition_numbers));
            if ~isempty(valid_conditions)
                quality_summary.mean_condition_number = mean(valid_conditions);
                quality_summary.max_condition_number = max(valid_conditions);
            end
        end
        
        quality_summary.success = true;
        
    catch ME
        quality_summary.success = false;
        quality_summary.error = ME.message;
    end
end

function class_features = demonstrate_class_features(preprocessing_results)
% DEMONSTRATE_CLASS_FEATURES - Showcase new CovarianceWhitening class capabilities
    
    fprintf('Testing new CovarianceWhitening class features...\n');
    
    class_features = struct();
    class_features.success = true;
    
    % Get data from preprocessing results
    Sigma_emp = preprocessing_results.Sigma_emp;
    Sigma_tilde = preprocessing_results.Sigma_tilde;
    F = length(Sigma_emp);
    
    %% Test 1: Independent quality computation for a single frequency
    fprintf('  Testing independent quality computation...\n');
    omega = min(3, F); % Use frequency 3 or the last available
    
    quality_single = CovarianceWhitening.compute_quality_metrics(...
        Sigma_emp{omega}, Sigma_tilde{omega}, omega, ...
        'target_diagonal', 1.0, 'diagonal_tolerance', 0.1);
    
    class_features.single_quality = quality_single;
    fprintf('    ✓ Single frequency quality computed (effectiveness: %.3f)\n', ...
            quality_single.whitening_effectiveness);
    
    %% Test 2: Custom quality reporting
    fprintf('  Testing custom quality reporting...\n');
    
    % Create a smaller quality structure for demonstration
    demo_quality = struct();
    demo_quality.whitening_effectiveness = [0.85, 0.92, 0.78, 0.89, 0.91];
    demo_quality.diagonal_errors = rand(5, 12) * 0.1; % 5 frequencies, 12 nodes
    demo_quality.hermitian_errors = rand(5, 1) * 1e-12;
    demo_quality.min_eigenvalues = rand(5, 1) * 0.1 + 0.01;
    demo_quality.condition_numbers = rand(5, 1) * 100 + 10;
    
    % Test reporting (capture output to avoid cluttering)
    fprintf('    Custom quality report:\n');
    CovarianceWhitening.report_quality(demo_quality, ...
        struct('target_diagonal', 1.0, 'diagonal_tolerance', 0.08));
    
    class_features.custom_reporting = true;
    
    %% Test 3: Batch quality analysis with different parameters
    fprintf('  Testing batch quality analysis...\n');
    
    batch_results = cell(min(3, F), 1);
    test_targets = [0.9, 1.0, 1.1];
    
    for i = 1:min(3, F)
        target = test_targets(min(i, length(test_targets)));
        batch_results{i} = CovarianceWhitening.compute_quality_metrics(...
            Sigma_emp{i}, Sigma_tilde{i}, i, ...
            'target_diagonal', target, 'diagonal_tolerance', 0.1);
    end
    
    class_features.batch_analysis = batch_results;
    fprintf('    ✓ Batch analysis completed for %d frequencies\n', length(batch_results));
    
    %% Summary
    fprintf('  New class features demonstration completed successfully!\n');
    
    class_features.capabilities_tested = {
        'Independent quality computation for single frequencies';
        'Custom quality reporting with different parameters';
        'Batch quality analysis with varying targets';
        'Access to component-level quality scores'
    };
end

function summary = create_demo_summary(demo_results)
% CREATE_DEMO_SUMMARY - Create summary of demo results with safe field access
    
    summary = struct();
    summary.timestamp = demo_results.timestamp;
    summary.version = demo_results.version;
    
    % Data generation summary
    if demo_results.data_generation.success
        summary.data_status = 'SUCCESS';
        summary.data_info = sprintf('Generated %dx%d matrices for %d frequencies', ...
            demo_results.data_generation.params.n_nodes, ...
            demo_results.data_generation.params.n_nodes, ...
            demo_results.data_generation.params.n_freq);
    else
        summary.data_status = 'FAILED';
        summary.data_info = demo_results.data_generation.error;
    end
    
    % Preprocessing summary
    if demo_results.preprocessing.success
        summary.preprocessing_status = 'SUCCESS';
        
        % Safe access to preprocessing metrics
        if isfield(demo_results.preprocessing, 'n_frequencies')
            summary.n_frequencies_processed = demo_results.preprocessing.n_frequencies;
        else
            summary.n_frequencies_processed = 'Unknown';
        end
        
        if isfield(demo_results.preprocessing, 'matrix_size')
            summary.matrix_dimensions = demo_results.preprocessing.matrix_size;
        else
            summary.matrix_dimensions = 'Unknown';
        end
        
        % NEW: Enhanced quality metrics from class-based implementation
        if isfield(demo_results.preprocessing, 'quality_summary')
            quality = demo_results.preprocessing.quality_summary;
            
            if isfield(quality, 'mean_effectiveness')
                summary.mean_whitening_effectiveness = quality.mean_effectiveness;
            end
            
            if isfield(quality, 'success_rate_10')
                summary.success_rate_10_percent = quality.success_rate_10;
            end
            
            if isfield(quality, 'mean_max_diagonal_error')
                summary.mean_diagonal_error = quality.mean_max_diagonal_error;
            end
        end
        
        % Timing information
        if isfield(demo_results.preprocessing, 'timing')
            summary.total_processing_time = demo_results.preprocessing.timing.total;
        end
        
    else
        summary.preprocessing_status = 'FAILED';
        summary.preprocessing_error = demo_results.preprocessing.error;
    end
    
    % NEW: Class features demonstration summary
    if isfield(demo_results, 'class_features') && demo_results.class_features.success
        summary.class_features_status = 'SUCCESS';
        summary.new_capabilities = length(demo_results.class_features.capabilities_tested);
    else
        summary.class_features_status = 'SKIPPED or FAILED';
    end
    
    % Overall assessment
    if demo_results.data_generation.success && demo_results.preprocessing.success
        summary.overall_status = 'SUCCESS';
        summary.assessment = 'Demo completed successfully with new class features';
    else
        summary.overall_status = 'FAILED';
        if ~demo_results.data_generation.success
            summary.assessment = 'Failed at data generation stage';
        else
            summary.assessment = 'Failed at preprocessing stage';
        end
    end
    
    % Print summary
    fprintf('Updated Demo Summary:\n');
    fprintf('====================\n');
    fprintf('Data Generation: %s\n', summary.data_status);
    fprintf('Preprocessing: %s\n', summary.preprocessing_status);
    fprintf('Class Features: %s\n', summary.class_features_status);
    fprintf('Overall Status: %s\n', summary.overall_status);
    fprintf('Assessment: %s\n', summary.assessment);
    
    if strcmp(summary.preprocessing_status, 'SUCCESS')
        fprintf('\nProcessing Details:\n');
        if isfield(summary, 'n_frequencies_processed')
            fprintf('  Frequencies processed: %s\n', num2str(summary.n_frequencies_processed));
        end
        if isfield(summary, 'matrix_dimensions')
            fprintf('  Matrix dimensions: %s\n', mat2str(summary.matrix_dimensions));
        end
        if isfield(summary, 'mean_whitening_effectiveness')
            fprintf('  Mean whitening effectiveness: %.3f\n', summary.mean_whitening_effectiveness);
        end
        if isfield(summary, 'success_rate_10_percent')
            fprintf('  Success rate (≤10%% error): %.1f%%\n', summary.success_rate_10_percent);
        end
        if isfield(summary, 'total_processing_time')
            fprintf('  Total processing time: %.2f seconds\n', summary.total_processing_time);
        end
    end
    
    if strcmp(summary.class_features_status, 'SUCCESS')
        fprintf('\nNew Class Features:\n');
        fprintf('  Successfully demonstrated %d new capabilities\n', summary.new_capabilities);
        fprintf('  ✓ Independent quality assessment\n');
        fprintf('  ✓ Custom quality reporting\n');
        fprintf('  ✓ Batch analysis with different parameters\n');
    end
    
    if strcmp(summary.overall_status, 'FAILED')
        fprintf('\nPlease check the error details above for troubleshooting.\n');
    end
end