function demo_results = demo_module1_preprocessing_updated()
% DEMO_MODULE1_PREPROCESSING_UPDATED - Updated demonstration with complex data support
%
% This function demonstrates the Module 1 preprocessing pipeline using
% the latest module7_simulation_improved_complex with complex Hermitian matrices.
% 
% MAJOR UPDATES:
% 1. Uses module7_simulation_improved_complex instead of module7_simulation
% 2. Added complex data validation and analysis
% 3. Enhanced error handling for complex matrix operations
% 4. Extended quality assessment for complex matrices
%
% Usage:
%   demo_results = demo_module1_preprocessing_updated()
%
% File location: examples/demo_module1_preprocessing_updated.m

    fprintf('========================================\n');
    fprintf('Module 1 Preprocessing Demo (Complex Data Support)\n');
    fprintf('========================================\n\n');
    
    % Initialize demo results
    demo_results = struct();
    demo_results.timestamp = datestr(now);
    demo_results.version = 'complex_support_v1.0';
    
    %% Generate simulation data using latest module7
    fprintf('=== Generating Complex Simulation Data ===\n');
    try
        % Use module7_simulation_improved_complex to generate complex Hermitian test data
        [true_prec, true_cov, emp_cov, sim_params] = module7_simulation_improved_complex(...
            'n_nodes', 12, ...
            'n_freq', 15, ...
            'n_samples', 100, ...
            'graph_type', 'chain', ...
            'complex_strength', 1.0, ...  % Enable complex components
            'sparsity_variation', 0.3, ...
            'edge_activation_smoothness', 0.8, ...
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
        
        % Analyze complex properties of generated data
        complex_analysis = analyze_complex_data(emp_cov, sim_params);
        
        demo_results.data_generation = struct();
        demo_results.data_generation.success = true;
        demo_results.data_generation.params = sim_params;
        demo_results.data_generation.complex_analysis = complex_analysis;
        
        fprintf('Complex data analysis:\n');
        fprintf('  - Matrices with complex entries: %d/%d\n', ...
                complex_analysis.matrices_with_complex, sim_params.n_freq);
        fprintf('  - Average complex fraction: %.3f\n', complex_analysis.avg_complex_fraction);
        fprintf('  - Max imaginary component: %.4f\n', complex_analysis.max_imag_component);
        fprintf('  - All matrices Hermitian: %s\n', logical_to_string(complex_analysis.all_hermitian));
        
    catch ME
        fprintf('Data generation failed: %s\n', ME.message);
        demo_results.data_generation = struct();
        demo_results.data_generation.success = false;
        demo_results.data_generation.error = ME.message;
        return;
    end
    
    %% Validate complex data compatibility
    fprintf('\n=== Validating Complex Data Compatibility ===\n');
    try
        compatibility_check = validate_complex_compatibility(input_data);
        demo_results.compatibility_check = compatibility_check;
        
        if ~compatibility_check.all_checks_passed
            fprintf('Warning: Some compatibility issues detected:\n');
            fields = fieldnames(compatibility_check);
            for i = 1:length(fields)
                field = fields{i};
                if islogical(compatibility_check.(field)) && ~compatibility_check.(field)
                    fprintf('  - %s: FAILED\n', strrep(field, '_', ' '));
                end
            end
        else
            fprintf('All compatibility checks passed!\n');
        end
        
    catch ME
        fprintf('Compatibility check failed: %s\n', ME.message);
        demo_results.compatibility_check = struct('success', false, 'error', ME.message);
    end
    
    %% Run preprocessing with enhanced error handling for complex data
    fprintf('\n=== Running Preprocessing (Complex Data Enhanced) ===\n');
    try
        % Run module1 preprocessing with parameters optimized for complex data
        preprocessing_results = module1_preprocessing_main(input_data, ...
            'smoothing_method', 'moving_average', ...
            'window_size', 5, ...
            'diagonal_loading', true, ...
            'loading_factor', 0.02, ...
            'target_diagonal', 1.0, ...
            'diagonal_tolerance', 0.1, ...
            'force_hermitian', true, ...      % Important for complex matrices
            'check_psd', true, ...            % Validate positive semi-definiteness
            'verbose', true);
        
        fprintf('Preprocessing completed successfully\n');
        
        % Store results with enhanced complex data analysis
        demo_results.preprocessing = struct();
        demo_results.preprocessing.success = true;
        demo_results.preprocessing.results = preprocessing_results;
        
        % Extract key metrics with complex data support
        if isfield(preprocessing_results, 'Sigma_tilde') && ...
           iscell(preprocessing_results.Sigma_tilde)
            demo_results.preprocessing.n_frequencies = length(preprocessing_results.Sigma_tilde);
            demo_results.preprocessing.matrix_size = size(preprocessing_results.Sigma_tilde{1});
            
            % Analyze complex properties of processed data
            processed_complex_analysis = analyze_complex_data(preprocessing_results.Sigma_tilde, sim_params);
            demo_results.preprocessing.processed_complex_analysis = processed_complex_analysis;
            
            fprintf('Processed data complex analysis:\n');
            fprintf('  - Matrices with complex entries: %d/%d\n', ...
                    processed_complex_analysis.matrices_with_complex, sim_params.n_freq);
            fprintf('  - Max imaginary component: %.4f\n', processed_complex_analysis.max_imag_component);
            fprintf('  - Hermitian preservation: %s\n', logical_to_string(processed_complex_analysis.all_hermitian));
        end
        
        % Extract whitening quality metrics
        if isfield(preprocessing_results, 'processing_stats') && ...
           isfield(preprocessing_results.processing_stats, 'whitening_quality')
            demo_results.preprocessing.whitening_quality = preprocessing_results.processing_stats.whitening_quality;
            
            % Compute quality summary with complex data considerations
            quality = preprocessing_results.processing_stats.whitening_quality;
            demo_results.preprocessing.quality_summary = compute_complex_quality_summary(quality);
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
        
        % Try to provide helpful debugging information
        if contains(ME.message, 'complex') || contains(ME.message, 'Hermitian')
            fprintf('\nSuggested fixes for complex data issues:\n');
            fprintf('1. Ensure CovarianceWhitening class supports complex matrices\n');
            fprintf('2. Check if diagonal smoothing handles complex data correctly\n');
            fprintf('3. Verify Hermitian matrix operations in whitening process\n');
        end
    end
    
    %% Generate comprehensive summary
    fprintf('\n=== Generating Enhanced Summary ===\n');
    try
        demo_results.summary = generate_enhanced_summary(demo_results);
        
        fprintf('Summary generated:\n');
        fprintf('  - Overall success: %s\n', logical_to_string(demo_results.summary.overall_success));
        fprintf('  - Complex data handling: %s\n', logical_to_string(demo_results.summary.complex_data_success));
        fprintf('  - Processing quality: %.1f%%\n', demo_results.summary.processing_quality_score);
        
    catch ME
        fprintf('Summary generation failed: %s\n', ME.message);
        demo_results.summary = struct('success', false, 'error', ME.message);
    end
    
    %% Visualization with complex data support
    fprintf('\n=== Creating Enhanced Visualization ===\n');
    try
        if demo_results.preprocessing.success
            % Use enhanced visualization function that supports complex data
            create_complex_data_visualization(demo_results);
            demo_results.visualization = struct('success', true, 'created', true);
        else
            fprintf('Skipping visualization due to preprocessing failure\n');
            demo_results.visualization = struct('success', false, 'reason', 'preprocessing_failed');
        end
        
    catch ME
        fprintf('Visualization failed: %s\n', ME.message);
        demo_results.visualization = struct('success', false, 'error', ME.message);
    end
    
    fprintf('\n========================================\n');
    fprintf('Module 1 Complex Data Demo Complete!\n');
    fprintf('========================================\n');
    
    % Final recommendations
    if demo_results.preprocessing.success
        fprintf('✓ Demo completed successfully with complex data support\n');
        fprintf('✓ Ready for integration with complex Module 7 data\n');
    else
        fprintf('❌ Issues detected - review error messages and compatibility checks\n');
        fprintf('→ Additional module1 modifications may be required\n');
    end

end

%% Helper functions for complex data analysis

function complex_analysis = analyze_complex_data(matrix_cell_array, params)
% Analyze complex properties of matrix data
    
    complex_analysis = struct();
    n_freq = length(matrix_cell_array);
    
    matrices_with_complex = 0;
    total_complex_fraction = 0;
    max_imag_component = 0;
    all_hermitian = true;
    
    for f = 1:n_freq
        matrix = matrix_cell_array{f};
        
        % Check for complex entries
        has_complex = any(abs(imag(matrix(:))) > 1e-12);
        if has_complex
            matrices_with_complex = matrices_with_complex + 1;
            complex_fraction = sum(abs(imag(matrix(:))) > 1e-12) / numel(matrix);
            total_complex_fraction = total_complex_fraction + complex_fraction;
        end
        
        % Track maximum imaginary component
        max_imag_component = max(max_imag_component, max(abs(imag(matrix(:)))));
        
        % Check Hermitian property
        hermitian_error = max(abs(matrix - matrix'));
        if hermitian_error > 1e-10
            all_hermitian = false;
        end
    end
    
    complex_analysis.matrices_with_complex = matrices_with_complex;
    complex_analysis.avg_complex_fraction = total_complex_fraction / max(matrices_with_complex, 1);
    complex_analysis.max_imag_component = max_imag_component;
    complex_analysis.all_hermitian = all_hermitian;
    complex_analysis.n_frequencies = n_freq;
end

function compatibility = validate_complex_compatibility(input_data)
% Validate that input data is compatible with complex processing
    
    compatibility = struct();
    
    try
        % Check if empirical covariance matrices exist
        compatibility.has_empirical_cov = isfield(input_data.sim_results, 'Sigma_emp') && ...
                                         iscell(input_data.sim_results.Sigma_emp);
        
        % Check matrix dimensions consistency
        if compatibility.has_empirical_cov
            Sigma_emp = input_data.sim_results.Sigma_emp;
            n_freq = length(Sigma_emp);
            
            compatibility.consistent_dimensions = true;
            for f = 1:n_freq
                if size(Sigma_emp{f}, 1) ~= size(Sigma_emp{f}, 2)
                    compatibility.consistent_dimensions = false;
                    break;
                end
            end
            
            % Check for complex data
            compatibility.contains_complex_data = false;
            for f = 1:n_freq
                if any(abs(imag(Sigma_emp{f}(:))) > 1e-12)
                    compatibility.contains_complex_data = true;
                    break;
                end
            end
            
            % Check Hermitian property
            compatibility.matrices_hermitian = true;
            for f = 1:n_freq
                hermitian_error = max(abs(Sigma_emp{f} - Sigma_emp{f}'));
                if hermitian_error > 1e-10
                    compatibility.matrices_hermitian = false;
                    break;
                end
            end
            
        else
            compatibility.consistent_dimensions = false;
            compatibility.contains_complex_data = false;
            compatibility.matrices_hermitian = false;
        end
        
        % Overall compatibility check
        compatibility.all_checks_passed = compatibility.has_empirical_cov && ...
                                         compatibility.consistent_dimensions && ...
                                         compatibility.matrices_hermitian;
        
    catch ME
        compatibility.validation_error = ME.message;
        compatibility.all_checks_passed = false;
    end
end

function quality_summary = compute_complex_quality_summary(quality_metrics)
% Compute quality summary with complex data considerations
    
    quality_summary = struct();
    
    try
        if isfield(quality_metrics, 'diagonal_errors')
            quality_summary.avg_diagonal_error = mean(quality_metrics.diagonal_errors);
            quality_summary.max_diagonal_error = max(quality_metrics.diagonal_errors);
        end
        
        if isfield(quality_metrics, 'hermitian_errors')
            quality_summary.avg_hermitian_error = mean(quality_metrics.hermitian_errors);
            quality_summary.max_hermitian_error = max(quality_metrics.hermitian_errors);
        end
        
        if isfield(quality_metrics, 'condition_numbers')
            quality_summary.avg_condition_number = mean(quality_metrics.condition_numbers);
            quality_summary.max_condition_number = max(quality_metrics.condition_numbers);
        end
        
        % Overall quality score (0-100)
        diagonal_score = max(0, 100 - quality_summary.max_diagonal_error * 1000);
        hermitian_score = max(0, 100 - quality_summary.max_hermitian_error * 1000);
        condition_score = max(0, 100 - log10(quality_summary.max_condition_number) * 10);
        
        quality_summary.overall_score = (diagonal_score + hermitian_score + condition_score) / 3;
        
    catch ME
        quality_summary.computation_error = ME.message;
        quality_summary.overall_score = 0;
    end
end

function summary = generate_enhanced_summary(demo_results)
% Generate comprehensive summary including complex data assessment
    
    summary = struct();
    
    % Overall success assessment
    summary.overall_success = demo_results.data_generation.success && ...
                             demo_results.preprocessing.success;
    
    % Complex data handling assessment
    summary.complex_data_success = false;
    if isfield(demo_results, 'data_generation') && demo_results.data_generation.success
        if isfield(demo_results.data_generation, 'complex_analysis')
            complex_analysis = demo_results.data_generation.complex_analysis;
            summary.complex_data_success = complex_analysis.matrices_with_complex > 0 && ...
                                          complex_analysis.all_hermitian;
        end
    end
    
    % Processing quality score
    summary.processing_quality_score = 0;
    if demo_results.preprocessing.success && ...
       isfield(demo_results.preprocessing, 'quality_summary')
        if isfield(demo_results.preprocessing.quality_summary, 'overall_score')
            summary.processing_quality_score = demo_results.preprocessing.quality_summary.overall_score;
        end
    end
    
    % Timing performance
    if isfield(demo_results.preprocessing, 'timing')
        summary.total_processing_time = demo_results.preprocessing.timing.total;
    end
    
    % Recommendations
    summary.recommendations = {};
    
    if ~summary.overall_success
        summary.recommendations{end+1} = 'Fix preprocessing pipeline errors';
    end
    
    if ~summary.complex_data_success
        summary.recommendations{end+1} = 'Enhance complex data handling capabilities';
    end
    
    if summary.processing_quality_score < 80
        summary.recommendations{end+1} = 'Improve numerical stability and accuracy';
    end
    
    if isempty(summary.recommendations)
        summary.recommendations{end+1} = 'System ready for production use with complex data';
    end
end

function str = logical_to_string(logical_value)
% Convert logical value to readable string
    if logical_value
        str = 'YES';
    else
        str = 'NO';
    end
end

function create_complex_data_visualization(demo_results)
% Create enhanced visualization for complex data using ComplexDataVisualizer class
    
    fprintf('Creating complex data visualization...\n');
    
    try
        % Use the ComplexDataVisualizer class
        visualizer = ComplexDataVisualizer();
        visualizer.visualize_results(demo_results);
        fprintf('Enhanced complex visualization completed\n');
        
    catch ME
        fprintf('Enhanced visualization failed: %s\n', ME.message);
        
        % Try fallback to original visualization
        try
            visualize_module1_results(demo_results);
            fprintf('Fallback to standard visualization completed\n');
        catch ME2
            fprintf('Standard visualization also failed: %s\n', ME2.message);
            fprintf('Visualization enhancement needed\n');
        end
    end
end