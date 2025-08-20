        % Analyze complex properties of generated data - FIXED
        complex_analysis = analyze_complex_data_safe(emp_cov);
        
        demo_results.data_generation = struct();
        demo_results.data_generation.success = true;
        demo_results.data_generation.params = sim_params;
        demo_results.data_generation.complex_analysis = complex_analysis;
        
        fprintf('Complex data analysis:\n');
        fprintf('  - Matrices with complex entries: %d/%d\n', ...
                complex_analysisfunction demo_results = demo_module1_preprocessing_updated()
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
% 5. FIXED: Quality calculation to prevent Inf% values
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
        
        % Analyze complex properties of generated data - FIXED
        complex_analysis = analyze_complex_data_safe(emp_cov);
        
        demo_results.data_generation = struct();
        demo_results.data_generation.success = true;
        demo_results.data_generation.params = sim_params;
        demo_results.data_generation.complex_analysis = complex_analysis;
        
        fprintf('Complex data analysis:\n');
        fprintf('  - Matrices with complex entries: %d/%d\n', ...
                complex_analysis.matrices_with_complex, length(emp_cov));
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
        compatibility = validate_complex_compatibility(input_data);
        
        if compatibility.all_checks_passed
            fprintf('All compatibility checks passed!\n');
        else
            fprintf('Some compatibility issues detected:\n');
            if ~compatibility.has_empirical_cov
                fprintf('  - Missing empirical covariance matrices\n');
            end
            if ~compatibility.consistent_dimensions
                fprintf('  - Inconsistent matrix dimensions\n');
            end
            if ~compatibility.matrices_hermitian
                fprintf('  - Non-Hermitian matrices detected\n');
            end
        end
        
    catch ME
        fprintf('Compatibility validation failed: %s\n', ME.message);
        compatibility = struct('all_checks_passed', false, 'error', ME.message);
    end
    
    %% Run preprocessing pipeline
    fprintf('\n=== Running Preprocessing (Complex Data Enhanced) ===\n');
    try
        % module1_preprocessing_main now uses CovarianceWhitening class internally
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
            
            % Analyze complex properties of processed data - FIXED
            processed_complex_analysis = analyze_complex_data_safe(preprocessing_results.Sigma_tilde);
            demo_results.preprocessing.processed_complex_analysis = processed_complex_analysis;
            
            fprintf('Processed data complex analysis:\n');
            fprintf('  - Matrices with complex entries: %d/%d\n', ...
                    processed_complex_analysis.matrices_with_complex, length(preprocessing_results.Sigma_tilde));
            fprintf('  - Max imaginary component: %.4f\n', processed_complex_analysis.max_imag_component);
            fprintf('  - Hermitian preservation: %s\n', logical_to_string(processed_complex_analysis.all_hermitian));
        end
        
        % NEW: Extract whitening quality metrics from the class-based implementation
        if isfield(preprocessing_results, 'processing_stats') && ...
           isfield(preprocessing_results.processing_stats, 'whitening_quality')
            demo_results.preprocessing.whitening_quality = preprocessing_results.processing_stats.whitening_quality;
            
            % Compute summary statistics
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
        fprintf('❌ Consider running test_module1_fix() to verify system integrity\n');
    end
end

function complex_analysis = analyze_complex_data(matrices, params)
% Analyze complex properties of matrix data
    
    complex_analysis = struct();
    n_freq = length(matrices);
    
    % Count matrices with complex entries
    matrices_with_complex = 0;
    total_complex_fraction = 0;
    max_imag_component = 0;
    hermitian_errors = zeros(n_freq, 1);
    
    for f = 1:n_freq
        matrix = matrices{f};
        
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
        hermitian_error = max(abs(matrix(:) - conj(matrix').'));
        hermitian_errors(f) = hermitian_error;
    end
    
    % Compile analysis results
    complex_analysis.matrices_with_complex = matrices_with_complex;
    complex_analysis.avg_complex_fraction = total_complex_fraction / max(matrices_with_complex, 1);
    complex_analysis.max_imag_component = max_imag_component;
    complex_analysis.hermitian_errors = hermitian_errors;
    complex_analysis.max_hermitian_error = max(hermitian_errors);
    complex_analysis.all_hermitian = all(hermitian_errors < 1e-10);
end

function compatibility = validate_complex_compatibility(input_data)
% Validate compatibility with complex data processing
    
    compatibility = struct();
    
    try
        % Check for empirical covariance matrices
        compatibility.has_empirical_cov = isfield(input_data, 'sim_results') && ...
                                         isfield(input_data.sim_results, 'Sigma_emp');
        
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
% FIXED: Compute quality summary with complex data considerations and proper Inf handling
    
    quality_summary = struct();
    
    try
        if isfield(quality_metrics, 'diagonal_errors')
            quality_summary.avg_diagonal_error = mean(quality_metrics.diagonal_errors);
            quality_summary.max_diagonal_error = max(quality_metrics.diagonal_errors);
        else
            quality_summary.avg_diagonal_error = NaN;
            quality_summary.max_diagonal_error = NaN;
        end
        
        if isfield(quality_metrics, 'hermitian_errors')
            quality_summary.avg_hermitian_error = mean(quality_metrics.hermitian_errors);
            quality_summary.max_hermitian_error = max(quality_metrics.hermitian_errors);
        else
            quality_summary.avg_hermitian_error = NaN;
            quality_summary.max_hermitian_error = NaN;
        end
        
        if isfield(quality_metrics, 'condition_numbers')
            quality_summary.avg_condition_number = mean(quality_metrics.condition_numbers);
            quality_summary.max_condition_number = max(quality_metrics.condition_numbers);
        else
            quality_summary.avg_condition_number = NaN;
            quality_summary.max_condition_number = NaN;
        end
        
        % FIXED: Overall quality score (0-100) with proper Inf/NaN handling
        
        % Diagonal score - handle NaN/Inf
        if isfinite(quality_summary.max_diagonal_error) && quality_summary.max_diagonal_error >= 0
            diagonal_score = max(0, 100 - quality_summary.max_diagonal_error * 1000);
        else
            diagonal_score = 0;
        end
        
        % Hermitian score - handle NaN/Inf  
        if isfinite(quality_summary.max_hermitian_error) && quality_summary.max_hermitian_error >= 0
            hermitian_score = max(0, 100 - quality_summary.max_hermitian_error * 1000);
        else
            hermitian_score = 0;
        end
        
        % Condition score - FIXED to handle Inf/NaN properly
        if isfinite(quality_summary.max_condition_number) && quality_summary.max_condition_number > 0
            % Use safe log10 calculation
            log_cond = log10(quality_summary.max_condition_number);
            if isfinite(log_cond)
                condition_score = max(0, 100 - log_cond * 10);
            else
                condition_score = 0;
            end
        else
            condition_score = 0;
        end
        
        % Final score - ensure it's finite
        quality_summary.overall_score = (diagonal_score + hermitian_score + condition_score) / 3;
        
        % Safety check - ensure the result is finite
        if ~isfinite(quality_summary.overall_score)
            quality_summary.overall_score = 0;
        end
        
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
% FIXED: Create enhanced visualization for complex data with better error handling
    
    fprintf('Creating complex data visualization...\n');
    
    try
        % Use the ComplexDataVisualizer class if available
        if exist('ComplexDataVisualizer', 'class')
            visualizer = ComplexDataVisualizer();
            visualizer.visualize_results(demo_results);
            fprintf('Enhanced complex visualization completed\n');
        else
            % Try fallback to original visualization
            visualize_module1_results(demo_results);
            fprintf('Fallback to standard visualization completed\n');
        end
        
    catch ME
        fprintf('Enhanced visualization failed: %s\n', ME.message);
        
        % Try fallback to original visualization
        try
            visualize_module1_results(demo_results);
            fprintf('Fallback to standard visualization completed\n');
        catch ME2
            fprintf('Standard visualization also failed: %s\n', ME2.message);
            
            % Final fallback - basic status report
            try
                create_fallback_visualization(demo_results);
                fprintf('Basic fallback visualization completed\n');
            catch ME3
                fprintf('All visualization methods failed: %s\n', ME3.message);
                create_basic_status_report(demo_results);
            end
        end
    end
    
    fprintf('Complex data visualization process completed\n');
end

function create_fallback_visualization(demo_results)
% Create minimal fallback visualization
    
    figure('Name', 'Module 1 Status Report', 'Position', [300, 300, 600, 400]);
    
    % Create status summary
    status_lines = {};
    status_lines{end+1} = 'Module 1 Demo Results';
    status_lines{end+1} = '=====================';
    status_lines{end+1} = '';
    
    if demo_results.data_generation.success
        status_lines{end+1} = '✓ Data Generation: SUCCESS';
    else
        status_lines{end+1} = '✗ Data Generation: FAILED';
    end
    
    if demo_results.preprocessing.success
        status_lines{end+1} = '✓ Preprocessing: SUCCESS';
    else
        status_lines{end+1} = '✗ Preprocessing: FAILED';
    end
    
    text(0.1, 0.9, status_lines, 'VerticalAlignment', 'top', 'FontSize', 12, ...
         'FontName', 'FixedWidth');
    axis off;
    title('Demo Status Report', 'FontSize', 14, 'FontWeight', 'bold');
end

function create_basic_status_report(demo_results)
% Create text-only status report
    
    fprintf('\n=== Basic Status Report ===\n');
    if demo_results.data_generation.success
        fprintf('Data Generation: ✓ SUCCESS\n');
    else
        fprintf('Data Generation: ✗ FAILED\n');
    end
    
    if demo_results.preprocessing.success
        fprintf('Preprocessing: ✓ SUCCESS\n');
    else
        fprintf('Preprocessing: ✗ FAILED\n');
    end
    fprintf('=== End Status Report ===\n\n');
end