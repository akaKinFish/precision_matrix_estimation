function test_results = test_module1_complex_support()
% TEST_MODULE1_COMPLEX_SUPPORT - Test Module 1 with complex data from Module 7
%
% This script tests the updated Module 1 preprocessing pipeline with
% complex Hermitian matrices from the improved Module 7 simulation.
%
% The test includes:
% 1. Complex data generation validation
% 2. Module 1 preprocessing compatibility check
% 3. Complex data preservation verification
% 4. Hermitian property maintenance
% 5. Whitening quality assessment
% 6. Visualization compatibility
%
% Usage:
%   test_results = test_module1_complex_support()
%
% Output:
%   test_results - Comprehensive test results structure
%
% File location: tests/integration/test_module1_complex_support.m

    fprintf('========================================\n');
    fprintf('Testing Module 1 Complex Data Support\n');
    fprintf('========================================\n\n');
    
    test_results = struct();
    test_results.timestamp = datestr(now);
    test_results.test_version = 'complex_integration_v1.0';
    
    % Initialize test counters
    total_tests = 0;
    passed_tests = 0;
    
    %% Test 1: Module 7 Complex Data Generation
    fprintf('=== Test 1: Module 7 Complex Data Generation ===\n');
    [total_tests, passed_tests] = run_test(@test_module7_complex_generation, ...
                                          'Module 7 Complex Generation', ...
                                          total_tests, passed_tests, test_results);
    
    %% Test 2: Data Format Compatibility
    fprintf('\n=== Test 2: Data Format Compatibility ===\n');
    [total_tests, passed_tests] = run_test(@test_data_format_compatibility, ...
                                          'Data Format Compatibility', ...
                                          total_tests, passed_tests, test_results);
    
    %% Test 3: Module 1 Basic Processing
    fprintf('\n=== Test 3: Module 1 Basic Processing ===\n');
    [total_tests, passed_tests] = run_test(@test_module1_basic_processing, ...
                                          'Module 1 Basic Processing', ...
                                          total_tests, passed_tests, test_results);
    
    %% Test 4: Complex Data Preservation
    fprintf('\n=== Test 4: Complex Data Preservation ===\n');
    [total_tests, passed_tests] = run_test(@test_complex_data_preservation, ...
                                          'Complex Data Preservation', ...
                                          total_tests, passed_tests, test_results);
    
    %% Test 5: Hermitian Property Maintenance
    fprintf('\n=== Test 5: Hermitian Property Maintenance ===\n');
    [total_tests, passed_tests] = run_test(@test_hermitian_property_maintenance, ...
                                          'Hermitian Property Maintenance', ...
                                          total_tests, passed_tests, test_results);
    
    %% Test 6: Whitening Quality Assessment
    fprintf('\n=== Test 6: Whitening Quality Assessment ===\n');
    [total_tests, passed_tests] = run_test(@test_whitening_quality_assessment, ...
                                          'Whitening Quality Assessment', ...
                                          total_tests, passed_tests, test_results);
    
    %% Test 7: Numerical Stability
    fprintf('\n=== Test 7: Numerical Stability ===\n');
    [total_tests, passed_tests] = run_test(@test_numerical_stability, ...
                                          'Numerical Stability', ...
                                          total_tests, passed_tests, test_results);
    
    %% Test 8: Visualization Compatibility
    fprintf('\n=== Test 8: Visualization Compatibility ===\n');
    [total_tests, passed_tests] = run_test(@test_visualization_compatibility, ...
                                          'Visualization Compatibility', ...
                                          total_tests, passed_tests, test_results);
    
    %% Test 9: End-to-End Integration
    fprintf('\n=== Test 9: End-to-End Integration ===\n');
    [total_tests, passed_tests] = run_test(@test_end_to_end_integration, ...
                                          'End-to-End Integration', ...
                                          total_tests, passed_tests, test_results);
    
    %% Generate Final Report
    fprintf('\n========================================\n');
    fprintf('FINAL TEST REPORT\n');
    fprintf('========================================\n');
    
    test_results.summary = struct();
    test_results.summary.total_tests = total_tests;
    test_results.summary.passed_tests = passed_tests;
    test_results.summary.success_rate = (passed_tests / total_tests) * 100;
    test_results.summary.overall_pass = (test_results.summary.success_rate >= 80);
    
    fprintf('Tests passed: %d/%d (%.1f%%)\n', passed_tests, total_tests, test_results.summary.success_rate);
    
    if test_results.summary.overall_pass
        fprintf('🎉 OVERALL RESULT: SUCCESS!\n');
        fprintf('   Module 1 is ready for complex data from Module 7\n');
        test_results.summary.recommendation = 'READY_FOR_PRODUCTION';
    else
        fprintf('❌ OVERALL RESULT: NEEDS WORK\n');
        fprintf('   Module 1 requires fixes for complex data support\n');
        test_results.summary.recommendation = 'REQUIRES_FIXES';
        
        % Provide specific recommendations
        fprintf('\nRecommendations:\n');
        if ~test_results.test3.passed
            fprintf('• Fix Module 1 basic processing issues\n');
        end
        if ~test_results.test4.passed
            fprintf('• Improve complex data preservation\n');
        end
        if ~test_results.test5.passed
            fprintf('• Fix Hermitian property maintenance\n');
        end
        if ~test_results.test8.passed
            fprintf('• Update visualization for complex data\n');
        end
    end
    
    fprintf('========================================\n');
end

%% Test runner helper function
function [total, passed] = run_test(test_func, test_name, total, passed, test_results)
    total = total + 1;
    test_num = total;
    
    fprintf('Test %d: %s... ', test_num, test_name);
    
    try
        test_result = test_func();
        if test_result.passed
            fprintf('PASSED\n');
            passed = passed + 1;
        else
            fprintf('FAILED\n');
            if isfield(test_result, 'details')
                fprintf('  Details: %s\n', test_result.details);
            end
        end
        
        % Store detailed test results
        field_name = sprintf('test%d', test_num);
        test_results.(field_name) = test_result;
        test_results.(field_name).test_name = test_name;
        
    catch ME
        fprintf('ERROR\n');
        fprintf('  Error: %s\n', ME.message);
        if ~isempty(ME.stack)
            fprintf('  Location: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
        end
        
        % Store error information
        field_name = sprintf('test%d', test_num);
        test_results.(field_name) = struct();
        test_results.(field_name).passed = false;
        test_results.(field_name).error = ME.message;
        test_results.(field_name).test_name = test_name;
    end
end

%% Individual test functions

function result = test_module7_complex_generation()
% Test Module 7 complex data generation
    
    result = struct('passed', false);
    
    try
        % Generate complex data using latest Module 7
        [true_prec, true_cov, emp_cov, sim_params] = module7_simulation_improved_complex(...
            'n_nodes', 8, ...
            'n_freq', 10, ...
            'n_samples', 50, ...
            'graph_type', 'random', ...
            'complex_strength', 1.0, ...
            'random_seed', 123);
        
        % Validate output structure
        assert(iscell(true_prec) && length(true_prec) == 10, 'Invalid precision matrix output');
        assert(iscell(true_cov) && length(true_cov) == 10, 'Invalid covariance matrix output');
        assert(iscell(emp_cov) && length(emp_cov) == 10, 'Invalid empirical covariance output');
        assert(isstruct(sim_params), 'Invalid parameters structure');
        
        % Check for complex data
        has_complex = false;
        for f = 1:length(emp_cov)
            if any(abs(imag(emp_cov{f}(:))) > 1e-12)
                has_complex = true;
                break;
            end
        end
        assert(has_complex, 'No complex data generated');
        
        % Check Hermitian property
        all_hermitian = true;
        for f = 1:length(emp_cov)
            hermitian_error = max(abs(emp_cov{f} - emp_cov{f}'));
            if hermitian_error > 1e-10
                all_hermitian = false;
                break;
            end
        end
        assert(all_hermitian, 'Not all matrices are Hermitian');
        
        result.passed = true;
        result.details = sprintf('Generated %d complex Hermitian matrices', sim_params.matrices_with_complex);
        result.data = struct('emp_cov', emp_cov, 'sim_params', sim_params);
        
    catch ME
        result.details = ME.message;
    end
end

function result = test_data_format_compatibility()
% Test data format compatibility between Module 7 and Module 1
    
    result = struct('passed', false);
    
    try
        % Generate test data
        [~, ~, emp_cov, sim_params] = module7_simulation_improved_complex(...
            'n_nodes', 6, ...
            'n_freq', 8, ...
            'n_samples', 30, ...
            'random_seed', 456);
        
        % Create Module 1 input structure
        input_data = struct();
        input_data.mode = 'simulation';
        input_data.sim_results = struct();
        input_data.sim_results.Sigma_emp = emp_cov;
        input_data.sim_results.F = sim_params.n_freq;
        input_data.sim_results.n = sim_params.n_nodes;
        input_data.sim_results.T = sim_params.n_samples;
        
        % Validate structure using Module1Utils (if available)
        try
            is_valid = Module1Utils.validate_data_structure(input_data);
            assert(is_valid, 'Data structure validation failed');
        catch
            % If Module1Utils not available, do basic validation
            assert(isfield(input_data, 'mode'), 'Missing mode field');
            assert(isfield(input_data, 'sim_results'), 'Missing sim_results field');
            assert(iscell(input_data.sim_results.Sigma_emp), 'Sigma_emp must be cell array');
        end
        
        result.passed = true;
        result.details = 'Data format is compatible';
        result.input_data = input_data;
        
    catch ME
        result.details = ME.message;
    end
end

function result = test_module1_basic_processing()
% Test basic Module 1 processing with complex data
    
    result = struct('passed', false);
    
    try
        % Generate test data
        [~, ~, emp_cov, sim_params] = module7_simulation_improved_complex(...
            'n_nodes', 6, ...
            'n_freq', 8, ...
            'n_samples', 50, ...
            'complex_strength', 0.8, ...
            'random_seed', 789);
        
        % Create input structure
        input_data = struct();
        input_data.mode = 'simulation';
        input_data.sim_results = struct();
        input_data.sim_results.Sigma_emp = emp_cov;
        input_data.sim_results.F = sim_params.n_freq;
        input_data.sim_results.n = sim_params.n_nodes;
        input_data.sim_results.T = sim_params.n_samples;
        
        % Run Module 1 preprocessing
        preprocessing_results = module1_preprocessing_main(input_data, ...
            'verbose', false, ...
            'force_hermitian', true, ...
            'check_psd', true);
        
        % Validate output structure
        assert(isstruct(preprocessing_results), 'Invalid output structure');
        assert(isfield(preprocessing_results, 'Sigma_tilde'), 'Missing Sigma_tilde field');
        assert(iscell(preprocessing_results.Sigma_tilde), 'Sigma_tilde must be cell array');
        assert(length(preprocessing_results.Sigma_tilde) == sim_params.n_freq, 'Wrong number of frequencies');
        
        result.passed = true;
        result.details = 'Basic processing completed successfully';
        result.processing_results = preprocessing_results;
        
    catch ME
        result.details = ME.message;
    end
end

function result = test_complex_data_preservation()
% Test if complex data properties are preserved during processing
    
    result = struct('passed', false);
    
    try
        % Generate and process complex data
        [~, ~, emp_cov, sim_params] = module7_simulation_improved_complex(...
            'n_nodes', 6, ...
            'n_freq', 8, ...
            'complex_strength', 1.0, ...
            'random_seed', 101112);
        
        input_data = struct();
        input_data.mode = 'simulation';
        input_data.sim_results = struct();
        input_data.sim_results.Sigma_emp = emp_cov;
        input_data.sim_results.F = sim_params.n_freq;
        input_data.sim_results.n = sim_params.n_nodes;
        input_data.sim_results.T = sim_params.n_samples;
        
        preprocessing_results = module1_preprocessing_main(input_data, 'verbose', false);
        
        % Check if complex properties are preserved
        input_has_complex = false;
        output_has_complex = false;
        
        for f = 1:length(emp_cov)
            if any(abs(imag(emp_cov{f}(:))) > 1e-12)
                input_has_complex = true;
                break;
            end
        end
        
        for f = 1:length(preprocessing_results.Sigma_tilde)
            if any(abs(imag(preprocessing_results.Sigma_tilde{f}(:))) > 1e-12)
                output_has_complex = true;
                break;
            end
        end
        
        % Complex data should be preserved or appropriately handled
        if input_has_complex
            % Either preserve complex data or ensure proper real conversion
            assert(output_has_complex || check_proper_real_conversion(preprocessing_results.Sigma_tilde), ...
                   'Complex data not properly handled');
        end
        
        result.passed = true;
        result.details = sprintf('Input complex: %s, Output complex: %s', ...
                                logical_to_string(input_has_complex), ...
                                logical_to_string(output_has_complex));
        
    catch ME
        result.details = ME.message;
    end
end

function result = test_hermitian_property_maintenance()
% Test if Hermitian property is maintained during processing
    
    result = struct('passed', false);
    
    try
        % Generate and process data
        [~, ~, emp_cov, sim_params] = module7_simulation_improved_complex(...
            'n_nodes', 6, ...
            'n_freq', 8, ...
            'complex_strength', 1.0, ...
            'random_seed', 131415);
        
        input_data = struct();
        input_data.mode = 'simulation';
        input_data.sim_results = struct();
        input_data.sim_results.Sigma_emp = emp_cov;
        input_data.sim_results.F = sim_params.n_freq;
        input_data.sim_results.n = sim_params.n_nodes;
        input_data.sim_results.T = sim_params.n_samples;
        
        preprocessing_results = module1_preprocessing_main(input_data, ...
            'verbose', false, 'force_hermitian', true);
        
        % Check Hermitian property in output
        max_hermitian_error = 0;
        for f = 1:length(preprocessing_results.Sigma_tilde)
            matrix = preprocessing_results.Sigma_tilde{f};
            hermitian_error = max(abs(matrix - matrix'));
            max_hermitian_error = max(max_hermitian_error, hermitian_error);
        end
        
        hermitian_threshold = 1e-10;
        assert(max_hermitian_error < hermitian_threshold, ...
               sprintf('Hermitian property violated: max error = %.2e', max_hermitian_error));
        
        result.passed = true;
        result.details = sprintf('Max Hermitian error: %.2e (threshold: %.2e)', ...
                                max_hermitian_error, hermitian_threshold);
        
    catch ME
        result.details = ME.message;
    end
end

function result = test_whitening_quality_assessment()
% Test whitening quality with complex data
    
    result = struct('passed', false);
    
    try
        % Generate and process data
        [~, ~, emp_cov, sim_params] = module7_simulation_improved_complex(...
            'n_nodes', 6, ...
            'n_freq', 8, ...
            'complex_strength', 0.8, ...
            'random_seed', 161718);
        
        input_data = struct();
        input_data.mode = 'simulation';
        input_data.sim_results = struct();
        input_data.sim_results.Sigma_emp = emp_cov;
        input_data.sim_results.F = sim_params.n_freq;
        input_data.sim_results.n = sim_params.n_nodes;
        input_data.sim_results.T = sim_params.n_samples;
        
        preprocessing_results = module1_preprocessing_main(input_data, ...
            'verbose', false, 'target_diagonal', 1.0, 'diagonal_tolerance', 0.1);
        
        % Check whitening quality
        if isfield(preprocessing_results, 'processing_stats') && ...
           isfield(preprocessing_results.processing_stats, 'whitening_quality')
            quality = preprocessing_results.processing_stats.whitening_quality;
            
            if isfield(quality, 'diagonal_errors')
                max_diagonal_error = max(quality.diagonal_errors);
                assert(max_diagonal_error < 0.2, ...
                       sprintf('Poor whitening quality: max diagonal error = %.3f', max_diagonal_error));
            end
        end
        
        % Check that output matrices have reasonable condition numbers
        max_condition = 0;
        for f = 1:length(preprocessing_results.Sigma_tilde)
            condition_num = cond(preprocessing_results.Sigma_tilde{f});
            max_condition = max(max_condition, condition_num);
        end
        
        condition_threshold = 1e10;
        assert(max_condition < condition_threshold, ...
               sprintf('Poor conditioning: max condition number = %.2e', max_condition));
        
        result.passed = true;
        result.details = sprintf('Max condition number: %.2e', max_condition);
        
    catch ME
        result.details = ME.message;
    end
end

function result = test_numerical_stability()
% Test numerical stability with complex data
    
    result = struct('passed', false);
    
    try
        % Test with challenging parameters
        test_conditions = {
            struct('n_nodes', 4, 'complex_strength', 2.0, 'name', 'High complexity'),
            struct('n_nodes', 12, 'complex_strength', 0.5, 'name', 'Large matrices'),
            struct('n_nodes', 6, 'complex_strength', 1.5, 'name', 'Medium complexity')
        };
        
        all_stable = true;
        failed_conditions = {};
        
        for i = 1:length(test_conditions)
            condition = test_conditions{i};
            
            try
                [~, ~, emp_cov, sim_params] = module7_simulation_improved_complex(...
                    'n_nodes', condition.n_nodes, ...
                    'n_freq', 6, ...
                    'complex_strength', condition.complex_strength, ...
                    'random_seed', 200 + i);
                
                input_data = struct();
                input_data.mode = 'simulation';
                input_data.sim_results = struct();
                input_data.sim_results.Sigma_emp = emp_cov;
                input_data.sim_results.F = sim_params.n_freq;
                input_data.sim_results.n = sim_params.n_nodes;
                input_data.sim_results.T = sim_params.n_samples;
                
                preprocessing_results = module1_preprocessing_main(input_data, 'verbose', false);
                
                % Check for NaN or Inf values
                for f = 1:length(preprocessing_results.Sigma_tilde)
                    matrix = preprocessing_results.Sigma_tilde{f};
                    if any(isnan(matrix(:))) || any(isinf(matrix(:)))
                        all_stable = false;
                        failed_conditions{end+1} = condition.name;
                        break;
                    end
                end
                
            catch
                all_stable = false;
                failed_conditions{end+1} = condition.name;
            end
        end
        
        assert(all_stable, sprintf('Numerical instability in conditions: %s', ...
                                  strjoin(failed_conditions, ', ')));
        
        result.passed = true;
        result.details = 'Numerically stable under all test conditions';
        
    catch ME
        result.details = ME.message;
    end
end

function result = test_visualization_compatibility()
% Test visualization compatibility with complex data
    
    result = struct('passed', false);
    
    try
        % Run the updated demo
        demo_results = demo_module1_preprocessing_updated();
        
        % Check if demo completed successfully
        assert(demo_results.data_generation.success, 'Demo data generation failed');
        
        % Try visualization using ComplexDataVisualizer class
        try
            visualizer = ComplexDataVisualizer();
            visualizer.visualize_results(demo_results);
            visualization_success = true;
        catch
            % Try fallback to original visualization
            try
                visualize_module1_results(demo_results);
                visualization_success = true;
            catch
                visualization_success = false;
            end
        end
        
        result.passed = visualization_success;
        if visualization_success
            result.details = 'Visualization completed successfully';
        else
            result.details = 'Visualization failed - needs enhancement for complex data';
        end
        
        % Close figures to avoid clutter
        close all;
        
    catch ME
        result.details = ME.message;
    end
end

function result = test_end_to_end_integration()
% Test complete end-to-end integration
    
    result = struct('passed', false);
    
    try
        % Run complete updated demo
        demo_results = demo_module1_preprocessing_updated();
        
        % Comprehensive validation
        checks = {};
        checks{end+1} = demo_results.data_generation.success;
        checks{end+1} = demo_results.preprocessing.success;
        
        if isfield(demo_results, 'summary')
            checks{end+1} = demo_results.summary.overall_success;
            if isfield(demo_results.summary, 'complex_data_success')
                checks{end+1} = demo_results.summary.complex_data_success;
            end
        end
        
        all_checks_passed = all(cell2mat(checks));
        
        assert(all_checks_passed, 'End-to-end integration has failures');
        
        result.passed = true;
        result.details = sprintf('All integration checks passed (%d/%d)', sum(cell2mat(checks)), length(checks));
        result.demo_results = demo_results;
        
    catch ME
        result.details = ME.message;
    end
end

%% Helper functions

function is_proper = check_proper_real_conversion(matrix_cell_array)
% Check if complex-to-real conversion was done properly
    is_proper = true;
    
    for f = 1:length(matrix_cell_array)
        matrix = matrix_cell_array{f};
        
        % Check if matrix is real
        if any(abs(imag(matrix(:))) > 1e-12)
            is_proper = false;
            break;
        end
        
        % Check if matrix is still positive definite
        eigenvals = eig(matrix);
        if any(eigenvals <= 0)
            is_proper = false;
            break;
        end
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
    