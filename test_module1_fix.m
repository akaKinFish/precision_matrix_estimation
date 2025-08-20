function test_result = test_module1_fix()
% TEST_MODULE1_FIX - Simple test to verify Module1 fixes are working
%
% This function tests the key repair points:
% 1. diagonal_smoothing function call fix
% 2. Basic complex data processing
% 3. End-to-end functionality verification
%
% Usage:
%   test_result = test_module1_fix()
%
% File location: test_module1_fix.m (project root directory)

    fprintf('========================================\n');
    fprintf('Module1 Fix Verification Test\n');
    fprintf('========================================\n\n');
    
    test_result = struct();
    test_result.timestamp = datestr(now);
    test_result.overall_success = false;
    test_result.details = {};
    
    %% Test 1: Check key function availability
    fprintf('Test 1: Check key function availability\n');
    fprintf('----------------------------------------\n');
    
    functions_to_check = {
        'module7_simulation_improved_complex',
        'module7_simulation',
        'module1_preprocessing_main',
        'diagonal_smoothing'
    };
    
    available_functions = {};
    for i = 1:length(functions_to_check)
        func_name = functions_to_check{i};
        if exist(func_name, 'file')
            fprintf('✅ %s available\n', func_name);
            available_functions{end+1} = func_name;
        else
            fprintf('❌ %s not available\n', func_name);
        end
    end
    
    test_result.available_functions = available_functions;
    
    %% Test 2: Test data generation
    fprintf('\nTest 2: Test data generation\n');
    fprintf('----------------------------\n');
    
    try
        % Try to use enhanced version
        if ismember('module7_simulation_improved_complex', available_functions)
            fprintf('Using enhanced complex simulation...\n');
            [~, ~, emp_cov, sim_params] = module7_simulation_improved_complex(...
                'n_nodes', 6, ...
                'n_freq', 8, ...
                'n_samples', 30, ...
                'complex_strength', 1.0, ...
                'random_seed', 123);
            test_result.data_generation_method = 'enhanced_complex';
        elseif ismember('module7_simulation', available_functions)
            fprintf('Using standard simulation...\n');
            [~, ~, emp_cov, sim_params] = module7_simulation(...
                'n_nodes', 6, ...
                'n_freq', 8, ...
                'n_samples', 30, ...
                'random_seed', 123);
            test_result.data_generation_method = 'standard';
        else
            error('No available data generation function');
        end
        
        fprintf('✅ Data generation successful\n');
        test_result.data_generation_success = true;
        test_result.details{end+1} = 'Data generation: SUCCESS';
        
        % Analyze generated data
        has_complex = false;
        for f = 1:length(emp_cov)
            if any(abs(imag(emp_cov{f}(:))) > 1e-12)
                has_complex = true;
                break;
            end
        end
        
        test_result.generated_complex_data = has_complex;
        if has_complex
            fprintf('ℹ️  Generated complex data\n');
        else
            fprintf('ℹ️  Generated real data\n');
        end
        
    catch ME
        fprintf('❌ Data generation failed: %s\n', ME.message);
        test_result.data_generation_success = false;
        test_result.data_generation_error = ME.message;
        test_result.details{end+1} = sprintf('Data generation: FAILED (%s)', ME.message);
        return;
    end
    
    %% Test 3: Test Module1 preprocessing (critical repair point)
    fprintf('\nTest 3: Test Module1 preprocessing\n');
    fprintf('----------------------------------\n');
    
    try
        % Create input data structure
        input_data = struct();
        input_data.mode = 'simulation';
        input_data.sim_results = struct();
        input_data.sim_results.Sigma_emp = emp_cov;
        input_data.sim_results.F = sim_params.n_freq;
        input_data.sim_results.n = sim_params.n_nodes;
        input_data.sim_results.T = sim_params.n_samples;
        
        % Call Module1 preprocessing (this will test the diagonal_smoothing fix)
        fprintf('Calling module1_preprocessing_main...\n');
        preprocessing_results = module1_preprocessing_main(input_data, ...
            'verbose', false, ...
            'force_hermitian', true, ...
            'check_psd', true);
        
        fprintf('✅ Module1 preprocessing successful\n');
        test_result.preprocessing_success = true;
        test_result.details{end+1} = 'Module1 preprocessing: SUCCESS';
        
        % Validate output structure
        required_outputs = {'Sigma_tilde', 'timing'};
        for i = 1:length(required_outputs)
            field = required_outputs{i};
            if isfield(preprocessing_results, field)
                fprintf('  ✅ Output contains %s\n', field);
            else
                fprintf('  ⚠️  Output missing %s\n', field);
            end
        end
        
        % Check output data quality
        if isfield(preprocessing_results, 'Sigma_tilde') && iscell(preprocessing_results.Sigma_tilde)
            output_matrices = preprocessing_results.Sigma_tilde;
            
            % Check matrix count
            if length(output_matrices) == sim_params.n_freq
                fprintf('  ✅ Output matrix count correct (%d)\n', length(output_matrices));
            else
                fprintf('  ❌ Output matrix count wrong: %d (expected %d)\n', ...
                        length(output_matrices), sim_params.n_freq);
            end
            
            % Check Hermitian property
            max_hermitian_error = 0;
            for f = 1:length(output_matrices)
                matrix = output_matrices{f};
                hermitian_error = max(abs(matrix - matrix'));
                max_hermitian_error = max(max_hermitian_error, hermitian_error);
            end
            
            if max_hermitian_error < 1e-10
                fprintf('  ✅ Hermitian property well preserved (error: %.2e)\n', max_hermitian_error);
            else
                fprintf('  ⚠️  Large Hermitian error: %.2e\n', max_hermitian_error);
            end
            
            test_result.hermitian_error = max_hermitian_error;
        end
        
    catch ME
        fprintf('❌ Module1 preprocessing failed: %s\n', ME.message);
        test_result.preprocessing_success = false;
        test_result.preprocessing_error = ME.message;
        test_result.details{end+1} = sprintf('Module1 preprocessing: FAILED (%s)', ME.message);
        
        % Check if this is the specific error we tried to fix
        if contains(ME.message, 'Too many output arguments')
            test_result.details{end+1} = '❌ diagonal_smoothing output argument error NOT FIXED';
            fprintf('  ❌ This is the diagonal_smoothing error we tried to fix\n');
            fprintf('  📝 Please apply the provided fix patch\n');
        end
        
        return;
    end
    
    %% Test 4: Test Demo execution
    fprintf('\nTest 4: Test Demo execution\n');
    fprintf('---------------------------\n');
    
    try
        demo_results = demo_module1_preprocessing();
        
        if demo_results.data_generation.success && demo_results.preprocessing.success
            fprintf('✅ Demo execution successful\n');
            test_result.demo_success = true;
            test_result.details{end+1} = 'Demo execution: SUCCESS';
        else
            fprintf('⚠️  Demo executed but with warnings\n');
            test_result.demo_success = false;
            test_result.details{end+1} = 'Demo execution: WARNING';
        end
        
    catch ME
        fprintf('❌ Demo execution failed: %s\n', ME.message);
        test_result.demo_success = false;
        test_result.demo_error = ME.message;
        test_result.details{end+1} = sprintf('Demo execution: FAILED (%s)', ME.message);
    end
    
    %% Test 5: Test visualization
    fprintf('\nTest 5: Test visualization\n');
    fprintf('--------------------------\n');
    
    try
        if exist('demo_results', 'var') && demo_results.preprocessing.success
            % Close figures to avoid interference
            close all;
            visualize_module1_results(demo_results);
            fprintf('✅ Visualization creation successful\n');
            test_result.visualization_success = true;
            test_result.details{end+1} = 'Visualization: SUCCESS';
            
            % Close created figures
            close all;
        else
            fprintf('⚠️  Skipping visualization test (Demo not successful)\n');
            test_result.visualization_success = false;
            test_result.details{end+1} = 'Visualization: SKIPPED';
        end
        
    catch ME
        fprintf('❌ Visualization failed: %s\n', ME.message);
        test_result.visualization_success = false;
        test_result.visualization_error = ME.message;
        test_result.details{end+1} = sprintf('Visualization: FAILED (%s)', ME.message);
    end
    
    %% Generate final results
    fprintf('\n========================================\n');
    fprintf('Test Results Summary\n');
    fprintf('========================================\n');
    
    % Calculate success rate
    successful_tests = 0;
    total_tests = 5;
    
    if test_result.data_generation_success
        successful_tests = successful_tests + 1;
    end
    if test_result.preprocessing_success
        successful_tests = successful_tests + 1;
    end
    if isfield(test_result, 'demo_success') && test_result.demo_success
        successful_tests = successful_tests + 1;
    end
    if isfield(test_result, 'visualization_success') && test_result.visualization_success
        successful_tests = successful_tests + 1;
    end
    
    test_result.success_rate = (successful_tests / total_tests) * 100;
    
    % Display detailed results
    for i = 1:length(test_result.details)
        fprintf('%s\n', test_result.details{i});
    end
    
    fprintf('\nSuccess rate: %.1f%% (%d/%d)\n', test_result.success_rate, successful_tests, total_tests);
    
    % Overall assessment
    if test_result.success_rate >= 80
        fprintf('🎉 Overall assessment: SUCCESS!\n');
        fprintf('   Module1 fix is effective, supports complex data processing\n');
        test_result.overall_success = true;
        test_result.recommendation = 'Fix successful, ready for normal use';
    elseif test_result.success_rate >= 60
        fprintf('⚠️  Overall assessment: PARTIAL SUCCESS\n');
        fprintf('   Core functionality working, but some features need refinement\n');
        test_result.overall_success = false;
        test_result.recommendation = 'Basic functionality working, recommend refining details';
    else
        fprintf('❌ Overall assessment: NEEDS REPAIR\n');
        fprintf('   Core issues not resolved, need to apply fix patches\n');
        test_result.overall_success = false;
        test_result.recommendation = 'Need to apply the provided fix patches';
    end
    
    fprintf('========================================\n');
    
end