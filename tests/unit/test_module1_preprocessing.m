function test_results = test_module1_preprocessing()
% TEST_MODULE1_PREPROCESSING - Comprehensive test suite for Module 1
%
% This function runs a complete test suite for all Module 1 components,
% including unit tests, integration tests, and validation tests.
%
% Usage:
%   test_results = test_module1_preprocessing()
%
% Outputs:
%   test_results - Structure containing detailed test results and statistics
%
% File location: tests/test_module1_preprocessing.m

    fprintf('========================================\n');
    fprintf('Module 1 Preprocessing Test Suite\n');
    fprintf('========================================\n\n');
    
    % Initialize test results
    test_results = struct();
    test_results.timestamp = datestr(now);
    test_results.total_tests = 0;
    test_results.passed_tests = 0;
    test_results.failed_tests = 0;
    test_results.test_details = {};
    
    %% Unit Tests
    fprintf('=== Unit Tests ===\n');
    test_results.unit_tests = run_unit_tests();
    update_test_counts(test_results, test_results.unit_tests);
    
    %% Integration Tests
    fprintf('\n=== Integration Tests ===\n');
    test_results.integration_tests = run_integration_tests();
    update_test_counts(test_results, test_results.integration_tests);
    
    %% Validation Tests
    fprintf('\n=== Validation Tests ===\n');
    test_results.validation_tests = run_validation_tests();
    update_test_counts(test_results, test_results.validation_tests);
    
    %% Performance Tests
    fprintf('\n=== Performance Tests ===\n');
    test_results.performance_tests = run_performance_tests();
    update_test_counts(test_results, test_results.performance_tests);
    
    %% Edge Case Tests
    fprintf('\n=== Edge Case Tests ===\n');
    test_results.edge_case_tests = run_edge_case_tests();
    update_test_counts(test_results, test_results.edge_case_tests);
    
    %% Generate Test Report
    fprintf('\n=== Test Summary ===\n');
    generate_test_report(test_results);
    
    fprintf('\nTest suite completed!\n');
end

function unit_test_results = run_unit_tests()
% Run unit tests for individual components
    
    unit_test_results = struct();
    unit_test_results.total = 0;
    unit_test_results.passed = 0;
    unit_test_results.failed = 0;
    unit_test_results.details = {};
    
    %% Test Data Acquisition
    fprintf('Testing data_acquisition...\n');
    [passed, details] = test_data_acquisition();
    unit_test_results = record_test_result(unit_test_results, 'data_acquisition', passed, details);
    
    %% Test Diagonal Smoothing
    fprintf('Testing diagonal_smoothing...\n');
    [passed, details] = test_diagonal_smoothing();
    unit_test_results = record_test_result(unit_test_results, 'diagonal_smoothing', passed, details);
    
    %% Test Whitening Matrix Construction
    fprintf('Testing whitening_matrix_construction...\n');
    [passed, details] = test_whitening_matrix_construction();
    unit_test_results = record_test_result(unit_test_results, 'whitening_matrix_construction', passed, details);
    
    %% Test Covariance Whitening
    fprintf('Testing covariance_whitening...\n');
    [passed, details] = test_covariance_whitening();
    unit_test_results = record_test_result(unit_test_results, 'covariance_whitening', passed, details);
    
    fprintf('Unit tests completed: %d/%d passed\n', unit_test_results.passed, unit_test_results.total);
end

function [passed, details] = test_data_acquisition()
% Test data_acquisition function
    
    passed = true;
    details = {};
    
    try
        % Test 1: Simple simulation data
        sim_results = create_simple_test_data(10, 5, 20);
        Sigma_emp = data_acquisition('simulation', sim_results);
        
        % Validate output
        assert(iscell(Sigma_emp), 'Output should be cell array');
        assert(length(Sigma_emp) == 5, 'Should have 5 frequencies');
        assert(all(cellfun(@(x) size(x, 1) == 10, Sigma_emp)), 'All matrices should be 10x10');
        
        details{end+1} = 'PASS: Simple simulation data test';
        
        % Test 2: Input validation
        try
            data_acquisition('invalid_mode', []);
            passed = false;
            details{end+1} = 'FAIL: Should reject invalid mode';
        catch
            details{end+1} = 'PASS: Correctly rejects invalid mode';
        end
        
        % Test 3: Missing arguments
        try
            data_acquisition('simulation');
            passed = false;
            details{end+1} = 'FAIL: Should require second argument';
        catch
            details{end+1} = 'PASS: Correctly requires second argument';
        end
        
    catch ME
        passed = false;
        details{end+1} = sprintf('FAIL: %s', ME.message);
    end
end

function [passed, details] = test_diagonal_smoothing()
% Test diagonal_smoothing function
    
    passed = true;
    details = {};
    
    try
        % Create test data
        Sigma_emp = create_test_covariances(8, 10);
        
        % Test 1: Basic functionality
        [g_smooth, Sigma_loaded] = diagonal_smoothing(Sigma_emp);
        
        assert(iscell(g_smooth), 'g_smooth should be cell array');
        assert(length(g_smooth) == length(Sigma_emp), 'Should have same number of frequencies');
        assert(all(cellfun(@(x) length(x) == 8, g_smooth)), 'All power vectors should have length 8');
        assert(all(cellfun(@(x) all(x > 0), g_smooth)), 'All powers should be positive');
        
        details{end+1} = 'PASS: Basic functionality test';
        
        % Test 2: Different smoothing methods
        methods = {'moving_average', 'lowpass', 'spline'};
        for i = 1:length(methods)
            method = methods{i};
            [g_smooth_method, ~] = diagonal_smoothing(Sigma_emp, 'smoothing_method', method);
            assert(iscell(g_smooth_method), sprintf('%s method should return cell array', method));
            assert(all(cellfun(@(x) all(isreal(x)), g_smooth_method)), ...
                   sprintf('%s method should return real values', method));
        end
        details{end+1} = 'PASS: Different smoothing methods test';
        
        % Test 3: Diagonal loading
        [~, Sigma_loaded] = diagonal_smoothing(Sigma_emp, 'diagonal_loading', true, 'loading_factor', 0.1);
        assert(iscell(Sigma_loaded), 'Sigma_loaded should be cell array');
        
        % Check that loading improved condition numbers
        orig_cond = cellfun(@cond, Sigma_emp);
        loaded_cond = cellfun(@cond, Sigma_loaded);
        assert(mean(loaded_cond) <= mean(orig_cond) * 2, 'Loading should not drastically worsen conditioning');
        
        details{end+1} = 'PASS: Diagonal loading test';
        
        % Test 4: Parameter validation
        try
            diagonal_smoothing(Sigma_emp, 'window_size', -1);
            passed = false;
            details{end+1} = 'FAIL: Should reject negative window size';
        catch
            details{end+1} = 'PASS: Correctly rejects negative window size';
        end
        
    catch ME
        passed = false;
        details{end+1} = sprintf('FAIL: %s', ME.message);
    end
end

function [passed, details] = test_whitening_matrix_construction()
% Test whitening_matrix_construction function
    
    passed = true;
    details = {};
    
    try
        % Create test data
        n = 6;
        F = 8;
        g_smooth = cell(F, 1);
        for omega = 1:F
            g_smooth{omega} = 0.5 + 2 * rand(n, 1); % Random powers between 0.5 and 2.5
        end
        
        % Test 1: Basic functionality
        [D, stats] = whitening_matrix_construction(g_smooth);
        
        assert(iscell(D), 'D should be cell array');
        assert(length(D) == F, 'Should have same number of frequencies');
        
        for omega = 1:F
            assert(size(D{omega}, 1) == n && size(D{omega}, 2) == n, 'Each D should be nxn');
            assert(isreal(D{omega}), 'D should be real');
            
            % Check diagonal structure
            off_diag_norm = norm(D{omega} - diag(diag(D{omega})), 'fro');
            assert(off_diag_norm < 1e-12, 'D should be diagonal');
            
            % Check positive diagonal
            assert(all(diag(D{omega}) > 0), 'Diagonal elements should be positive');
            
            % Check relationship to input powers
            expected_diag = 1 ./ sqrt(g_smooth{omega});
            actual_diag = diag(D{omega});
            relative_error = max(abs(actual_diag - expected_diag) ./ expected_diag);
            assert(relative_error < 0.1, 'Diagonal should match 1/sqrt(g) within 10%');
        end
        
        details{end+1} = 'PASS: Basic functionality test';
        
        % Test 2: Regularization
        % Create data with very small powers
        g_bad = g_smooth;
        g_bad{1}(1) = 1e-15; % Very small power
        
        [D_reg, ~] = whitening_matrix_construction(g_bad, 'regularization', 'floor', 'min_power', 1e-10);
        assert(all(diag(D_reg{1}) < 1e5), 'Regularization should prevent extreme values');
        
        details{end+1} = 'PASS: Regularization test';
        
        % Test 3: Statistics structure
        assert(isstruct(stats), 'Should return statistics structure');
        required_fields = {'condition_numbers', 'power_ranges', 'regularized_count', 'final_power_ranges'};
        for i = 1:length(required_fields)
            assert(isfield(stats, required_fields{i}), sprintf('Missing field: %s', required_fields{i}));
        end
        
        details{end+1} = 'PASS: Statistics structure test';
        
    catch ME
        passed = false;
        details{end+1} = sprintf('FAIL: %s', ME.message);
    end
end

function [passed, details] = test_covariance_whitening()
% Test covariance_whitening function
    
    passed = true;
    details = {};
    
    try
        % Create test data
        n = 5;
        F = 6;
        Sigma_emp = create_test_covariances(n, F);
        
        % Create corresponding whitening matrices
        D = cell(F, 1);
        for omega = 1:F
            diag_powers = real(diag(Sigma_emp{omega}));
            D{omega} = diag(1 ./ sqrt(diag_powers));
        end
        
        % Test 1: Basic functionality
        [Sigma_tilde, quality] = covariance_whitening(Sigma_emp, D);
        
        assert(iscell(Sigma_tilde), 'Sigma_tilde should be cell array');
        assert(length(Sigma_tilde) == F, 'Should have same number of frequencies');
        
        for omega = 1:F
            assert(size(Sigma_tilde{omega}, 1) == n && size(Sigma_tilde{omega}, 2) == n, ...
                   'Each matrix should be nxn');
            
            % Check Hermitian property
            hermitian_error = norm(Sigma_tilde{omega} - Sigma_tilde{omega}', 'fro') / ...
                             norm(Sigma_tilde{omega}, 'fro');
            assert(hermitian_error < 1e-10, 'Should be Hermitian');
            
            % Check diagonal values are close to 1
            diag_vals = real(diag(Sigma_tilde{omega}));
            max_diag_error = max(abs(diag_vals - 1));
            assert(max_diag_error < 0.2, 'Diagonal values should be close to 1');
        end
        
        details{end+1} = 'PASS: Basic functionality test';
        
        % Test 2: Quality metrics
        assert(isstruct(quality), 'Should return quality structure');
        required_fields = {'diagonal_errors', 'hermitian_errors', 'min_eigenvalues', ...
                          'condition_numbers', 'whitening_effectiveness'};
        for i = 1:length(required_fields)
            assert(isfield(quality, required_fields{i}), sprintf('Missing field: %s', required_fields{i}));
        end
        
        assert(size(quality.diagonal_errors, 1) == F, 'Diagonal errors should have F rows');
        assert(size(quality.diagonal_errors, 2) == n, 'Diagonal errors should have n columns');
        assert(length(quality.whitening_effectiveness) == F, 'Effectiveness should have F elements');
        assert(all(quality.whitening_effectiveness >= 0 & quality.whitening_effectiveness <= 1), ...
               'Effectiveness should be between 0 and 1');
        
        details{end+1} = 'PASS: Quality metrics test';
        
        % Test 3: Input validation
        try
            covariance_whitening(Sigma_emp(1:end-1), D); % Mismatched lengths
            passed = false;
            details{end+1} = 'FAIL: Should reject mismatched input lengths';
        catch
            details{end+1} = 'PASS: Correctly rejects mismatched input lengths';
        end
        
    catch ME
        passed = false;
        details{end+1} = sprintf('FAIL: %s', ME.message);
    end
end

function integration_test_results = run_integration_tests()
% Run integration tests for complete pipeline
    
    integration_test_results = struct();
    integration_test_results.total = 0;
    integration_test_results.passed = 0;
    integration_test_results.failed = 0;
    integration_test_results.details = {};
    
    %% Test Complete Pipeline
    fprintf('Testing complete pipeline...\n');
    [passed, details] = test_complete_pipeline();
    integration_test_results = record_test_result(integration_test_results, 'complete_pipeline', passed, details);
    
    %% Test Pipeline Consistency
    fprintf('Testing pipeline consistency...\n');
    [passed, details] = test_pipeline_consistency();
    integration_test_results = record_test_result(integration_test_results, 'pipeline_consistency', passed, details);
    
    %% Test Parameter Propagation
    fprintf('Testing parameter propagation...\n');
    [passed, details] = test_parameter_propagation();
    integration_test_results = record_test_result(integration_test_results, 'parameter_propagation', passed, details);
    
    fprintf('Integration tests completed: %d/%d passed\n', ...
            integration_test_results.passed, integration_test_results.total);
end

function [passed, details] = test_complete_pipeline()
% Test complete preprocessing pipeline
    
    passed = true;
    details = {};
    
    try
        % Create test input
        input_data = struct();
        input_data.mode = 'simulation';
        input_data.sim_results = create_simple_test_data(12, 8, 25);
        
        % Run complete pipeline
        results = module1_preprocessing_main(input_data, 'verbose', false);
        
        % Validate output structure
        required_fields = {'Sigma_emp', 'g_smooth', 'D', 'Sigma_tilde', ...
                          'processing_stats', 'parameters', 'timing'};
        for i = 1:length(required_fields)
            assert(isfield(results, required_fields{i}), sprintf('Missing field: %s', required_fields{i}));
        end
        
        % Validate data consistency
        F = length(results.Sigma_emp);
        assert(length(results.g_smooth) == F, 'Inconsistent number of frequencies');
        assert(length(results.D) == F, 'Inconsistent number of frequencies');
        assert(length(results.Sigma_tilde) == F, 'Inconsistent number of frequencies');
        
        % Validate final output quality
        overall_stats = results.processing_stats.overall;
        assert(overall_stats.success.completed_all_steps, 'Should complete all steps');
        assert(overall_stats.success.all_matrices_valid, 'All matrices should be valid');
        
        details{end+1} = 'PASS: Complete pipeline test';
        
    catch ME
        passed = false;
        details{end+1} = sprintf('FAIL: %s', ME.message);
    end
end

function [passed, details] = test_pipeline_consistency()
% Test consistency of pipeline steps
    
    passed = true;
    details = {};
    
    try
        % Create test data
        input_data = struct();
        input_data.mode = 'simulation';
        input_data.sim_results = create_simple_test_data(10, 6, 30);
        
        % Run pipeline
        results = module1_preprocessing_main(input_data, 'verbose', false);
        
        F = length(results.Sigma_emp);
        n = size(results.Sigma_emp{1}, 1);
        
        % Test 1: Whitening transformation consistency
        for omega = 1:F
            % Manual whitening
            manual_whitened = results.D{omega} * results.Sigma_emp{omega} * results.D{omega};
            pipeline_whitened = results.Sigma_tilde{omega};
            
            relative_error = norm(manual_whitened - pipeline_whitened, 'fro') / ...
                            norm(pipeline_whitened, 'fro');
            assert(relative_error < 1e-10, sprintf('Whitening inconsistency at frequency %d', omega));
        end
        
        details{end+1} = 'PASS: Whitening transformation consistency';
        
        % Test 2: Power smoothing consistency
        for omega = 1:F
            g_vals = results.g_smooth{omega};
            D_vals = diag(results.D{omega});
            expected_D = 1 ./ sqrt(g_vals);
            
            relative_error = max(abs(D_vals - expected_D) ./ expected_D);
            assert(relative_error < 0.1, sprintf('Power-whitening inconsistency at frequency %d', omega));
        end
        
        details{end+1} = 'PASS: Power-whitening consistency';
        
    catch ME
        passed = false;
        details{end+1} = sprintf('FAIL: %s', ME.message);
    end
end

function [passed, details] = test_parameter_propagation()
% Test that parameters are properly propagated through pipeline
    
    passed = true;
    details = {};
    
    try
        % Create test data
        input_data = struct();
        input_data.mode = 'simulation';
        input_data.sim_results = create_simple_test_data(8, 5, 20);
        
        % Test specific parameters
        test_params = struct();
        test_params.smoothing_method = 'spline';
        test_params.window_size = 7;
        test_params.loading_factor = 0.05;
        test_params.target_diagonal = 1.5;
        
        results = module1_preprocessing_main(input_data, ...
            'smoothing_method', test_params.smoothing_method, ...
            'window_size', test_params.window_size, ...
            'loading_factor', test_params.loading_factor, ...
            'target_diagonal', test_params.target_diagonal, ...
            'verbose', false);
        
        % Check parameter storage
        stored_params = results.parameters;
        assert(strcmp(stored_params.smoothing_method, test_params.smoothing_method), ...
               'Smoothing method not stored correctly');
        assert(stored_params.window_size == test_params.window_size, ...
               'Window size not stored correctly');
        assert(stored_params.loading_factor == test_params.loading_factor, ...
               'Loading factor not stored correctly');
        assert(stored_params.target_diagonal == test_params.target_diagonal, ...
               'Target diagonal not stored correctly');
        
        details{end+1} = 'PASS: Parameter storage test';
        
        # Check parameter effects
        # Target diagonal should affect final whitened matrices
        for omega = 1:length(results.Sigma_tilde)
            diag_vals = real(diag(results.Sigma_tilde{omega}));
            mean_diag = mean(diag_vals);
            assert(abs(mean_diag - test_params.target_diagonal) < 0.3, ...
                   'Target diagonal parameter not effective');
        end
        
        details{end+1} = 'PASS: Parameter effectiveness test';
        
    catch ME
        passed = false;
        details{end+1} = sprintf('FAIL: %s', ME.message);
    end
end

function validation_test_results = run_validation_tests()
% Run validation tests with known ground truth
    
    validation_test_results = struct();
    validation_test_results.total = 0;
    validation_test_results.passed = 0;
    validation_test_results.failed = 0;
    validation_test_results.details = {};
    
    %% Test Known Ground Truth
    fprintf('Testing with known ground truth...\n');
    [passed, details] = test_known_ground_truth();
    validation_test_results = record_test_result(validation_test_results, 'known_ground_truth', passed, details);
    
    %% Test Mathematical Properties
    fprintf('Testing mathematical properties...\n');
    [passed, details] = test_mathematical_properties();
    validation_test_results = record_test_result(validation_test_results, 'mathematical_properties', passed, details);
    
    fprintf('Validation tests completed: %d/%d passed\n', ...
            validation_test_results.passed, validation_test_results.total);
end

function [passed, details] = test_known_ground_truth()
% Test with analytically known results
    
    passed = true;
    details = {};
    
    try
        % Create identity covariance matrices (already whitened)
        n = 8;
        F = 5;
        Sigma_emp = cell(F, 1);
        for omega = 1:F
            Sigma_emp{omega} = eye(n);
        end
        
        # Create sim_results structure
        sim_results = struct();
        sim_results.Sigma_emp = Sigma_emp;
        sim_results.F = F;
        sim_results.n = n;
        sim_results.T = 50;
        
        input_data = struct();
        input_data.mode = 'simulation';
        input_data.sim_results = sim_results;
        
        # Run preprocessing
        results = module1_preprocessing_main(input_data, 'verbose', false);
        
        # For identity matrices, whitening should preserve identity (approximately)
        for omega = 1:F
            whitened = results.Sigma_tilde{omega};
            identity_error = norm(whitened - eye(n), 'fro') / norm(eye(n), 'fro');
            assert(identity_error < 0.1, sprintf('Identity not preserved at frequency %d', omega));
        end
        
        details{end+1} = 'PASS: Identity preservation test';
        
        # Test with known scaling
        scale_factor = 4.0;
        for omega = 1:F
            Sigma_emp{omega} = scale_factor * eye(n);
        end
        sim_results.Sigma_emp = Sigma_emp;
        input_data.sim_results = sim_results;
        
        results2 = module1_preprocessing_main(input_data, 'verbose', false);
        
        # Whitening should remove the scaling
        for omega = 1:F
            whitened = results2.Sigma_tilde{omega};
            scaling_error = norm(whitened - eye(n), 'fro') / norm(eye(n), 'fro');
            assert(scaling_error < 0.1, sprintf('Scaling not removed at frequency %d', omega));
        end
        
        details{end+1} = 'PASS: Scaling removal test';
        
    catch ME
        passed = false;
        details{end+1} = sprintf('FAIL: %s', ME.message);
    end
end

function [passed, details] = test_mathematical_properties()
% Test fundamental mathematical properties
    
    passed = true;
    details = {};
    
    try
        # Create test data
        input_data = struct();
        input_data.mode = 'simulation';
        input_data.sim_results = create_simple_test_data(10, 6, 40);
        
        results = module1_preprocessing_main(input_data, 'verbose', false);
        
        F = length(results.Sigma_tilde);
        n = size(results.Sigma_tilde{1}, 1);
        
        # Test 1: Positive semi-definiteness
        for omega = 1:F
            eigenvals = eig(results.Sigma_tilde{omega});
            min_eigenval = min(real(eigenvals));
            assert(min_eigenval > -1e-10, sprintf('Not PSD at frequency %d (min eig: %.2e)', omega, min_eigenval));
        end
        
        details{end+1} = 'PASS: Positive semi-definiteness test';
        
        # Test 2: Hermitian symmetry
        for omega = 1:F
            matrix = results.Sigma_tilde{omega};
            hermitian_error = norm(matrix - matrix', 'fro') / norm(matrix, 'fro');
            assert(hermitian_error < 1e-12, sprintf('Not Hermitian at frequency %d', omega));
        end
        
        details{end+1} = 'PASS: Hermitian symmetry test';
        
        # Test 3: Determinant bounds
        for omega = 1:F
            det_val = det(results.Sigma_tilde{omega});
            assert(det_val > 0, sprintf('Non-positive determinant at frequency %d', omega));
            assert(isfinite(det_val), sprintf('Infinite determinant at frequency %d', omega));
        end
        
        details{end+1} = 'PASS: Determinant bounds test';
        
        # Test 4: Trace bounds (for whitened matrices, trace should be around n)
        for omega = 1:F
            trace_val = trace(results.Sigma_tilde{omega});
            expected_trace = n; % Approximately n for well-whitened matrices
            relative_error = abs(trace_val - expected_trace) / expected_trace;
            assert(relative_error < 0.5, sprintf('Unexpected trace at frequency %d', omega));
        end
        
        details{end+1} = 'PASS: Trace bounds test';
        
    catch ME
        passed = false;
        details{end+1} = sprintf('FAIL: %s', ME.message);
    end
end

function performance_test_results = run_performance_tests()
% Run performance and scalability tests
    
    performance_test_results = struct();
    performance_test_results.total = 0;
    performance_test_results.passed = 0;
    performance_test_results.failed = 0;
    performance_test_results.details = {};
    
    %% Test Scalability
    fprintf('Testing scalability...\n');
    [passed, details] = test_scalability();
    performance_test_results = record_test_result(performance_test_results, 'scalability', passed, details);
    
    %% Test Memory Usage
    fprintf('Testing memory usage...\n');
    [passed, details] = test_memory_usage();
    performance_test_results = record_test_result(performance_test_results, 'memory_usage', passed, details);
    
    fprintf('Performance tests completed: %d/%d passed\n', ...
            performance_test_results.passed, performance_test_results.total);
end

function [passed, details] = test_scalability()
% Test performance scaling with problem size
    
    passed = true;
    details = {};
    
    try
        # Test different problem sizes
        sizes = [10, 20, 30];
        times = zeros(size(sizes));
        
        for i = 1:length(sizes)
            n = sizes(i);
            F = 10;
            T = 50;
            
            input_data = struct();
            input_data.mode = 'simulation';
            input_data.sim_results = create_simple_test_data(n, F, T);
            
            tic;
            results = module1_preprocessing_main(input_data, 'verbose', false);
            times(i) = toc;
            
            # Verify successful completion
            assert(results.processing_stats.overall.success.overall, ...
                   sprintf('Failed to complete for size %d', n));
        end
        
        # Check that scaling is reasonable (should be roughly O(n^3) or better)
        if length(times) >= 2
            time_ratio = times(end) / times(1);
            size_ratio = sizes(end) / sizes(1);
            scaling_factor = time_ratio / (size_ratio^3);
            
            assert(scaling_factor < 5, 'Scaling worse than expected (>5x cubic)');
        end
        
        details{end+1} = sprintf('PASS: Scalability test (times: %s)', mat2str(times, 3));
        
    catch ME
        passed = false;
        details{end+1} = sprintf('FAIL: %s', ME.message);
    end
end

function [passed, details] = test_memory_usage()
% Test memory usage is reasonable
    
    passed = true;
    details = {};
    
    try
        # Test moderate size problem
        n = 25;
        F = 20;
        T = 100;
        
        input_data = struct();
        input_data.mode = 'simulation';
        input_data.sim_results = create_simple_test_data(n, F, T);
        
        # Monitor memory before and after
        mem_before = memory;
        results = module1_preprocessing_main(input_data, 'verbose', false);
        mem_after = memory;
        
        # Estimate expected memory usage
        bytes_per_complex = 16; # 8 bytes real + 8 bytes imaginary
        expected_memory = F * n * n * bytes_per_complex * 4; # Factor of 4 for intermediate storage
        expected_memory_mb = expected_memory / (1024^2);
        
        # Memory increase should be reasonable
        if isfield(mem_after, 'MemUsedMATLAB') && isfield(mem_before, 'MemUsedMATLAB')
            memory_increase = (mem_after.MemUsedMATLAB - mem_before.MemUsedMATLAB) / (1024^2);
            assert(memory_increase < expected_memory_mb * 10, 'Memory usage too high');
        end
        
        details{end+1} = sprintf('PASS: Memory usage test (expected ~%.1f MB)', expected_memory_mb);
        
    catch ME
        passed = false;
        details{end+1} = sprintf('FAIL: %s', ME.message);
    end
end

function edge_case_test_results = run_edge_case_tests()
% Run edge case and robustness tests
    
    edge_case_test_results = struct();
    edge_case_test_results.total = 0;
    edge_case_test_results.passed = 0;
    edge_case_test_results.failed = 0;
    edge_case_test_results.details = {};
    
    %% Test Extreme Parameters
    fprintf('Testing extreme parameters...\n');
    [passed, details] = test_extreme_parameters();
    edge_case_test_results = record_test_result(edge_case_test_results, 'extreme_parameters', passed, details);
    
    %% Test Degenerate Cases
    fprintf('Testing degenerate cases...\n');
    [passed, details] = test_degenerate_cases();
    edge_case_test_results = record_test_result(edge_case_test_results, 'degenerate_cases', passed, details);
    
    %% Test Numerical Stability
    fprintf('Testing numerical stability...\n');
    [passed, details] = test_numerical_stability();
    edge_case_test_results = record_test_result(edge_case_test_results, 'numerical_stability', passed, details);
    
    fprintf('Edge case tests completed: %d/%d passed\n', ...
            edge_case_test_results.passed, edge_case_test_results.total);
end

function [passed, details] = test_extreme_parameters()
% Test with extreme parameter values
    
    passed = true;
    details = {};
    
    try
        input_data = struct();
        input_data.mode = 'simulation';
        input_data.sim_results = create_simple_test_data(8, 5, 30);
        
        # Test 1: Very large loading factor
        results = module1_preprocessing_main(input_data, 'loading_factor', 1.0, 'verbose', false);
        assert(results.processing_stats.overall.success.completed_all_steps, ...
               'Should handle large loading factor');
        
        details{end+1} = 'PASS: Large loading factor test';
        
        # Test 2: Very small loading factor
        results = module1_preprocessing_main(input_data, 'loading_factor', 1e-10, 'verbose', false);
        assert(results.processing_stats.overall.success.completed_all_steps, ...
               'Should handle small loading factor');
        
        details{end+1} = 'PASS: Small loading factor test';
        
        # Test 3: Large window size
        results = module1_preprocessing_main(input_data, 'window_size', 99, 'verbose', false);
        assert(results.processing_stats.overall.success.completed_all_steps, ...
               'Should handle large window size');
        
        details{end+1} = 'PASS: Large window size test';
        
    catch ME
        passed = false;
        details{end+1} = sprintf('FAIL: %s', ME.message);
    end
end

function [passed, details] = test_degenerate_cases()
% Test degenerate input cases
    
    passed = true;
    details = {};
    
    try
        # Test 1: Single frequency
        input_data = struct();
        input_data.mode = 'simulation';
        input_data.sim_results = create_simple_test_data(5, 1, 20); # F = 1
        
        results = module1_preprocessing_main(input_data, 'verbose', false);
        assert(length(results.Sigma_tilde) == 1, 'Should handle single frequency');
        
        details{end+1} = 'PASS: Single frequency test';
        
        # Test 2: Small matrix size
        input_data.sim_results = create_simple_test_data(2, 3, 10); # n = 2
        
        results = module1_preprocessing_main(input_data, 'verbose', false);
        assert(size(results.Sigma_tilde{1}, 1) == 2, 'Should handle small matrices');
        
        details{end+1} = 'PASS: Small matrix test';
        
        # Test 3: Rank deficient matrices
        n = 5;
        F = 3;
        Sigma_emp = cell(F, 1);
        for omega = 1:F
            # Create rank-deficient matrix
            A = randn(n, 2); # Rank 2
            Sigma_emp{omega} = A * A' + 1e-6 * eye(n); # Small regularization
        end
        
        sim_results = struct();
        sim_results.Sigma_emp = Sigma_emp;
        sim_results.F = F;
        sim_results.n = n;
        sim_results.T = 20;
        
        input_data.sim_results = sim_results;
        
        results = module1_preprocessing_main(input_data, 'verbose', false);
        assert(results.processing_stats.overall.success.completed_all_steps, ...
               'Should handle rank-deficient matrices');
        
        details{end+1} = 'PASS: Rank-deficient matrices test';
        
    catch ME
        passed = false;
        details{end+1} = sprintf('FAIL: %s', ME.message);
    end
end

function [passed, details] = test_numerical_stability()
% Test numerical stability with challenging conditions
    
    passed = true;
    details = {};
    
    try
        % Test 1: Ill-conditioned matrices
        n = 6;
        F = 4;
        Sigma_emp = cell(F, 1);
        
        for omega = 1:F
            % Create ill-conditioned matrix
            [U, ~, V] = svd(randn(n));
            singular_vals = [1e6, 1e3, 1e0, 1e-3, 1e-6, 1e-9]; % Wide range
            Sigma_emp{omega} = U * diag(singular_vals) * V';
            Sigma_emp{omega} = (Sigma_emp{omega} + Sigma_emp{omega}') / 2; % Ensure Hermitian
        end
        
        sim_results = struct();
        sim_results.Sigma_emp = Sigma_emp;
        sim_results.F = F;
        sim_results.n = n;
        sim_results.T = 25;
        
        input_data = struct();
        input_data.mode = 'simulation';
        input_data.sim_results = sim_results;
        
        results = module1_preprocessing_main(input_data, 'verbose', false);
        
        % Check that all outputs are finite
        for omega = 1:F
            assert(all(isfinite(results.Sigma_tilde{omega}(:))), ...
                   sprintf('Non-finite values at frequency %d', omega));
        end
        
        details{end+1} = 'PASS: Ill-conditioned matrices test';
        
        % Test 2: Very large and very small values
        for omega = 1:F
            diagonal_vals = [1e10, 1e5, 1e0, 1e-5, 1e-10, 1e-12];
            Sigma_emp{omega} = diag(diagonal_vals);
        end
        
        input_data.sim_results.Sigma_emp = Sigma_emp;
        
        results = module1_preprocessing_main(input_data, 'verbose', false);
        
        % Check numerical stability
        for omega = 1:F
            matrix = results.Sigma_tilde{omega};
            assert(all(isfinite(matrix(:))), 'Should handle extreme values');
            assert(cond(matrix) < 1e12, 'Condition number should be controlled');
        end
        
        details{end+1} = 'PASS: Extreme values test';
        
    catch ME
        passed = false;
        details{end+1} = sprintf('FAIL: %s', ME.message);
    end
end

% Helper functions
function sim_results = create_simple_test_data(n, F, T)
% Create simple test data for testing
    
    Sigma_emp = cell(F, 1);
    
    for omega = 1:F
        % Create structured covariance matrix
        freq_factor = omega / F;
        
        % Base correlation with exponential decay
        base_corr = zeros(n, n);
        for i = 1:n
            for j = 1:n
                base_corr(i, j) = exp(-abs(i-j) / (n/4)) * 0.5;
            end
        end
        base_corr = base_corr + eye(n);
        
        % Frequency-dependent scaling
        scaling = 1 + 0.3 * sin(2 * pi * freq_factor);
        
        % Add some random variation
        diagonal_powers = scaling * (0.5 + rand(n, 1));
        
        Sigma_emp{omega} = diag(sqrt(diagonal_powers)) * base_corr * diag(sqrt(diagonal_powers));
        
        % Ensure positive definiteness
        [V, D] = eig(Sigma_emp{omega});
        D = diag(max(real(diag(D)), 1e-8));
        Sigma_emp{omega} = V * D * V';
    end
    
    sim_results = struct();
    sim_results.Sigma_emp = Sigma_emp;
    sim_results.F = F;
    sim_results.n = n;
    sim_results.T = T;
end

function Sigma_emp = create_test_covariances(n, F)
% Create test covariance matrices
    
    Sigma_emp = cell(F, 1);
    
    for omega = 1:F
        % Create random positive definite matrix
        A = randn(n, n) + 1i * randn(n, n);
        Sigma_emp{omega} = A * A' + 0.1 * eye(n);
        
        % Ensure Hermitian
        Sigma_emp{omega} = (Sigma_emp{omega} + Sigma_emp{omega}') / 2;
    end
end

function test_results = record_test_result(test_results, test_name, passed, details)
% Record test result in test_results structure
    
    test_results.total = test_results.total + 1;
    
    if passed
        test_results.passed = test_results.passed + 1;
        status = 'PASS';
    else
        test_results.failed = test_results.failed + 1;
        status = 'FAIL';
    end
    
    test_results.details{end+1} = struct('name', test_name, 'status', status, 'details', {details});
    
    fprintf('  %s: %s\n', test_name, status);
end

function update_test_counts(test_results, sub_results)
% Update overall test counts from sub-test results
    
    test_results.total_tests = test_results.total_tests + sub_results.total;
    test_results.passed_tests = test_results.passed_tests + sub_results.passed;
    test_results.failed_tests = test_results.failed_tests + sub_results.failed;
    
    % Accumulate detailed results
    for i = 1:length(sub_results.details)
        test_results.test_details{end+1} = sub_results.details{i};
    end
end

function generate_test_report(test_results)
% Generate comprehensive test report
    
    fprintf('========================================\n');
    fprintf('COMPREHENSIVE TEST REPORT\n');
    fprintf('========================================\n');
    fprintf('Total Tests: %d\n', test_results.total_tests);
    fprintf('Passed: %d (%.1f%%)\n', test_results.passed_tests, ...
            test_results.passed_tests / test_results.total_tests * 100);
    fprintf('Failed: %d (%.1f%%)\n', test_results.failed_tests, ...
            test_results.failed_tests / test_results.total_tests * 100);
    fprintf('========================================\n');
    
    % Test category breakdown
    categories = {'unit_tests', 'integration_tests', 'validation_tests', ...
                  'performance_tests', 'edge_case_tests'};
    category_names = {'Unit Tests', 'Integration Tests', 'Validation Tests', ...
                      'Performance Tests', 'Edge Case Tests'};
    
    fprintf('\nTest Category Breakdown:\n');
    fprintf('------------------------\n');
    
    for i = 1:length(categories)
        if isfield(test_results, categories{i})
            cat_results = test_results.(categories{i});
            fprintf('%s: %d/%d passed\n', category_names{i}, ...
                    cat_results.passed, cat_results.total);
        end
    end
    
    % Failed test details
    if test_results.failed_tests > 0
        fprintf('\nFailed Test Details:\n');
        fprintf('-------------------\n');
        
        for i = 1:length(test_results.test_details)
            test_detail = test_results.test_details{i};
            if strcmp(test_detail.status, 'FAIL')
                fprintf('%s:\n', test_detail.name);
                for j = 1:length(test_detail.details)
                    if contains(test_detail.details{j}, 'FAIL')
                        fprintf('  %s\n', test_detail.details{j});
                    end
                end
            end
        end
    end
    
    % Overall assessment
    fprintf('\nOverall Assessment:\n');
    fprintf('------------------\n');
    
    pass_rate = test_results.passed_tests / test_results.total_tests;
    if pass_rate == 1.0
        assessment = 'EXCELLENT - All tests passed!';
    elseif pass_rate >= 0.95
        assessment = 'VERY GOOD - Minor issues detected';
    elseif pass_rate >= 0.85
        assessment = 'GOOD - Some issues need attention';
    elseif pass_rate >= 0.70
        assessment = 'FAIR - Multiple issues detected';
    else
        assessment = 'POOR - Significant problems detected';
    end
    
    fprintf('%s\n', assessment);
    fprintf('Pass Rate: %.1f%%\n', pass_rate * 100);
    
    % Recommendations
    if test_results.failed_tests > 0
        fprintf('\nRecommendations:\n');
        fprintf('---------------\n');
        
        % Check specific failure patterns
        unit_failed = test_results.unit_tests.failed > 0;
        integration_failed = test_results.integration_tests.failed > 0;
        validation_failed = test_results.validation_tests.failed > 0;
        
        if unit_failed
            fprintf('- Review individual component implementations\n');
        end
        if integration_failed
            fprintf('- Check component interaction and data flow\n');
        end
        if validation_failed
            fprintf('- Verify mathematical correctness and algorithm implementation\n');
        end
        if test_results.performance_tests.failed > 0
            fprintf('- Optimize performance and memory usage\n');
        end
        if test_results.edge_case_tests.failed > 0
            fprintf('- Improve robustness and error handling\n');
        end
    end
    
    fprintf('\n========================================\n');
end