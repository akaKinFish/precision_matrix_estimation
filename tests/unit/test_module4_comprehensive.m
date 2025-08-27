function test_results = test_module4_comprehensive()
% TEST_MODULE4_COMPREHENSIVE - Comprehensive test suite for Module 4
%
% Syntax:
%   test_results = test_module4_comprehensive()
%
% Description:
%   Runs a complete test suite for Module 4 (Objective and Gradient Computation),
%   including mathematical verification, numerical stability tests, performance
%   benchmarks, and cross-method validation.
%   
%   Test Categories:
%   1. Mathematical Verification - Zero smoothing test, gradient symmetry
%   2. Numerical Differentiation - Finite difference gradient validation  
%   3. Method Comparison - Direct vs Laplacian implementation consistency
%   4. Boundary Conditions - Singular matrices, extreme parameters
%   5. Performance Analysis - Timing and memory usage assessment
%   6. Integration Tests - End-to-end workflow validation
%
% Output Arguments:
%   test_results - (struct) Contains:
%     .total_tests               - (integer) Number of tests run
%     .passed_tests              - (integer) Number of tests passed
%     .failed_tests              - (integer) Number of tests failed
%     .test_details              - (cell array) Detailed results per test
%     .performance_metrics       - (struct) Performance benchmarking results
%     .overall_success          - (logical) Whether all tests passed
%
% Examples:
%   % Run all tests
%   results = test_module4_comprehensive();
%   
%   % Check overall success
%   if results.overall_success
%       fprintf('All tests passed!\n');
%   else
%       fprintf('%d/%d tests failed\n', results.failed_tests, results.total_tests);
%   end
%   
%   % Review performance metrics
%   perf = results.performance_metrics;
%   fprintf('Average gradient computation: %.3fs\n', perf.avg_gradient_time);
%
% See also: MODULE4_OBJECTIVE_GRADIENT_MAIN, MODULE4_OBJECTIVE_EVALUATION
%
% Author: [Your Name]
% Date: [Current Date]
% Version: 1.0

fprintf('\n=== Comprehensive Module 4 Test Suite ===\n');
fprintf('Testing objective and gradient computation functionality\n\n');

% Initialize results structure
test_results = struct();
test_results.total_tests = 0;
test_results.passed_tests = 0;
test_results.failed_tests = 0;
test_results.test_details = {};
test_results.performance_metrics = struct();
test_results.overall_success = false;

% ==================== Test 1: Mathematical Verification ====================
fprintf('Test 1: Mathematical Verification Tests\n');
[passed, details] = test_mathematical_verification();
test_results = record_test_result(test_results, 'Mathematical Verification', passed, details);

% ==================== Test 2: Numerical Differentiation Validation ====================
fprintf('\nTest 2: Numerical Differentiation Validation\n');
[passed, details] = test_numerical_differentiation();
test_results = record_test_result(test_results, 'Numerical Differentiation', passed, details);

% ==================== Test 3: Method Comparison ====================
fprintf('\nTest 3: Direct vs Laplacian Method Comparison\n');
[passed, details] = test_method_comparison();
test_results = record_test_result(test_results, 'Method Comparison', passed, details);

% ==================== Test 4: Boundary Conditions ====================
fprintf('\nTest 4: Boundary Conditions and Edge Cases\n');
[passed, details] = test_boundary_conditions();
test_results = record_test_result(test_results, 'Boundary Conditions', passed, details);

% ==================== Test 5: Hermitian Properties ====================
fprintf('\nTest 5: Hermitian Symmetry Enforcement\n');
[passed, details] = test_hermitian_properties();
test_results = record_test_result(test_results, 'Hermitian Properties', passed, details);

% ==================== Test 6: Performance Benchmarking ====================
fprintf('\nTest 6: Performance Benchmarking\n');
[passed, details, perf_metrics] = test_performance_benchmarking();
test_results = record_test_result(test_results, 'Performance Benchmarking', passed, details);
test_results.performance_metrics = perf_metrics;

% ==================== Test 7: Integration Tests ====================
fprintf('\nTest 7: Integration Tests\n');
[passed, details] = test_integration_workflow();
test_results = record_test_result(test_results, 'Integration Workflow', passed, details);

% ==================== Final Results ====================
test_results.overall_success = (test_results.failed_tests == 0);

fprintf('\n============================================\n');
fprintf('Module 4 Test Suite Summary:\n');
fprintf('  Total tests: %d\n', test_results.total_tests);
fprintf('  Passed: %d\n', test_results.passed_tests);
fprintf('  Failed: %d\n', test_results.failed_tests);
fprintf('  Success rate: %.1f%%\n', 100 * test_results.passed_tests / test_results.total_tests);
if test_results.overall_success
    fprintf('  Overall success: YES\n');
else
    fprintf('  Overall success: NO\n');
end
fprintf('============================================\n');

if ~test_results.overall_success
    fprintf('\nFailed test details:\n');
    for i = 1:length(test_results.test_details)
        detail = test_results.test_details{i};
        if ~detail.passed
            if isfield(detail.details, 'error_message')
                fprintf('  - %s: %s\n', detail.test_name, detail.details.error_message);
            else
                fprintf('  - %s: Test failed\n', detail.test_name);
            end
        end
    end
end

end

% ==================== Individual Test Functions ====================

function [passed, details] = test_mathematical_verification()
% Test mathematical properties: zero smoothing, analytical solutions

passed = true;
details = struct();
details.subtests = {};

try
    % Create test problem
    [input_data, true_solution] = create_test_problem_small();
    
    % ========== Subtest 1: Zero Smoothing Test ==========
    fprintf('  1.1 Zero smoothing test (λ₁=λ₂=0)... ');
    
    % Set parameters for zero smoothing
    params = struct();
    params.lambda1 = 0;  % No smoothing
    params.lambda2 = 0;  % No L1 penalty
    params.verbose = false;
    params.force_hermitian = true;
    
    % Set precision matrices to analytical optimum: Γ = Σ^(-1)
    for f = 1:length(input_data.precision_matrices)
        input_data.precision_matrices{f} = inv(input_data.whitened_covariances{f});
    end
    
    % Compute gradients - should be near zero
    results = module4_objective_gradient_main(input_data, params);
    
    % Check gradient magnitudes
    max_gradient_norm = 0;
    for f = 1:length(results.smooth_gradients)
        grad_norm = norm(results.smooth_gradients{f}, 'fro');
        max_gradient_norm = max(max_gradient_norm, grad_norm);
    end
    
    tolerance = 1e-10;
    subtest1_passed = max_gradient_norm < tolerance;
    
    if subtest1_passed
        fprintf('PASSED (max gradient norm: %.2e)\n', max_gradient_norm);
    else
        fprintf('FAILED (max gradient norm: %.2e > %.2e)\n', max_gradient_norm, tolerance);
        passed = false;
    end
    
    details.subtests{end+1} = struct('name', 'Zero Smoothing', 'passed', subtest1_passed, ...
                                    'max_gradient_norm', max_gradient_norm);
    
    % ========== Subtest 2: Gradient Symmetry Test ==========
    fprintf('  1.2 Gradient Hermitian symmetry... ');
    
    % Use non-trivial parameters
    params.lambda1 = 0.05;
    params.lambda2 = 0.02;
    
    % Reset to test problem
    [input_data, ~] = create_test_problem_small();
    
    results = module4_objective_gradient_main(input_data, params);
    
    % Check Hermitian property of gradients
    max_hermitian_error = 0;
    for f = 1:length(results.smooth_gradients)
        G = results.smooth_gradients{f};
        hermitian_error = norm(G - G', 'fro');
        max_hermitian_error = max(max_hermitian_error, hermitian_error);
    end
    
    hermitian_tolerance = 1e-12;
    subtest2_passed = max_hermitian_error < hermitian_tolerance;
    
    if subtest2_passed
        fprintf('PASSED (max Hermitian error: %.2e)\n', max_hermitian_error);
    else
        fprintf('FAILED (max Hermitian error: %.2e > %.2e)\n', max_hermitian_error, hermitian_tolerance);
        passed = false;
    end
    
    details.subtests{end+1} = struct('name', 'Hermitian Symmetry', 'passed', subtest2_passed, ...
                                    'max_hermitian_error', max_hermitian_error);
    
    % ========== Subtest 3: Single Frequency Reduction ==========
    fprintf('  1.3 Single frequency reduction... ');
    
    % Create single-frequency problem
    single_freq_data = input_data;
    single_freq_data.precision_matrices = single_freq_data.precision_matrices(1);
    single_freq_data.whitened_covariances = single_freq_data.whitened_covariances(1);
    single_freq_data.smoothing_kernel = 0;  % 1x1 matrix with zero
    
    params_single = params;
    params_single.lambda1 = 0;  % No smoothing for single frequency
    
    results_single = module4_objective_gradient_main(single_freq_data, params_single);
    
    % Compare with expected result for standard graphical lasso
    expected_grad = -inv(single_freq_data.precision_matrices{1}) + single_freq_data.whitened_covariances{1};
    actual_grad = results_single.smooth_gradients{1};
    
    single_freq_error = norm(actual_grad - expected_grad, 'fro');
    single_freq_tolerance = 1e-10;
    subtest3_passed = single_freq_error < single_freq_tolerance;
    
    if subtest3_passed
        fprintf('PASSED (error: %.2e)\n', single_freq_error);
    else
        fprintf('FAILED (error: %.2e > %.2e)\n', single_freq_error, single_freq_tolerance);
        passed = false;
    end
    
    details.subtests{end+1} = struct('name', 'Single Frequency', 'passed', subtest3_passed, ...
                                    'reduction_error', single_freq_error);
    
catch ME
    passed = false;
    details.error_message = sprintf('Mathematical verification failed: %s', ME.message);
    fprintf('ERROR: %s\n', ME.message);
end

end

function [passed, details] = test_numerical_differentiation()
% Validate gradients using finite differences

passed = true;
details = struct();

try
    fprintf('  2.1 Finite difference validation... ');
    
    % Create small test problem for numerical differentiation
    [input_data, ~] = create_test_problem_tiny();  % Very small for finite diff
    
    params = struct();
    params.lambda1 = 0.03;
    params.lambda2 = 0;  % Only test smooth part
    params.verbose = false;
    
    % Compute analytical gradients
    results = module4_objective_gradient_main(input_data, params);
    analytical_gradients = results.smooth_gradients;
    
    % Finite difference parameters
    h = 1e-8;  % Step size
    max_relative_error = 0;
    
    % Test first frequency only (full test would be too expensive)
    f = 1;
    G_analytical = analytical_gradients{f};
    G_numerical = zeros(size(G_analytical));
    
    p = size(G_analytical, 1);
    
    % Compute numerical gradients using central differences
    for i = 1:min(p, 3)  % Limit to first few entries for speed
        for j = i:min(p, 3)  % Upper triangular only due to Hermitian
            
            % Generate Hermitian perturbation direction
            E_ij = zeros(p, p);
            if i == j
                E_ij(i, j) = 1;  % Real perturbation for diagonal
            else
                E_ij(i, j) = 0.5;
                E_ij(j, i) = 0.5;  % Hermitian perturbation for off-diagonal
            end
            
            % Forward step
            input_plus = input_data;
            input_plus.precision_matrices{f} = input_data.precision_matrices{f} + h * E_ij;
            [obj_plus, ~] = module4_objective_evaluation(input_plus, params);
            
            % Backward step
            input_minus = input_data;
            input_minus.precision_matrices{f} = input_data.precision_matrices{f} - h * E_ij;
            [obj_minus, ~] = module4_objective_evaluation(input_minus, params);
            
            % Central difference
            numerical_derivative = (obj_plus.smooth_objective - obj_minus.smooth_objective) / (2 * h);
            
            % Analytical directional derivative: <G, E_ij>
            analytical_derivative = real(trace(G_analytical' * E_ij));
            
            % Compare
            if abs(analytical_derivative) > 1e-12
                relative_error = abs(numerical_derivative - analytical_derivative) / abs(analytical_derivative);
                max_relative_error = max(max_relative_error, relative_error);
            end
        end
    end
    
    % Check tolerance
    finite_diff_tolerance = 1e-2;  % More reasonable for finite differences
    numerical_test_passed = max_relative_error < finite_diff_tolerance;
    
    if numerical_test_passed
        fprintf('PASSED (max relative error: %.2e)\n', max_relative_error);
    else
        fprintf('FAILED (max relative error: %.2e > %.2e)\n', max_relative_error, finite_diff_tolerance);
        passed = false;
    end
    
    details.max_relative_error = max_relative_error;
    details.tolerance = finite_diff_tolerance;
    
catch ME
    passed = false;
    details.error_message = sprintf('Numerical differentiation test failed: %s', ME.message);
    fprintf('ERROR: %s\n', ME.message);
end

end

function [passed, details] = test_method_comparison()
% Compare direct vs Laplacian methods

passed = true;
details = struct();

try
    fprintf('  3.1 Direct vs Laplacian consistency... ');
    
    [input_data, ~] = create_test_problem_medium();
    
    params_base = struct();
    params_base.lambda1 = 0.05;
    params_base.verbose = false;
    params_base.force_hermitian = true;
    
    % Compute with direct method
    params_direct = params_base;
    params_direct.use_graph_laplacian = false;
    results_direct = module4_objective_gradient_main(input_data, params_direct);
    
    % Compute with Laplacian method  
    params_laplacian = params_base;
    params_laplacian.use_graph_laplacian = true;
    results_laplacian = module4_objective_gradient_main(input_data, params_laplacian);
    
    % Compare results
    F = length(results_direct.smooth_gradients);
    max_difference = 0;
    
    for f = 1:F
        diff_f = norm(results_direct.smooth_gradients{f} - results_laplacian.smooth_gradients{f}, 'fro');
        max_difference = max(max_difference, diff_f);
    end
    
    % Check consistency
    consistency_tolerance = 1e-12;
    methods_consistent = max_difference < consistency_tolerance;
    
    if methods_consistent
        fprintf('PASSED (max difference: %.2e)\n', max_difference);
    else
        fprintf('FAILED (max difference: %.2e > %.2e)\n', max_difference, consistency_tolerance);
        passed = false;
    end
    
    details.max_difference = max_difference;
    details.direct_time = results_direct.computation_stats.total_computation_time;
    details.laplacian_time = results_laplacian.computation_stats.total_computation_time;
    details.speedup_ratio = details.direct_time / details.laplacian_time;
    
    fprintf('  3.2 Performance comparison - Direct: %.3fs, Laplacian: %.3fs (%.2fx speedup)\n', ...
            details.direct_time, details.laplacian_time, details.speedup_ratio);
    
catch ME
    passed = false;
    details.error_message = sprintf('Method comparison failed: %s', ME.message);
    fprintf('ERROR: %s\n', ME.message);
end

end

function [passed, details] = test_boundary_conditions()
% Test edge cases and boundary conditions

passed = true;
details = struct();
details.boundary_tests = {};

try
    % ========== Test 1: Near-singular matrices ==========
    fprintf('  4.1 Near-singular precision matrices... ');
    
    [input_data, ~] = create_test_problem_small();
    
    % Make one matrix nearly singular by reducing smallest eigenvalue
    [V, D] = eig(input_data.precision_matrices{1});
    eigenvals = diag(D);
    eigenvals(end) = 1e-12;  % Very small but not zero eigenvalue
    input_data.precision_matrices{1} = V * diag(eigenvals) * V';
    input_data.precision_matrices{1} = (input_data.precision_matrices{1} + input_data.precision_matrices{1}') / 2;
    
    params = struct('lambda1', 0.01, 'verbose', false, 'chol_tolerance', 1e-10);
    
    try
        results = module4_objective_gradient_main(input_data, params);
        
        % Check if gradients are finite
        all_finite = true;
        for f = 1:length(results.smooth_gradients)
            if any(~isfinite(results.smooth_gradients{f}(:)))
                all_finite = false;
                break;
            end
        end
        
        boundary_test1_passed = results.success && all_finite;
        
        if boundary_test1_passed
            fprintf('PASSED (handled with regularization)\n');
        else
            fprintf('FAILED (computation unsuccessful or non-finite results)\n');
            passed = false;
        end
    catch ME
        % Check if it's a controlled failure (acceptable)
        if contains(ME.message, 'positive definite') || contains(ME.message, 'singular')
            fprintf('PASSED (properly detected and handled singularity)\n');
            boundary_test1_passed = true;  % This is actually good behavior
        else
            fprintf('FAILED (unexpected exception: %s)\n', ME.message);
            passed = false;
            boundary_test1_passed = false;
        end
    end
    
    details.boundary_tests{end+1} = struct('name', 'Near-singular matrices', ...
                                          'passed', boundary_test1_passed);
    
    % ========== Test 2: Extreme parameters ==========
    fprintf('  4.2 Extreme parameter values... ');
    
    [input_data, ~] = create_test_problem_small();
    
    % Test very large lambda1
    params_extreme = struct('lambda1', 1e6, 'verbose', false);
    
    try
        results_extreme = module4_objective_gradient_main(input_data, params_extreme);
        
        % Check gradients are finite
        all_finite = true;
        for f = 1:length(results_extreme.smooth_gradients)
            if any(~isfinite(results_extreme.smooth_gradients{f}(:)))
                all_finite = false;
                break;
            end
        end
        
        boundary_test2_passed = results_extreme.success && all_finite;
        
        if boundary_test2_passed
            fprintf('PASSED\n');
        else
            fprintf('FAILED (non-finite gradients)\n');
            passed = false;
        end
    catch
        fprintf('FAILED (exception with extreme parameters)\n');
        passed = false;
        boundary_test2_passed = false;
    end
    
    details.boundary_tests{end+1} = struct('name', 'Extreme parameters', ...
                                          'passed', boundary_test2_passed);
    
catch ME
    passed = false;
    details.error_message = sprintf('Boundary condition tests failed: %s', ME.message);
    fprintf('ERROR: %s\n', ME.message);
end

end

function [passed, details] = test_hermitian_properties()
% Test Hermitian symmetry enforcement

passed = true;
details = struct();

try
    fprintf('  5.1 Hermitian symmetry enforcement... ');
    
    [input_data, ~] = create_test_problem_small();
    
    % Intentionally create non-Hermitian input (will be corrected internally)
    for f = 1:length(input_data.whitened_covariances)
        Sigma = input_data.whitened_covariances{f};
        % Add small non-Hermitian part
        non_hermitian_part = 1e-6 * randn(size(Sigma));
        non_hermitian_part = non_hermitian_part - non_hermitian_part';  % Make skew-Hermitian
        input_data.whitened_covariances{f} = Sigma + non_hermitian_part;
    end
    
    params = struct();
    params.lambda1 = 0.02;
    params.force_hermitian = true;
    params.symmetrization_tolerance = 1e-12;
    params.verbose = false;
    
    results = module4_objective_gradient_main(input_data, params);
    
    % Check all gradients are Hermitian
    max_hermitian_violation = max(results.hermitian_violations);
    hermitian_test_passed = max_hermitian_violation < 1e-10;
    
    if hermitian_test_passed
        fprintf('PASSED (max violation: %.2e)\n', max_hermitian_violation);
    else
        fprintf('FAILED (max violation: %.2e)\n', max_hermitian_violation);
        passed = false;
    end
    
    details.max_hermitian_violation = max_hermitian_violation;
    details.enforcement_count = results.computation_stats.hermitian_enforcement_count;
    
catch ME
    passed = false;
    details.error_message = sprintf('Hermitian property test failed: %s', ME.message);
    fprintf('ERROR: %s\n', ME.message);
end

end

function [passed, details, perf_metrics] = test_performance_benchmarking()
% Performance benchmarking tests

passed = true;
details = struct();
perf_metrics = struct();

try
    fprintf('  6.1 Performance benchmarking... ');
    
    % Test different problem sizes
    sizes = [10, 25, 50];
    frequencies = [5, 10, 20];
    
    perf_metrics.gradient_times = zeros(length(sizes), length(frequencies));
    perf_metrics.objective_times = zeros(length(sizes), length(frequencies));
    perf_metrics.memory_usage = zeros(length(sizes), length(frequencies));
    
    for i = 1:length(sizes)
        for j = 1:length(frequencies)
            p = sizes(i);
            F = frequencies(j);
            
            % Create test problem of specified size
            input_data = create_test_problem_size(p, F);
            params = struct('lambda1', 0.02, 'lambda2', 0.01, 'verbose', false);
            
            % Benchmark gradient computation
            tic_grad = tic;
            results_grad = module4_objective_gradient_main(input_data, params);
            time_grad = toc(tic_grad);
            
            % Benchmark objective evaluation
            tic_obj = tic;
            [results_obj, ~] = module4_objective_evaluation(input_data, params);
            time_obj = toc(tic_obj);
            
            perf_metrics.gradient_times(i, j) = time_grad;
            perf_metrics.objective_times(i, j) = time_obj;
            
            % Simple memory usage estimate (rough)
            mem_estimate = 8 * p^2 * F * 10;  % Bytes, rough estimate
            perf_metrics.memory_usage(i, j) = mem_estimate / 1024^2;  % MB
        end
    end
    
    % Compute averages
    perf_metrics.avg_gradient_time = mean(perf_metrics.gradient_times(:));
    perf_metrics.avg_objective_time = mean(perf_metrics.objective_times(:));
    perf_metrics.max_gradient_time = max(perf_metrics.gradient_times(:));
    perf_metrics.max_objective_time = max(perf_metrics.objective_times(:));
    
    % Simple performance check - should complete reasonably fast
    performance_acceptable = perf_metrics.max_gradient_time < 10.0;  % 10 seconds max
    
    if performance_acceptable
        fprintf('PASSED (avg gradient: %.3fs, max: %.3fs)\n', ...
                perf_metrics.avg_gradient_time, perf_metrics.max_gradient_time);
    else
        fprintf('WARNING (max gradient time: %.3fs > 10s)\n', perf_metrics.max_gradient_time);
        % Don't fail the test, just warn
    end
    
    details.performance_acceptable = performance_acceptable;
    
catch ME
    passed = false;
    details.error_message = sprintf('Performance benchmarking failed: %s', ME.message);
    perf_metrics = struct();
    fprintf('ERROR: %s\n', ME.message);
end

end

function [passed, details] = test_integration_workflow()
% Test complete integration workflow

passed = true;
details = struct();

try
    fprintf('  7.1 End-to-end integration workflow... ');
    
    % Create realistic test problem
    [input_data, ground_truth] = create_test_problem_medium();
    
    % Step 1: Compute initial objective
    params = struct('lambda1', 0.03, 'lambda2', 0.01, 'verbose', false);
    [initial_objective, ~] = module4_objective_evaluation(input_data, params);
    
    % Step 2: Compute gradients
    gradient_results = module4_objective_gradient_main(input_data, params);
    
    % Step 3: Test ObjectiveGradientComputer class
    computer = ObjectiveGradientComputer('lambda1', params.lambda1, 'lambda2', params.lambda2);
    computer.computation_options.verbose = false;
    
    % Test gradient computation through class
    class_gradients = computer.compute_smooth_gradients(input_data);
    
    % Test objective computation through class
    class_objective = computer.compute_full_objective(input_data);
    
    % Step 4: Compare results
    gradients_match = true;
    F = length(gradient_results.smooth_gradients);
    
    for f = 1:F
        diff_f = norm(gradient_results.smooth_gradients{f} - class_gradients{f}, 'fro');
        if diff_f > 1e-12
            gradients_match = false;
            break;
        end
    end
    
    objective_match = abs(initial_objective.total_objective - class_objective.total_objective) < 1e-12;
    
    % Step 5: Test method comparison through class
    comparison_results = computer.compare_gradient_methods(input_data);
    methods_consistent = comparison_results.methods_consistent;
    
    % Overall integration test
    integration_passed = gradient_results.success && gradients_match && ...
                        objective_match && methods_consistent;
    
    if integration_passed
        fprintf('PASSED\n');
    else
        fprintf('FAILED (gradients_match: %d, objective_match: %d, methods_consistent: %d)\n', ...
                gradients_match, objective_match, methods_consistent);
        passed = false;
    end
    
    details.gradients_match = gradients_match;
    details.objective_match = objective_match;
    details.methods_consistent = methods_consistent;
    details.comparison_results = comparison_results;
    
    % Additional integration checks
    fprintf('  7.2 Wrapper function compatibility... ');
    
    % Test wrapper function
    wrapper_results = module4_gradient(input_data, params);
    
    wrapper_compatible = isstruct(wrapper_results) && ...
                        isfield(wrapper_results, 'smooth_gradients') && ...
                        wrapper_results.success;
    
    if wrapper_compatible
        fprintf('PASSED\n');
    else
        fprintf('FAILED\n');
        passed = false;
    end
    
    details.wrapper_compatible = wrapper_compatible;
    
catch ME
    passed = false;
    details.error_message = sprintf('Integration workflow test failed: %s', ME.message);
    fprintf('ERROR: %s\n', ME.message);
end

end

% ==================== Helper Functions ====================

function test_results = record_test_result(test_results, test_name, passed, details)
% Record individual test result

test_results.total_tests = test_results.total_tests + 1;

if passed
    test_results.passed_tests = test_results.passed_tests + 1;
    fprintf('  ✓ %s: PASSED\n', test_name);
else
    test_results.failed_tests = test_results.failed_tests + 1;
    fprintf('  ✗ %s: FAILED\n', test_name);
end

% Store detailed results
detail_entry = struct();
detail_entry.test_name = test_name;
detail_entry.passed = passed;
detail_entry.details = details;

test_results.test_details{end+1} = detail_entry;

end

function [input_data, ground_truth] = create_test_problem_tiny()
% Create tiny test problem for numerical differentiation (3x3, 2 frequencies)

p = 3;  % Very small for finite differences
F = 2;

% Create ground truth sparse precision matrices
ground_truth = cell(F, 1);
for f = 1:F
    % Simple sparse structure
    Omega_true = eye(p) + 0.3 * (rand(p, p) - 0.5);
    Omega_true = (Omega_true + Omega_true') / 2;  % Make symmetric
    Omega_true = Omega_true + 0.5 * eye(p);      % Ensure positive definite
    ground_truth{f} = Omega_true;
end

% Generate empirical covariances
Sigmas = cell(F, 1);
for f = 1:F
    Sigmas{f} = inv(ground_truth{f}) + 0.01 * eye(p);  % Add noise
    Sigmas{f} = (Sigmas{f} + Sigmas{f}') / 2;  % Ensure Hermitian
end

% Initial precision estimates (slightly perturbed)
Gammas = cell(F, 1);
for f = 1:F
    Gammas{f} = ground_truth{f} + 0.05 * randn(p, p);
    Gammas{f} = (Gammas{f} + Gammas{f}') / 2;
    Gammas{f} = Gammas{f} + 1.0 * eye(p);  % Ensure positive definite
end

% Smoothing kernel (simple)
K = 0.5 * ones(F, F) - 0.5 * eye(F);  % Off-diagonal coupling

% Weight matrix (identity)
W = eye(p);

input_data = struct();
input_data.precision_matrices = Gammas;
input_data.whitened_covariances = Sigmas;
input_data.smoothing_kernel = K;
input_data.weight_matrix = W;

end

function [input_data, ground_truth] = create_test_problem_small()
% Create small test problem (5x5, 3 frequencies)

p = 5;
F = 3;

% Create ground truth precision matrices with sparse structure
ground_truth = cell(F, 1);
for f = 1:F
    % Start with identity
    Omega_true = eye(p);
    
    % Add some off-diagonal entries
    if f == 1
        Omega_true(1, 2) = 0.4; Omega_true(2, 1) = 0.4;
        Omega_true(2, 3) = -0.3; Omega_true(3, 2) = -0.3;
    elseif f == 2
        Omega_true(1, 2) = 0.3; Omega_true(2, 1) = 0.3;
        Omega_true(3, 4) = 0.2; Omega_true(4, 3) = 0.2;
    else
        Omega_true(2, 3) = -0.2; Omega_true(3, 2) = -0.2;
        Omega_true(4, 5) = 0.3; Omega_true(5, 4) = 0.3;
    end
    
    ground_truth{f} = Omega_true;
end

% Generate empirical covariances
Sigmas = cell(F, 1);
for f = 1:F
    Sigma_true = inv(ground_truth{f});
    noise = 0.02 * randn(p, p);
    noise = (noise + noise') / 2;
    Sigmas{f} = Sigma_true + noise;
    Sigmas{f} = (Sigmas{f} + Sigmas{f}') / 2;
    
    % Ensure positive definite
    [V, D] = eig(Sigmas{f});
    D = max(D, 0.01 * eye(p));
    Sigmas{f} = V * D * V';
end

% Initial precision estimates
Gammas = cell(F, 1);
for f = 1:F
    Gammas{f} = inv(Sigmas{f}) + 0.05 * randn(p, p);
    Gammas{f} = (Gammas{f} + Gammas{f}') / 2;
    
    % Ensure positive definite
    [V, D] = eig(Gammas{f});
    D = max(D, 0.1 * eye(p));
    Gammas{f} = V * D * V';
end

% Smoothing kernel (chain graph)
K = zeros(F, F);
for f = 1:F-1
    K(f, f+1) = 0.3;
    K(f+1, f) = 0.3;
end

% Weight matrix
W = eye(p) + 0.1 * ones(p, p);
W = W / norm(W);  % Normalize

input_data = struct();
input_data.precision_matrices = Gammas;
input_data.whitened_covariances = Sigmas;
input_data.smoothing_kernel = K;
input_data.weight_matrix = W;

end

function [input_data, ground_truth] = create_test_problem_medium()
% Create medium test problem (10x10, 5 frequencies)

p = 10;
F = 5;

% Create ground truth with block structure
ground_truth = cell(F, 1);
for f = 1:F
    Omega_true = 2 * eye(p);  % Strong diagonal
    
    % Add block structure
    block_size = 3;
    for block = 1:floor(p/block_size)
        start_idx = (block-1) * block_size + 1;
        end_idx = min(block * block_size, p);
        
        % Random connections within block
        block_strength = 0.2 + 0.1 * f / F;  % Vary across frequencies
        for i = start_idx:end_idx
            for j = (i+1):end_idx
                if rand() < 0.6  % 60% connection probability
                    val = block_strength * (rand() - 0.5) * 2;
                    Omega_true(i, j) = val;
                    Omega_true(j, i) = val;
                end
            end
        end
    end
    
    ground_truth{f} = Omega_true;
end

% Generate realistic empirical covariances
Sigmas = cell(F, 1);
for f = 1:F
    Sigma_true = inv(ground_truth{f});
    
    % Add sampling noise
    noise_level = 0.05;
    noise = noise_level * randn(p, p);
    noise = (noise + noise') / 2;
    
    Sigmas{f} = Sigma_true + noise;
    
    % Ensure positive definite with reasonable condition number
    [V, D] = eig(Sigmas{f});
    eigenvals = diag(D);
    eigenvals = max(eigenvals, 0.01);
    eigenvals = min(eigenvals, 100 * min(eigenvals));  % Limit condition number
    Sigmas{f} = V * diag(eigenvals) * V';
end

% Initial precision estimates (reasonable starting point)
Gammas = cell(F, 1);
for f = 1:F
    % Start with inverse of empirical covariance
    Gammas{f} = inv(Sigmas{f});
    
    % Add small perturbation
    perturbation = 0.02 * randn(p, p);
    perturbation = (perturbation + perturbation') / 2;
    Gammas{f} = Gammas{f} + perturbation;
    
    % Ensure positive definite
    [V, D] = eig(Gammas{f});
    eigenvals = diag(D);
    eigenvals = max(eigenvals, 0.05);
    Gammas{f} = V * diag(eigenvals) * V';
end

% Smoothing kernel (temporal chain with some long-range connections)
K = zeros(F, F);
for f = 1:F-1
    K(f, f+1) = 0.4;  % Temporal neighbors
    K(f+1, f) = 0.4;
end
for f = 1:F-2
    K(f, f+2) = 0.1;  % Long-range connections
    K(f+2, f) = 0.1;
end

% Weight matrix with some structure
W = eye(p);
for i = 1:p-1
    W(i, i+1) = 0.3;
    W(i+1, i) = 0.3;
end
W = W / max(eig(W)) * 0.8;  % Scale to be well-conditioned

input_data = struct();
input_data.precision_matrices = Gammas;
input_data.whitened_covariances = Sigmas;
input_data.smoothing_kernel = K;
input_data.weight_matrix = W;

end

function input_data = create_test_problem_size(p, F)
% Create test problem of specified size

% Simple but realistic structure
ground_truth = cell(F, 1);
for f = 1:F
    Omega_true = eye(p) + 0.1 * randn(p, p);
    Omega_true = (Omega_true + Omega_true') / 2;
    
    % Ensure positive definite with reasonable condition number
    [V, D] = eig(Omega_true);
    eigenvals = diag(D);
    eigenvals = max(eigenvals, 0.1);
    eigenvals = min(eigenvals, 10 * min(eigenvals));
    Omega_true = V * diag(eigenvals) * V';
    
    ground_truth{f} = Omega_true;
end

% Generate data
Sigmas = cell(F, 1);
Gammas = cell(F, 1);
for f = 1:F
    Sigma_true = inv(ground_truth{f});
    noise = 0.03 * randn(p, p);
    noise = (noise + noise') / 2;
    
    Sigmas{f} = Sigma_true + noise;
    [V, D] = eig(Sigmas{f});
    D = max(D, 0.01 * eye(p));
    Sigmas{f} = V * D * V';
    
    Gammas{f} = inv(Sigmas{f}) + 0.01 * randn(p, p);
    Gammas{f} = (Gammas{f} + Gammas{f}') / 2;
    [V, D] = eig(Gammas{f});
    D = max(D, 0.05 * eye(p));
    Gammas{f} = V * D * V';
end

% Simple smoothing kernel
K = 0.2 * (ones(F, F) - eye(F));

% Simple weight matrix
W = eye(p);

input_data = struct();
input_data.precision_matrices = Gammas;
input_data.whitened_covariances = Sigmas;
input_data.smoothing_kernel = K;
input_data.weight_matrix = W;

end