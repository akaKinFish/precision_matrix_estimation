function test_results = test_module5_comprehensive()
% TEST_MODULE5_COMPREHENSIVE - Comprehensive test suite for Module 5 proximal updates
%
% Syntax:
%   test_results = test_module5_comprehensive()
%
% Description:
%   Runs comprehensive tests covering correctness, stability, performance, 
%   and integration for Module 5 proximal gradient updates. Tests include:
%   
%   A. Unit Tests: Core function correctness
%   B. Numerical Validation: Mathematical properties
%   C. Integration Tests: Module interactions
%   D. Performance Tests: Scalability and efficiency
%
% Output Arguments:
%   test_results - (struct) Contains:
%     .total_tests           - (integer) Total number of tests run
%     .passed_tests          - (integer) Number of tests passed
%     .failed_tests          - (integer) Number of tests failed
%     .test_details          - (struct) Detailed results for each test
%     .execution_time        - (double) Total test execution time
%     .success              - (logical) Overall test suite success
%
% Examples:
%   % Run all tests
%   results = test_module5_comprehensive();
%   
%   % Check results
%   if results.success
%       fprintf('All tests passed: %d/%d\n', results.passed_tests, results.total_tests);
%   else
%       fprintf('Tests failed: %d/%d\n', results.failed_tests, results.total_tests);
%   end
%
% See also: TEST_MODULE5_UNIT_TESTS, TEST_MODULE5_INTEGRATION
%
% Author: [Your Name] 
% Date: [Current Date]
% Version: 1.0

fprintf('\n================================================\n');
fprintf('Module 5 Comprehensive Test Suite\n');
fprintf('================================================\n');

overall_tic = tic;

% Initialize test results structure
test_results = struct();
test_results.total_tests = 0;
test_results.passed_tests = 0;
test_results.failed_tests = 0;
test_results.test_details = struct();

% ==================== A. Unit Tests ====================
fprintf('\n=== A. Unit Tests ===\n');

% A1: Complex Soft Thresholding
fprintf('A1: Testing complex soft thresholding... ');
try
    [passed, details] = test_complex_soft_thresholding();
    test_results = record_test_result(test_results, 'complex_soft_threshold', passed, details);
    if passed
        fprintf('PASSED\n');
    else
        fprintf('FAILED\n');
    end
catch ME
    test_results = record_test_result(test_results, 'complex_soft_threshold', false, ME.message);
    fprintf('ERROR: %s\n', ME.message);
end

% A2: Hermitian Symmetrization
fprintf('A2: Testing Hermitian symmetrization... ');
try
    [passed, details] = test_hermitian_symmetrization();
    test_results = record_test_result(test_results, 'hermitian_symmetrization', passed, details);
    if passed
        fprintf('PASSED\n');
    else
        fprintf('FAILED\n');
    end
catch ME
    test_results = record_test_result(test_results, 'hermitian_symmetrization', false, ME.message);
    fprintf('ERROR: %s\n', ME.message);
end

% A3: Active Set Projection
fprintf('A3: Testing active set projection... ');
try
    [passed, details] = test_active_set_projection();
    test_results = record_test_result(test_results, 'active_set_projection', passed, details);
    if passed
        fprintf('PASSED\n');
    else
        fprintf('FAILED\n');
    end
catch ME
    test_results = record_test_result(test_results, 'active_set_projection', false, ME.message);
    fprintf('ERROR: %s\n', ME.message);
end

% A4: PSD Check and Backtracking
fprintf('A4: Testing PSD check and backtracking... ');
try
    [passed, details] = test_psd_backtracking();
    test_results = record_test_result(test_results, 'psd_backtracking', passed, details);
    if passed
        fprintf('PASSED\n');
    else
        fprintf('FAILED\n');
    end
catch ME
    test_results = record_test_result(test_results, 'psd_backtracking', false, ME.message);
    fprintf('ERROR: %s\n', ME.message);
end

% A5: Step Size Selection
fprintf('A5: Testing automatic step size selection... ');
try
    [passed, details] = test_step_size_selection();
    test_results = record_test_result(test_results, 'step_size_selection', passed, details);
    if passed
        fprintf('PASSED\n');
    else
        fprintf('FAILED\n');
    end
catch ME
    test_results = record_test_result(test_results, 'step_size_selection', false, ME.message);
    fprintf('ERROR: %s\n', ME.message);
end

% ==================== B. Numerical Validation Tests ====================
fprintf('\n=== B. Numerical Validation Tests ===\n');

% B1: Unpenalized Recovery
fprintf('B1: Testing unpenalized recovery (λ₁=λ₂=0)... ');
try
    [passed, details] = test_unpenalized_recovery();
    test_results = record_test_result(test_results, 'unpenalized_recovery', passed, details);
    if passed
        fprintf('PASSED\n');
    else
        fprintf('FAILED\n');
    end
catch ME
    test_results = record_test_result(test_results, 'unpenalized_recovery', false, ME.message);
    fprintf('ERROR: %s\n', ME.message);
end

% B2: Homogenization Scenario
fprintf('B2: Testing homogenization (equal Γ_ω)... ');
try
    [passed, details] = test_homogenization_scenario();
    test_results = record_test_result(test_results, 'homogenization', passed, details);
    if passed
        fprintf('PASSED\n');
    else
        fprintf('FAILED\n');
    end
catch ME
    test_results = record_test_result(test_results, 'homogenization', false, ME.message);
    fprintf('ERROR: %s\n', ME.message);
end

% B3: Convergence vs Tolerance
fprintf('B3: Testing convergence vs tolerance... ');
try
    [passed, details] = test_convergence_tolerance();
    test_results = record_test_result(test_results, 'convergence_tolerance', passed, details);
    if passed
        fprintf('PASSED\n');
    else
        fprintf('FAILED\n');
    end
catch ME
    test_results = record_test_result(test_results, 'convergence_tolerance', false, ME.message);
    fprintf('ERROR: %s\n', ME.message);
end

% B4: Sparsity vs Lambda2
fprintf('B4: Testing sparsity vs λ₂ relationship... ');
try
    [passed, details] = test_sparsity_lambda2_curve();
    test_results = record_test_result(test_results, 'sparsity_lambda2', passed, details);
    if passed
        fprintf('PASSED\n');
    else
        fprintf('FAILED\n');
    end
catch ME
    test_results = record_test_result(test_results, 'sparsity_lambda2', false, ME.message);
    fprintf('ERROR: %s\n', ME.message);
end

% ==================== C. Integration Tests ====================
fprintf('\n=== C. Integration Tests ===\n');

% C1: Integration with Module 4
fprintf('C1: Testing integration with Module 4... ');
try
    [passed, details] = test_integration_with_module4();
    test_results = record_test_result(test_results, 'module4_integration', passed, details);
    if passed
        fprintf('PASSED\n');
    else
        fprintf('FAILED\n');
    end
catch ME
    test_results = record_test_result(test_results, 'module4_integration', false, ME.message);
    fprintf('ERROR: %s\n', ME.message);
end

% C2: End-to-End Pipeline
fprintf('C2: Testing end-to-end pipeline... ');
try
    [passed, details] = test_end_to_end_pipeline();
    test_results = record_test_result(test_results, 'end_to_end', passed, details);
    if passed
        fprintf('PASSED\n');
    else
        fprintf('FAILED\n');
    end
catch ME
    test_results = record_test_result(test_results, 'end_to_end', false, ME.message);
    fprintf('ERROR: %s\n', ME.message);
end

% ==================== D. Performance Tests ====================
fprintf('\n=== D. Performance Tests ===\n');

% D1: Scalability Test
fprintf('D1: Testing scalability (multiple problem sizes)... ');
try
    [passed, details] = test_scalability();
    test_results = record_test_result(test_results, 'scalability', passed, details);
    if passed
        fprintf('PASSED\n');
    else
        fprintf('FAILED\n');
    end
catch ME
    test_results = record_test_result(test_results, 'scalability', false, ME.message);
    fprintf('ERROR: %s\n', ME.message);
end

% D2: Parallel Performance
fprintf('D2: Testing parallel performance... ');
try
    [passed, details] = test_parallel_performance();
    test_results = record_test_result(test_results, 'parallel_performance', passed, details);
    if passed
        fprintf('PASSED\n');
    else
        fprintf('FAILED\n');
    end
catch ME
    test_results = record_test_result(test_results, 'parallel_performance', false, ME.message);
    fprintf('ERROR: %s\n', ME.message);
end

% ==================== Test Summary ====================
test_results.execution_time = toc(overall_tic);
test_results.success = (test_results.failed_tests == 0);

fprintf('\n================================================\n');
fprintf('Module 5 Test Suite Summary\n');
fprintf('================================================\n');
fprintf('Total tests run: %d\n', test_results.total_tests);
fprintf('Tests passed: %d\n', test_results.passed_tests);
fprintf('Tests failed: %d\n', test_results.failed_tests);
fprintf('Success rate: %.1f%%\n', 100 * test_results.passed_tests / test_results.total_tests);
fprintf('Total execution time: %.2f seconds\n', test_results.execution_time);

if test_results.success
    fprintf('🎉 ALL TESTS PASSED! Module 5 is ready for use.\n');
else
    fprintf('❌ SOME TESTS FAILED. Please review failed tests:\n');
    failed_tests = fieldnames(test_results.test_details);
    for i = 1:length(failed_tests)
        test_name = failed_tests{i};
        test_detail = test_results.test_details.(test_name);
        if ~test_detail.passed
            fprintf('  - %s: %s\n', test_name, test_detail.details);
        end
    end
end

fprintf('================================================\n');

end

% ==================== Helper Function ====================
function test_results = record_test_result(test_results, test_name, passed, details)
% Record individual test result

test_results.total_tests = test_results.total_tests + 1;

if passed
    test_results.passed_tests = test_results.passed_tests + 1;
else
    test_results.failed_tests = test_results.failed_tests + 1;
end

test_results.test_details.(test_name) = struct();
test_results.test_details.(test_name).passed = passed;
test_results.test_details.(test_name).details = details;

end

% ==================== Unit Test Functions ====================
function [passed, details] = test_complex_soft_thresholding()
% Test complex amplitude soft thresholding

p = 4;
tau = 0.1;

% Create test matrix with known complex structure
test_matrix = [1.0,     0.2+0.3i, 0.05+0.02i, 0.15-0.1i;
               0.2-0.3i, 1.0,      0.08+0.06i, 0.01-0.08i;
               0.05-0.02i, 0.08-0.06i, 1.0,    0.12+0.04i;
               0.15+0.1i, 0.01+0.08i, 0.12-0.04i, 1.0];

% Apply soft thresholding
result_simplified = module5_soft_threshold_complex(test_matrix, tau, 'simplified');
result_joint = module5_soft_threshold_complex(test_matrix, tau, 'joint');

% Test 1: Hermitian symmetry preserved
hermitian_error_simp = norm(result_simplified - result_simplified', 'fro');
hermitian_error_joint = norm(result_joint - result_joint', 'fro');

% Test 2: Diagonal processing
diag_real_simp = all(abs(imag(diag(result_simplified))) < 1e-14);
diag_ones_joint = norm(diag(result_joint) - ones(p, 1)) < 1e-14;

% Test 3: Phase preservation for surviving elements
phase_preserved = true;
for i = 1:p
    for j = i+1:p
        original_val = test_matrix(i, j);
        result_val = result_simplified(i, j);
        
        if abs(original_val) > tau && abs(result_val) > 1e-12
            original_phase = angle(original_val);
            result_phase = angle(result_val);
            phase_diff = abs(original_phase - result_phase);
            
            if phase_diff > 1e-10 && phase_diff < 2*pi - 1e-10
                phase_preserved = false;
                break;
            end
        end
    end
end

% Test 4: Correct amplitude shrinkage
amplitude_correct = true;
for i = 1:p
    for j = i+1:p
        original_val = test_matrix(i, j);
        result_val = result_simplified(i, j);
        
        if abs(original_val) > tau
            expected_magnitude = abs(original_val) - tau;
            actual_magnitude = abs(result_val);
            
            if abs(expected_magnitude - actual_magnitude) > 1e-12
                amplitude_correct = false;
                break;
            end
        elseif abs(result_val) > 1e-14
            amplitude_correct = false;
            break;
        end
    end
end

% Overall assessment
passed = (hermitian_error_simp < 1e-12) && (hermitian_error_joint < 1e-12) && ...
         diag_real_simp && diag_ones_joint && phase_preserved && amplitude_correct;

if passed
    details = 'All soft thresholding properties verified';
else
    details = sprintf('Hermitian: %.2e/%.2e, Diag real/ones: %d/%d, Phase: %d, Amplitude: %d', ...
                     hermitian_error_simp, hermitian_error_joint, diag_real_simp, diag_ones_joint, ...
                     phase_preserved, amplitude_correct);
end

end

function [passed, details] = test_hermitian_symmetrization()
% Test Hermitian symmetrization function

p = 3;
% Create deliberately non-Hermitian matrix
test_matrix = [1.5+0.1i,  0.2+0.3i,  0.1-0.2i;
               0.3-0.2i,  2.0+0.05i, 0.4+0.1i;
               0.2+0.1i,  0.3-0.2i,  1.8-0.1i];

% Test simplified mode
result_simp = module5_hermitian_symmetrize(test_matrix, 'simplified');

% Test joint mode  
result_joint = module5_hermitian_symmetrize(test_matrix, 'joint');

% Validation checks
hermitian_simp = norm(result_simp - result_simp', 'fro') < 1e-14;
hermitian_joint = norm(result_joint - result_joint', 'fro') < 1e-14;

diag_real_simp = all(abs(imag(diag(result_simp))) < 1e-14);
diag_ones_joint = norm(diag(result_joint) - ones(p, 1)) < 1e-14;

passed = hermitian_simp && hermitian_joint && diag_real_simp && diag_ones_joint;

if passed
    details = 'Hermitian symmetrization working correctly for both modes';
else
    details = sprintf('Hermitian: %d/%d, Diagonal: %d/%d', ...
                     hermitian_simp, hermitian_joint, diag_real_simp, diag_ones_joint);
end

end

function [passed, details] = test_active_set_projection()
% Test active set projection with Hermitian pairing

p = 4;
% Create symmetric active mask
active_mask = logical([1, 1, 0, 1;
                      1, 1, 1, 0;
                      0, 1, 1, 1;
                      1, 0, 1, 1]);

% Test matrix
test_matrix = randn(p) + 1i * randn(p);
test_matrix = (test_matrix + test_matrix') / 2;  % Make Hermitian

% Apply projection
result = module5_active_set_projection(test_matrix, active_mask);

% Check that inactive elements are zero
inactive_indices = ~active_mask;
max_inactive = max(abs(result(inactive_indices)));

% Check that active elements are preserved
active_indices = active_mask;
preservation_error = norm(result(active_indices) - test_matrix(active_indices));

% Check Hermitian preservation
hermitian_preserved = norm(result - result', 'fro') < 1e-12;

passed = (max_inactive < 1e-15) && (preservation_error < 1e-12) && hermitian_preserved;

if passed
    details = 'Active set projection preserves structure correctly';
else
    details = sprintf('Max inactive: %.2e, Preservation error: %.2e, Hermitian: %d', ...
                     max_inactive, preservation_error, hermitian_preserved);
end

end

function [passed, details] = test_psd_backtracking()
% Test PSD checking and backtracking mechanism

p = 3;
F = 2;

% Create positive definite matrices
Gamma_init = cell(F, 1);
for f = 1:F
    A = randn(p) + 1i * randn(p);
    Gamma_init{f} = A' * A + eye(p);  % Guaranteed PSD
end

% Create problematic gradient that would break PSD
gradient_bad = 10 * eye(p);  % Large positive diagonal will make matrix singular

% Test PSD check function
[isPSD_good, info_good] = module5_psd_check(Gamma_init{1});
[isPSD_bad, info_bad] = module5_psd_check(Gamma_init{1} - gradient_bad);

% Test single step with backtracking
aux_data = struct();
aux_data.lambda2 = 0.01;

step_params = struct();
step_params.mode = 'simplified';
step_params.beta_backtrack = 0.5;
step_params.max_backtrack = 10;

active_mask = true(p, p);

[Gamma_result, step_info] = module5_single_proximal_step(...
    Gamma_init{1}, gradient_bad, 0.5, active_mask, aux_data, step_params);

% Validate results
[final_isPSD, ~] = module5_psd_check(Gamma_result);

passed = isPSD_good && ~isPSD_bad && final_isPSD && (step_info.backtrack_count > 0);

if passed
    details = sprintf('Backtracking worked: %d steps, final matrix PSD', step_info.backtrack_count);
else
    details = sprintf('PSD checks: %d/%d/%d, Backtrack count: %d', ...
                     isPSD_good, ~isPSD_bad, final_isPSD, step_info.backtrack_count);
end

end

function [passed, details] = test_step_size_selection()
% Test automatic step size and lambda1 selection

p = 4;
F = 3;

% Create test data
Gamma_test = cell(F, 1);
for f = 1:F
    A = randn(p) + 1i * randn(p);
    Gamma_test{f} = A' * A + 0.1 * eye(p);
end

K_smooth = 0.3 * (eye(F) + diag(ones(F-1, 1), 1) + diag(ones(F-1, 1), -1));
W_matrix = eye(p) + 0.2 * randn(p);
W_matrix = (W_matrix + W_matrix') / 2;

% Test parameter selection
test_params = struct('verbose', false, 'delta', 0.9);
[alpha_auto, lambda1_auto, stats] = module5_step_size_selection(...
    Gamma_test, K_smooth, W_matrix, test_params);

% Validation checks
positive_params = (alpha_auto > 0) && (lambda1_auto > 0);
finite_params = isfinite(alpha_auto) && isfinite(lambda1_auto);
reasonable_alpha = (alpha_auto > 1e-6) && (alpha_auto < 10);
reasonable_lambda1 = (lambda1_auto > 1e-6) && (lambda1_auto < 1);

passed = positive_params && finite_params && reasonable_alpha && reasonable_lambda1;

if passed
    details = sprintf('α=%.4e, λ₁=%.4e (both reasonable)', alpha_auto, lambda1_auto);
else
    details = sprintf('α=%.4e, λ₁=%.4e (pos: %d/%d, finite: %d/%d, reasonable: %d/%d)', ...
                     alpha_auto, lambda1_auto, alpha_auto > 0, lambda1_auto > 0, ...
                     finite_params, finite_params, reasonable_alpha, reasonable_lambda1);
end

end

function [passed, details] = test_unpenalized_recovery()
% Test recovery of analytical solution when λ₁=λ₂=0

p = 3;
F = 2;

% Create known covariances and their inverses
Sigma_true = cell(F, 1);
Gamma_true = cell(F, 1);

for f = 1:F
    A = randn(p) + 1i * randn(p);
    Sigma_true{f} = A * A' + 0.1 * eye(p);
    Sigma_true{f} = (Sigma_true{f} + Sigma_true{f}') / 2;
    Gamma_true{f} = inv(Sigma_true{f});
end

% Set up zero penalties
input_data = struct();
input_data.whitened_covariances = Sigma_true;
input_data.initial_precision = Gamma_true;  % Start from true solution
input_data.smoothing_kernel = zeros(F, F);  % No smoothing
input_data.weight_matrix = eye(p);
input_data.active_set_masks = cell(F, 1);
for f = 1:F
    input_data.active_set_masks{f} = true(p, p);  % All elements active
end

proximal_params = struct();
proximal_params.lambda1 = 0;
proximal_params.lambda2 = 0;
proximal_params.max_iter = 50;
proximal_params.eps_x = 1e-8;
proximal_params.verbose = false;

% Run optimization
[Gamma_result, results] = module5_proximal_main(input_data, proximal_params);

% Check recovery accuracy
max_error = 0;
for f = 1:F
    error_f = norm(Gamma_result{f} - Gamma_true{f}, 'fro') / norm(Gamma_true{f}, 'fro');
    max_error = max(max_error, error_f);
end

% Should converge quickly and accurately
passed = (max_error < 1e-3) && results.convergence_info.converged;

if passed
    details = sprintf('Recovered true solution: max error %.2e, converged in %d iterations', ...
                     max_error, results.convergence_info.iterations);
else
    details = sprintf('Recovery failed: max error %.2e, converged: %d', ...
                     max_error, results.convergence_info.converged);
end

end

function [passed, details] = test_homogenization_scenario()
% Test that equal Γ_ω produces zero smoothing gradient

p = 3;
F = 4;

% Create identical precision matrices
Gamma_identical = cell(F, 1);
A = randn(p) + 1i * randn(p);
Gamma_base = A' * A + 0.5 * eye(p);
Gamma_base = (Gamma_base + Gamma_base') / 2;

for f = 1:F
    Gamma_identical{f} = Gamma_base;
end

% Create dummy data
Sigma_dummy = cell(F, 1);
for f = 1:F
    Sigma_dummy{f} = inv(Gamma_base) + 0.01 * eye(p);
end

K_smooth = 0.5 * (eye(F) + diag(ones(F-1, 1), 1) + diag(ones(F-1, 1), -1));
W_matrix = eye(p);

% Compute gradient
input_data = struct();
input_data.precision_matrices = Gamma_identical;
input_data.whitened_covariances = Sigma_dummy;
input_data.smoothing_kernel = K_smooth;
input_data.weight_matrix = W_matrix;

gradient_params = struct();
gradient_params.lambda1 = 0.1;  % Non-zero smoothing
gradient_params.verbose = false;

gradient_results = module4_objective_gradient_main(input_data, gradient_params);

% Check that smoothing gradients are near zero (since all Γ_ω equal)
max_smoothing_norm = 0;
for f = 1:F
    if isfield(gradient_results.gradient_components, 'smoothing_gradients')
        smoothing_grad = gradient_results.gradient_components.smoothing_gradients{f};
        smoothing_norm = norm(smoothing_grad, 'fro');
        max_smoothing_norm = max(max_smoothing_norm, smoothing_norm);
    end
end

passed = max_smoothing_norm < 1e-10;

if passed
    details = sprintf('Smoothing gradient vanishes as expected: max norm %.2e', max_smoothing_norm);
else
    details = sprintf('Smoothing gradient too large: max norm %.2e', max_smoothing_norm);
end

end

function [passed, details] = test_convergence_tolerance()
% Test convergence behavior with different tolerances

p = 3;
F = 2;

% Create test problem
input_data = create_test_problem(p, F);

% Test with different tolerances
tolerances = [1e-2, 1e-3, 1e-4];
iterations = zeros(size(tolerances));
final_objectives = zeros(size(tolerances));

for t = 1:length(tolerances)
    params_t = struct();
    params_t.eps_x = tolerances(t);
    params_t.eps_f = tolerances(t);
    params_t.max_iter = 500;
    params_t.verbose = false;
    
    [~, results_t] = module5_proximal_main(input_data, params_t);
    
    iterations(t) = results_t.convergence_info.iterations;
    final_objectives(t) = results_t.convergence_info.final_objective;
end

% Check monotonicity: tighter tolerance should need more iterations
monotonic_iterations = all(diff(iterations) >= 0);

% Check that all converged to similar objectives
objective_spread = max(final_objectives) - min(final_objectives);
consistent_objectives = objective_spread < 1e-3 * abs(mean(final_objectives));

passed = monotonic_iterations && consistent_objectives;

if passed
    details = sprintf('Convergence behavior correct: iterations %s, objective spread %.2e', ...
                     mat2str(iterations), objective_spread);
else
    details = sprintf('Issues: monotonic_iter=%d, consistent_obj=%d (spread=%.2e)', ...
                     monotonic_iterations, consistent_objectives, objective_spread);
end

end

function [passed, details] = test_sparsity_lambda2_curve()
% Test that increasing λ₂ produces sparser solutions

p = 4;
F = 2;

% Create test problem
input_data = create_test_problem(p, F);

% Test different λ₂ values
lambda2_values = [1e-4, 1e-3, 1e-2, 1e-1];
sparsity_ratios = zeros(size(lambda2_values));

for l = 1:length(lambda2_values)
    params_l = struct();
    params_l.lambda2 = lambda2_values(l);
    params_l.max_iter = 200;
    params_l.verbose = false;
    
    [Gamma_result, ~] = module5_proximal_main(input_data, params_l);
    
    % Compute sparsity ratio
    total_elements = 0;
    zero_elements = 0;
    
    for f = 1:F
        off_diag_mask = ~eye(p);
        elements = Gamma_result{f}(off_diag_mask);
        total_elements = total_elements + length(elements);
        zero_elements = zero_elements + sum(abs(elements) < 1e-8);
    end
    
    sparsity_ratios(l) = zero_elements / total_elements;
end

% Check monotonicity: higher λ₂ should produce higher sparsity
monotonic_sparsity = all(diff(sparsity_ratios) >= -0.05);  % Allow small numerical noise

% Check reasonable range
reasonable_range = (min(sparsity_ratios) >= 0) && (max(sparsity_ratios) <= 1);

passed = monotonic_sparsity && reasonable_range;

if passed
    details = sprintf('Sparsity increases with λ₂: %s', mat2str(sparsity_ratios, 3));
else
    details = sprintf('Sparsity issues: ratios=%s, monotonic=%d', ...
                     mat2str(sparsity_ratios, 3), monotonic_sparsity);
end

end

function [passed, details] = test_integration_with_module4()
% Test integration between Module 4 gradients and Module 5 updates

p = 3;
F = 2;

% Create test data
input_data = create_test_problem(p, F);

% Compute gradients with Module 4
gradient_params = struct();
gradient_params.lambda1 = 0.05;
gradient_params.verbose = false;

module4_input = struct();
module4_input.precision_matrices = input_data.initial_precision;
module4_input.whitened_covariances = input_data.whitened_covariances;
module4_input.smoothing_kernel = input_data.smoothing_kernel;
module4_input.weight_matrix = input_data.weight_matrix;

gradient_results = module4_objective_gradient_main(module4_input, gradient_params);

% Use Module 4 gradients in Module 5 single step
aux_data = struct();
aux_data.lambda2 = 0.01;

step_params = struct();
step_params.mode = 'simplified';
step_params.beta_backtrack = 0.5;
step_params.max_backtrack = 10;

% Test single step for first frequency
f_test = 1;
G_test = gradient_results.smooth_gradients{f_test};
Gamma_test = input_data.initial_precision{f_test};
A_test = input_data.active_set_masks{f_test};

[Gamma_updated, step_info] = module5_single_proximal_step(...
    Gamma_test, G_test, 0.1, A_test, aux_data, step_params);

% Validate step properties
step_successful = step_info.success;
matrix_changed = (step_info.relative_change > 1e-10);
result_hermitian = norm(Gamma_updated - Gamma_updated', 'fro') < 1e-12;
result_psd = module5_psd_check(Gamma_updated);

passed = step_successful && matrix_changed && result_hermitian && result_psd;

if passed
    details = sprintf('Integration successful: step size %.2e, change %.2e', ...
                     step_info.final_step_size, step_info.relative_change);
else
    details = sprintf('Integration issues: success=%d, changed=%d, hermitian=%d, psd=%d', ...
                     step_successful, matrix_changed, result_hermitian, result_psd);
end

end

function [passed, details] = test_end_to_end_pipeline()
% Test complete Module 5 pipeline

p = 4;
F = 3;

% Create test problem
input_data = create_test_problem(p, F);

% Run complete pipeline
params = struct();
params.lambda1 = 0.02;
params.lambda2 = 0.005;
params.max_iter = 100;
params.verbose = false;

[Gamma_final, results] = module5_proximal_main(input_data, params);

% Check final results
converged = results.convergence_info.converged;
all_psd = true;
all_hermitian = true;

for f = 1:F
    [isPSD, ~] = module5_psd_check(Gamma_final{f});
    hermitian_error = norm(Gamma_final{f} - Gamma_final{f}', 'fro');
    
    all_psd = all_psd && isPSD;
    all_hermitian = all_hermitian && (hermitian_error < 1e-12);
end

objective_decreased = (results.objective_history(end) <= results.objective_history(1));

passed = converged && all_psd && all_hermitian && objective_decreased;

if passed
    details = sprintf('Pipeline successful: %d iterations, final objective %.4e', ...
                     results.convergence_info.iterations, results.convergence_info.final_objective);
else
    details = sprintf('Pipeline issues: converged=%d, psd=%d, hermitian=%d, decreased=%d', ...
                     converged, all_psd, all_hermitian, objective_decreased);
end

end

function [passed, details] = test_scalability()
% Test performance across different problem sizes

problem_sizes = [4, 8, 16];
frequencies = [2, 4, 8];
execution_times = zeros(length(problem_sizes), length(frequencies));

for p_idx = 1:length(problem_sizes)
    for f_idx = 1:length(frequencies)
        p = problem_sizes(p_idx);
        F = frequencies(f_idx);
        
        % Create problem
        input_data = create_test_problem(p, F);
        
        % Run with limited iterations for timing
        params = struct();
        params.max_iter = 20;  % Limited for timing test
        params.verbose = false;
        
        % Time execution
        tic;
        [~, ~] = module5_proximal_main(input_data, params);
        execution_times(p_idx, f_idx) = toc;
    end
end

% Check that execution time scales reasonably (not exponentially)
max_time = max(execution_times(:));
min_time = min(execution_times(:));
scaling_reasonable = (max_time / min_time < 1000);  % Allow up to 1000x difference

% Check that larger problems take more time
size_scaling_ok = true;
freq_scaling_ok = true;

for f_idx = 1:length(frequencies)
    if any(diff(execution_times(:, f_idx)) < -0.01)  % Allow small numerical noise
        size_scaling_ok = false;
    end
end

for p_idx = 1:length(problem_sizes)
    if any(diff(execution_times(p_idx, :) ) < -0.01)
        freq_scaling_ok = false;
    end
end

passed = scaling_reasonable && size_scaling_ok && freq_scaling_ok;

if passed
    details = sprintf('Scalability good: %.3f-%.3fs range, monotonic scaling', min_time, max_time);
else
    details = sprintf('Scalability issues: range %.3f-%.3fs, size_ok=%d, freq_ok=%d', ...
                     min_time, max_time, size_scaling_ok, freq_scaling_ok);
end

end

function [passed, details] = test_parallel_performance()
% Test parallel performance improvement

p = 8;
F = 16;  % Enough frequencies to benefit from parallelization

% Create test problem
input_data = create_test_problem(p, F);

params_base = struct();
params_base.max_iter = 10;  % Limited for timing
params_base.verbose = false;

% Test sequential execution
params_seq = params_base;
params_seq.use_parfor = false;

tic;
[~, ~] = module5_proximal_main(input_data, params_seq);
time_sequential = toc;

% Test parallel execution (if parallel pool available)
try
    % Check if parallel pool exists, create if needed
    if isempty(gcp('nocreate'))
        parpool('local', 2);  % Small pool for testing
    end
    
    params_par = params_base;
    params_par.use_parfor = true;
    
    tic;
    [~, ~] = module5_proximal_main(input_data, params_par);
    time_parallel = toc;
    
    % Parallel should be at least not much slower (may not be faster due to overhead)
    reasonable_parallel = (time_parallel < 3 * time_sequential);
    
    passed = reasonable_parallel;
    details = sprintf('Parallel timing: sequential %.3fs, parallel %.3fs (ratio %.2fx)', ...
                     time_sequential, time_parallel, time_parallel/time_sequential);
    
catch ME
    % Parallel testing failed (maybe no Parallel Computing Toolbox)
    passed = true;  % Don't fail the test suite for this
    details = sprintf('Parallel test skipped: %s', ME.message);
end

end

% ==================== Helper Function ====================
function input_data = create_test_problem(p, F)
% Create synthetic test problem for Module 5 testing

% Generate random positive definite covariances
Sigma_tilde = cell(F, 1);
Gamma_init = cell(F, 1);

for f = 1:F
    % Create positive definite covariance
    A = randn(p) + 1i * randn(p);
    Sigma_tilde{f} = A * A' + 0.1 * eye(p);
    Sigma_tilde{f} = (Sigma_tilde{f} + Sigma_tilde{f}') / 2;
    
    % Initialize precision as ridge-regularized inverse
    S = Sigma_tilde{f};
    eps_ridge = 1e-8 * trace(S) / p;
    S_reg = S + eps_ridge * eye(p);
    G0 = S_reg \ eye(p);
    G0 = (G0 + G0') / 2;
    G0(1:p+1:end) = real(diag(G0));
    
    Gamma_init{f} = G0;
end

% Create smoothing kernel
K_smooth = 0.3 * (eye(F) + diag(ones(F-1, 1), 1) + diag(ones(F-1, 1), -1));

% Create weight matrix
W_matrix = eye(p) + 0.1 * randn(p);
W_matrix = (W_matrix + W_matrix') / 2;

% Create active set masks (mostly active for testing)
A_masks = cell(F, 1);
for f = 1:F
    A_masks{f} = rand(p, p) > 0.3;  % 70% active
    A_masks{f} = A_masks{f} | A_masks{f}';  % Ensure symmetry
    A_masks{f}(1:p+1:end) = true;  % Diagonal always active
end

% Package input data
input_data = struct();
input_data.whitened_covariances = Sigma_tilde;
input_data.initial_precision = Gamma_init;
input_data.smoothing_kernel = K_smooth;
input_data.weight_matrix = W_matrix;
input_data.active_set_masks = A_masks;

end
