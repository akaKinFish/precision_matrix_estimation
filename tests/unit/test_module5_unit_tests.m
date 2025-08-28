function test_module5_unit_tests()
% TEST_MODULE5_UNIT_TESTS - Essential unit tests for Module 5 components
%
% Syntax:
%   test_module5_unit_tests()
%
% Description:
%   Runs essential unit tests for Module 5 proximal updates. Focuses on
%   core functionality verification with minimal dependencies.
%
% Tests included:
%   1. Complex soft thresholding correctness
%   2. Hermitian symmetrization
%   3. Active set projection
%   4. PSD checking
%   5. Basic proximal step
%
% See also: TEST_MODULE5_COMPREHENSIVE, MODULE5_PROXIMAL_MAIN
%
% Author: [Your Name]
% Date: [Current Date]
% Version: 1.0

fprintf('\n=== Module 5 Unit Tests ===\n');

total_tests = 0;
passed_tests = 0;

% Test 1: Complex Soft Thresholding
fprintf('Test 1: Complex soft thresholding... ');
total_tests = total_tests + 1;
try
    test_passed = test_soft_thresholding_basic();
    if test_passed
        fprintf('PASSED\n');
        passed_tests = passed_tests + 1;
    else
        fprintf('FAILED\n');
    end
catch ME
    fprintf('ERROR: %s\n', ME.message);
end

% Test 2: Hermitian Symmetrization
fprintf('Test 2: Hermitian symmetrization... ');
total_tests = total_tests + 1;
try
    test_passed = test_hermitian_basic();
    if test_passed
        fprintf('PASSED\n');
        passed_tests = passed_tests + 1;
    else
        fprintf('FAILED\n');
    end
catch ME
    fprintf('ERROR: %s\n', ME.message);
end

% Test 3: Active Set Projection
fprintf('Test 3: Active set projection... ');
total_tests = total_tests + 1;
try
    test_passed = test_active_set_basic();
    if test_passed
        fprintf('PASSED\n');
        passed_tests = passed_tests + 1;
    else
        fprintf('FAILED\n');
    end
catch ME
    fprintf('ERROR: %s\n', ME.message);
end

% Test 4: PSD Check
fprintf('Test 4: PSD checking... ');
total_tests = total_tests + 1;
try
    test_passed = test_psd_check_basic();
    if test_passed
        fprintf('PASSED\n');
        passed_tests = passed_tests + 1;
    else
        fprintf('FAILED\n');
    end
catch ME
    fprintf('ERROR: %s\n', ME.message);
end

% Test 5: Single Proximal Step
fprintf('Test 5: Single proximal step... ');
total_tests = total_tests + 1;
try
    test_passed = test_single_step_basic();
    if test_passed
        fprintf('PASSED\n');
        passed_tests = passed_tests + 1;
    else
        fprintf('FAILED\n');
    end
catch ME
    fprintf('ERROR: %s\n', ME.message);
end

% Summary
fprintf('\n=== Test Summary ===\n');
fprintf('Tests passed: %d/%d\n', passed_tests, total_tests);
if passed_tests == total_tests
    fprintf('🎉 All unit tests passed!\n');
else
    fprintf('❌ Some tests failed. Please check the implementations.\n');
end
fprintf('===================\n');

end

% ==================== Individual Test Functions ====================

function passed = test_soft_thresholding_basic()
% Test basic soft thresholding functionality

% Create test matrix
p = 3;
test_matrix = [1.0,        0.15+0.2i,  0.05+0.03i;
               0.15-0.2i,  1.0,        0.12-0.1i;
               0.05-0.03i, 0.12+0.1i,  1.0];

tau = 0.1;

% Test simplified mode
result = module5_soft_threshold_complex(test_matrix, tau, 'simplified');

% Check Hermitian symmetry
hermitian_error = norm(result - result', 'fro');

% Check diagonal is real
diag_imag_error = max(abs(imag(diag(result))));

% Check thresholding worked (element (1,2) should be thresholded)
original_magnitude = abs(test_matrix(1,2));
result_magnitude = abs(result(1,2));

if original_magnitude > tau
    expected_magnitude = original_magnitude - tau;
    threshold_correct = abs(result_magnitude - expected_magnitude) < 1e-12;
else
    threshold_correct = result_magnitude < 1e-12;
end

% Check phase preservation for surviving elements
phase_preserved = true;
if result_magnitude > 1e-12
    original_phase = angle(test_matrix(1,2));
    result_phase = angle(result(1,2));
    phase_diff = abs(original_phase - result_phase);
    phase_preserved = (phase_diff < 1e-10) || (abs(phase_diff - 2*pi) < 1e-10);
end

passed = (hermitian_error < 1e-12) && (diag_imag_error < 1e-12) && ...
         threshold_correct && phase_preserved;

end

function passed = test_hermitian_basic()
% Test basic Hermitian symmetrization

p = 2;
% Create non-Hermitian matrix
test_matrix = [1.5+0.1i,  0.2+0.3i;
               0.3-0.4i,  2.0+0.05i];

% Test simplified mode
result_simp = module5_hermitian_symmetrize(test_matrix, 'simplified');

% Test joint mode
result_joint = module5_hermitian_symmetrize(test_matrix, 'joint');

% Check Hermitian properties
hermitian_simp = norm(result_simp - result_simp', 'fro') < 1e-14;
hermitian_joint = norm(result_joint - result_joint', 'fro') < 1e-14;

% Check diagonal properties
diag_real_simp = all(abs(imag(diag(result_simp))) < 1e-14);
diag_ones_joint = norm(diag(result_joint) - ones(p, 1)) < 1e-14;

passed = hermitian_simp && hermitian_joint && diag_real_simp && diag_ones_joint;

end

function passed = test_active_set_basic()
% Test basic active set projection

p = 3;
% Create test matrix and mask
test_matrix = randn(p) + 1i*randn(p);
test_matrix = (test_matrix + test_matrix')/2;  % Make Hermitian

% Create symmetric active mask
active_mask = logical([1, 1, 0; 1, 1, 1; 0, 1, 1]);

% Apply projection
result = module5_active_set_projection(test_matrix, active_mask);

% Check inactive elements are zero
inactive_positions = ~active_mask;
max_inactive = max(abs(result(inactive_positions)));

% Check active elements preserved
active_positions = active_mask;
preservation_error = norm(result(active_positions) - test_matrix(active_positions));

% Check Hermitian symmetry preserved
hermitian_preserved = norm(result - result', 'fro') < 1e-12;

passed = (max_inactive < 1e-15) && (preservation_error < 1e-12) && hermitian_preserved;

end

function passed = test_psd_check_basic()
% Test basic PSD checking

p = 2;

% Create positive definite matrix
A = randn(p) + 1i*randn(p);
psd_matrix = A'*A + eye(p);

% Create non-PSD matrix
non_psd_matrix = [1, 2; 2, 1];  % Eigenvalues: 3, -1

% Test PSD check
[isPSD_good, ~] = module5_psd_check(psd_matrix);
[isPSD_bad, ~] = module5_psd_check(non_psd_matrix);

passed = isPSD_good && ~isPSD_bad;

end

function passed = test_single_step_basic()
% Test basic single proximal step

p = 2;

% Create test data
Gamma_current = [2, 0.1+0.2i; 0.1-0.2i, 2];
gradient = [0.1, 0.05-0.1i; 0.05+0.1i, 0.1];
alpha = 0.1;
active_mask = true(p, p);

aux_data = struct();
aux_data.lambda2 = 0.01;

params = struct();
params.mode = 'simplified';
params.beta_backtrack = 0.5;
params.max_backtrack = 10;

% Perform single step
[Gamma_new, step_info] = module5_single_proximal_step(...
    Gamma_current, gradient, alpha, active_mask, aux_data, params);

% Check basic properties
step_successful = step_info.success;
result_hermitian = norm(Gamma_new - Gamma_new', 'fro') < 1e-12;
result_psd = module5_psd_check(Gamma_new);
matrix_changed = (step_info.relative_change > 1e-12);

passed = step_successful && result_hermitian && result_psd && matrix_changed;

end