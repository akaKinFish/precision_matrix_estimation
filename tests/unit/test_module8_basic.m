function test_results = test_module8_basic()
% TEST_MODULE8_BASIC - Basic unit tests for Module 8 recoloring functionality
%
% Coverage:
%  1) Formula correctness
%  2) Cycle consistency (Omega * Sigma ≈ I)
%  3) Hermitian + SPD preservation
%  4) Sparsity pattern preservation (with correct threshold scaling by d_i d_j)
%  5) Complex Hermitian support
%  6) Numerical stability (condition number sweep)
%  7) Edge cases

fprintf('========================================\n');
fprintf('Module 8 Basic Unit Test Suite\n');
fprintf('========================================\n\n');

test_results = struct();
test_results.timestamp = datestr(now);
test_results.total_tests = 0;
test_results.passed_tests = 0;
test_results.failed_tests = 0;
test_results.test_details = {};

% U1
fprintf('=== Test U1: Mathematical Formula Correctness ===\n');
[test_results.total_tests, test_results.passed_tests] = run_test(@test_formula_correctness, ...
    'Mathematical Formula', test_results.total_tests, test_results.passed_tests, test_results);

% U2
fprintf('\n=== Test U2: Cycle Consistency ===\n');
[test_results.total_tests, test_results.passed_tests] = run_test(@test_cycle_consistency, ...
    'Cycle Consistency', test_results.total_tests, test_results.passed_tests, test_results);

% U3
fprintf('\n=== Test U3: Hermitian and SPD Properties ===\n');
[test_results.total_tests, test_results.passed_tests] = run_test(@test_hermitian_spd_properties, ...
    'Hermitian and SPD', test_results.total_tests, test_results.passed_tests, test_results);

% U4
fprintf('\n=== Test U4: Sparsity Pattern Preservation ===\n');
[test_results.total_tests, test_results.passed_tests] = run_test(@test_sparsity_preservation, ...
    'Sparsity Preservation', test_results.total_tests, test_results.passed_tests, test_results);

% U5
fprintf('\n=== Test U5: Complex Hermitian Matrix Support ===\n');
[test_results.total_tests, test_results.passed_tests] = run_test(@test_complex_hermitian, ...
    'Complex Hermitian', test_results.total_tests, test_results.passed_tests, test_results);

% U6
fprintf('\n=== Test U6: Numerical Stability ===\n');
[test_results.total_tests, test_results.passed_tests] = run_test(@test_numerical_stability, ...
    'Numerical Stability', test_results.total_tests, test_results.passed_tests, test_results);

% U7
fprintf('\n=== Test U7: Edge Cases and Boundary Conditions ===\n');
[test_results.total_tests, test_results.passed_tests] = run_test(@test_edge_cases, ...
    'Edge Cases', test_results.total_tests, test_results.passed_tests, test_results);

% Summary
test_results.failed_tests = test_results.total_tests - test_results.passed_tests;
test_results.success_rate = test_results.passed_tests / test_results.total_tests;

fprintf('\n========================================\n');
fprintf('Module 8 Basic Test Summary\n');
fprintf('========================================\n');
fprintf('Total tests:   %d\n', test_results.total_tests);
fprintf('Passed tests:  %d\n', test_results.passed_tests);
fprintf('Failed tests:  %d\n', test_results.failed_tests);
fprintf('Success rate:  %.1f%%\n', 100 * test_results.success_rate);
if test_results.success_rate >= 0.9
    fprintf('Result: ✓ EXCELLENT - Module 8 basic functionality verified\n');
elseif test_results.success_rate >= 0.8
    fprintf('Result: ⚠ GOOD - Minor issues may exist\n');
else
    fprintf('Result: ✗ NEEDS WORK - Significant issues detected\n');
end
fprintf('\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Individual Test Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [passed, details] = test_formula_correctness()
passed = true; details = 'Formula correctness validation';
try
    p = 6; F = 3;
    for omega = 1:F
        A = randn(p, p);
        Sigma_original = (A * A') + 0.1 * eye(p);
        g_vals = 0.5 + 2 * rand(p, 1);
        D_omega = diag(1 ./ sqrt(g_vals));
        Sigma_whitened = D_omega * Sigma_original * D_omega;
        Gamma_tilde_theoretical = inv(Sigma_whitened);

        input_data.whitened_precision_matrices = {Gamma_tilde_theoretical};
        input_data.whitening_matrices = {D_omega};
        results = module8_recoloring_main(input_data, struct('verbose', false));
        Omega_computed = results.recolored_precision_matrices{1};

        Omega_theoretical = inv(Sigma_original);
        formula_error = norm(Omega_computed - Omega_theoretical, 'fro') / norm(Omega_theoretical, 'fro');

        machine_eps = eps('double');
        condition_D = cond(D_omega);
        tolerance = 1e3 * machine_eps * condition_D^2;

        if formula_error > tolerance
            passed = false;
            details = sprintf('Formula error %.2e > tolerance %.2e at ω=%d', ...
                formula_error, tolerance, omega);
            return;
        end
    end
    details = sprintf('Formula correctness verified for %d frequencies', F);
catch ME
    passed = false;
    details = sprintf('Formula test failed: %s', ME.message);
end
end

function [passed, details] = test_cycle_consistency()
passed = true; details = 'Cycle consistency validation';
try
    p = 8; F = 4; max_cycle_error = 0;
    for omega = 1:F
        [Sigma_original, D_omega, Gamma_tilde] = create_test_matrices(p);
        input_data.whitened_precision_matrices = {Gamma_tilde};
        input_data.whitening_matrices = {D_omega};
        input_data.original_covariances = {Sigma_original};
        params = struct('verbose', false, 'compute_quality_metrics', true);
        results = module8_recoloring_main(input_data, params);
        quality = results.recoloring_quality{1};
        cycle_error = quality.inv_error;
        max_cycle_error = max(max_cycle_error, cycle_error);

        machine_eps = eps('double');
        condition_D = cond(D_omega);
        tolerance = 1e3 * machine_eps * condition_D^2;

        if cycle_error > tolerance
            passed = false;
            details = sprintf('Cycle error %.2e > tolerance %.2e at ω=%d', ...
                cycle_error, tolerance, omega);
            return;
        end
    end
    details = sprintf('Cycle consistency verified: max error %.2e across %d frequencies', ...
        max_cycle_error, F);
catch ME
    passed = false;
    details = sprintf('Cycle consistency test failed: %s', ME.message);
end
end

function [passed, details] = test_hermitian_spd_properties()
passed = true; details = 'Hermitian and SPD property validation';
try
    p = 5; F = 3;
    for omega = 1:F
        [Sigma_original, D_omega, Gamma_tilde] = create_test_matrices(p);
        Gamma_tilde = (Gamma_tilde + Gamma_tilde') / 2;
        eigenvals = eig(Gamma_tilde);
        if any(real(eigenvals) <= 0)
            Gamma_tilde = Gamma_tilde + (abs(min(real(eigenvals))) + 0.01) * eye(p);
        end
        input_data.whitened_precision_matrices = {Gamma_tilde};
        input_data.whitening_matrices = {D_omega};
        params = struct('verbose', false, 'force_hermitian', true, 'validate_spd', true);
        results = module8_recoloring_main(input_data, params);
        quality = results.recoloring_quality{1};

        hermitian_tolerance = 1e-12;
        if quality.hermitian_error > hermitian_tolerance
            passed = false;
            details = sprintf('Hermitian error %.2e > %.2e at ω=%d', ...
                quality.hermitian_error, hermitian_tolerance, omega);
            return;
        end
        if ~quality.spd_status.is_spd
            passed = false;
            details = sprintf('SPD property lost at ω=%d: %s', omega, quality.spd_status.details.message);
            return;
        end
    end
    details = sprintf('Hermitian and SPD properties preserved across %d frequencies', F);
catch ME
    passed = false;
    details = sprintf('Hermitian/SPD test failed: %s', ME.message);
end
end

function [passed, details] = test_sparsity_preservation()
passed = true; details = 'Sparsity pattern preservation validation';
try
    p = 6; F = 2;
    for omega = 1:F
        [Sigma_original, D_omega] = create_sparse_test_matrices(p, 0.3);
        Sigma_whitened = D_omega * Sigma_original * D_omega;
        Gamma_tilde = inv(Sigma_whitened);

        threshold_gamma = 0.1 * max(abs(Gamma_tilde(:)));
        Gamma_tilde_sparse = Gamma_tilde;
        Gamma_tilde_sparse(abs(Gamma_tilde) < threshold_gamma) = 0;
        Gamma_tilde_sparse = (Gamma_tilde_sparse + Gamma_tilde_sparse') / 2;
        eigenvals = eig(Gamma_tilde_sparse);
        if any(real(eigenvals) <= 0)
            Gamma_tilde_sparse = Gamma_tilde_sparse + (abs(min(real(eigenvals))) + 0.01) * eye(p);
        end

        input_data.whitened_precision_matrices = {Gamma_tilde_sparse};
        input_data.whitening_matrices = {D_omega};
        results = module8_recoloring_main(input_data, struct('verbose', false));
        Omega_recolored = results.recolored_precision_matrices{1};

        d_vec = diag(D_omega);
        original_nonzeros = (abs(Gamma_tilde_sparse) > threshold_gamma);

        % Correct threshold scaling: multiply by (d_i d_j)
        threshold_matrix = threshold_gamma * (d_vec * d_vec.');
        recolored_nonzeros = (abs(Omega_recolored) > threshold_matrix);

        total_elements = p * p;
        consistent_elements = sum((original_nonzeros(:) & recolored_nonzeros(:)) | ...
                                  (~original_nonzeros(:) & ~recolored_nonzeros(:)));
        consistency_rate = consistent_elements / total_elements;

        min_consistency = 0.95;
        if consistency_rate < min_consistency
            passed = false;
            details = sprintf('Sparsity consistency %.3f < %.3f at ω=%d', ...
                consistency_rate, min_consistency, omega);
            return;
        end
    end
    details = sprintf('Sparsity patterns preserved across %d frequencies', F);
catch ME
    passed = false;
    details = sprintf('Sparsity test failed: %s', ME.message);
end
end

function [passed, details] = test_complex_hermitian()
passed = true; details = 'Complex Hermitian matrix support validation';
try
    p = 4; F = 2;
    for omega = 1:F
        A = randn(p, p) + 1i * randn(p, p);
        lambda = 0.1;
        Gamma_tilde_complex = A * A' + lambda * eye(p);
        Gamma_tilde_complex = (Gamma_tilde_complex + Gamma_tilde_complex') / 2;

        d_vals = 0.5 + rand(p, 1);
        D_omega = diag(1 ./ sqrt(d_vals));

        input_data.whitened_precision_matrices = {Gamma_tilde_complex};
        input_data.whitening_matrices = {D_omega};
        params = struct('verbose', false, 'force_hermitian', true, 'validate_spd', true);
        results = module8_recoloring_main(input_data, params);
        Omega_complex = results.recolored_precision_matrices{1};
        quality = results.recoloring_quality{1};

        if quality.hermitian_error > 1e-12
            passed = false;
            details = sprintf('Complex Hermitian error %.2e > 1e-12 at ω=%d', ...
                quality.hermitian_error, omega); return;
        end
        if ~quality.spd_status.is_spd
            passed = false;
            details = sprintf('Complex matrix lost SPD property at ω=%d', omega); return;
        end
        if isreal(Omega_complex) && ~isreal(Gamma_tilde_complex)
            passed = false;
            details = sprintf('Complex input produced real output at ω=%d', omega); return;
        end
    end
    details = sprintf('Complex Hermitian matrices processed correctly across %d frequencies', F);
catch ME
    passed = false;
    details = sprintf('Complex Hermitian test failed: %s', ME.message);
end
end

function [passed, details] = test_numerical_stability()
passed = true; details = 'Numerical stability validation';
try
    p = 5; condition_numbers = [1e2, 1e4, 1e6, 1e8]; max_observed_error = 0;
    for cond_target = condition_numbers
        [Sigma_original, D_omega, Gamma_tilde] = create_test_matrices_with_condition(p, cond_target);
        input_data.whitened_precision_matrices = {Gamma_tilde};
        input_data.whitening_matrices = {D_omega};
        input_data.original_covariances = {Sigma_original};
        results = module8_recoloring_main(input_data, struct('verbose', false));
        quality = results.recoloring_quality{1};

        inv_error = quality.inv_error; max_observed_error = max(max_observed_error, inv_error);
        machine_eps = eps('double'); condition_D = cond(D_omega);
        expected_tolerance = 1e3 * machine_eps * condition_D^2;

        if inv_error > expected_tolerance
            passed = false;
            details = sprintf('Error %.2e > tolerance %.2e for κ=%.1e', ...
                inv_error, expected_tolerance, cond_target);
            return;
        end
    end
    details = sprintf('Numerical stability verified: max error %.2e across κ ∈ [1e2, 1e8]', max_observed_error);
catch ME
    passed = false;
    details = sprintf('Numerical stability test failed: %s', ME.message);
end
end

function [passed, details] = test_edge_cases()
passed = true; details = 'Edge cases validation';
try
    test_cases = {'minimal_size', 'diagonal_matrix', 'near_singular', 'high_sparsity'};
    for i = 1:numel(test_cases)
        case_name = test_cases{i};
        [input_data, expected_behavior] = create_edge_case_data(case_name);
        results = module8_recoloring_main(input_data, struct('verbose', false));
        if ~validate_edge_case_result(results, expected_behavior, case_name)
            passed = false;
            details = sprintf('Edge case "%s" failed validation', case_name);
            return;
        end
    end
    details = sprintf('All %d edge cases passed validation', numel(test_cases));
catch ME
    passed = false;
    details = sprintf('Edge cases test failed: %s', ME.message);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helpers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Sigma_original, D_omega, Gamma_tilde] = create_test_matrices(p)
A = randn(p, p);
Sigma_original = (A * A') + 0.1 * eye(p);
g_vals = 0.2 + 1.8 * rand(p, 1);
D_omega = diag(1 ./ sqrt(g_vals));
Sigma_whitened = D_omega * Sigma_original * D_omega;
Gamma_tilde = inv(Sigma_whitened);
end

function [Sigma_original, D_omega] = create_sparse_test_matrices(p, sparsity_level)
Omega_true = 2 * eye(p);
n_connections = round(sparsity_level * p * (p-1) / 2);
[I, J] = find(triu(ones(p), 1));
idxs = randperm(numel(I), n_connections);
for k = 1:n_connections
    i = I(idxs(k)); j = J(idxs(k));
    val = 0.2 * (rand() - 0.5);
    Omega_true(i, j) = val; Omega_true(j, i) = val;
end
Sigma_original = inv(Omega_true);
g_vals = 0.5 + rand(p, 1);
D_omega = diag(1 ./ sqrt(g_vals));
end

function [Sigma_original, D_omega, Gamma_tilde] = create_test_matrices_with_condition(p, target_condition)
% Create matrices with a controlled condition number. Ensure all vectors are column vectors.

U = orth(randn(p, p));
lambda_max = 1.0; lambda_min = lambda_max / target_condition;
lambdas = linspace(lambda_min, lambda_max, p).';   % force column vector (p x 1)
Sigma_original = U * diag(lambdas) * U';

% g_vals must be (p x 1); avoid broadcasting to (p x p)
g_vals = sqrt(lambdas) .* (0.8 + 0.4 * rand(p, 1));  % (p x 1) .* (p x 1) -> (p x 1)
D_omega = diag(1 ./ sqrt(g_vals));

Sigma_whitened = D_omega * Sigma_original * D_omega;
Gamma_tilde = inv(Sigma_whitened);
end

function [input_data, expected_behavior] = create_edge_case_data(case_name)
switch case_name
    case 'minimal_size'
        Gamma_tilde = {2.5}; D = {diag(0.8)};
        expected_behavior.should_succeed = true;
        expected_behavior.expected_result = 2.5 * 0.8^2;
    case 'diagonal_matrix'
        Gamma_tilde = {diag([1, 2, 3, 4])}; D = {diag([0.5, 0.8, 1.2, 0.9])};
        expected_behavior.should_succeed = true; expected_behavior.should_be_diagonal = true;
    case 'near_singular'
        p = 3; g_vals = [1, 1, 1e-10];
        D = {diag(1 ./ sqrt(g_vals))};
        X = randn(p); Gamma_tilde = {eye(p) + 0.1 * (X + X')/2};
        expected_behavior.should_succeed = true; expected_behavior.may_have_warnings = true;
    case 'high_sparsity'
        p = 8; Gamma_sparse = 3 * eye(p);
        Gamma_tilde = {Gamma_sparse}; D = {diag(0.5 + rand(p, 1))};
        expected_behavior.should_succeed = true; expected_behavior.should_be_sparse = true;
    otherwise
        error('Unknown edge case: %s', case_name);
end
input_data.whitened_precision_matrices = Gamma_tilde;
input_data.whitening_matrices = D;
end

function is_valid = validate_edge_case_result(results, expected_behavior, case_name)
is_valid = true;
try
    if expected_behavior.should_succeed && ~results.success, is_valid = false; return; end
    Omega_result = results.recolored_precision_matrices{1};
    switch case_name
        case 'minimal_size'
            actual = Omega_result(1,1); expected = expected_behavior.expected_result;
            err = abs(actual - expected) / max(abs(expected), eps);
            is_valid = (err < 1e-12);
        case 'diagonal_matrix'
            off = norm(Omega_result - diag(diag(Omega_result)), 'fro');
            is_valid = (off < 1e-12);
        case 'high_sparsity'
            dp = norm(diag(Omega_result))^2; tp = norm(Omega_result, 'fro')^2;
            is_valid = (dp / max(tp, eps) > 0.9);
        otherwise
            is_valid = results.success;
    end
catch
    is_valid = false;
end
end

function [total_tests, passed_tests] = run_test(test_function, test_name, total_tests, passed_tests, test_results)
total_tests = total_tests + 1;
try
    [passed, details] = test_function();
    if passed
        passed_tests = passed_tests + 1;
        fprintf('✓ %s: PASSED - %s\n', test_name, details); status = 'PASSED';
    else
        fprintf('✗ %s: FAILED - %s\n', test_name, details); status = 'FAILED';
    end
catch ME
    fprintf('✗ %s: ERROR - %s\n', test_name, ME.message);
    details = sprintf('Test execution error: %s', ME.message); status = 'ERROR';
end
rec.name = test_name; rec.status = status; rec.details = details; rec.timestamp = datetime('now');
test_results.test_details{end+1} = rec;
end
