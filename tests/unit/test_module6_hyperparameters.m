function test_results = test_module6_hyperparameters()
% TEST_MODULE6_HYPERPARAMETERS - Basic test suite for Module 6

rng(42);  % reproducible

fprintf('============================================\n');
fprintf('Module 6 Hyperparameter Configuration Tests\n');
fprintf('============================================\n\n');

test_results = struct();
test_results.total_tests = 0;
test_results.passed_tests = 0;
test_results.failed_tests = 0;
test_results.test_details = {};
test_results.success_rate = 0;

%% Test 1: Basic functionality
fprintf('=== Test 1: Basic Functionality ===\n');
[test_results.total_tests, test_results.passed_tests] = ...
    run_test(@test_basic_functionality, 'Basic Functionality', ...
             test_results.total_tests, test_results.passed_tests, test_results);

%% Test 2: Input validation
fprintf('\n=== Test 2: Input Validation ===\n');
[test_results.total_tests, test_results.passed_tests] = ...
    run_test(@test_input_validation, 'Input Validation', ...
             test_results.total_tests, test_results.passed_tests, test_results);

%% Test 3: Parameter ranges
fprintf('\n=== Test 3: Parameter Ranges ===\n');
[test_results.total_tests, test_results.passed_tests] = ...
    run_test(@test_parameter_ranges, 'Parameter Ranges', ...
             test_results.total_tests, test_results.passed_tests, test_results);

%% Test 4: Method comparison
fprintf('\n=== Test 4: Method Comparison ===\n');
[test_results.total_tests, test_results.passed_tests] = ...
    run_test(@test_method_comparison, 'Method Comparison', ...
             test_results.total_tests, test_results.passed_tests, test_results);

%% Summary
test_results.failed_tests = test_results.total_tests - test_results.passed_tests;
test_results.success_rate = test_results.passed_tests / max(1, test_results.total_tests);

fprintf('\n============================================\n');
fprintf('Module 6 Test Summary\n');
fprintf('============================================\n');
fprintf('Total tests:   %d\n', test_results.total_tests);
fprintf('Passed tests:  %d\n', test_results.passed_tests);
fprintf('Failed tests:  %d\n', test_results.failed_tests);
fprintf('Success rate:  %.1f%%\n', 100 * test_results.success_rate);

if test_results.success_rate >= 0.9
    fprintf('Result: ✓ EXCELLENT - Module 6 is production ready\n');
elseif test_results.success_rate >= 0.8
    fprintf('Result: ⚠ GOOD - Minor issues may exist\n');
else
    fprintf('Result: ✗ NEEDS WORK - Significant issues detected\n');
end

test_results.summary = sprintf('Module 6: %d/%d tests passed (%.1f%%)', ...
    test_results.passed_tests, test_results.total_tests, 100 * test_results.success_rate);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Individual Test Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [passed, details] = test_basic_functionality()
try
    test_input = create_simple_test_input(4, 3);
    output = module6_hyperparameters(test_input, 'verbose', false);

    required_fields = {'lambda1', 'alpha', 'lambda2_suggested'};
    for i = 1:numel(required_fields)
        if ~isfield(output, required_fields{i})
            passed = false; details = sprintf('Missing field: %s', required_fields{i}); return;
        end
    end

    params_positive = (output.lambda1 > 0) && (output.alpha > 0) && (output.lambda2_suggested > 0);
    params_finite   = isfinite(output.lambda1) && isfinite(output.alpha) && isfinite(output.lambda2_suggested);

    if ~params_positive
        passed = false;
        details = sprintf('Parameters not positive: λ₁=%.2e, α=%.2e, λ₂=%.2e', ...
            output.lambda1, output.alpha, output.lambda2_suggested);
        return;
    end
    if ~params_finite
        passed = false; details = 'Parameters are not finite'; return;
    end

    passed = true;
    details = sprintf('OK: λ₁=%.3e, α=%.3e, λ₂=%.3e', ...
        output.lambda1, output.alpha, output.lambda2_suggested);
catch ME
    passed = false;
    details = sprintf('Basic functionality test failed: %s', ME.message);
end
end

function [passed, details] = test_input_validation()
try
    try
        module6_hyperparameters([]); %#ok<*NASGU>
        passed = false; details = 'Should reject empty input'; return;
    catch ME
        if ~contains(ME.identifier, 'invalid_input')
            passed = false; details = 'Wrong error type for empty input'; return;
        end
    end

    incomplete = struct(); incomplete.whitened_covariances = {};
    try
        module6_hyperparameters(incomplete);
        passed = false; details = 'Should reject incomplete input structure'; return;
    catch ME
        if ~contains(ME.identifier, 'missing_field')
            passed = false; details = 'Wrong error type for missing fields'; return;
        end
    end

    valid = create_simple_test_input(4, 3);
    out = module6_hyperparameters(valid, 'verbose', false);
    if ~isstruct(out)
        passed = false; details = 'Valid input should produce struct output'; return;
    end

    passed = true; details = 'Input validation tests passed';
catch ME
    passed = false;
    details = sprintf('Input validation test failed: %s', ME.message);
end
end

function [passed, details] = test_parameter_ranges()
try
    test_input = create_simple_test_input(5, 4);

    safety_margins = [0.5, 0.9, 0.95];
    lambda1_values = zeros(size(safety_margins));
    alpha_values   = zeros(size(safety_margins));

    for i = 1:numel(safety_margins)
        cfg = module6_hyperparameters(test_input, 'safety_margin', safety_margins(i), 'verbose', false);
        lambda1_values(i) = cfg.lambda1;
        alpha_values(i)   = cfg.alpha;
    end

    lambda1_decreasing = all(diff(lambda1_values) <= 1e-10);
    alpha_decreasing   = all(diff(alpha_values)   <= 1e-10);

    lambda1_reasonable = all(lambda1_values > 1e-8) && all(lambda1_values < 1);
    alpha_reasonable   = all(alpha_values   > 1e-6) && all(alpha_values   < 1);

    passed = lambda1_decreasing && alpha_decreasing && lambda1_reasonable && alpha_reasonable;

    if passed
        details = sprintf('Ranges OK: λ₁∈[%.2e,%.2e], α∈[%.2e,%.2e]', ...
            min(lambda1_values), max(lambda1_values), min(alpha_values), max(alpha_values));
    else
        details = sprintf('Range issues: λ₁_mono=%s, α_mono=%s, λ₁_range=%s, α_range=%s', ...
            tern(lambda1_decreasing, 'OK', 'FAIL'), tern(alpha_decreasing, 'OK', 'FAIL'), ...
            tern(lambda1_reasonable, 'OK', 'FAIL'), tern(alpha_reasonable, 'OK', 'FAIL'));
    end
catch ME
    passed = false;
    details = sprintf('Parameter range test failed: %s', ME.message);
end
end

function [passed, details] = test_method_comparison()
try
    test_input = create_simple_test_input(6, 3);

    cfg_exact = module6_hyperparameters(test_input, 'use_gershgorin', false, 'verbose', false);
    cfg_gersh = module6_hyperparameters(test_input, 'use_gershgorin', true,  'verbose', false);

    exact_valid = isfield(cfg_exact, 'lambda1') && isfield(cfg_exact, 'alpha');
    gersh_valid = isfield(cfg_gersh, 'lambda1') && isfield(cfg_gersh, 'alpha');
    if ~exact_valid || ~gersh_valid
        passed = false;
        details = sprintf('Method validation: exact=%s, gershgorin=%s', ...
            tern(exact_valid, 'OK', 'FAIL'), tern(gersh_valid, 'OK', 'FAIL'));
        return;
    end

    exact_positive = (cfg_exact.lambda1 > 0) && (cfg_exact.alpha > 0);
    gersh_positive = (cfg_gersh.lambda1 > 0) && (cfg_gersh.alpha > 0);

    L_exact = cfg_exact.diagnostics.L_logdet;
    L_gersh = cfg_gersh.diagnostics.L_logdet;
    bound_reasonable = L_gersh >= L_exact - 1e-10;  % Gershgorin is conservative

    passed = exact_positive && gersh_positive && bound_reasonable;

    if passed
        details = sprintf('Method comparison OK: L_ratio=%.2f', L_gersh / max(1e-12, L_exact));
    else
        details = sprintf('Method issues: exact_pos=%s, gersh_pos=%s, bound_ok=%s', ...
            tern(exact_positive, 'OK', 'FAIL'), tern(gersh_positive, 'OK', 'FAIL'), ...
            tern(bound_reasonable, 'OK', 'FAIL'));
    end
catch ME
    passed = false;
    details = sprintf('Method comparison test failed: %s', ME.message);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helpers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function input_struct = create_simple_test_input(n_nodes, n_freq)
input_struct = struct();
input_struct.whitened_covariances = cell(n_freq, 1);
for f = 1:n_freq
    A = randn(n_nodes); A = A + A';
    A = A + (abs(min(eig(A))) + 0.1) * eye(n_nodes);
    input_struct.whitened_covariances{f} = (A + A')/2;
end
input_struct.kernel_matrix = eye(n_freq);
input_struct.weight_matrix = eye(n_nodes); % diag ignored in R_max computation
end

function [total_tests, passed_tests] = run_test(fun, name, total_tests, passed_tests, test_results)
total_tests = total_tests + 1;
try
    [ok, details] = fun();
    if ok
        passed_tests = passed_tests + 1;
        fprintf('✓ %s: PASSED - %s\n', name, details);
        status = 'PASSED';
    else
        fprintf('✗ %s: FAILED - %s\n', name, details);
        status = 'FAILED';
    end
catch ME
    fprintf('✗ %s: ERROR - %s\n', name, ME.message);
    details = sprintf('Test execution error: %s', ME.message);
    status = 'ERROR';
end

rec = struct(); rec.name = name; rec.status = status; rec.details = details; rec.timestamp = datetime("now");
test_results.test_details{end+1} = rec; 
end

function s = tern(cond, a, b)
if cond, s = a; else, s = b; end
end
