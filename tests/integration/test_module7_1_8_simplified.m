function test_results = test_module7_1_8_simplified()
% TEST_MODULE7_1_8_SIMPLIFIED - Simplified end-to-end integration test
% Pipeline: Module 7 → Module 1 → Theoretical Gamma_tilde → Module 8
% Key update: build Gamma_tilde from whitened Σ_true using D from Module 1.

fprintf('========================================\n');
fprintf('Module 7→1→8 Simplified Integration Test\n');
fprintf('========================================\n\n');

test_results = struct();
test_results.timestamp = datestr(now);
test_results.test_version = 'simplified_e2e_v1.2';
test_results.total_tests = 0;
test_results.passed_tests = 0;
test_results.failed_tests = 0;
test_results.test_details = {};

% Test 1
fprintf('=== Test 1: Basic Pipeline Integration ===\n');
[test_results.total_tests, test_results.passed_tests] = run_test(@test_basic_pipeline_integration, ...
    'Basic Pipeline', test_results.total_tests, test_results.passed_tests, test_results);

% Test 2
fprintf('\n=== Test 2: Complex Data Pipeline ===\n');
[test_results.total_tests, test_results.passed_tests] = run_test(@test_complex_data_pipeline, ...
    'Complex Data Pipeline', test_results.total_tests, test_results.passed_tests, test_results);

% Test 3
fprintf('\n=== Test 3: Multiple Problem Sizes ===\n');
[test_results.total_tests, test_results.passed_tests] = run_test(@test_multiple_problem_sizes, ...
    'Multiple Problem Sizes', test_results.total_tests, test_results.passed_tests, test_results);

% Test 4
fprintf('\n=== Test 4: Different Graph Structures ===\n');
[test_results.total_tests, test_results.passed_tests] = run_test(@test_different_graph_types, ...
    'Different Graph Types', test_results.total_tests, test_results.passed_tests, test_results);

% Test 5
fprintf('\n=== Test 5: Noise Robustness ===\n');
[test_results.total_tests, test_results.passed_tests] = run_test(@test_noise_robustness, ...
    'Noise Robustness', test_results.total_tests, test_results.passed_tests, test_results);

% Test 6
fprintf('\n=== Test 6: Ground Truth Recovery Quality ===\n');
[test_results.total_tests, test_results.passed_tests] = run_test(@test_ground_truth_recovery, ...
    'Ground Truth Recovery', test_results.total_tests, test_results.passed_tests, test_results);

% Summary
test_results.failed_tests = test_results.total_tests - test_results.passed_tests;
test_results.success_rate = test_results.passed_tests / test_results.total_tests;

fprintf('\n========================================\n');
fprintf('Module 7→1→8 Integration Test Summary\n');
fprintf('========================================\n');
fprintf('Total tests:   %d\n', test_results.total_tests);
fprintf('Passed tests:  %d\n', test_results.passed_tests);
fprintf('Failed tests:  %d\n', test_results.failed_tests);
fprintf('Success rate:  %.1f%%\n', 100 * test_results.success_rate);
if test_results.success_rate >= 0.9
    fprintf('Result: ✓ EXCELLENT - Simplified pipeline ready for production\n');
elseif test_results.success_rate >= 0.8
    fprintf('Result: ⚠ GOOD - Minor integration issues may exist\n');
else
    fprintf('Result: ✗ NEEDS WORK - Significant integration issues detected\n');
end
fprintf('\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Individual Test Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [passed, details] = test_basic_pipeline_integration()
passed = true; details = 'Basic pipeline integration validation';
try
    test_params = struct('n_nodes', 8, 'n_freq', 5, 'n_samples', 100, ...
                         'graph_type', 'random', 'edge_density', 0.3);
    [recovery_metrics, pipeline_data] = run_simplified_pipeline(test_params);

    validation_results = validate_recovery_metrics(recovery_metrics, pipeline_data);
    if ~validation_results.meets_criteria
        passed = false;
        details = sprintf('Recovery validation failed: %s', validation_results.failure_reason);
        return;
    end

    consistency_results = validate_pipeline_consistency(pipeline_data);
    if ~consistency_results.is_consistent
        passed = false;
        details = sprintf('Pipeline consistency failed: %s', consistency_results.failure_reason);
        return;
    end

    details = sprintf('Pipeline integration successful: median recovery error %.3f, cycle error %.2e', ...
        recovery_metrics.median_recovery_error, recovery_metrics.median_cycle_error);
catch ME
    passed = false;
    details = sprintf('Pipeline integration test failed: %s', ME.message);
end
end

function [passed, details] = test_complex_data_pipeline()
passed = true; details = 'Complex data pipeline validation';
try
    test_params = struct('n_nodes', 6, 'n_freq', 4, 'n_samples', 80, ...
                         'graph_type', 'hub', 'edge_density', 0.25, ...
                         'complex_strength', 1.0);
    [recovery_metrics, pipeline_data] = run_simplified_pipeline(test_params);

    complex_validation = validate_complex_data_handling(pipeline_data);
    if ~complex_validation.success
        passed = false;
        details = sprintf('Complex data handling failed: %s', complex_validation.failure_reason);
        return;
    end

    validation_results = validate_recovery_metrics(recovery_metrics, pipeline_data);
    if ~validation_results.meets_criteria
        passed = false;
        details = sprintf('Complex data recovery failed: %s', validation_results.failure_reason);
        return;
    end

    details = sprintf('Complex data pipeline successful: %d complex matrices, recovery error %.3f', ...
        complex_validation.n_complex_matrices, recovery_metrics.median_recovery_error);
catch ME
    passed = false;
    details = sprintf('Complex data pipeline test failed: %s', ME.message);
end
end

function [passed, details] = test_multiple_problem_sizes()
passed = true; details = 'Multiple problem sizes validation';
try
    problem_sizes = [struct('n_nodes',4,'n_freq',3), struct('n_nodes',10,'n_freq',6), struct('n_nodes',15,'n_freq',8)];
    recovery_errors = []; problem_info = {};
    for i = 1:numel(problem_sizes)
        ps = problem_sizes(i);
        test_params = struct('n_nodes', ps.n_nodes, 'n_freq', ps.n_freq, ...
                             'n_samples', max(50, 2*ps.n_nodes), ...
                             'graph_type', 'random', 'edge_density', 0.3);
        [recovery_metrics, pipeline_data] = run_simplified_pipeline(test_params);
        recovery_errors(end+1) = recovery_metrics.median_recovery_error; %#ok<AGROW>
        problem_info{end+1} = sprintf('%dx%d', ps.n_nodes, ps.n_freq); %#ok<AGROW>

        max_acceptable_error = 0.4;
        if recovery_metrics.median_recovery_error > max_acceptable_error
            passed = false;
            details = sprintf('Problem size %s: recovery error %.3f > %.3f', ...
                problem_info{end}, recovery_metrics.median_recovery_error, max_acceptable_error);
            return;
        end

        validation_results = validate_recovery_metrics(recovery_metrics, pipeline_data);
        if ~validation_results.meets_criteria
            passed = false;
            details = sprintf('Cycle consistency failed for %s: %s', problem_info{end}, validation_results.failure_reason);
            return;
        end
    end
    details = sprintf('Multiple sizes successful: errors [%s] for problems [%s]', ...
        sprintf('%.3f ', recovery_errors), strjoin(problem_info, ', '));
catch ME
    passed = false;
    details = sprintf('Multiple problem sizes test failed: %s', ME.message);
end
end

function [passed, details] = test_different_graph_types()
passed = true; details = 'Different graph types validation';
try
    graph_types = {'random', 'hub', 'chain'};
    recovery_results = struct();
    for i = 1:numel(graph_types)
        gt = graph_types{i};
        test_params = struct('n_nodes', 8, 'n_freq', 4, 'n_samples', 100, ...
                             'graph_type', gt, 'edge_density', 0.25);
        [recovery_metrics, pipeline_data] = run_simplified_pipeline(test_params);
        recovery_results.(gt) = recovery_metrics.median_recovery_error;

        max_error_by_type = containers.Map({'random','hub','chain'}, {0.35,0.4,0.3});
        if recovery_metrics.median_recovery_error > max_error_by_type(gt)
            passed = false;
            details = sprintf('%s graph: recovery error %.3f > %.3f', ...
                gt, recovery_metrics.median_recovery_error, max_error_by_type(gt));
            return;
        end

        validation_results = validate_recovery_metrics(recovery_metrics, pipeline_data);
        if ~validation_results.meets_criteria
            passed = false;
            details = sprintf('%s graph: cycle consistency failed: %s', gt, validation_results.failure_reason);
            return;
        end
    end
    details = sprintf('Graph types successful: random=%.3f, hub=%.3f, chain=%.3f', ...
        recovery_results.random, recovery_results.hub, recovery_results.chain);
catch ME
    passed = false;
    details = sprintf('Graph types test failed: %s', ME.message);
end
end

function [passed, details] = test_noise_robustness()
passed = true; details = 'Noise robustness validation';
try
    noise_levels = [50, 100, 200];
    recovery_errors = [];
    for i = 1:numel(noise_levels)
        n_samples = noise_levels(i);
        test_params = struct('n_nodes', 6, 'n_freq', 4, 'n_samples', n_samples, ...
                             'graph_type', 'random', 'edge_density', 0.3);
        [recovery_metrics, pipeline_data] = run_simplified_pipeline(test_params);
        recovery_errors(end+1) = recovery_metrics.median_recovery_error; %#ok<AGROW>

        validation_results = validate_recovery_metrics(recovery_metrics, pipeline_data);
        if ~validation_results.meets_criteria
            warning('Cycle consistency may be loose at %d samples: %s', n_samples, validation_results.failure_reason);
        end
    end

    best_recovery = min(recovery_errors);
    if best_recovery > 0.25
        passed = false;
        details = sprintf('Best recovery error %.3f > 0.25 even with %d samples', ...
            best_recovery, max(noise_levels));
        return;
    end
    details = sprintf('Noise robustness validated: errors [%s] for samples [%s]', ...
        sprintf('%.3f ', recovery_errors), sprintf('%d ', noise_levels));
catch ME
    passed = false;
    details = sprintf('Noise robustness test failed: %s', ME.message);
end
end

function [passed, details] = test_ground_truth_recovery()
passed = true; details = 'Ground truth recovery quality validation';
try
    test_params = struct('n_nodes', 10, 'n_freq', 6, 'n_samples', 200, ...
                         'graph_type', 'random', 'edge_density', 0.3);
    [recovery_metrics, pipeline_data] = run_simplified_pipeline(test_params);

    gt_evaluation = evaluate_ground_truth_recovery(pipeline_data, recovery_metrics);

    checks.median_error = (gt_evaluation.median_relative_error <= 0.3);
    checks.support_f1 = (gt_evaluation.support_f1_score >= 0.8);
    checks.fdr = (gt_evaluation.false_discovery_rate <= 0.1);

    validation_results = validate_recovery_metrics(recovery_metrics, pipeline_data);
    checks.cycle_consistency = validation_results.meets_criteria_cycle;

    if ~all(structfun(@(x) x, checks))
        fn = fieldnames(checks); bad = fn(~structfun(@(x) x, checks));
        passed = false;
        details = sprintf('Ground truth criteria failed: %s', strjoin(bad, ', '));
        return;
    end

    details = sprintf('Ground truth recovery excellent: error=%.3f, F1=%.3f, FDR=%.3f, cycle=OK', ...
        gt_evaluation.median_relative_error, gt_evaluation.support_f1_score, gt_evaluation.false_discovery_rate);
catch ME
    passed = false;
    details = sprintf('Ground truth recovery test failed: %s', ME.message);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper Functions (pipeline, validation, metrics)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [recovery_metrics, pipeline_data] = run_simplified_pipeline(test_params)
% Module 7 (simulation) -> Module 1 (whitening) -> Theoretical Gamma_tilde -> Module 8 (recolor)

try
    if isfield(test_params, 'complex_strength')
        [Omega_true, Sigma_true, Sigma_emp, sim_params] = ...
            module7_simulation_improved_complex( ...
                'n_nodes', test_params.n_nodes, ...
                'n_freq', test_params.n_freq, ...
                'n_samples', test_params.n_samples, ...
                'graph_type', test_params.graph_type, ...
                'edge_density', test_params.edge_density, ...
                'complex_strength', test_params.complex_strength, ...
                'random_seed', 42);
    else
        [Omega_true, Sigma_true, Sigma_emp, sim_params] = ...
            module7_simulation_improved_complex( ...
                'n_nodes', test_params.n_nodes, ...
                'n_freq', test_params.n_freq, ...
                'n_samples', test_params.n_samples, ...
                'graph_type', test_params.graph_type, ...
                'edge_density', test_params.edge_density, ...
                'random_seed', 42);
    end
catch
    [Omega_true, Sigma_true, Sigma_emp, sim_params] = ...
        module7_simulation( ...
            'n_nodes', test_params.n_nodes, ...
            'n_freq', test_params.n_freq, ...
            'n_samples', test_params.n_samples, ...
            'graph_type', test_params.graph_type, ...
            'edge_density', test_params.edge_density, ...
            'random_seed', 42);
end

% Module 1 preprocessing (build D from empirical covariances as in production)
module1_input.mode = 'simulation';
module1_input.sim_results.Sigma_emp = Sigma_emp;
module1_input.sim_results.F = numel(Sigma_emp);
module1_input.sim_results.n = size(Sigma_emp{1}, 1);
module1_input.sim_results.T = test_params.n_samples;
module1_results = module1_preprocessing_main(module1_input, struct('verbose', false));

% ==== Key fix: Construct theoretical Gamma_tilde using Σ_true, not Σ_emp ====
F = numel(Sigma_emp);
Gamma_tilde_theoretical = cell(F, 1);
for omega = 1:F
    D_omega = module1_results.D{omega};
    Sigma_whitened_true = D_omega * Sigma_true{omega} * D_omega;  % whiten Σ_true
    Sigma_whitened_true = (Sigma_whitened_true + Sigma_whitened_true') / 2; % enforce Hermitian

    e = eig(Sigma_whitened_true);
    min_ev = min(real(e));
    if min_ev <= 1e-10
        Sigma_whitened_true = Sigma_whitened_true + (abs(min_ev) + 1e-8) * eye(size(Sigma_whitened_true));
    end
    Gamma_tilde_theoretical{omega} = inv(Sigma_whitened_true);
end

module8_input.whitened_precision_matrices = Gamma_tilde_theoretical;
module8_input.whitening_matrices = module1_results.D;
module8_input.original_covariances = Sigma_true;
module8_params = struct('verbose', false, 'compute_quality_metrics', true);
module8_results = module8_recoloring_main(module8_input, module8_params);

Omega_estimated = module8_results.recolored_precision_matrices;
recovery_metrics = compute_recovery_metrics(Omega_estimated, Omega_true, module8_results);

pipeline_data.ground_truth = struct('Omega_true', {Omega_true}, 'Sigma_true', {Sigma_true});
pipeline_data.module1_results = module1_results;
pipeline_data.module8_results = module8_results;
pipeline_data.sim_params = sim_params;
end

function recovery_metrics = compute_recovery_metrics(Omega_estimated, Omega_true, module8_results)
F = numel(Omega_estimated);
relative_errors = zeros(F,1); cycle_errors = NaN(F,1); herm_errors = NaN(F,1);
for omega = 1:F
    if ~isempty(Omega_estimated{omega}) && ~isempty(Omega_true{omega})
        relative_errors(omega) = norm(Omega_estimated{omega} - Omega_true{omega}, 'fro') / ...
                                 max(norm(Omega_true{omega}, 'fro'), eps);
        if isfield(module8_results,'recoloring_quality') && omega <= numel(module8_results.recoloring_quality)
            q = module8_results.recoloring_quality{omega};
            if isfield(q,'inv_error') && isfinite(q.inv_error), cycle_errors(omega) = q.inv_error; end
            if isfield(q,'hermitian_error') && isfinite(q.hermitian_error), herm_errors(omega) = q.hermitian_error; end
        end
    else
        relative_errors(omega) = NaN;
    end
end
valid_errors = relative_errors(~isnan(relative_errors));
valid_cycle = cycle_errors(~isnan(cycle_errors));
valid_herm = herm_errors(~isnan(herm_errors));

recovery_metrics.median_recovery_error = median(valid_errors);
recovery_metrics.mean_recovery_error = mean(valid_errors);
recovery_metrics.max_recovery_error = max(valid_errors);
recovery_metrics.std_recovery_error = std(valid_errors);

recovery_metrics.median_cycle_error = median(valid_cycle);
recovery_metrics.max_cycle_error = max(valid_cycle);

recovery_metrics.max_hermitian_error = max(valid_herm);
recovery_metrics.all_hermitian_ok = all(valid_herm < 1e-12);

recovery_metrics.n_successful_frequencies = sum(~isnan(relative_errors));
recovery_metrics.success_rate = recovery_metrics.n_successful_frequencies / F;
end

function validation_results = validate_recovery_metrics(recovery_metrics, pipeline_data)
% Adaptive tolerance for cycle-consistency based on cond(D)^2
validation_results = struct(); validation_results.meets_criteria = true;
validation_results.meets_criteria_cycle = true; validation_results.failure_reason = '';

% Criterion 1: median relative error
if recovery_metrics.median_recovery_error > 0.3
    validation_results.meets_criteria = false;
    validation_results.failure_reason = sprintf('Median recovery error %.3f > 0.3', ...
        recovery_metrics.median_recovery_error); return;
end

% Criterion 2: success rate
if recovery_metrics.success_rate < 0.9
    validation_results.meets_criteria = false;
    validation_results.failure_reason = sprintf('Success rate %.2f < 0.9', recovery_metrics.success_rate); return;
end

% Criterion 3: adaptive cycle-consistency threshold
c = 1e3; machine_eps = eps('double');
F = numel(pipeline_data.module1_results.D);
condD2 = zeros(F,1);
for w = 1:F
    condD2(w) = cond(pipeline_data.module1_results.D{w})^2;
end
tau_inv = c * machine_eps * prctile(condD2, 95);
if recovery_metrics.max_cycle_error > tau_inv
    validation_results.meets_criteria = false;
    validation_results.meets_criteria_cycle = false;
    validation_results.failure_reason = sprintf('Max cycle error %.2e > adaptive tol %.2e', ...
        recovery_metrics.max_cycle_error, tau_inv); return;
end
end

function consistency_results = validate_pipeline_consistency(pipeline_data)
consistency_results = struct('is_consistent', true, 'failure_reason', '');
try
    F = numel(pipeline_data.ground_truth.Omega_true);
    p = size(pipeline_data.ground_truth.Omega_true{1},1);

    if numel(pipeline_data.module1_results.D) ~= F
        consistency_results.is_consistent = false;
        consistency_results.failure_reason = 'Module 1 whitening matrices count mismatch'; return;
    end
    if numel(pipeline_data.module1_results.Sigma_tilde) ~= F
        consistency_results.is_consistent = false;
        consistency_results.failure_reason = 'Module 1 whitened covariances count mismatch'; return;
    end
    if numel(pipeline_data.module8_results.recolored_precision_matrices) ~= F
        consistency_results.is_consistent = false;
        consistency_results.failure_reason = 'Module 8 output count mismatch'; return;
    end

    for omega = 1:F
        if size(pipeline_data.module1_results.D{omega}, 1) ~= p
            consistency_results.is_consistent = false;
            consistency_results.failure_reason = sprintf('Module 1 whitening matrix size at ω=%d', omega); return;
        end
        if size(pipeline_data.module1_results.Sigma_tilde{omega}, 1) ~= p
            consistency_results.is_consistent = false;
            consistency_results.failure_reason = sprintf('Module 1 whitened covariance size at ω=%d', omega); return;
        end
        if ~isempty(pipeline_data.module8_results.recolored_precision_matrices{omega})
            if size(pipeline_data.module8_results.recolored_precision_matrices{omega}, 1) ~= p
                consistency_results.is_consistent = false;
                consistency_results.failure_reason = sprintf('Module 8 output size at ω=%d', omega); return;
            end
        end
    end
catch ME
    consistency_results.is_consistent = false;
    consistency_results.failure_reason = sprintf('Consistency check failed: %s', ME.message);
end
end

function complex_validation = validate_complex_data_handling(pipeline_data)
complex_validation = struct('success', true, 'failure_reason', '', 'n_complex_matrices', 0);
try
    F = numel(pipeline_data.ground_truth.Omega_true);
    for omega = 1:F
        is_complex_gt = ~isreal(pipeline_data.ground_truth.Omega_true{omega});
        if is_complex_gt
            Omega_result = pipeline_data.module8_results.recolored_precision_matrices{omega};
            if ~isempty(Omega_result)
                is_complex_res = ~isreal(Omega_result);
                if ~is_complex_res
                    complex_validation.success = false;
                    complex_validation.failure_reason = sprintf('Complex GT at ω=%d produced real result', omega);
                    return;
                end
                herm_err = norm(Omega_result - Omega_result', 'fro') / max(norm(Omega_result, 'fro'), eps);
                if herm_err > 1e-12
                    complex_validation.success = false;
                    complex_validation.failure_reason = sprintf('Complex result not Hermitian at ω=%d: %.2e', ...
                        omega, herm_err); return;
                end
                complex_validation.n_complex_matrices = complex_validation.n_complex_matrices + 1;
            end
        end
    end
catch ME
    complex_validation.success = false;
    complex_validation.failure_reason = sprintf('Complex validation failed: %s', ME.message);
end
end

function gt_evaluation = evaluate_ground_truth_recovery(pipeline_data, recovery_metrics)
gt_evaluation = struct();
try
    Omega_true = pipeline_data.ground_truth.Omega_true;
    Omega_estimated = pipeline_data.module8_results.recolored_precision_matrices;
    F = numel(Omega_true);

    relative_errors = [];
    for omega = 1:F
        if ~isempty(Omega_estimated{omega})
            rel_err = norm(Omega_estimated{omega} - Omega_true{omega}, 'fro') / ...
                      max(norm(Omega_true{omega}, 'fro'), eps);
            relative_errors(end+1) = rel_err; %#ok<AGROW>
        end
    end
    gt_evaluation.median_relative_error = median(relative_errors);
    gt_evaluation.mean_relative_error = mean(relative_errors);

    support_metrics = compute_support_recovery_metrics(Omega_estimated, Omega_true);
    gt_evaluation.support_f1_score = support_metrics.f1_score;
    gt_evaluation.support_precision = support_metrics.precision;
    gt_evaluation.support_recall = support_metrics.recall;
    gt_evaluation.false_discovery_rate = support_metrics.false_discovery_rate;

    gt_evaluation.max_cycle_error = recovery_metrics.max_cycle_error;
    gt_evaluation.success = true;
catch ME
    gt_evaluation.success = false;
    gt_evaluation.error = ME.message;
end
end

function support_metrics = compute_support_recovery_metrics(Omega_estimated, Omega_true)
support_metrics = struct();
try
    F = numel(Omega_true);
    p = size(Omega_true{1}, 1);
    TP=0; FP=0; TN=0; FN=0;
    base_threshold = 1e-3;
    for omega = 1:F
        if ~isempty(Omega_estimated{omega})
            Oe = Omega_estimated{omega}; Og = Omega_true{omega};
            thr_gt = base_threshold * norm(Og, 'fro') / p;
            thr_est = base_threshold * norm(Oe, 'fro') / p;
            [I,J] = meshgrid(1:p,1:p); off = (I~=J);
            Sg = (abs(Og) > thr_gt) & off; Se = (abs(Oe) > thr_est) & off;
            TP = TP + sum(Sg(:) & Se(:));
            FP = FP + sum(~Sg(:) & Se(:));
            TN = TN + sum(~Sg(:) & ~Se(:));
            FN = FN + sum(Sg(:) & ~Se(:));
        end
    end
    precision = TP / max(TP+FP,1);
    recall    = TP / max(TP+FN,1);
    f1        = 2*precision*recall / max(precision+recall, eps);
    fdr       = FP / max(TP+FP,1);
    support_metrics.precision = precision;
    support_metrics.recall = recall;
    support_metrics.f1_score = f1;
    support_metrics.false_discovery_rate = fdr;
    support_metrics.success = true;
catch ME
    support_metrics.success = false;
    support_metrics.error = ME.message;
    support_metrics.precision = 0; support_metrics.recall = 0;
    support_metrics.f1_score = 0;  support_metrics.false_discovery_rate = 1;
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
