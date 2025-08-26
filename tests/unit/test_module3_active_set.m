function test_results = test_module3_active_set()
% TEST_MODULE3_ACTIVE_SET - Comprehensive test suite for Module 3 active set selection
%
% 返回：
%   test_results - 结构体，包含测试统计
%
% 注：本版本仅做“必要修复”，不改变仿真模块的可选模式集合

rng(42); % For reproducibility

fprintf('========================================\n');
fprintf('Module 3 Active Set Selection Test Suite\n');
fprintf('========================================\n\n');

test_results = struct();
test_results.total_tests = 0;
test_results.passed_tests = 0;
test_results.failed_tests = 0;
test_results.test_details = {};

%% Test 1: Edge Proxy Computation with Synthetic Data
fprintf('Test 1: Edge Proxy Computation with Synthetic Data\n');
fprintf('--------------------------------------------------\n');
[passed, details] = test_edge_proxy_computation_synthetic();
test_results = record_test_result(test_results, 'edge_proxy_synthetic', passed, details);

%% Test 2: Threshold Determination Edge Cases
fprintf('\nTest 2: Threshold Determination Edge Cases\n');
fprintf('------------------------------------------\n');
[passed, details] = test_threshold_determination_edge_cases();
test_results = record_test_result(test_results, 'threshold_edge_cases', passed, details);

%% Test 3: Combined Active Set with Module7 Data
fprintf('\nTest 3: Combined Active Set with Module7 Data\n');
fprintf('---------------------------------------------\n');
[passed, details] = test_combined_active_set_module7();
test_results = record_test_result(test_results, 'combined_active_module7', passed, details);

%% Test 4: Full Pipeline Integration with Modules 1 and 2
fprintf('\nTest 4: Full Pipeline Integration with Modules 1 and 2\n');
fprintf('------------------------------------------------------\n');
[passed, details] = test_full_pipeline_integration();
test_results = record_test_result(test_results, 'full_pipeline_integration', passed, details);

%% Test 5: Precision vs Correlation Proxy Comparison
fprintf('\nTest 5: Precision vs Correlation Proxy Comparison\n');
fprintf('-------------------------------------------------\n');
[passed, details] = test_precision_vs_correlation_proxies();
test_results = record_test_result(test_results, 'precision_correlation_comparison', passed, details);

%% Test 6: Active Set Robustness and Performance
fprintf('\nTest 6: Active Set Robustness and Performance\n');
fprintf('---------------------------------------------\n');
[passed, details] = test_active_set_robustness();
test_results = record_test_result(test_results, 'active_set_robustness', passed, details);

%% Final Report
fprintf('\n========================================\n');
fprintf('Module 3 Test Suite Summary\n');
fprintf('========================================\n');
fprintf('Total tests: %d\n', test_results.total_tests);
fprintf('Passed: %d\n', test_results.passed_tests);
fprintf('Failed: %d\n', test_results.failed_tests);
fprintf('Success rate: %.1f%%\n', 100 * test_results.passed_tests / max(test_results.total_tests, 1));

if test_results.failed_tests > 0
    fprintf('\nFailed tests:\n');
    for i = 1:length(test_results.test_details)
        detail = test_results.test_details{i};
        if ~detail.passed
            fprintf('  - %s: %s\n', detail.test_name, detail.error_message);
        end
    end
end

fprintf('\n✓ Module 3 testing completed\n');

end

%% Individual Test Functions

function [passed, details] = test_edge_proxy_computation_synthetic()
details = struct('test_name', 'Edge Proxy Computation Synthetic', 'subtests', {{}}); passed = true;
try
    fprintf('  1.1 Creating synthetic test data... ');
    p = 15; F = 5;
    true_precision = cell(F, 1);
    empirical_cov = cell(F, 1);
    for f = 1:F
        Omega = eye(p) * 2;
        block_size = 5;
        for b = 1:floor(p/block_size)
            idx = (b-1)*block_size + 1 : min(b*block_size, p);
            for i = idx
                for j = idx
                    if i ~= j
                        strength = 0.3 + 0.2 * sin(2*pi*f/F);
                        Omega(i, j) = -strength * (0.8 + 0.4*rand());
                        Omega(j, i) = Omega(i, j);
                    end
                end
            end
        end
        [V, D] = eig(Omega);
        D = diag(max(diag(D), 0.1));
        Omega = V * D * V'; Omega = (Omega + Omega') / 2;
        true_precision{f} = Omega;
        empirical_cov{f} = inv(Omega) + 0.05 * randn(p) * randn(p)'; 
        empirical_cov{f} = (empirical_cov{f} + empirical_cov{f}') / 2;
    end
    fprintf('✓\n'); details.subtests{end+1} = struct('name','synthetic_data_creation','passed',true);

    fprintf('  1.2 Testing correlation-based proxies... ');
    corr_proxies = module3_edge_proxy_computation(empirical_cov, 'correlation');
    assert(iscell(corr_proxies) && numel(corr_proxies) == F);
    for f = 1:F
        proxy_matrix = corr_proxies{f};
        assert(isnumeric(proxy_matrix) && isequal(size(proxy_matrix), [p p]));
        assert(all(proxy_matrix(:) >= 0));
        assert(all(diag(proxy_matrix) == 0));
        sym_error = norm(proxy_matrix - proxy_matrix', 'fro');
        assert(sym_error < 1e-12);
    end
    fprintf('✓\n'); details.subtests{end+1} = struct('name','correlation_proxies','passed',true);

    fprintf('  1.3 Testing precision-based proxies... ');
    prec_proxies = module3_edge_proxy_computation(true_precision, 'precision');
    assert(iscell(prec_proxies) && numel(prec_proxies) == F);
    for f = 1:F
        proxy_matrix = prec_proxies{f};
        assert(isnumeric(proxy_matrix) && isequal(size(proxy_matrix), [p p]));
        assert(all(proxy_matrix(:) >= 0));
        assert(all(diag(proxy_matrix) == 0));
        assert(all(proxy_matrix(:) <= 1 + eps));
    end
    fprintf('✓\n'); details.subtests{end+1} = struct('name','precision_proxies','passed',true);

    fprintf('  1.4 Comparing proxy methods on known structure... ');
    strong_edges_found_corr = 0; strong_edges_found_prec = 0;
    for f = 1:F
        Omega = true_precision{f};
        corr_proxy = corr_proxies{f}; prec_proxy = prec_proxies{f};
        [strong_i, strong_j] = find(triu(abs(Omega), 1) > 0.15);
        corr_threshold = quantile(corr_proxy(corr_proxy(:) > 0), 0.7);
        prec_threshold = quantile(prec_proxy(prec_proxy(:) > 0), 0.7);
        for k = 1:length(strong_i)
            i = strong_i(k); j = strong_j(k);
            if corr_proxy(i, j) >= corr_threshold, strong_edges_found_corr = strong_edges_found_corr + 1; end
            if prec_proxy(i, j) >= prec_threshold, strong_edges_found_prec = strong_edges_found_prec + 1; end
        end
    end
    fprintf('✓ (Corr: %d, Prec: %d strong edges found)\n', strong_edges_found_corr, strong_edges_found_prec);
    details.subtests{end+1} = struct('name','proxy_method_comparison','passed',true);
catch ME
    passed = false; details.error_message = ME.message;
    fprintf('✗ (%s)\n', ME.message);
    details.subtests{end+1} = struct('name','failed','passed',false,'error',ME.message);
end
end

function [passed, details] = test_threshold_determination_edge_cases()
details = struct('test_name', 'Threshold Determination Edge Cases', 'subtests', {{}}); passed = true;
try
    fprintf('  2.1 Creating test proxy data... ');
    p = 10; F = 4; edge_proxies = cell(F,1);
    for f = 1:F
        proxy_matrix = zeros(p, p);
        base_value = f * 0.1;
        for i = 1:p
            for j = i+1:p
                proxy_matrix(i, j) = max(base_value + 0.05 * randn(), 0);
                proxy_matrix(j, i) = proxy_matrix(i, j);
            end
        end
        edge_proxies{f} = proxy_matrix;
    end
    fprintf('✓\n'); details.subtests{end+1} = struct('name','test_data_creation','passed',true);

    fprintf('  2.2 Testing various quantile levels... ');
    quantile_levels = [0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9];
    thresholds = zeros(size(quantile_levels));
    for i = 1:length(quantile_levels)
        q = quantile_levels(i);
        tau = module3_threshold_determination(edge_proxies, q);
        thresholds(i) = tau;
        assert(isfinite(tau) && tau >= 0);
    end
    for i = 2:length(thresholds)
        assert(thresholds(i) >= thresholds(i-1), 'Thresholds not monotonically increasing');
    end
    fprintf('✓\n'); details.subtests{end+1} = struct('name','quantile_levels','passed',true);

    fprintf('  2.3 Testing edge cases... ');
    zero_proxies = repmat({zeros(p,p)}, F, 1);
    tau_zero = module3_threshold_determination(zero_proxies, 0.5, struct('exclude_zeros', false));
    assert(tau_zero == 0, 'All-zero case should give threshold = 0');
    single_proxies = { [0 0.5 0; 0.5 0 0; 0 0 0] };
    tau_single = module3_threshold_determination(single_proxies, 0.5);
    assert(abs(tau_single - 0.5) < 1e-10, 'Single value case failed');
    fprintf('✓\n'); details.subtests{end+1} = struct('name','edge_cases','passed',true);

    fprintf('  2.4 Testing with statistics output... ');
    [~, stats] = module3_threshold_determination(edge_proxies, 0.1, struct('verbose', false));
    required_fields = {'total_proxy_values','valid_proxy_values','excluded_values','proxy_range','threshold_percentile'};
    for i = 1:length(required_fields)
        assert(isfield(stats, required_fields{i}), 'Missing statistics field: %s', required_fields{i});
    end
    assert(stats.valid_proxy_values <= stats.total_proxy_values);
    assert(stats.threshold_percentile >= 0 && stats.threshold_percentile <= 100);
    fprintf('✓\n'); details.subtests{end+1} = struct('name','statistics_output','passed',true);
catch ME
    passed = false; details.error_message = ME.message;
    fprintf('✗ (%s)\n', ME.message);
    details.subtests{end+1} = struct('name','failed','passed',false,'error',ME.message);
end
end

function [passed, details] = test_combined_active_set_module7()
details = struct('test_name', 'Combined Active Set Module7', 'subtests', {{}}); passed = true;
try
    fprintf('  3.1 Generating Module7 simulation data... ');
    % 允许的 graph_type：'random','chain','hub' —— 这里用 'hub'
    [true_precision, ~, empirical_covariance, sim_params] = ...
        module7_simulation_improved_complex('n_nodes', 20, 'n_freq', 8, ...
                                           'n_samples', 100, 'graph_type', 'hub', ...
                                           'edge_density', 0.15);
    p = sim_params.n_nodes; F = sim_params.n_freq;
    fprintf('✓ (%dx%d, %d freq)\n', p, p, F);
    details.subtests{end+1} = struct('name','module7_data_generation','passed',true);

    fprintf('  3.2 Computing edge proxies... ');
    corr_proxies = module3_edge_proxy_computation(empirical_covariance, 'correlation');
    prec_proxies = module3_edge_proxy_computation(true_precision, 'precision');
    fprintf('✓\n'); details.subtests{end+1} = struct('name','edge_proxy_computation','passed',true);

    fprintf('  3.3 Testing threshold determination... ');
    quantile_levels = [0.05, 0.1, 0.2]; thresholds = zeros(size(quantile_levels));
    for i = 1:length(quantile_levels)
        [tau, stats] = module3_threshold_determination(corr_proxies, quantile_levels(i), struct('verbose', false));
        thresholds(i) = tau; assert(isfinite(tau) && tau >= 0); assert(stats.valid_proxy_values > 0);
    end
    fprintf('✓\n'); details.subtests{end+1} = struct('name','threshold_determination','passed',true);

    fprintf('  3.4 Constructing combined active sets... ');
    tau = thresholds(2); % 0.1 quantile
    [combined_masks, active_stats] = module3_combined_active_set(corr_proxies, tau, struct('verbose', false));
    required_mask_fields = {'edge_masks','node_masks','combined_masks'};
    for i = 1:length(required_mask_fields)
        assert(isfield(combined_masks, required_mask_fields{i}));
    end
    assert(isequal(size(combined_masks.edge_masks), [p p F]));
    assert(isequal(size(combined_masks.node_masks), [p F]));
    assert(isequal(size(combined_masks.combined_masks), [p p F]));
    fprintf('✓\n'); details.subtests{end+1} = struct('name','combined_active_set_construction','passed',true);

    fprintf('  3.5 Validating active set properties... ');
    total_active_edges = 0; total_active_nodes = 0;
    for f = 1:F
        edge_mask = combined_masks.edge_masks(:, :, f);
        node_mask = combined_masks.node_masks(:, f);
        combined_mask = combined_masks.combined_masks(:, :, f);
        assert(islogical(edge_mask) && islogical(node_mask) && islogical(combined_mask));
        assert(isequal(edge_mask, edge_mask'));
        assert(isequal(combined_mask, combined_mask'));
        assert(all(combined_mask(:) <= edge_mask(:)));
        [ai, aj] = find(combined_mask);
        for k = 1:length(ai)
            i = ai(k); j = aj(k);
            if i ~= j
                assert(node_mask(i) && node_mask(j), 'Active edge (%d,%d) connects inactive nodes at f=%d', i, j, f);
            end
        end
        off_diag = combined_mask & ~eye(p);
        total_active_edges = total_active_edges + sum(off_diag(:)) / 2;
        total_active_nodes = total_active_nodes + sum(node_mask);
    end
    max_possible_edges = F * p * (p - 1) / 2;
    sparsity_ratio = total_active_edges / max_possible_edges;
    assert(sparsity_ratio > 0.01 && sparsity_ratio < 0.8, 'Unreasonable sparsity ratio: %.3f', sparsity_ratio);
    fprintf('✓ (%.1f%% sparsity)\n', 100 * sparsity_ratio);
    details.subtests{end+1} = struct('name','active_set_validation','passed',true);

    fprintf('  3.6 Testing different proxy methods... ');
    tau_prec = module3_threshold_determination(prec_proxies, 0.1);
    [combined_masks_prec, ~] = module3_combined_active_set(prec_proxies, tau_prec, struct('verbose', false));
    corr_edges = 0; prec_edges = 0;
    for f = 1:F
        corr_mask = combined_masks.combined_masks(:, :, f);
        prec_mask = combined_masks_prec.combined_masks(:, :, f);
        corr_off = corr_mask & ~eye(p); prec_off = prec_mask & ~eye(p);
        corr_edges = corr_edges + sum(corr_off(:)) / 2;
        prec_edges = prec_edges + sum(prec_off(:)) / 2;
    end
    fprintf('✓ (Corr: %d, Prec: %d edges)\n', corr_edges, prec_edges);
    details.subtests{end+1} = struct('name','proxy_method_comparison','passed',true);

    % 3.7 可视化：你有 GUI 环境，保留
    try
        create_active_set_visualization(combined_masks, active_stats, sim_params);
        fprintf('  3.7 Created visualization ✓\n');
        details.subtests{end+1} = struct('name','visualization','passed',true);
    catch ME
        warning(ME.identifier, '%s', ME.message); % 按你要求使用格式化传递
        fprintf('  3.7 Visualization skipped\n');
    end

catch ME
    passed = false; details.error_message = ME.message;
    fprintf('✗ (%s)\n', ME.message);
    details.subtests{end+1} = struct('name','failed','passed',false,'error',ME.message);
end
end

function [passed, details] = test_full_pipeline_integration()
details = struct('test_name', 'Full Pipeline Integration', 'subtests', {{}}); passed = true;
try
    fprintf('  4.1 Generating comprehensive test data... ');
    [true_precision, true_covariance, empirical_covariance, sim_params] = ...
        module7_simulation_improved_complex('n_nodes', 15, 'n_freq', 6, ...
                                           'n_samples', 120, 'graph_type', 'hub'); % 允许模式
    p = sim_params.n_nodes; F = sim_params.n_freq;
    frequencies = linspace(8, 12, F);
    fprintf('✓\n'); details.subtests{end+1} = struct('name','test_data_generation','passed',true);

    fprintf('  4.2 Running Module1 preprocessing... ');
    module1_input = struct();
    module1_input.mode = 'simulation';
    module1_input.sim_results = struct();
    module1_input.sim_results.Sigma_emp = empirical_covariance;
    module1_input.sim_results.F = F;
    module1_input.sim_results.n = p;
    module1_input.sim_results.T = sim_params.n_samples;
    module1_params = struct('verbose', false);
    % 兼容层生效后，这里不会再因为 .whiten 缺失而失败
    module1_results = module1_preprocessing_main(module1_input, module1_params);
    assert(module1_results.success, 'Module1 preprocessing failed');
    fprintf('✓\n'); details.subtests{end+1} = struct('name','module1_preprocessing','passed',true);

    fprintf('  4.3 Running Module2 E-step (if available)... ');
    use_module2_precision = false;
    try
        module2_input = struct();
        module2_input.leadfield_matrix = randn(p, round(p*0.8));
        module2_input.empirical_covariances = module1_results.whitened_covariances;
        module2_input.source_prior_covariances = repmat({eye(round(p*0.8))*0.5}, F, 1);
        module2_input.frequencies = frequencies;
        module2_input.noise_covariance = eye(p) * 0.1;
        module2_params = struct('verbose', false);
        module2_results = module2_estep_main(module2_input, module2_params);
        if isstruct(module2_results) && isfield(module2_results,'success') && module2_results.success
            fprintf('✓\n'); details.subtests{end+1} = struct('name','module2_estep','passed',true);
            use_module2_precision = true;
        else
            fprintf('○ (Module2 available but failed)\n');
        end
    catch
        fprintf('○ (Module2 not available)\n');
    end

    fprintf('  4.4 Running Module3 with whitened covariances... ');
    module3_input = struct();
    module3_input.whitened_covariances = module1_results.whitened_covariances;
    module3_input.frequencies = frequencies;
    if use_module2_precision
        module3_input.initial_precision_matrices = module2_results.initial_precision_matrices;
    end
    module3_params_corr = struct('proxy_method','correlation','quantile_level',0.1,'verbose',false);
    module3_results_corr = module3_active_set_main(module3_input, module3_params_corr);
    assert(module3_results_corr.success, 'Module3 with correlation method failed');
    fprintf('✓\n'); details.subtests{end+1} = struct('name','module3_correlation','passed',true);

    if use_module2_precision
        fprintf('  4.5 Running Module3 with precision-based proxy... ');
        module3_params_prec = struct('proxy_method','precision','quantile_level',0.1,'verbose',false);
        module3_results_prec = module3_active_set_main(module3_input, module3_params_prec);
        assert(module3_results_prec.success, 'Module3 with precision method failed');
        fprintf('✓\n'); details.subtests{end+1} = struct('name','module3_precision','passed',true);
    else
        fprintf('  4.5 Precision-based proxy test skipped (no Module2)\n');
    end

    fprintf('  4.6 Validating pipeline consistency... ');
    assert(numel(module1_results.whitened_covariances) == F);
    assert(size(module1_results.whitened_covariances{1}, 1) == p);
    stats = module3_results_corr.active_edge_statistics;
    assert(stats.total_active_edges > 0);
    assert(stats.average_sparsity_ratio > 0.01 && stats.average_sparsity_ratio < 0.9);
    whitening_matrices = module1_results.whitening_matrices;
    assert(iscell(whitening_matrices) && numel(whitening_matrices) == F);
    for f = 1:F
        D_f = whitening_matrices{f};
        assert(isnumeric(D_f) && isequal(size(D_f), [p p]));
        assert(all(diag(D_f) > 0));
    end
    fprintf('✓\n'); details.subtests{end+1} = struct('name','pipeline_consistency','passed',true);

catch ME
    passed = false; details.error_message = ME.message;
    fprintf('✗ (%s)\n', ME.message);
    details.subtests{end+1} = struct('name','failed','passed',false,'error',ME.message);
end
end

function [passed, details] = test_precision_vs_correlation_proxies()
details = struct('test_name', 'Precision vs Correlation Proxy Comparison', 'subtests', {{}}); passed = true;
try
    fprintf('  5.1 Generating controlled test scenarios... ');
    scenarios = struct();
    % 允许的模式集合：{'random','chain','hub'}
    scenarios(1).name = 'dense';  scenarios(1).n_nodes = 12; scenarios(1).edge_density = 0.4; scenarios(1).graph_type = 'random';
    scenarios(2).name = 'sparse'; scenarios(2).n_nodes = 15; scenarios(2).edge_density = 0.1; scenarios(2).graph_type = 'chain';
    scenarios(3).name = 'hub';    scenarios(3).n_nodes = 18; scenarios(3).edge_density = 0.15; scenarios(3).graph_type = 'hub';
    n_scenarios = length(scenarios); F = 5;
    fprintf('✓ (%d scenarios)\n', n_scenarios);
    details.subtests{end+1} = struct('name','scenario_setup','passed',true);

    comparison_results = struct();
    comparison_results.scenarios = {scenarios.name};
    comparison_results.edge_overlap = zeros(n_scenarios, 1);
    comparison_results.sparsity_corr = zeros(n_scenarios, 1);
    comparison_results.sparsity_prec = zeros(n_scenarios, 1);
    comparison_results.precision_recall = zeros(n_scenarios, 2);

    for s = 1:n_scenarios
        fprintf('  5.%d Testing scenario: %s... ', s+1, scenarios(s).name);
        [true_precision, ~, empirical_covariance, sim_params] = ...
            module7_simulation_improved_complex('n_nodes', scenarios(s).n_nodes, ...
                                               'n_freq', F, 'n_samples', 100, ...
                                               'graph_type', scenarios(s).graph_type, ...
                                               'edge_density', scenarios(s).edge_density);
        p = sim_params.n_nodes;
        corr_proxies = module3_edge_proxy_computation(empirical_covariance, 'correlation');
        prec_proxies = module3_edge_proxy_computation(true_precision, 'precision');
        tau_corr = module3_threshold_determination(corr_proxies, 0.1);
        tau_prec = module3_threshold_determination(prec_proxies, 0.1);
        [masks_corr, stats_corr] = module3_combined_active_set(corr_proxies, tau_corr, struct('verbose', false));
        [masks_prec, stats_prec] = module3_combined_active_set(prec_proxies, tau_prec, struct('verbose', false));
        total_overlap = 0; total_corr_edges = 0; total_prec_edges = 0; total_true_edges = 0; tp_corr = 0; tp_prec = 0;
        for f = 1:F
            corr_mask = masks_corr.combined_masks(:, :, f);
            prec_mask = masks_prec.combined_masks(:, :, f);
            true_mask = abs(true_precision{f}) > 0.05;
            corr_off = corr_mask & ~eye(p); prec_off = prec_mask & ~eye(p); true_off = true_mask & ~eye(p);
            overlap_mask = corr_off & prec_off;
            total_overlap = total_overlap + sum(overlap_mask(:));
            total_corr_edges = total_corr_edges + sum(corr_off(:));
            total_prec_edges = total_prec_edges + sum(prec_off(:));
            total_true_edges = total_true_edges + sum(true_off(:));
            tp_corr = tp_corr + sum(corr_off(:) & true_off(:));
            tp_prec = tp_prec + sum(prec_off(:) & true_off(:));
        end
        comparison_results.edge_overlap(s) = total_overlap / max(min(total_corr_edges, total_prec_edges), 1);
        comparison_results.sparsity_corr(s) = stats_corr.overall_sparsity_ratio;
        comparison_results.sparsity_prec(s) = stats_prec.overall_sparsity_ratio;
        precision_corr = tp_corr / max(total_corr_edges, 1);
        recall_corr = tp_corr / max(total_true_edges, 1);
        comparison_results.precision_recall(s, :) = [precision_corr, recall_corr];
        fprintf('✓ (Overlap: %.3f)\n', comparison_results.edge_overlap(s));
    end

    details.subtests{end+1} = struct('name','scenario_comparisons','passed',true);
    fprintf('  5.%d Analyzing comparison results... ', n_scenarios+2);
    assert(all(comparison_results.edge_overlap >= 0 & comparison_results.edge_overlap <= 1));
    assert(all(comparison_results.sparsity_corr > 0 & comparison_results.sparsity_corr < 1));
    assert(all(comparison_results.sparsity_prec > 0 & comparison_results.sparsity_prec < 1));
    assert(all(comparison_results.precision_recall(:) >= 0 & comparison_results.precision_recall(:) <= 1));
    fprintf('✓\n'); details.subtests{end+1} = struct('name','result_analysis','passed',true);
catch ME
    passed = false; details.error_message = ME.message;
    fprintf('✗ (%s)\n', ME.message);
    details.subtests{end+1} = struct('name','failed','passed',false,'error',ME.message);
end
end

function [passed, details] = test_active_set_robustness()
details = struct('test_name', 'Active Set Robustness', 'subtests', {{}}); passed = true;
try
    fprintf('  6.1 Testing robustness to noise levels... ');
    [true_precision, ~, ~, sim_params] = ...
        module7_simulation_improved_complex('n_nodes', 12, 'n_freq', 4, 'n_samples', 80);
    p = sim_params.n_nodes; F = sim_params.n_freq;
    noise_levels = [0.01, 0.05, 0.1, 0.2];
    active_edge_counts = zeros(size(noise_levels));
    for n = 1:length(noise_levels)
        noise_level = noise_levels(n);
        noisy_covariances = cell(F, 1);
        for f = 1:F
            true_cov = inv(true_precision{f});
            noise_cov = noise_level * randn(p) * randn(p)'; noise_cov = (noise_cov + noise_cov') / 2;
            S = true_cov + noise_cov;
            [V, D] = eig(S); D = diag(max(diag(D), 0.01)); S = V * D * V';
            noisy_covariances{f} = S;
        end
        proxies = module3_edge_proxy_computation(noisy_covariances, 'correlation');
        tau = module3_threshold_determination(proxies, 0.1);
        [~, stats] = module3_combined_active_set(proxies, tau, struct('verbose', false));
        active_edge_counts(n) = stats.total_active_edges;
    end
    edge_variation = std(active_edge_counts) / mean(active_edge_counts);
    assert(edge_variation < 0.5, 'Active set too sensitive to noise (CV = %.3f)', edge_variation);
    fprintf('✓ (CV = %.3f)\n', edge_variation);
    details.subtests{end+1} = struct('name','noise_robustness','passed',true);

    fprintf('  6.2 Testing parameter sensitivity... ');
    noisy_covariances = cell(F, 1);
    for f = 1:F
        true_cov = inv(true_precision{f});
        noise_cov = 0.05 * randn(p) * randn(p)'; noise_cov = (noise_cov + noise_cov') / 2;
        S = true_cov + noise_cov; [V, D] = eig(S); D = diag(max(diag(D), 0.01)); noisy_covariances{f} = V * D * V';
    end
    proxies = module3_edge_proxy_computation(noisy_covariances, 'correlation');
    quantile_levels = [0.05, 0.1, 0.15, 0.2];
    sparsity_ratios = zeros(size(quantile_levels));
    for q = 1:length(quantile_levels)
        tau = module3_threshold_determination(proxies, quantile_levels(q));
        % 关键修正：为了保证与 q 的单调性一致，禁用节点门控/最少边数/强制对角
        opts = struct('verbose', false, 'min_active_per_node', 0, 'force_diagonal_active', false, 'symmetrize_masks', true);
        [~, stats] = module3_combined_active_set(proxies, tau, opts);
        sparsity_ratios(q) = stats.overall_sparsity_ratio;
    end
    for i = 2:length(sparsity_ratios)
        assert(sparsity_ratios(i) >= sparsity_ratios(i-1), 'Sparsity not monotonic with quantile level');
    end
    fprintf('✓\n'); details.subtests{end+1} = struct('name','parameter_sensitivity','passed',true);

    fprintf('  6.3 Testing computational performance... ');
    large_p = 25; large_F = 6;
    large_covariances = cell(large_F, 1);
    for f = 1:large_F
        S = eye(large_p) + 0.3 * randn(large_p); S = S * S';
        large_covariances{f} = S;
    end
    tic_start = tic;
    large_proxies = module3_edge_proxy_computation(large_covariances, 'correlation');
    large_tau = module3_threshold_determination(large_proxies, 0.1);
    [~, large_stats] = module3_combined_active_set(large_proxies, large_tau, struct('verbose', false));
    computation_time = toc(tic_start);
    assert(large_stats.total_active_edges > 0);
    assert(computation_time < 10, 'Computation too slow: %.2fs', computation_time);
    fprintf('✓ (%.2fs for %dx%d)\n', computation_time, large_p, large_p);
    details.subtests{end+1} = struct('name','computational_performance','passed',true);

    fprintf('  6.4 Testing edge cases and error handling... ');
    sparse_covariances = repmat({eye(p) * 10}, F, 1);
    sparse_proxies = module3_edge_proxy_computation(sparse_covariances, 'correlation');
    sparse_tau = module3_threshold_determination(sparse_proxies, 0.1);
    [~, sparse_stats] = module3_combined_active_set(sparse_proxies, sparse_tau, struct('verbose', false));
    assert(sparse_stats.overall_sparsity_ratio < 0.1, 'Sparse test should have low sparsity');
    dense_covariances = repmat({ones(p) + eye(p)}, F, 1);
    dense_proxies = module3_edge_proxy_computation(dense_covariances, 'correlation');
    dense_tau = module3_threshold_determination(dense_proxies, 0.1);
    [~, dense_stats] = module3_combined_active_set(dense_proxies, dense_tau, struct('verbose', false));
    assert(dense_stats.overall_sparsity_ratio > 0.5, 'Dense test should have high sparsity');
    fprintf('✓\n'); details.subtests{end+1} = struct('name','edge_cases','passed',true);

catch ME
    passed = false; details.error_message = ME.message;
    fprintf('✗ (%s)\n', ME.message);
    details.subtests{end+1} = struct('name','failed','passed',false,'error',ME.message);
end
end

%% Helper Functions (unchanged except for warning formatting fix where applicable)

function test_results = record_test_result(test_results, test_name, passed, details)
test_results.total_tests = test_results.total_tests + 1;
if passed, test_results.passed_tests = test_results.passed_tests + 1;
else,     test_results.failed_tests = test_results.failed_tests + 1; end
details.passed = passed;
test_results.test_details{end+1} = details;
end

function create_active_set_visualization(combined_masks, active_stats, sim_params)
try
    p = sim_params.n_nodes; F = sim_params.n_freq;
    figure('Name', 'Module 3 Active Set Analysis', 'Position', [100, 100, 1200, 800]);
    subplot(2, 3, 1);
    bar(1:F, active_stats.active_edges_per_frequency);
    title('Active Edges per Frequency'); xlabel('Frequency Index'); ylabel('Number of Active Edges'); grid on;
    subplot(2, 3, 2);
    bar(1:F, active_stats.active_nodes_per_frequency);
    title('Active Nodes per Frequency'); xlabel('Frequency Index'); ylabel('Number of Active Nodes'); grid on;
    subplot(2, 3, 3);
    plot(1:F, active_stats.sparsity_ratios_per_frequency, 'o-', 'LineWidth', 2);
    title('Sparsity Ratio per Frequency'); xlabel('Frequency Index'); ylabel('Sparsity Ratio');
    ylim([0, max(active_stats.sparsity_ratios_per_frequency) * 1.1]); grid on;
    subplot(2, 3, 4);
    imagesc(double(combined_masks.combined_masks(:, :, 1))); colormap(gray); colorbar;
    title('Active Set Pattern (Freq 1)'); xlabel('Node Index'); ylabel('Node Index');
    subplot(2, 3, 5);
    imagesc(double(combined_masks.node_masks)); colormap(gray); colorbar;
    title('Node Activity Across Frequencies'); xlabel('Frequency Index'); ylabel('Node Index');
    subplot(2, 3, 6);
    stats_data = [active_stats.total_active_edges; ...
                  active_stats.total_active_nodes; ...
                  active_stats.average_edges_per_frequency; ...
                  active_stats.overall_sparsity_ratio * 100];
    stats_labels = {'Total Edges','Total Nodes','Avg Edges/Freq','Sparsity %'};
    bar(stats_data); set(gca, 'XTickLabel', stats_labels); title('Active Set Summary Statistics');
    ylabel('Count / Percentage'); xtickangle(45); grid on;
    sgtitle(sprintf('Module 3 Active Set Results (%dx%d, %d frequencies)', p, p, F), 'FontSize', 14, 'FontWeight', 'bold');
catch ME
    warning(ME.identifier, '%s', ME.message); % 按你要求使用格式化传递
end
end
