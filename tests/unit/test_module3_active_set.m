function test_results = test_module3_active_set()
% TEST_MODULE3_ACTIVE_SET - Tests for Module 3 (Active-set deactivation)
%
% GOAL: Verify that Module 3 deactivates edges/nodes that SHOULD be
% deactivated given the chosen threshold and rules. We DO NOT require
% sparsity to increase if there is nothing to deactivate.
%
% What we check:
%   1) Edge proxy computation sanity (nonnegative, symmetric, zero diagonal)
%   2) Threshold determination sanity and monotonicity across quantiles
%   3) Deactivation correctness: any edge with proxy < tau MUST NOT be active
%   4) If nothing is below tau, the module may keep everything active (OK)
%   5) Monotonicity: with larger quantile (higher tau), active set is
%      non-increasing (i.e., A(q1) ⊇ A(q2) for q1<q2)
%   6) Integration smoke test with Module 1/2 (no dependency on their
%      'success' flags or internals)
%
% Notes:
%   - All comments are in English as requested.
%   - We do NOT assert global sparsity targets or “overall density” fields.
%   - We avoid relying on Module 1/2 custom fields; only use what we need.

rng(42);

fprintf('========================================\n');
fprintf('Module 3 Active-Set Deactivation Test Suite\n');
fprintf('========================================\n\n');

test_results = struct();
test_results.total_tests = 0;
test_results.passed_tests = 0;
test_results.failed_tests = 0;
test_results.test_details = {};

%% Test 1: Edge Proxy Computation (synthetic)
fprintf('Test 1: Edge Proxy Computation (synthetic)\n');
fprintf('------------------------------------------\n');
[passed, details] = t1_edge_proxy_sanity();
test_results = record(test_results, 'edge_proxy_sanity', passed, details);

%% Test 2: Threshold Determination Sanity
fprintf('\nTest 2: Threshold Determination Sanity\n');
fprintf('--------------------------------------\n');
[passed, details] = t2_threshold_sanity();
test_results = record(test_results, 'threshold_sanity', passed, details);

%% Test 3: Deactivation correctness (Module 7 data; correlation proxies)
fprintf('\nTest 3: Deactivation correctness (module7, correlation)\n');
fprintf('-------------------------------------------------------\n');
[passed, details] = t3_deactivation_correctness_corr();
test_results = record(test_results, 'deactivation_corr', passed, details);

%% Test 4: Monotonicity of active set vs. quantile
fprintf('\nTest 4: Monotonicity of active set vs. quantile\n');
fprintf('-----------------------------------------------\n');
[passed, details] = t4_monotonicity_active_set();
test_results = record(test_results, 'active_set_monotonicity', passed, details);

%% Test 5: Integration smoke test with Module 1 & optional Module 2
fprintf('\nTest 5: Integration smoke test (Module 1 & optional Module 2)\n');
fprintf('-------------------------------------------------------------\n');
[passed, details, viz] = t5_integration_smoke();
test_results = record(test_results, 'integration_smoke', passed, details);

% Move the visualization block here (and delete the old broken one):
if passed && ~isempty(viz)
    try
        create_active_set_visualization(viz.combined_masks, viz.active_stats, viz.sim_params);
        fprintf('  ✓ visualization created\n');
    catch ME
        warning(ME.identifier, 'Visualization failed: %s', ME.message);
    end
end
endtest_results = record(test_results, 'integration_smoke', passed, details);

%% Summary
fprintf('\n========================================\n');
fprintf('Summary\n');
fprintf('========================================\n');
fprintf('Total tests: %d\n', test_results.total_tests);
fprintf('Passed: %d\n', test_results.passed_tests);
fprintf('Failed: %d\n', test_results.failed_tests);
fprintf('Success rate: %.1f%%\n', 100 * test_results.passed_tests / max(test_results.total_tests, 1));

if test_results.failed_tests > 0
    fprintf('\nFailed tests:\n');
    for i = 1:numel(test_results.test_details)
        d = test_results.test_details{i};
        if ~d.passed
            fprintf('  - %s: %s\n', d.test_name, d.error_message);
        end
    end
end

 % === Visualization of final active set ===


fprintf('\n✓ Module 3 testing completed\n');

end

% -------------------------------------------------------------------------
% Test 1
% -------------------------------------------------------------------------
function [passed, details] = t1_edge_proxy_sanity()
details = struct('test_name','Edge Proxy Sanity','subtests',{{}}); passed = true;
try
    p = 12; F = 4;
    S = cell(F,1);
    for f = 1:F
        A = randn(p); S{f} = A*A';  % PSD-ish real symmetric
    end
    proxies = module3_edge_proxy_computation(S, 'correlation');
    assert(iscell(proxies) && numel(proxies)==F);
    for f = 1:F
        P = proxies{f};
        assert(isnumeric(P) && isequal(size(P),[p p]));
        assert(all(P(:) >= 0));
        assert(norm(P - P.', 'fro') < 1e-12);
        assert(all(diag(P)==0));
    end
    details.subtests{end+1} = struct('name','correlation_sanity','passed',true);
    fprintf('  ✓ correlation proxies sanity\n');
catch ME
    passed = false; details.error_message = ME.message;
    fprintf('  ✗ %s\n', ME.message);
end
end

% -------------------------------------------------------------------------
% Test 2
% -------------------------------------------------------------------------
function [passed, details] = t2_threshold_sanity()
details = struct('test_name','Threshold Sanity','subtests',{{}}); passed = true;
try
    p=10; F=3; edge_proxies = cell(F,1);
    for f=1:F
        M = triu(abs(randn(p)),1); M = M + M.';   % nonnegative, symmetric, zero diag
        M(1:p+1:end) = 0;
        edge_proxies{f} = M;
    end
    qs = [0.05 0.10 0.25 0.50 0.75 0.90];
    taus = zeros(size(qs));
    for i=1:numel(qs)
        taus(i) = module3_threshold_determination(edge_proxies, qs(i));
        assert(isfinite(taus(i)) && taus(i) >= 0);
    end
    % Monotone nondecreasing with q
    assert(all(diff(taus) >= -1e-12));
    details.subtests{end+1} = struct('name','quantile_monotonicity','passed',true);
    fprintf('  ✓ quantile monotonicity\n');
catch ME
    passed = false; details.error_message = ME.message;
    fprintf('  ✗ %s\n', ME.message);
end
end

% -------------------------------------------------------------------------
% Test 3
% -------------------------------------------------------------------------
function [passed, details] = t3_deactivation_correctness_corr()
details = struct('test_name','Deactivation Correctness (Correlation)','subtests',{{}}); passed = true;
try
    % Generate simulation via Module 7 (any allowed graph_type is fine)
    [~, ~, Sigma_emp, sim_params] = module7_simulation_improved_complex( ...
        'n_nodes', 20, 'n_freq', 6, 'n_samples', 120, 'graph_type', 'hub', 'edge_density', 0.15);
    p = sim_params.n_nodes; F = sim_params.n_freq;

    % Compute proxies and threshold
    proxies = module3_edge_proxy_computation(Sigma_emp, 'correlation');
    tau = module3_threshold_determination(proxies, 0.10);

    % Build combined active set with default options
    [masks, stats] = module3_combined_active_set(proxies, tau, struct('verbose', false));

    % Deactivation rule check:
    % Any off-diagonal edge with proxy < tau MUST NOT appear in combined mask.
    violations = 0; total_checked = 0;
    for f=1:F
        P = proxies{f};
        CM = masks.combined_masks(:,:,f);
        for i=1:p
            for j=i+1:p
                if isfinite(P(i,j))
                    total_checked = total_checked + 1;
                    if P(i,j) < tau
                        if CM(i,j) || CM(j,i)
                            violations = violations + 1;
                        end
                    end
                end
            end
        end
    end
    assert(violations == 0, 'Found %d edges below tau that remained active', violations);
    details.subtests{end+1} = struct('name','below_tau_deactivated','passed',true);
    fprintf('  ✓ edges below tau are deactivated (checked %d edges)\n', total_checked);

    % If none are below tau, this test is trivially satisfied (we do not
    % require sparsity change).
    % Optional node-level sanity: active edges must connect active nodes
    for f=1:F
        CM = masks.combined_masks(:,:,f);
        NM = masks.node_masks(:,f);
        [ii,jj] = find(triu(CM,1));
        for k=1:numel(ii)
            assert(NM(ii(k)) && NM(jj(k)), 'Active edge connects an inactive node');
        end
    end
    details.subtests{end+1} = struct('name','edge_node_consistency','passed',true);
    fprintf('  ✓ combined edges connect active nodes\n');

    % Stats existence smoke
    must = {'active_edges_per_frequency','active_nodes_per_frequency','sparsity_ratios_per_frequency','overall_sparsity_ratio'};
    for k=1:numel(must)
        assert(isfield(stats, must{k}), 'Missing stats field: %s', must{k});
    end
    details.subtests{end+1} = struct('name','stats_presence','passed',true);
    fprintf('  ✓ stats presence\n');

catch ME
    passed = false; details.error_message = ME.message;
    fprintf('  ✗ %s\n', ME.message);
end
end

% -------------------------------------------------------------------------
% Test 4
% -------------------------------------------------------------------------
function [passed, details] = t4_monotonicity_active_set()
details = struct('test_name','Active-set Monotonicity','subtests',{{}}); passed = true;
try
    % Create a controlled proxy set to ensure a spread of values
    p=18; F=5; proxies = cell(F,1);
    for f=1:F
        X = abs(randn(p)); X = (X + X.')/2; X(1:p+1:end)=0;
        proxies{f} = X;
    end

    % Increasing quantiles -> increasing tau -> active set must NOT grow
    quantiles = [0.05 0.10 0.20 0.30 0.40];
    active_counts = zeros(numel(quantiles),1);

    opts = struct('verbose', false, 'force_diagonal_active', false, 'symmetrize_masks', true, 'min_active_per_node', 0);

    for q = 1:numel(quantiles)
        tau = module3_threshold_determination(proxies, quantiles(q));
        [masks, ~] = module3_combined_active_set(proxies, tau, opts);
        c = 0;
        for f=1:F
            M = masks.combined_masks(:,:,f);
            c = c + nnz(triu(M,1)); % count undirected edges once
        end
        active_counts(q) = c;
    end

    % Non-increasing sequence
    assert(all(diff(active_counts) <= 1e-12), 'Active edges are not non-increasing with higher quantiles');
    details.subtests{end+1} = struct('name','non_increasing_active_edges','passed',true);
    fprintf('  ✓ active edges non-increasing with higher quantiles\n');

catch ME
    passed = false; details.error_message = ME.message;
    fprintf('  ✗ %s\n', ME.message);
end
end

% -------------------------------------------------------------------------
% Test 5
% -------------------------------------------------------------------------
function [passed, details, viz] = t5_integration_smoke()
viz = [];  % default
details = struct('test_name','Integration Smoke','subtests',{{}}); passed = true;
try
    % Generate simulation (any allowed type)
    [~, ~, Sigma_emp, sim_params] = module7_simulation_improved_complex( ...
        'n_nodes', 15, 'n_freq', 6, 'n_samples', 120, 'graph_type', 'hub');
    p = sim_params.n_nodes; F = sim_params.n_freq;
    freqs = linspace(8,12,F);

    % Run Module 1 preprocessing to get whitened covariances
    m1_in = struct();
    m1_in.mode = 'simulation';
    m1_in.sim_results = struct('Sigma_emp', {Sigma_emp}, 'F', F, 'n', p, 'T', sim_params.n_samples);
    m1_out = module1_preprocessing_main(m1_in, 'verbose', false);

    % Use whitened covariances from Module 1 (or fallback to Sigma_emp if absent)
    if isfield(m1_out, 'Sigma_tilde') && ~isempty(m1_out.Sigma_tilde)
        S_use = m1_out.Sigma_tilde;
    else
        S_use = Sigma_emp;
    end

    % Prepare Module 3 input
    in3 = struct();
    in3.whitened_covariances = S_use;
    in3.frequencies = freqs;

    % Correlation-based active set
    p3 = struct('proxy_method','correlation','quantile_level',0.10,'verbose',false);
    res3 = module3_active_set_main(in3, p3);
    assert(isstruct(res3) && isfield(res3,'success') && res3.success);

    % Basic invariants
    assert(numel(res3.edge_proxies)==F);
    assert(isequal(size(res3.edge_active_mask), [p p F]));
    assert(isequal(size(res3.node_active_mask), [p F]));
    assert(isequal(size(res3.combined_active_mask), [p p F]));
    fprintf('  ✓ module3 pipeline with module1 output\n');

    % Optional: try Module 2 to supply precision-based proxies (if available)
    try
        L = randn(p, round(0.7*p));
        pri = repmat({eye(size(L,2))}, F, 1);
        est_in = struct('leadfield_matrix', L, ...
            'empirical_covariances', {S_use}, ...
            'source_prior_covariances', {pri}, ...
            'noise_covariance', eye(p)*0.1, ...
            'frequencies', freqs);
        est_out = module2_estep_main(est_in, struct('verbose', false));
        if isstruct(est_out) && isfield(est_out, 'initial_precision_matrices')
            in3.initial_precision_matrices = est_out.initial_precision_matrices;
            p3b = struct('proxy_method','precision','quantile_level',0.10,'verbose',false);
            res3b = module3_active_set_main(in3, p3b);
            assert(isstruct(res3b) && isfield(res3b,'success') && res3b.success);
            fprintf('  ✓ module3 pipeline with module2 precision (if available)\n');
        else
            fprintf('  ○ module2 precision not used\n');
        end
        % Build a minimal payload the visualizer expects
viz = struct();
viz.combined_masks = struct( ...
    'combined_masks', res3.combined_active_mask, ...
    'node_masks',     res3.node_active_mask);
viz.active_stats = res3.active_edge_statistics;
viz.sim_params  = sim_params;
    catch
        fprintf('  ○ module2 not available or failed (ignored)\n');
    end

catch ME
    passed = false; details.error_message = ME.message;
    fprintf('  ✗ %s\n', ME.message);
end
end

function create_active_set_visualization(combined_masks, active_stats, sim_params)
try p = sim_params.n_nodes;
    F = sim_params.n_freq;
    figure('Name', 'Module 3 Active Set Analysis', 'Position', [100, 100, 1200, 800]);
    subplot(2, 3, 1); bar(1:F, active_stats.active_edges_per_frequency);
    title('Active Edges per Frequency');
    xlabel('Frequency Index');
    ylabel('Number of Active Edges');
    grid on; subplot(2, 3, 2);
    bar(1:F, active_stats.active_nodes_per_frequency);
    title('Active Nodes per Frequency');
    xlabel('Frequency Index');
    ylabel('Number of Active Nodes');
    grid on;
    subplot(2, 3, 3);
    plot(1:F, active_stats.sparsity_ratios_per_frequency, 'o-', 'LineWidth', 2);
    title('Sparsity Ratio per Frequency');
    xlabel('Frequency Index'); ylabel('Sparsity Ratio');
    ylim([0, max(active_stats.sparsity_ratios_per_frequency) * 1.1]);
    grid on;
    subplot(2, 3, 4);
    imagesc(double(combined_masks.combined_masks(:, :, 1)));
    colormap(gray);
    colorbar;
    title('Active Set Pattern (Freq 1)');
    xlabel('Node Index');
    ylabel('Node Index');
    subplot(2, 3, 5);
    imagesc(double(combined_masks.node_masks));
    colormap(gray);
    colorbar;
    title('Node Activity Across Frequencies');
    xlabel('Frequency Index');
    ylabel('Node Index');
    subplot(2, 3, 6);
    stats_data = [active_stats.total_active_edges; ...
        active_stats.total_active_nodes; ...
        active_stats.average_edges_per_frequency; ...
        active_stats.overall_sparsity_ratio * 100];
    stats_labels = {'Total Edges','Total Nodes','Avg Edges/Freq','Sparsity %'};
    bar(stats_data);
    set(gca, 'XTickLabel', stats_labels);
    title('Active Set Summary Statistics');
    ylabel('Count / Percentage');
    xtickangle(45);
    grid on;
    sgtitle(sprintf('Module 3 Active Set Results (%dx%d, %d frequencies)', p, p, F), 'FontSize', 14, 'FontWeight', 'bold');
catch ME warning(ME.identifier, '%s', ME.message); % 按你要求使用格式化传递
end
end

% -------------------------------------------------------------------------
% Helpers
% -------------------------------------------------------------------------
function test_results = record(test_results, test_name, passed, details)
test_results.total_tests = test_results.total_tests + 1;
if passed
    test_results.passed_tests = test_results.passed_tests + 1;
else
    test_results.failed_tests = test_results.failed_tests + 1;
end
details.passed = passed;
details.test_name = test_name;
test_results.test_details{end+1} = details;
end
