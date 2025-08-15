function test_module7_simulation()
% TEST_MODULE7_SIMULATION - Unit tests for simulation data generation module
%
% This function runs comprehensive tests for the module7_simulation function
% to ensure correct generation of synthetic precision matrices.
%
% Usage:
%   test_module7_simulation()  % Run all tests
%
% Tests included:
%   1. Basic functionality and output structure
%   2. Positive definiteness of all matrices
%   3. Sparsity pattern consistency
%   4. Frequency smoothness of entries
%   5. Wishart sampling statistical properties
%   6. Different graph types
%   7. Reproducibility with random seeds
%   8. Edge cases and error handling

    fprintf('\n========================================\n');
    fprintf('Testing Module 7: Simulation Data Generation\n');
    fprintf('========================================\n\n');
    
    % Initialize test counters
    total_tests = 0;
    passed_tests = 0;
    
    % Run all tests
    [total_tests, passed_tests] = run_test(@test_basic_generation, 'Basic functionality', total_tests, passed_tests);
    [total_tests, passed_tests] = run_test(@test_positive_definiteness, 'Positive definiteness', total_tests, passed_tests);
    [total_tests, passed_tests] = run_test(@test_sparsity_pattern, 'Sparsity pattern', total_tests, passed_tests);
    [total_tests, passed_tests] = run_test(@test_frequency_smoothness, 'Frequency smoothness', total_tests, passed_tests);
    [total_tests, passed_tests] = run_test(@test_wishart_properties, 'Wishart properties', total_tests, passed_tests);
    [total_tests, passed_tests] = run_test(@test_graph_types, 'Graph types', total_tests, passed_tests);
    [total_tests, passed_tests] = run_test(@test_reproducibility, 'Reproducibility', total_tests, passed_tests);
    [total_tests, passed_tests] = run_test(@test_edge_cases, 'Edge cases', total_tests, passed_tests);
    
    % Summary
    fprintf('\n========================================\n');
    fprintf('Test Summary: %d/%d tests passed\n', passed_tests, total_tests);
    if passed_tests == total_tests
        fprintf('All tests PASSED!\n');
    else
        fprintf('Some tests FAILED!\n');
    end
    fprintf('========================================\n');
end

%% Test runner helper
function [total, passed] = run_test(test_func, test_name, total, passed)
    total = total + 1;
    fprintf('Test %d: %s... ', total, test_name);
    try
        test_func();
        fprintf('PASSED\n');
        passed = passed + 1;
    catch ME
        fprintf('FAILED\n');
        fprintf('  Error: %s\n', ME.message);
        fprintf('  Location: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
    end
end

%% Test 1: Basic functionality
function test_basic_generation()
    % Test with small parameters
    [true_prec, true_cov, emp_cov, params] = module7_simulation_improved(...
        'n_nodes', 5, ...
        'n_freq', 10, ...
        'n_samples', 50);
    
    % Check dimensions
    assert(length(true_prec) == 10, 'Wrong number of frequency points');
    assert(all(size(true_prec{1}) == [5, 5]), 'Wrong matrix dimensions');
    
    % Check all outputs are provided
    assert(~isempty(true_cov), 'True covariance is empty');
    assert(~isempty(emp_cov), 'Empirical covariance is empty');
    assert(~isempty(params), 'Parameters structure is empty');
    
    % Check params structure
    assert(isfield(params, 'n_nodes'), 'Missing n_nodes in params');
    assert(isfield(params, 'n_freq'), 'Missing n_freq in params');
    assert(isfield(params, 'edge_list'), 'Missing edge_list in params');
    assert(params.n_nodes == 5, 'Incorrect n_nodes in params');
end

%% Test 2: Positive definiteness
function test_positive_definiteness()
    [true_prec, true_cov, emp_cov, ~] = module7_simulation_improved(...
        'n_nodes', 10, ...
        'n_freq', 20, ...
        'n_samples', 100);
    
    % Check all precision matrices
    for f = 1:length(true_prec)
        eigvals_prec = eig(true_prec{f});
        assert(all(eigvals_prec > 0), ...
            sprintf('Precision matrix at freq %d is not positive definite', f));
        
        % Check symmetry
        assert(norm(true_prec{f} - true_prec{f}', 'fro') < 1e-12, ...
            sprintf('Precision matrix at freq %d is not symmetric', f));
    end
    
    % Check all covariance matrices
    for f = 1:length(true_cov)
        eigvals_cov = eig(true_cov{f});
        assert(all(eigvals_cov > 0), ...
            sprintf('Covariance matrix at freq %d is not positive definite', f));
    end
    
    % Check empirical covariances (allow small numerical errors)
    for f = 1:length(emp_cov)
        eigvals_emp = eig(emp_cov{f});
        assert(all(eigvals_emp > -1e-10), ...
            sprintf('Empirical covariance at freq %d has large negative eigenvalues', f));
    end
end

%% Test 3: Sparsity pattern
function test_sparsity_pattern()
    % Test chain graph - should have exactly n-1 edges
    n_nodes = 10;
    [true_prec, ~, ~, params] = module7_simulation(...
        'n_nodes', n_nodes, ...
        'n_freq', 5, ...
        'graph_type', 'chain');
    
    % Check correct number of edges
    expected_edges = n_nodes - 1;
    assert(params.n_edges == expected_edges, ...
        sprintf('Chain graph should have %d edges, got %d', expected_edges, params.n_edges));
    
    % Check sparsity level (considering L*L' fills in some entries)
    Omega = true_prec{1};
    threshold = max(abs(Omega(:))) * 1e-3;
    
    % Count significant off-diagonal entries
    n_nonzero = 0;
    for i = 1:n_nodes
        for j = 1:i-1
            if abs(Omega(i,j)) > threshold
                n_nonzero = n_nonzero + 1;
            end
        end
    end
    
    % For chain, expect limited fill-in
    max_expected = 3 * expected_edges; % Allow some fill-in from L*L'
    assert(n_nonzero <= max_expected, ...
        sprintf('Too many non-zeros: %d (expected <= %d)', n_nonzero, max_expected));
end

%% Test 4: Frequency smoothness
function test_frequency_smoothness()
    [true_prec, ~, ~, params] = module7_simulation(...
        'n_nodes', 5, ...
        'n_freq', 50, ...
        'n_basis', 5, ...
        'graph_type', 'hub');
    
    % Track a specific edge across frequencies
    edge_values = zeros(50, 1);
    for f = 1:50
        edge_values(f) = true_prec{f}(2, 1); % Edge (2,1)
    end
    
    % Compute roughness (sum of squared second differences)
    if length(edge_values) > 2
        second_diff = diff(diff(edge_values));
        roughness = sum(second_diff.^2);
        
        % Compare with random values
        random_values = randn(50, 1) * std(edge_values);
        random_second_diff = diff(diff(random_values));
        random_roughness = sum(random_second_diff.^2);
        
        % Smooth values should have much lower roughness
        smoothness_ratio = roughness / (random_roughness + eps);
        assert(smoothness_ratio < 0.1, ...
            sprintf('Frequency variation not smooth enough (ratio: %.3f)', smoothness_ratio));
    end
end

%% Test 5: Wishart sampling properties
function test_wishart_properties()
    % Use larger sample size for better statistics
    n_samples = 1000;
    n_nodes = 5;
    [~, true_cov, emp_cov, ~] = module7_simulation(...
        'n_nodes', n_nodes, ...
        'n_freq', 1, ...
        'n_samples', n_samples);
    
    Sigma_true = true_cov{1};
    Sigma_emp = emp_cov{1};
    
    % Test 1: Unbiasedness - empirical should approximate true covariance
    relative_error = norm(Sigma_emp - Sigma_true, 'fro') / norm(Sigma_true, 'fro');
    assert(relative_error < 0.15, ...
        sprintf('Empirical covariance error too large: %.3f', relative_error));
    
    % Test 2: Trace test - E[tr(Sigma_emp)] = tr(Sigma_true)
    trace_ratio = trace(Sigma_emp) / trace(Sigma_true);
    assert(abs(trace_ratio - 1) < 0.1, ...
        sprintf('Trace ratio is %.3f, should be close to 1', trace_ratio));
    
    % Test 3: Diagonal variance - theoretical vs empirical
    % For Wishart, Var(S_ii) = 2*Sigma_ii^2/T
    for i = 1:n_nodes
        theoretical_std = sqrt(2) * Sigma_true(i,i) / sqrt(n_samples);
        empirical_error = abs(Sigma_emp(i,i) - Sigma_true(i,i));
        % Check if error is within 3 standard deviations
        assert(empirical_error < 3 * theoretical_std, ...
            sprintf('Diagonal element %d variance inconsistent with Wishart', i));
    end
end

%% Test 6: Different graph types
function test_graph_types()
    graph_types = {'random', 'chain', 'grid', 'hub'};
    expected_properties = struct();
    
    % Define expected properties for each graph type
    n_nodes = 16; % 4x4 for grid
    expected_properties.chain = n_nodes - 1; % n-1 edges
    expected_properties.hub = n_nodes - 1;   % n-1 edges
    expected_properties.grid = 2*(4*3); % 24 edges for 4x4 grid
    expected_properties.random = -1; % Variable, check later
    
    for i = 1:length(graph_types)
        graph_type = graph_types{i};
        
        [~, ~, ~, params] = module7_simulation(...
            'n_nodes', n_nodes, ...
            'n_freq', 2, ...
            'graph_type', graph_type, ...
            'edge_density', 0.2);
        
        % Check edges were generated
        assert(params.n_edges > 0, ...
            sprintf('No edges generated for %s graph', graph_type));
        
        % Check specific properties
        if ~strcmp(graph_type, 'random')
            assert(params.n_edges == expected_properties.(graph_type), ...
                sprintf('%s graph has wrong number of edges: %d (expected %d)', ...
                graph_type, params.n_edges, expected_properties.(graph_type)));
        else
            % For random graph, check density
            max_edges = n_nodes * (n_nodes - 1) / 2;
            actual_density = params.n_edges / max_edges;
            assert(abs(actual_density - 0.2) < 0.1, ...
                'Random graph density too far from specified');
        end
    end
end

%% Test 7: Reproducibility
function test_reproducibility()
    % Generate data twice with same seed
    [prec1, cov1, emp1, ~] = module7_simulation(...
        'n_nodes', 5, ...
        'n_freq', 10, ...
        'random_seed', 123);
    
    [prec2, cov2, emp2, ~] = module7_simulation(...
        'n_nodes', 5, ...
        'n_freq', 10, ...
        'random_seed', 123);
    
    % Should be identical
    for f = 1:10
        assert(isequal(prec1{f}, prec2{f}), ...
            'Precision matrices not reproducible with same seed');
        assert(isequal(cov1{f}, cov2{f}), ...
            'Covariance matrices not reproducible with same seed');
        assert(isequal(emp1{f}, emp2{f}), ...
            'Empirical covariances not reproducible with same seed');
    end
    
    % Generate with different seed
    [prec3, ~, ~, ~] = module7_simulation(...
        'n_nodes', 5, ...
        'n_freq', 10, ...
        'random_seed', 456);
    
    % Should be different
    different = false;
    for f = 1:10
        if ~isequal(prec1{f}, prec3{f})
            different = true;
            break;
        end
    end
    assert(different, 'Results should differ with different seeds');
end

%% Test 8: Edge cases and input validation
function test_edge_cases()
    % Test 1: Minimum size
    [~, ~, ~, params] = module7_simulation(...
        'n_nodes', 2, ...
        'n_freq', 1, ...
        'n_samples', 10);
    assert(params.n_nodes == 2, 'Failed with minimum nodes');
    
    % Test 2: Large epsilon_reg ensures strong diagonal dominance
    [prec, ~, ~, ~] = module7_simulation(...
        'n_nodes', 5, ...
        'n_freq', 1, ...
        'epsilon_reg', 10);
    
    % Check strong diagonal dominance
    P = prec{1};
    for i = 1:5
        off_diag_sum = sum(abs(P(i,:))) - abs(P(i,i));
        assert(abs(P(i,i)) > off_diag_sum, ...
            'Strong diagonal dominance not achieved with large epsilon');
    end
    
    % Test 3: High edge density
    [~, ~, ~, params] = module7_simulation(...
        'n_nodes', 6, ...
        'graph_type', 'random', ...
        'edge_density', 0.9);
    
    max_edges = 6 * 5 / 2; % Complete graph has n(n-1)/2 edges
    assert(params.n_edges >= 0.8 * max_edges, ...
        'High density graph has too few edges');
end