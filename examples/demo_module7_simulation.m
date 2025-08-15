%% Demo: Module 7 - Simulation Data Generation
% This script demonstrates various use cases of the module7_simulation function
% Run sections individually or the entire script

clear; close all; clc;

% Add paths if needed
if exist('startup.m', 'file')
    run('startup.m');
end

%% Demo 1: Basic Usage with Default Parameters
fprintf('==========================================\n');
fprintf('Demo 1: Basic Usage with Default Parameters\n');
fprintf('==========================================\n\n');

% Generate data with default settings
[true_prec, true_cov, emp_cov, params] = module7_simulation();

% Display summary
fprintf('Generated data summary:\n');
fprintf('  Number of nodes: %d\n', params.n_nodes);
fprintf('  Number of frequencies: %d\n', params.n_freq);
fprintf('  Number of samples: %d\n', params.n_samples);
fprintf('  Graph type: %s\n', params.graph_type);
fprintf('  Number of edges: %d\n', params.n_edges);
fprintf('  Matrix size: %d x %d\n', size(true_prec{1}));
fprintf('\n');

% Basic visualization
visualize_simulation_data(true_prec, true_cov, emp_cov, params);

%% Demo 2: Chain Graph Structure
fprintf('\n==========================================\n');
fprintf('Demo 2: Chain Graph Structure\n');
fprintf('==========================================\n\n');

% Generate chain graph
[true_prec_chain, true_cov_chain, emp_cov_chain, params_chain] = module7_simulation(...
    'n_nodes', 10, ...
    'n_freq', 30, ...
    'n_samples', 200, ...
    'graph_type', 'chain', ...
    'n_basis', 3);

fprintf('Chain graph properties:\n');
fprintf('  Expected edges: %d\n', params_chain.n_nodes - 1);
fprintf('  Actual edges: %d\n', params_chain.n_edges);

% Show sparsity pattern
figure('Name', 'Chain Graph Sparsity Pattern');
subplot(1, 2, 1);
imagesc(abs(true_prec_chain{1}) > 1e-10);
colormap(gray);
title('Sparsity Pattern (Frequency 1)');
axis square;

subplot(1, 2, 2);
imagesc(abs(true_prec_chain{end}) > 1e-10);
colormap(gray);
title(sprintf('Sparsity Pattern (Frequency %d)', params_chain.n_freq));
axis square;

%% Demo 3: Different Graph Types Comparison
fprintf('\n==========================================\n');
fprintf('Demo 3: Comparing Different Graph Types\n');
fprintf('==========================================\n\n');

graph_types = {'random', 'chain', 'grid', 'hub'};
n_nodes = 16; % Perfect square for grid

figure('Name', 'Graph Types Comparison', 'Position', [100, 100, 1000, 800]);

for i = 1:length(graph_types)
    % Generate data for each graph type
    [prec, ~, ~, params_temp] = module7_simulation(...
        'n_nodes', n_nodes, ...
        'n_freq', 1, ...
        'graph_type', graph_types{i}, ...
        'edge_density', 0.2);
    
    % Plot sparsity pattern
    subplot(2, 2, i);
    imagesc(abs(prec{1}) > 1e-10);
    colormap(gray);
    title(sprintf('%s Graph (edges: %d)', ...
                  graph_types{i}, params_temp.n_edges));
    axis square;
    colorbar;
end

%% Demo 4: Effect of Sample Size on Estimation Error
fprintf('\n==========================================\n');
fprintf('Demo 4: Sample Size Effect on Estimation\n');
fprintf('==========================================\n\n');

sample_sizes = [50, 100, 200, 500, 1000];
errors = zeros(length(sample_sizes), 1);

for i = 1:length(sample_sizes)
    % Generate data with different sample sizes
    [~, true_cov_temp, emp_cov_temp, ~] = module7_simulation(...
        'n_nodes', 10, ...
        'n_freq', 5, ...
        'n_samples', sample_sizes(i), ...
        'random_seed', 123); % Fixed seed for comparison
    
    % Compute average relative error
    total_error = 0;
    for f = 1:5
        diff = emp_cov_temp{f} - true_cov_temp{f};
        total_error = total_error + norm(diff, 'fro') / norm(true_cov_temp{f}, 'fro');
    end
    errors(i) = total_error / 5;
    
    fprintf('  Samples: %d, Mean relative error: %.4f\n', sample_sizes(i), errors(i));
end

% Plot error vs sample size
figure('Name', 'Sample Size Effect');
loglog(sample_sizes, errors, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
% Theoretical 1/sqrt(n) rate
theoretical = errors(1) * sqrt(sample_sizes(1)) ./ sqrt(sample_sizes);
loglog(sample_sizes, theoretical, 'r--', 'LineWidth', 2);
hold off;
xlabel('Number of Samples');
ylabel('Mean Relative Error');
legend('Empirical', '1/√n theoretical', 'Location', 'best');
grid on;
title('Convergence of Empirical Covariance');

%% Demo 5: Frequency Smoothness with Different Basis Numbers
fprintf('\n==========================================\n');
fprintf('Demo 5: Effect of Basis Functions on Smoothness\n');
fprintf('==========================================\n\n');

basis_numbers = [2, 5, 10];
figure('Name', 'Basis Function Effect', 'Position', [100, 100, 1200, 400]);

for i = 1:length(basis_numbers)
    % Generate data with different numbers of basis functions
    [prec_temp, ~, ~, params_temp] = module7_simulation(...
        'n_nodes', 5, ...
        'n_freq', 50, ...
        'n_basis', basis_numbers(i), ...
        'graph_type', 'hub', ...
        'random_seed', 456);
    
    % Plot edge trajectory
    subplot(1, 3, i);
    edge_trajectory = zeros(50, 1);
    for f = 1:50
        edge_trajectory(f) = prec_temp{f}(2, 1); % Track edge (2,1)
    end
    
    plot(1:50, edge_trajectory, 'LineWidth', 2);
    xlabel('Frequency');
    ylabel('Edge (2,1) value');
    title(sprintf('%d Basis Functions', basis_numbers(i)));
    grid on;
    
    % Compute roughness (check length first)
    if length(edge_trajectory) > 2
        second_diff = diff(diff(edge_trajectory));
        roughness = sum(second_diff.^2);
        text(0.5, 0.1, sprintf('Roughness: %.2e', roughness), ...
             'Units', 'normalized', 'HorizontalAlignment', 'center');
    else
        text(0.5, 0.1, 'Too few points for roughness', ...
             'Units', 'normalized', 'HorizontalAlignment', 'center');
    end
end

%% Demo 6: Performance Testing
fprintf('\n==========================================\n');
fprintf('Demo 6: Performance Testing\n');
fprintf('==========================================\n\n');

% Test generation time for different problem sizes
node_sizes = [10, 20, 50, 100];
times = zeros(length(node_sizes), 1);

for i = 1:length(node_sizes)
    tic;
    [~, ~, ~, ~] = module7_simulation(...
        'n_nodes', node_sizes(i), ...
        'n_freq', 20, ...
        'n_samples', 100);
    times(i) = toc;
    
    fprintf('  Nodes: %d, Generation time: %.3f seconds\n', node_sizes(i), times(i));
end

% Plot scaling
figure('Name', 'Performance Scaling');
subplot(1, 2, 1);
plot(node_sizes, times, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Number of Nodes');
ylabel('Generation Time (seconds)');
title('Computation Time Scaling');
grid on;

subplot(1, 2, 2);
loglog(node_sizes, times, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
% Fit power law
p = polyfit(log(node_sizes), log(times), 1);
times_fit = exp(polyval(p, log(node_sizes)));
loglog(node_sizes, times_fit, 'r--', 'LineWidth', 2);
hold off;
xlabel('Number of Nodes');
ylabel('Generation Time (seconds)');
title(sprintf('Log-Log Scale (slope=%.2f)', p(1)));
legend('Actual', sprintf('O(n^{%.2f})', p(1)), 'Location', 'best');
grid on;

%% Demo 7: Save and Load Generated Data
fprintf('\n==========================================\n');
fprintf('Demo 7: Saving Generated Data\n');
fprintf('==========================================\n\n');

% Generate data to save
[prec_save, cov_save, emp_save, params_save] = module7_simulation(...
    'n_nodes', 15, ...
    'n_freq', 25, ...
    'graph_type', 'grid');

% Create filename with timestamp
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
filename = sprintf('data/simulated/simulation_data_%s.mat', timestamp);

% Ensure directory exists
if ~exist('data/simulated', 'dir')
    mkdir('data/simulated');
end

% Save data
save(filename, 'prec_save', 'cov_save', 'emp_save', 'params_save');
fprintf('Data saved to: %s\n', filename);
fprintf('File size: %.2f MB\n', dir(filename).bytes / 1024^2);

% Clear and reload to verify
clear prec_save cov_save emp_save params_save;
loaded_data = load(filename);
fprintf('Successfully reloaded data with fields: %s\n', ...
        strjoin(fieldnames(loaded_data), ', '));

%% Summary
fprintf('\n==========================================\n');
fprintf('Demo Complete!\n');
fprintf('==========================================\n\n');
fprintf('Key takeaways:\n');
fprintf('1. module7_simulation generates synthetic data with controllable properties\n');
fprintf('2. Different graph types produce different sparsity patterns\n');
fprintf('3. Sample size affects empirical covariance accuracy (~1/sqrt(n))\n');
fprintf('4. Number of basis functions controls frequency smoothness\n');
fprintf('5. Computation scales well with problem size\n');
fprintf('\n');

%% Optional: Run unit tests
fprintf('\nRun unit tests? (y/n): ');
user_input = input('', 's');
if strcmpi(user_input, 'y')
    fprintf('\nRunning unit tests...\n');
    test_module7_simulation();
end