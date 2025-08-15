%% Demo: Improved Simulation with Dynamic Sparsity
% This script demonstrates the enhanced simulation capabilities
% Location: examples/demo_module7_improved.m

clear; close all; clc;

%% Example 1: Compare Original vs Improved Version
fprintf('=========================================\n');
fprintf('Example 1: Original vs Improved Simulation\n');
fprintf('=========================================\n\n');

% Common parameters - defined as a proper cell array
common_params = {...
    'n_nodes', 15, ...
    'n_freq', 40, ...
    'n_samples', 200, ...
    'graph_type', 'random', ...
    'edge_density', 0.15, ...
    'random_seed', 123 ...
};

% Generate with original method
fprintf('Generating with original method...\n');
[prec_orig, cov_orig, emp_orig, params_orig] = module7_simulation(common_params{:});

% Generate with improved method
fprintf('Generating with improved method...\n');
[prec_imp, cov_imp, emp_imp, params_imp] = module7_simulation_improved(...
    common_params{:}, ...
    'sparsity_variation', 0.3, ...
    'edge_activation_smoothness', 0.7);

% Compare sparsity patterns
figure('Name', 'Original vs Improved Comparison', 'Position', [100, 100, 1200, 600]);

% Original method visualization
subplot(2, 3, 1);
imagesc(abs(prec_orig{1}) > 1e-10);
colormap(gray);
title('Original: Freq 1');
axis square;

subplot(2, 3, 2);
imagesc(abs(prec_orig{20}) > 1e-10);
colormap(gray);
title('Original: Freq 20');
axis square;

subplot(2, 3, 3);
imagesc(abs(prec_orig{40}) > 1e-10);
colormap(gray);
title('Original: Freq 40');
axis square;

% Improved method visualization
subplot(2, 3, 4);
imagesc(abs(prec_imp{1}) > 1e-10);
colormap(gray);
title('Improved: Freq 1');
axis square;

subplot(2, 3, 5);
imagesc(abs(prec_imp{20}) > 1e-10);
colormap(gray);
title('Improved: Freq 20');
axis square;

subplot(2, 3, 6);
imagesc(abs(prec_imp{40}) > 1e-10);
colormap(gray);
title('Improved: Freq 40');
axis square;

% Analyze differences
fprintf('\nSparsity Analysis:\n');
if isfield(params_orig, 'sparsity_changes')
    fprintf('Original - pattern changes: %d\n', params_orig.sparsity_changes);
else
    % Calculate sparsity changes for original method
    n_changes = 0;
    threshold = 1e-10;
    for f = 2:length(prec_orig)
        prev_pattern = abs(prec_orig{f-1}) > threshold;
        curr_pattern = abs(prec_orig{f}) > threshold;
        % Count edge changes (upper triangular only)
        for i = 1:params_orig.n_nodes
            for j = i+1:params_orig.n_nodes
                if prev_pattern(i,j) ~= curr_pattern(i,j)
                    n_changes = n_changes + 1;
                end
            end
        end
    end
    fprintf('Original - pattern changes: %d (calculated)\n', n_changes);
end

if isfield(params_imp, 'sparsity_changes')
    fprintf('Improved - pattern changes: %d\n', params_imp.sparsity_changes);
    fprintf('Improved - base edges: %d\n', params_imp.n_base_edges);
    fprintf('Improved - variable edges: %d\n', params_imp.n_variable_edges);
else
    fprintf('Improved - pattern changes: Not available\n');
end

%% Example 2: Effect of Sparsity Variation Parameter
fprintf('\n=========================================\n');
fprintf('Example 2: Effect of Sparsity Variation\n');
fprintf('=========================================\n\n');

sparsity_vars = [0, 0.2, 0.5, 0.8];
pattern_changes = zeros(length(sparsity_vars), 1);

figure('Name', 'Sparsity Variation Effect', 'Position', [100, 100, 1000, 800]);

for i = 1:length(sparsity_vars)
    % Generate with different sparsity variation
    [prec_temp, ~, ~, params_temp] = module7_simulation_improved(...
        'n_nodes', 10, ...
        'n_freq', 30, ...
        'graph_type', 'chain', ...
        'sparsity_variation', sparsity_vars(i), ...
        'random_seed', 456);
    
    pattern_changes(i) = params_temp.sparsity_changes;
    
    % Plot middle frequency
    subplot(2, 2, i);
    imagesc(abs(prec_temp{15}));
    colormap(hot);
    colorbar;
    title(sprintf('Sparsity Var = %.1f (changes: %d)', ...
                  sparsity_vars(i), pattern_changes(i)));
    axis square;
end

%% Example 3: Visualize Dynamic Sparsity Evolution
fprintf('\n=========================================\n');
fprintf('Example 3: Dynamic Sparsity Evolution\n');
fprintf('=========================================\n\n');

% Generate data with clear dynamic behavior
[prec_dyn, cov_dyn, emp_dyn, params_dyn] = module7_simulation_improved(...
    'n_nodes', 12, ...
    'n_freq', 50, ...
    'graph_type', 'hub', ...
    'sparsity_variation', 0.4, ...
    'edge_activation_smoothness', 0.6, ...
    'n_basis', 4);

% Use the new visualization function
visualize_dynamic_sparsity(prec_dyn, params_dyn);

%% Example 4: Smoothness vs Dynamics Trade-off
fprintf('\n=========================================\n');
fprintf('Example 4: Smoothness vs Dynamics\n');
fprintf('=========================================\n\n');

% Test different smoothness levels
smoothness_levels = [0.2, 0.5, 0.8, 0.95];
figure('Name', 'Smoothness Levels', 'Position', [100, 100, 1200, 400]);

for i = 1:length(smoothness_levels)
    [prec_smooth, ~, ~, params_smooth] = module7_simulation_improved(...
        'n_nodes', 8, ...
        'n_freq', 40, ...
        'graph_type', 'grid', ...
        'sparsity_variation', 0.3, ...
        'edge_activation_smoothness', smoothness_levels(i), ...
        'random_seed', 789);
    
    % Track a variable edge
    if ~isempty(params_smooth.variable_edge_list)
        edge = params_smooth.variable_edge_list(1, :);
        values = zeros(40, 1);
        for f = 1:40
            values(f) = prec_smooth{f}(edge(1), edge(2));
        end
        
        subplot(1, 4, i);
        plot(1:40, values, 'LineWidth', 2);
        title(sprintf('Smoothness = %.2f', smoothness_levels(i)));
        xlabel('Frequency');
        ylabel('Edge Value');
        grid on;
        ylim([-0.5, 0.5]);
    end
end

%% Example 5: Practical Scenario - Realistic EEG-like Structure
fprintf('\n=========================================\n');
fprintf('Example 5: Realistic Scenario\n');
fprintf('=========================================\n\n');

% Parameters mimicking EEG connectivity
[prec_eeg, cov_eeg, emp_eeg, params_eeg] = module7_simulation_improved(...
    'n_nodes', 32, ... % Like 32 EEG channels
    'n_freq', 100, ... % 100 frequency bins
    'n_samples', 500, ... % Reasonable sample size
    'graph_type', 'random', ...
    'edge_density', 0.1, ... % Sparse connectivity
    'sparsity_variation', 0.25, ... % Some edges appear/disappear
    'edge_activation_smoothness', 0.85, ... % Smooth transitions
    'n_basis', 8, ... % More basis for complex patterns
    'sigma_coef', 0.3); % Moderate edge strengths

% Analyze the results
stats = Module7Utilities.analyze_simulation_quality(prec_eeg, cov_eeg, emp_eeg, params_eeg);

% Create summary visualization
figure('Name', 'Realistic EEG-like Simulation', 'Position', [100, 100, 1000, 600]);

% Show connectivity at different frequency bands
freq_bands = [1, 25, 50, 75, 100]; % Different frequency bands
band_names = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma'};

for i = 1:5
    subplot(2, 3, i);
    f = freq_bands(i);
    
    % Create adjacency matrix
    P = prec_eeg{f};
    threshold = 0.01 * max(abs(P(:)));
    A = abs(P) > threshold;
    A = A - diag(diag(A)); % Remove diagonal
    
    % Plot as graph if possible
    if exist('graph', 'class') == 8
        G = graph(A);
        plot(G, 'Layout', 'circle', 'NodeColor', [0.3, 0.5, 0.8], ...
             'MarkerSize', 6, 'LineWidth', 1.5);
        title(sprintf('%s Band (f=%d)', band_names{i}, f));
    else
        imagesc(A);
        colormap(gray);
        title(sprintf('%s Band (f=%d)', band_names{i}, f));
    end
end

% Summary statistics
subplot(2, 3, 6);
text(0.1, 0.9, 'Simulation Summary:', 'FontWeight', 'bold', 'FontSize', 12);
text(0.1, 0.7, sprintf('Nodes: %d', params_eeg.n_nodes));
text(0.1, 0.6, sprintf('Frequencies: %d', params_eeg.n_freq));
text(0.1, 0.5, sprintf('Base edges: %d', params_eeg.n_base_edges));
text(0.1, 0.4, sprintf('Variable edges: %d', params_eeg.n_variable_edges));
text(0.1, 0.3, sprintf('Pattern changes: %d', params_eeg.sparsity_changes));
text(0.1, 0.2, sprintf('Estimation error: %.3f', stats.cov_errors.relative_mean));
text(0.1, 0.1, sprintf('Min condition #: %.1f', stats.conditioning.precision_mean));
axis off;

%% Summary and Recommendations
fprintf('\n=========================================\n');
fprintf('Summary and Recommendations\n');
fprintf('=========================================\n\n');

fprintf('Key insights from improved simulation:\n');
fprintf('1. Dynamic sparsity creates more realistic patterns\n');
fprintf('2. Sparsity variation 0.2-0.3 provides good balance\n');
fprintf('3. Edge activation smoothness 0.7-0.9 prevents abrupt changes\n');
fprintf('4. More basis functions allow complex frequency patterns\n');
fprintf('\n');

fprintf('Recommended parameter ranges:\n');
fprintf('- sparsity_variation: 0.2-0.4 for moderate dynamics\n');
fprintf('- edge_activation_smoothness: 0.7-0.9 for smooth transitions\n');
fprintf('- n_basis: 5-10 depending on frequency resolution\n');
fprintf('- sigma_coef: 0.3-0.5 for reasonable SNR\n');

%% Save an example dataset
if ~exist('data/simulated', 'dir')
    mkdir('data/simulated');
end

timestamp = datestr(now, 'yyyymmdd_HHMMSS');
filename = sprintf('data/simulated/improved_simulation_%s.mat', timestamp);
save(filename, 'prec_eeg', 'cov_eeg', 'emp_eeg', 'params_eeg');
fprintf('\nExample dataset saved to: %s\n', filename);