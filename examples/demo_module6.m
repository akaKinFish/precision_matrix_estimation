%% Module 6 Demonstration: Hyperparameter Configuration
% This script demonstrates Module 6 (Hyperparameter Configuration).
% It auto-computes λ₁ and α using safe bounds (Gershgorin for W^Γ, optional for L_logdet).
%
% File: examples/demo_module6.m

clear; clc;
fprintf('========================================\n');
fprintf('Module 6: Hyperparameter Configuration Demo\n');
fprintf('========================================\n\n');

%% Step 1: Basic hyperparameter computation
fprintf('Step 1: Basic hyperparameter computation\n');
fprintf('----------------------------------------\n');

n_nodes = 6;
n_freq  = 4;
test_input = create_demo_input(n_nodes, n_freq);

tic;
config = module6_hyperparameters(test_input);
computation_time = toc;

fprintf('Computed hyperparameters:\n');
fprintf('  λ₁ = %.6e (smoothing strength)\n', config.lambda1);
fprintf('  α  = %.6e (step size)\n', config.alpha);
fprintf('  λ₂ = %.6e (suggested sparsity)\n', config.lambda2_suggested);
fprintf('  Computation time: %.3f seconds\n', computation_time);

%% Step 2: Safety margin sensitivity analysis
fprintf('\nStep 2: Safety margin sensitivity analysis\n');
fprintf('-----------------------------------------\n');

safety_margins = [0.5, 0.7, 0.9, 0.95];
lambda1_values = zeros(size(safety_margins));
alpha_values   = zeros(size(safety_margins));

fprintf('Safety Margin | λ₁        | α         \n');
fprintf('-------------|-----------|----------\n');

for i = 1:numel(safety_margins)
    margin = safety_margins(i);
    config_temp = module6_hyperparameters(test_input, 'safety_margin', margin, 'verbose', false);

    lambda1_values(i) = config_temp.lambda1;
    alpha_values(i)   = config_temp.alpha;

    fprintf('    %.2f     | %.3e | %.3e\n', margin, lambda1_values(i), alpha_values(i));
end

lambda1_monotonic = all(diff(lambda1_values) <= 0);
alpha_monotonic   = all(diff(alpha_values)   <= 0);

fprintf('\nMonotonicity check:\n');
fprintf('  λ₁ decreases with safety margin: %s\n', tern(lambda1_monotonic, 'YES', 'NO'));
fprintf('  α  decreases with safety margin: %s\n', tern(alpha_monotonic,   'YES', 'NO'));

%% Step 3: Problem size scaling behavior
fprintf('\nStep 3: Problem size scaling behavior\n');
fprintf('------------------------------------\n');

problem_sizes = [4, 6, 8, 10];
scaling_results = struct();
scaling_results.n_nodes = problem_sizes;
scaling_results.lambda1 = zeros(size(problem_sizes));
scaling_results.alpha   = zeros(size(problem_sizes));
scaling_results.times   = zeros(size(problem_sizes));

fprintf('Problem Size | λ₁        | α         | Time (s)\n');
fprintf('------------|-----------|-----------|--------\n');

for i = 1:numel(problem_sizes)
    p = problem_sizes(i);
    scale_input = create_demo_input(p, 3);  % fixed frequency count

    t0 = tic;
    scale_config = module6_hyperparameters(scale_input, 'verbose', false);
    exec_time = toc(t0);

    scaling_results.lambda1(i) = scale_config.lambda1;
    scaling_results.alpha(i)   = scale_config.alpha;
    scaling_results.times(i)   = exec_time;

    fprintf('     %2d     | %.3e | %.3e | %.4f\n', p, scaling_results.lambda1(i), scaling_results.alpha(i), exec_time);
end

fprintf('\nScaling analysis:\n');
fprintf('  Time complexity: %.2fx increase from 4 to 10 nodes\n', ...
        scaling_results.times(end) / max(1e-12, scaling_results.times(1)));

%% Step 4: Method comparison (Gershgorin vs exact)
fprintf('\nStep 4: Method comparison (Gershgorin vs exact)\n');
fprintf('----------------------------------------------\n');

comparison_input = create_demo_input(8, 4);

t0 = tic;
config_exact = module6_hyperparameters(comparison_input, 'use_gershgorin', false, 'verbose', false);
time_exact = toc(t0);

t0 = tic;
config_gersh = module6_hyperparameters(comparison_input, 'use_gershgorin', true,  'verbose', false);
time_gersh = toc(t0);

L_exact = config_exact.diagnostics.L_logdet;
L_gersh = config_gersh.diagnostics.L_logdet;
bound_ratio = L_gersh / max(1e-12, L_exact);

fprintf('Method comparison:\n');
fprintf('                  | Exact      | Gershgorin | Ratio\n');
fprintf('------------------|------------|------------|------\n');
fprintf('L_logdet          | %.3e | %.3e | %.2f\n', L_exact, L_gersh, bound_ratio);
fprintf('λ₁                | %.3e | %.3e | %.2f\n', config_exact.lambda1, config_gersh.lambda1, config_gersh.lambda1/max(1e-12, config_exact.lambda1));
fprintf('α                 | %.3e | %.3e | %.2f\n', config_exact.alpha,   config_gersh.alpha,   config_gersh.alpha/max(1e-12, config_exact.alpha));
fprintf('Computation time  | %.4f s  | %.4f s  | %.2fx\n', time_exact, time_gersh, max(1e-12, time_gersh)/max(1e-12, time_exact));

fprintf('\nGershgorin efficiency: %.1fx speedup, %.1fx conservative bound\n', ...
        max(1e-12, time_exact)/max(1e-12, time_gersh), bound_ratio);

%% Step 5: Integration with Module 7 data (if available)
fprintf('\nStep 5: Integration test with Module 7\n');
fprintf('-------------------------------------\n');

try
    [~, ~, Sigma_emp, ~] = module7_simulation_improved_complex( ...
        'n_nodes', 6, 'n_freq', 4, 'n_samples', 80, 'verbose', false);
    
    integration_input = struct();
    integration_input.whitened_covariances = Sigma_emp;
    integration_input.kernel_matrix = create_rbf_kernel(4, 0.8);
    integration_input.weight_matrix = create_correlation_weights(6);

    realistic_config = module6_hyperparameters(integration_input, 'verbose', false);

    fprintf('Module 7 integration successful:\n');
    fprintf('  λ₁ = %.6e\n', realistic_config.lambda1);
    fprintf('  α  = %.6e\n', realistic_config.alpha);
    fprintf('  λ₂ = %.6e\n', realistic_config.lambda2_suggested);

    integration_success = true;
catch ME
    fprintf('Module 7 integration failed: %s\n', ME.message);
    fprintf('Using synthetic data instead...\n');

    synthetic_config = module6_hyperparameters(create_demo_input(6, 4), 'verbose', false);
    fprintf('Synthetic data results:\n');
    fprintf('  λ₁ = %.6e\n', synthetic_config.lambda1);
    fprintf('  α  = %.6e\n', synthetic_config.alpha);

    integration_success = false;
end

%% Step 6: Parameter quality assessment
fprintf('\nStep 6: Parameter quality assessment\n');
fprintf('-----------------------------------\n');

quality_input = create_demo_input(8, 5);
quality_config = module6_hyperparameters(quality_input, 'verbose', false);

if isfield(quality_config, 'quality_assessment')
    q = quality_config.quality_assessment;
    fprintf('Quality assessment:\n');
    fprintf('  λ₁ reasonable: %s\n', tern(q.lambda1_reasonable, 'YES', 'NO'));
    fprintf('  α  reasonable: %s\n', tern(q.alpha_reasonable,   'YES', 'NO'));
    fprintf('  Condition acceptable: %s\n', tern(q.condition_acceptable, 'YES', 'NO'));
    fprintf('  Overall quality: %s\n', tern(q.overall_quality, 'GOOD', 'NEEDS ATTENTION'));
else
    fprintf('Quality assessment not available in output\n');
end

%% Step 7: Diagnostic information
fprintf('\nStep 7: Diagnostic information\n');
fprintf('-----------------------------\n');

if isfield(quality_config, 'diagnostics')
    d = quality_config.diagnostics;
    fprintf('Mathematical bounds:\n');
    fprintf('  L_logdet = %.4e\n', d.L_logdet);
    fprintf('  K_max    = %.4f\n', d.K_max);
    fprintf('  R_max    = %.4f\n', d.R_max);
    if isfield(d, 'condition_numbers')
        fprintf('  Mean cond(Σ): %.2e\n', mean(d.condition_numbers));
        fprintf('  Max  cond(Σ): %.2e\n', max(d.condition_numbers));
    end
else
    fprintf('Diagnostics not available in output\n');
end

%% Step 8: Summary and recommendations
fprintf('\nStep 8: Summary and recommendations\n');
fprintf('----------------------------------\n');

fprintf('\n=== Module 6 Configuration Summary ===\n');
fprintf('✓ Automatic λ₁ and α computation using safe bounds\n');
fprintf('✓ Safety margin control for convergence\n');
fprintf('✓ Empirical λ₂ suggestion with diagnostics\n');

if integration_success
    fprintf('✓ Successfully integrated with Module 7 data\n');
else
    fprintf('⚠ Module 7 integration needs attention\n');
end

fprintf('\nRecommended usage:\n');
fprintf('• Use default safety_margin=0.9 for robust convergence\n');
fprintf('• Start with suggested λ₂, then tune based on sparsity needs\n');
fprintf('• Monitor diagnostics and quality assessment\n');
fprintf('• Use Gershgorin (default) for computational efficiency\n');

final_config = module6_hyperparameters(create_demo_input(6, 4), 'verbose', false);
fprintf('\nTypical parameter values (6 nodes, 4 frequencies):\n');
fprintf('• λ₁ = %.3e\n', final_config.lambda1);
fprintf('• α  = %.3e\n', final_config.alpha);
fprintf('• λ₂ = %.3e\n', final_config.lambda2_suggested);

fprintf('\n=== Module 6 Demo Complete ===\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local helpers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function input_struct = create_demo_input(n_nodes, n_freq)
% Create demo input structure with realistic SPD covariances.

input_struct = struct();
input_struct.whitened_covariances = cell(n_freq, 1);
for f = 1:n_freq
    A = randn(n_nodes); A = A + A';           % symmetric
    min_eig = min(eig(A));
    if min_eig <= 0
        A = A + (abs(min_eig) + 0.1) * eye(n_nodes);
    end
    input_struct.whitened_covariances{f} = (A + A')/2;
end

input_struct.kernel_matrix = create_rbf_kernel(n_freq, 1.0);
input_struct.weight_matrix = create_correlation_weights(n_nodes);
end

function K = create_rbf_kernel(n_freq, bandwidth)
% Create RBF kernel matrix and row-normalize.

freq_grid = linspace(0, 1, n_freq);
K = zeros(n_freq);
for i = 1:n_freq
    for j = 1:n_freq
        K(i, j) = exp(-0.5 * ((freq_grid(i) - freq_grid(j)) / bandwidth)^2);
    end
end
K = K ./ sum(K, 2);
end

function W = create_correlation_weights(n_nodes)
% Create a nonnegative weight matrix with zero diagonal (as in theory).

B = randn(n_nodes); B = B + B';    % symmetric
W = abs(B);
W = W ./ max(W(:));                % normalize to [0,1]
W(1:n_nodes+1:end) = 0;            % zero diagonal for theory consistency
end

function s = tern(cond, a, b)
if cond, s = a; else, s = b; end
end
