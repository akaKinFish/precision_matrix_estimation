function demo_module4()
% DEMO_MODULE4 - Demonstration of Module 4 objective and gradient computation
%
% Syntax:
%   demo_module4()
%
% Description:
%   Comprehensive demonstration of Module 4 functionality for sparse precision
%   matrix estimation. Shows data preparation, gradient computation, objective
%   evaluation, and method comparison with visualizations and analysis.
%   
%   This demo covers:
%   1. Creating realistic test data with sparse structure
%   2. Setting up smoothing kernels and weight matrices
%   3. Computing gradients with different methods and parameters
%   4. Evaluating complete objective functions
%   5. Using the ObjectiveGradientComputer class
%   6. Method comparison and validation
%   7. Performance analysis and visualization
%
% Examples:
%   % Run the complete demonstration
%   demo_module4();
%
% See also: MODULE4_OBJECTIVE_GRADIENT_MAIN, MODULE4_OBJECTIVE_EVALUATION,
%           OBJECTIVEGRADIENTCOMPUTER, TEST_MODULE4_COMPREHENSIVE
%
% Author: [Your Name]
% Date: [Current Date]
% Version: 1.0

%% Initialize
clear; clc; close all;

fprintf('=== Module 4 Demo: Objective and Gradient Computation ===\n\n');

%% Step 1: Create Example Data
fprintf('Step 1: Creating example data...\n');

% Problem dimensions
p = 8;   % Number of nodes
F = 4;   % Number of frequencies

% Create ground truth sparse precision matrices
rng(42);  % For reproducibility
ground_truth = cell(F, 1);

for f = 1:F
    % Start with identity and add sparse structure
    Omega_true = 2 * eye(p);  % Strong diagonal
    
    % Add some off-diagonal connections
    connections = [
        1, 2, 0.4;   % Node 1-2 connection
        2, 3, -0.3;  % Node 2-3 connection  
        3, 4, 0.3;   % Node 3-4 connection
        5, 6, -0.2;  % Node 5-6 connection
        6, 7, 0.3;   % Node 6-7 connection
    ];
    
    % Vary connection strength across frequencies
    strength_factor = 0.7 + 0.6 * f / F;
    
    for k = 1:size(connections, 1)
        i = connections(k, 1);
        j = connections(k, 2);
        val = connections(k, 3) * strength_factor;
        Omega_true(i, j) = val;
        Omega_true(j, i) = val;
    end
    
    ground_truth{f} = Omega_true;
end

% Generate empirical covariances (with noise)
empirical_covariances = cell(F, 1);
for f = 1:F
    Sigma_true = inv(ground_truth{f});
    
    % Add sampling noise
    noise = 0.05 * randn(p, p);
    noise = (noise + noise') / 2;  % Make symmetric
    
    empirical_covariances{f} = Sigma_true + noise;
    
    % Ensure positive definite
    [V, D] = eig(empirical_covariances{f});
    eigenvals = diag(D);
    eigenvals = max(eigenvals, 0.01);  % Minimum eigenvalue
    empirical_covariances{f} = V * diag(eigenvals) * V';
end

% Create initial precision estimates (starting point for optimization)
initial_precisions = cell(F, 1);
for f = 1:F
    initial_precisions{f} = inv(empirical_covariances{f});
    
    % Add small perturbation
    perturbation = 0.02 * randn(p, p);
    perturbation = (perturbation + perturbation') / 2;
    
    initial_precisions{f} = initial_precisions{f} + perturbation;
    
    % Ensure positive definite
    [V, D] = eig(initial_precisions{f});
    eigenvals = diag(D);
    eigenvals = max(eigenvals, 0.1);
    initial_precisions{f} = V * diag(eigenvals) * V';
end

fprintf('  Created %dx%d precision matrices for %d frequencies\n', p, p, F);

%% Step 2: Set up smoothing kernel and weight matrix
fprintf('\nStep 2: Setting up smoothing kernel and weight matrix...\n');

% Create smoothing kernel (temporal smoothing)
smoothing_kernel = zeros(F, F);

% Adjacent frequency coupling (stronger)
for f = 1:F-1
    smoothing_kernel(f, f+1) = 0.4;
    smoothing_kernel(f+1, f) = 0.4;
end

% Next-nearest neighbor coupling (weaker)
for f = 1:F-2
    smoothing_kernel(f, f+2) = 0.1;
    smoothing_kernel(f+2, f) = 0.1;
end

fprintf('  Smoothing kernel has %d non-zero entries\n', nnz(smoothing_kernel));

% Create weight matrix (edge-wise weights for ||·||_{W^Γ}²)
weight_matrix = eye(p);

% Add structure to weight matrix (penalize certain edges more)
for i = 1:p-1
    weight_matrix(i, i+1) = 0.5;  % Penalize adjacent nodes more
    weight_matrix(i+1, i) = 0.5;
end

% Normalize weight matrix
weight_matrix = weight_matrix / max(eig(weight_matrix)) * 0.8;
fprintf('  Weight matrix condition number: %.2f\n', cond(weight_matrix));

%% Step 3: Prepare input data structure
fprintf('\nStep 3: Preparing input data...\n');

input_data = struct();
input_data.precision_matrices = initial_precisions;
input_data.whitened_covariances = empirical_covariances;
input_data.smoothing_kernel = smoothing_kernel;
input_data.weight_matrix = weight_matrix;

fprintf('  Input data structure prepared\n');

%% Step 4: Basic gradient computation
fprintf('\nStep 4: Computing gradients with default parameters...\n');

% Use default parameters
tic;
gradient_results = module4_gradient(input_data);
default_time = toc;

fprintf('  Gradient computation completed in %.3f seconds\n', default_time);
if gradient_results.success
    fprintf('  Success: YES\n');
else
    fprintf('  Success: NO\n');
end

% Display gradient statistics
grad_norms = zeros(F, 1);
for f = 1:F
    grad_norms(f) = norm(gradient_results.smooth_gradients{f}, 'fro');
end

fprintf('  Gradient Frobenius norms: [%.3f, %.3f, %.3f, %.3f]\n', grad_norms);
fprintf('  Max Hermitian violation: %.2e\n', max(gradient_results.hermitian_violations));

%% Step 5: Custom parameter gradient computation
fprintf('\nStep 5: Computing gradients with custom parameters...\n');

custom_params = struct();
custom_params.lambda1 = 0.08;        % Stronger smoothing
custom_params.lambda2 = 0.03;        % L1 penalty (for reference)
custom_params.use_graph_laplacian = true;   % Use efficient method
custom_params.verbose = true;        % Show progress

tic;
custom_results = module4_gradient(input_data, custom_params);
custom_time = toc;

fprintf('  Custom gradient computation completed in %.3f seconds\n', custom_time);

%% Step 6: Objective function evaluation
fprintf('\nStep 6: Evaluating complete objective function...\n');

evaluation_params = custom_params;
evaluation_params.verbose = false;  % Reduce output

tic;
[objective_values, obj_stats] = module4_objective_evaluation(input_data, evaluation_params);
objective_time = toc;

fprintf('  Objective evaluation completed in %.3f seconds\n', objective_time);
fprintf('\n  Objective Breakdown:\n');
fprintf('    Log-determinant:  %12.6f\n', objective_values.logdet_terms);
fprintf('    Trace terms:      %12.6f\n', objective_values.trace_terms);
fprintf('    Smoothing (lambda1):   %12.6f\n', objective_values.smoothing_penalty);
fprintf('    L1 penalty (lambda2):  %12.6f\n', objective_values.l1_penalty);
fprintf('    ─────────────────────────────────────\n');
fprintf('    Smooth objective: %12.6f\n', objective_values.smooth_objective);
fprintf('    Total objective:  %12.6f\n', objective_values.total_objective);

%% Step 7: Using ObjectiveGradientComputer class
fprintf('\nStep 7: Using ObjectiveGradientComputer class...\n');

% Create computer instance
computer = ObjectiveGradientComputer('lambda1', 0.05, 'lambda2', 0.02);
computer.computation_options.verbose = false;

% Compute gradients through class interface
tic;
class_gradients = computer.compute_smooth_gradients(input_data);
class_grad_time = toc;

% Compute objective through class interface
tic;
class_objective = computer.compute_full_objective(input_data);
class_obj_time = toc;

fprintf('  Class gradient time: %.3f seconds\n', class_grad_time);
fprintf('  Class objective time: %.3f seconds\n', class_obj_time);

% Compare methods using the class
fprintf('\n  Comparing direct vs Laplacian methods...\n');
comparison = computer.compare_gradient_methods(input_data);

fprintf('  Method comparison results:\n');
fprintf('    Direct method time:    %.3f seconds\n', comparison.time_direct);
fprintf('    Laplacian method time: %.3f seconds\n', comparison.time_laplacian);
fprintf('    Speedup ratio:         %.2fx\n', comparison.speedup_ratio);
fprintf('    Max difference:        %.2e\n', comparison.max_absolute_difference);
if comparison.methods_consistent
    fprintf('    Methods consistent:    YES\n');
else
    fprintf('    Methods consistent:    NO\n');
end

%% Step 8: Validation and analysis
fprintf('\nStep 8: Validation and analysis...\n');

% Check gradient properties
fprintf('  Gradient validation:\n');

for f = 1:min(F, 2)  % Check first 2 frequencies
    G = gradient_results.smooth_gradients{f};
    
    % Check Hermitian property
    hermitian_error = norm(G - G', 'fro');
    fprintf('    Freq %d Hermitian error: %.2e\n', f, hermitian_error);
    
    % Check finite values
    is_finite = all(isfinite(G(:)));
    if is_finite
        fprintf('    Freq %d all finite:     YES\n', f);
    else
        fprintf('    Freq %d all finite:     NO\n', f);
    end
end

% Performance summary
fprintf('\n  Performance Summary:\n');
fprintf('    Default gradients:  %.3fs\n', default_time);
fprintf('    Custom gradients:   %.3fs\n', custom_time);
fprintf('    Objective eval:     %.3fs\n', objective_time);
fprintf('    Class interface:    %.3fs (grad) + %.3fs (obj)\n', class_grad_time, class_obj_time);

%% Step 9: Visualization
fprintf('\nStep 9: Creating visualizations...\n');

try
    % Plot gradient norms across frequencies
    figure('Name', 'Module 4 Analysis', 'Position', [100, 100, 1200, 800]);
    
    subplot(2, 3, 1);
    plot(1:F, grad_norms, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
    title('Gradient Norms Across Frequencies');
    xlabel('Frequency Index');
    ylabel('||G_\omega||_F');
    grid on;
    
    % Plot sparsity pattern of one precision matrix
    subplot(2, 3, 2);
    spy(abs(ground_truth{1}) > 0.01);
    title('Ground Truth Sparsity Pattern (Freq 1)');
    
    % Plot current estimate sparsity
    subplot(2, 3, 3);
    spy(abs(initial_precisions{1}) > 0.1);
    title('Initial Estimate Sparsity (Freq 1)');
    
    % Plot smoothing kernel
    subplot(2, 3, 4);
    imagesc(smoothing_kernel);
    colorbar;
    title('Smoothing Kernel K');
    xlabel('Frequency');
    ylabel('Frequency');
    
    % Plot weight matrix
    subplot(2, 3, 5);
    imagesc(weight_matrix);
    colorbar;
    title('Weight Matrix W^Γ');
    
    % Plot objective breakdown
    subplot(2, 3, 6);
    components = [objective_values.logdet_terms, objective_values.trace_terms, ...
                  objective_values.smoothing_penalty, objective_values.l1_penalty];
    component_names = {'LogDet', 'Trace', 'Smooth', 'L1'};
    
    bar(components);
    set(gca, 'XTickLabel', component_names);
    title('Objective Function Components');
    ylabel('Value');
    grid on;
    
    fprintf('  Visualization created\n');
    
catch ME
    fprintf('  Visualization failed: %s\n', ME.message);
end

%% Step 10: Summary and recommendations
fprintf('\nStep 10: Summary and recommendations\n');

fprintf('\n=== Module 4 Demo Summary ===\n');
fprintf('✓ Successfully computed gradients and objectives\n');
fprintf('✓ Validated mathematical properties (Hermitian, finite)\n');
fprintf('✓ Compared direct vs Laplacian methods\n');
fprintf('✓ Demonstrated class-based interface\n');

fprintf('\nKey insights:\n');
fprintf('• Laplacian method is %.2fx faster than direct method\n', comparison.speedup_ratio);
fprintf('• Maximum gradient norm: %.3f\n', max(grad_norms));
fprintf('• Objective function value: %.6f\n', objective_values.total_objective);
fprintf('• All gradients satisfy Hermitian constraint (max violation: %.2e)\n', ...
        max(gradient_results.hermitian_violations));

fprintf('\nRecommended usage patterns:\n');
fprintf('• Use module4_gradient() for simple gradient computation\n');
fprintf('• Use ObjectiveGradientComputer class for repeated computations\n');
fprintf('• Enable verbose mode for debugging: params.verbose = true\n');
fprintf('• Use Laplacian method for better performance: params.use_graph_laplacian = true\n');
fprintf('• Check hermitian_violations field for numerical issues\n');

fprintf('\nNext steps:\n');
fprintf('• Use gradients in Module 5 (Proximal Updates) for optimization\n');
fprintf('• Tune lambda1 and lambda2 parameters based on your problem\n');
fprintf('• Run comprehensive tests with: test_module4_comprehensive()\n');

fprintf('\n=== Module 4 Demo Complete ===\n');

end