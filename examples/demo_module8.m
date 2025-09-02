function demo_results = demo_module8()
% DEMO_MODULE8 - Comprehensive demonstration of Module 8 recoloring functionality
%
% Returns:
%   demo_results - Structure containing demo results, metrics, and visualizations
%
% Description:
%   This demo shows the complete Module 8 recoloring functionality using
%   the simplified pipeline (Module 7 → Module 1 → Module 8) without
%   complex optimization. It demonstrates:
%
%   1. Integration with Module 7 simulation data
%   2. Integration with Module 1 whitening matrices
%   3. Theoretical precision matrix recovery
%   4. Quality assessment and validation
%   5. Visualization of results and error analysis
%   6. Performance metrics and recommendations
%
% Features:
%   - Uses realistic EEG-like parameters from Module 7
%   - Tests both real and complex Hermitian matrices
%   - Provides comprehensive error analysis
%   - Creates detailed visualizations
%   - Benchmarks against ground truth
%
% Usage:
%   demo_results = demo_module8();
%
% File location: examples/demo_module8.m

fprintf('========================================\n');
fprintf('Module 8 Recoloring Demonstration\n');
fprintf('Simplified Pipeline: Module 7→1→8\n');
fprintf('========================================\n\n');

% Initialize demo results
demo_results = struct();
demo_results.timestamp = datestr(now);
demo_results.matlab_version = version();

%% Demo 1: Basic Recoloring with Real Data
fprintf('=== Demo 1: Basic Recoloring with Real Matrices ===\n');
fprintf('----------------------------------------------------\n');

try
    demo1_params = struct();
    demo1_params.n_nodes = 8;
    demo1_params.n_freq = 6;
    demo1_params.n_samples = 120;
    demo1_params.graph_type = 'random';
    demo1_params.edge_density = 0.3;
    
    fprintf('Generating Module 7 simulation data:\n');
    fprintf('  - Nodes: %d\n', demo1_params.n_nodes);
    fprintf('  - Frequencies: %d\n', demo1_params.n_freq);
    fprintf('  - Samples: %d\n', demo1_params.n_samples);
    fprintf('  - Graph: %s (density: %.2f)\n', demo1_params.graph_type, demo1_params.edge_density);
    
    [demo1_results, demo1_data] = run_basic_recoloring_demo(demo1_params);
    demo_results.demo1 = demo1_results;
    
    if demo1_results.success
        fprintf('\n✓ Demo 1 successful!\n');
        fprintf('  - Recovery error: %.4f\n', demo1_results.recovery_metrics.median_recovery_error);
        fprintf('  - Cycle error: %.2e\n', demo1_results.recovery_metrics.median_cycle_error);
        fprintf('  - Processing time: %.3f seconds\n', demo1_results.processing_time);
    else
        fprintf('\n✗ Demo 1 failed: %s\n', demo1_results.error_message);
    end
    
catch ME
    fprintf('\n✗ Demo 1 failed: %s\n', ME.message);
    demo_results.demo1.success = false;
    demo_results.demo1.error_message = ME.message;
end

%% Demo 2: Complex Hermitian Data
fprintf('\n=== Demo 2: Complex Hermitian Matrices ===\n');
fprintf('-------------------------------------------\n');

try
    demo2_params = struct();
    demo2_params.n_nodes = 6;
    demo2_params.n_freq = 4;
    demo2_params.n_samples = 100;
    demo2_params.graph_type = 'hub';
    demo2_params.edge_density = 0.25;
    demo2_params.complex_strength = 1.2;  % Enable complex data
    
    fprintf('Testing complex Hermitian matrix support:\n');
    fprintf('  - Complex strength: %.1f\n', demo2_params.complex_strength);
    
    [demo2_results, demo2_data] = run_complex_recoloring_demo(demo2_params);
    demo_results.demo2 = demo2_results;
    
    if demo2_results.success
        fprintf('\n✓ Demo 2 successful!\n');
        fprintf('  - Complex matrices processed: %d/%d\n', ...
                demo2_results.n_complex_matrices, demo2_params.n_freq);
        fprintf('  - Recovery error: %.4f\n', demo2_results.recovery_metrics.median_recovery_error);
        fprintf('  - All Hermitian: %s\n', logical_to_string(demo2_results.all_hermitian_preserved));
    else
        fprintf('\n✗ Demo 2 failed: %s\n', demo2_results.error_message);
    end
    
catch ME
    fprintf('\n✗ Demo 2 failed: %s\n', ME.message);
    demo_results.demo2.success = false;
    demo_results.demo2.error_message = ME.message;
end

%% Demo 3: Performance Scaling Analysis
fprintf('\n=== Demo 3: Performance Scaling Analysis ===\n');
fprintf('---------------------------------------------\n');

try
    fprintf('Testing performance across different problem sizes...\n');
    
    [demo3_results, demo3_data] = run_performance_scaling_demo();
    demo_results.demo3 = demo3_results;
    
    if demo3_results.success
        fprintf('\n✓ Demo 3 successful!\n');
        fprintf('  - Tested problem sizes: %s\n', demo3_results.size_summary);
        fprintf('  - Performance scaling: %s\n', demo3_results.scaling_assessment);
        fprintf('  - Average time per op: %.2e s/op\n', demo3_results.avg_time_per_operation);
    else
        fprintf('\n✗ Demo 3 failed: %s\n', demo3_results.error_message);
    end
    
catch ME
    fprintf('\n✗ Demo 3 failed: %s\n', ME.message);
    demo_results.demo3.success = false;
    demo_results.demo3.error_message = ME.message;
end

%% Demo 4: Error Analysis and Quality Assessment
fprintf('\n=== Demo 4: Error Analysis and Quality Assessment ===\n');
fprintf('-----------------------------------------------------\n');

try
    fprintf('Performing detailed error analysis...\n');
    
    [demo4_results, demo4_data] = run_error_analysis_demo();
    demo_results.demo4 = demo4_results;
    
    if demo4_results.success
        fprintf('\n✓ Demo 4 successful!\n');
        fprintf('  - Error components analyzed: %d\n', length(demo4_results.error_components));
        fprintf('  - Quality score: %.2f/1.0\n', demo4_results.overall_quality_score);
        fprintf('  - Robustness assessment: %s\n', demo4_results.robustness_level);
    else
        fprintf('\n✗ Demo 4 failed: %s\n', demo4_results.error_message);
    end
    
catch ME
    fprintf('\n✗ Demo 4 failed: %s\n', ME.message);
    demo_results.demo4.success = false;
    demo_results.demo4.error_message = ME.message;
end

%% Summary and Visualization
% Generate overall assessment FIRST (it is used by visualization)
demo_results.overall_assessment = generate_overall_assessment(demo_results);

fprintf('\n=== Creating Comprehensive Visualization ===\n');

try
    if can_make_plots()
        create_module8_comprehensive_visualization(demo_results);
        demo_results.visualization_created = true;
        fprintf('✓ Comprehensive visualization created\n');
    else
        fprintf('⚠ Plotting disabled - skipping visualization\n');
        demo_results.visualization_created = false;
    end
    
catch ME
    fprintf('✗ Visualization failed: %s\n', ME.message);
    demo_results.visualization_created = false;
    demo_results.visualization_error = ME.message;
end

fprintf('\n========================================\n');
fprintf('Module 8 Demonstration Summary\n');
fprintf('========================================\n');

successful_demos = sum([
    safe_get_field(demo_results, 'demo1', 'success', false), ...
    safe_get_field(demo_results, 'demo2', 'success', false), ...
    safe_get_field(demo_results, 'demo3', 'success', false), ...
    safe_get_field(demo_results, 'demo4', 'success', false) ...
]);

fprintf('Successful demos: %d/4\n', successful_demos);
fprintf('Overall assessment: %s\n', demo_results.overall_assessment.summary);

if demo_results.overall_assessment.ready_for_production
    fprintf('\n🎉 MODULE 8 IS READY FOR PRODUCTION USE!\n');
else
    fprintf('\n⚠️  Module 8 needs attention before production use.\n');
end

fprintf('\nRecommendations:\n');
for i = 1:length(demo_results.overall_assessment.recommendations)
    fprintf('  • %s\n', demo_results.overall_assessment.recommendations{i});
end

fprintf('\n========================================\n');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo Implementation Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [demo1_results, demo1_data] = run_basic_recoloring_demo(demo_params)
% Run basic recoloring demonstration (uses Σ_emp for "realistic" error)

demo1_results = struct();
demo1_results.success = false;

tic;

% Generate Module 7 data
[Omega_true, Sigma_true, Sigma_emp, sim_params] = ...
    module7_simulation_improved_complex( ...
        'n_nodes', demo_params.n_nodes, ...
        'n_freq', demo_params.n_freq, ...
        'n_samples', demo_params.n_samples, ...
        'graph_type', demo_params.graph_type, ...
        'edge_density', demo_params.edge_density, ...
        'random_seed', 123);

fprintf('  ✓ Module 7 data generated: %d SPD matrices\n', length(Omega_true));

% Apply Module 1 preprocessing (uses Σ_emp to compute D and Σ̃)
module1_input = struct();
module1_input.mode = 'simulation';
module1_input.sim_results = struct();
module1_input.sim_results.Sigma_emp = Sigma_emp;
module1_input.sim_results.F = length(Sigma_emp);
module1_input.sim_results.n = size(Sigma_emp{1}, 1);
module1_input.sim_results.T = demo_params.n_samples;

module1_results = module1_preprocessing_main(module1_input, struct('verbose', false));
fprintf('  ✓ Module 1 preprocessing completed\n');

% Create theoretical whitened precision matrices
F = length(Sigma_emp);
Gamma_tilde_theoretical = cell(F, 1);

for omega = 1:F
    Sigma_whitened = module1_results.Sigma_tilde{omega};
    Sigma_whitened = (Sigma_whitened + Sigma_whitened') / 2;  % Ensure Hermitian
    
    % Add regularization if needed
    min_eigenval = min(real(eig(Sigma_whitened)));
    if min_eigenval <= 1e-10
        reg_factor = abs(min_eigenval) + 1e-8;
        Sigma_whitened = Sigma_whitened + reg_factor * eye(size(Sigma_whitened));
    end
    
    Gamma_tilde_theoretical{omega} = inv(Sigma_whitened);
end

fprintf('  ✓ Theoretical whitened precision matrices computed\n');

% Apply Module 8 recoloring (uses Σ_emp here to show realistic recovery error)
module8_input = struct();
module8_input.whitened_precision_matrices = Gamma_tilde_theoretical;
module8_input.whitening_matrices = module1_results.D;
module8_input.original_covariances = Sigma_emp;  % keep empirical here on purpose

module8_params = struct('verbose', true, 'compute_quality_metrics', true);
module8_results = module8_recoloring_main(module8_input, module8_params);

fprintf('  ✓ Module 8 recoloring completed\n');

% Compute recovery metrics
Omega_estimated = module8_results.recolored_precision_matrices;
demo1_results.recovery_metrics = compute_recovery_metrics_demo(Omega_estimated, Omega_true, module8_results);

demo1_results.processing_time = toc;
demo1_results.success = module8_results.success && demo1_results.recovery_metrics.success;

% Store data for visualization
demo1_data = struct();
demo1_data.ground_truth = struct('Omega_true', {Omega_true}, 'Sigma_true', {Sigma_true});
demo1_data.module1_results = module1_results;
demo1_data.module8_results = module8_results;
demo1_data.demo_params = demo_params;

end

function [demo2_results, demo2_data] = run_complex_recoloring_demo(demo_params)
% Run complex Hermitian recoloring demonstration

demo2_results = struct();
demo2_results.success = false;

% Generate complex Module 7 data
[Omega_true, Sigma_true, Sigma_emp, sim_params] = ...
    module7_simulation_improved_complex( ...
        'n_nodes', demo_params.n_nodes, ...
        'n_freq', demo_params.n_freq, ...
        'n_samples', demo_params.n_samples, ...
        'graph_type', demo_params.graph_type, ...
        'edge_density', demo_params.edge_density, ...
        'complex_strength', demo_params.complex_strength, ...
        'random_seed', 456);

% Count complex matrices
n_complex = sum(cellfun(@(x) ~isreal(x), Omega_true));
fprintf('  ✓ Generated %d complex matrices out of %d\n', n_complex, demo_params.n_freq);

% Apply Module 1 preprocessing
module1_input = struct();
module1_input.mode = 'simulation';
module1_input.sim_results = struct();
module1_input.sim_results.Sigma_emp = Sigma_emp;
module1_input.sim_results.F = length(Sigma_emp);
module1_input.sim_results.n = size(Sigma_emp{1}, 1);
module1_input.sim_results.T = demo_params.n_samples;

module1_results = module1_preprocessing_main(module1_input, struct('verbose', false));

% Create theoretical whitened precision matrices (complex-aware)
F = length(Sigma_emp);
Gamma_tilde_theoretical = cell(F, 1);

for omega = 1:F
    Sigma_whitened = module1_results.Sigma_tilde{omega};
    Sigma_whitened = (Sigma_whitened + Sigma_whitened') / 2;  % Force Hermitian
    
    % Regularization for complex matrices
    eigenvals = eig(Sigma_whitened);
    min_real_eigenval = min(real(eigenvals));
    if min_real_eigenval <= 1e-10
        reg_factor = abs(min_real_eigenval) + 1e-8;
        Sigma_whitened = Sigma_whitened + reg_factor * eye(size(Sigma_whitened));
    end
    
    Gamma_tilde_theoretical{omega} = inv(Sigma_whitened);
end

% Apply Module 8 recoloring (empirical covariances to illustrate realistic error)
module8_input = struct();
module8_input.whitened_precision_matrices = Gamma_tilde_theoretical;
module8_input.whitening_matrices = module1_results.D;
module8_input.original_covariances = Sigma_emp;

module8_results = module8_recoloring_main(module8_input, struct('verbose', false));

% Analyze complex data handling
demo2_results.n_complex_matrices = n_complex;
demo2_results.recovery_metrics = compute_recovery_metrics_demo( ...
    module8_results.recolored_precision_matrices, Omega_true, module8_results);

% Check Hermitian preservation for complex matrices
hermitian_preserved = true;
for omega = 1:F
    if ~isreal(Omega_true{omega}) && ~isempty(module8_results.recolored_precision_matrices{omega})
        Omega_result = module8_results.recolored_precision_matrices{omega};
        hermitian_error = norm(Omega_result - Omega_result', 'fro') / norm(Omega_result, 'fro');
        if hermitian_error > 1e-12
            hermitian_preserved = false;
            break;
        end
    end
end

demo2_results.all_hermitian_preserved = hermitian_preserved;
demo2_results.success = module8_results.success && demo2_results.recovery_metrics.success && hermitian_preserved;

% Store data
demo2_data = struct();
demo2_data.ground_truth = struct('Omega_true', {Omega_true}, 'Sigma_true', {Sigma_true});
demo2_data.module8_results = module8_results;
demo2_data.demo_params = demo_params;

end

function [demo3_results, demo3_data] = run_performance_scaling_demo()
% Run performance scaling demonstration

demo3_results = struct();
demo3_results.success = false;

problem_configs = [
    struct('n_nodes', 4, 'n_freq', 3, 'name', '4x3');
    struct('n_nodes', 8, 'n_freq', 5, 'name', '8x5');
    struct('n_nodes', 12, 'n_freq', 8, 'name', '12x8');
    struct('n_nodes', 16, 'n_freq', 10, 'name', '16x10');
];

performance_data = struct();
performance_data.problem_sizes = {};
performance_data.processing_times = [];
performance_data.matrix_operations = []; % use p^2*F as proxy for work

for i = 1:length(problem_configs)
    config = problem_configs(i);
    
    fprintf('  Testing %s... ', config.name);
    
    test_params = struct();
    test_params.n_nodes = config.n_nodes;
    test_params.n_freq = config.n_freq;
    test_params.n_samples = max(50, config.n_nodes * 2);
    test_params.graph_type = 'random';
    test_params.edge_density = 0.3;
    
    t0 = tic;
    [~, ~] = run_basic_recoloring_demo(test_params);
    processing_time = toc(t0);
    
    performance_data.problem_sizes{end+1} = config.name; %#ok<AGROW>
    performance_data.processing_times(end+1) = processing_time; %#ok<AGROW>
    performance_data.matrix_operations(end+1) = (config.n_nodes^2) * config.n_freq; %#ok<AGROW>
    
    fprintf('%.3fs\n', processing_time);
end

% Analyze scaling behavior
ops = performance_data.matrix_operations(:);
times = performance_data.processing_times(:);

% Fit a simple linear model time = a*ops + b
p = polyfit(ops, times, 1);
slope = p(1);
scaling_quality = 'good';
if slope <= 0
    scaling_quality = 'poor';
end

% Average time per operation (seconds per "p^2F" unit)
avg_time_per_operation = mean(times ./ ops);

demo3_results.size_summary = strjoin(performance_data.problem_sizes, ', ');
demo3_results.avg_time_per_operation = avg_time_per_operation;
demo3_results.scaling_assessment = scaling_quality;
demo3_results.performance_data = performance_data;
demo3_results.success = true;

demo3_data = performance_data;

end

function [demo4_results, demo4_data] = run_error_analysis_demo()
% Run detailed error analysis demonstration

demo4_results = struct();
demo4_results.success = false;

% Create test case with known characteristics
test_params = struct();
test_params.n_nodes = 10;
test_params.n_freq = 6;
test_params.n_samples = 150;
test_params.graph_type = 'random';
test_params.edge_density = 0.3;

% Run pipeline
[basic_results, basic_data] = run_basic_recoloring_demo(test_params);

if ~basic_results.success
    demo4_results.error_message = 'Basic pipeline failed';
    demo4_data = [];
    return;
end

% Analyze error components
error_analysis = analyze_error_components(basic_data);

demo4_results.error_components = error_analysis.component_names;
demo4_results.error_magnitudes = error_analysis.component_magnitudes;
demo4_results.dominant_error_source = error_analysis.dominant_source;

% Compute overall quality score
quality_components = [
    error_analysis.transformation_quality, ...
    error_analysis.numerical_quality, ...
    error_analysis.consistency_quality ...
];

demo4_results.overall_quality_score = mean(quality_components);

% Assess robustness level
if demo4_results.overall_quality_score > 0.9
    demo4_results.robustness_level = 'excellent';
elseif demo4_results.overall_quality_score > 0.8
    demo4_results.robustness_level = 'good';
elseif demo4_results.overall_quality_score > 0.6
    demo4_results.robustness_level = 'acceptable';
else
    demo4_results.robustness_level = 'needs_improvement';
end

demo4_results.success = true;
demo4_data = error_analysis;

end

function recovery_metrics = compute_recovery_metrics_demo(Omega_estimated, Omega_true, module8_results)
% Compute recovery metrics for demo purposes

recovery_metrics = struct();
recovery_metrics.success = false;

try
    F = length(Omega_estimated);
    
    % Relative recovery errors
    relative_errors = [];
    cycle_errors = [];
    
    for omega = 1:F
        if ~isempty(Omega_estimated{omega}) && ~isempty(Omega_true{omega})
            % Recovery error
            rel_err = norm(Omega_estimated{omega} - Omega_true{omega}, 'fro') / ...
                     norm(Omega_true{omega}, 'fro');
            relative_errors(end+1) = rel_err;
            
            % Cycle error from Module 8 quality
            if isfield(module8_results, 'recoloring_quality') && ...
               omega <= length(module8_results.recoloring_quality) && ...
               isfield(module8_results.recoloring_quality{omega}, 'inv_error')
                cycle_errors(end+1) = module8_results.recoloring_quality{omega}.inv_error;
            end
        end
    end
    
    % Summary statistics
    recovery_metrics.median_recovery_error = median(relative_errors);
    recovery_metrics.mean_recovery_error = mean(relative_errors);
    recovery_metrics.max_recovery_error = max(relative_errors);
    recovery_metrics.std_recovery_error = std(relative_errors);
    
    if ~isempty(cycle_errors)
        recovery_metrics.median_cycle_error = median(cycle_errors);
        recovery_metrics.max_cycle_error = max(cycle_errors);
    else
        recovery_metrics.median_cycle_error = NaN;
        recovery_metrics.max_cycle_error = NaN;
    end
    
    recovery_metrics.n_successful_frequencies = length(relative_errors);
    recovery_metrics.success_rate = recovery_metrics.n_successful_frequencies / F;
    recovery_metrics.success = true;
    
catch ME
    recovery_metrics.error = ME.message;
end

end

function error_analysis = analyze_error_components(pipeline_data)
% Analyze different sources of error in the pipeline

error_analysis = struct();

try
    ground_truth = pipeline_data.ground_truth;
    module8_results = pipeline_data.module8_results;
    
    F = length(ground_truth.Omega_true);
    
    % Component 1: Transformation accuracy
    transformation_errors = [];
    for omega = 1:F
        if isfield(module8_results, 'recoloring_quality') && ...
           omega <= length(module8_results.recoloring_quality)
            quality = module8_results.recoloring_quality{omega};
            if isfield(quality, 'scaling_metrics') && isfield(quality.scaling_metrics, 'success') && quality.scaling_metrics.success
                diag_scaling_error = quality.scaling_metrics.diagonal_scaling_error;
                transformation_errors(end+1) = diag_scaling_error;
            end
        end
    end
    
    % Component 2: Numerical precision
    numerical_errors = [];
    for omega = 1:F
        if isfield(module8_results, 'recoloring_quality') && ...
           omega <= length(module8_results.recoloring_quality)
            quality = module8_results.recoloring_quality{omega};
            if isfield(quality, 'hermitian_error')
                numerical_errors(end+1) = quality.hermitian_error;
            end
        end
    end
    
    % Component 3: Cycle consistency
    cycle_consistency_errors = [];
    for omega = 1:F
        if isfield(module8_results, 'recoloring_quality') && ...
           omega <= length(module8_results.recoloring_quality)
            quality = module8_results.recoloring_quality{omega};
            if isfield(quality, 'inv_error')
                cycle_consistency_errors(end+1) = quality.inv_error;
            end
        end
    end
    
    % Summarize components
    error_analysis.component_names = {'transformation', 'numerical', 'cycle_consistency'};
    error_analysis.component_magnitudes = [
        median_safe(transformation_errors), ...
        median_safe(numerical_errors), ...
        median_safe(cycle_consistency_errors) ...
    ];
    
    % Identify dominant error source
    [~, dominant_idx] = max(error_analysis.component_magnitudes);
    error_analysis.dominant_source = error_analysis.component_names{dominant_idx};
    
    % Quality scores (higher is better)
    error_analysis.transformation_quality = 1 - min(median_safe(transformation_errors) / 0.1, 1);
    error_analysis.numerical_quality = 1 - min(median_safe(numerical_errors) / 1e-10, 1);
    error_analysis.consistency_quality = 1 - min(median_safe(cycle_consistency_errors) / 0.01, 1);
    
catch ME
    error_analysis.error = ME.message;
end

end

function overall_assessment = generate_overall_assessment(demo_results)
% Generate overall assessment and recommendations

overall_assessment = struct();

% Count successful demos
successful_demos = [
    safe_get_field(demo_results, 'demo1', 'success', false), ...
    safe_get_field(demo_results, 'demo2', 'success', false), ...
    safe_get_field(demo_results, 'demo3', 'success', false), ...
    safe_get_field(demo_results, 'demo4', 'success', false) ...
];

success_count = sum(successful_demos);
overall_assessment.success_count = success_count;
overall_assessment.success_rate = success_count / 4;

% Generate summary
if success_count == 4
    overall_assessment.summary = 'All demos successful - Module 8 fully functional';
    overall_assessment.ready_for_production = true;
elseif success_count >= 3
    overall_assessment.summary = 'Most demos successful - Minor issues detected';
    overall_assessment.ready_for_production = true;
elseif success_count >= 2
    overall_assessment.summary = 'Partial success - Significant issues need attention';
    overall_assessment.ready_for_production = false;
else
    overall_assessment.summary = 'Major issues detected - Module needs debugging';
    overall_assessment.ready_for_production = false;
end

% Generate recommendations
recommendations = {};

if ~successful_demos(1)
    recommendations{end+1} = 'Fix basic pipeline integration issues';
end
if ~successful_demos(2)
    recommendations{end+1} = 'Address complex matrix handling problems';
end
if ~successful_demos(3)
    recommendations{end+1} = 'Investigate performance scaling issues';
end
if ~successful_demos(4)
    recommendations{end+1} = 'Review error analysis and numerical stability';
end

if isempty(recommendations)
    recommendations{end+1} = 'Module 8 is ready for integration with complete pipeline';
    recommendations{end+1} = 'Consider running full pipeline tests with Modules 2-5';
    recommendations{end+1} = 'Monitor performance in production environments';
else
    recommendations{end+1} = 'Run individual unit tests to identify specific issues';
    recommendations{end+1} = 'Check MATLAB path and required dependencies';
end

overall_assessment.recommendations = recommendations;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualization Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function create_module8_comprehensive_visualization(demo_results)
% Create comprehensive visualization of Module 8 results

figure('Name', 'Module 8 Comprehensive Results', 'Position', [100, 100, 1600, 1200]);

% Plot 1: Recovery Error Comparison
subplot(2, 3, 1);
if safe_get_field(demo_results, 'demo1', 'success', false)
    recovery_data = demo_results.demo1.recovery_metrics;
    if isfield(recovery_data, 'median_recovery_error')
        bar([recovery_data.median_recovery_error, recovery_data.mean_recovery_error, recovery_data.max_recovery_error]);
        set(gca, 'XTickLabel', {'Median', 'Mean', 'Max'});
        ylabel('Recovery Error');
        title('Recovery Error Statistics');
        grid on;
    end
else
    text(0.5, 0.5, 'Demo 1 Failed', 'HorizontalAlignment', 'center');
    title('Recovery Error (N/A)');
end

% Plot 2: Complex vs Real Data Performance
subplot(2, 3, 2);
demo1_success = safe_get_field(demo_results, 'demo1', 'success', false);
demo2_success = safe_get_field(demo_results, 'demo2', 'success', false);

performance_comparison = [double(demo1_success), double(demo2_success)];
bar(performance_comparison);
set(gca, 'XTickLabel', {'Real Data', 'Complex Data'});
ylabel('Success (0/1)');
title('Data Type Performance');
ylim([0, 1.2]);
grid on;

% Plot 3: Problem Size Scaling
subplot(2, 3, 3);
if safe_get_field(demo_results, 'demo3', 'success', false)
    perf_data = demo_results.demo3.performance_data;
    if isfield(perf_data, 'processing_times') && ~isempty(perf_data.processing_times)
        plot(perf_data.matrix_operations, perf_data.processing_times, 'o-', 'LineWidth', 2);
        xlabel('Matrix Operations (p^2 F)');
        ylabel('Processing Time (s)');
        title('Performance Scaling');
        grid on;
    end
else
    text(0.5, 0.5, 'Demo 3 Failed', 'HorizontalAlignment', 'center');
    title('Performance Scaling (N/A)');
end

% Plot 4: Error Analysis Components
subplot(2, 3, 4);
if safe_get_field(demo_results, 'demo4', 'success', false)
    error_data = demo_results.demo4;
    if isfield(error_data, 'error_components') && isfield(error_data, 'error_magnitudes')
        bar(error_data.error_magnitudes);
        set(gca, 'XTickLabel', error_data.error_components);
        ylabel('Error Magnitude');
        title('Error Component Analysis');
        set(gca, 'YScale', 'log'); % MATLAB way to set log scale
        grid on;
    end
else
    text(0.5, 0.5, 'Demo 4 Failed', 'HorizontalAlignment', 'center');
    title('Error Analysis (N/A)');
end

% Plot 5: Overall Quality Assessment
subplot(2, 3, 5);
quality_scores = [];
quality_labels = {};

if safe_get_field(demo_results, 'demo1', 'success', false)
    quality_scores(end+1) = 1.0; %#ok<AGROW>
else
    quality_scores(end+1) = 0.0; %#ok<AGROW>
end
quality_labels{end+1} = 'Basic'; %#ok<AGROW>

if safe_get_field(demo_results, 'demo2', 'success', false)
    quality_scores(end+1) = 1.0; %#ok<AGROW>
else
    quality_scores(end+1) = 0.0; %#ok<AGROW>
end
quality_labels{end+1} = 'Complex'; %#ok<AGROW>

if safe_get_field(demo_results, 'demo3', 'success', false)
    quality_scores(end+1) = 1.0; %#ok<AGROW>
else
    quality_scores(end+1) = 0.0; %#ok<AGROW>
end
quality_labels{end+1} = 'Scaling'; %#ok<AGROW>

if safe_get_field(demo_results, 'demo4', 'success', false)
    if isfield(demo_results.demo4, 'overall_quality_score')
        quality_scores(end+1) = demo_results.demo4.overall_quality_score; %#ok<AGROW>
    else
        quality_scores(end+1) = 1.0; %#ok<AGROW>
    end
else
    quality_scores(end+1) = 0.0; %#ok<AGROW>
end
quality_labels{end+1} = 'Quality'; %#ok<AGROW>

bar(quality_scores);
set(gca, 'XTickLabel', quality_labels);
ylabel('Quality Score');
title('Demo Quality Assessment');
ylim([0, 1.2]);
grid on;

% Plot 6: Summary Status
subplot(2, 3, 6);
axis off;

% Create status text
status_text = {};
status_text{end+1} = 'Module 8 Demonstration Results';
status_text{end+1} = '================================';
status_text{end+1} = '';

success_rate = demo_results.overall_assessment.success_rate;
status_text{end+1} = sprintf('Success Rate: %.1f%% (%d/4)', success_rate * 100, demo_results.overall_assessment.success_count);
status_text{end+1} = '';
status_text{end+1} = sprintf('Status: %s', demo_results.overall_assessment.summary);
status_text{end+1} = '';

if demo_results.overall_assessment.ready_for_production
    status_text{end+1} = '✓ Ready for Production';
    text_color = 'green';
else
    status_text{end+1} = '⚠ Needs Attention';
    text_color = [0.8500 0.3250 0.0980]; % orange-ish
end

text(0.05, 0.95, status_text, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
     'FontSize', 10, 'FontName', 'FixedWidth', 'Color', text_color);

sgtitle('Module 8 Recoloring - Comprehensive Demo Results', 'FontSize', 16, 'FontWeight', 'bold');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Utility Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = safe_get_field(struct_data, field1, field2, default_value)
% Safely get nested struct field with default

if isfield(struct_data, field1) && isfield(struct_data.(field1), field2)
    result = struct_data.(field1).(field2);
else
    result = default_value;
end

end

function median_val = median_safe(data_array)
% Compute median safely, handling empty arrays

if isempty(data_array)
    median_val = NaN;
else
    median_val = median(data_array);
end

end

function str_out = logical_to_string(logical_val)
% Convert logical to readable string

if logical_val
    str_out = 'YES';
else
    str_out = 'NO';
end

end

function can_plot = can_make_plots()
% Check if plotting is available

can_plot = ~isempty(get(groot, 'Children')) || usejava('desktop');

end
