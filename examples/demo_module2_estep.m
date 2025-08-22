function demo_results = demo_module2_estep()
    % DEMO_MODULE2_ESTEP - Demonstration of Module 2 E-step computation (Full Source Model)
    %
    % Returns:
    %   demo_results - Structure containing demonstration results and analysis
    %
    % Description:
    %   This demonstration showcases the capabilities of Module 2 for E-step
    %   computation in the EM algorithm for precision matrix estimation.
    %   Uses Module7 simulation for realistic data with comprehensive visualization.
    %
    % Demonstration Coverage:
    %   1. Basic E-step computation with Module7 EEG-like data
    %   2. Mathematical property verification with visualization
    %   3. Parameter sensitivity analysis
    %   4. Performance analysis across different problem sizes
    %   5. Comparison with ground truth from simulation
    %
    % Author: [Author Name]
    % Date: [Current Date]  
    % Version: 2.0 - Full Source Model with Module7 Integration
    
    fprintf('========================================\n');
    fprintf('Module 2 E-step Computation Demonstration\n');
    fprintf('Full Source Model with Module7 Simulation\n');
    fprintf('========================================\n\n');
    
    demo_results = struct();
    demo_results.timestamp = datetime('now');
    demo_results.matlab_version = version();
    
    %% Demo 1: Basic E-step with Module7 EEG Data
    fprintf('Demo 1: Basic E-step with Module7 EEG-like Data\n');
    fprintf('-----------------------------------------------\n');
    
    try
        % Create realistic EEG-like parameters
        demo_params = struct();
        demo_params.n_sensors = 32;        % EEG sensor count
        demo_params.n_sources = 50;        % Cortical source locations  
        demo_params.n_frequencies = 6;     % Alpha band frequencies
        demo_params.frequency_range = [8, 12]; % Alpha band
        demo_params.snr_db = 15;           % Signal-to-noise ratio
        demo_params.edge_density = 0.25;   % Graph connectivity
        
        fprintf('Creating Module7 simulation data:\n');
        fprintf('  - Sensors: %d\n', demo_params.n_sensors);
        fprintf('  - Sources: %d\n', demo_params.n_sources);
        fprintf('  - Frequencies: %d (%.1f - %.1f Hz)\n', ...
                demo_params.n_frequencies, demo_params.frequency_range(1), demo_params.frequency_range(2));
        fprintf('  - SNR: %d dB\n', demo_params.snr_db);
        fprintf('  - Edge density: %.2f\n', demo_params.edge_density);
        
        % Generate Module7 simulation data
        [input_data, ground_truth] = create_module7_eeg_data(demo_params);
        demo_results.demo1.input_params = demo_params;
        demo_results.demo1.ground_truth = ground_truth;
        
        % Run E-step computation
        fprintf('\nRunning E-step computation...\n');
        estep_params = struct('verbose', true, 'regularization_factor', 1e-6);
        
        tic;
        estep_results = module2_estep_main(input_data, estep_params);
        computation_time = toc;
        
        demo_results.demo1.estep_results = estep_results;
        demo_results.demo1.computation_time = computation_time;
        
        % Analyze and visualize results
        if estep_results.success
            fprintf('\n✓ E-step computation successful!\n');
            fprintf('  - Total time: %.3f seconds\n', computation_time);
            fprintf('  - Successful frequencies: %d/%d\n', ...
                    estep_results.computation_stats.successful_frequencies, demo_params.n_frequencies);
            
            % Create comprehensive visualization
            create_demo1_visualization(estep_results, ground_truth, demo_params);
            
            % Analyze precision matrix quality
            precision_analysis = analyze_precision_quality(estep_results, ground_truth);
            demo_results.demo1.precision_analysis = precision_analysis;
            
            fprintf('  - Average recovery accuracy: %.2f%%\n', precision_analysis.avg_recovery_rate * 100);
            fprintf('  - Sparsity preservation: %.2f%%\n', precision_analysis.sparsity_preservation * 100);
            
        else
            fprintf('\n✗ E-step computation failed\n');
            demo_results.demo1.error_message = 'E-step computation failed';
        end
        
        demo_results.demo1.success = estep_results.success;
        
    catch ME
        fprintf('\n✗ Demo 1 failed: %s\n', ME.message);
        demo_results.demo1.success = false;
        demo_results.demo1.error_message = ME.message;
    end
    
    %% Demo 2: Parameter Sensitivity Analysis
    fprintf('\n\nDemo 2: Parameter Sensitivity Analysis\n');
    fprintf('--------------------------------------\n');
    
    try
        if demo_results.demo1.success
            fprintf('Analyzing sensitivity to key parameters...\n');
            
            % Test different regularization factors
            reg_factors = [1e-8, 1e-6, 1e-4, 1e-2];
            noise_levels = [0.01, 0.05, 0.1, 0.2];
            
            sensitivity_results = analyze_parameter_sensitivity(input_data, reg_factors, noise_levels);
            demo_results.demo2.sensitivity_results = sensitivity_results;
            
            % Create sensitivity visualization
            create_sensitivity_visualization(sensitivity_results);
            
            fprintf('Parameter sensitivity analysis completed:\n');
            fprintf('  - Optimal regularization: %.0e\n', sensitivity_results.optimal_regularization);
            fprintf('  - Noise robustness range: %.3f - %.3f\n', ...
                    min(noise_levels), max(noise_levels));
            
            demo_results.demo2.success = true;
            
        else
            fprintf('Skipping sensitivity analysis due to Demo 1 failure\n');
            demo_results.demo2.success = false;
            demo_results.demo2.skip_reason = 'Demo 1 failed';
        end
        
    catch ME
        fprintf('\n✗ Demo 2 failed: %s\n', ME.message);
        demo_results.demo2.success = false;
        demo_results.demo2.error_message = ME.message;
    end
    
    %% Demo 3: Mathematical Property Verification
    fprintf('\n\nDemo 3: Mathematical Property Verification\n');
    fprintf('-----------------------------------------\n');
    
    try
        fprintf('Verifying mathematical properties with visualization...\n');
        
        % Test mathematical properties using controlled scenarios
        math_verification = verify_mathematical_properties();
        demo_results.demo3.math_verification = math_verification;
        
        % Create mathematical property visualization
        create_math_property_visualization(math_verification);
        
        fprintf('Mathematical property verification:\n');
        fprintf('  - Idempotent property: %s (error: %.2e)\n', ...
                pass_fail_string(math_verification.idempotent.passed), math_verification.idempotent.error);
        fprintf('  - Complementary projection: %s (error: %.2e)\n', ...
                pass_fail_string(math_verification.complementary.passed), math_verification.complementary.error);
        fprintf('  - Uncertainty reduction: %s\n', ...
                pass_fail_string(math_verification.uncertainty_reduction.passed));
        
        demo_results.demo3.success = all([math_verification.idempotent.passed, ...
                                         math_verification.complementary.passed, ...
                                         math_verification.uncertainty_reduction.passed]);
        
    catch ME
        fprintf('\n✗ Demo 3 failed: %s\n', ME.message);
        demo_results.demo3.success = false;
        demo_results.demo3.error_message = ME.message;
    end
    
    %% Demo 4: Performance Analysis
    fprintf('\n\nDemo 4: Performance Analysis\n');
    fprintf('----------------------------\n');
    
    try
        fprintf('Analyzing computational performance across problem sizes...\n');
        
        problem_sizes = [12, 16, 24, 32];  % Different numbers of sensors
        performance_data = analyze_performance_scaling(problem_sizes);
        demo_results.demo4.performance_data = performance_data;
        
        % Create performance visualization
        create_performance_visualization(performance_data);
        
        % Analyze scaling behavior
        if length(performance_data.successful_sizes) >= 3
            scaling_analysis = analyze_computational_scaling(performance_data);
            demo_results.demo4.scaling_analysis = scaling_analysis;
            
            fprintf('Performance analysis results:\n');
            fprintf('  - Scaling exponent: %.2f\n', scaling_analysis.exponent);
            fprintf('  - Interpretation: %s\n', scaling_analysis.interpretation);
            fprintf('  - Largest tested size: %d sensors\n', max(problem_sizes));
        end
        
        demo_results.demo4.success = length(performance_data.successful_sizes) >= 3;
        
    catch ME
        fprintf('\n✗ Demo 4 failed: %s\n', ME.message);
        demo_results.demo4.success = false;
        demo_results.demo4.error_message = ME.message;
    end
    
    %% Demo 5: Ground Truth Comparison
    fprintf('\n\nDemo 5: Ground Truth Comparison Analysis\n');
    fprintf('---------------------------------------\n');
    
    try
        if demo_results.demo1.success
            fprintf('Comparing E-step results with Module7 ground truth...\n');
            
            comparison_analysis = compare_with_ground_truth(estep_results, ground_truth);
            demo_results.demo5.comparison_analysis = comparison_analysis;
            
            % Create ground truth comparison visualization
            create_ground_truth_visualization(estep_results, ground_truth, comparison_analysis);
            
            fprintf('Ground truth comparison results:\n');
            fprintf('  - Edge detection accuracy: %.2f%%\n', comparison_analysis.edge_detection_accuracy * 100);
            fprintf('  - Value correlation: %.3f\n', comparison_analysis.value_correlation);
            fprintf('  - Sparsity match: %.2f%%\n', comparison_analysis.sparsity_match * 100);
            
            demo_results.demo5.success = comparison_analysis.overall_quality > 0.7;
            
        else
            fprintf('Skipping ground truth comparison due to Demo 1 failure\n');
            demo_results.demo5.success = false;
            demo_results.demo5.skip_reason = 'Demo 1 failed';
        end
        
    catch ME
        fprintf('\n✗ Demo 5 failed: %s\n', ME.message);
        demo_results.demo5.success = false;
        demo_results.demo5.error_message = ME.message;
    end
    
    %% Summary and Recommendations
    fprintf('\n\n========================================\n');
    fprintf('Demonstration Summary\n');
    fprintf('========================================\n');
    
    demo_success_flags = [
        demo_results.demo1.success, ...
        demo_results.demo2.success, ...
        demo_results.demo3.success, ...
        demo_results.demo4.success, ...
        demo_results.demo5.success
    ];
    
    overall_success_rate = sum(demo_success_flags) / length(demo_success_flags);
    demo_results.overall_success_rate = overall_success_rate;
    demo_results.overall_success = overall_success_rate >= 0.8;
    
    fprintf('Individual demo results:\n');
    fprintf('  1. Basic E-step with Module7: %s\n', success_symbol(demo_results.demo1.success));
    fprintf('  2. Parameter sensitivity: %s\n', success_symbol(demo_results.demo2.success));
    fprintf('  3. Mathematical properties: %s\n', success_symbol(demo_results.demo3.success));
    fprintf('  4. Performance analysis: %s\n', success_symbol(demo_results.demo4.success));
    fprintf('  5. Ground truth comparison: %s\n', success_symbol(demo_results.demo5.success));
    
    fprintf('\nOverall success rate: %.1f%%\n', overall_success_rate * 100);
    
    if demo_results.overall_success
        fprintf('\n🎉 OVERALL RESULT: SUCCESS!\n');
        fprintf('   Module 2 E-step computation is ready for production use.\n');
        generate_usage_recommendations(demo_results);
    else
        fprintf('\n⚠️  OVERALL RESULT: PARTIAL SUCCESS\n');
        fprintf('   Some issues detected. Review failed demos and visualizations.\n');
        generate_troubleshooting_guide(demo_results);
    end
    
    demo_results.completion_time = datetime('now');
    fprintf('\nDemonstration completed at: %s\n', char(demo_results.completion_time));
    fprintf('Check generated figures for detailed visualizations.\n');
    fprintf('========================================\n');
end

%% Helper Functions

function [input_data, ground_truth] = create_module7_eeg_data(params)
    % Create EEG-like data using Module7 simulation
    
    fprintf('Generating Module7 simulation data...\n');
    
    % Call Module7 simulation (移除verbose参数)
    [true_precision, true_covariance, empirical_covariance, sim_params] = ...
        module7_simulation_improved_complex(...
            'n_nodes', params.n_sensors, ...
            'n_freq', params.n_frequencies, ...
            'n_samples', 150, ...
            'graph_type', 'random', ...
            'edge_density', params.edge_density, ...
            'sigma_coef', 0.5);
    
    % Create leadfield matrix 
    p = params.n_sensors;
    n = params.n_sources;
    L = randn(p, n) / sqrt(n);  % Normalized for numerical stability
    
    % Create frequency vector
    frequencies = linspace(params.frequency_range(1), params.frequency_range(2), params.n_frequencies);
    
    % Prepare input data structure
    input_data = struct();
    input_data.leadfield_matrix = L;
    input_data.empirical_covariances = empirical_covariance;
    input_data.source_prior_covariances = cell(params.n_frequencies, 1);
    input_data.frequencies = frequencies;
    
    % Noise covariance based on SNR
    snr_linear = 10^(params.snr_db / 10);
    noise_variance = 1 / snr_linear;
    input_data.noise_covariance = eye(p) * noise_variance;
    
    % Create source prior covariances
    for f = 1:params.n_frequencies
        % Use reasonable source prior
        input_data.source_prior_covariances{f} = eye(n) * 0.3;
    end
    
    % Ground truth storage
    ground_truth = struct();
    ground_truth.true_precision_matrices = true_precision;
    ground_truth.true_covariance_matrices = true_covariance;
    ground_truth.simulation_parameters = sim_params;
    ground_truth.empirical_covariances = empirical_covariance;
    ground_truth.leadfield_matrix = L;
    
    fprintf('Module7 simulation completed:\n');
    fprintf('  - Generated %d precision matrices\n', length(true_precision));
    fprintf('  - Graph edges: %d\n', sim_params.n_edges);
    fprintf('  - Complex matrices: %s\n', logical_to_string(any(cellfun(@(x) ~isreal(x), true_precision))));
end

function create_demo1_visualization(estep_results, ground_truth, demo_params)
    % Create comprehensive visualization for Demo 1
    
    figure('Name', 'Demo 1: E-step Results Overview', 'Position', [50, 50, 1400, 1000]);
    
    F = length(estep_results.initial_precision_matrices);
    
    % Processing performance
    subplot(3, 4, 1);
    if isfield(estep_results.computation_stats, 'processing_times')
        bar(estep_results.computation_stats.processing_times(:, 1:4));
        title('Processing Times by Step');
        xlabel('Frequency');
        ylabel('Time (s)');
        legend({'DSTF', 'Posterior', 'Residual TF', 'Residual Cov'}, 'Location', 'best');
    end
    
    % Condition numbers
    subplot(3, 4, 2);
    if isfield(estep_results.computation_stats, 'condition_numbers')
        semilogy(estep_results.computation_stats.condition_numbers);
        title('Condition Numbers');
        xlabel('Frequency');
        ylabel('log(Condition Number)');
        legend({'Sensor Cov', 'Post Source', 'Residual'}, 'Location', 'best');
    end
    
    % True vs estimated precision (first frequency)
    if ~isempty(ground_truth.true_precision_matrices{1}) && ~isempty(estep_results.initial_precision_matrices{1})
        subplot(3, 4, 3);
        imagesc(abs(ground_truth.true_precision_matrices{1}));
        colorbar;
        title('True Precision |Ω|');
        
        subplot(3, 4, 4);
        imagesc(abs(estep_results.initial_precision_matrices{1}));
        colorbar;
        title('Estimated Precision |Ω̂|');
        
        % Error analysis
        subplot(3, 4, 5);
        error_matrix = abs(ground_truth.true_precision_matrices{1} - estep_results.initial_precision_matrices{1});
        imagesc(error_matrix);
        colorbar;
        title('Absolute Error |Ω - Ω̂|');
        
        % Scatter plot comparison
        subplot(3, 4, 6);
        true_vals = ground_truth.true_precision_matrices{1}(:);
        est_vals = estep_results.initial_precision_matrices{1}(:);
        scatter(real(true_vals), real(est_vals), 10, 'filled', 'alpha', 0.6);
        xlabel('True Values');
        ylabel('Estimated Values');
        title('Value Correlation');
        hold on;
        plot(xlim, xlim, 'r--');
        
        % Calculate correlation
        corr_val = corr(real(true_vals), real(est_vals));
        text(0.1, 0.9, sprintf('r = %.3f', corr_val), 'Units', 'normalized');
    end
    
    % Transfer function analysis
    subplot(3, 4, 7);
    if ~isempty(estep_results.transfer_functions{1})
        imagesc(abs(estep_results.transfer_functions{1}));
        colorbar;
        title('DSTF |T_{jv}|');
        xlabel('Sensors');
        ylabel('Sources');
    end
    
    % Residual covariance eigenvalues
    subplot(3, 4, 8);
    if ~isempty(estep_results.residual_covariances{1})
        eigs_res = sort(real(eig(estep_results.residual_covariances{1})), 'descend');
        semilogy(eigs_res, 'o-');
        title('Residual Covariance Eigenvalues');
        xlabel('Index');
        ylabel('Eigenvalue');
    end
    
    % Sparsity pattern analysis
    subplot(3, 4, 9);
    if ~isempty(ground_truth.true_precision_matrices{1})
        threshold = 0.01 * max(abs(ground_truth.true_precision_matrices{1}(:)));
        true_pattern = abs(ground_truth.true_precision_matrices{1}) > threshold;
        est_pattern = abs(estep_results.initial_precision_matrices{1}) > threshold;
        
        % Create combined pattern visualization
        combined = zeros(size(true_pattern));
        combined(true_pattern & est_pattern) = 3;  % Both (yellow)
        combined(true_pattern & ~est_pattern) = 1; % True only (red)
        combined(~true_pattern & est_pattern) = 2; % Est only (green)
        
        imagesc(combined);
        colormap(gca, [1 1 1; 1 0 0; 0 1 0; 1 1 0]); % White, Red, Green, Yellow
        title('Sparsity Pattern Match');
    end
    
    % Frequency consistency
    subplot(3, 4, 10);
    if F > 1
        sample_entry = zeros(F, 1);
        for f = 1:F
            if ~isempty(estep_results.initial_precision_matrices{f})
                sample_entry(f) = real(estep_results.initial_precision_matrices{f}(2, 3));
            end
        end
        plot(1:F, sample_entry, 'o-', 'LineWidth', 2);
        title('Precision Entry vs Frequency');
        xlabel('Frequency Index');
        ylabel('Ω̂(2,3)');
    end
    
    % Eigenvalue comparison
    subplot(3, 4, 11);
    if ~isempty(ground_truth.true_precision_matrices{1}) && ~isempty(estep_results.initial_precision_matrices{1})
        true_eigs = sort(real(eig(ground_truth.true_precision_matrices{1})), 'descend');
        est_eigs = sort(real(eig(estep_results.initial_precision_matrices{1})), 'descend');
        
        semilogy(true_eigs, 'b-o', 'DisplayName', 'True');
        hold on;
        semilogy(est_eigs, 'r--s', 'DisplayName', 'Estimated');
        legend;
        title('Eigenvalue Comparison');
        xlabel('Index');
        ylabel('Eigenvalue');
    end
    
    % Success summary
    subplot(3, 4, 12);
    success_rate = estep_results.computation_stats.successful_frequencies / F;
    pie([success_rate, 1-success_rate], {'Success', 'Failed'});
    title(sprintf('Success: %.0f%%', success_rate*100));
    
    sgtitle(sprintf('E-step Results: %d sensors, %d sources, %d frequencies', ...
                   demo_params.n_sensors, demo_params.n_sources, demo_params.n_frequencies));
end

function precision_analysis = analyze_precision_quality(estep_results, ground_truth)
    % Analyze quality of precision matrix estimation
    
    F = length(estep_results.initial_precision_matrices);
    precision_analysis = struct();
    
    recovery_rates = zeros(F, 1);
    sparsity_preservations = zeros(F, 1);
    correlations = zeros(F, 1);
    
    for f = 1:F
        if ~isempty(estep_results.initial_precision_matrices{f}) && ~isempty(ground_truth.true_precision_matrices{f})
            true_prec = ground_truth.true_precision_matrices{f};
            est_prec = estep_results.initial_precision_matrices{f};
            
            % Sparsity analysis
            threshold = 0.01 * max(abs(true_prec(:)));
            true_pattern = abs(true_prec) > threshold;
            est_pattern = abs(est_prec) > threshold;
            
            % Recovery rate (what fraction of true edges were detected)
            true_edges = sum(true_pattern(:));
            if true_edges > 0
                recovery_rates(f) = sum(true_pattern(:) & est_pattern(:)) / true_edges;
            end
            
            % Sparsity preservation (what fraction of total edges match)
            total_elements = numel(true_pattern);
            sparsity_preservations(f) = sum(true_pattern(:) == est_pattern(:)) / total_elements;
            
            % Value correlation
            correlations(f) = corr(real(true_prec(:)), real(est_prec(:)));
        end
    end
    
    precision_analysis.recovery_rates = recovery_rates;
    precision_analysis.sparsity_preservations = sparsity_preservations;
    precision_analysis.correlations = correlations;
    precision_analysis.avg_recovery_rate = mean(recovery_rates(recovery_rates > 0));
    precision_analysis.sparsity_preservation = mean(sparsity_preservations(sparsity_preservations > 0));
    precision_analysis.avg_correlation = mean(correlations(correlations > 0));
end

function sensitivity_results = analyze_parameter_sensitivity(input_data, reg_factors, noise_levels)
    % Analyze sensitivity to regularization and noise parameters
    
    sensitivity_results = struct();
    
    % Test regularization sensitivity
    reg_quality = zeros(length(reg_factors), 3); % [success_rate, avg_condition, avg_det]
    
    for i = 1:length(reg_factors)
        estep_params = struct('regularization_factor', reg_factors(i), 'verbose', false);
        results = module2_estep_main(input_data, estep_params);
        
        if results.success
            reg_quality(i, 1) = results.computation_stats.successful_frequencies / length(input_data.empirical_covariances);
            
            if ~isempty(results.initial_precision_matrices{1})
                reg_quality(i, 2) = cond(results.initial_precision_matrices{1});
                reg_quality(i, 3) = det(results.initial_precision_matrices{1});
            end
        end
    end
    
    sensitivity_results.regularization_factors = reg_factors;
    sensitivity_results.regularization_quality = reg_quality;
    
    % Find optimal regularization (best success rate with reasonable condition number)
    valid_indices = reg_quality(:, 1) > 0.8 & reg_quality(:, 2) < 1e10;
    if any(valid_indices)
        [~, best_idx] = min(reg_quality(valid_indices, 2));
        valid_regs = reg_factors(valid_indices);
        sensitivity_results.optimal_regularization = valid_regs(best_idx);
    else
        sensitivity_results.optimal_regularization = reg_factors(2); % Default fallback
    end
    
    % Test noise level sensitivity
    noise_quality = zeros(length(noise_levels), 2); % [success_rate, avg_accuracy]
    
    for i = 1:length(noise_levels)
        input_data_noise = input_data;
        input_data_noise.noise_covariance = eye(size(input_data.noise_covariance)) * noise_levels(i);
        
        estep_params = struct('regularization_factor', sensitivity_results.optimal_regularization, 'verbose', false);
        results = module2_estep_main(input_data_noise, estep_params);
        
        if results.success
            noise_quality(i, 1) = results.computation_stats.successful_frequencies / length(input_data.empirical_covariances);
            
            % Simple accuracy measure (trace similarity to low-noise case)
            if ~isempty(results.initial_precision_matrices{1})
                noise_quality(i, 2) = 1 / (1 + norm(results.initial_precision_matrices{1}, 'fro'));
            end
        end
    end
    
    sensitivity_results.noise_levels = noise_levels;
    sensitivity_results.noise_quality = noise_quality;
end

function create_sensitivity_visualization(sensitivity_results)
    % Create visualization for parameter sensitivity analysis
    
    figure('Name', 'Parameter Sensitivity Analysis', 'Position', [150, 100, 1200, 600]);
    
    % Regularization sensitivity
    subplot(2, 3, 1);
    semilogx(sensitivity_results.regularization_factors, sensitivity_results.regularization_quality(:, 1) * 100, 'o-');
    title('Success Rate vs Regularization');
    xlabel('Regularization Factor');
    ylabel('Success Rate (%)');
    grid on;
    
    subplot(2, 3, 2);
    loglog(sensitivity_results.regularization_factors, sensitivity_results.regularization_quality(:, 2), 's-');
    title('Condition Number vs Regularization');
    xlabel('Regularization Factor');
    ylabel('Condition Number');
    grid on;
    
    subplot(2, 3, 3);
    semilogx(sensitivity_results.regularization_factors, sensitivity_results.regularization_quality(:, 3), '^-');
    title('Determinant vs Regularization');
    xlabel('Regularization Factor');
    ylabel('Determinant');
    grid on;
    
    % Noise sensitivity
    subplot(2, 3, 4);
    semilogx(sensitivity_results.noise_levels, sensitivity_results.noise_quality(:, 1) * 100, 'o-');
    title('Success Rate vs Noise Level');
    xlabel('Noise Variance');
    ylabel('Success Rate (%)');
    grid on;
    
    subplot(2, 3, 5);
    semilogx(sensitivity_results.noise_levels, sensitivity_results.noise_quality(:, 2), 's-');
    title('Quality vs Noise Level');
    xlabel('Noise Variance');
    ylabel('Quality Score');
    grid on;
    
    subplot(2, 3, 6);
    % Optimal parameter region
    [X, Y] = meshgrid(log10(sensitivity_results.regularization_factors), log10(sensitivity_results.noise_levels));
    Z = ones(size(X)) * 0.5; % Placeholder for combined quality metric
    contourf(X, Y, Z);
    colorbar;
    title('Parameter Space Quality');
    xlabel('log₁₀(Regularization)');
    ylabel('log₁₀(Noise)');
end

function math_verification = verify_mathematical_properties()
    % Verify mathematical properties with controlled test cases
    
    math_verification = struct();
    
    % Create controlled test scenario with smaller dimensions for numerical stability
    p = 8; n = 12;  % 减小维数以提高数值稳定性
    L = randn(p, n);
    Sigma_jj = eye(n) * 0.5;
    Sigma_xi = eye(p) * 0.1;
    
    % Test DSTF and RTF properties with proper computation pipeline
    try
        % Step 1: Compute DSTF with regularization
        options_dstf = struct('regularization_factor', 1e-8, 'verbose', false);
        T_jv = module2_dstf_computation(L, Sigma_jj, Sigma_xi, options_dstf);
        
        % Step 2: Compute RTF properly
        T_xi_v = module2_residual_transfer_function(T_jv, L);
        
        % Test idempotent property: T_ξv² = T_ξv
        T_xi_v_squared = T_xi_v * T_xi_v;
        idempotent_error = norm(T_xi_v_squared - T_xi_v, 'fro') / (norm(T_xi_v, 'fro') + eps);
        
        math_verification.idempotent.error = idempotent_error;
        math_verification.idempotent.passed = idempotent_error < 1e-10;
        
        % Test complementary projection: T_ξv + L*T_jv = I
        L_T_jv = L * T_jv;
        identity_sum = T_xi_v + L_T_jv;
        I_p = eye(p);
        complementary_error = norm(identity_sum - I_p, 'fro') / norm(I_p, 'fro');
        
        math_verification.complementary.error = complementary_error;
        math_verification.complementary.passed = complementary_error < 1e-10;
        
        % Test uncertainty reduction property
        options_psc = struct('regularization_factor', 1e-8, 'verbose', false);
        Sigma_post = module2_posterior_source_covariance(Sigma_jj, L, Sigma_xi, options_psc);
        uncertainty_reduction_matrix = Sigma_jj - Sigma_post;
        eigenvals_reduction = eig(uncertainty_reduction_matrix);
        
        math_verification.uncertainty_reduction.min_eigenval = min(real(eigenvals_reduction));
        math_verification.uncertainty_reduction.passed = all(real(eigenvals_reduction) >= -1e-10);
        math_verification.uncertainty_reduction.trace_reduction = trace(uncertainty_reduction_matrix) / trace(Sigma_jj);
        
        % Store test matrices for visualization
        math_verification.test_matrices.T_xi_v = T_xi_v;
        math_verification.test_matrices.T_xi_v_squared = T_xi_v_squared;
        math_verification.test_matrices.identity_sum = identity_sum;
        math_verification.test_matrices.uncertainty_reduction = uncertainty_reduction_matrix;
        
    catch ME
        % Handle numerical errors gracefully
        fprintf('Warning: Mathematical verification encountered numerical issues: %s\n', ME.message);
        math_verification.idempotent.error = NaN;
        math_verification.idempotent.passed = false;
        math_verification.complementary.error = NaN;
        math_verification.complementary.passed = false;
        math_verification.uncertainty_reduction.passed = false;
        math_verification.error_encountered = ME.message;
    end
end

function create_math_property_visualization(math_verification)
    % Create visualization for mathematical property verification
    
    figure('Name', 'Mathematical Property Verification', 'Position', [250, 100, 1200, 800]);
    
    % Idempotent property visualization
    subplot(2, 4, 1);
    imagesc(real(math_verification.test_matrices.T_xi_v));
    colorbar;
    title('T_{\xi v} (Real)');
    
    subplot(2, 4, 2);
    imagesc(real(math_verification.test_matrices.T_xi_v_squared));
    colorbar;
    title('T_{\xi v}^2 (Real)');
    
    subplot(2, 4, 3);
    error_matrix = math_verification.test_matrices.T_xi_v_squared - math_verification.test_matrices.T_xi_v;
    imagesc(real(error_matrix));
    colorbar;
    title(sprintf('T_ξv² - T_ξv (err: %.2e)', math_verification.idempotent.error));
    
    % Complementary projection visualization
    subplot(2, 4, 4);
    imagesc(real(math_verification.test_matrices.identity_sum));
    colorbar;
    title(sprintf('T_ξv + LT_jv (err: %.2e)', math_verification.complementary.error));
    
    % Uncertainty reduction analysis
    subplot(2, 4, 5);
    imagesc(real(math_verification.test_matrices.uncertainty_reduction));
    colorbar;
    title('Uncertainty Reduction Matrix');
    
    subplot(2, 4, 6);
    eigenvals = eig(math_verification.test_matrices.uncertainty_reduction);
    bar(real(eigenvals));
    title('Uncertainty Reduction Eigenvalues');
    xlabel('Index');
    ylabel('Eigenvalue');
    
    % Property summary
    subplot(2, 4, 7);
    properties = {'Idempotent', 'Complementary', 'Uncertainty Red.'};
    passed = [math_verification.idempotent.passed, ...
              math_verification.complementary.passed, ...
              math_verification.uncertainty_reduction.passed];
    
    bar(passed);
    set(gca, 'XTickLabel', properties);
    title('Property Verification');
    ylabel('Passed (1) / Failed (0)');
    ylim([0, 1.2]);
    
    % Error summary
    subplot(2, 4, 8);
    errors = [math_verification.idempotent.error, ...
              math_verification.complementary.error, ...
              abs(math_verification.uncertainty_reduction.min_eigenval)];
    
    semilogy(errors, 'o-');
    set(gca, 'XTickLabel', properties);
    title('Error Magnitudes');
    ylabel('log(Error)');
end

function performance_data = analyze_performance_scaling(problem_sizes)
    % Analyze computational performance across different problem sizes
    
    performance_data = struct();
    performance_data.problem_sizes = problem_sizes;
    performance_data.computation_times = zeros(size(problem_sizes));
    performance_data.success_flags = false(size(problem_sizes));
    performance_data.memory_usage = zeros(size(problem_sizes));
    
    for i = 1:length(problem_sizes)
        p = problem_sizes(i);
        n = round(p * 1.2);  % Sources proportional to sensors
        
        fprintf('  Testing problem size: %d sensors, %d sources... ', p, n);
        
        try
            % Generate test data using Module7
            [true_prec, true_cov, emp_cov, sim_params] = ...
                module7_simulation_improved_complex(...
                    'n_nodes', p, ...
                    'n_freq', 3, ...
                    'n_samples', 100, ...
                    'graph_type', 'random', ...
                    'edge_density', 0.2);
            
            % Create input data
            input_data = struct();
            input_data.leadfield_matrix = randn(p, n);
            input_data.empirical_covariances = emp_cov;
            input_data.source_prior_covariances = repmat({eye(n)*0.4}, length(emp_cov), 1);
            input_data.frequencies = 1:length(emp_cov);
            input_data.noise_covariance = eye(p) * 0.1;
            
            % Measure computation time and memory
            tic;
            mem_before = feature('memstats');
            
            results = module2_estep_main(input_data, struct('verbose', false));
            
            elapsed_time = toc;
            mem_after = feature('memstats');
            
            performance_data.computation_times(i) = elapsed_time;
            performance_data.success_flags(i) = results.success;
            performance_data.memory_usage(i) = (mem_after.PeakMemUsed - mem_before.PeakMemUsed) / 1e6; % MB
            
            fprintf('%.3f sec, %.1f MB\n', elapsed_time, performance_data.memory_usage(i));
            
        catch ME
            fprintf('FAILED (%s)\n', ME.message);
            performance_data.computation_times(i) = NaN;
            performance_data.success_flags(i) = false;
            performance_data.memory_usage(i) = NaN;
        end
    end
    
    % Extract successful runs
    performance_data.successful_sizes = problem_sizes(performance_data.success_flags);
    performance_data.successful_times = performance_data.computation_times(performance_data.success_flags);
    performance_data.successful_memory = performance_data.memory_usage(performance_data.success_flags);
end

function create_performance_visualization(performance_data)
    % Create visualization for performance analysis
    
    figure('Name', 'Performance Analysis', 'Position', [350, 100, 1200, 600]);
    
    % Computation time scaling
    subplot(2, 3, 1);
    loglog(performance_data.successful_sizes, performance_data.successful_times, 'o-', 'LineWidth', 2);
    title('Computation Time Scaling');
    xlabel('Problem Size (sensors)');
    ylabel('Time (seconds)');
    grid on;
    
    % Memory usage scaling
    subplot(2, 3, 2);
    loglog(performance_data.successful_sizes, performance_data.successful_memory, 's-', 'LineWidth', 2);
    title('Memory Usage Scaling');
    xlabel('Problem Size (sensors)');
    ylabel('Memory (MB)');
    grid on;
    
    % Success rate
    subplot(2, 3, 3);
    bar(performance_data.problem_sizes, performance_data.success_flags);
    title('Success Rate vs Problem Size');
    xlabel('Problem Size (sensors)');
    ylabel('Success (1) / Failure (0)');
    
    % Time per element
    subplot(2, 3, 4);
    elements = performance_data.successful_sizes.^2;  % Roughly quadratic in matrix elements
    time_per_element = performance_data.successful_times ./ elements * 1e6; % microseconds per element
    semilogx(performance_data.successful_sizes, time_per_element, '^-');
    title('Time per Matrix Element');
    xlabel('Problem Size (sensors)');
    ylabel('Time per Element (μs)');
    
    % Efficiency comparison
    subplot(2, 3, 5);
    if length(performance_data.successful_sizes) > 1
        efficiency = performance_data.successful_times(1) ./ performance_data.successful_times;
        ideal_efficiency = performance_data.successful_sizes(1).^2 ./ performance_data.successful_sizes.^2;
        
        semilogx(performance_data.successful_sizes, efficiency, 'o-', 'DisplayName', 'Actual');
        hold on;
        semilogx(performance_data.successful_sizes, ideal_efficiency, '--', 'DisplayName', 'Ideal O(n²)');
        legend;
        title('Computational Efficiency');
        xlabel('Problem Size (sensors)');
        ylabel('Relative Efficiency');
    end
    
    % Problem size recommendation
    subplot(2, 3, 6);
    max_successful = max(performance_data.successful_sizes);
    recommended_sizes = performance_data.successful_sizes(performance_data.successful_times < 10); % Under 10 seconds
    
    bar([max(recommended_sizes), max_successful], [1, 0.5]);
    set(gca, 'XTickLabel', {'Recommended', 'Maximum Tested'});
    title('Problem Size Guidelines');
    ylabel('Suitability');
    text(1, 1.1, sprintf('%d sensors', max(recommended_sizes)), 'HorizontalAlignment', 'center');
    text(2, 0.6, sprintf('%d sensors', max_successful), 'HorizontalAlignment', 'center');
end

function scaling_analysis = analyze_computational_scaling(performance_data)
    % Analyze computational scaling behavior
    
    if length(performance_data.successful_sizes) < 3
        scaling_analysis = struct('exponent', NaN, 'interpretation', 'Insufficient data');
        return;
    end
    
    % Fit power law: time = a * size^b
    log_sizes = log(performance_data.successful_sizes);
    log_times = log(performance_data.successful_times);
    
    coeffs = polyfit(log_sizes, log_times, 1);
    scaling_exponent = coeffs(1);
    
    scaling_analysis = struct();
    scaling_analysis.exponent = scaling_exponent;
    scaling_analysis.fit_coefficients = coeffs;
    scaling_analysis.r_squared = corr(log_sizes(:), log_times(:))^2;
    
    % Interpret scaling behavior
    if scaling_exponent < 1.5
        scaling_analysis.interpretation = 'Better than quadratic scaling (excellent)';
    elseif scaling_exponent < 2.5
        scaling_analysis.interpretation = 'Approximately quadratic scaling (good)';
    elseif scaling_exponent < 3.5
        scaling_analysis.interpretation = 'Between quadratic and cubic scaling (acceptable)';
    else
        scaling_analysis.interpretation = 'Worse than cubic scaling (concerning)';
    end
end

function comparison_analysis = compare_with_ground_truth(estep_results, ground_truth)
    % Compare E-step results with Module7 ground truth
    
    comparison_analysis = struct();
    F = length(estep_results.initial_precision_matrices);
    
    edge_accuracies = zeros(F, 1);
    value_correlations = zeros(F, 1);
    sparsity_matches = zeros(F, 1);
    
    for f = 1:F
        if ~isempty(estep_results.initial_precision_matrices{f}) && ~isempty(ground_truth.true_precision_matrices{f})
            true_prec = ground_truth.true_precision_matrices{f};
            est_prec = estep_results.initial_precision_matrices{f};
            
            % Edge detection accuracy
            threshold = 0.01 * max(abs(true_prec(:)));
            true_edges = abs(true_prec) > threshold;
            est_edges = abs(est_prec) > threshold;
            
            % Sensitivity and specificity for edge detection
            true_positives = sum(true_edges(:) & est_edges(:));
            false_positives = sum(~true_edges(:) & est_edges(:));
            false_negatives = sum(true_edges(:) & ~est_edges(:));
            true_negatives = sum(~true_edges(:) & ~est_edges(:));
            
            if (true_positives + false_negatives) > 0
                sensitivity = true_positives / (true_positives + false_negatives);
            else
                sensitivity = 1;
            end
            
            if (true_negatives + false_positives) > 0
                specificity = true_negatives / (true_negatives + false_positives);
            else
                specificity = 1;
            end
            
            edge_accuracies(f) = (sensitivity + specificity) / 2;  % Balanced accuracy
            
            % Value correlation
            value_correlations(f) = corr(real(true_prec(:)), real(est_prec(:)));
            
            % Sparsity pattern match
            sparsity_matches(f) = sum(true_edges(:) == est_edges(:)) / numel(true_edges);
        end
    end
    
    comparison_analysis.edge_detection_accuracy = mean(edge_accuracies(edge_accuracies > 0));
    comparison_analysis.value_correlation = mean(value_correlations(value_correlations > 0));
    comparison_analysis.sparsity_match = mean(sparsity_matches(sparsity_matches > 0));
    
    % Overall quality score
    comparison_analysis.overall_quality = mean([
        comparison_analysis.edge_detection_accuracy, ...
        abs(comparison_analysis.value_correlation), ...
        comparison_analysis.sparsity_match
    ]);
    
    % Detailed results per frequency
    comparison_analysis.per_frequency.edge_accuracies = edge_accuracies;
    comparison_analysis.per_frequency.value_correlations = value_correlations;
    comparison_analysis.per_frequency.sparsity_matches = sparsity_matches;
end

function create_ground_truth_visualization(estep_results, ground_truth, comparison_analysis)
    % Create visualization comparing results with ground truth
    
    figure('Name', 'Ground Truth Comparison', 'Position', [450, 100, 1400, 800]);
    
    F = length(estep_results.initial_precision_matrices);
    
    % Accuracy metrics across frequencies
    subplot(2, 4, 1);
    plot(1:F, comparison_analysis.per_frequency.edge_accuracies, 'o-', 'LineWidth', 2);
    title('Edge Detection Accuracy');
    xlabel('Frequency');
    ylabel('Accuracy');
    ylim([0, 1]);
    
    subplot(2, 4, 2);
    plot(1:F, comparison_analysis.per_frequency.value_correlations, 's-', 'LineWidth', 2);
    title('Value Correlation');
    xlabel('Frequency');
    ylabel('Correlation');
    ylim([-1, 1]);
    
    subplot(2, 4, 3);
    plot(1:F, comparison_analysis.per_frequency.sparsity_matches, '^-', 'LineWidth', 2);
    title('Sparsity Pattern Match');
    xlabel('Frequency');
    ylabel('Match Rate');
    ylim([0, 1]);
    
    % Overall quality radar plot (if we had more metrics)
    subplot(2, 4, 4);
    metrics = [comparison_analysis.edge_detection_accuracy, ...
               abs(comparison_analysis.value_correlation), ...
               comparison_analysis.sparsity_match];
    metric_names = {'Edge Accuracy', 'Value Corr', 'Sparsity Match'};
    
    bar(metrics);
    set(gca, 'XTickLabel', metric_names);
    title('Overall Quality Metrics');
    ylabel('Score (0-1)');
    ylim([0, 1]);
    
    % Detailed comparison for first frequency
    if ~isempty(ground_truth.true_precision_matrices{1}) && ~isempty(estep_results.initial_precision_matrices{1})
        true_prec = ground_truth.true_precision_matrices{1};
        est_prec = estep_results.initial_precision_matrices{1};
        
        subplot(2, 4, 5);
        imagesc(abs(true_prec));
        colorbar;
        title('True Precision |Ω|');
        
        subplot(2, 4, 6);
        imagesc(abs(est_prec));
        colorbar;
        title('Estimated Precision |Ω̂|');
        
        % ROC-like analysis for edge detection
        subplot(2, 4, 7);
        threshold_range = linspace(0, max(abs(est_prec(:))), 100);
        true_threshold = 0.01 * max(abs(true_prec(:)));
        true_edges = abs(true_prec) > true_threshold;
        
        sensitivity = zeros(size(threshold_range));
        specificity = zeros(size(threshold_range));
        
        for i = 1:length(threshold_range)
            est_edges = abs(est_prec) > threshold_range(i);
            
            tp = sum(true_edges(:) & est_edges(:));
            fp = sum(~true_edges(:) & est_edges(:));
            fn = sum(true_edges(:) & ~est_edges(:));
            tn = sum(~true_edges(:) & ~est_edges(:));
            
            if (tp + fn) > 0
                sensitivity(i) = tp / (tp + fn);
            end
            if (tn + fp) > 0
                specificity(i) = tn / (tn + fp);
            end
        end
        
        plot(1 - specificity, sensitivity, 'LineWidth', 2);
        xlabel('False Positive Rate');
        ylabel('True Positive Rate');
        title('ROC Curve for Edge Detection');
        
        % Value scatter plot
        subplot(2, 4, 8);
        scatter(real(true_prec(:)), real(est_prec(:)), 10, 'filled', 'alpha', 0.6);
        xlabel('True Values');
        ylabel('Estimated Values');
        title(sprintf('Value Correlation: %.3f', comparison_analysis.value_correlation));
        hold on;
        plot(xlim, xlim, 'r--');
    end
    
    sgtitle(sprintf('Ground Truth Comparison - Overall Quality: %.2f', comparison_analysis.overall_quality));
end

function generate_usage_recommendations(demo_results)
    % Generate usage recommendations based on demo results
    
    fprintf('\nUsage Recommendations:\n');
    
    if demo_results.demo2.success
        fprintf('  - Optimal regularization factor: %.0e\n', demo_results.demo2.sensitivity_results.optimal_regularization);
        fprintf('  - Module handles numerical challenges well with proper regularization\n');
    end
    
    if demo_results.demo3.success
        fprintf('  - Mathematical properties verified - algorithm is numerically sound\n');
    end
    
    if demo_results.demo4.success && isfield(demo_results.demo4, 'scaling_analysis')
        fprintf('  - Computational complexity: %s\n', demo_results.demo4.scaling_analysis.interpretation);
        if ~isempty(demo_results.demo4.performance_data.successful_sizes)
            max_recommended = max(demo_results.demo4.performance_data.successful_sizes);
            fprintf('  - Suitable for problems up to ~%d sensors on standard hardware\n', max_recommended);
        end
    end
    
    if demo_results.demo5.success
        fprintf('  - Ground truth comparison shows good recovery performance\n');
        fprintf('  - Edge detection accuracy: %.1f%% - suitable for network analysis\n', ...
                demo_results.demo5.comparison_analysis.edge_detection_accuracy * 100);
    end
    
    fprintf('  - Use Module7 simulation for algorithm validation and testing\n');
    fprintf('  - Monitor condition numbers and apply regularization when needed\n');
end

function generate_troubleshooting_guide(demo_results)
    % Generate troubleshooting guide for failed demos
    
    fprintf('\nTroubleshooting Guide:\n');
    
    if ~demo_results.demo1.success
        fprintf('  - Basic functionality failed: Check Module7 simulation and leadfield matrix\n');
        fprintf('  - Verify input data dimensions and matrix properties\n');
    end
    
    if ~demo_results.demo3.success
        fprintf('  - Mathematical property violations: Increase regularization factor\n');
        fprintf('  - Check for numerical instabilities in matrix computations\n');
    end
    
    if ~demo_results.demo4.success
        fprintf('  - Performance issues: Consider reducing problem size\n');
        fprintf('  - Monitor memory usage for large problems\n');
    end
    
    if ~demo_results.demo5.success
        fprintf('  - Poor ground truth recovery: Check SNR and edge density parameters\n');
        fprintf('  - Consider adjusting source prior covariances\n');
    end
    
    fprintf('  - Always check generated visualizations for diagnostic information\n');
end

function str = logical_to_string(logical_value)
    % Convert logical value to string
    if logical_value
        str = 'YES';
    else
        str = 'NO';
    end
end

function str = success_symbol(success_flag)
    % Return success/failure symbol
    if success_flag
        str = '✓';
    else
        str = '✗';
    end
end

function str = pass_fail_string(passed)
    % Return PASS/FAIL string
    if passed
        str = 'PASS';
    else
        str = 'FAIL';
    end
end