function demo_results = demo_module2_estep()
% DEMO_MODULE2_ESTEP - Demonstration of Module 2 E-step computation (Full Source Model)
%
% Returns:
%   demo_results - structure with results, metrics, and (optional) figures
%
% Highlights:
%   - Uses Module7 simulation to create realistic EEG-like data
%   - Runs the full Module2 E-step pipeline
%   - Verifies key mathematical properties (corrected: RTF is a contraction)
%   - Performs parameter sensitivity and performance scaling analyses
%
% NOTE:
%   Plotting is automatically skipped in headless environments.

    rng(2025); % reproducibility across runs
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
        % Realistic EEG-like parameters
        demo_params = struct();
        demo_params.n_sensors = 32;               % EEG sensors
        demo_params.n_sources = 50;               % Cortical sources
        demo_params.n_frequencies = 6;            % e.g. alpha band bins
        demo_params.frequency_range = [8, 12];    % alpha band
        demo_params.snr_db = 15;                  % in dB
        demo_params.edge_density = 0.25;          % graph connectivity

        fprintf('Creating Module7 simulation data:\n');
        fprintf('  - Sensors: %d\n', demo_params.n_sensors);
        fprintf('  - Sources: %d\n', demo_params.n_sources);
        fprintf('  - Frequencies: %d (%.1f - %.1f Hz)\n', ...
                demo_params.n_frequencies, demo_params.frequency_range(1), demo_params.frequency_range(2));
        fprintf('  - SNR: %d dB\n', demo_params.snr_db);
        fprintf('  - Edge density: %.2f\n', demo_params.edge_density);

        % Generate Module7 data and pack input for Module2
        [input_data, ground_truth] = create_module7_eeg_data(demo_params);
        demo_results.demo1.input_params = demo_params;
        demo_results.demo1.ground_truth = ground_truth;

        % Run Module2 E-step
        fprintf('\nRunning E-step computation...\n');
        estep_params = struct('verbose', true, 'regularization_factor', 1e-6);
        t0 = tic;
        estep_results = module2_estep_main(input_data, estep_params);
        computation_time = toc(t0);

        demo_results.demo1.estep_results = estep_results;
        demo_results.demo1.computation_time = computation_time;

        if estep_results.success
            fprintf('\n✓ E-step computation successful!\n');
            fprintf('  - Total time: %.3f seconds\n', computation_time);
            fprintf('  - Successful frequencies: %d/%d\n', ...
                    estep_results.computation_stats.successful_frequencies, demo_params.n_frequencies);

            % Comprehensive visualization
            if can_make_plots(), create_demo1_visualization(estep_results, ground_truth, demo_params); end

            % Precision quality analysis
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
        if isfield(demo_results,'demo1') && demo_results.demo1.success
            fprintf('Analyzing sensitivity to key parameters...\n');

            reg_factors  = [1e-8, 1e-6, 1e-4, 1e-2];
            noise_levels = [0.01, 0.05, 0.10, 0.20];

            sensitivity_results = analyze_parameter_sensitivity(input_data, reg_factors, noise_levels);
            demo_results.demo2.sensitivity_results = sensitivity_results;

            if can_make_plots(), create_sensitivity_visualization(sensitivity_results); end

            fprintf('Parameter sensitivity analysis completed:\n');
            fprintf('  - Optimal regularization: %.0e\n', sensitivity_results.optimal_regularization);
            fprintf('  - Noise robustness range tested: %.3f - %.3f\n', min(noise_levels), max(noise_levels));

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

    %% Demo 3: Mathematical Property Verification (corrected)
    fprintf('\n\nDemo 3: Mathematical Property Verification\n');
    fprintf('-----------------------------------------\n');

    try
        fprintf('Verifying mathematical properties...\n');
        math_verification = verify_mathematical_properties();
        demo_results.demo3.math_verification = math_verification;

        if can_make_plots(), create_math_property_visualization(math_verification); end

        fprintf('Summary:\n');
        fprintf('  - Complementary identity: %s (error: %.2e)\n', ...
                pass_fail_string(math_verification.complementary.passed), math_verification.complementary.error);
        fprintf('  - Contraction (whitened spectrum): %s (min=%.2e, max=%.2e)\n', ...
                pass_fail_string(math_verification.contraction.passed), ...
                math_verification.contraction.min_eval, math_verification.contraction.max_eval);
        fprintf('  - Uncertainty reduction: %s\n', ...
                pass_fail_string(math_verification.uncertainty_reduction.passed));

        demo_results.demo3.success = all([ ...
            math_verification.complementary.passed, ...
            math_verification.contraction.passed, ...
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
        problem_sizes = [12, 16, 24, 32];  % different numbers of sensors
        performance_data = analyze_performance_scaling(problem_sizes);
        demo_results.demo4.performance_data = performance_data;

        if can_make_plots(), create_performance_visualization(performance_data); end

        if numel(performance_data.successful_sizes) >= 3
            scaling_analysis = analyze_computational_scaling(performance_data);
            demo_results.demo4.scaling_analysis = scaling_analysis;

            fprintf('Performance analysis results:\n');
            fprintf('  - Scaling exponent: %.2f (R^2=%.2f)\n', ...
                    scaling_analysis.exponent, scaling_analysis.r_squared);
            fprintf('  - Interpretation: %s\n', scaling_analysis.interpretation);
            fprintf('  - Largest tested size: %d sensors\n', max(problem_sizes));
        end

        demo_results.demo4.success = numel(performance_data.successful_sizes) >= 3;

    catch ME
        fprintf('\n✗ Demo 4 failed: %s\n', ME.message);
        demo_results.demo4.success = false;
        demo_results.demo4.error_message = ME.message;
    end

    %% Demo 5: Ground Truth Comparison
    fprintf('\n\nDemo 5: Ground Truth Comparison Analysis\n');
    fprintf('---------------------------------------\n');

    try
        if isfield(demo_results,'demo1') && demo_results.demo1.success
            fprintf('Comparing E-step results with Module7 ground truth...\n');
            comparison_analysis = compare_with_ground_truth(estep_results, ground_truth);
            demo_results.demo5.comparison_analysis = comparison_analysis;

            if can_make_plots(), create_ground_truth_visualization(estep_results, ground_truth, comparison_analysis); end

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

    flags = [safe_flag(demo_results,'demo1'), safe_flag(demo_results,'demo2'), ...
             safe_flag(demo_results,'demo3'), safe_flag(demo_results,'demo4'), ...
             safe_flag(demo_results,'demo5')];

    overall_success_rate = sum(flags) / numel(flags);
    demo_results.overall_success_rate = overall_success_rate;
    demo_results.overall_success = overall_success_rate >= 0.8;

    fprintf('Individual demo results:\n');
    fprintf('  1. Basic E-step with Module7: %s\n', success_symbol(flags(1)));
    fprintf('  2. Parameter sensitivity: %s\n',   success_symbol(flags(2)));
    fprintf('  3. Mathematical properties: %s\n', success_symbol(flags(3)));
    fprintf('  4. Performance analysis: %s\n',    success_symbol(flags(4)));
    fprintf('  5. Ground truth comparison: %s\n', success_symbol(flags(5)));

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

%% ======================= Helper Functions =======================

function [input_data, ground_truth] = create_module7_eeg_data(params)
% Create EEG-like Module7 data and Module2 inputs

    fprintf('Generating Module7 simulation data...\n');
    [true_precision, true_covariance, empirical_covariance, sim_params] = ...
        module7_simulation_improved_complex( ...
            'n_nodes',   params.n_sensors, ...
            'n_freq',    params.n_frequencies, ...
            'n_samples', 150, ...
            'graph_type','random', ...
            'edge_density', params.edge_density, ...
            'sigma_coef', 0.5);

    % Leadfield
    p = params.n_sensors;
    n = params.n_sources;
    L = randn(p, n) / sqrt(n); % normalization helps conditioning

    % Frequency grid
    frequencies = linspace(params.frequency_range(1), params.frequency_range(2), params.n_frequencies);

    % Pack Module2 input
    input_data = struct();
    input_data.leadfield_matrix = L;
    input_data.empirical_covariances = empirical_covariance;
    input_data.source_prior_covariances = repmat({eye(n)*0.3}, params.n_frequencies, 1);
    input_data.frequencies = frequencies;

    % Noise covariance from SNR
    snr_linear   = 10^(params.snr_db/10);
    noise_variance = 1 / max(1e-12, snr_linear);
    input_data.noise_covariance = eye(p) * noise_variance;

    % Ground truth
    ground_truth = struct();
    ground_truth.true_precision_matrices  = true_precision;
    ground_truth.true_covariance_matrices = true_covariance;
    ground_truth.empirical_covariances    = empirical_covariance;
    ground_truth.simulation_parameters     = sim_params;
    ground_truth.leadfield_matrix         = L;

    fprintf('Module7 simulation completed:\n');
    fprintf('  - Generated %d precision matrices\n', numel(true_precision));
    fprintf('  - Graph edges: %d\n', sim_params.n_edges);
    fprintf('  - Complex matrices: %s\n', logical_to_string(any(cellfun(@(x) ~isreal(x), true_precision))));
end

function create_demo1_visualization(estep_results, ground_truth, demo_params)
% Comprehensive visualization for Demo 1 (skips if headless)
    if ~can_make_plots(), return; end

    figure('Name','Demo 1: E-step Results Overview','Position',[50,50,1400,1000]);
    F = numel(estep_results.initial_precision_matrices);

    % Processing times (first four steps)
    subplot(3,4,1);
    if isfield(estep_results.computation_stats, 'processing_times')
        bar(estep_results.computation_stats.processing_times(:,1:4));
        title('Processing Times by Step');
        xlabel('Frequency'); ylabel('Time (s)');
        legend({'DSTF','Posterior','Residual TF','Residual Cov'}, 'Location','best');
    else
        axis off; text(0.1,0.5,'No timing stats','FontWeight','bold');
    end

    % Condition numbers (only if available)
    subplot(3,4,2);
    if isfield(estep_results.computation_stats, 'condition_numbers')
        semilogy(estep_results.computation_stats.condition_numbers);
        title('Condition Numbers'); xlabel('Frequency'); ylabel('Condition');
        legend({'Sensor Cov','Post Source','Residual'}, 'Location','best');
    else
        axis off; text(0.1,0.5,'No condition stats','FontWeight','bold');
    end

    % True vs estimated precision (f=1)
    if ~isempty(ground_truth.true_precision_matrices{1}) && ~isempty(estep_results.initial_precision_matrices{1})
        subplot(3,4,3); imagesc(abs(ground_truth.true_precision_matrices{1})); colorbar; title('True |Ω|');
        subplot(3,4,4); imagesc(abs(estep_results.initial_precision_matrices{1})); colorbar; title('Estimated |Ω̂|');

        subplot(3,4,5);
        errM = abs(ground_truth.true_precision_matrices{1} - estep_results.initial_precision_matrices{1});
        imagesc(errM); colorbar; title('|Ω - Ω̂|');

        subplot(3,4,6);
        tv = ground_truth.true_precision_matrices{1}(:);
        ev = estep_results.initial_precision_matrices{1}(:);
        scatter(real(tv), real(ev), 10, 'filled', 'MarkerFaceAlpha', 0.6);
        hold on; xl = xlim; plot(xl, xl, 'r--');
        xlabel('True'); ylabel('Estimated'); title('Value Correlation');
        try
            c = corr(real(tv), real(ev), 'rows','complete'); txt = sprintf('r = %.3f', c);
        catch, txt = 'r = NaN'; end
        text(0.1, 0.9, txt, 'Units', 'normalized');
    end

    % DSTF magnitude
    subplot(3,4,7);
    if ~isempty(estep_results.transfer_functions{1})
        imagesc(abs(estep_results.transfer_functions{1})); colorbar; title('|T_{jv}|');
        xlabel('Sensors'); ylabel('Sources');
    end

    % Residual covariance eigenvalues
    subplot(3,4,8);
    if ~isempty(estep_results.residual_covariances{1})
        ev = sort(real(eig(estep_results.residual_covariances{1})), 'descend');
        semilogy(ev, 'o-'); title('Residual Covariance Eigenvalues'); xlabel('Index'); ylabel('Value');
    end

    % Sparsity match (f=1)
    subplot(3,4,9);
    if ~isempty(ground_truth.true_precision_matrices{1}) && ~isempty(estep_results.initial_precision_matrices{1})
        thr = 0.01 * max(abs(ground_truth.true_precision_matrices{1}(:)));
        tp = abs(ground_truth.true_precision_matrices{1}) > thr;
        ep = abs(estep_results.initial_precision_matrices{1}) > thr;
        M  = tp + 2*ep; % 1:true-only, 2:est-only, 3:both
        imagesc(M); colormap(gca, [1 1 1; 1 0 0; 0 1 0; 1 1 0]);
        title('Sparsity Pattern Match');
    end

    % Frequency consistency of one entry
    subplot(3,4,10);
    if F > 1
        v = zeros(F,1);
        for f = 1:F
            if ~isempty(estep_results.initial_precision_matrices{f})
                v(f) = real(estep_results.initial_precision_matrices{f}(2,3));
            end
        end
        plot(1:F, v, 'o-','LineWidth',2); title('Ω̂(2,3) vs Frequency'); xlabel('f'); ylabel('Value');
    end

    % Eigenvalue comparison (f=1)
    subplot(3,4,11);
    if ~isempty(ground_truth.true_precision_matrices{1}) && ~isempty(estep_results.initial_precision_matrices{1})
        e1 = sort(real(eig(ground_truth.true_precision_matrices{1})),'descend');
        e2 = sort(real(eig(estep_results.initial_precision_matrices{1})),'descend');
        semilogy(e1, 'b-o','DisplayName','True'); hold on;
        semilogy(e2, 'r--s','DisplayName','Estimated'); legend; title('Precision Eigenvalues');
        xlabel('Index'); ylabel('Value');
    end

    % Success summary
    subplot(3,4,12);
    success_rate = estep_results.computation_stats.successful_frequencies / max(1,F);
    pie([success_rate, 1-success_rate], {'Success','Failed'});
    title(sprintf('Success: %.0f%%', success_rate*100));

    sgtitle(sprintf('E-step Results: %d sensors, %d sources, %d frequencies', ...
                   demo_params.n_sensors, demo_params.n_sources, demo_params.n_frequencies));
end

function precision_analysis = analyze_precision_quality(estep_results, ground_truth)
% Basic structural/metric quality analysis for the initial precision matrices
    F = numel(estep_results.initial_precision_matrices);
    precision_analysis = struct();

    rec  = zeros(F,1);
    spars = zeros(F,1);
    corrv = zeros(F,1);

    for f = 1:F
        Om = estep_results.initial_precision_matrices{f};
        Tr = ground_truth.true_precision_matrices{f};
        if ~isempty(Om) && ~isempty(Tr)
            thr = 0.01 * max(abs(Tr(:)));
            tp  = abs(Tr) > thr; ep = abs(Om) > thr;

            true_edges = sum(tp(:));
            if true_edges > 0
                rec(f) = sum(tp(:) & ep(:)) / true_edges; % recovery rate
            end
            spars(f) = sum(tp(:) == ep(:)) / numel(tp);   % sparsity match
            try
                c = corr(real(Tr(:)), real(Om(:)), 'rows','complete');
            catch
                c = NaN;
            end
            corrv(f) = c;
        end
    end

    precision_analysis.recovery_rates        = rec;
    precision_analysis.sparsity_preservations= spars;
    precision_analysis.correlations          = corrv;
    precision_analysis.avg_recovery_rate     = mean(rec(rec>0));
    precision_analysis.sparsity_preservation = mean(spars(spars>0));
    precision_analysis.avg_correlation       = mean(corrv(isfinite(corrv)));
end

function sensitivity_results = analyze_parameter_sensitivity(input_data, reg_factors, noise_levels)
% Evaluate sensitivity w.r.t regularization and noise parameters
    sensitivity_results = struct();

    % --- Regularization sweep ---
    R = numel(reg_factors);
    reg_quality = zeros(R, 3); % [success_rate, logcond(Ω), logdet(Ω)]
    for i = 1:R
        results = module2_estep_main(input_data, struct('regularization_factor', reg_factors(i), 'verbose', false));
        if results.success
            reg_quality(i,1) = results.computation_stats.successful_frequencies / numel(input_data.empirical_covariances);
            if ~isempty(results.initial_precision_matrices{1})
                Omega = (results.initial_precision_matrices{1}+results.initial_precision_matrices{1}')/2;
                ev = real(eig(Omega)); ev(ev<=0) = eps;
                reg_quality(i,2) = log(max(ev)) - log(min(ev)); % log condition
                reg_quality(i,3) = safe_logdet_chol(Omega);     % robust logdet
            end
        else
            reg_quality(i,:) = [0, NaN, NaN];
        end
    end
    sensitivity_results.regularization_factors = reg_factors;
    sensitivity_results.regularization_quality = reg_quality;

    % Pick "optimal" reg: highest success rate; tie-break by smallest logcond
    best = find(reg_quality(:,1) == max(reg_quality(:,1)));
    [~,k] = min(reg_quality(best,2)); sensitivity_results.optimal_regularization = reg_factors(best(k));

    % --- Noise sweep ---
    N = numel(noise_levels);
    noise_quality = zeros(N, 2); % [success_rate, quality_score]
    for i = 1:N
        inp = input_data;
        p   = size(input_data.noise_covariance,1);
        inp.noise_covariance = eye(p) * noise_levels(i);
        results = module2_estep_main(inp, struct('regularization_factor', sensitivity_results.optimal_regularization, 'verbose', false));
        if results.success
            noise_quality(i,1) = results.computation_stats.successful_frequencies / numel(inp.empirical_covariances);
            if ~isempty(results.initial_precision_matrices{1})
                Om = (results.initial_precision_matrices{1}+results.initial_precision_matrices{1}')/2;
                ev = real(eig(Om)); ev(ev<=0) = eps;
                % quality proxy: inverse of log-condition (bounded to [0,1] roughly)
                q = 1./(1 + (log(max(ev)) - log(min(ev))));
                noise_quality(i,2) = q;
            end
        end
    end
    sensitivity_results.noise_levels = noise_levels;
    sensitivity_results.noise_quality = noise_quality;
end

function create_sensitivity_visualization(S)
% Visual summaries for parameter sensitivity (skips if headless)
    if ~can_make_plots(), return; end

    figure('Name','Parameter Sensitivity Analysis','Position',[150,100,1200,600]);

    % Regularization
    subplot(2,3,1);
    semilogx(S.regularization_factors, S.regularization_quality(:,1)*100, 'o-');
    title('Success Rate vs Regularization'); xlabel('Reg'); ylabel('Success (%)'); grid on;

    subplot(2,3,2);
    semilogx(S.regularization_factors, exp(S.regularization_quality(:,2)), 's-');
    set(gca,'YScale','log'); title('Condition Number vs Regularization'); xlabel('Reg'); ylabel('Cond'); grid on;

    subplot(2,3,3);
    semilogx(S.regularization_factors, S.regularization_quality(:,3), '^-');
    title('logdet(Ω) vs Regularization'); xlabel('Reg'); ylabel('logdet'); grid on;

    % Noise
    subplot(2,3,4);
    semilogx(S.noise_levels, S.noise_quality(:,1)*100, 'o-');
    title('Success Rate vs Noise'); xlabel('Noise Var'); ylabel('Success (%)'); grid on;

    subplot(2,3,5);
    semilogx(S.noise_levels, S.noise_quality(:,2), 's-');
    title('Quality vs Noise'); xlabel('Noise Var'); ylabel('Quality Score'); grid on;

    subplot(2,3,6);
    % Placeholder for combined landscape (optional)
    axis off; text(0.1,0.5,'Parameter landscape omitted','FontWeight','bold');
end

function math_verification = verify_mathematical_properties()
% Verify corrected math: complementary identity, contraction, uncertainty reduction
    math_verification = struct();

    % Small controlled setup for stability
    p = 8; n = 12;
    L = randn(p, n);
    Sigma_jj = eye(n) * 0.5;
    Sigma_xi = eye(p) * 0.1;

    try
        % DSTF
        T_jv = module2_dstf_computation(L, Sigma_jj, Sigma_xi, struct('regularization_factor',1e-8,'verbose',false));
        % RTF (standard path)
        T_xi_v = module2_residual_transfer_function(T_jv, L, 'validate_properties', false, 'verbose', false);

        % 1) Complementary identity: T_xi_v + L*T_jv ≈ I
        I_p = eye(p);
        comp_err = rel_fro(T_xi_v + L*T_jv - I_p, I_p);
        math_verification.complementary.error  = comp_err;
        math_verification.complementary.passed = comp_err < 1e-10;

        % 2) Contraction (whitened): eig( S^{-1/2} T S^{1/2} ) ∈ [0,1]
        [Rc,flag] = chol((Sigma_xi+Sigma_xi')/2);
        if flag ~= 0, error('Sigma_xi must be PD for whitening'); end
        M = (Rc \ T_xi_v) * Rc'; M = (M+M')/2;
        eM = eig(M);
        math_verification.contraction.min_eval = min(real(eM));
        math_verification.contraction.max_eval = max(real(eM));
        math_verification.contraction.passed   = (math_verification.contraction.max_eval <= 1+1e-10) ...
                                               && (math_verification.contraction.min_eval >= -1e-10);

        % 3) Uncertainty reduction: Sigma_post ≼ Sigma_prior
        Sigma_post = module2_posterior_source_covariance(Sigma_jj, L, Sigma_xi, ...
                        struct('regularization_factor',1e-8,'verbose',false));
        D = (Sigma_jj + Sigma_jj')/2 - (Sigma_post + Sigma_post')/2;
        ev = eig(D);
        math_verification.uncertainty_reduction.min_eigenval = min(real(ev));
        math_verification.uncertainty_reduction.passed = all(real(ev) >= -1e-10);
        math_verification.uncertainty_reduction.trace_reduction = trace(D) / max(1,trace(Sigma_jj));

        % store matrices for plots
        math_verification.test_matrices.T_xi_v = T_xi_v;
        math_verification.test_matrices.identity_sum = T_xi_v + L*T_jv;
        math_verification.test_matrices.uncertainty_reduction = D;

    catch ME
        fprintf('Warning: Math verification issue: %s\n', ME.message);
        math_verification.complementary = struct('error',NaN,'passed',false);
        math_verification.contraction   = struct('min_eval',NaN,'max_eval',NaN,'passed',false);
        math_verification.uncertainty_reduction = struct('min_eigenval',NaN,'passed',false,'trace_reduction',NaN);
        math_verification.error_encountered = ME.message;
    end
end

function create_math_property_visualization(V)
% Visualize corrected properties (skips if headless)
    if ~can_make_plots(), return; end

    figure('Name','Mathematical Property Verification','Position',[250,100,1200,700]);

    % Complementary identity
    subplot(2,3,1);
    imagesc(real(V.test_matrices.identity_sum)); colorbar;
    title(sprintf('T_{\\xi v} + L T_{jv} (err=%.2e)', V.complementary.error));

    % Residual transform (real part)
    subplot(2,3,2);
    imagesc(real(V.test_matrices.T_xi_v)); colorbar; title('T_{\\xi v} (Real)');

    % Uncertainty reduction matrix
    subplot(2,3,3);
    imagesc(real(V.test_matrices.uncertainty_reduction)); colorbar;
    title('Uncertainty Reduction (Prior - Post)');

    % Spectrum placeholder for contraction (values printed in title)
    subplot(2,3,4); axis off;
    text(0.1,0.6, sprintf('Contraction: min=%.2e, max=%.2e', V.contraction.min_eval, V.contraction.max_eval), 'FontSize',12);

    % Eigenvalues of uncertainty reduction
    subplot(2,3,5);
    ev = eig(V.test_matrices.uncertainty_reduction);
    bar(real(ev)); title('Uncertainty Reduction Eigenvalues'); xlabel('Index'); ylabel('Eigenvalue');

    % Pass/Fail summary
    subplot(2,3,6); props = {'Complementary','Contraction','Uncertainty Red.'};
    passed = [V.complementary.passed, V.contraction.passed, V.uncertainty_reduction.passed];
    bar(passed); set(gca,'XTickLabel',props); ylim([0,1.2]); ylabel('Pass (1) / Fail (0)'); title('Property Summary');
end

function performance_data = analyze_performance_scaling(problem_sizes)
% Measure runtime/memory across problem sizes (best-effort memory stats)
    performance_data = struct();
    performance_data.problem_sizes   = problem_sizes;
    performance_data.computation_times = zeros(size(problem_sizes));
    performance_data.success_flags   = false(size(problem_sizes));
    performance_data.memory_usage    = zeros(size(problem_sizes));

    for i = 1:numel(problem_sizes)
        p = problem_sizes(i);
        n = round(p * 1.2);
        fprintf('  Testing problem size: %d sensors, %d sources... ', p, n);

        try
            [~, emp_true_cov, emp_cov, sim_params] = ...
                module7_simulation_improved_complex('n_nodes', p, 'n_freq', 3, ...
                                                    'n_samples', 100, 'graph_type','random', ...
                                                    'edge_density', 0.2); %#ok<ASGLU>

            input_data = struct();
            input_data.leadfield_matrix = randn(p, n);
            input_data.empirical_covariances = emp_cov;
            input_data.source_prior_covariances = repmat({eye(n)*0.4}, numel(emp_cov), 1);
            input_data.frequencies = 1:numel(emp_cov);
            input_data.noise_covariance = eye(p) * 0.1;

            % Timing & memory (best effort, platform dependent)
            try, mem_before = get_mem_used_mb(); catch, mem_before = NaN; end
            t0 = tic;
            results = module2_estep_main(input_data, struct('verbose', false));
            elapsed = toc(t0);
            try, mem_after = get_mem_used_mb(); catch, mem_after = NaN; end

            performance_data.computation_times(i) = elapsed;
            performance_data.success_flags(i)     = results.success;
            performance_data.memory_usage(i)      = mem_after - mem_before;

            fprintf('%.3f sec, %s MB\n', elapsed, num2str(performance_data.memory_usage(i)));

        catch ME
            fprintf('FAILED (%s)\n', ME.message);
            performance_data.computation_times(i) = NaN;
            performance_data.success_flags(i)     = false;
            performance_data.memory_usage(i)      = NaN;
        end
    end

    performance_data.successful_sizes  = problem_sizes(performance_data.success_flags);
    performance_data.successful_times  = performance_data.computation_times(performance_data.success_flags);
    performance_data.successful_memory = performance_data.memory_usage(performance_data.success_flags);
end

function create_performance_visualization(P)
% Visualize performance scaling (skips if headless)
    if ~can_make_plots(), return; end
    figure('Name','Performance Analysis','Position',[350,100,1200,600]);

    % Time scaling
    subplot(2,3,1);
    if ~isempty(P.successful_sizes)
        loglog(P.successful_sizes, P.successful_times, 'o-','LineWidth',2);
        title('Computation Time Scaling'); xlabel('Sensors'); ylabel('Time (s)'); grid on;
    else
        axis off; text(0.1,0.5,'No successful runs','FontWeight','bold');
    end

    % Memory scaling
    subplot(2,3,2);
    if ~isempty(P.successful_sizes)
        loglog(P.successful_sizes, abs(P.successful_memory), 's-','LineWidth',2);
        title('Memory Usage Scaling'); xlabel('Sensors'); ylabel('Memory (MB)'); grid on;
    else
        axis off;
    end

    % Success vs size
    subplot(2,3,3);
    bar(P.problem_sizes, P.success_flags);
    title('Success by Size'); xlabel('Sensors'); ylabel('Success (1/0)');

    % Time per element
    subplot(2,3,4);
    if ~isempty(P.successful_sizes)
        elems = P.successful_sizes.^2;
        tpe = P.successful_times ./ elems * 1e6; % μs per element
        semilogx(P.successful_sizes, tpe, '^-'); title('Time per Matrix Element');
        xlabel('Sensors'); ylabel('μs/element');
    else
        axis off;
    end

    % Efficiency
    subplot(2,3,5);
    if numel(P.successful_sizes) > 1
        eff = P.successful_times(1) ./ P.successful_times;
        ideal = (P.successful_sizes(1).^2) ./ (P.successful_sizes.^2);
        semilogx(P.successful_sizes, eff, 'o-','DisplayName','Actual'); hold on;
        semilogx(P.successful_sizes, ideal, '--','DisplayName','Ideal O(n^2)'); legend;
        title('Computational Efficiency'); xlabel('Sensors'); ylabel('Rel Efficiency');
    else
        axis off;
    end

    % Size guideline (illustrative)
    subplot(2,3,6);
    if ~isempty(P.successful_sizes)
        under10 = P.successful_sizes(P.successful_times < 10);
        rec_size = max(under10(:), [], 'omitnan');
        max_size = max(P.successful_sizes);
        bar([max(rec_size,0), max_size], [1, 0.5]);
        set(gca, 'XTickLabel', {'Recommended','Maximum Tested'});
        title('Problem Size Guidelines'); ylabel('Suitability');
        text(1, 1.05, sprintf('%d sensors', max(rec_size,0)), 'HorizontalAlignment','center');
        text(2, 0.55, sprintf('%d sensors', max_size), 'HorizontalAlignment','center');
    else
        axis off; text(0.1,0.5,'No guideline available','FontWeight','bold');
    end
end

function scaling_analysis = analyze_computational_scaling(P)
% Fit power law time ~ a * size^b (log-log regression)
    if numel(P.successful_sizes) < 3
        scaling_analysis = struct('exponent', NaN, 'interpretation', 'Insufficient data', 'r_squared', NaN);
        return;
    end
    xs = log(P.successful_sizes(:)); ys = log(P.successful_times(:));
    C = polyfit(xs, ys, 1); b = C(1);
    scaling_analysis = struct();
    scaling_analysis.exponent = b;
    scaling_analysis.fit_coefficients = C;
    scaling_analysis.r_squared = corr(xs, ys)^2;

    if b < 1.5,     txt = 'Better than quadratic scaling (excellent)';
    elseif b < 2.5, txt = 'Approximately quadratic scaling (good)';
    elseif b < 3.5, txt = 'Between quadratic and cubic scaling (acceptable)';
    else            txt = 'Worse than cubic scaling (concerning)';
    end
    scaling_analysis.interpretation = txt;
end

function comparison_analysis = compare_with_ground_truth(estep_results, ground_truth)
% Compare estimated vs true precision across frequencies
    F = numel(estep_results.initial_precision_matrices);
    edge_accuracies = zeros(F,1);
    value_correlations = zeros(F,1);
    sparsity_matches = zeros(F,1);

    for f = 1:F
        Om = estep_results.initial_precision_matrices{f};
        Tr = ground_truth.true_precision_matrices{f};
        if ~isempty(Om) && ~isempty(Tr)
            thr = 0.01 * max(abs(Tr(:)));
            tp = abs(Tr) > thr; ep = abs(Om) > thr;

            tp_count = sum(tp(:));
            tn_count = sum(~tp(:));
            tp_pos   = sum(tp(:) & ep(:));
            fp_pos   = sum(~tp(:) & ep(:));
            fn_pos   = sum(tp(:) & ~ep(:));
            tn_pos   = sum(~tp(:) & ~ep(:));

            if (tp_count) > 0, sens = tp_pos / (tp_count); else, sens = 1; end
            if (tn_count) > 0, spec = tn_pos / (tn_count); else, spec = 1; end

            edge_accuracies(f) = (sens + spec)/2; % balanced accuracy

            try
                value_correlations(f) = corr(real(Tr(:)), real(Om(:)), 'rows','complete');
            catch
                value_correlations(f) = NaN;
            end

            sparsity_matches(f) = sum(tp(:) == ep(:)) / numel(tp);
        end
    end

    comparison_analysis = struct();
    comparison_analysis.edge_detection_accuracy = mean(edge_accuracies(edge_accuracies>0));
    comparison_analysis.value_correlation       = mean(value_correlations(isfinite(value_correlations)));
    comparison_analysis.sparsity_match          = mean(sparsity_matches(sparsity_matches>0));

    comparison_analysis.overall_quality = mean([ ...
        comparison_analysis.edge_detection_accuracy, ...
        abs(comparison_analysis.value_correlation), ...
        comparison_analysis.sparsity_match]);

    comparison_analysis.per_frequency.edge_accuracies   = edge_accuracies;
    comparison_analysis.per_frequency.value_correlations= value_correlations;
    comparison_analysis.per_frequency.sparsity_matches  = sparsity_matches;
end

function create_ground_truth_visualization(estep_results, ground_truth, C)
% Visualize comparison to ground truth (skips if headless)
    if ~can_make_plots(), return; end

    figure('Name','Ground Truth Comparison','Position',[450,100,1400,800]);
    F = numel(estep_results.initial_precision_matrices);

    % Metrics across frequencies
    subplot(2,4,1);
    plot(1:F, C.per_frequency.edge_accuracies, 'o-','LineWidth',2);
    title('Edge Detection Accuracy'); xlabel('Frequency'); ylabel('Accuracy'); ylim([0,1]);

    subplot(2,4,2);
    plot(1:F, C.per_frequency.value_correlations, 's-','LineWidth',2);
    title('Value Correlation'); xlabel('Frequency'); ylabel('Correlation'); ylim([-1,1]);

    subplot(2,4,3);
    plot(1:F, C.per_frequency.sparsity_matches, '^-','LineWidth',2);
    title('Sparsity Pattern Match'); xlabel('Frequency'); ylabel('Match'); ylim([0,1]);

    subplot(2,4,4);
    bar([C.edge_detection_accuracy, abs(C.value_correlation), C.sparsity_match]);
    set(gca,'XTickLabel',{'Edge Acc','Value Corr','Sparsity'}); ylim([0,1]); title('Overall Metrics');

    % Detailed for first frequency
    if ~isempty(ground_truth.true_precision_matrices{1}) && ~isempty(estep_results.initial_precision_matrices{1})
        Tr = ground_truth.true_precision_matrices{1};
        Om = estep_results.initial_precision_matrices{1};

        subplot(2,4,5); imagesc(abs(Tr)); colorbar; title('True |Ω|');
        subplot(2,4,6); imagesc(abs(Om)); colorbar; title('Estimated |Ω̂|');

        % ROC-like curve (threshold sweep on estimated matrix)
        subplot(2,4,7);
        thr_true = 0.01 * max(abs(Tr(:)));
        tp = abs(Tr) > thr_true;
        thrs = linspace(0, max(abs(Om(:))), 100);
        sens = zeros(size(thrs)); spec = zeros(size(thrs));
        for i = 1:numel(thrs)
            ep = abs(Om) > thrs(i);
            tp_pos = sum(tp(:) & ep(:));
            fp_pos = sum(~tp(:) & ep(:));
            fn_pos = sum(tp(:) & ~ep(:));
            tn_pos = sum(~tp(:) & ~ep(:));
            sens(i) = tp_pos / max(1,(tp_pos+fn_pos));
            spec(i) = tn_pos / max(1,(tn_pos+fp_pos));
        end
        plot(1-spec, sens, 'LineWidth',2); xlabel('FPR'); ylabel('TPR'); title('ROC Curve');

        % Value scatter
        subplot(2,4,8);
        tv = Tr(:); ev = Om(:);
        scatter(real(tv), real(ev), 10, 'filled', 'MarkerFaceAlpha', 0.6);
        hold on; xl = xlim; plot(xl, xl, 'r--');
        xlabel('True'); ylabel('Estimated'); title(sprintf('Value Corr: %.3f', C.value_correlation));
    end

    sgtitle(sprintf('Ground Truth Comparison - Overall Quality: %.2f', C.overall_quality));
end

function generate_usage_recommendations(D)
% Print practical guidance derived from demo results
    fprintf('\nUsage Recommendations:\n');
    if safe_flag(D,'demo2')
        fprintf('  - Optimal regularization factor: %.0e\n', D.demo2.sensitivity_results.optimal_regularization);
        fprintf('  - Apply mild ridge (ε≈1e-6~1e-4) if conditioning degrades.\n');
    end
    if safe_flag(D,'demo3')
        fprintf('  - Mathematical checks passed (complementary + contraction + uncertainty).\n');
    end
    if safe_flag(D,'demo4') && isfield(D.demo4,'scaling_analysis')
        fprintf('  - Computational complexity: %s\n', D.demo4.scaling_analysis.interpretation);
        if isfield(D.demo4,'performance_data') && ~isempty(D.demo4.performance_data.successful_sizes)
            fprintf('  - Suitable up to ~%d sensors on standard hardware.\n', ...
                    max(D.demo4.performance_data.successful_sizes));
        end
    end
    if safe_flag(D,'demo5')
        fprintf('  - Ground-truth comparison is satisfactory for network analysis (edges/values/sparsity).\n');
    end
    fprintf('  - Monitor conditioning of sensor/residual covariances; use regularization when needed.\n');
end

function generate_troubleshooting_guide(D)
% Print troubleshooting guidance for failed demos
    fprintf('\nTroubleshooting Guide:\n');
    if ~safe_flag(D,'demo1')
        fprintf('  - Check Module7 data generation and leadfield dimensions.\n');
    end
    if ~safe_flag(D,'demo3')
        fprintf('  - If math properties fail, increase regularization and verify covariance Hermitianity.\n');
    end
    if ~safe_flag(D,'demo4')
        fprintf('  - Performance issues: reduce problem size or enable parallel BLAS if available.\n');
    end
    if ~safe_flag(D,'demo5')
        fprintf('  - Poor ground truth recovery: adjust SNR/edge density or source priors.\n');
    end
    fprintf('  - Inspect figures for hints (if available).\n');
end

function str = logical_to_string(v)
    if v, str = 'YES'; else, str = 'NO'; end
end

function str = success_symbol(flag)
    if flag, str = '✓'; else, str = '✗'; end
end

function str = pass_fail_string(flag)
    if flag, str = 'PASS'; else, str = 'FAIL'; end
end

function ok = safe_flag(S, fieldname)
    ok = isfield(S, fieldname) && isfield(S.(fieldname),'success') && S.(fieldname).success;
end

%% ---------- Small Utilities ----------
function tf = can_make_plots()
    tf = usejava('jvm') && feature('ShowFigureWindows');
end

function e = rel_fro(A, B)
% Relative Frobenius norm ||A||_F / ||B||_F (absolute if B omitted)
    if nargin < 2 || isempty(B), e = norm(A,'fro'); else, e = norm(A,'fro') / max(1, norm(B,'fro')); end
end

function v = safe_logdet_chol(A)
% Robust log(det(A)) using Cholesky with small jitter
    A = (A + A')/2; p = size(A,1);
    base = 1e-12 * max(1, trace(A)/max(1,p));
    for k = 0:6
        [R,flag] = chol(A + (base*10^k)*eye(p,class(A)));
        if flag == 0, v = 2*sum(log(abs(diag(R)))); return; end
    end
    ev = real(eig(A)); ev(ev<=0) = eps; v = sum(log(ev));
end

function mb = get_mem_used_mb()
% Best-effort memory usage (platform dependent). Falls back to NaN.
    mb = NaN;
    try
        ms = memory; mb = ms.MemUsedMATLAB/1e6;
    catch
        try
            s = feature('memstats'); mb = s.PeakMemUsed/1e6;
        catch
            % leave NaN
        end
    end
end
