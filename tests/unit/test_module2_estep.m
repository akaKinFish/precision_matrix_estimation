function test_results = test_module2_estep()
% TEST_MODULE2_ESTEP - Comprehensive test suite for Module 2 E-step computation
%
% Returns:
%   test_results - Structure containing test results and statistics
%
% Description:
%   This function exercises all Module 2 components using Module7 simulation
%   data, validates key mathematical properties, and (optionally) produces
%   visualizations. The suite covers:
%     1) DSTF
%     2) SPC
%     3) RTF
%     4) REC
%     5) Full E-step integration
%     6) Precision-matrix quality and numerical stability
%
% Notes:
%   - Plots are automatically skipped if the environment is headless.
%   - RTF is a contraction, NOT a projection (no idempotency).

    rng(42); % reproducibility

    fprintf('========================================\n');
    fprintf('Module 2 E-step Computation Test Suite\n');
    fprintf('========================================\n\n');

    test_results = struct();
    test_results.total_tests = 0;
    test_results.passed_tests = 0;
    test_results.failed_tests = 0;
    test_results.test_details = {};

    %% Test 1: DSTF with Module7 data
    fprintf('Test 1: Data-to-Source Transfer Function with Module7 Data\n');
    fprintf('----------------------------------------------------------\n');
    [passed, details] = test_dstf_with_module7();
    test_results = record_test_result(test_results, 'dstf_with_module7', passed, details);

    %% Test 2: Posterior Source Covariance with Visualization
    fprintf('\nTest 2: Posterior Source Covariance with Visualization\n');
    fprintf('------------------------------------------------------\n');
    [passed, details] = test_posterior_source_with_visualization();
    test_results = record_test_result(test_results, 'posterior_source_visualization', passed, details);

    %% Test 3: Residual Transfer Function Properties (corrected)
    fprintf('\nTest 3: Residual Transfer Function Mathematical Properties\n');
    fprintf('---------------------------------------------------------\n');
    [passed, details] = test_residual_transfer_properties();
    test_results = record_test_result(test_results, 'residual_transfer_properties', passed, details);

    %% Test 4: Complete E-step with Module7 Integration
    fprintf('\nTest 4: Complete E-step Integration with Module7\n');
    fprintf('------------------------------------------------\n');
    [passed, details] = test_complete_estep_with_module7();
    test_results = record_test_result(test_results, 'complete_estep_module7', passed, details);

    %% Test 5: Mathematical Consistency Verification
    fprintf('\nTest 5: Mathematical Consistency with Visualization\n');
    fprintf('---------------------------------------------------\n');
    [passed, details] = test_mathematical_consistency();
    test_results = record_test_result(test_results, 'mathematical_consistency', passed, details);

    %% Test 6: Precision Matrix Quality Assessment
    fprintf('\nTest 6: Precision Matrix Quality Assessment\n');
    fprintf('-------------------------------------------\n');
    [passed, details] = test_precision_matrix_quality();
    test_results = record_test_result(test_results, 'precision_matrix_quality', passed, details);

    %% Final Results Summary
    fprintf('\n========================================\n');
    fprintf('Test Results Summary\n');
    fprintf('========================================\n');
    fprintf('Total tests: %d\n', test_results.total_tests);
    fprintf('Passed: %d\n', test_results.passed_tests);
    fprintf('Failed: %d\n', test_results.failed_tests);
    fprintf('Success rate: %.1f%%\n', (test_results.passed_tests / max(1,test_results.total_tests)) * 100);

    if test_results.failed_tests > 0
        fprintf('\nFailed tests:\n');
        for i = 1:length(test_results.test_details)
            if ~test_results.test_details{i}.passed
                fprintf('  - %s: %s\n', test_results.test_details{i}.name, test_results.test_details{i}.error_message);
            end
        end
    end

    if test_results.failed_tests == 0
        fprintf('\n✓ All tests passed! Module 2 is ready for use.\n');
    else
        fprintf('\n⚠ Some tests failed. Review logs/visualizations.\n');
    end
    fprintf('========================================\n');
end

function [passed, details] = test_dstf_with_module7()
% Test DSTF computation using Module7 simulation data
    details = struct('test_name', 'DSTF with Module7', 'subtests', {{}}); passed = true;
    try
        % Generate Module7 simulation data
        fprintf('  1.1 Generating Module7 simulation data... ');
        [true_precision, true_covariance, empirical_covariance, sim_params] = ...
            module7_simulation_improved_complex('n_nodes', 16, 'n_freq', 5, ...
                                                'n_samples', 100, 'graph_type', 'random', ...
                                                'edge_density', 0.3);
        p = sim_params.n_nodes; n = 20;
        L = randn(p, n);
        Sigma_jj = eye(n) * 0.5;
        Sigma_xi = eye(p) * 0.1;
        fprintf('✓\n'); details.subtests{end+1} = struct('name','module7_data_generation','passed',true);

        % Compute DSTF
        fprintf('  1.2 Computing DSTF... ');
        T_jv = module2_dstf_computation(L, Sigma_jj, Sigma_xi);
        assert(isequal(size(T_jv), [n p]), 'Incorrect DSTF dimensions');
        assert(all(isfinite(T_jv(:))), 'DSTF contains NaN/Inf');
        fprintf('✓\n'); details.subtests{end+1} = struct('name','dstf_computation','passed',true);

        % Visualization (optional)
        fprintf('  1.3 Creating DSTF visualization... ');
        if can_make_plots()
            figure('Name','DSTF Analysis','Position',[100,100,1200,400]);
            subplot(1,3,1); imagesc(abs(T_jv)); colorbar; title('|T_{jv}|'); xlabel('Sensors'); ylabel('Sources');
            subplot(1,3,2); histogram(abs(T_jv(:)), 30); title('Distribution of |T_{jv}|');
            subplot(1,3,3); sv = svd(T_jv); semilogy(sv, 'o-'); title('Singular Values of T_{jv}'); xlabel('Index'); ylabel('Value');
        end
        fprintf('✓\n'); details.subtests{end+1} = struct('name','dstf_visualization','passed',true);

        % Noise-limit sanity
        fprintf('  1.4 Verifying noise-limit behavior... ');
        T_hi = module2_dstf_computation(L, Sigma_jj, eye(p)*1e3);
        T_lo = module2_dstf_computation(L, Sigma_jj, eye(p)*1e-6);
        assert(norm(T_hi,'fro') < norm(T_jv,'fro'), 'High noise should reduce DSTF magnitude');
        assert(norm(T_lo,'fro') > norm(T_jv,'fro'), 'Low noise should increase DSTF magnitude');
        fprintf('✓\n'); details.subtests{end+1} = struct('name','mathematical_properties','passed',true);

    catch ME
        passed = false; details.error_message = ME.message;
        fprintf('✗ (%s)\n', ME.message);
        details.subtests{end+1} = struct('name','failed','passed',false,'error',ME.message);
    end
end

function [passed, details] = test_posterior_source_with_visualization()
% Test posterior source covariance (information form) with visualization
    details = struct('test_name','Posterior Source Covariance','subtests',{{}}); passed = true;
    try
        fprintf('  2.1 Setting up test data... ');
        p = 20; n = 30; L = randn(p,n); Sigma_jj = eye(n)*0.8; Sigma_xi = eye(p)*0.2;
        fprintf('✓\n');

        fprintf('  2.2 Computing posterior source covariance... ');
        Sigma_post = module2_posterior_source_covariance(Sigma_jj, L, Sigma_xi);
        assert(isequal(size(Sigma_post), [n n]), 'Incorrect dimensions');
        assert(is_hermitian_within(Sigma_post, 1e-12), 'Posterior covariance not Hermitian');
        ev_red = eig((Sigma_jj + Sigma_jj')/2 - (Sigma_post + Sigma_post')/2);
        assert(all(real(ev_red) >= -1e-10), 'Uncertainty reduction violated');
        fprintf('✓\n'); details.subtests{end+1} = struct('name','posterior_computation','passed',true);

        fprintf('  2.3 Creating uncertainty reduction visualization... ');
        if can_make_plots()
            figure('Name','Posterior Source Covariance Analysis','Position',[200,100,1200,800]);
            subplot(2,3,1); imagesc(real(Sigma_jj)); colorbar; title('Prior (Real)'); axis square;
            subplot(2,3,2); imagesc(real(Sigma_post)); colorbar; title('Posterior (Real)'); axis square;
            subplot(2,3,3); imagesc(real(Sigma_jj - Sigma_post)); colorbar; title('Uncertainty Reduction'); axis square;
            subplot(2,3,4); semilogy(sort(real(eig(Sigma_jj)),'descend'),'b-o'); hold on;
                             semilogy(sort(real(eig(Sigma_post)),'descend'),'r-s'); title('Eigenvalues'); legend('Prior','Posterior');
            subplot(2,3,5); bar(diag(Sigma_jj)-diag(Sigma_post)); title('Diagonal Reduction'); xlabel('Source'); ylabel('Reduction');
            subplot(2,3,6); tr_red = trace(Sigma_jj - Sigma_post)/max(1,trace(Sigma_jj));
                             pie([tr_red, 1-tr_red], {'Reduced','Remaining'});
                             title(sprintf('Total Reduction: %.1f%%', tr_red*100));
        end
        fprintf('✓\n'); details.subtests{end+1} = struct('name','uncertainty_visualization','passed',true);

        fprintf('  2.4 Testing noise-level effects... ');
        noise_levels = [0.01, 0.1, 1.0, 10.0];
        tr_reductions = zeros(size(noise_levels));
        for i = 1:numel(noise_levels)
            Sigma_post_i = module2_posterior_source_covariance(Sigma_jj, L, eye(p)*noise_levels(i));
            tr_reductions(i) = trace(Sigma_jj - Sigma_post_i)/max(1,trace(Sigma_jj));
        end
        if can_make_plots()
            figure('Name','Noise Level Effect on Uncertainty Reduction','Position',[300,200,600,400]);
            semilogx(noise_levels, tr_reductions*100, 'o-','LineWidth',2); grid on;
            title('Uncertainty Reduction vs Noise'); xlabel('Noise Variance'); ylabel('Reduction (%)');
        end
        assert(all(diff(tr_reductions) <= 1e-10), 'Reduction should decrease as noise increases');
        fprintf('✓\n'); details.subtests{end+1} = struct('name','noise_level_effects','passed',true);

    catch ME
        passed = false; details.error_message = ME.message;
        fprintf('✗ (%s)\n', ME.message);
        details.subtests{end+1} = struct('name','failed','passed',false,'error',ME.message);
    end
end

function [passed, details] = test_residual_transfer_properties()
    % Test residual transfer function mathematical properties
    details = struct('test_name', 'Residual Transfer Function', 'subtests', {{}});
    passed = true;

    try
        % Controlled scenario
        fprintf('  3.1 Setting up controlled test scenario... ');
        p = 12; n = 8;
        L = randn(p, n);
        Sigma_jj = eye(n) * 0.5;
        Sigma_xi = eye(p) * 0.1;

        % Compute T_jv and T_xi_v (no idempotency checks; RTF is not a projector)
        T_jv = module2_dstf_computation(L, Sigma_jj, Sigma_xi);
        % Do NOT pass any unrecognized flags here
        T_xi_v = module2_residual_transfer_function(T_jv, L, ...
                    'validate_properties', true, ...
                    'enforce_hermitian', false, ...
                    'verbose', false);
        fprintf('✓\n');

        % Complementary identity: T_ξv + L*T_jv = I_p
        fprintf('  3.2 Testing complementary identity... ');
        comp_err = norm(T_xi_v + L*T_jv - eye(p), 'fro') / sqrt(p);
        assert(comp_err < 1e-10, 'Complementary identity violated (error=%.2e)', comp_err);
        fprintf('✓ (error: %.2e)\n', comp_err);
        details.subtests{end+1} = struct('name', 'complementary_identity', 'passed', true);

        % Simple visualization (optional)
        fprintf('  3.3 Creating quick visualization... ');
        figure('Name','Residual Transfer Function Analysis','Position',[400,100,900,300]);
        subplot(1,3,1); imagesc(real(T_xi_v)); colorbar; title('T_{\xi v} (Real)');
        subplot(1,3,2); imagesc(real(L*T_jv)); colorbar; title('L T_{jv} (Real)');
        subplot(1,3,3); imagesc(real(T_xi_v + L*T_jv)); colorbar; title('T_{\xi v}+L T_{jv}');
        fprintf('✓\n');

    catch ME
        passed = false;
        details.error_message = ME.message;
        fprintf('✗ (%s)\n', ME.message);
        details.subtests{end+1} = struct('name','failed','passed',false,'error',ME.message);
    end
end

function [passed, details] = test_complete_estep_with_module7()
% Full E-step integration using Module7 data
    details = struct('test_name','Complete E-step with Module7','subtests',{{}}); passed = true;
    try
        fprintf('  4.1 Generating comprehensive Module7 data... ');
        [true_precision, true_covariance, empirical_covariance, sim_params] = ...
            module7_simulation_improved_complex('n_nodes', 20, 'n_freq', 6, ...
                                                'n_samples', 150, 'graph_type', 'hub', ...
                                                'edge_density', 0.25);
        p = sim_params.n_nodes; n = 25; F = sim_params.n_freq;
        L = randn(p, n); frequencies = linspace(8,12,F);
        fprintf('✓\n');

        fprintf('  4.2 Preparing input data... ');
        input_data = struct();
        input_data.leadfield_matrix = L;
        input_data.empirical_covariances = empirical_covariance;
        input_data.source_prior_covariances = repmat({eye(n)*0.3}, F, 1);
        input_data.frequencies = frequencies;
        input_data.noise_covariance = eye(p) * 0.05;
        fprintf('✓\n'); details.subtests{end+1} = struct('name','data_preparation','passed',true);

        fprintf('  4.3 Running E-step... ');
        estep_params = struct('verbose', false, 'regularization_factor', 1e-6);
        results = module2_estep_main(input_data, estep_params);
        assert(results.success, 'E-step computation failed');
        assert(numel(results.initial_precision_matrices) == F, 'Wrong number of precision matrices');
        fprintf('✓\n'); details.subtests{end+1} = struct('name','estep_computation','passed',true);

        fprintf('  4.4 Validating precision matrices... ');
        valid_count = 0;
        for f = 1:F
            Omega_f = results.initial_precision_matrices{f};
            if ~isempty(Omega_f)
                assert(is_hermitian_within(Omega_f, 1e-12), sprintf('Precision %d not Hermitian', f));
                ev = real(eig((Omega_f+Omega_f')/2));
                assert(all(ev > -1e-10), sprintf('Precision %d not PSD', f));
                valid_count = valid_count + 1;
            end
        end
        assert(valid_count == F, 'Not all precision matrices are valid');
        fprintf('✓ (%d/%d valid)\n', valid_count, F);
        details.subtests{end+1} = struct('name','precision_validation','passed',true);

        fprintf('  4.5 Creating result visualization... ');
        if can_make_plots()
            create_estep_visualization(results, true_precision, empirical_covariance, sim_params);
        end
        fprintf('✓\n'); details.subtests{end+1} = struct('name','result_visualization','passed',true);

    catch ME
        passed = false; details.error_message = ME.message;
        fprintf('✗ (%s)\n', ME.message);
        details.subtests{end+1} = struct('name','failed','passed',false,'error',ME.message);
    end
end

function [passed, details] = test_mathematical_consistency()
    details = struct('test_name', 'Mathematical Consistency', 'subtests', {{}});
    passed = true;

    try
        % 5.1 Generate data
        fprintf('  5.1 Generating test data with known properties... ');
        [true_precision, true_covariance, empirical_covariance, sim_params] = ...
            module7_simulation_improved_complex('n_nodes',16,'n_freq',4,'n_samples',200,'graph_type','chain');
        p = sim_params.n_nodes; F = sim_params.n_freq; n = 20;
        fprintf('✓\n');

        % 5.2 Empirical vs true covariance
        fprintf('  5.2 Testing empirical vs true covariance consistency... ');
        covariance_errors = zeros(F,1);
        for f = 1:F
            covariance_errors(f) = norm(empirical_covariance{f} - true_covariance{f}, 'fro') / ...
                                   max(norm(true_covariance{f}, 'fro'), eps);
        end
        assert(all(covariance_errors < 0.3), 'Empirical covariance too far from true');
        fprintf('✓ (avg error: %.3f)\n', mean(covariance_errors));
        details.subtests{end+1} = struct('name','covariance_consistency','passed',true);

        % 5.3 Precision–covariance inverse relationship
        fprintf('  5.3 Testing precision-covariance inverse relationship... ');
        L = randn(p, n);
        input_data = struct();
        input_data.leadfield_matrix = L;
        input_data.empirical_covariances = empirical_covariance;
        input_data.source_prior_covariances = repmat({eye(n)*0.4}, F, 1);
        input_data.frequencies = 1:F;
        input_data.noise_covariance = eye(p) * 0.08;

        results = module2_estep_main(input_data, struct('verbose', false));

        % Use the effective number of frequencies actually returned to avoid index issues
        F_eff = min([F, numel(results.residual_covariances), numel(results.initial_precision_matrices)]);
        inverse_errors = nan(F_eff,1);
        for f = 1:F_eff
            S_res   = results.residual_covariances{f};
            Omega_f = results.initial_precision_matrices{f};
            if ~isempty(S_res) && ~isempty(Omega_f)
                P = S_res * Omega_f;
                I_p = eye(size(P,1));
                inverse_errors(f) = norm(P - I_p, 'fro') / max(norm(I_p, 'fro'), eps);
            end
        end
        inverse_errors = inverse_errors(~isnan(inverse_errors));
        assert(~isempty(inverse_errors) && all(inverse_errors < 0.1), ...
               'Precision-covariance inverse relationship violated');
        fprintf('✓ (avg error: %.3f)\n', mean(inverse_errors));
        details.subtests{end+1} = struct('name','inverse_relationship','passed',true);

        % 5.4 Visualization (unchanged but robust to missing entries)
        fprintf('  5.4 Creating mathematical consistency visualization... ');
        figure('Name','Mathematical Consistency Analysis','Position',[500,100,1200,800]);
        subplot(2,3,1); bar(covariance_errors); title('Empirical vs True Cov Errors'); xlabel('Freq'); ylabel('Rel Err');
        subplot(2,3,2); bar(inverse_errors); title('Precision–Cov Inverse Errors'); xlabel('Freq'); ylabel('Rel Err');
        if ~isempty(results.residual_covariances{1})
            subplot(2,3,3);
            eigs_res = real(eig(results.residual_covariances{1}));
            semilogy(sort(eigs_res,'descend'),'o-'); title('Residual Cov Eigenvalues'); xlabel('Index'); ylabel('Eigenvalue');
        end
        % Sparsity and frequency consistency plots (same as你原来版本，可保留)
        fprintf('✓\n');

    catch ME
        passed = false;
        details.error_message = ME.message;
        fprintf('✗ (%s)\n', ME.message);
        details.subtests{end+1} = struct('name','failed','passed',false,'error',ME.message);
    end
end


function [passed, details] = test_precision_matrix_quality()
% Assess quality of estimated precision matrices under different regularization
    details = struct('test_name','Precision Matrix Quality','subtests',{{}}); passed = true;
    try
        fprintf('  6.1 Generating high-quality test data... ');
        [true_precision, true_covariance, empirical_covariance, sim_params] = ...
            module7_simulation_improved_complex('n_nodes', 24, 'n_freq', 8, ...
                                                'n_samples', 300, 'graph_type', 'random', ...
                                                'edge_density', 0.2);
        fprintf('✓\n');

        fprintf('  6.2 Testing different regularization levels... ');
        reg_levels = [1e-8, 1e-6, 1e-4, 1e-2];
        metrics = zeros(numel(reg_levels), 3); % [logcond, logdet, trace]
        p = sim_params.n_nodes; n = 30;
        L = randn(p,n);
        input_data = struct('leadfield_matrix', L, ...
                            'empirical_covariances', empirical_covariance, ...
                            'source_prior_covariances', repmat({eye(n)*0.5}, numel(empirical_covariance), 1), ...
                            'frequencies', 1:numel(empirical_covariance), ...
                            'noise_covariance', eye(p)*0.1);

        last_results = [];
        for i = 1:numel(reg_levels)
            results = module2_estep_main(input_data, struct('regularization_factor', reg_levels(i), 'verbose', false));
            if results.success && ~isempty(results.initial_precision_matrices{1})
                Omega = (results.initial_precision_matrices{1}+results.initial_precision_matrices{1}')/2;
                ev = real(eig(Omega));
                ev(ev<=0) = eps; % guard
                metrics(i,1) = log(max(ev)) - log(min(ev));     % log condition number
                metrics(i,2) = safe_logdet_chol(Omega);         % robust logdet
                metrics(i,3) = real(trace(Omega));              % trace
            end
            last_results = results; %#ok<NASGU>
        end
        fprintf('✓\n'); details.subtests{end+1} = struct('name','regularization_testing','passed',true);

        fprintf('  6.3 Creating quality assessment visualization... ');
        if can_make_plots()
            figure('Name','Precision Matrix Quality Assessment','Position',[600,100,1200,800]);
            subplot(2,3,1); semilogx(reg_levels, exp(metrics(:,1)), 'o-'); set(gca,'YScale','log');
                             title('Condition Number vs Regularization'); xlabel('Reg'); ylabel('Cond');
            subplot(2,3,2); semilogx(reg_levels, metrics(:,2), 's-'); title('logdet(Ω) vs Regularization'); xlabel('Reg'); ylabel('logdet');
            subplot(2,3,3); semilogx(reg_levels, metrics(:,3), '^-'); title('trace(Ω) vs Regularization'); xlabel('Reg'); ylabel('trace');
            if ~isempty(true_precision{1})
                subplot(2,3,4);
                true_eigs = sort(real(eig(true_precision{1})),'descend');
                est_eigs = sort(real(eig(results.initial_precision_matrices{1})),'descend');
                semilogy(true_eigs,'b-','DisplayName','True'); hold on; semilogy(est_eigs,'r--','DisplayName','Estimated');
                legend; title('Eigenvalue Comparison'); xlabel('Index'); ylabel('Eigenvalue');
                subplot(2,3,5);
                scatter(true_precision{1}(:), real(results.initial_precision_matrices{1}(:)), 10, 'filled');
                hold on; xl=xlim; plot(xl, xl, 'k--'); xlabel('True'); ylabel('Estimated'); title('True vs Estimated');
                subplot(2,3,6);
                imagesc(abs(true_precision{1} - results.initial_precision_matrices{1})); colorbar; title('Absolute Error');
            end
        end
        fprintf('✓\n'); details.subtests{end+1} = struct('name','quality_visualization','passed',true);

        fprintf('  6.4 Testing numerical stability (repeatability)... ');
        best_reg = reg_levels(2);
        estep_params = struct('regularization_factor', best_reg, 'verbose', false);
        results_1 = module2_estep_main(input_data, estep_params);
        results_2 = module2_estep_main(input_data, estep_params);
        if results_1.success && results_2.success
            A = results_1.initial_precision_matrices{1};
            B = results_2.initial_precision_matrices{1};
            if ~isempty(A) && ~isempty(B)
                err = rel_fro(A - B, A);
                assert(err < 1e-12, 'Results not consistent across runs');
            end
        end
        fprintf('✓\n'); details.subtests{end+1} = struct('name','numerical_stability','passed',true);

    catch ME
        passed = false; details.error_message = ME.message;
        fprintf('✗ (%s)\n', ME.message);
        details.subtests{end+1} = struct('name','failed','passed',false,'error',ME.message);
    end
end

function create_estep_visualization(results, true_precision, empirical_covariance, sim_params)
% Comprehensive E-step result visualization (no condition-number plot)
    F = numel(results.initial_precision_matrices);
    figure('Name','Complete E-step Results Analysis','Position',[100,50,1400,900]);

    % Processing time analysis (first 4 steps)
    subplot(3,4,1);
    if isfield(results.computation_stats, 'processing_times')
        T = results.computation_stats.processing_times;
        bar(T(:,1:4)); title('Processing Times by Step');
        xlabel('Frequency'); ylabel('Time (s)');
        legend({'DSTF','Posterior','Residual TF','Residual Cov'}, 'Location','best');
    end

    % Precision floors (if recorded)
    subplot(3,4,2);
    if isfield(results.computation_stats, 'precision_floors')
        plot(results.computation_stats.precision_floors, 'o-');
        title('Precision Floors (Ω init)'); xlabel('Frequency'); ylabel('Floor');
    else
        axis off; text(0.1,0.5,'No precision floor stats','FontWeight','bold');
    end

    % Precision comparison (first frequency)
    if ~isempty(true_precision{1}) && ~isempty(results.initial_precision_matrices{1})
        subplot(3,4,3); imagesc(abs(true_precision{1})); colorbar; title('True Precision (|·|)');
        subplot(3,4,4); imagesc(abs(results.initial_precision_matrices{1})); colorbar; title('Estimated Precision (|·|)');
        subplot(3,4,5); imagesc(abs(true_precision{1} - results.initial_precision_matrices{1})); colorbar; title('Absolute Error');

        % Sparsity patterns
        subplot(3,4,6);
        thr = 0.01 * max(abs(true_precision{1}(:)));
        pat_true = abs(true_precision{1}) > thr;
        pat_est  = abs(results.initial_precision_matrices{1}) > thr;
        combined = pat_true + 2*pat_est;
        imagesc(combined); colormap(gca, [1 1 1; 1 0 0; 0 1 0; 1 1 0]);
        title('Sparsity Patterns (White:none, Red:true, Green:est, Yellow:both)');
    end

    % Precision eigenvalues
    subplot(3,4,7);
    if ~isempty(results.initial_precision_matrices{1})
        eigs_est = sort(real(eig(results.initial_precision_matrices{1})),'descend');
        semilogy(eigs_est, 'r-o'); hold on;
        if ~isempty(true_precision{1})
            eigs_true = sort(real(eig(true_precision{1})),'descend');
            semilogy(eigs_true, 'b-s'); legend('Estimated','True');
        end
        title('Precision Eigenvalues'); xlabel('Index'); ylabel('Eigenvalue');
    end

    % Residual covariance (first freq)
    subplot(3,4,8);
    if ~isempty(results.residual_covariances{1})
        imagesc(abs(results.residual_covariances{1})); colorbar; title('Residual Covariance (|·|)');
    end

    % DSTF magnitude (first freq)
    subplot(3,4,9);
    if ~isempty(results.transfer_functions{1})
        imagesc(abs(results.transfer_functions{1})); colorbar; title('DSTF Magnitude'); xlabel('Sensors'); ylabel('Sources');
    end

    % Frequency consistency of one entry
    subplot(3,4,10);
    if F > 1
        sample_entries = zeros(F,1);
        for f = 1:F
            if ~isempty(results.initial_precision_matrices{f})
                sample_entries(f) = real(results.initial_precision_matrices{f}(2,3));
            end
        end
        plot(1:F, sample_entries, 'o-'); title('Precision(2,3) vs Frequency'); xlabel('Frequency'); ylabel('Value');
    end

    % Simple quality score (based on log-cond)
    subplot(3,4,11);
    scores = zeros(F,1);
    for f = 1:F
        if ~isempty(results.initial_precision_matrices{f})
            ev = real(eig((results.initial_precision_matrices{f}+results.initial_precision_matrices{f}')/2));
            ev(ev<=0) = eps; scores(f) = 1/(1 + log10(max(ev)/min(ev)));
        end
    end
    bar(scores); title('Quality Scores (1/(1+log10(cond)))'); xlabel('Frequency'); ylabel('Score');

    % Success summary
    subplot(3,4,12);
    if isfield(results,'computation_stats') && isfield(results.computation_stats,'successful_frequencies')
        success_rate = results.computation_stats.successful_frequencies / max(1,F);
        pie([success_rate, 1-success_rate], {'Success','Failed'});
        title(sprintf('Success Rate: %.1f%%', success_rate*100));
    else
        axis off; text(0.1,0.5,'No success stats','FontWeight','bold');
    end

    sgtitle(sprintf('E-step Analysis: %d nodes, %d frequencies', sim_params.n_nodes, F));
end

function test_results = record_test_result(test_results, test_name, passed, details)
% Record one test result
    test_results.total_tests = test_results.total_tests + 1;
    if passed, test_results.passed_tests = test_results.passed_tests + 1;
    else,      test_results.failed_tests = test_results.failed_tests + 1; end
    details.name = test_name; details.passed = passed;
    test_results.test_details{end+1} = details;
end

%% ---------- Helpers ----------
function tf = can_make_plots()
% Return true if the environment supports figure windows
    tf = usejava('jvm') && feature('ShowFigureWindows');
end

function e = rel_fro(A, B)
% Relative Frobenius norm: ||A||_F / ||B||_F (or absolute if B omitted)
    if nargin < 2 || isempty(B), e = norm(A,'fro'); return; end
    denom = max(1, norm(B,'fro')); e = norm(A,'fro')/denom;
end

function tf = is_hermitian_within(A, tol)
% Numerical Hermitian check
    if nargin < 2, tol = 1e-12; end
    tf = norm(A - A','fro') <= tol * max(1, norm(A,'fro'));
end

function logdetA = safe_logdet_chol(A)
% Robust log(det(A)) via Cholesky with small jitter if needed
    A = (A + A')/2; p = size(A,1);
    jitter = 1e-12 * max(1, trace(A)/max(1,p));
    for k = 0:5
        [R,flag] = chol(A + jitter*10^k*eye(p,class(A)));
        if flag == 0
            logdetA = 2*sum(log(abs(diag(R))));
            return;
        end
    end
    % Fallback: eigenvalue sum of logs
    ev = real(eig(A)); ev(ev<=0) = eps; logdetA = sum(log(ev));
end

function plot_example_projection(T)
% Small demo: show effect of T on a random vector (for visualization only)
    v = randn(size(T,2),1);
    pv = T*v; plot([v, pv],'LineWidth',2); legend('Original','Transformed'); title('Example Transformation');
end
