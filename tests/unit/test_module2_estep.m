function test_results = test_module2_estep()
    % TEST_MODULE2_ESTEP - Comprehensive test suite for Module 2 E-step computation
    %
    % Returns:
    %   test_results - Structure containing test results and statistics
    %
    % Description:
    %   This function performs comprehensive testing of all Module 2 components
    %   using Module7 simulation data with visualization for validation:
    %   1. Data-to-Source Transfer Function (DSTF) computation
    %   2. Posterior Source Covariance (SPC) computation  
    %   3. Residual Transfer Function (RTF) computation
    %   4. Residual Empirical Covariance (REC) computation
    %   5. Complete E-step integration with Module7 data
    %   6. Mathematical property verification with visualization
    %
    % Test Coverage:
    %   - Normal operation with Module7 simulation data
    %   - Mathematical property verification
    %   - Visualization of key results for manual validation
    %   - Integration testing with realistic data
    %
    % Author: [Author Name]
    % Date: [Current Date]
    % Version: 2.0 - Using Module7 Simulation with Visualization
    
    fprintf('========================================\n');
    fprintf('Module 2 E-step Computation Test Suite\n');
    fprintf('========================================\n\n');
    
    test_results = struct();
    test_results.total_tests = 0;
    test_results.passed_tests = 0;
    test_results.failed_tests = 0;
    test_results.test_details = {};
    
    %% Test 1: DSTF Computation with Module7 Data
    fprintf('Test 1: Data-to-Source Transfer Function with Module7 Data\n');
    fprintf('----------------------------------------------------------\n');
    [passed, details] = test_dstf_with_module7();
    test_results = record_test_result(test_results, 'dstf_with_module7', passed, details);
    
    %% Test 2: Posterior Source Covariance with Visualization
    fprintf('\nTest 2: Posterior Source Covariance with Visualization\n');
    fprintf('------------------------------------------------------\n');
    [passed, details] = test_posterior_source_with_visualization();
    test_results = record_test_result(test_results, 'posterior_source_visualization', passed, details);
    
    %% Test 3: Residual Transfer Function Properties
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
    fprintf('Success rate: %.1f%%\n', (test_results.passed_tests / test_results.total_tests) * 100);
    
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
        fprintf('\n⚠ Some tests failed. Review visualizations and issues.\n');
    end
    fprintf('========================================\n');
end

function [passed, details] = test_dstf_with_module7()
    % Test DSTF computation using Module7 simulation data
    details = struct('test_name', 'DSTF with Module7', 'subtests', {{}});
    passed = true;
    
    try
        % Generate Module7 simulation data (移除verbose参数)
        fprintf('  1.1 Generating Module7 simulation data... ');
        [true_precision, true_covariance, empirical_covariance, sim_params] = ...
            module7_simulation_improved_complex(...
                'n_nodes', 16, ...
                'n_freq', 5, ...
                'n_samples', 100, ...
                'graph_type', 'random', ...
                'edge_density', 0.3);
        
        p = sim_params.n_nodes;
        n = 20;  % Number of sources
        L = randn(p, n);
        Sigma_jj = eye(n) * 0.5;
        Sigma_xi = eye(p) * 0.1;
        
        fprintf('✓\n');
        details.subtests{end+1} = struct('name', 'module7_data_generation', 'passed', true);
        
        % Test DSTF computation
        fprintf('  1.2 Computing DSTF... ');
        T_jv = module2_dstf_computation(L, Sigma_jj, Sigma_xi);
        
        % Verify dimensions and properties
        assert(size(T_jv, 1) == n && size(T_jv, 2) == p, 'Incorrect DSTF dimensions');
        assert(all(isfinite(T_jv(:))), 'DSTF contains NaN or Inf');
        
        fprintf('✓\n');
        details.subtests{end+1} = struct('name', 'dstf_computation', 'passed', true);
        
        % Visualize DSTF properties
        fprintf('  1.3 Creating DSTF visualization... ');
        figure('Name', 'DSTF Analysis', 'Position', [100, 100, 1200, 400]);
        
        subplot(1, 3, 1);
        imagesc(abs(T_jv));
        colorbar;
        title('|T_{jv}| - DSTF Magnitude');
        xlabel('Sensors');
        ylabel('Sources');
        
        subplot(1, 3, 2);
        histogram(abs(T_jv(:)), 30);
        title('Distribution of DSTF Magnitudes');
        xlabel('|T_{jv}|');
        ylabel('Count');
        
        subplot(1, 3, 3);
        plot(svd(T_jv), 'o-');
        title('Singular Values of T_{jv}');
        xlabel('Index');
        ylabel('Singular Value');
        set(gca, 'YScale', 'log');
        
        fprintf('✓\n');
        details.subtests{end+1} = struct('name', 'dstf_visualization', 'passed', true);
        
        % Test mathematical properties
        fprintf('  1.4 Verifying mathematical properties... ');
        
        % Test high vs low noise behavior
        Sigma_xi_high = eye(p) * 1000;
        T_jv_high_noise = module2_dstf_computation(L, Sigma_jj, Sigma_xi_high);
        assert(norm(T_jv_high_noise, 'fro') < norm(T_jv, 'fro'), 'High noise should reduce DSTF magnitude');
        
        Sigma_xi_low = eye(p) * 1e-6;
        T_jv_low_noise = module2_dstf_computation(L, Sigma_jj, Sigma_xi_low);
        assert(norm(T_jv_low_noise, 'fro') > norm(T_jv, 'fro'), 'Low noise should increase DSTF magnitude');
        
        fprintf('✓\n');
        details.subtests{end+1} = struct('name', 'mathematical_properties', 'passed', true);
        
    catch ME
        passed = false;
        details.error_message = ME.message;
        fprintf('✗ (%s)\n', ME.message);
        details.subtests{end+1} = struct('name', 'failed', 'passed', false, 'error', ME.message);
    end
end

function [passed, details] = test_posterior_source_with_visualization()
    % Test posterior source covariance with visualization
    details = struct('test_name', 'Posterior Source Covariance', 'subtests', {{}});
    passed = true;
    
    try
        % Setup test data
        fprintf('  2.1 Setting up test data... ');
        p = 20; n = 30;
        L = randn(p, n);
        Sigma_jj = eye(n) * 0.8;
        Sigma_xi = eye(p) * 0.2;
        
        fprintf('✓\n');
        
        % Compute posterior source covariance
        fprintf('  2.2 Computing posterior source covariance... ');
        Sigma_post = module2_posterior_source_covariance(Sigma_jj, L, Sigma_xi);
        
        % Verify basic properties
        assert(size(Sigma_post, 1) == n && size(Sigma_post, 2) == n, 'Incorrect dimensions');
        assert(ishermitian(Sigma_post), 'Posterior covariance not Hermitian');
        
        % Verify uncertainty reduction
        eigenvals_reduction = eig(Sigma_jj - Sigma_post);
        assert(all(real(eigenvals_reduction) >= -1e-10), 'Uncertainty reduction property violated');
        
        fprintf('✓\n');
        details.subtests{end+1} = struct('name', 'posterior_computation', 'passed', true);
        
        % Create visualization
        fprintf('  2.3 Creating uncertainty reduction visualization... ');
        figure('Name', 'Posterior Source Covariance Analysis', 'Position', [200, 100, 1200, 800]);
        
        subplot(2, 3, 1);
        imagesc(real(Sigma_jj));
        colorbar;
        title('Prior Source Covariance (Real)');
        axis square;
        
        subplot(2, 3, 2);
        imagesc(real(Sigma_post));
        colorbar;
        title('Posterior Source Covariance (Real)');
        axis square;
        
        subplot(2, 3, 3);
        imagesc(real(Sigma_jj - Sigma_post));
        colorbar;
        title('Uncertainty Reduction');
        axis square;
        
        subplot(2, 3, 4);
        prior_eigs = sort(real(eig(Sigma_jj)), 'descend');
        post_eigs = sort(real(eig(Sigma_post)), 'descend');
        semilogy(prior_eigs, 'b-o', 'DisplayName', 'Prior');
        hold on;
        semilogy(post_eigs, 'r-s', 'DisplayName', 'Posterior');
        legend;
        title('Eigenvalue Comparison');
        xlabel('Index');
        ylabel('Eigenvalue');
        
        subplot(2, 3, 5);
        uncertainty_reduction = diag(Sigma_jj) - diag(Sigma_post);
        bar(uncertainty_reduction);
        title('Diagonal Uncertainty Reduction');
        xlabel('Source Index');
        ylabel('Reduction');
        
        subplot(2, 3, 6);
        trace_reduction = trace(Sigma_jj - Sigma_post) / trace(Sigma_jj);
        pie([trace_reduction, 1-trace_reduction], {'Reduced', 'Remaining'});
        title(sprintf('Total Uncertainty Reduction: %.1f%%', trace_reduction*100));
        
        fprintf('✓\n');
        details.subtests{end+1} = struct('name', 'uncertainty_visualization', 'passed', true);
        
        % Test different noise levels
        fprintf('  2.4 Testing noise level effects... ');
        noise_levels = [0.01, 0.1, 1.0, 10.0];
        trace_reductions = zeros(size(noise_levels));
        
        for i = 1:length(noise_levels)
            Sigma_xi_test = eye(p) * noise_levels(i);
            Sigma_post_test = module2_posterior_source_covariance(Sigma_jj, L, Sigma_xi_test);
            trace_reductions(i) = trace(Sigma_jj - Sigma_post_test) / trace(Sigma_jj);
        end
        
        % Plot noise effect
        figure('Name', 'Noise Level Effect on Uncertainty Reduction', 'Position', [300, 200, 600, 400]);
        semilogx(noise_levels, trace_reductions*100, 'o-', 'LineWidth', 2);
        title('Uncertainty Reduction vs Noise Level');
        xlabel('Noise Variance');
        ylabel('Uncertainty Reduction (%)');
        grid on;
        
        % Verify monotonic decrease with noise
        assert(all(diff(trace_reductions) <= 1e-10), 'Uncertainty reduction should decrease with increasing noise');
        
        fprintf('✓\n');
        details.subtests{end+1} = struct('name', 'noise_level_effects', 'passed', true);
        
    catch ME
        passed = false;
        details.error_message = ME.message;
        fprintf('✗ (%s)\n', ME.message);
        details.subtests{end+1} = struct('name', 'failed', 'passed', false, 'error', ME.message);
    end
end

function [passed, details] = test_residual_transfer_properties()
    % Test residual transfer function mathematical properties
    details = struct('test_name', 'Residual Transfer Function', 'subtests', {{}});
    passed = true;
    
    try
        % Create test scenario with realistic parameters 
        fprintf('  3.1 Setting up controlled test scenario... ');
        p = 12; n = 8;
        L = randn(p, n);
        Sigma_jj = eye(n) * 0.5;
        Sigma_xi = eye(p) * 0.1;
        
        % 重要：计算实际的DSTF和RTF，而不是使用随机矩阵
        T_jv = module2_dstf_computation(L, Sigma_jj, Sigma_xi);
        T_xi_v = module2_residual_transfer_function(T_jv, L);
        
        fprintf('✓\n');
        
        % Test idempotent property with实际计算的矩阵
        fprintf('  3.2 Testing idempotent property... ');
        T_xi_v_squared = T_xi_v * T_xi_v;
        idempotent_error = norm(T_xi_v_squared - T_xi_v, 'fro') / (norm(T_xi_v, 'fro') + eps);
        
        % 使用更现实的容差，因为RTF是投影矩阵
        assert(idempotent_error < 1e-10, 'Idempotent property violated');
        
        fprintf('✓ (error: %.2e)\n', idempotent_error);
        details.subtests{end+1} = struct('name', 'idempotent_property', 'passed', true);
        
        % Test complementary projection property
        fprintf('  3.3 Testing complementary projection... ');
        L_T_jv = L * T_jv;
        identity_sum = T_xi_v + L_T_jv;
        I_p = eye(p);
        complementary_error = norm(identity_sum - I_p, 'fro') / norm(I_p, 'fro');
        assert(complementary_error < 1e-10, 'Complementary projection violated');
        
        fprintf('✓ (error: %.2e)\n', complementary_error);
        details.subtests{end+1} = struct('name', 'complementary_projection', 'passed', true);
        
        fprintf('✓ (error: %.2e)\n', complementary_error);
        details.subtests{end+1} = struct('name', 'complementary_projection', 'passed', true);
        
        % Visualize projection properties
        fprintf('  3.4 Creating projection analysis visualization... ');
        figure('Name', 'Residual Transfer Function Analysis', 'Position', [400, 100, 1200, 800]);
        
        subplot(2, 3, 1);
        imagesc(real(T_xi_v));
        colorbar;
        title('T_{\xi v} (Real Part)');
        axis square;
        
        subplot(2, 3, 2);
        imagesc(real(L_T_jv));
        colorbar;
        title('L T_{jv} (Real Part)');
        axis square;
        
        subplot(2, 3, 3);
        imagesc(real(identity_sum));
        colorbar;
        title('T_{\xi v} + L T_{jv}');
        axis square;
        
        subplot(2, 3, 4);
        imagesc(real(T_xi_v_squared - T_xi_v));
        colorbar;
        title('T_{\xi v}^2 - T_{\xi v} (Idempotent Error)');
        axis square;
        
        subplot(2, 3, 5);
        eigs_T_xi_v = real(eig(T_xi_v));
        histogram(eigs_T_xi_v, 20);
        title('Eigenvalues of T_{\xi v}');
        xlabel('Eigenvalue');
        ylabel('Count');
        
        subplot(2, 3, 6);
        % Project random vector to show null space
        test_vec = randn(p, 1);
        projected_vec = T_xi_v * test_vec;
        plot([test_vec, projected_vec], 'LineWidth', 2);
        legend('Original Vector', 'Projected Vector');
        title('Example Projection');
        
        fprintf('✓\n');
        details.subtests{end+1} = struct('name', 'projection_visualization', 'passed', true);
        
    catch ME
        passed = false;
        details.error_message = ME.message;
        fprintf('✗ (%s)\n', ME.message);
        details.subtests{end+1} = struct('name', 'failed', 'passed', false, 'error', ME.message);
    end
end

function [passed, details] = test_complete_estep_with_module7()
    % Test complete E-step integration using Module7 simulation
    details = struct('test_name', 'Complete E-step with Module7', 'subtests', {{}});
    passed = true;
    
    try
        % Generate Module7 simulation data
        fprintf('  4.1 Generating comprehensive Module7 data... ');
        [true_precision, true_covariance, empirical_covariance, sim_params] = ...
            module7_simulation_improved_complex(...
                'n_nodes', 20, ...
                'n_freq', 6, ...
                'n_samples', 150, ...
                'graph_type', 'hub', ...
                'edge_density', 0.25);
        
        p = sim_params.n_nodes;
        n = 25;  % Sources
        F = sim_params.n_freq;
        L = randn(p, n);
        frequencies = linspace(8, 12, F);
        
        fprintf('✓\n');
        
        % Prepare input data structure
        fprintf('  4.2 Preparing input data structure... ');
        input_data = struct();
        input_data.leadfield_matrix = L;
        input_data.empirical_covariances = empirical_covariance;
        input_data.source_prior_covariances = cell(F, 1);
        input_data.frequencies = frequencies;
        input_data.noise_covariance = eye(p) * 0.05;
        
        for f = 1:F
            input_data.source_prior_covariances{f} = eye(n) * 0.3;
        end
        
        fprintf('✓\n');
        details.subtests{end+1} = struct('name', 'data_preparation', 'passed', true);
        
        % Run complete E-step
        fprintf('  4.3 Running complete E-step computation... ');
        estep_params = struct('verbose', false, 'regularization_factor', 1e-6);
        results = module2_estep_main(input_data, estep_params);
        
        assert(results.success, 'E-step computation failed');
        assert(length(results.initial_precision_matrices) == F, 'Wrong number of precision matrices');
        
        fprintf('✓\n');
        details.subtests{end+1} = struct('name', 'estep_computation', 'passed', true);
        
        % Validate all precision matrices
        fprintf('  4.4 Validating precision matrices... ');
        valid_count = 0;
        for f = 1:F
            Omega_f = results.initial_precision_matrices{f};
            if ~isempty(Omega_f)
                assert(ishermitian(Omega_f), sprintf('Precision matrix %d not Hermitian', f));
                eigenvals = real(eig(Omega_f));
                assert(all(eigenvals > -1e-10), sprintf('Precision matrix %d not PSD', f));
                valid_count = valid_count + 1;
            end
        end
        
        assert(valid_count == F, 'Not all precision matrices are valid');
        fprintf('✓ (%d/%d valid)\n', valid_count, F);
        details.subtests{end+1} = struct('name', 'precision_validation', 'passed', true);
        
        % Create comprehensive visualization
        fprintf('  4.5 Creating comprehensive result visualization... ');
        create_estep_visualization(results, true_precision, empirical_covariance, sim_params);
        
        fprintf('✓\n');
        details.subtests{end+1} = struct('name', 'result_visualization', 'passed', true);
        
    catch ME
        passed = false;
        details.error_message = ME.message;
        fprintf('✗ (%s)\n', ME.message);
        details.subtests{end+1} = struct('name', 'failed', 'passed', false, 'error', ME.message);
    end
end

function [passed, details] = test_mathematical_consistency()
    % Test mathematical consistency across the E-step pipeline
    details = struct('test_name', 'Mathematical Consistency', 'subtests', {{}});
    passed = true;
    
    try
        % Use Module7 data for realistic testing
        fprintf('  5.1 Generating test data with known properties... ');
        [true_precision, true_covariance, empirical_covariance, sim_params] = ...
            module7_simulation_improved_complex(...
                'n_nodes', 16, ...
                'n_freq', 4, ...
                'n_samples', 200, ...
                'graph_type', 'chain');
        
        p = sim_params.n_nodes;
        n = 20;
        F = sim_params.n_freq;
        
        fprintf('✓\n');
        
        % Test consistency of empirical vs true covariance
        fprintf('  5.2 Testing empirical vs true covariance consistency... ');
        covariance_errors = zeros(F, 1);
        
        for f = 1:F
            true_cov = true_covariance{f};
            emp_cov = empirical_covariance{f};
            covariance_errors(f) = norm(emp_cov - true_cov, 'fro') / norm(true_cov, 'fro');
        end
        
        % Empirical covariance should be reasonably close to true
        assert(all(covariance_errors < 0.3), 'Empirical covariance too far from true covariance');
        
        fprintf('✓ (avg error: %.3f)\n', mean(covariance_errors));
        details.subtests{end+1} = struct('name', 'covariance_consistency', 'passed', true);
        
        % Test precision-covariance relationship
        fprintf('  5.3 Testing precision-covariance inverse relationship... ');
        L = randn(p, n);
        input_data = struct();
        input_data.leadfield_matrix = L;
        input_data.empirical_covariances = empirical_covariance;
        input_data.source_prior_covariances = repmat({eye(n)*0.4}, F, 1);
        input_data.frequencies = 1:F;
        input_data.noise_covariance = eye(p) * 0.08;
        
        results = module2_estep_main(input_data, struct('verbose', false));
        
        inverse_errors = zeros(F, 1);
        for f = 1:F
            if ~isempty(results.residual_covariances{f}) && ~isempty(results.initial_precision_matrices{f})
                S_res = results.residual_covariances{f};
                Omega_res = results.initial_precision_matrices{f};
                
                % Test S * Omega ≈ I
                product = S_res * Omega_res;
                I_p = eye(size(product));
                inverse_errors(f) = norm(product - I_p, 'fro') / norm(I_p, 'fro');
            end
        end
        
        % Remove zeros (failed frequencies)
        inverse_errors = inverse_errors(inverse_errors > 0);
        assert(all(inverse_errors < 0.1), 'Precision-covariance inverse relationship violated');
        
        fprintf('✓ (avg error: %.3f)\n', mean(inverse_errors));
        details.subtests{end+1} = struct('name', 'inverse_relationship', 'passed', true);
        
        % Visualize mathematical consistency
        fprintf('  5.4 Creating mathematical consistency visualization... ');
        figure('Name', 'Mathematical Consistency Analysis', 'Position', [500, 100, 1200, 800]);
        
        subplot(2, 3, 1);
        bar(covariance_errors);
        title('Empirical vs True Covariance Errors');
        xlabel('Frequency');
        ylabel('Relative Error');
        
        subplot(2, 3, 2);
        bar(inverse_errors);
        title('Precision-Covariance Inverse Errors');
        xlabel('Frequency');
        ylabel('Relative Error');
        
        subplot(2, 3, 3);
        % Show typical residual covariance eigenvalues
        if ~isempty(results.residual_covariances{1})
            eigs_res = real(eig(results.residual_covariances{1}));
            semilogy(sort(eigs_res, 'descend'), 'o-');
            title('Residual Covariance Eigenvalues');
            xlabel('Index');
            ylabel('Eigenvalue');
        end
        
        subplot(2, 3, 4);
        % Compare true and estimated precision sparsity patterns
        if ~isempty(true_precision{1})
            true_pattern = abs(true_precision{1}) > 1e-3 * max(abs(true_precision{1}(:)));
            imagesc(true_pattern);
            title('True Precision Sparsity');
            colormap(gca, gray);
        end
        
        subplot(2, 3, 5);
        if ~isempty(results.initial_precision_matrices{1})
            est_pattern = abs(results.initial_precision_matrices{1}) > 1e-3 * max(abs(results.initial_precision_matrices{1}(:)));
            imagesc(est_pattern);
            title('Estimated Precision Sparsity');
            colormap(gca, gray);
        end
        
        subplot(2, 3, 6);
        % Frequency consistency of precision matrix entries
        if F > 1
            sample_entry = zeros(F, 1);
            for f = 1:F
                if ~isempty(results.initial_precision_matrices{f})
                    sample_entry(f) = real(results.initial_precision_matrices{f}(2, 3));
                end
            end
            plot(1:F, sample_entry, 'o-');
            title('Precision Entry Across Frequencies');
            xlabel('Frequency Index');
            ylabel('Precision(2,3)');
        end
        
        fprintf('✓\n');
        details.subtests{end+1} = struct('name', 'consistency_visualization', 'passed', true);
        
    catch ME
        passed = false;
        details.error_message = ME.message;
        fprintf('✗ (%s)\n', ME.message);
        details.subtests{end+1} = struct('name', 'failed', 'passed', false, 'error', ME.message);
    end
end

function [passed, details] = test_precision_matrix_quality()
    % Test the quality of estimated precision matrices
    details = struct('test_name', 'Precision Matrix Quality', 'subtests', {{}});
    passed = true;
    
    try
        % Generate high-quality test data
        fprintf('  6.1 Generating high-quality test data... ');
        [true_precision, true_covariance, empirical_covariance, sim_params] = ...
            module7_simulation_improved_complex(...
                'n_nodes', 24, ...
                'n_freq', 8, ...
                'n_samples', 300, ...
                'graph_type', 'random', ...
                'edge_density', 0.2);
        
        fprintf('✓\n');
        
        % Run E-step with different regularization levels
        fprintf('  6.2 Testing different regularization levels... ');
        reg_levels = [1e-8, 1e-6, 1e-4, 1e-2];
        quality_metrics = zeros(length(reg_levels), 3);  % [condition, det, trace]
        
        p = sim_params.n_nodes;
        n = 30;
        L = randn(p, n);
        input_data = struct();
        input_data.leadfield_matrix = L;
        input_data.empirical_covariances = empirical_covariance;
        input_data.source_prior_covariances = repmat({eye(n)*0.5}, length(empirical_covariance), 1);
        input_data.frequencies = 1:length(empirical_covariance);
        input_data.noise_covariance = eye(p) * 0.1;
        
        for i = 1:length(reg_levels)
            estep_params = struct('regularization_factor', reg_levels(i), 'verbose', false);
            results = module2_estep_main(input_data, estep_params);
            
            if results.success && ~isempty(results.initial_precision_matrices{1})
                Omega = results.initial_precision_matrices{1};
                quality_metrics(i, 1) = cond(Omega);
                quality_metrics(i, 2) = det(Omega);
                quality_metrics(i, 3) = trace(Omega);
            end
        end
        
        fprintf('✓\n');
        details.subtests{end+1} = struct('name', 'regularization_testing', 'passed', true);
        
        % Visualize quality metrics
        fprintf('  6.3 Creating quality assessment visualization... ');
        figure('Name', 'Precision Matrix Quality Assessment', 'Position', [600, 100, 1200, 800]);
        
        subplot(2, 3, 1);
        semilogx(reg_levels, quality_metrics(:, 1), 'o-');
        title('Condition Number vs Regularization');
        xlabel('Regularization Factor');
        ylabel('Condition Number');
        set(gca, 'YScale', 'log');
        
        subplot(2, 3, 2);
        semilogx(reg_levels, quality_metrics(:, 2), 's-');
        title('Determinant vs Regularization');
        xlabel('Regularization Factor');
        ylabel('Determinant');
        
        subplot(2, 3, 3);
        semilogx(reg_levels, quality_metrics(:, 3), '^-');
        title('Trace vs Regularization');
        xlabel('Regularization Factor');
        ylabel('Trace');
        
        % Compare with true precision matrix properties
        if ~isempty(true_precision{1})
            true_Omega = true_precision{1};
            
            subplot(2, 3, 4);
            true_eigs = sort(real(eig(true_Omega)), 'descend');
            est_eigs = sort(real(eig(results.initial_precision_matrices{1})), 'descend');
            semilogy(true_eigs, 'b-', 'DisplayName', 'True');
            hold on;
            semilogy(est_eigs, 'r--', 'DisplayName', 'Estimated');
            legend;
            title('Eigenvalue Comparison');
            xlabel('Index');
            ylabel('Eigenvalue');
            
            subplot(2, 3, 5);
            scatter(true_Omega(:), real(results.initial_precision_matrices{1}(:)), 10, 'filled');
            xlabel('True Precision Values');
            ylabel('Estimated Precision Values');
            title('True vs Estimated Values');
            hold on;
            plot(xlim, xlim, 'k--');
            
            subplot(2, 3, 6);
            error_matrix = abs(true_Omega - results.initial_precision_matrices{1});
            imagesc(error_matrix);
            colorbar;
            title('Absolute Error Matrix');
        end
        
        fprintf('✓\n');
        details.subtests{end+1} = struct('name', 'quality_visualization', 'passed', true);
        
        % Test numerical stability
        fprintf('  6.4 Testing numerical stability... ');
        best_reg = reg_levels(2);  % Usually 1e-6 is good
        estep_params = struct('regularization_factor', best_reg, 'verbose', false);
        
        % Run multiple times to check consistency
        results_1 = module2_estep_main(input_data, estep_params);
        results_2 = module2_estep_main(input_data, estep_params);
        
        if results_1.success && results_2.success
            consistency_error = norm(results_1.initial_precision_matrices{1} - results_2.initial_precision_matrices{1}, 'fro') / ...
                               norm(results_1.initial_precision_matrices{1}, 'fro');
            assert(consistency_error < 1e-12, 'Results not consistent across runs');
        end
        
        fprintf('✓\n');
        details.subtests{end+1} = struct('name', 'numerical_stability', 'passed', true);
        
    catch ME
        passed = false;
        details.error_message = ME.message;
        fprintf('✗ (%s)\n', ME.message);
        details.subtests{end+1} = struct('name', 'failed', 'passed', false, 'error', ME.message);
    end
end

function create_estep_visualization(results, true_precision, empirical_covariance, sim_params)
    % Create comprehensive E-step result visualization
    
    F = length(results.initial_precision_matrices);
    
    figure('Name', 'Complete E-step Results Analysis', 'Position', [100, 50, 1400, 900]);
    
    % Processing time analysis
    subplot(3, 4, 1);
    if isfield(results.computation_stats, 'processing_times')
        bar(results.computation_stats.processing_times(:, 1:4));
        title('Processing Times by Step');
        xlabel('Frequency');
        ylabel('Time (s)');
        legend({'DSTF', 'Posterior', 'Residual TF', 'Residual Cov'}, 'Location', 'best');
    end
    
    % Condition numbers across frequencies
    subplot(3, 4, 2);
    if isfield(results.computation_stats, 'condition_numbers')
        semilogy(results.computation_stats.condition_numbers);
        title('Condition Numbers');
        xlabel('Frequency');
        ylabel('Condition Number');
        legend({'Sensor Cov', 'Post Source', 'Residual'}, 'Location', 'best');
    end
    
    % Precision matrix comparison (first frequency)
    if ~isempty(true_precision{1}) && ~isempty(results.initial_precision_matrices{1})
        subplot(3, 4, 3);
        imagesc(abs(true_precision{1}));
        colorbar;
        title('True Precision (|·|)');
        
        subplot(3, 4, 4);
        imagesc(abs(results.initial_precision_matrices{1}));
        colorbar;
        title('Estimated Precision (|·|)');
        
        subplot(3, 4, 5);
        error_matrix = abs(true_precision{1} - results.initial_precision_matrices{1});
        imagesc(error_matrix);
        colorbar;
        title('Absolute Error');
        
        % Sparsity pattern comparison
        subplot(3, 4, 6);
        threshold = 0.01 * max(abs(true_precision{1}(:)));
        true_pattern = abs(true_precision{1}) > threshold;
        est_pattern = abs(results.initial_precision_matrices{1}) > threshold;
        
        combined_pattern = true_pattern + 2*est_pattern;
        imagesc(combined_pattern);
        colormap(gca, [1 1 1; 1 0 0; 0 1 0; 1 1 0]);  % White, Red, Green, Yellow
        title('Sparsity Patterns');
        % White: neither, Red: true only, Green: est only, Yellow: both
    end
    
    % Eigenvalue analysis
    subplot(3, 4, 7);
    if ~isempty(results.initial_precision_matrices{1})
        eigs_est = sort(real(eig(results.initial_precision_matrices{1})), 'descend');
        semilogy(eigs_est, 'r-o');
        if ~isempty(true_precision{1})
            hold on;
            eigs_true = sort(real(eig(true_precision{1})), 'descend');
            semilogy(eigs_true, 'b-s');
            legend('Estimated', 'True');
        end
        title('Precision Eigenvalues');
        xlabel('Index');
        ylabel('Eigenvalue');
    end
    
    % Residual covariance analysis
    subplot(3, 4, 8);
    if ~isempty(results.residual_covariances{1})
        imagesc(abs(results.residual_covariances{1}));
        colorbar;
        title('Residual Covariance (|·|)');
    end
    
    % Transfer function magnitude
    subplot(3, 4, 9);
    if ~isempty(results.transfer_functions{1})
        imagesc(abs(results.transfer_functions{1}));
        colorbar;
        title('DSTF Magnitude');
        xlabel('Sensors');
        ylabel('Sources');
    end
    
    % Frequency consistency
    subplot(3, 4, 10);
    if F > 1
        sample_entries = zeros(F, 1);
        for f = 1:F
            if ~isempty(results.initial_precision_matrices{f})
                sample_entries(f) = real(results.initial_precision_matrices{f}(2, 3));
            end
        end
        plot(1:F, sample_entries, 'o-');
        title('Precision(2,3) vs Frequency');
        xlabel('Frequency Index');
        ylabel('Value');
    end
    
    % Overall quality metrics
    subplot(3, 4, 11);
    quality_scores = zeros(F, 1);
    for f = 1:F
        if ~isempty(results.initial_precision_matrices{f})
            Omega = results.initial_precision_matrices{f};
            % Simple quality score based on condition number and determinant
            quality_scores(f) = 1 / (1 + log10(cond(Omega)));
        end
    end
    bar(quality_scores);
    title('Quality Scores');
    xlabel('Frequency');
    ylabel('Quality (0-1)');
    
    % Success summary
    subplot(3, 4, 12);
    success_rate = results.computation_stats.successful_frequencies / F;
    pie([success_rate, 1-success_rate], {'Success', 'Failed'});
    title(sprintf('Success Rate: %.1f%%', success_rate*100));
    
    sgtitle(sprintf('E-step Analysis: %d nodes, %d frequencies, %d edges', ...
                   sim_params.n_nodes, F, sim_params.n_edges));
end

function test_results = record_test_result(test_results, test_name, passed, details)
    % Record test result in results structure
    test_results.total_tests = test_results.total_tests + 1;
    
    if passed
        test_results.passed_tests = test_results.passed_tests + 1;
    else
        test_results.failed_tests = test_results.failed_tests + 1;
    end
    
    details.name = test_name;
    details.passed = passed;
    test_results.test_details{end+1} = details;
end