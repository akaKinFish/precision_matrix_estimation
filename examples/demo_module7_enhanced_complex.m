function demo_results = demo_module7_enhanced_complex()
% DEMO_MODULE7_ENHANCED_COMPLEX - Test enhanced complex simulation functionality
%
% This demo specifically tests the module7_simulation_improved_complex features:
% 1. Complex Hermitian matrix generation
% 2. Dynamic sparsity variation
% 3. Smooth frequency evolution  
% 4. Numerical stability
%
% Usage:
%   demo_results = demo_module7_enhanced_complex()

    fprintf('========================================\n');
    fprintf('Module 7 Enhanced Complex Simulation Demo\n');
    fprintf('========================================\n\n');
    
    demo_results = struct();
    demo_results.timestamp = datestr(now);
    
    %% Test 1: Basic functionality test
    fprintf('=== Test 1: Basic Functionality ===\n');
    demo_results.test1 = test_basic_functionality();
    
    %% Test 2: Dynamic sparsity test
    fprintf('\n=== Test 2: Dynamic Sparsity Testing ===\n');
    demo_results.test2 = test_dynamic_sparsity();
    
    %% Test 3: Complex properties test
    fprintf('\n=== Test 3: Complex Properties Testing ===\n');
    demo_results.test3 = test_complex_properties();
    
    %% Test 4: Parameter sensitivity test
    fprintf('\n=== Test 4: Parameter Sensitivity Testing ===\n');
    demo_results.test4 = test_parameter_sensitivity();
    
    %% Test 5: Numerical stability test
    fprintf('\n=== Test 5: Numerical Stability Testing ===\n');
    demo_results.test5 = test_numerical_stability();
    
    %% Test 6: Visualization demo
    fprintf('\n=== Test 6: Visualization Demo ===\n');
    demo_results.test6 = create_comprehensive_visualization();
    
    %% Generate summary report
    fprintf('\n=== Final Assessment ===\n');
    demo_results.summary = generate_final_assessment(demo_results);
    
    fprintf('\nModule 7 Enhanced Complex Demo Completed!\n');
    fprintf('Check demo_results for detailed analysis.\n');
end

function results = test_basic_functionality()
% Test basic functionality
    
    fprintf('Testing basic functionality with standard parameters...\n');
    
    results = struct();
    results.success = true;
    
    try
        % Generate test data with improved parameters
        [prec, cov, emp, params] = module7_simulation_improved_complex(...
            'n_nodes', 8, ...
            'n_freq', 12, ...
            'n_samples', 100, ...
            'graph_type', 'random', ...
            'edge_density', 0.3, ...
            'sparsity_variation', 0.4, ...  % Increased for more changes
            'edge_activation_smoothness', 0.6, ... % Reduced for more variation
            'complex_strength', 1.0, ...
            'random_seed', 42);
        
        % Basic validation
        results.n_matrices = length(prec);
        results.matrix_size = size(prec{1});
        
        % Check Hermitian property
        hermitian_check = true;
        for f = 1:params.n_freq
            if ~ishermitian(prec{f})
                hermitian_check = false;
                break;
            end
        end
        results.all_hermitian = hermitian_check;
        
        % Check positive definiteness
        pd_check = true;
        for f = 1:params.n_freq
            eigvals = eig(prec{f});
            if any(real(eigvals) <= 1e-10)
                pd_check = false;
                break;
            end
        end
        results.all_positive_definite = pd_check;
        
        % Check for complex values
        has_complex = false;
        total_complex_elements = 0;
        for f = 1:params.n_freq
            complex_elements = sum(abs(imag(prec{f}(:))) > 1e-12);
            total_complex_elements = total_complex_elements + complex_elements;
            if complex_elements > 0
                has_complex = true;
            end
        end
        results.has_complex = has_complex;
        results.avg_complex_elements = total_complex_elements / params.n_freq;
        
        % Check sparsity changes
        results.sparsity_changes = params.sparsity_changes;
        results.base_edges = params.n_base_edges;
        results.variable_edges = params.n_variable_edges;
        
        % Output summary
        fprintf('  ✓ Generated %d matrices of size %dx%d\n', ...
                results.n_matrices, results.matrix_size(1), results.matrix_size(2));
        fprintf('  ✓ All matrices Hermitian: %s\n', yesno(results.all_hermitian));
        fprintf('  ✓ All matrices positive definite: %s\n', yesno(results.all_positive_definite));
        fprintf('  ✓ Contains complex values: %s\n', yesno(results.has_complex));
        fprintf('  ✓ Sparsity pattern changes: %d\n', results.sparsity_changes);
        
        if results.all_hermitian && results.all_positive_definite && ...
           results.has_complex && results.sparsity_changes > 0
            fprintf('  🎉 PASS: All basic requirements met!\n');
            results.overall_pass = true;
        else
            fprintf('  ❌ FAIL: Some requirements not met\n');
            results.overall_pass = false;
        end
        
    catch ME
        fprintf('  ❌ ERROR: %s\n', ME.message);
        results.success = false;
        results.error = ME.message;
        results.overall_pass = false;
    end
end

function results = test_dynamic_sparsity()
% Test different sparsity variation parameters
    
    fprintf('Testing different sparsity variation levels...\n');
    
    results = struct();
    sparsity_levels = [0, 0.2, 0.5, 0.8];
    
    for i = 1:length(sparsity_levels)
        sparsity_var = sparsity_levels(i);
        
        try
            [prec, ~, ~, params] = module7_simulation_improved_complex(...
                'n_nodes', 6, ...
                'n_freq', 10, ...
                'n_samples', 30, ...
                'sparsity_variation', sparsity_var, ...
                'edge_activation_smoothness', 0.5, ... % Lower for more variation
                'random_seed', 100 + i);
            
            % Analyze sparsity
            edge_counts = zeros(params.n_freq, 1);
            for f = 1:params.n_freq
                P = prec{f};
                P_upper = triu(P, 1);
                edge_counts(f) = sum(abs(P_upper(:)) > 1e-10);
            end
            
            results.(['level_' num2str(i)]) = struct();
            results.(['level_' num2str(i)]).sparsity_variation = sparsity_var;
            results.(['level_' num2str(i)]).pattern_changes = params.sparsity_changes;
            results.(['level_' num2str(i)]).edge_counts = edge_counts;
            results.(['level_' num2str(i)]).edge_std = std(edge_counts);
            results.(['level_' num2str(i)]).edge_range = [min(edge_counts), max(edge_counts)];
            
            fprintf('  Sparsity variation %.1f: %d changes, edge std %.1f, range [%d,%d]\n', ...
                    sparsity_var, params.sparsity_changes, std(edge_counts), ...
                    min(edge_counts), max(edge_counts));
            
        catch ME
            fprintf('  ❌ Failed for sparsity_variation = %.1f: %s\n', sparsity_var, ME.message);
            results.(['level_' num2str(i)]).error = ME.message;
        end
    end
    
    % Verify trend: sparsity changes should increase with parameter
    fprintf('  Sparsity variation trend analysis:\n');
    changes = [];
    for i = 1:length(sparsity_levels)
        if isfield(results, ['level_' num2str(i)]) && ...
           isfield(results.(['level_' num2str(i)]), 'pattern_changes')
            changes(end+1) = results.(['level_' num2str(i)]).pattern_changes;
        end
    end
    
    if length(changes) >= 3 && changes(end) > changes(1)
        fprintf('  ✓ PASS: Sparsity changes increase with parameter\n');
        results.trend_correct = true;
    else
        fprintf('  ⚠ WARNING: Sparsity trend may be inconsistent\n');
        results.trend_correct = false;
    end
end

function results = test_complex_properties()
% Test complex properties with different strengths
    
    fprintf('Testing complex properties with different strengths...\n');
    
    results = struct();
    complex_strengths = [0.2, 0.5, 1.0, 1.5];
    
    for i = 1:length(complex_strengths)
        strength = complex_strengths(i);
        
        try
            [prec, ~, ~, params] = module7_simulation_improved_complex(...
                'n_nodes', 6, ...
                'n_freq', 8, ...
                'n_samples', 40, ...
                'complex_strength', strength, ...
                'coefficient_complex_fraction', 1.0, ...
                'random_seed', 200 + i);
            
            % Analyze complex properties
            complex_fraction = params.avg_complex_fraction;
            max_imag = params.max_imag_component;
            avg_phase = params.avg_phase_magnitude;
            
            results.(['strength_' num2str(i)]) = struct();
            results.(['strength_' num2str(i)]).complex_strength = strength;
            results.(['strength_' num2str(i)]).complex_fraction = complex_fraction;
            results.(['strength_' num2str(i)]).max_imag = max_imag;
            results.(['strength_' num2str(i)]).avg_phase = avg_phase;
            
            fprintf('  Strength %.1f: complex fraction %.2f, max imag %.3f, avg phase %.3f\n', ...
                    strength, complex_fraction, max_imag, avg_phase);
            
        catch ME
            fprintf('  ❌ Failed for complex_strength = %.1f: %s\n', strength, ME.message);
            results.(['strength_' num2str(i)]).error = ME.message;
        end
    end
    
    % Check if max imaginary component generally increases with strength
    max_imags = [];
    for i = 1:length(complex_strengths)
        if isfield(results, ['strength_' num2str(i)]) && ...
           isfield(results.(['strength_' num2str(i)]), 'max_imag')
            max_imags(end+1) = results.(['strength_' num2str(i)]).max_imag;
        end
    end
    
    if length(max_imags) >= 3
        % Check if there's generally increasing trend (allowing some variation)
        trend_score = corr((1:length(max_imags))', max_imags');
        if trend_score > 0.3
            fprintf('  ✓ PASS: Complex strength shows reasonable trend\n');
            results.strength_trend_ok = true;
        else
            fprintf('  ⚠ WARNING: Complex strength trend inconsistent\n');
            results.strength_trend_ok = false;
        end
    end
end

function results = test_parameter_sensitivity()
% Test parameter sensitivity
    
    fprintf('Testing parameter sensitivity...\n');
    
    results = struct();
    
    % Test edge activation smoothness
    smoothness_levels = [0.3, 0.6, 0.9];
    
    fprintf('  Testing edge activation smoothness:\n');
    for i = 1:length(smoothness_levels)
        smoothness = smoothness_levels(i);
        
        try
            [prec, ~, ~, params] = module7_simulation_improved_complex(...
                'n_nodes', 8, ...
                'n_freq', 15, ...
                'n_samples', 60, ...
                'sparsity_variation', 0.3, ...
                'edge_activation_smoothness', smoothness, ...
                'random_seed', 300 + i);
            
            % Calculate roughness metric for edge activation
            if params.n_variable_edges > 0
                activations = params.variable_activations;
                roughness = 0;
                for e = 1:params.n_variable_edges
                    activation_curve = activations(:, e);
                    diff_curve = diff(activation_curve);
                    roughness = roughness + mean(abs(diff_curve));
                end
                roughness = roughness / params.n_variable_edges;
            else
                roughness = 0;
            end
            
            results.(['smoothness_' num2str(i)]) = struct();
            results.(['smoothness_' num2str(i)]).smoothness_param = smoothness;
            results.(['smoothness_' num2str(i)]).edge_roughness = roughness;
            
            fprintf('    Smoothness %.1f: roughness %.2e\n', smoothness, roughness);
            
        catch ME
            fprintf('    ❌ Failed for smoothness = %.1f: %s\n', smoothness, ME.message);
        end
    end
    
    % Test different graph types
    graph_types = {'random', 'chain', 'hub'};
    
    fprintf('  Testing different graph types:\n');
    for i = 1:length(graph_types)
        graph_type = graph_types{i};
        
        try
            [prec, ~, ~, params] = module7_simulation_improved_complex(...
                'n_nodes', 8, ...
                'n_freq', 10, ...
                'graph_type', graph_type, ...
                'edge_density', 0.3, ...
                'random_seed', 400 + i);
            
            results.(['graph_' num2str(i)]) = struct();
            results.(['graph_' num2str(i)]).graph_type = graph_type;
            results.(['graph_' num2str(i)]).n_base_edges = params.n_base_edges;
            results.(['graph_' num2str(i)]).n_variable_edges = params.n_variable_edges;
            results.(['graph_' num2str(i)]).sparsity_changes = params.sparsity_changes;
            
            fprintf('    %s: %d base edges, %d variable, %d changes\n', ...
                    graph_type, params.n_base_edges, params.n_variable_edges, ...
                    params.sparsity_changes);
            
        catch ME
            fprintf('    ❌ Failed for graph_type = %s: %s\n', graph_type, ME.message);
        end
    end
    
    results.test_completed = true;
end

function results = test_numerical_stability()
% Test numerical stability
    
    fprintf('Testing numerical stability under various conditions...\n');
    
    results = struct();
    test_conditions = {
        struct('name', 'Small matrices', 'n_nodes', 4, 'n_freq', 5),
        struct('name', 'Large matrices', 'n_nodes', 20, 'n_freq', 30),
        struct('name', 'High sparsity', 'sparsity_variation', 0.8, 'edge_density', 0.8),
        struct('name', 'High complexity', 'complex_strength', 2.0, 'n_basis', 10)
    };
    
    for i = 1:length(test_conditions)
        condition = test_conditions{i};
        
        fprintf('  Testing %s...\n', condition.name);
        
        try
            % Set base parameters
            base_params = {...
                'n_nodes', 10, ...
                'n_freq', 12, ...
                'n_samples', 50, ...
                'random_seed', 500 + i ...
            };
            
            % Add specific condition parameters
            all_params = base_params;
            fields = fieldnames(condition);
            for j = 1:length(fields)
                if ~strcmp(fields{j}, 'name')
                    param_name = fields{j};
                    param_value = condition.(fields{j});
                    
                    % Find and replace existing parameter, or add new parameter
                    param_found = false;
                    for k = 1:2:length(all_params)
                        if strcmp(all_params{k}, param_name)
                            all_params{k+1} = param_value;
                            param_found = true;
                            break;
                        end
                    end
                    
                    if ~param_found
                        all_params{end+1} = param_name;
                        all_params{end+1} = param_value;
                    end
                end
            end
            
            [prec, cov, emp, params] = module7_simulation_improved_complex(all_params{:});
            
            % Stability checks
            stability_check = struct();
            stability_check.all_finite = all(cellfun(@(x) all(isfinite(x(:))), prec));
            stability_check.all_hermitian = all(cellfun(@ishermitian, prec));
            stability_check.all_psd = all(cellfun(@(x) all(real(eig(x)) > -1e-10), prec));
            stability_check.reasonable_condition = all(cellfun(@(x) cond(x) < 1e12, prec));
            
            % Covariance checks
            stability_check.cov_all_finite = all(cellfun(@(x) all(isfinite(x(:))), cov));
            stability_check.cov_all_psd = all(cellfun(@(x) all(real(eig(x)) > -1e-10), cov));
            
            % Empirical covariance checks
            stability_check.emp_all_finite = all(cellfun(@(x) all(isfinite(x(:))), emp));
            
            results.(['condition_' num2str(i)]) = stability_check;
            results.(['condition_' num2str(i)]).name = condition.name;
            
            if all(struct2array(stability_check))
                fprintf('    ✓ PASS: All stability checks passed\n');
            else
                fprintf('    ⚠ WARNING: Some stability issues detected\n');
                failed_checks = fieldnames(stability_check);
                for k = 1:length(failed_checks)
                    if ~stability_check.(failed_checks{k})
                        fprintf('      - Failed: %s\n', failed_checks{k});
                    end
                end
            end
            
        catch ME
            fprintf('    ❌ ERROR: %s\n', ME.message);
            results.(['condition_' num2str(i)]).error = ME.message;
        end
    end
end

function results = create_comprehensive_visualization()
% Create comprehensive visualization
    
    fprintf('Creating comprehensive visualization...\n');
    
    results = struct();
    
    try
        % Generate demo data
        [prec, cov, emp, params] = module7_simulation_improved_complex(...
            'n_nodes', 12, ...
            'n_freq', 20, ...
            'n_samples', 80, ...
            'graph_type', 'random', ...
            'edge_density', 0.25, ...
            'complex_strength', 1.0, ...
            'sparsity_variation', 0.4, ...
            'edge_activation_smoothness', 0.75, ...
            'random_seed', 999);
        
        % Create comprehensive figure
        figure('Name', 'Module 7 Enhanced Complex Simulation - Comprehensive Analysis', ...
               'Position', [50, 50, 1400, 1000]);
        
        % 1. Sparsity pattern evolution
        frequencies_to_show = [1, 7, 13, 20];
        for i = 1:4
            f = frequencies_to_show(i);
            subplot(3, 4, i);
            pattern = abs(prec{f}) > 1e-10;
            imagesc(pattern);
            colormap(gca, gray);
            title(sprintf('Sparsity Pattern - Freq %d', f));
            axis square;
        end
        
        % 5. Edge count evolution
        subplot(3, 4, 5);
        edge_counts = zeros(params.n_freq, 1);
        for f = 1:params.n_freq
            P_upper = triu(prec{f}, 1);
            edge_counts(f) = sum(abs(P_upper(:)) > 1e-10);
        end
        plot(1:params.n_freq, edge_counts, 'o-', 'LineWidth', 2);
        xlabel('Frequency');
        ylabel('Number of Edges');
        title('Edge Count Evolution');
        grid on;
        
        % 6. Complex content analysis
        subplot(3, 4, 6);
        complex_fractions = zeros(params.n_freq, 1);
        for f = 1:params.n_freq
            P = prec{f};
            complex_fractions(f) = sum(abs(imag(P(:))) > 1e-12) / numel(P);
        end
        plot(1:params.n_freq, complex_fractions, 's-', 'LineWidth', 2, 'Color', [0.8, 0.2, 0.2]);
        xlabel('Frequency');
        ylabel('Complex Element Fraction');
        title('Complex Content Evolution');
        grid on;
        
        % 7. Numerical conditioning
        subplot(3, 4, 7);
        condition_numbers = zeros(params.n_freq, 1);
        for f = 1:params.n_freq
            condition_numbers(f) = cond(prec{f});
        end
        semilogy(1:params.n_freq, condition_numbers, '^-', 'LineWidth', 2, 'Color', [0.2, 0.6, 0.8]);
        xlabel('Frequency');
        ylabel('Condition Number (log scale)');
        title('Numerical Conditioning');
        grid on;
        
        % 8. Eigenvalue distribution for one frequency
        subplot(3, 4, 8);
        f_sample = 10;
        eigvals = eig(prec{f_sample});
        plot(real(eigvals), imag(eigvals), 'o', 'MarkerSize', 8, 'LineWidth', 2);
        xlabel('Real Part');
        ylabel('Imaginary Part');
        title(sprintf('Eigenvalue Distribution - Freq %d', f_sample));
        grid on;
        axis equal;
        
        % 9. Phase analysis for one frequency
        subplot(3, 4, 9);
        P_sample = prec{f_sample};
        phase_matrix = angle(P_sample);
        imagesc(phase_matrix);
        colorbar;
        colormap(gca, hsv);
        title(sprintf('Phase - Freq %d', f_sample));
        axis square;
        caxis([-pi, pi]);
        
        % 10. Magnitude for same frequency
        subplot(3, 4, 10);
        imagesc(abs(P_sample));
        colorbar;
        title(sprintf('Magnitude - Freq %d', f_sample));
        axis square;
        
        % 11. Variable edge activation
        subplot(3, 4, 11);
        if params.n_variable_edges > 0
            for e = 1:min(3, params.n_variable_edges)
                plot(1:params.n_freq, params.variable_activations(:, e), ...
                     'LineWidth', 2, 'DisplayName', sprintf('Edge %d', e));
                hold on;
            end
            xlabel('Frequency');
            ylabel('Activation Level');
            title('Variable Edge Activation');
            legend('show');
            grid on;
        else
            text(0.5, 0.5, 'No Variable Edges', 'HorizontalAlignment', 'center');
            title('Variable Edge Activation');
        end
        
        % 12. Summary statistics
        subplot(3, 4, 12);
        axis off;
        summary_text = {
            'Simulation Summary:',
            sprintf('Nodes: %d', params.n_nodes),
            sprintf('Frequencies: %d', params.n_freq),
            sprintf('Base edges: %d', params.n_base_edges),
            sprintf('Variable edges: %d', params.n_variable_edges),
            sprintf('Pattern changes: %d', params.sparsity_changes),
            sprintf('Complex matrices: %d/%d', params.matrices_with_complex, params.n_freq),
            sprintf('Avg complex fraction: %.2f', params.avg_complex_fraction),
            sprintf('Max imaginary: %.3f', params.max_imag_component)
        };
        
        for i = 1:length(summary_text)
            text(0.1, 0.9 - (i-1)*0.08, summary_text{i}, 'FontSize', 10, ...
                 'VerticalAlignment', 'top', 'FontWeight', 'bold');
        end
        
        results.visualization_created = true;
        fprintf('  ✓ Comprehensive visualization created\n');
        
    catch ME
        fprintf('  ❌ Visualization failed: %s\n', ME.message);
        results.error = ME.message;
        results.visualization_created = false;
    end
end

function summary = generate_final_assessment(demo_results)
% Generate final assessment
    
    fprintf('Generating final assessment...\n');
    
    summary = struct();
    
    % Count passed tests
    passed_tests = 0;
    total_tests = 6;
    
    test_names = {'test1', 'test2', 'test3', 'test4', 'test5', 'test6'};
    
    for i = 1:length(test_names)
        test_name = test_names{i};
        if isfield(demo_results, test_name)
            test_result = demo_results.(test_name);
            if isfield(test_result, 'overall_pass') && test_result.overall_pass
                passed_tests = passed_tests + 1;
            elseif isfield(test_result, 'test_completed') && test_result.test_completed
                passed_tests = passed_tests + 1;
            elseif isfield(test_result, 'visualization_created') && test_result.visualization_created
                passed_tests = passed_tests + 1;
            end
        end
    end
    
    summary.passed_tests = passed_tests;
    summary.total_tests = total_tests;
    summary.success_rate = (passed_tests / total_tests) * 100;
    
    % Check core requirements
    summary.requirements_met = struct();
    
    % Complex Hermitian requirement
    if isfield(demo_results, 'test1') && isfield(demo_results.test1, 'all_hermitian') && ...
       isfield(demo_results.test1, 'has_complex')
        summary.requirements_met.complex_hermitian = demo_results.test1.all_hermitian && ...
                                                     demo_results.test1.has_complex;
    else
        summary.requirements_met.complex_hermitian = false;
    end
    
    % Dynamic sparsity requirement
    if isfield(demo_results, 'test1') && isfield(demo_results.test1, 'sparsity_changes')
        summary.requirements_met.dynamic_sparsity = demo_results.test1.sparsity_changes > 0;
    else
        summary.requirements_met.dynamic_sparsity = false;
    end
    
    if isfield(demo_results, 'test1') && isfield(demo_results.test1, 'all_positive_definite')
        summary.requirements_met.numerical_stability = demo_results.test1.all_positive_definite;
    else
        summary.requirements_met.numerical_stability = false;
    end
    
    % Generate recommendations
    all_requirements_met = all(struct2array(summary.requirements_met));
    
    fprintf('\n========================================\n');
    fprintf('FINAL ASSESSMENT REPORT\n');
    fprintf('========================================\n');
    fprintf('Tests passed: %d/%d (%.1f%%)\n', passed_tests, total_tests, summary.success_rate);
    fprintf('\nCore Requirements:\n');
    fprintf('  ✓ Complex Hermitian matrices: %s\n', yesno(summary.requirements_met.complex_hermitian));
    fprintf('  ✓ Dynamic sparsity patterns: %s\n', yesno(summary.requirements_met.dynamic_sparsity));
    fprintf('  ✓ Numerical stability: %s\n', yesno(summary.requirements_met.numerical_stability));
    
    if all_requirements_met && summary.success_rate >= 80
        fprintf('\n🎉 OVERALL RESULT: SUCCESS!\n');
        fprintf('   Module7 enhanced complex simulation is ready for use.\n');
        fprintf('   All core requirements are satisfied.\n');
        summary.overall_recommendation = 'APPROVED';
    elseif all_requirements_met
        fprintf('\n⚠️  OVERALL RESULT: CONDITIONAL SUCCESS\n');
        fprintf('   Core requirements met but some tests had issues.\n');
        fprintf('   Recommend reviewing failed tests before production use.\n');
        summary.overall_recommendation = 'CONDITIONAL';
    else
        fprintf('\n❌ OVERALL RESULT: NEEDS WORK\n');
        fprintf('   Core requirements not fully satisfied.\n');
        fprintf('   Module needs fixes before integration with Module1.\n');
        summary.overall_recommendation = 'NEEDS_FIXES';
    end
    
    summary.all_requirements_met = all_requirements_met;
    
    fprintf('========================================\n');
end

function str = yesno(logical_value)
% Helper function: convert logical value to yes/no string
    if logical_value
        str = 'YES';
    else
        str = 'NO';
    end
end
                