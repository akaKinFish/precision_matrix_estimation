function demo_results = demo_module7_enhanced_complex()
% DEMO_MODULE7_ENHANCED_COMPLEX - 测试增强版复数仿真功能
%
% 这个demo专门测试 module7_simulation_improved_complex 的新功能：
% 1. 复数 Hermitian 矩阵生成
% 2. 动态稀疏性变化
% 3. 平滑频率演变
% 4. 数值稳定性
%
% Usage:
%   demo_results = demo_module7_enhanced_complex()

    fprintf('========================================\n');
    fprintf('Module 7 Enhanced Complex Simulation Demo\n');
    fprintf('========================================\n\n');
    
    demo_results = struct();
    demo_results.timestamp = datestr(now);
    
    %% Test 1: 基本功能测试
    fprintf('=== Test 1: Basic Functionality ===\n');
    demo_results.test1 = test_basic_functionality();
    
    %% Test 2: 动态稀疏性测试
    fprintf('\n=== Test 2: Dynamic Sparsity Testing ===\n');
    demo_results.test2 = test_dynamic_sparsity();
    
    %% Test 3: 复数特性测试
    fprintf('\n=== Test 3: Complex Properties Testing ===\n');
    demo_results.test3 = test_complex_properties();
    
    %% Test 4: 参数敏感性测试
    fprintf('\n=== Test 4: Parameter Sensitivity Testing ===\n');
    demo_results.test4 = test_parameter_sensitivity();
    
    %% Test 5: 数值稳定性测试
    fprintf('\n=== Test 5: Numerical Stability Testing ===\n');
    demo_results.test5 = test_numerical_stability();
    
    %% Test 6: 可视化演示
    fprintf('\n=== Test 6: Visualization Demo ===\n');
    demo_results.test6 = create_comprehensive_visualization();
    
    %% 生成总结报告
    fprintf('\n=== Final Assessment ===\n');
    demo_results.summary = generate_final_assessment(demo_results);
    
    fprintf('\nModule 7 Enhanced Complex Demo Completed!\n');
    fprintf('Check demo_results for detailed analysis.\n');
end

function results = test_basic_functionality()
% 测试基本功能
    
    fprintf('Testing basic functionality with standard parameters...\n');
    
    results = struct();
    
    try
        [prec, cov, emp, params] = module7_simulation_improved_complex(...
            'n_nodes', 8, ...
            'n_freq', 12, ...
            'n_samples', 50, ...
            'graph_type', 'random', ...
            'edge_density', 0.4, ...
            'complex_strength', 0.8, ...
            'sparsity_variation', 0.8, ...
            'edge_activation_smoothness', 0.7, ...
            'random_seed', 123);
        
        results.success = true;
        results.params = params;
        
        % 基本验证
        results.n_frequencies = length(prec);
        results.matrix_size = size(prec{1});
        results.all_hermitian = all(cellfun(@ishermitian, prec));
        results.all_positive_definite = all(cellfun(@(x) all(real(eig(x)) > 1e-12), prec));
        results.has_complex = any(cellfun(@(x) any(~isreal(x(:))), prec));
        results.sparsity_changes = params.sparsity_changes;
        
        fprintf('  ✓ Generated %d matrices of size %dx%d\n', ...
                results.n_frequencies, results.matrix_size(1), results.matrix_size(2));
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
% 测试不同稀疏性变化参数
    
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
                'edge_activation_smoothness', 0.8, ...
                'random_seed', 100 + i);
            
            % 分析稀疏性
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
    
    % 验证趋势：稀疏性变化应该随参数增加
    fprintf('  Sparsity variation trend analysis:\n');
    changes = [];
    for i = 1:length(sparsity_levels)
        if isfield(results, ['level_' num2str(i)]) && ...
           isfield(results.(['level_' num2str(i)]), 'pattern_changes')
            changes(end+1) = results.(['level_' num2str(i)]).pattern_changes;
        end
    end
    
    if length(changes) >= 3 && all(diff(changes) >= 0)
        fprintf('  ✓ PASS: Sparsity changes increase with parameter\n');
        results.trend_correct = true;
    else
        fprintf('  ⚠ WARNING: Sparsity trend may not be monotonic\n');
        results.trend_correct = false;
    end
end

function results = test_complex_properties()
% 测试复数特性
    
    fprintf('Testing complex properties with different strengths...\n');
    
    results = struct();
    complex_strengths = [0.2, 0.5, 1.0, 1.5];
    
    for i = 1:length(complex_strengths)
        strength = complex_strengths(i);
        
        try
            [prec, ~, ~, params] = module7_simulation_improved_complex(...
                'n_nodes', 6, ...
                'n_freq', 8, ...
                'complex_strength', strength, ...
                'random_seed', 200 + i);
            
            % 分析复数特性
            complex_fractions = zeros(params.n_freq, 1);
            max_imaginary = zeros(params.n_freq, 1);
            phase_magnitudes = zeros(params.n_freq, 1);
            
            for f = 1:params.n_freq
                P = prec{f};
                
                % 复数元素比例
                complex_fractions(f) = sum(abs(imag(P(:))) > 1e-12) / numel(P);
                
                % 最大虚部
                max_imaginary(f) = max(abs(imag(P(:))));
                
                % 平均相位大小
                nonzero_elements = P(abs(P) > 1e-12);
                if ~isempty(nonzero_elements)
                    phases = angle(nonzero_elements);
                    phase_magnitudes(f) = mean(abs(phases));
                end
            end
            
            results.(['strength_' num2str(i)]) = struct();
            results.(['strength_' num2str(i)]).complex_strength = strength;
            results.(['strength_' num2str(i)]).avg_complex_fraction = mean(complex_fractions);
            results.(['strength_' num2str(i)]).avg_max_imaginary = mean(max_imaginary);
            results.(['strength_' num2str(i)]).avg_phase_magnitude = mean(phase_magnitudes);
            results.(['strength_' num2str(i)]).all_hermitian = all(cellfun(@ishermitian, prec));
            
            fprintf('  Strength %.1f: complex fraction %.2f, max imag %.3f, avg phase %.3f\n', ...
                    strength, mean(complex_fractions), mean(max_imaginary), mean(phase_magnitudes));
            
        catch ME
            fprintf('  ❌ Failed for complex_strength = %.1f: %s\n', strength, ME.message);
            results.(['strength_' num2str(i)]).error = ME.message;
        end
    end
    
    % 验证复数强度趋势
    complex_fracs = [];
    for i = 1:length(complex_strengths)
        if isfield(results, ['strength_' num2str(i)]) && ...
           isfield(results.(['strength_' num2str(i)]), 'avg_complex_fraction')
            complex_fracs(end+1) = results.(['strength_' num2str(i)]).avg_complex_fraction;
        end
    end
    
    if length(complex_fracs) >= 3 && all(diff(complex_fracs) >= -0.05) % Allow small decrease due to randomness
        fprintf('  ✓ PASS: Complex properties scale with strength parameter\n');
        results.strength_trend_correct = true;
    else
        fprintf('  ⚠ WARNING: Complex strength trend inconsistent\n');
        results.strength_trend_correct = false;
    end
end

function results = test_parameter_sensitivity()
% 测试参数敏感性
    
    fprintf('Testing parameter sensitivity...\n');
    
    results = struct();
    
    % 测试边激活平滑性
    smoothness_levels = [0.3, 0.6, 0.9];
    
    fprintf('  Testing edge activation smoothness:\n');
    for i = 1:length(smoothness_levels)
        smoothness = smoothness_levels(i);
        
        try
            [prec, ~, ~, params] = module7_simulation_improved_complex(...
                'n_nodes', 8, ...
                'n_freq', 15, ...
                'sparsity_variation', 0.4, ...
                'edge_activation_smoothness', smoothness, ...
                'random_seed', 300 + i);
            
            % 分析边值的平滑性
            if ~isempty(params.variable_edge_list)
                edge_idx = params.variable_edge_list(1, :);
                edge_values = zeros(params.n_freq, 1);
                for f = 1:params.n_freq
                    edge_values(f) = abs(prec{f}(edge_idx(1), edge_idx(2)));
                end
                
                % 计算二阶差分作为平滑性指标
                if length(edge_values) > 2
                    second_diffs = diff(diff(edge_values));
                    roughness = sum(second_diffs.^2);
                else
                    roughness = 0;
                end
                
                results.(['smoothness_' num2str(i)]) = struct();
                results.(['smoothness_' num2str(i)]).smoothness_param = smoothness;
                results.(['smoothness_' num2str(i)]).edge_roughness = roughness;
                
                fprintf('    Smoothness %.1f: roughness %.2e\n', smoothness, roughness);
            end
            
        catch ME
            fprintf('    ❌ Failed for smoothness = %.1f: %s\n', smoothness, ME.message);
        end
    end
    
    % 测试不同图类型
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
% 测试数值稳定性
    
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
            % 设置基础参数
            base_params = {...
                'n_nodes', 10, ...
                'n_freq', 12, ...
                'n_samples', 50, ...
                'random_seed', 500 + i ...
            };
            
            % 添加特定条件参数 - FIXED: 正确处理参数覆盖
            all_params = base_params;
            fields = fieldnames(condition);
            for j = 1:length(fields)
                if ~strcmp(fields{j}, 'name')
                    param_name = fields{j};
                    param_value = condition.(fields{j});
                    
                    % 查找并替换现有参数，或添加新参数
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
            
            % 稳定性检查
            stability_check = struct();
            stability_check.all_finite = all(cellfun(@(x) all(isfinite(x(:))), prec));
            stability_check.all_hermitian = all(cellfun(@ishermitian, prec));
            stability_check.all_psd = all(cellfun(@(x) all(real(eig(x)) > -1e-10), prec));
            stability_check.reasonable_condition = all(cellfun(@(x) cond(x) < 1e12, prec));
            
            % 协方差检查
            stability_check.cov_all_finite = all(cellfun(@(x) all(isfinite(x(:))), cov));
            stability_check.cov_all_psd = all(cellfun(@(x) all(real(eig(x)) > -1e-10), cov));
            
            % 经验协方差检查
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
% 创建综合可视化
    
    fprintf('Creating comprehensive visualization...\n');
    
    results = struct();
    
    try
        % 生成演示数据
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
        
        % 创建大型图表
        figure('Name', 'Module 7 Enhanced Complex Simulation', 'Position', [50, 50, 1400, 1000]);
        
        % 1. 稀疏模式演变
        subplot(3, 4, 1:2);
        frequencies_to_show = [1, 7, 13, 20];
        for i = 1:4
            f = frequencies_to_show(i);
            subplot(3, 4, i);
            pattern = abs(prec{f}) > 1e-10;
            imagesc(pattern);
            colormap(gray);
            title(sprintf('Sparsity Pattern - Freq %d', f));
            axis square;
        end
        
        % 5. 边数变化
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
        
        % 6. 复数成分分析
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
        
        % 7. 条件数分析
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
        
        % 8. 幅度和相位示例
        subplot(3, 4, 8);
        mid_freq = round(params.n_freq / 2);
        P_mid = prec{mid_freq};
        imagesc(abs(P_mid));
        colorbar;
        title(sprintf('Magnitude - Freq %d', mid_freq));
        axis square;
        
        % 9-12. 更多分析图
        subplot(3, 4, 9);
        imagesc(angle(P_mid));
        colorbar;
        colormap(gca, hsv);
        title(sprintf('Phase - Freq %d', mid_freq));
        axis square;
        
        % 特征值分布
        subplot(3, 4, 10);
        eigenvals = eig(P_mid);
        plot(real(eigenvals), imag(eigenvals), 'o', 'MarkerSize', 8, 'LineWidth', 2);
        xlabel('Real Part');
        ylabel('Imaginary Part');
        title('Eigenvalue Distribution');
        grid on;
        axis equal;
        
        % 边激活模式
        subplot(3, 4, 11);
        if ~isempty(params.variable_edge_list) && size(params.variable_activations, 2) > 0
            plot(1:params.n_freq, params.variable_activations(:, 1:min(3, end)), 'LineWidth', 2);
            xlabel('Frequency');
            ylabel('Activation Level');
            title('Variable Edge Activation');
            legend('Edge 1', 'Edge 2', 'Edge 3', 'Location', 'best');
            grid on;
        else
            text(0.5, 0.5, 'No Variable Edges', 'HorizontalAlignment', 'center');
            title('Edge Activation');
        end
        
        % 总结统计
        subplot(3, 4, 12);
        axis off;
        text(0.1, 0.9, 'Simulation Summary:', 'FontWeight', 'bold', 'FontSize', 12);
        text(0.1, 0.8, sprintf('Nodes: %d', params.n_nodes));
        text(0.1, 0.7, sprintf('Frequencies: %d', params.n_freq));
        text(0.1, 0.6, sprintf('Base edges: %d', params.n_base_edges));
        text(0.1, 0.5, sprintf('Variable edges: %d', params.n_variable_edges));
        text(0.1, 0.4, sprintf('Pattern changes: %d', params.sparsity_changes));
        text(0.1, 0.3, sprintf('Complex matrices: %d/%d', ...
                               params.complex_analysis.matrices_with_complex, params.n_freq));
        text(0.1, 0.2, sprintf('Avg complex fraction: %.2f', params.complex_analysis.avg_complex_fraction));
        text(0.1, 0.1, sprintf('Max imaginary: %.3f', params.complex_analysis.max_imaginary));
        
        sgtitle('Module 7 Enhanced Complex Simulation - Comprehensive Analysis');
        
        results.visualization_created = true;
        results.params = params;
        
        fprintf('  ✓ Comprehensive visualization created\n');
        
    catch ME
        fprintf('  ❌ Visualization failed: %s\n', ME.message);
        results.visualization_created = false;
        results.error = ME.message;
    end
end

function summary = generate_final_assessment(demo_results)
% 生成最终评估报告
    
    summary = struct();
    
    fprintf('Generating final assessment...\n');
    
    % 收集所有测试结果
    tests = {'test1', 'test2', 'test3', 'test4', 'test5', 'test6'};
    test_names = {'Basic Functionality', 'Dynamic Sparsity', 'Complex Properties', ...
                  'Parameter Sensitivity', 'Numerical Stability', 'Visualization'};
    
    passed_tests = 0;
    total_tests = 0;
    
    for i = 1:length(tests)
        test_name = tests{i};
        if isfield(demo_results, test_name)
            total_tests = total_tests + 1;
            
            % 判断测试是否通过
            test_result = demo_results.(test_name);
            if isfield(test_result, 'success') && test_result.success
                passed_tests = passed_tests + 1;
            elseif isfield(test_result, 'overall_pass') && test_result.overall_pass
                passed_tests = passed_tests + 1;
            elseif isfield(test_result, 'visualization_created') && test_result.visualization_created
                passed_tests = passed_tests + 1;
            elseif isfield(test_result, 'test_completed') && test_result.test_completed
                passed_tests = passed_tests + 1;
            end
        end
    end
    
    summary.total_tests = total_tests;
    summary.passed_tests = passed_tests;
    summary.success_rate = passed_tests / total_tests * 100;
    
    % 检查关键功能
    summary.requirements_met = struct();
    
    if isfield(demo_results, 'test1') && isfield(demo_results.test1, 'has_complex')
        summary.requirements_met.complex_hermitian = demo_results.test1.has_complex && demo_results.test1.all_hermitian;
    else
        summary.requirements_met.complex_hermitian = false;
    end
    
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
    
    % 生成推荐
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
% 辅助函数：将逻辑值转换为是/否字符串
    if logical_value
        str = 'YES';
    else
        str = 'NO';
    end
end