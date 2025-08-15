function test_results = test_simulation_functions()
% TEST_SIMULATION_FUNCTIONS - 验证三个simulation函数的特性
%
% 这个脚本测试三个simulation函数，帮助确定哪个最符合需求：
% 1. 生成Hermitian矩阵
% 2. 稀疏模式有变化
%
% Usage:
%   test_results = test_simulation_functions()

    fprintf('========================================\n');
    fprintf('Testing All Simulation Functions\n');
    fprintf('========================================\n\n');
    
    test_results = struct();
    
    % 通用测试参数
    test_params = {
        'n_nodes', 10, ...
        'n_freq', 15, ...
        'n_samples', 50, ...
        'graph_type', 'random', ...
        'edge_density', 0.3, ...
        'random_seed', 123
    };
    
    %% 测试 1: module7_simulation (原始版本)
    fprintf('=== Testing module7_simulation (Original) ===\n');
    try
        [prec1, cov1, emp1, params1] = module7_simulation(test_params{:});
        
        test_results.original = analyze_simulation_output(prec1, cov1, emp1, params1, 'Original');
        fprintf('✓ Original version completed\n\n');
        
    catch ME
        fprintf('✗ Original version failed: %s\n\n', ME.message);
        test_results.original = struct('success', false, 'error', ME.message);
    end
    
    %% 测试 2: module7_simulation_improved (改进版本)
    fprintf('=== Testing module7_simulation_improved (Improved) ===\n');
    try
        [prec2, cov2, emp2, params2] = module7_simulation_improved(...
            test_params{:}, ...
            'sparsity_variation', 0.3, ...
            'edge_activation_smoothness', 0.7);
        
        test_results.improved = analyze_simulation_output(prec2, cov2, emp2, params2, 'Improved');
        fprintf('✓ Improved version completed\n\n');
        
    catch ME
        fprintf('✗ Improved version failed: %s\n\n', ME.message);
        test_results.improved = struct('success', false, 'error', ME.message);
    end
    
    %% 测试 3: module7_simulation_improved_complex (复数版本)
    fprintf('=== Testing module7_simulation_improved_complex (Complex) ===\n');
    try
        [prec3, cov3, emp3, params3] = module7_simulation_improved_complex(...
            test_params{:}, ...
            'complex_strength', 1.0);
        
        test_results.complex = analyze_simulation_output(prec3, cov3, emp3, params3, 'Complex');
        fprintf('✓ Complex version completed\n\n');
        
    catch ME
        fprintf('✗ Complex version failed: %s\n\n', ME.message);
        test_results.complex = struct('success', false, 'error', ME.message);
    end
    
    %% 生成对比报告
    fprintf('=== Comparison Summary ===\n');
    test_results.recommendation = generate_recommendation(test_results);
    
    %% 可视化对比（如果有成功的结果）
    create_comparison_visualization(test_results);
    
    fprintf('\nTest completed! Check test_results for detailed analysis.\n');
end

function analysis = analyze_simulation_output(prec, cov, emp, params, version_name)
% 分析simulation输出的特性
    
    analysis = struct();
    analysis.success = true;
    analysis.version = version_name;
    
    F = length(prec);
    n = size(prec{1}, 1);
    
    fprintf('Analyzing %s version:\n', version_name);
    
    %% 1. 基本信息
    analysis.n_nodes = n;
    analysis.n_frequencies = F;
    analysis.matrix_type = 'unknown';
    
    %% 2. 检查是否为Hermitian
    hermitian_count = 0;
    complex_count = 0;
    
    for f = 1:F
        % 检查Hermitian性质
        if ishermitian(prec{f})
            hermitian_count = hermitian_count + 1;
        end
        
        % 检查是否包含复数
        if ~isreal(prec{f})
            complex_count = complex_count + 1;
        end
    end
    
    analysis.hermitian_matrices = hermitian_count;
    analysis.hermitian_percentage = hermitian_count / F * 100;
    analysis.complex_matrices = complex_count;
    analysis.complex_percentage = complex_count / F * 100;
    
    if complex_count > 0
        analysis.matrix_type = 'complex';
    else
        analysis.matrix_type = 'real';
    end
    
    fprintf('  Matrix type: %s\n', analysis.matrix_type);
    fprintf('  Hermitian matrices: %d/%d (%.1f%%)\n', hermitian_count, F, analysis.hermitian_percentage);
    
    %% 3. 稀疏性分析
    sparsity_patterns = zeros(F, 1);
    edge_counts = zeros(F, 1);
    
    for f = 1:F
        % 计算非零元素（排除对角线）
        P_f = prec{f};
        P_f(1:n+1:end) = 0; % 移除对角线
        edge_counts(f) = sum(abs(P_f(:)) > 1e-10) / 2; % 上三角的边数
        sparsity_patterns(f) = edge_counts(f);
    end
    
    analysis.edge_counts = edge_counts;
    analysis.mean_edges = mean(edge_counts);
    analysis.std_edges = std(edge_counts);
    analysis.min_edges = min(edge_counts);
    analysis.max_edges = max(edge_counts);
    analysis.sparsity_variation = analysis.std_edges / analysis.mean_edges;
    
    fprintf('  Average edges: %.1f ± %.1f\n', analysis.mean_edges, analysis.std_edges);
    fprintf('  Edge range: [%d, %d]\n', analysis.min_edges, analysis.max_edges);
    fprintf('  Sparsity variation: %.3f\n', analysis.sparsity_variation);
    
    %% 4. 正定性检查
    min_eigenvals = zeros(F, 1);
    condition_numbers = zeros(F, 1);
    
    for f = 1:F
        try
            eigenvals = eig(prec{f});
            min_eigenvals(f) = min(real(eigenvals));
            condition_numbers(f) = max(real(eigenvals)) / min(real(eigenvals));
        catch
            min_eigenvals(f) = NaN;
            condition_numbers(f) = Inf;
        end
    end
    
    analysis.min_eigenvalues = min_eigenvals;
    analysis.condition_numbers = condition_numbers;
    analysis.positive_definite_count = sum(min_eigenvals > 1e-12);
    analysis.well_conditioned_count = sum(condition_numbers < 1e6);
    
    fprintf('  Positive definite: %d/%d\n', analysis.positive_definite_count, F);
    fprintf('  Well conditioned: %d/%d\n', analysis.well_conditioned_count, F);
    
    %% 5. 检查频率平滑性
    if F > 2
        % 追踪一个off-diagonal元素的变化
        element_trajectory = zeros(F, 1);
        for f = 1:F
            element_trajectory(f) = prec{f}(2, 1); % 元素(2,1)
        end
        
        % 计算二阶差分来衡量平滑性
        if length(element_trajectory) > 2
            second_diff = diff(diff(element_trajectory));
            analysis.smoothness_metric = sum(second_diff.^2);
        else
            analysis.smoothness_metric = 0;
        end
        
        fprintf('  Frequency smoothness metric: %.2e\n', analysis.smoothness_metric);
    end
end

function recommendation = generate_recommendation(test_results)
% 基于测试结果生成推荐
    
    recommendation = struct();
    recommendation.recommended_function = 'none';
    recommendation.reasons = {};
    recommendation.scores = struct();
    
    % 评分标准
    functions = {'original', 'improved', 'complex'};
    scores = zeros(length(functions), 1);
    
    for i = 1:length(functions)
        func_name = functions{i};
        if ~isfield(test_results, func_name) || ~test_results.(func_name).success
            scores(i) = 0;
            continue;
        end
        
        result = test_results.(func_name);
        score = 0;
        
        % 1. Hermitian矩阵 (30分)
        if result.hermitian_percentage >= 95
            score = score + 30;
        elseif result.hermitian_percentage >= 80
            score = score + 20;
        elseif result.hermitian_percentage >= 50
            score = score + 10;
        end
        
        % 2. 稀疏性变化 (25分)
        if result.sparsity_variation > 0.1
            score = score + 25;
        elseif result.sparsity_variation > 0.05
            score = score + 15;
        elseif result.sparsity_variation > 0.01
            score = score + 10;
        end
        
        % 3. 数值稳定性 (25分)
        if result.positive_definite_count == result.n_frequencies
            score = score + 15;
        end
        if result.well_conditioned_count >= result.n_frequencies * 0.9
            score = score + 10;
        end
        
        % 4. 复数支持 (20分) - 如果需要复数的话
        if strcmp(result.matrix_type, 'complex') && result.complex_percentage > 50
            score = score + 20;
        elseif strcmp(result.matrix_type, 'real')
            score = score + 10; % 实数也是可接受的
        end
        
        scores(i) = score;
        recommendation.scores.(func_name) = score;
    end
    
    % 选择最高分的函数
    [max_score, best_idx] = max(scores);
    
    if max_score > 50
        recommendation.recommended_function = functions{best_idx};
        
        % 生成推荐理由
        best_result = test_results.(functions{best_idx});
        
        if best_result.hermitian_percentage >= 95
            recommendation.reasons{end+1} = sprintf('Generates %.1f%% Hermitian matrices', best_result.hermitian_percentage);
        end
        
        if best_result.sparsity_variation > 0.1
            recommendation.reasons{end+1} = sprintf('Good sparsity variation (%.3f)', best_result.sparsity_variation);
        end
        
        if best_result.positive_definite_count == best_result.n_frequencies
            recommendation.reasons{end+1} = 'All matrices are positive definite';
        end
        
        if strcmp(best_result.matrix_type, 'complex')
            recommendation.reasons{end+1} = 'Supports complex Hermitian matrices';
        end
        
    else
        recommendation.recommended_function = 'none';
        recommendation.reasons{end+1} = 'No function meets the requirements well enough';
    end
    
    % 打印推荐
    fprintf('RECOMMENDATION:\n');
    fprintf('===============\n');
    if strcmp(recommendation.recommended_function, 'none')
        fprintf('❌ No clear recommendation - all functions have issues\n');
    else
        fprintf('✅ Recommended: %s\n', recommendation.recommended_function);
        fprintf('Score: %d/100\n', max_score);
        fprintf('Reasons:\n');
        for i = 1:length(recommendation.reasons)
            fprintf('  • %s\n', recommendation.reasons{i});
        end
    end
    
    fprintf('\nAll scores:\n');
    for i = 1:length(functions)
        if isfield(recommendation.scores, functions{i})
            fprintf('  %s: %d/100\n', functions{i}, recommendation.scores.(functions{i}));
        else
            fprintf('  %s: Failed\n', functions{i});
        end
    end
end

function create_comparison_visualization(test_results)
% 创建可视化对比
    
    % 检查是否有成功的结果
    successful_results = {};
    functions = {'original', 'improved', 'complex'};
    
    for i = 1:length(functions)
        if isfield(test_results, functions{i}) && test_results.(functions{i}).success
            successful_results{end+1} = functions{i};
        end
    end
    
    if isempty(successful_results)
        fprintf('No successful results to visualize.\n');
        return;
    end
    
    figure('Name', 'Simulation Functions Comparison', 'Position', [100, 100, 1200, 800]);
    
    n_plots = length(successful_results);
    
    for i = 1:n_plots
        func_name = successful_results{i};
        result = test_results.(func_name);
        
        % 稀疏性变化图
        subplot(2, n_plots, i);
        plot(1:length(result.edge_counts), result.edge_counts, 'o-', 'LineWidth', 2);
        title(sprintf('%s: Edge Count vs Frequency', func_name));
        xlabel('Frequency');
        ylabel('Number of Edges');
        grid on;
        
        % 条件数图
        subplot(2, n_plots, i + n_plots);
        valid_conditions = result.condition_numbers(isfinite(result.condition_numbers));
        if ~isempty(valid_conditions)
            semilogy(1:length(valid_conditions), valid_conditions, 's-', 'LineWidth', 2);
            title(sprintf('%s: Condition Numbers', func_name));
            xlabel('Frequency');
            ylabel('Condition Number (log scale)');
            grid on;
        else
            text(0.5, 0.5, 'No valid condition numbers', 'HorizontalAlignment', 'center');
            title(sprintf('%s: Condition Numbers', func_name));
        end
    end
    
    sgtitle('Simulation Functions Comparison');
end