%% 修正的复数Hermitian演示脚本
clear; close all; clc;

% 添加路径
addpath('utils');

fprintf('=========================================\n');
fprintf('Fixed Complex Hermitian Simulation Demo\n');
fprintf('=========================================\n\n');

%% 使用调试版本生成复数矩阵
fprintf('Generating complex Hermitian matrices with debug version...\n');

[prec_complex, cov_complex, emp_complex, params_complex] = module7_simulation_improved_complex(...
    'n_nodes', 12, ...
    'n_freq', 25, ...
    'n_samples', 150, ...
    'graph_type', 'random', ...
    'edge_density', 0.3, ...  % 增加边密度确保有足够的边
    'complex_strength', 1.2, ... % 增强复数特性
    'random_seed', 123);

%% 示例1：可视化复数矩阵属性
fprintf('\nExample: Complex Matrix Visualization\n');
fprintf('------------------------------------\n');

f_idx = 10;
P_complex = prec_complex{f_idx};

figure('Name', 'Fixed Complex Hermitian Properties', 'Position', [100, 100, 1200, 800]);

% 幅度
subplot(2, 3, 1);
imagesc(abs(P_complex));
colorbar;
title(sprintf('Amplitude |Ω(%d)|', f_idx));
axis square;

% 相位
subplot(2, 3, 2);
phase_matrix = angle(P_complex);
imagesc(phase_matrix);
colorbar;
colormap(gca, hsv);
title(sprintf('Phase ∠Ω(%d)', f_idx));
axis square;
caxis([-pi, pi]);  % 确保相位范围正确

% 实部
subplot(2, 3, 3);
imagesc(real(P_complex));
colorbar;
title(sprintf('Real(Ω(%d))', f_idx));
axis square;

% 虚部
subplot(2, 3, 4);
imagesc(imag(P_complex));
colorbar;
title(sprintf('Imag(Ω(%d))', f_idx));
axis square;

% 特征值谱
subplot(2, 3, 5);
eigvals = eig(P_complex);
plot(real(eigvals), imag(eigvals), 'bo', 'MarkerSize', 8, 'LineWidth', 2);
xlabel('Real Part');
ylabel('Imaginary Part');
title('Eigenvalue Spectrum');
grid on;

% 稀疏模式
subplot(2, 3, 6);
spy(abs(P_complex) > 1e-6);
title('Sparsity Pattern');

%% 示例2：修正的相位演化分析
fprintf('\nExample: Fixed Phase Evolution Analysis\n');
fprintf('--------------------------------------\n');

% 智能选择边：从实际存在的边中选择
actual_edges = [];
threshold = 1e-6;

% 找到实际存在的边
for i = 1:params_complex.n_nodes
    for j = i+1:params_complex.n_nodes
        % 检查这条边在任何频率是否存在
        edge_exists = false;
        for f = 1:params_complex.n_freq
            if abs(prec_complex{f}(i, j)) > threshold
                edge_exists = true;
                break;
            end
        end
        if edge_exists
            actual_edges = [actual_edges; i, j];
        end
    end
end

fprintf('Found %d actual edges\n', size(actual_edges, 1));

if size(actual_edges, 1) >= 3
    % 选择前3条边进行分析
    selected_edges = actual_edges(1:3, :);
else
    % 如果实际边不够，使用基础边列表
    if isfield(params_complex, 'base_edge_list') && size(params_complex.base_edge_list, 1) >= 3
        selected_edges = params_complex.base_edge_list(1:3, :);
    else
        % 最后的备选方案
        selected_edges = [2, 1; 3, 1; 4, 2];
    end
end

figure('Name', 'Fixed Phase Evolution Analysis', 'Position', [100, 100, 1200, 400]);

for e = 1:size(selected_edges, 1)
    i = selected_edges(e, 1);
    j = selected_edges(e, 2);
    
    % 提取值
    amplitudes = zeros(params_complex.n_freq, 1);
    phases = zeros(params_complex.n_freq, 1);
    
    for f = 1:params_complex.n_freq
        val = prec_complex{f}(i, j);
        amplitudes(f) = abs(val);
        phases(f) = angle(val);
    end
    
    % 检查是否有复数值
    has_complex = any(abs(imag([prec_complex{:}])) > 1e-10, 'all');
    fprintf('Edge (%d,%d): Max amplitude = %.4f, Phase range = [%.4f, %.4f]\n', ...
            i, j, max(amplitudes), min(phases), max(phases));
    
    % 幅度演化
    subplot(2, 3, e);
    plot(1:params_complex.n_freq, amplitudes, 'b-', 'LineWidth', 2);
    xlabel('Frequency');
    ylabel('Amplitude');
    title(sprintf('Edge (%d,%d) Amplitude', i, j));
    grid on;
    
    % 相位演化
    subplot(2, 3, e + 3);
    plot(1:params_complex.n_freq, phases, 'r-', 'LineWidth', 2);
    xlabel('Frequency');
    ylabel('Phase (radians)');
    title(sprintf('Edge (%d,%d) Phase', i, j));
    grid on;
    ylim([-pi, pi]);
    
    % 添加相位变化统计
    phase_variation = std(phases);
    text(0.05, 0.95, sprintf('σ=%.3f', phase_variation), ...
         'Units', 'normalized', 'BackgroundColor', 'white');
end

%% 示例3：矩阵复数特性验证
fprintf('\nExample: Matrix Complex Properties Verification\n');
fprintf('----------------------------------------------\n');

% 分析所有频率的复数特性
complex_fractions = zeros(params_complex.n_freq, 1);
phase_variations = zeros(params_complex.n_freq, 1);
max_imaginary = zeros(params_complex.n_freq, 1);

for f = 1:params_complex.n_freq
    P = prec_complex{f};
    n = size(P, 1);
    
    % 统计非对角元素
    off_diag_elements = [];
    for i = 1:n
        for j = 1:n
            if i ~= j
                off_diag_elements = [off_diag_elements; P(i, j)];
            end
        end
    end
    
    % 复数分数
    complex_count = sum(abs(imag(off_diag_elements)) > 1e-10);
    complex_fractions(f) = complex_count / length(off_diag_elements);
    
    % 相位变化
    phases = angle(off_diag_elements);
    phases = phases(abs(imag(off_diag_elements)) > 1e-10);
    if ~isempty(phases)
        phase_variations(f) = std(phases);
    end
    
    % 最大虚部
    max_imaginary(f) = max(abs(imag(off_diag_elements)));
end

figure('Name', 'Complex Properties Analysis', 'Position', [100, 100, 1200, 400]);

subplot(1, 3, 1);
plot(1:params_complex.n_freq, complex_fractions, 'g-', 'LineWidth', 2);
xlabel('Frequency');
ylabel('Complex Fraction');
title('Fraction of Complex Elements');
grid on;

subplot(1, 3, 2);
plot(1:params_complex.n_freq, phase_variations, 'm-', 'LineWidth', 2);
xlabel('Frequency');
ylabel('Phase Variation (std)');
title('Phase Diversity');
grid on;

subplot(1, 3, 3);
plot(1:params_complex.n_freq, max_imaginary, 'c-', 'LineWidth', 2);
xlabel('Frequency');
ylabel('Max |Imaginary|');
title('Maximum Imaginary Component');
grid on;

%% 总结报告
fprintf('\n=== Complex Simulation Summary ===\n');
fprintf('Overall complex fraction: %.2f\n', mean(complex_fractions));
fprintf('Average phase variation: %.4f radians\n', mean(phase_variations));
fprintf('Maximum imaginary component: %.4f\n', max(max_imaginary));

if mean(complex_fractions) > 0.1
    fprintf('✓ Successfully generated complex matrices!\n');
else
    fprintf('✗ Matrices are still predominantly real.\n');
    fprintf('  Consider increasing complex_strength parameter.\n');
end

fprintf('================================\n');