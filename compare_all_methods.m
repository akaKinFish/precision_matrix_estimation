function results = compare_all_methods(cfg)
% COMPARE_ALL_METHODS - 统一比较你的算法与 BC-V 方法
%
% 用法:
%   results = compare_all_methods();           % 用默认配置
%   results = compare_all_methods(my_cfg);     % 自定义参数
%
% 产出:
%   results 结构体，含各方法的 Ω_src 估计、指标等

if nargin < 1
    cfg = struct('n_nodes', 10, 'n_sensors', 3, 'n_freq', 3, ...
                 'n_samples', 4096, 'leadfield_type', 'simple', ...
                 'random_seed', 42);
end

ts = datestr(now,'yyyymmdd_HHMMSS');
outdir = fullfile('results', ['exp_', ts]);
if ~exist(outdir,'dir'), mkdir(outdir); end

%% ===== 1. Simulation =====
fprintf('==[1/5] Simulation with Module 7 ==\n');
[Omega_true, Sigma_true, emp_covariance, sim] = module7_simulation_improved_complex( ...
    'n_nodes', cfg.n_nodes, 'n_sensors', cfg.n_sensors, ...
    'n_freq', cfg.n_freq, 'n_samples', cfg.n_samples, ...
    'generate_leadfield', true, 'leadfield_type', cfg.leadfield_type, ...
    'random_seed', cfg.random_seed);

L = sim.leadfield_matrix;                          % p×n
emp_cov_cell = coerce_cov_cell(emp_covariance);    % {F×1}
F = numel(emp_cov_cell);
n = size(L, 2);
T = cfg.n_samples;

save(fullfile(outdir,'sim_data.mat'), 'Omega_true', 'Sigma_true', ...
     'emp_covariance', 'sim', 'L', 'cfg');

%% ===== 2. Your algorithm (修正后的调用) =====
fprintf('==[2/5] Run YOUR algorithm (test2 adapter) ==\n');

% 配置参数（对齐 test2 脚本）
cfg_adapter = struct();
cfg_adapter.Sigma_xixi = sim.Sigma_xixi;           % 初始噪声协方差
cfg_adapter.verbose = false;
cfg_adapter.warm_method = 'ssblpp';                % 或 'eloreta'
cfg_adapter.enable_refit = true;                   % 开启 support-refit
cfg_adapter.max_iter = 30;

% ⭐ 关键：统一接口，不传 Omega_true（真实场景不可用）
[Omega_my, Dsrc_my, Gamma_my] = my_test2_adapter( ...
    emp_cov_cell, L, T, cfg_adapter);

%% ===== 3. BC-V methods =====
fprintf('==[3/5] Run BC-V methods on the SAME data ==\n');
opts_bcv = struct();
bcv = run_bcv_methods(emp_cov_cell, L, T, opts_bcv);

%% ===== 4. 构造 Omega_all（逐字段赋值，避免 struct array）=====
fprintf('==[4/5] Evaluate & summarize ==\n');

Omega_all = struct();
Omega_all.my_test2 = Omega_my;
Omega_all.bcv_higgs = bcv.higgs;
Omega_all.bcv_eloreta = bcv.eloreta_hglasso;
Omega_all.bcv_lcmv = bcv.lcmv_hglasso;
Omega_all.bcv_vareta = bcv.vareta;

methods = {'my_test2', 'bcv_higgs', 'bcv_eloreta', 'bcv_lcmv', 'bcv_vareta'};

% ⭐ 评估：传入 Dsrc 和 Gamma 用于白化域指标（如需要）
metrics = evaluate_all_methods(Omega_true, Omega_all, ...
                                struct('f_view', 1), Dsrc_my, Gamma_my);

%% ===== 5. 可视化 =====
fprintf('==[5/5] Plot partial coherence comparison ==\n');

% ⭐ 使用你已有的 plot_pcoh_panel 函数
opt_compare = struct('mode', 'match_sparsity', 'f_view', 1, ...
                     'show_pr', true, 'show_roc', true);
opt_compare.Gamma = Omega_all.my_test2;  % 可选：传白化域的 Gamma

% 调用你已有的可视化函数
plot_pcoh_panel(Omega_true, Omega_all, 1);
saveas(gcf, fullfile(outdir, 'pcoh_compare_f1.png'));

% 如果有更详细的对比绘图函数，也可以调用
if exist('gt_compare_and_plot', 'file')
    figure;
    gt_compare_and_plot(Omega_true, Omega_all.my_test2, opt_compare);
    saveas(gcf, fullfile(outdir, 'detailed_my_test2.png'));
end

%% ===== 6. 保存结果 =====
results = struct();
results.cfg = cfg;
results.outdir = outdir;
results.Omega_true = Omega_true;
results.Omega_all = Omega_all;
results.metrics = metrics;
results.Dsrc_my = Dsrc_my;         % 保存白化矩阵（可选）
results.Gamma_my = Gamma_my;       % 保存白化域 precision（可选）

save(fullfile(outdir, 'results_all.mat'), '-struct', 'results');
fprintf('\n[OK] Done. All artifacts saved under: %s\n', outdir);
end

%% ==================== Helper Functions ====================

function C = coerce_cov_cell(X)
    if iscell(X), C = X(:); return; end
    if isnumeric(X) && ndims(X)==3 && size(X,1)==size(X,2)
        F = size(X,3); C = cell(F,1);
        for f = 1:F, C{f} = (X(:,:,f) + X(:,:,f)') / 2; end
        return;
    end
    if isnumeric(X) && ismatrix(X) && size(X,1)==size(X,2)
        C = {(X + X') / 2}; return;
    end
    error('emp_covariance: unsupported format');
end