function test_results = test_module7_simulation_bayesian()
% TEST_MODULE7_SIMULATION_BAYESIAN
% 单元测试：调用 "module7_simulation_bayesian" 生成 Bayes 真值与样本，
% 检查 SPD、样本协方差与真协方差的相对误差、超参数记录，并可选调用
% 你的 my_test2_adapter 验证整条流水线是否跑通（若该函数在路径上）。

fprintf('============================================\n');
fprintf('Testing Bayes simulator: module7_simulation_bayesian\n');
fprintf('============================================\n');
rng(2025);

% ---------- 配置 ----------
cfg = struct();
cfg.p = 10; cfg.m = 3; cfg.F = 3; cfg.T = 200;
cfg.edge_density = 0.15;
cfg.modeN = 'diff';            % 'diff' 或 'kernel'（需要 cfg.K）
cfg.lambda1_star = 0.35;       % 固定真超参（推荐）
cfg.lambda2_star = 1.0;
cfg.sigma_xi2 = 0.2;
cfg.delta_eps = [0.1, 1e-6];
cfg.complex_samples = true;    % 复频域样本

% ---------- 生成 ----------
[Omega_true, Sigma_true, emp_cov, L, T, bayes_truth, sim] = module7_simulation_bayesian(cfg);

% ---------- 统计/校验 ----------
ok_spd = true;
min_eigs = zeros(sim.F,1);
avg_rel_cov_err = 0;
for f=1:sim.F
    lam = eig( (Omega_true(:,:,f)+Omega_true(:,:,f)')/2 );
    min_eigs(f) = min(real(lam));
    ok_spd = ok_spd && all(real(lam) > 0);
    relE = norm(emp_cov{f} - Sigma_true{f}, 'fro')/max(norm(Sigma_true{f},'fro'), eps);
    avg_rel_cov_err = avg_rel_cov_err + relE/sim.F;
end

fprintf('\n=== Ground-truth hyperparameters ===\n');
fprintf('lambda1_true = %.4f\n', bayes_truth.lambda1_true);
fprintf('lambda2_true = %.4f\n', bayes_truth.lambda2_true);
fprintf('SPD check (all freqs): %d (min eig range [%.4g, %.4g])\n', ok_spd, min(min_eigs), max(min_eigs));
fprintf('Average rel. error of sample cov vs true cov: %.3f\n', avg_rel_cov_err);

% 简单成功判据（可按需放宽/收紧）
success = ok_spd && (avg_rel_cov_err < 0.25);
fprintf('Success flag (simulator self-check): %d\n', success);
% 统计每个频点的活跃边数（不含对角）
act_per_f = zeros(sim.F,1);
for f=1:sim.F
    Of = Omega_true(:,:,f);
    supf = triu(abs(Of)>1e-12, 1);   % off-diagonal support
    act_per_f(f) = nnz(supf);
end
fprintf('Active off-diagonal edges per freq: [');
fprintf('%d ', act_per_f);
fprintf(']\n');
S = bayes_truth.support_mask;   % d × F
figure('Name','Support across frequency','Color','w');
imagesc(S'); colorbar; xlabel('edge index'); ylabel('frequency'); title('support (1=active)');
% ---------- 可视化 ----------
visualize_sim_results(Omega_true, Sigma_true, emp_cov, L, bayes_truth);

% ---------- （可选）调用你的算法适配器，串通后续流程 ----------
% if exist('my_test2_adapter','file') == 2
%     try
%         fprintf('\nRunning my_test2_adapter on simulated data...\n');
%         % 你的适配器接口通常为：my_test2_adapter(emp_cov, L, T, cfg, Omega_true)
%         cfg_for_adapter = struct();  % 需要什么就填什么；这里给个空壳即可
%         [Omega_est, Dsrc, Gamma_tilde, outs] = my_test2_adapter(emp_cov, L, T, cfg_for_adapter, Omega_true); %#ok<ASGLU>
%         fprintf('my_test2_adapter finished. (Not scoring here; use your existing evaluators.)\n');
%     catch ME
%     end
% else
%     fprintf('\n(my_test2_adapter not found on path; skipped downstream run.)\n');
% end

% ---------- 汇总 ----------
test_results = struct();
test_results.success = success;
test_results.avg_rel_cov_err = avg_rel_cov_err;
test_results.ok_spd = ok_spd;
test_results.min_eigs = min_eigs;
test_results.lambda1_true = bayes_truth.lambda1_true;
test_results.lambda2_true = bayes_truth.lambda2_true;
test_results.sim = sim;

fprintf('\nSummary: success=%d, avg_rel_cov_err=%.3f\n', success, avg_rel_cov_err);
end

% ================= helpers for test =================
function visualize_sim_results(Omega_true, Sigma_true, emp_cov, L, bayes_truth)
[p,~,F] = size(Omega_true);

% 1) Ω_true：每个频率的 |Ω| 热图
fig1 = figure('Name','Omega_true (per frequency)','Color','w');
tiledlayout(fig1,1,F,'Padding','compact','TileSpacing','compact');
for f=1:F
    nexttile;
    imagesc(abs(Omega_true(:,:,f)));
    axis image; colorbar;
    title(sprintf('|\\Omega|, f=%d', f));
    xlabel('node'); ylabel('node');
end

% 2) Σ_true vs emp_cov：显示第1个频率
fig2 = figure('Name','Sigma_true vs emp_cov (f=1)','Color','w');
tiledlayout(fig2,1,3,'Padding','compact','TileSpacing','compact');
f = 1;
nexttile; imagesc(abs(Sigma_true{f})); axis image; colorbar; title('|Σ\_true|');
nexttile; imagesc(abs(emp_cov{f}));    axis image; colorbar; title('|Σ\_emp|');
nexttile; imagesc(abs(emp_cov{f}-Sigma_true{f})); axis image; colorbar; title('|ΔΣ|');

% 3) L 的简单诊断
fig3 = figure('Name','Leadfield diagnostics','Color','w');
subplot(1,2,1); imagesc(L); colorbar; title('L entries'); xlabel('source'); ylabel('sensor');
subplot(1,2,2); bar([norm(L,'fro'), cond(L)]); set(gca,'XTickLabel',{'||L||_F','cond(L)'}); grid on;

% 4) 打印超参
% 计算“至少一个频点激活”的边数（不含对角）
if isfield(bayes_truth, 'support_idx')
    n_unique = numel(bayes_truth.support_idx);
elseif isfield(bayes_truth, 'support_mask')
    n_unique = nnz(any(bayes_truth.support_mask, 2));
else
    % 兜底：从 Omega_true 现算
    [p,~,F] = size(Omega_true);
    sup_any = false(p*(p-1)/2,1);
    [iu,ju] = nchoosek(1:p,2); iu=iu(:); ju=ju(:);
    for f=1:F
        Of = Omega_true(:,:,f);
        sup_f = abs(triu(Of,1)) > 1e-12;
        sup_any = sup_any | sup_f(sub2ind([p,p],iu,ju));
    end
    n_unique = nnz(sup_any);
end
fprintf('\n[viz] lambda1_true=%.4f, lambda2_true=%.4f, #unique active edges=%d\n', ...
    bayes_truth.lambda1_true, bayes_truth.lambda2_true, n_unique);
end
