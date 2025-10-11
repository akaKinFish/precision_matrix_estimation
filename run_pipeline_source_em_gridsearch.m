function run_pipeline_source_em_gridsearch()
% RUN_PIPELINE_SOURCE_EM_GRIDSEARCH
% - 手动网格搜索 (lambda1, lambda2, lambda3_ratio)
% - 固定 λ2（禁用内层退火），alpha0 由 Lipschitz 候选给出（更稳）
% - 与 module6 的自动建议对比并作图
% 
% 主要改进：
% 1. 频率内核参数自适应调整
% 2. EM收敛判据收紧
% 3. 参数冲突清理
% 4. 评分函数阈值调整

%% 0) Simulation (与主脚本一致)
n  = 10;  p  = 3;  F  = 3;  T  = 4096;
[Omega_true, Sigma_true, emp_covariance, sim] = module7_simulation_improved_complex( ...
    'n_nodes', n, 'n_sensors', p, 'n_freq', F, 'n_samples', T, ...
    'generate_leadfield', true, 'leadfield_type', 'simple', 'random_seed', 42);

L          = sim.leadfield_matrix;
Sigma_xixi = sim.Sigma_xixi;

% 传感器互谱集合 {F×1}
emp_cov_cell = coerce_cov_cell(emp_covariance);
F            = numel(emp_cov_cell);

% eLORETA warm start（一次性，供每个 grid 组合复用，公平）
gamma_grid = logspace(-4, 1, 30);
opts_el    = struct('maxit', 50, 'tol', 1e-6, 'verbose', false);
warm       = eloreta_warmstart_from_covs(emp_cov_cell, L, gamma_grid, opts_el);

%% 1) 固定资源：K / W / Lsp（改进：自适应频率内核）
% 改进1：根据频率数量自适应调整sigma参数
sigma_adaptive = max(1.0, F/3.0);  % F=3时sigma=1.0，更大F时适当增加
K = make_frequency_kernel(F, sigma_adaptive);   % 使用自适应sigma
fprintf('使用自适应频率内核参数: sigma=%.2f (F=%d)\n', sigma_adaptive, F);

W  = make_uniform_weight(n);

% 图拉普拉斯矩阵处理
if isfield(sim,'source_graph') && isfield(sim.source_graph,'L')
    Lsp_raw = sim.source_graph.L;
else
    Lsp_raw = laplacian_placeholder(n);
end
[Lsp, ~] = normalize_graph_laplacian(Lsp_raw, 'spectral');

% K的谱归一化，控制 ||K||_2 ≤ 1 （更稳）
K = (K + K')/2;
try
    smax = eigs(K,1,'lm');
catch
    smax = max(abs(eig(full(K))));
end
if isfinite(smax) && real(smax) > 1
    K = K / real(smax);
    fprintf('频率内核谱归一化: ||K||_2 从 %.3f 缩放到 %.3f\n', real(smax), 1.0);
end

%% 2) EM 外层固定设置（改进：收紧收敛判据）
em = struct();
em.max_em_iter      = 20;    % 改进2：增加最大迭代次数 (10->20)
em.tol_Omega        = 5e-4;  % 改进2：收紧Omega收敛容差 (1e-3->5e-4)
em.tol_S            = 2e-4;  % 改进2：收紧S收敛容差 (5e-4->2e-4)
em.min_hold_iters   = 3;     % 改进2：增加稳定轮次 (2->3)
em.update_act_every = 2;

fprintf('EM设置: max_iter=%d, tol_Omega=%.1e, tol_S=%.1e, min_hold=%d\n', ...
    em.max_em_iter, em.tol_Omega, em.tol_S, em.min_hold_iters);

% 内层 prox/line-search 控制（基础值；alpha0 会按每次 E-step 的 L 候选更新）
inner = struct( ...
    'alpha0',  0.1, ...
    'alpha_max', 0.6, ...
    'alpha_up', 1.2, ...
    'alpha_down', 0.6, ...
    'alpha_grow_patience', 1, ...
    'max_iter', 30, ...
    'obj_improve_tol', 5e-6, ...
    'weight_mode','hadamard', ...
    'use_graph_laplacian', true, ...
    'diag', struct('enable', false));

%% 3) 定义网格（可按需调整）
lambda1_grid       = logspace(-3, -1, 9);   % 1e-3 ... 1e-1
lambda2_grid       = logspace(-4, -1.5, 8); % 1e-4 ... ~3e-2
lambda3_ratio_grid = [0, 0.1, 0.3];         % λ3 = ratio * λ1

fprintf('Grid sizes: |λ1|=%d, |λ2|=%d, |ratio|=%d ⇒ total=%d\n', ...
    numel(lambda1_grid), numel(lambda2_grid), numel(lambda3_ratio_grid), ...
    numel(lambda1_grid)*numel(lambda2_grid)*numel(lambda3_ratio_grid));

%% 4) 先跑一遍 module6（仅用于"对比点"）
prior_cell   = repmat({eye(n)}, F, 1); 
F            = numel(emp_cov_cell);
estep_in0 = struct();                                  % scalar struct
estep_in0.leadfield_matrix         = L;                % p×n
estep_in0.empirical_covariances    = emp_cov_cell;     % {F×1}, p×p
estep_in0.source_prior_covariances = prior_cell;       % {F×1}, n×n
estep_in0.noise_covariance         = Sigma_xixi;       % p×p
estep_in0.frequencies              = (1:F); 
estep_out0 = module2_estep(estep_in0, struct('ensure_hermitian',true,'ensure_real_diag',true,'ensure_psd',true,'psd_tol',1e-10,'diag_loading',1e-10));
pre0 = module1_preproc_from_covset(estep_out0.source_second_moments, struct('smoothing_method','moving_average','loading_factor',1e-6,'min_power',1e-10,'verbose',false));
input_data_m6 = struct();
input_data_m6.whitened_covariances = pre0.Sigma_tilde;
input_data_m6.kernel_matrix = K;
input_data_m6.weight_matrix = W;
input_data_m6.active_set_mask = [];
hp6 = module6_hyperparameter_config(input_data_m6, struct('use_gershgorin', true));
lambda1_m6 = hp6.lambda1;
lambda2_m6 = hp6.lambda2_suggested;
lambda3_ratio_m6 = 0.3;  % 与主脚本一致
fprintf('[module6] λ1=%.3g, λ2=%.3g, ratio(λ3/λ1)=%.3g\n', lambda1_m6, lambda2_m6, lambda3_ratio_m6);

%% 5) 评分设置（改进：调整阈值）
score_mode = 'AUPR';   % 'AUPR' | 'AUROC' | 'combo'
use_abs    = true;
edge_threshold = 1e-8; % 改进4：真边判定阈值从1e-12调整到1e-8
fprintf('评分设置: mode=%s, use_abs=%s, edge_threshold=%.1e\n', ...
    score_mode, mat2str(use_abs), edge_threshold);

%% 6) Grid search 主循环
S = struct();
best = struct('score',-Inf,'l1',NaN,'l2',NaN,'r3',NaN,'Omega',[]);
S.scores = zeros(numel(lambda1_grid), numel(lambda2_grid), numel(lambda3_ratio_grid));

total_combinations = numel(lambda1_grid) * numel(lambda2_grid) * numel(lambda3_ratio_grid);
combination_count = 0;

fprintf('\n开始网格搜索...\n');
tic;

for k3 = 1:numel(lambda3_ratio_grid)
    r3 = lambda3_ratio_grid(k3);
    fprintf('\n--- 处理 λ3_ratio = %.2g ---\n', r3);
    
    for i1 = 1:numel(lambda1_grid)
        l1 = lambda1_grid(i1);
        for j2 = 1:numel(lambda2_grid)
            l2 = lambda2_grid(j2);
            combination_count = combination_count + 1;

            % 跑一遍 EM（固定 λ）
            out = em_run_fixed_lambdas( ...
                emp_cov_cell, L, Sigma_xixi, warm, ...
                K, W, Lsp, ...
                l1, l2, r3, ...
                em, inner);

            % 打分（使用改进的阈值）
            sc = score_against_gt(out.Omega_src, Omega_true, use_abs, edge_threshold);
            switch upper(score_mode)
                case 'AUPR',  score = sc.aupr;
                case 'AUROC', score = sc.auroc_pos;
                case 'COMBO', score = 0.7*sc.aupr + 0.3*sc.auroc_pos;
                otherwise,    score = sc.aupr;
            end

            S.scores(i1,j2,k3) = score;

            if score > best.score
                best.score = score; best.l1=l1; best.l2=l2; best.r3=r3; best.Omega = out.Omega_src;
            end

            % 进度显示
            if mod(combination_count, 20) == 1 || combination_count <= 10
                elapsed = toc;
                eta = elapsed * (total_combinations - combination_count) / combination_count;
                fprintf('[%3d/%3d] ratio=%.2g | λ1=%.3g λ2=%.3g => %s=%.4f (ΔΩ=%.3g, ΔS=%.3g) [ETA: %.1fs]\n', ...
                        combination_count, total_combinations, r3, l1, l2, score_mode, score, ...
                        out.delta_Omega_end, out.delta_S_end, eta);
            end
        end
    end
end

total_time = toc;
fprintf('\n网格搜索完成! 总时间: %.1f秒 (平均 %.2f秒/组合)\n', ...
    total_time, total_time/total_combinations);

%% 7) 作图：每个 ratio 画一张 (λ1, λ2) 热力图，并标注 module6 点与最优点
for k3 = 1:numel(lambda3_ratio_grid)
    r3 = lambda3_ratio_grid(k3);
    figure('Name',sprintf('Grid heatmap ratio=%.2g (%s)', r3, score_mode));
    imagesc(log10(lambda2_grid), log10(lambda1_grid), S.scores(:,:,k3));
    set(gca,'YDir','normal'); 
    colorbar;
    xlabel('log_{10}(\lambda_2)'); 
    ylabel('log_{10}(\lambda_1)');
    title(sprintf('%s heatmap (\\lambda_3=%.2g\\lambda_1, \\sigma_{freq}=%.2f)', ...
        score_mode, r3, sigma_adaptive));

    hold on;
    % 标注module6建议点
    [~,i1m] = min(abs(lambda1_grid - lambda1_m6));
    [~,j2m] = min(abs(lambda2_grid - lambda2_m6));
    plot(log10(lambda2_grid(j2m)), log10(lambda1_grid(i1m)), 'wx', 'markersize', 12, 'linewidth', 2);
    text(log10(lambda2_grid(j2m)), log10(lambda1_grid(i1m)), '  module6', 'color','w', 'fontweight','bold');

    % 标注最优点
    [mx,idx] = max(S.scores(:,:,k3), [], 'all','linear'); %#ok<ASGLU>
    [i1b,j2b] = ind2sub(size(S.scores(:,:,k3)), idx);
    plot(log10(lambda2_grid(j2b)), log10(lambda1_grid(i1b)), 'wo', 'markersize', 10, 'linewidth', 2);
    text(log10(lambda2_grid(j2b)), log10(lambda1_grid(i1b)), '  best', 'color','w', 'fontweight','bold');
    hold off;
    
    % 添加分数范围信息
    score_range = [min(S.scores(:,:,k3),[],'all'), max(S.scores(:,:,k3),[],'all')];
    subtitle(sprintf('Score range: [%.4f, %.4f]', score_range(1), score_range(2)));
end

%% 8) 汇总
fprintf('\n===== BEST (by %s) =====\n', score_mode);
fprintf('ratio=%.3g | λ1=%.3g λ2=%.3g | score=%.4f\n', best.r3, best.l1, best.l2, best.score);
fprintf('module6: λ1=%.3g λ2=%.3g ratio=%.3g\n', lambda1_m6, lambda2_m6, lambda3_ratio_m6);

% 计算改进幅度
module6_score = NaN;
try
    [~,i1m] = min(abs(lambda1_grid - lambda1_m6));
    [~,j2m] = min(abs(lambda2_grid - lambda2_m6));
    [~,k3m] = min(abs(lambda3_ratio_grid - lambda3_ratio_m6));
    module6_score = S.scores(i1m, j2m, k3m);
    improvement = (best.score - module6_score) / max(module6_score, 1e-6) * 100;
    fprintf('相对于module6的改进: %.2f%%\n', improvement);
catch
    fprintf('无法计算module6性能（参数超出网格范围）\n');
end

opt_compare = struct('mode','match_sparsity','f_view',1,'show_pr',true,'show_roc',true);
opt_compare.Gamma = best.Omega;
gt_compare_and_plot(Omega_true, best.Omega, opt_compare);

end % run_pipeline_source_em_gridsearch


%% ===== 子程序：固定 λ 的 EM（改进版） =====
function out = em_run_fixed_lambdas(Svv_cell, L, Sigma_xixi, warm, K, W, Lsp, lambda1, lambda2, lambda3_ratio, em, inner)
    F  = numel(Svv_cell);
    n  = size(L,2);

    estep_in = struct();
    estep_in.leadfield_matrix         = L;
    estep_in.empirical_covariances    = Svv_cell;
    estep_in.source_prior_covariances = repmat({eye(n)}, F, 1);
    estep_in.noise_covariance         = Sigma_xixi;
    estep_in.frequencies              = 1:F;

    Omega_prev = warm.Omega_init;  % 每个 grid 组合均以相同 warm start 起步
    Sjj_prev   = warm.Sjj_e;

    A_masks = [];
    hold_counter = 0;

    for t = 1:em.max_em_iter
        % E-step
        estep_out = module2_estep(estep_in, struct('ensure_hermitian',true,'ensure_real_diag',true,'ensure_psd',true,'psd_tol',1e-10,'diag_loading',1e-10));
        Sjj_hat   = estep_out.source_second_moments;

        % 预处理白化
        pre = module1_preproc_from_covset(Sjj_hat, struct( ...
            'smoothing_method','moving_average','loading_factor',1e-6,'min_power',1e-10,'verbose',false));
        D_src     = pre.D;
        Sjj_tilde = pre.Sigma_tilde;

        % Active set
        if t<=2 || isempty(A_masks) || mod(t, em.update_act_every)==0
            input_data_m3 = struct();
            input_data_m3.whitened_covariances = Sjj_tilde;
            input_data_m3.frequencies = 1:F;
            act = module3_active_set(input_data_m3, struct('proxy_method','correlation','quantile_level',0.25,'force_diagonal_active',true,'verbose',false));
            A_masks = arrayfun(@(f) logical(act.combined_active_mask(:,:,f)), 1:F, 'uni', 0);
        end

        % α_rec
        [L_byS, L_byG] = estimate_L_candidates_local(Omega_prev, D_src, Sjj_tilde);
        alpha_rec = 0.8 / max(min([L_byS, L_byG]), 1e-12);
        alpha0_used = min(alpha_rec, inner.alpha_max);

        % M-step（白化域）
        Gamma_init = transport_init_local(Omega_prev, D_src, Sjj_tilde);

        input_data_m5 = struct();
        input_data_m5.whitened_covariances = Sjj_tilde;
        input_data_m5.initial_precision = Gamma_init;
        input_data_m5.smoothing_kernel = K;
        input_data_m5.weight_matrix = W;
        input_data_m5.active_set_mask = {A_masks{:}};
        input_data_m5.whitening_matrices = D_src;
        
        % 改进3：清理参数冲突，明确禁用λ2退火
        params5 = struct( ...
            'lambda1', lambda1, ...
            'lambda2', lambda2, ...
            'lambda2_suggested', lambda2, ...
            'alpha0',  alpha0_used, ...
            'max_iter', inner.max_iter, 'verbose', false, ...
            'active_set_update_freq', 10, ...
            'alpha_max', inner.alpha_max, 'alpha_up', inner.alpha_up, ...
            'alpha_down', inner.alpha_down, 'alpha_grow_patience', inner.alpha_grow_patience, ...
            'obj_improve_tol', inner.obj_improve_tol, ...
            'weight_mode', inner.weight_mode, 'use_graph_laplacian', inner.use_graph_laplacian, ...
            'diag', inner.diag);

        % 改进3：明确禁用λ2退火机制
        params5.backtrack_beta              = 0.5;
        params5.armijo_c1                   = 1e-4;
        params5.max_backtrack_per_iter      = 20;
        params5.backtrack_patience          = Inf;       % 不触发降 λ2
        params5.lambda2_decay_factor        = 1.0;       % 不衰减
        params5.lambda2_min                 = lambda2;   % 保持原值
        params5.alpha_min                   = 1e-4;
        params5.penalize_diagonal           = false;
        params5.enable_lambda2_annealing    = false;     % 明确禁用标志

        % 空间平滑
        lambda3 = lambda3_ratio * lambda1;
        if lambda3 > 0
            params5.lambda3                    = lambda3;
            params5.spatial_graph_matrix       = Lsp;
            params5.spatial_graph_is_laplacian = true;
            params5.spatial_weight_mode        = 'node';
        end

        [Gamma_tilde_star, prox_res] = module5_proximal(input_data_m5, params5); %#ok<NASGU>
        
        input_data_m8 = struct();
        input_data_m8.whitened_precision_matrices = Gamma_tilde_star;
        input_data_m8.whitening_matrices = D_src;
        % recolor
        recol = module8_recoloring(input_data_m8, struct());
        
        Omega_src = recol.recolored_precision_matrices;

        % 收敛统计
        [delta_Omega, delta_S] = compute_deltas_local(Omega_src, Omega_prev, Sjj_hat, Sjj_prev);

        % 停止判据
        if delta_Omega < em.tol_Omega && delta_S < em.tol_S
            hold_counter = hold_counter + 1;
        else
            hold_counter = 0;
        end
        if hold_counter >= em.min_hold_iters
            break;
        end

        % 回写先验 + 缓存
        estep_in.source_prior_covariances = invert_and_fix_local(Omega_src, 1e-10);
        Omega_prev = Omega_src; Sjj_prev = Sjj_hat;
    end

    out = struct('Omega_src',{Omega_src}, 'delta_Omega_end',delta_Omega, 'delta_S_end',delta_S);
end


%% ===== 评分：与真值对比（AUPR/AUROC）- 改进版 =====
function sc = score_against_gt(Om_est_cell, Om_true_cell, use_abs, edge_threshold)
    if nargin<3, use_abs=true; end
    if nargin<4, edge_threshold=1e-8; end  % 改进4：默认阈值调整
    
    F = numel(Om_est_cell);
    n = size(Om_est_cell{1},1);
    idxUT = triu(true(n),1);
    s_all=[]; y_all=[];
    
    for f=1:F
        E = Om_est_cell{f}; T = Om_true_cell{min(f,numel(Om_true_cell))};
        s = E(idxUT); t = T(idxUT);
        if use_abs, s=abs(s); t=abs(t); end
        y = t > edge_threshold;  % 改进4：使用调整后的阈值
        s_all = [s_all; s];
        y_all = [y_all; y];
    end
    
    % 添加边缘情况检查
    if all(y_all == y_all(1))  % 所有标签相同
        warning('所有边的标签相同，AUROC/AUPR可能不可靠');
        sc = struct('aupr', 0, 'auroc_pos', 0.5);
        return;
    end
    
    [~,~,~,auroc_pos] = perfcurve(y_all, s_all, true);
    [rec,prec,~,aupr] = perfcurve(y_all, s_all, true, 'xCrit','reca','yCrit','prec'); %#ok<ASGLU>
    sc = struct('aupr',aupr,'auroc_pos',auroc_pos);
end


%% ====== 轻量本地工具（保持不变） ======
function Gamma_init = transport_init_local(Omega_prev, D_src, Sjj_tilde)
% 将上一轮源域精度搬运到当前白化域：Γ_init = D^{-1} Ω_prev D^{-1}
% 若上一轮为空，则用 “稳健逆 S̃” 作为启动（Hermitian → 地板 → 特征逆 → 对称化）
    F = numel(D_src); 
    Gamma_init = cell(F,1);

    if isempty(Omega_prev)
        for f = 1:F
            St = (Sjj_tilde{f} + Sjj_tilde{f}')/2;   % Hermitian
            [U, D] = eig(full(St), 'vector');
            d = real(D);
            d = max(d, 1e-10);                       % 地板
            G = U * diag(1./d) * U';                 % 稳健逆
            G = (G + G')/2;                          % 数值对称化
            Gamma_init{f} = G;                       % ←← 必须是 G；启动阶段**不要**访问 Omega_prev
        end
        return;
    end

    for f = 1:F
        D = D_src{f};
        Gamma_init{f} = (D \ Omega_prev{f}) / D;     % D^{-1} * Ω_prev * D^{-1}
        Gamma_init{f} = (Gamma_init{f} + Gamma_init{f}')/2;
    end
end


function [L_byS, L_byG] = estimate_L_candidates_local(Omega_prev, D_src, Sjj_tilde)
    F = numel(Sjj_tilde); Ls=0; Lg=0;
    % 用 S̃ 的谱范数估计一个上界
    for f=1:F
        St = (Sjj_tilde{f} + Sjj_tilde{f}')/2;
        s  = svds(St,1);
        Ls = max(Ls, s^2);
    end
    % 用 Γ_init 的最小特征值估计另一个上界
    Gam = transport_init_local(Omega_prev, D_src, Sjj_tilde);
    for f=1:F
        G = (Gam{f} + Gam{f}')/2;
        try
            lam_min = min(real(eig(G)));
        catch
            lam_min = NaN;
        end
        if ~isfinite(lam_min) || lam_min <= 0
            invnorm = svds(pinv(G), 1);
        else
            invnorm = 1/lam_min;
        end
        Lg = max(Lg, invnorm^2);
    end
    L_byS=Ls; L_byG=Lg;
end

function [dOmega, dS] = compute_deltas_local(Omega, Omega_prev, Sjj, Sjj_prev)
    if isempty(Omega_prev)
        dOmega = inf;
    else
        num=0; den=0; 
        for f=1:numel(Omega)
            num=num+norm(Omega{f}-Omega_prev{f},'fro'); 
            den=den+norm(Omega_prev{f},'fro'); 
        end
        dOmega = num / max(1,den);
    end
    if isempty(Sjj_prev)
        dS = inf;
    else
        num=0; den=0; 
        for f=1:numel(Sjj)
            num=num+norm(Sjj{f}-Sjj_prev{f},'fro'); 
            den=den+norm(Sjj_prev{f},'fro'); 
        end
        dS = num / max(1,den);
    end
end

function Sigma_prior = invert_and_fix_local(Omega_cell, eps_ld)
    if nargin<2, eps_ld=1e-10; end
    F = numel(Omega_cell); n = size(Omega_cell{1},1);
    Sigma_prior = cell(F,1);
    for f=1:F
        Om = (Omega_cell{f}+Omega_cell{f}')/2; Om(~isfinite(Om))=0;
        d  = max(real(diag(Om)), eps_ld); Om(1:n+1:end) = d;
        [U,S] = eig(full(Om),'vector'); S = max(real(S), 2*eps_ld);
        Sigma = U*diag(1./S)*U'; Sigma = (Sigma+Sigma')/2; Sigma(1:n+1:end)=real(diag(Sigma));
        Sigma_prior{f}=Sigma;
    end
end

%% ===== Local helpers =====
function C = coerce_cov_cell(X, F_hint)
% Turn X into {F×1} cell of square matrices.
    if isa(X,'cell'), C = X(:); return; end
    if isnumeric(X) && ndims(X)==3 && size(X,1)==size(X,2)
        F = size(X,3); C = cell(F,1);
        for f = 1:F, C{f} = X(:,:,f); end, return;
    end
    if isnumeric(X) && ismatrix(X) && size(X,1)==size(X,2)
        if nargin>=2 && ~isempty(F_hint), C = repmat({X}, F_hint, 1);
        else, C = {X}; end, return;
    end
    error('coerce_cov_cell:unsupported','Expect cell{F,1} | p×p×F | single p×p.');
end

function K = make_frequency_kernel(F, sigma)
% Gaussian kernel along frequency index (1..F), row-normalized
    if nargin < 2, sigma = 3.0; end
    [I,J] = ndgrid(1:F,1:F);
    K = exp(-((I-J).^2)/(2*sigma^2));
    K = (K + K')/2;               % 保证对称
    K = K / max(sum(K,2));        % 统一尺度（不破坏对称）
end

function W = make_uniform_weight(n)
% Ones off-diagonal, zeros on diagonal
    W = ones(n); W(1:n+1:end) = 0;
end

function L = laplacian_placeholder(n)
% 极简占位：全连通均匀图的归一化拉普拉斯（仅为维持接口）
    A = ones(n) - eye(n);
    d = sum(A,2);
    L = diag(d) - A;
end

function [Lnorm, info] = normalize_graph_laplacian(L, mode)
% 光谱归一化：把谱搬到 [0,2]（或单位谱半径）
    if nargin<2, mode='spectral'; end
    L = (L+L')/2;
    ev = eig(full(L));
    info.min_eig_before = min(real(ev));
    info.max_eig_before = max(real(ev));
    switch lower(mode)
        case 'spectral'
            s = max(1, info.max_eig_before);
            Lnorm = L / s;
        otherwise
            Lnorm = L;
    end
    ev2 = eig(full(Lnorm));
    info.min_eig_after = min(real(ev2));
    info.max_eig_after = max(real(ev2));
end

function Gamma_init = transport_init(Omega_prev, D_src, Sjj_tilde)
% 将上一轮源域精度搬运到当前白化域：Γ_init = D^{-1} Ω_prev D^{-1}
% 若上一轮为空，则用 "稳健逆 S̃" 作为启动（替代原来的 diag(1./diag(S̃))）。
    F = numel(D_src);
    Gamma_init = cell(F,1);
    if isempty(Omega_prev)
        for f = 1:F
            St = (Sjj_tilde{f} + Sjj_tilde{f}')/2;   % Hermitian
            % —— SPD 投影 + 特征值地板，避免病态/负特征值
            [U, D] = eig(full(St), 'vector');
            d = real(D);
            d = max(d, 1e-10);                       % 小特征值抬到阈值
            G = U * diag(1./d) * U';                 % 稳健逆
            G = (G + G')/2;                          % 数值对称化
            Gamma_init{f} = G;
        end
        return;
    end
    for f=1:F
        D = D_src{f};
        Gamma_init{f} = D \ (Omega_prev{f} / D);     % D^{-1} * Ω_prev * D^{-1}
        Gamma_init{f} = (Gamma_init{f} + Gamma_init{f}')/2;
    end
end

function Sigma_prior = invert_and_fix(Omega_cell, eps_ld)
% Robust inversion for next-round priors:
% - Hermitianize, force real diagonals
% - SPD projection via eigenvalue flooring
% - Invert by eigen decomposition (no chol)
    if nargin < 2, eps_ld = 1e-10; end
    F = numel(Omega_cell);
    n = size(Omega_cell{1},1);
    Sigma_prior = cell(F,1);

    for f = 1:F
        Om = Omega_cell{f};

        % --- Hermitianize & clean ---
        Om  = (Om + Om')/2;                         % Hermitian
        Om(~isfinite(Om)) = 0;                      % drop NaN/Inf safely
        d   = real(diag(Om));                       % real diagonal
        d   = max(d, eps_ld);
        Om(1:n+1:end) = d;

        % --- SPD projection via eig ---
        [U, S] = eig(full(Om), 'vector');           % S: eigenvalues (should be ~real)
        S  = real(S);
        S  = max(S, eps_ld);                        % floor small/neg eigenvalues
        if any(S < 2*eps_ld)
            S = max(S, 2*eps_ld);
        end

        % Invert via eig
        Sinv = 1 ./ S;
        Sigma = U * diag(Sinv) * U';

        % --- Hermitianize & real diagonal again ---
        Sigma = (Sigma + Sigma')/2;
        Sigma(1:n+1:end) = real(diag(Sigma));

        Sigma_prior{f} = Sigma;
    end
end

function [dOmega, dS] = compute_deltas(Omega, Omega_prev, Sjj, Sjj_prev)
% 计算外层相对变化（Frobenius）
    if isempty(Omega_prev), dOmega = inf; else
        num=0; den=0;
        for f=1:numel(Omega)
            num = num + norm(Omega{f}-Omega_prev{f}, 'fro');
            den = den + norm(Omega_prev{f}, 'fro');
        end
        dOmega = num / max(1, den);
    end
    if isempty(Sjj_prev), dS = inf; else
        num=0; den=0;
        for f=1:numel(Sjj)
            num = num + norm(Sjj{f}-Sjj_prev{f}, 'fro');
            den = den + norm(Sjj_prev{f}, 'fro');
        end
        dS = num / max(1, den);
    end
end

function g = read_grad_norm_(prox_res)
% 尝试从 prox_res 读取梯度范数（不同实现字段名不一）
    g = NaN;
    try
        if isfield(prox_res,'grad_norm'), g = prox_res.grad_norm; end
        if isnumeric(g) && ~isnan(g), return; end
        if isfield(prox_res,'history') && isfield(prox_res.history,'grad_norm')
            h = prox_res.history.grad_norm;
            if ~isempty(h), g = h(end); end
        end
        if isnumeric(g) && ~isnan(g), return; end
        if isfield(prox_res,'grad_norm_mean'), g = prox_res.grad_norm_mean; end
        if isnan(g), g = -1; end
    catch
        g = -1;
    end
end

function [pe, se] = edge_corr_abs(Acell, Bcell)
% 计算 |A|-|B| 的边级别相关（只用上三角）
    F = numel(Acell); ps = zeros(F,1); ss = zeros(F,1);
    for f=1:F
        A = abs(Acell{f}); B = abs(Bcell{f});
        idx = triu(true(size(A)),1);
        a = A(idx); b = B(idx);
        if isempty(a) || isempty(b) || all(a==a(1)) || all(b==b(1))
            ps(f)=0; ss(f)=0;
        else
            ps(f) = corr(a, b, 'type','Pearson','rows','complete');
            ss(f) = corr(a, b, 'type','Spearman','rows','complete');
        end
    end
    pe = mean(ps); se = mean(ss);
end

%% ========= Self-check helpers (add-on) =========
function diag_whitening_sanity(Sjj_hat, Sjj_tilde, D_src, iter_id)
% 目的：确认白化是否"真白化了"，以及 D*S*D 是否等于传给下游的 S̃
    F = numel(Sjj_hat);
    md_hat = zeros(F,1); sd_hat = zeros(F,1);
    md_til = zeros(F,1); sd_til = zeros(F,1);
    rel_err = zeros(F,1);
    is_spd  = true;

    for f = 1:F
        Sh = (Sjj_hat{f} + Sjj_hat{f}')/2;
        St = (Sjj_tilde{f} + Sjj_tilde{f}')/2;
        D  = D_src{f};

        md_hat(f) = mean(real(diag(Sh)));
        sd_hat(f) = std (real(diag(Sh)));
        md_til(f) = mean(real(diag(St)));
        sd_til(f) = std (real(diag(St)));

        % 验证 D*S*D 与 S̃ 的一致性（相对 Fro 误差）
        R = D*Sh*D - St;
        rel_err(f) = norm(R,'fro') / max(1, norm(St,'fro'));

        % 粗略 SPD 检查
        try
            ev = eig(St);
            is_spd = is_spd && all(real(ev) > -1e-10);
        catch
            is_spd = false;
        end
    end

    fprintf(['[DIAG][t=%d] <Whitening>\n' ...
             '    mean(diag S_hat)   = %.3g ± %.3g\n' ...
             '    mean(diag S_tilde) = %.3g ± %.3g  (理想≈1)\n' ...
             '    max relErr ||D*S*D - S_tilde||_F / ||S_tilde||_F = %.2e\n' ...
             '    S_tilde SPD (eig>0)? %s\n'], ...
             iter_id, mean(md_hat), mean(sd_hat), ...
             mean(md_til), mean(sd_til), ...
             max(rel_err), ternary(is_spd,'YES','NO'));

    if mean(md_til) > 2 || mean(md_til) < 0.5
        warning('[DIAG] whitened diagonals far from 1 (mean=%.3g). Whitening may be wrong.', mean(md_til));
    end
    if max(rel_err) > 1e-6
        warning('[DIAG] D*S*D and S_tilde mismatch (max rel err=%.2e). Check Module 1 output and pass-through.', max(rel_err));
    end
end

function diag_alpha_sanity(Gamma_init, Sjj_tilde, hp, params5, iter_id)
% 目的：用两个常见上界估计 Lipschitz L，并给出推荐 α，与 hp/alpha_max 对照
% L_byS ：近似把 Γ≈S̃^{-1} 时的上界 → max ||S̃||_2^2
% L_byG ：直接用当前 Γ 的上界 → max ||Γ^{-1}||_2^2
    [L_byS, L_byG] = estimate_L_candidates(Gamma_init, Sjj_tilde);
    L_est = min([L_byS, L_byG]);           % 取更紧的那个
    alpha_rec = 0.9 / max(L_est, 1e-12);   % 推荐步长（留 10% 裕度）

    % 读出"天花板"
    amax = [];
    if isfield(params5,'alpha_max') && ~isempty(params5.alpha_max)
        amax = params5.alpha_max;
    end

    fprintf(['[DIAG][t=%d] <Step-size>\n' ...
             '    L_byS ≈ %.3g, L_byG ≈ %.3g  ->  α_rec ≈ %.3g\n' ...
             '    hp.alpha = %.3g, alpha_max = %s\n'], ...
            iter_id, L_byS, L_byG, alpha_rec, hp.alpha, ...
            iif(isempty(amax),'[]',num2str(amax)));

    if ~isempty(amax) && amax < 0.1*alpha_rec
        warning('[DIAG] alpha_max (%.3g) << recommended α (%.3g). Consider increasing alpha_max or using line-search.', amax, alpha_rec);
    end
end

function [L_byS, L_byG] = estimate_L_candidates(Gamma_init, Sjj_tilde)
% 两种常用上界：
%   L_byS = max_f ||S̃_f||_2^2；
%   L_byG = max_f ||Γ_f^{-1}||_2^2（用最小特征值的倒数近似）。
    F = numel(Sjj_tilde);
    Ls = 0; Lg = 0;
    for f = 1:F
        St = (Sjj_tilde{f} + Sjj_tilde{f}')/2;
        s  = svds(St, 1);              % 最大奇异值 = 谱范数
        Ls = max(Ls, s^2);

        G  = (Gamma_init{f} + Gamma_init{f}')/2;
        try
            lam_min = min(real(eig(G)));
            if ~isfinite(lam_min) || lam_min <= 0
                invnorm = svds(pinv(G), 1);
            else
                invnorm = 1 / lam_min;
            end
        catch
            invnorm = svds(pinv(G), 1);
        end
        Lg = max(Lg, invnorm^2);
    end
    L_byS = Ls; L_byG = Lg;
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
function out = iif(cond, a, b)
    if cond, out = a; else, out = b; end
end

function warm = eloreta_warmstart_from_covs(Svv_cell, L, gamma_grid, opts)
% ELORETA_WARMSTART_FROM_COVS  — 多频 eLORETA → Sjj_e → Ω_init
% Inputs:
%   Svv_cell   : {F×1}, each p×p  (复厄米传感器互谱)
%   L          : p×n leadfield    (实/复均可)
%   gamma_grid : 含岭参数比例的扫描向量（如 logspace(-4,1,30)）
%   opts.maxit(50), opts.tol(1e-6), opts.verbose(false)
% Output (全部 {F×1}):
%   .Omega_init{f}  = inv_psd_robust(Sjj_e{f})    % 作为 M-step warm start
%   .Sjj_e{f}       = Tjv*Svv*Tjv'               % eLORETA 源协方差
%   .Tjv{f}         = n×p                         % eLORETA 反演算子
%   .gamma_opt{f}, .gcv_curve{f}

    if nargin < 4, opts = struct(); end
    maxit = getf(opts,'maxit',50);
    tol   = getf(opts,'tol',1e-6);
    verb  = getf(opts,'verbose',false);

    F = numel(Svv_cell);
    Tjv_cell = cell(F,1); Sjj_cell = cell(F,1);
    Om0_cell = cell(F,1); gcv_cell = cell(F,1); gopt_cell = cell(F,1);

    for f = 1:F
        Svv = Svv_cell{f};
        [Tjv, Sjj, ~, gamma_opt, gcv] = eloreta_simple_(Svv, L, gamma_grid, maxit, tol, verb);

        Sjj = psd_project_(Sjj);                          % 厄米化 + PSD 投影
        Om0 = inv_psd_robust_(Sjj, 1e-8, 1e-12);          % 稳健逆作为 Ω^(0)

        Tjv_cell{f}  = Tjv;
        Sjj_cell{f}  = Sjj;
        Om0_cell{f}  = Om0;
        gcv_cell{f}  = gcv;
        gopt_cell{f} = gamma_opt;
    end

    warm = struct('Omega_init',{Om0_cell}, ...
                  'Sjj_e',{Sjj_cell}, ...
                  'Tjv',{Tjv_cell}, ...
                  'gamma_opt',{gopt_cell}, ...
                  'gcv_curve',{gcv_cell});
end

function [Tjv, Sjj, W, gamma_opt, gcv] = eloreta_simple_(Svv, L, gamma_grid, maxit, tol, verb)
% 单频 eLORETA + GCV（单方向版；复数安全）

    [p,n] = size(L);
    gcv = zeros(numel(gamma_grid),1);
    best.T = []; best.W = []; best.gamma = NaN; best.score = Inf;

    for k = 1:numel(gamma_grid)
        gamma = gamma_grid(k);
        % --- eLORETA 内循环：更新 W
        w = ones(n,1);                             % W = diag(w)
        for it=1:maxit
            Winv = diag(1./w);
            A    = hermi_( L*Winv*L' );            % A = L W^{-1} L^H
            alpha= gamma * trace(A)/p;             % 归一化岭
            M    = inv_psd_( A + alpha*eye(p) );

            w_old = w;
            for i=1:n
                li  = L(:,i);
                mii = real(li' * M * li);         % Mb = l_i^H M l_i
                w(i)= sqrt(max(mii, eps));
            end
            if norm(w-w_old)/max(1,norm(w_old)) < tol, break; end
        end
        Winv = diag(1./w);
        T    = Winv * L' * M;                      % n×p: W^{-1} L^H (A+αI)^{-1}

        Txiv = eye(p) - L*T;                       % I - L T
        num  = real(trace( hermi_(Txiv*Svv*Txiv') ))/p;
        den  = ( real(trace(Txiv))/p )^2 + eps;
        gcv(k) = num / den;

        if gcv(k) < best.score
            best.T = T; best.W = diag(w); best.gamma = gamma; best.score = gcv(k);
        end

        if verb && (mod(k,10)==1)
            fprintf('[eLORETA] gamma=%.3g, GCV=%.3g\n', gamma, gcv(k));
        end
    end

    Tjv       = best.T;
    gamma_opt = best.gamma;
    W         = best.W;
    Sjj       = Tjv * Svv * Tjv';
end

% ---------- tiny utils ----------
function A = hermi_(A), A = (A + A')/2; end
function X = inv_psd_(X)
    X = hermi_(X);
    [U,S] = eig(X,'vector'); S = max(real(S), eps);
    X = U*diag(1./S)*U'; X = hermi_(X);
end
function S = psd_project_(S)
    S = hermi_(S);
    [U,d] = eig(S,'vector'); d = max(real(d), 0);
    S = U*diag(d)*U'; S = hermi_(S);
end
function v = getf(s, f, d), if isfield(s,f) && ~isempty(s.(f)), v=s.(f); else, v=d; end, end
function Om = inv_psd_robust_(A, eps_reg, min_ratio)
    A = hermi_(A);
    [V,D] = eig(A); d = real(diag(D)); dmax = max(d);
    floor_val = max(min_ratio * max(dmax, eps), 0);
    d(d<floor_val) = floor_val;
    if eps_reg > 0, d = (d + eps_reg * dmax) / (1 + eps_reg); end
    Om = V * diag(1./d) * V'; Om = hermi_(Om);
end