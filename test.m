%% PIPELINE_SOURCE_EM  (True EM loop; scalar estep_in version)
% Build estep_in as a *scalar struct* with per-frequency matrices stored as {F×1} cells.

%% 0) Simulation (Module 7)
n  = 10;  p  = 3;  F  = 3;  T  = 4096;

[Omega_true, Sigma_true, emp_covariance, sim] = module7_simulation_improved_complex( ...
    'n_nodes', n, 'n_sensors', p, 'n_freq', F, 'n_samples', T, ...
    'generate_leadfield', true, ...
    'leadfield_type', 'simple', ...
    'random_seed', 42);

L          = sim.leadfield_matrix;   % p×n
Sigma_xixi = sim.Sigma_xixi;         % p×p

%% 1) Build E-step input as a *scalar struct* (NOT a struct array)
emp_cov_cell = coerce_cov_cell(emp_covariance);     % supports cell / 3D / single
F            = numel(emp_cov_cell);
prior_cell   = repmat({eye(n)}, F, 1);              % Σ_j^(0)

estep_in = struct();                                  % scalar struct
estep_in.leadfield_matrix         = L;                % p×n
estep_in.empirical_covariances    = emp_cov_cell;     % {F×1}, p×p
estep_in.source_prior_covariances = prior_cell;       % {F×1}, n×n
estep_in.noise_covariance         = Sigma_xixi;       % p×p
estep_in.frequencies              = (1:F);            % 1×F

assert(isstruct(estep_in) && isscalar(estep_in), 'estep_in must be a scalar struct');
assert(isa(estep_in.empirical_covariances,'cell') && numel(estep_in.empirical_covariances)==F);
assert(all(cellfun(@(A) isequal(size(A),[p p]), estep_in.empirical_covariances)));
assert(isa(estep_in.source_prior_covariances,'cell') && numel(estep_in.source_prior_covariances)==F);
assert(all(cellfun(@(S) isequal(size(S),[n n]), estep_in.source_prior_covariances)));

%% 1.5) eLORETA warm start (Plan A): build Ω^(0) from {Svv_f}
gamma_grid = logspace(-4, 1, 30);                 % 可按需调整扫描范围与密度
opts_el    = struct('maxit', 50, 'tol', 1e-6, 'verbose', false);

warm = eloreta_warmstart_from_covs(emp_cov_cell, L, gamma_grid, opts_el);

% 作为 M-step 的 warm start：上一轮源域精度 Ω_prev
Omega_prev = warm.Omega_init;                 % {F×1}, 每个 n×n

% （可选）也把 eLORETA 的源协方差作为上一轮二阶矩，便于 ΔS 统计更稳定
Sjj_prev   = warm.Sjj_e;                          % {F×1}, 每个 n×n

% （可选）若希望第一轮 E-step 的先验不再是 I，可用下行把 prior 设为 inv(Ω_init)
% estep_in.source_prior_covariances = invert_and_fix(Omega_prev, 1e-10);


%% 2) Fixed resources (K, W, spatial L) —— once for all EM iterations
K = make_frequency_kernel(F, 3.0);    % single source of truth
W = make_uniform_weight(n);

% Spatial Laplacian (spectrally normalized) — if available in sim, use it.
if isfield(sim, 'source_graph') && isfield(sim.source_graph, 'L')
    L_raw = sim.source_graph.L;
else
    % Fallback: harmless placeholder (不会用于推理，仅保接口完整)
    L_raw = laplacian_placeholder(n);
end
[Lsp, Linfo] = normalize_graph_laplacian(L_raw, 'spectral');
fprintf('Spatial L normalized: eig[min,max]=[%.2e, %.2e] -> [%.2e, %.2e]\n', ...
    Linfo.min_eig_before, Linfo.max_eig_before, Linfo.min_eig_after, Linfo.max_eig_after);

%% 3) EM loop settings
%% 3) EM loop settings  —— 10轮外层 + 快速内层
em = struct();
em.max_em_iter      = 10;      % 外层 EM 10 轮
em.tol_Omega        = 5e-3;
em.tol_S            = 1e-3;
em.min_hold_iters   = 2;
em.update_act_every = 2;
em.recfg_hp_every   = 1;

% λ2 不再持续增大（避免后期过度平滑拖慢进度）
em.lambda2_anneal   = struct('enable', false, 'grow', 1.02, 'max_scale', 1.3);
em.lambda3_ratio    = 0.3;

% —— 内层(Prox)关键：放开步长上限、减少最大迭代数、加速自适应
em.inner = struct( ...
    'alpha0',  0.1, ...          % 初始步长：给一个“可迈得动”的值
    'alpha_max', 1.0, ...        % 允许快速放大到 ~O(1) 的水平
    'alpha_up', 1.5, ...         % 放大倍率更激进
    'alpha_down', 0.6, ...       % 退火时稍保守
    'alpha_grow_patience', 1, ...% 每 1 步就允许尝试放大
    'max_iter', 30, ...          % 每次 M-step 最多 30 步（原来100）
    'obj_improve_tol', 5e-6, ...
    'weight_mode','hadamard', ...
    'use_graph_laplacian', true, ...
    'diag', struct('enable', true, 'update_every', 1, 'metrics_every', 1, ...
                   'print_every', 1, 'f_view', 1, 'log_csv', 'prox_trace.csv', ...
                   'keep_fig_open', true));


%% 4) State for EM
A_masks = [];                     % 活跃集（滞后更新）

% >>> 不要覆盖 1.5 节 warm 的结果；仅在未设置时才置空/默认
if ~exist('Omega_prev','var') || isempty(Omega_prev)
    Omega_prev = [];              % 若你想强制不用 warm，再手动清空
end
if ~exist('Sjj_prev','var') || isempty(Sjj_prev)
    Sjj_prev   = [];
end

hold_counter = 0;
lambda2_base = [];
hp_last      = [];


% 可视化设置
live_plot_cfg = struct('enable', true, 'f_view', 1, 'plot_every', 5, ...
    'value_mode','abs', 'ground_truth_domain','source');
live_plot_cfg.ground_truth_precision = Omega_true;

fprintf('\n===== EM starts (max %d iters) =====\n', em.max_em_iter);
for t = 1:em.max_em_iter
    fprintf('\n--- EM iter %d ---\n', t);

    %% E-step (Module 2 + wrapper) —— 数值保障：开启 PSD & 统一加载
    estep_out = module2_estep(estep_in, struct( ...
        'ensure_hermitian', true, ...
        'ensure_real_diag', true, ...
        'ensure_psd',       true, ...
        'psd_tol',          1e-10, ...
        'diag_loading',     1e-10 ));
    Sjj_hat = estep_out.source_second_moments;   % {F×1}, n×n

    %% ===== Noise M-step (Plan B: full SPD sensor noise) =====
    % Σ_{ξξ,new} = (1/F) * Σ_f  S_res(f)
    % S_res(f) = (I - L*T_jv(f)) * Svv(f) * (I - L*T_jv(f))' + L*Σ_{i|v}(f)*L'
    I_p = eye(p);
    S_xixi_accum = zeros(p,p,'like',estep_in.noise_covariance);

    for f = 1:F
        % 后验协方差 Σ_{i|v}(f)
        Sigma_post_f = module2_posterior_source_covariance( ...
            estep_in.source_prior_covariances{f}, ...
            L, ...
            estep_in.noise_covariance, ...
            struct('regularization_factor',1e-8,'verbose',false));

        % 传递算子 T_{jv}(f)
        T_jv_f  = Sigma_post_f * (L' / estep_in.noise_covariance);

        % 残差传递 T_{ξv}(f) = I - L*T_{jv}(f)
        T_xi_v_f = I_p - L * T_jv_f;

        % 经验残差协方差
        S_res_f = T_xi_v_f * emp_cov_cell{f} * T_xi_v_f' ...
                + L * Sigma_post_f * L';

        % 累加并对称化
        S_res_f = 0.5*(S_res_f + S_res_f');
        S_xixi_accum = S_xixi_accum + S_res_f;
    end

    % 频点平均（MLE 闭式）
    Sigma_xixi_new = S_xixi_accum / F;

    % 稳健化：对称化 + 轻微对角加载（按迹归一）
    Sigma_xixi_new = 0.5*(Sigma_xixi_new + Sigma_xixi_new');
    ridge = 1e-10 * trace(Sigma_xixi_new)/p;
    Sigma_xixi_new = Sigma_xixi_new + ridge * eye(p,'like',Sigma_xixi_new);

    % （可选）EM 指数平滑，抑制抖动
    eta = 0.3;  % 0~1，小=更平滑
    Sigma_xixi = (1-eta)*estep_in.noise_covariance + eta*Sigma_xixi_new;

    % 确保 SPD（必要时再微量加载）
    [~,chol_flag] = chol(0.5*(Sigma_xixi + Sigma_xixi'),'lower');
    if chol_flag ~= 0
        Sigma_xixi = Sigma_xixi + 1e-6*eye(p,'like',Sigma_xixi);
    end

    % 回写：下一轮 E 步使用更新后的噪声
    estep_in.noise_covariance = Sigma_xixi;
    %% ===== End Noise M-step =====

    %% Preprocessing / Whitening in source domain (Module 1)
    pre = module1_preproc_from_covset(Sjj_hat, struct( ...
        'smoothing_method','moving_average', ...
        'loading_factor',  1e-6, ...
        'min_power',       1e-10, ...
        'verbose',         false ));
    D_src     = pre.D;            % {F×1}, n×n
    Sjj_tilde = pre.Sigma_tilde;  % {F×1}, n×n

    %% ======================= [SELF-CHECK-1] Whitening sanities =======================
    diag_whitening_sanity(Sjj_hat, Sjj_tilde, D_src, t);
    % ================================================================================

    %% Active set（前2轮必更新，否则每 update_act_every 轮更新）
    if t <= 2 || isempty(A_masks) || mod(t, em.update_act_every) == 0
        input_data_m3 = struct();
        input_data_m3.whitened_covariances = Sjj_tilde;
        input_data_m3.frequencies          = 1:F;
        act = module3_active_set(input_data_m3, struct( ...
            'proxy_method','correlation', 'quantile_level',0.25, ...
            'force_diagonal_active', true, 'verbose', false));
        A_masks = arrayfun(@(f) logical(act.combined_active_mask(:,:,f)), 1:F, 'uni', 0);
        act_mark = '*updated*';
    else
        act_mark = 'kept';
    end
    fprintf('[EM] Active set %s.\n', act_mark);

    %% Hyperparameters（每轮重估；λ2 可缓慢增大）
    if mod(t, em.recfg_hp_every) == 0 || isempty(hp_last)
        input_data_m6 = struct();
        input_data_m6.whitened_covariances = Sjj_tilde;
        input_data_m6.kernel_matrix        = K;
        input_data_m6.weight_matrix        = W;
        input_data_m6.active_set_mask      = {A_masks{:}};
        hp = module6_hyperparameter_config(input_data_m6, struct('use_gershgorin', true));
        hp_last = hp;
        if isempty(lambda2_base), lambda2_base = hp.lambda2_suggested; end
    else
        hp = hp_last;
    end
    lambda1 = hp.lambda1;
    % λ2: 建议值基础上做 anneal
    if em.lambda2_anneal.enable
        grow = em.lambda2_anneal.grow^(max(0,t-1));
        grow = min(grow, em.lambda2_anneal.max_scale);
        lambda2_effective = hp.lambda2_suggested * grow;
    else
        lambda2_effective = hp.lambda2_suggested;
    end
    lambda3 = em.lambda3_ratio * lambda1;
    alpha0  = hp.alpha;

    fprintf('[HP ] λ1=%.3g, λ2(sugg)=%.3g -> λ2(eff)=%.3g, λ3=%.3g, α=%.3g\n', ...
        lambda1, hp.lambda2_suggested, lambda2_effective, lambda3, alpha0);

    %% M-step（Module 5，白化域）
    % Warm start transport: Γ_init = D^{-1} Ω_prev D^{-1}
    % Warm start transport: Γ_init = D^{-1} Ω_prev D^{-1}   (首轮亦如此，Ω_prev 来自 eLORETA)
    Gamma_init = transport_init(Omega_prev, D_src, Sjj_tilde);



    input_data_m5 = struct();
    input_data_m5.whitened_covariances = Sjj_tilde;
    input_data_m5.initial_precision    = Gamma_init;
    input_data_m5.smoothing_kernel     = K;
    input_data_m5.weight_matrix        = W;
    input_data_m5.active_set_mask      = {A_masks{:}};
    input_data_m5.whitening_matrices   = D_src;                   % 仅供 Live 可视化
alpha0_override = max(0.05, min(0.4, hp.alpha * 200));   % 经验兜底：把 hp.alpha 提升 ~200 倍但不超过 0.4
em.inner.alpha0 = alpha0_override;                       % 回写到 em.inner 供后面使用
fprintf('[STEP] alpha0(hp)=%.3g -> alpha0(used)=%.3g, alpha_max=%.2f\n', ...
        hp.alpha, em.inner.alpha0, em.inner.alpha_max);
    % params5 = struct( ...
    %     'lambda1', lambda1, ...
    %     'lambda2', lambda2_effective, ...
    %     'lambda2_suggested', hp.lambda2_suggested, ...
    %     'alpha0',  alpha0, ...
    %     'max_iter', em.inner.max_iter, 'verbose', true, ...
    %     'active_set_update_freq', 10, ...
    %     'alpha_max', em.inner.alpha_max, 'alpha_up', em.inner.alpha_up, ...
    %     'alpha_down', em.inner.alpha_down, 'alpha_grow_patience', em.inner.alpha_grow_patience, ...
    %     'obj_improve_tol', em.inner.obj_improve_tol, ...
    %     'weight_mode', em.inner.weight_mode, 'use_graph_laplacian', em.inner.use_graph_laplacian, ...
    %     'diag', em.inner.diag );
    %% 组装 M-step 参数（保持你现有的其它字段不变）
params5 = struct( ...
    'lambda1', lambda1, ...
    'lambda2', lambda2_effective, ...
    'lambda2_suggested', hp.lambda2_suggested, ...
    'alpha0',  em.inner.alpha0, ...      % <== 使用上面的 override
    'max_iter', em.inner.max_iter, 'verbose', true, ...
    'active_set_update_freq', 10, ...
    'alpha_max', em.inner.alpha_max, 'alpha_up', em.inner.alpha_up, ...
    'alpha_down', em.inner.alpha_down, 'alpha_grow_patience', em.inner.alpha_grow_patience, ...
    'obj_improve_tol', em.inner.obj_improve_tol, ...
    'weight_mode', em.inner.weight_mode, 'use_graph_laplacian', em.inner.use_graph_laplacian, ...
    'diag', em.inner.diag );

    % 在线可视化（Live）
    params5.diag.live_plot = live_plot_cfg;

    params5.alpha_min              = 1e-4;   % 防止 α 被回溯到数值噪声
params5.armijo_c1              = 1e-4;   % Armijo 常数
params5.backtrack_beta         = 0.5;    % 回溯时 α 缩放
params5.max_backtrack_per_iter = 20;     % 每步最多回溯次数

% 回溯累计触发“优先降 λ2”（退火），避免一直只缩 α
params5.backtrack_patience     = 10;     % 连续回溯阈值
params5.lambda2_decay_factor   = 0.9;    % 退火倍率
params5.lambda2_min            = 0.01 * params5.lambda2;  % 不把 λ2 降到 0

% 是否对角也做 L1（保持和你现在一致：不惩罚对角）
params5.penalize_diagonal      = false;


    % 空间平滑
    if lambda3 > 0
        params5.lambda3                    = lambda3;
        params5.spatial_graph_matrix       = Lsp;
        params5.spatial_graph_is_laplacian = true;
        if ~isfield(params5,'spatial_weight_mode') || isempty(params5.spatial_weight_mode)
            params5.spatial_weight_mode = 'node';      % 'node' or 'hadamard'
        end
    end

    %% ======================= [SELF-CHECK-2] Step-size / Lipschitz ====================
    diag_alpha_sanity(Gamma_init, Sjj_tilde, hp, params5, t);
    % ================================================================================

    [Gamma_tilde_star, prox_res] = module5_proximal(input_data_m5, params5);

    %% Recolor to source scale (Module 8)
    input_data_m8 = [];
    input_data_m8.whitened_precision_matrices = Gamma_tilde_star;
    input_data_m8.whitening_matrices = D_src;
    recol = module8_recoloring(input_data_m8, struct());
    Omega_src = recol.recolored_precision_matrices;

    %% 外层收敛统计
    [delta_Omega, delta_S] = compute_deltas(Omega_src, Omega_prev, Sjj_hat, Sjj_prev);
    gnorm = read_grad_norm_(prox_res);
    fprintf('[EM] ΔΩ=%.3e, ΔS=%.3e, ||grad||≈%.3g\n', delta_Omega, delta_S, gnorm);

    % 终止条件：两项同时小于阈值则累计一次；反之清零
    if delta_Omega < em.tol_Omega && delta_S < em.tol_S
        hold_counter = hold_counter + 1;
    else
        hold_counter = 0;
    end
    if hold_counter >= em.min_hold_iters
        fprintf('>>> EM converged by (ΔΩ,ΔS) for %d consecutive iters at t=%d.\n', hold_counter, t);
        break;
    end

    %% 回写先验进入下一轮
    estep_in.source_prior_covariances = invert_and_fix(Omega_src, 1e-10);

    %% 缓存进入下一轮
    Omega_prev = Omega_src;
    Sjj_prev   = Sjj_hat;
end

fprintf('===== EM finished at iter %d =====\n', t);

%% 5) 评估/可视化（最终一次）
f_view = 1;
Om     = Omega_src{f_view};
pcorr  = abs(-Om) ./ sqrt((abs(diag(Om))+eps) * (abs(diag(Om))+eps)');
pcorr(1:n+1:end) = 0;
fprintf('Done. Example partial coherence at f=%d: max=%g, median=%g\n', ...
    f_view, max(pcorr(:)), median(pcorr(pcorr>0)));

viz_pipeline_summary(prox_res, Gamma_tilde_star, Sjj_tilde, K, A_masks, Omega_src, Omega_true);

% ====== 使用“源域 Ω_src”做打分（而不是白化域 Γ）======
opt_compare = struct('mode','match_sparsity', 'f_view', 1, ...
                     'show_pr', true, 'show_roc', true);
opt_compare.Gamma = Omega_src;                % 与 Ω_true 同域同尺度
metrics = gt_compare_and_plot(Omega_true, Omega_src, opt_compare);

%% ====== SELF-CHECK: ranking and metric consistency ======
% Inputs assumed ready in workspace:
%   - Omega_src  : {F×1} estimated precisions (source domain)
%   - Omega_true : {F×1} ground-truth precisions (source domain)
% (如果 Omega_true 只有一个频率，也可以传 {F×1} 的复本)

opts_sc = struct('zero_tol',1e-12, 'use_abs',true, 'per_freq',true);

sc1 = sanity_edge_ranking(Omega_src, Omega_true, opts_sc, 'raw_abs');    % |Ω|
sc2 = sanity_edge_ranking(Omega_src, Omega_true, opts_sc, 'pcorr_abs');  % |partial correlation|

fprintf('\n[SELF-CHECK] ===== Summary =====\n');
fprintf('RAW |Ω|     : AUROC=%.3f (vs -score: %.3f) | AUPR=%.3f | ρ(|Ω_est|,|Ω_true|)=%.3f\n', ...
    sc1.auroc_pos, sc1.auroc_neg, sc1.aupr, sc1.corr_pos);
fprintf('PCorr |Ω|   : AUROC=%.3f (vs -score: %.3f) | AUPR=%.3f | ρ(|Ω_est|,|Ω_true|)=%.3f\n', ...
    sc2.auroc_pos, sc2.auroc_neg, sc2.aupr, sc2.corr_pos);

if sc1.auroc_neg > sc1.auroc_pos + 0.03 || sc2.auroc_neg > sc2.auroc_pos + 0.03
    fprintf('>>> [ALERT] 排序方向疑似反了：用 -score 的 AUROC 更大。请检查评分是否取了负号/倒数。\n');
end

if isfield(sc1,'perF') && ~isempty(sc1.perF)
    fprintf('[Per-freq AUROC] raw |Ω|:'); fprintf(' %.3f', sc1.perF); fprintf('\n');
    fprintf('[Per-freq AUROC] pcorr  :'); fprintf(' %.3f', sc2.perF); fprintf('\n');
end

fprintf('[Density] est>0 (raw) ratio=%.3f, pos_edges_ratio=%.3f\n', sc1.est_density, sc1.pos_ratio);

%% ===== helper function =====
function sc = sanity_edge_ranking(Om_est_cell, Om_true_cell, opts, mode)
    if nargin<4, mode='raw_abs'; end
    if nargin<3 || isempty(opts), opts = struct(); end
    ztol = getfield_with_default(opts,'zero_tol',1e-12);
    use_abs = getfield_with_default(opts,'use_abs',true);
    want_perF = getfield_with_default(opts,'per_freq',true);

    F = numel(Om_est_cell);
    n = size(Om_est_cell{1},1);
    idxUT = triu(true(n),1);

    s_all = []; y_all = []; sgt_all = []; perF = nan(F,1);
    for f=1:F
        E = Om_est_cell{f}; T = Om_true_cell{min(f, numel(Om_true_cell))};
        switch mode
            case 'raw_abs'
                s = E(idxUT);
                if use_abs, s = abs(s); end
            case 'pcorr_abs'
                d = sqrt( max(real(diag(E)), ztol) );
                D = d*d.';
                pc = -E ./ (D + eps);
                s = pc(idxUT);
                if use_abs, s = abs(s); end
            otherwise
                error('unknown mode');
        end
        t = T(idxUT);
        y = abs(t) > ztol;                % 1=有边, 0=无边
        sgt = abs(t);                      % 真值强度（只用于相关性）

        s_all   = [s_all;   s];
        y_all   = [y_all;   y];
        sgt_all = [sgt_all; sgt];

        if want_perF
            try
                [~,~,~,perF(f)] = perfcurve(y, s, true);
            catch
                perF(f) = NaN;
            end
        end
    end

    % AUROC/AUPR（正向/反向）
    [~,~,~,auroc_pos] = perfcurve(y_all,  s_all,  true);
    [~,~,~,auroc_neg] = perfcurve(y_all, -s_all, true);
    [rec,prec,~,aupr] = perfcurve(y_all, s_all, true, 'xCrit','reca','yCrit','prec'); %#ok<ASGLU>

    % 只在正样本上看强度相关（是否把强边排在前面）
    sel_pos = (y_all==1);
    if any(sel_pos)
        try
            c = corr(s_all(sel_pos), sgt_all(sel_pos), 'type','Spearman','rows','complete');
        catch
            c = NaN;
        end
    else
        c = NaN;
    end

    sc = struct('auroc_pos',auroc_pos,'auroc_neg',auroc_neg, ...
                'aupr',aupr,'corr_pos',c,'perF',perF);

    % 粗略密度信息
    est_density = mean(s_all > max(ztol, 1e-9));
    pos_ratio   = mean(y_all==1);
    sc.est_density = est_density; sc.pos_ratio = pos_ratio;
end

function v = getfield_with_default(s, name, def)
    if isfield(s,name) && ~isempty(s.(name)), v = s.(name); else, v = def; end
end

% —— 额外快速一致性/相关性检查（建议保留）——
[pe,se] = edge_corr_abs(Omega_src, Omega_true);
fprintf('[CorrΩ ] |Omega_est|-|Omega_true| : Pearson=%.3f | Spearman=%.3f (mean over F)\n', pe, se);

Sigma_hat_test   = ensure_sigma_collection(Sjj_hat);
Sigma_tildeteest = ensure_sigma_collection(Sjj_tilde);

% 让诊断工具与主循环 K/W 对齐；若老版本不支持透传参数，就回退到原调用
try
    debug_em_convergence(Sigma_hat_test,   'weight_mode','hadamard','use_laplacian',true, ...
                         'kernel_matrix',K,'weight_matrix',W);
    debug_em_convergence(Sigma_tildeteest, 'weight_mode','hadamard','use_laplacian',true, ...
                         'kernel_matrix',K,'weight_matrix',W);
catch
    debug_em_convergence(Sigma_hat_test,   'weight_mode','hadamard','use_laplacian',true);
    debug_em_convergence(Sigma_tildeteest, 'weight_mode','hadamard','use_laplacian',true);
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
% 若上一轮为空，则用 “稳健逆 S̃” 作为启动（替代原来的 diag(1./diag(S̃))）。
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
            Gamma_init{f} = D \ (Omega_prev{f} / D);     % D^{-1} * Ω_prev * D^{-1)
        end
        return;
    end
    for f=1:F
        D = D_src{f};
        Gamma_init{f} = (D \ Omega_prev{f}) / D;     % D^{-1} * Ω_prev * D^{-1}
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
% 目的：确认白化是否“真白化了”，以及 D*S*D 是否等于传给下游的 S̃
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

    % 读出“天花板”
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
