%% PIPELINE_SOURCE_EM  (True EM loop; scalar estep_in version, strengthened prior + frozen whitening)
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
em = struct();
em.max_em_iter      = 50;
em.tol_Omega        = 5e-3;     % 外层 Ω 相对变化阈值（0.5%）
em.tol_S            = 1e-3;     % 外层 S 相对变化阈值
em.min_hold_iters   = 2;        % 连续满足两轮再停
em.update_act_every = 2;        % 活跃集更新频率（前两轮必更）
em.recfg_hp_every   = 1;        % 超参数每轮重估
em.lambda2_anneal   = struct('enable', true, 'grow', 1.05, 'max_scale', 2.0); % 缓慢增大
em.lambda3_ratio    = 0.3;      % λ3 = ratio * λ1
% —— 强化先验 + 前几轮冻结白化 ——
em.prior_floor      = 1e-3;     % << 从 1e-10 提高，避免先验被“放大成巨大协方差”
em.prior_strength   = 5;        % 先验精度的强度（等价于协方差 /5）
em.prior_strength_decay = 0.85; % 每轮衰减，最终回到 1 左右
em.freeze_whiten_iters = 2;     % 前两轮冻结白化矩阵

em.inner = struct( ...          % M-step(Prox) 内层参数
    'alpha0', 5e-3, 'alpha_max', 2e-3, 'alpha_up', 1.05, 'alpha_down', 0.7, ...
    'alpha_grow_patience', 2, 'max_iter', 100, 'obj_improve_tol', 1e-6, ...
    'weight_mode','hadamard', 'use_graph_laplacian', true, ...
    'diag', struct('enable', true, 'update_every', 1, 'metrics_every', 1, ...
                   'print_every', 1, 'f_view', 1, 'log_csv', 'prox_trace.csv', ...
                   'keep_fig_open', true));

%% 4) State for EM
A_masks = [];                     % 活跃集（滞后更新）
Omega_prev   = [];                % 上一轮源域精度
Sjj_prev     = [];                % 上一轮 E 步二阶矩
hold_counter = 0;                 % 连续满足计数
lambda2_base = [];                % 记录初始 λ2 建议值（用于 anneal 上限）
hp_last      = [];                % 记录上一轮超参数
D_src_fixed  = [];                % 冻结白化所用

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

    %% Preprocessing / Whitening in source domain (Module 1)
    pre = module1_preproc_from_covset(Sjj_hat, struct( ...
        'smoothing_method','moving_average', ...
        'loading_factor',  1e-6, ...
        'min_power',       1e-10, ...
        'verbose',         false ));
    D_src     = pre.D;            % {F×1}, n×n
    Sjj_tilde = pre.Sigma_tilde;  % {F×1}, n×n

    % —— 冻结白化：前 em.freeze_whiten_iters 轮使用首轮的 D —— 
    if t == 1
        D_src_fixed = D_src;
    elseif t <= em.freeze_whiten_iters
        D_src = D_src_fixed;
        Sjj_tilde = cellfun(@(S,D) (D \ S) / D, Sjj_hat, D_src, 'uni', 0);
        fprintf('[EM] Whitening frozen (use t=1 D) at iter %d.\n', t);
    end

    %% Active set（前2轮必更新，否则每 update_act_every 轮更新）
    if t <= 5 || isempty(A_masks) || mod(t, em.update_act_every) == 0
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
    Gamma_init = transport_init(Omega_prev, D_src, Sjj_tilde);

    input_data_m5 = struct();
    input_data_m5.whitened_covariances = Sjj_tilde;
    input_data_m5.initial_precision    = Gamma_init;
    input_data_m5.smoothing_kernel     = K;
    input_data_m5.weight_matrix        = W;
    input_data_m5.active_set_mask      = {A_masks{:}};
    input_data_m5.whitening_matrices   = D_src;                   % 仅供 Live 可视化

    params5 = struct( ...
        'lambda1', lambda1, ...
        'lambda2', lambda2_effective, ...
        'lambda2_suggested', hp.lambda2_suggested, ...
        'alpha0',  alpha0, ...
        'max_iter', em.inner.max_iter, 'verbose', true, ...
        'active_set_update_freq', 10, ...
        'alpha_max', em.inner.alpha_max, 'alpha_up', em.inner.alpha_up, ...
        'alpha_down', em.inner.alpha_down, 'alpha_grow_patience', em.inner.alpha_grow_patience, ...
        'obj_improve_tol', em.inner.obj_improve_tol, ...
        'weight_mode', em.inner.weight_mode, 'use_graph_laplacian', em.inner.use_graph_laplacian, ...
        'diag', em.inner.diag );

    % 在线可视化（Live）
    params5.diag.live_plot = live_plot_cfg;

    % 空间平滑
    if lambda3 > 0
        params5.lambda3                    = lambda3;
        params5.spatial_graph_matrix       = Lsp;
        params5.spatial_graph_is_laplacian = true;
        if ~isfield(params5,'spatial_weight_mode') || isempty(params5.spatial_weight_mode)
            params5.spatial_weight_mode = 'node';      % 'node' or 'hadamard'
        end
    end

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

    % —— 量级自检：先验精度 vs 观测项 ——
    A_like = (L' / Sigma_xixi) * L;   % == L' * inv(Sigma_xixi) * L
    nl = norm(A_like, 'fro');
    np = mean(cellfun(@(Om) norm(Om, 'fro'), Omega_src));
    fprintf('[scale] ||L^H Σ^{-1} L||_F≈%.2f, mean||Ω||_F≈%.2f, prior_strength=%.2f\n', nl, np, em.prior_strength);

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

    %% 回写先验进入下一轮（强化 + 衰减）
    Sigma_prior = invert_and_fix(Omega_src, em.prior_floor);           % 更高地板
    Sigma_prior = cellfun(@(S) S / em.prior_strength, Sigma_prior, 'uni', 0); % 强化先验（等价于提升精度）
    estep_in.source_prior_covariances = Sigma_prior;
    em.prior_strength = max(1, em.prior_strength * em.prior_strength_decay);

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

opt_compare = struct('mode','match_sparsity', 'f_view', 1, 'show_pr', true, 'show_roc', true);
opt_compare.Gamma = Gamma_tilde_star;
metrics = gt_compare_and_plot(Omega_true, Omega_src, opt_compare);

Sigma_hat_test   = ensure_sigma_collection(Sjj_hat);
Sigma_tildeteest = ensure_sigma_collection(Sjj_tilde);
debug_em_convergence(Sigma_hat_test,   'weight_mode','hadamard', 'use_laplacian',true);
debug_em_convergence(Sigma_tildeteest, 'weight_mode','hadamard', 'use_laplacian',true);

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

function Gamma_init = transport_init(Omega_prev, D_src, Sjj_tilde)
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


function Sigma_prior = invert_and_fix(Omega_cell, eps_ld)
% Robust inversion for next-round priors:
% - Hermitianize, force real diagonals
% - SPD projection via eigenvalue flooring
% - Invert by eigen decomposition (no chol)
    if nargin < 2, eps_ld = 1e-3; end  % << 默认就更高一些
    F = numel(Omega_cell);
    n = size(Omega_cell{1},1);
    Sigma_prior = cell(F,1);

    for f = 1:F
        Om = Omega_cell{f};
        Om  = (Om + Om')/2;            % Hermitian
        Om(~isfinite(Om)) = 0;
        d   = real(diag(Om));
        d   = max(d, eps_ld);
        Om(1:n+1:end) = d;

        [U, S] = eig(full(Om), 'vector');
        S  = real(S);
        S  = max(S, eps_ld);
        if any(S < 2*eps_ld), S = max(S, 2*eps_ld); end

        Sinv = 1 ./ S;
        Sigma = U * diag(Sinv) * U';

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
