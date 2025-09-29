%% PIPELINE_SOURCE_EM  (Source-domain EM pipeline → extended to True EM loop)
% 保留现有 Source 版设定，仅扩展为外层EM=10轮；第1轮仍用 Omega_true 做 M步初值以便可视化对照

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

%% 2) Fixed resources: frequency kernel, weights, spatial Laplacian (spectrally normalized)
K = make_frequency_kernel(F, 3.0);    % single source of truth
W = make_uniform_weight(n);

if isfield(sim, 'source_graph') && isfield(sim.source_graph, 'L')
    L_raw = sim.source_graph.L;
else
    % 简单兜底（不用于推理，仅为接口完整）
    A = ones(n) - eye(n);
    d = sum(A,2);
    L_raw = diag(d) - A;
end
[Lsp, Linfo] = normalize_graph_laplacian(L_raw, 'spectral');
fprintf('Spatial L normalized: eig[min,max]=[%.2e, %.2e] -> [%.2e, %.2e]\n', ...
    Linfo.min_eig_before, Linfo.max_eig_before, Linfo.min_eig_after, Linfo.max_eig_after);

%% 3) EM loop settings（只做外层包装，保持你现有 prox/HP 设定不变）
em = struct();
em.max_em_iter      = 10;
em.tol_Omega        = 5e-3;
em.tol_S            = 1e-3;
em.min_hold_iters   = 2;          % ΔΩ与ΔS连续满足两轮才停
em.update_act_every = 2;          % 活跃集更新频率
em.recfg_hp_every   = 1;          % 每轮重估超参

%% 4) State for EM
A_masks   = [];     % 活跃集（滞后更新）
Omega_prev= [];     % 上一轮源域精度
Sjj_prev  = [];     % 上一轮 E 步二阶矩
hold_counter = 0;   % 连续满足次数

% Live 诊断（保持你的 Source 版 ground truth 对照）
live_plot_cfg = struct('enable', true, 'f_view', 1, 'plot_every', 5, ...
    'value_mode','abs', 'ground_truth_domain','source');
live_plot_cfg.ground_truth_precision = Omega_true;

fprintf('\n===== EM starts (max %d iters) =====\n', em.max_em_iter);
for t = 1:em.max_em_iter
    fprintf('\n--- EM iter %d ---\n', t);

    %% 5) E-step (Module 2 + wrapper) —— PSD & 统一加载（保持你的设置）
    estep_out = module2_estep(estep_in, struct( ...
        'ensure_hermitian', true, ...
        'ensure_real_diag', true, ...
        'ensure_psd',       true, ...
        'psd_tol',          1e-10, ...
        'diag_loading',     1e-10 ));
    Sjj_hat = estep_out.source_second_moments;   % {F×1}, n×n

    %% 6) Preprocessing / Whitening in source domain (Module 1)
    pre = module1_preproc_from_covset(Sjj_hat, struct( ...
        'smoothing_method','moving_average', ...
        'loading_factor',  1e-6, ...
        'min_power',       1e-10, ...
        'verbose',         false ));
    D_src     = pre.D;            % {F×1}, n×n
    Sjj_tilde = pre.Sigma_tilde;  % {F×1}, n×n

    %% 7) Active set（前两轮必更新；之后每 update_act_every 轮更新一次）
    if t <= 2 || isempty(A_masks) || mod(t, em.update_act_every) == 0
        input_data_m3 = struct();
        input_data_m3.whitened_covariances = Sjj_tilde;
        input_data_m3.frequencies          = 1:F;
        act = module3_active_set(input_data_m3, struct( ...
                'proxy_method','correlation', 'quantile_level',0.25, ...
                'force_diagonal_active', true, 'verbose', false));
        A_masks = arrayfun(@(f) logical(act.combined_active_mask(:,:,f)), 1:F, 'uni', 0);
        fprintf('[EM] Active set *updated*.\n');
    else
        fprintf('[EM] Active set kept.\n');
    end

    %% 8) Hyperparameters（每轮重估；保持你当前 Source 版的使用方式）
    input_data_m6 = struct();
    input_data_m6.whitened_covariances = Sjj_tilde;
    input_data_m6.kernel_matrix        = K;
    input_data_m6.weight_matrix        = W;
    input_data_m6.active_set_mask      = {A_masks{:}};
    hp = module6_hyperparameter_config(input_data_m6, struct('use_gershgorin', true));

    lambda1 = hp.lambda1;
    lambda2_suggested = hp.lambda2_suggested;
    alpha   = hp.alpha;

    % 你当前脚本中“生效值”仍然是硬编码（这里保持不变）
    lambda2_effective = lambda2_suggested; %#ok<NASGU> % 仅记录，不改变下方 params5 的设定

    %% 9) M-step（Module 5，白化域）
    % —— 初始化：首轮保留你的 trick（用 Omega_true）；其后轮用 transport_init 搬运上一轮 Ω_prev → Γ_init
    if t == 1
        initial_precision_for_mstep = Omega_true;   % 源域真值（仅用于可视化对照的 trick）
    else
        initial_precision_for_mstep = transport_init(Omega_prev, D_src, Sjj_tilde);
    end

    input_data_m5 = struct();
    input_data_m5.whitened_covariances = Sjj_tilde;
    input_data_m5.initial_precision    = initial_precision_for_mstep;
    input_data_m5.smoothing_kernel     = K;
    input_data_m5.weight_matrix        = W;
    input_data_m5.active_set_mask      = {A_masks{:}};
    input_data_m5.whitening_matrices   = D_src;      % 为 Live 可视化白化 GT

    % —— 空间平滑（按你当前 Source 版的建议：λ3 = 0.3*λ1；但 params5 仍保持原硬编码 λ1/λ2，不改）
    lambda3 = 0.3 * lambda1;

    % —— Prox 参数：完全保留你当前 Source 版设定（硬编码生效值）
    params5 = struct( ...
        'lambda1', 1e-3, ...
        'lambda2', 0.02, ...
        'lambda2_suggested', lambda2_suggested, ...  % 仅记录来源
        'alpha0',  5e-3, ...
        'max_iter', 30, 'verbose', true, ...
        'active_set_update_freq', 10, ...
        'alpha_max', 2e-3, 'alpha_up', 1.05, ...
        'alpha_down', 0.7, 'alpha_grow_patience', 2, ...
        'obj_improve_tol', 1e-6, ...
        'weight_mode','hadamard', 'use_graph_laplacian', true, ...
        'diag', struct( 'enable', true, 'update_every', 1, 'metrics_every', 1, ...
                        'print_every', 1, 'f_view', 1, 'log_csv', 'prox_trace.csv', ...
                        'keep_fig_open', true, 'weight_mode','hadamard', ...
                        'use_graph_laplacian', true ) );

    % —— 在线可视化（Live）：保留你的 Source 版 ground truth 对照
    params5.diag.live_plot = live_plot_cfg;

    if lambda3 > 0
        params5.lambda3                    = lambda3;
        params5.spatial_graph_matrix       = Lsp;
        params5.spatial_graph_is_laplacian = true;
        if ~isfield(params5,'spatial_weight_mode') || isempty(params5.spatial_weight_mode)
            params5.spatial_weight_mode = 'node';      % 'node' or 'hadamard'
        end
    end

    [Gamma_tilde_star, prox_res] = module5_proximal(input_data_m5, params5);

    %% 10) Recolor to source scale (Module 8)
    input_data_m8 = struct();
    input_data_m8.whitened_precision_matrices = Gamma_tilde_star;
    input_data_m8.whitening_matrices          = D_src;
    recol = module8_recoloring(input_data_m8, struct());
    Omega_src = recol.recolored_precision_matrices;

    %% 11) 外层收敛统计（保持简单的 Fro 比例变化）
    [delta_Omega, delta_S] = compute_deltas(Omega_src, Omega_prev, Sjj_hat, Sjj_prev);
    fprintf('[EM] ΔΩ=%.3e, ΔS=%.3e\n', delta_Omega, delta_S);

    if delta_Omega < em.tol_Omega && delta_S < em.tol_S
        hold_counter = hold_counter + 1;
    else
        hold_counter = 0;
    end
    if hold_counter >= em.min_hold_iters
        fprintf('>>> EM converged by (ΔΩ,ΔS) for %d consecutive iters at t=%d.\n', hold_counter, t);
        break;
    end

    %% 12) 回写先验进入下一轮（保持稳健逆 + 小地板）
    estep_in.source_prior_covariances = invert_and_fix(Omega_src, 1e-10);

    %% 13) 缓存进入下一轮
    Omega_prev = Omega_src;
    Sjj_prev   = Sjj_hat;
end

fprintf('===== EM finished at iter %d =====\n', t);

%% 14) 评估/可视化（最终一次）—— 保持你的 Source 版写法
f_view = 1;
Om     = Omega_src{f_view};
pcorr  = abs(-Om) ./ sqrt((abs(diag(Om))+eps) * (abs(diag(Om))+eps)');
pcorr(1:n+1:end) = 0;
fprintf('Done. Example partial coherence at f=%d: max=%g, median=%g\n', ...
    f_view, max(pcorr(:)), median(pcorr(pcorr>0)));

viz_pipeline_summary(prox_res, Gamma_tilde_star, Sjj_tilde, K, A_masks, Omega_src, Omega_true);
opt_compare = struct('mode','match_sparsity', 'f_view', 1, 'show_pr', true, 'show_roc', true);
opt_compare.Gamma = Gamma_tilde_star;   % 保持你原来的设定（如需源域评估，可改为 Omega_src）
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

function [Lnorm, info] = normalize_graph_laplacian(L, mode)
% 光谱归一化：把谱搬到 [0,2] 或单位谱半径
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
% 若上一轮为空，则用 “稳健逆 S̃” 作为启动
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
            Gamma_init{f} = G;
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
% 稳健逆：Hermitian→地板→特征逆→Hermitian&实对角
    if nargin < 2, eps_ld = 1e-10; end
    F = numel(Omega_cell);
    n = size(Omega_cell{1},1);
    Sigma_prior = cell(F,1);
    for f = 1:F
        Om = Omega_cell{f};
        Om  = (Om + Om')/2;
        Om(~isfinite(Om)) = 0;
        d   = real(diag(Om)); d = max(d, eps_ld);
        Om(1:n+1:end) = d;
        [U, S] = eig(full(Om), 'vector');
        S  = real(S); S = max(S, eps_ld);
        Sigma = U * diag(1./S) * U';
        Sigma = (Sigma + Sigma')/2;
        Sigma(1:n+1:end) = real(diag(Sigma));
        Sigma_prior{f} = Sigma;
    end
end

function [dOmega, dS] = compute_deltas(Omega, Omega_prev, Sjj, Sjj_prev)
% 外层相对变化（Fro）
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
