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

%% === [ADD-1] sSSBL风格：尺度统一（Leadfield + Data） ===
% 1) 统一 leadfield 尺度
scaleLvj = sqrt(trace(L*L')/size(L,1));
if isfinite(scaleLvj) && scaleLvj>0
    L = L / scaleLvj;
    estep_in.leadfield_matrix = L;
    fprintf('[SCALE] Leadfield scaled by 1/%.3g\n', scaleLvj);
end

% 2) 用归一化后的L估计源功率尺度，统一 Svv 的量纲
gamma_grid = logspace(-4, 1, 30);
opts_el    = struct('maxit', 50, 'tol', 1e-6, 'verbose', false);
warm_tmp   = eloreta_warmstart_from_covs(emp_cov_cell, L, gamma_grid, opts_el);
scaleJ = zeros(F,1);
for f=1:F
    Sjj_e = warm_tmp.Sjj_e{f};
    scaleJ(f) = mean(abs(diag(Sjj_e)));
    if ~isfinite(scaleJ(f)) || scaleJ(f)<=0, scaleJ(f)=1; end
    emp_cov_cell{f} = emp_cov_cell{f} / scaleJ(f);
end
estep_in.empirical_covariances = emp_cov_cell;
fprintf('[SCALE] Data scaled per-freq by 1/mean(diag(Sjj_e)), median factor=%.3g\n', median(scaleJ));

%% 1.5) Warm start：eLORETA 或 sSSBLpp
warm_method = 'ssblpp';   % 'ssblpp' | 'eloreta'
switch lower(warm_method)
    case 'ssblpp'
        opts_ws = struct(); opts_ws.m = T;
        warm = warmstart_ssblpp_from_covs(emp_cov_cell, L, opts_ws);
    case 'eloreta'
        warm = warm_tmp;  % 复用前面已计算的
    otherwise
        error('Unknown warm_method.');
end
Omega_prev = warm.Omega_init; 
Sjj_prev   = warm.Sjj_e;

%% 2) Fixed resources (K, W, spatial L)
K = make_frequency_kernel(F, 3.0);
K = real(0.5*(K + K')); K = max(K, 0);
row_sums_before = sum(K,2); maxrow_before = max(row_sums_before);
if maxrow_before > 0, K = K / maxrow_before; end
rho_after = max(abs(eig(K)));
fprintf('[DIAG] K normalized: maxrow %.3f -> %.3f, rho=%.3f\n', maxrow_before, max(sum(K,2)), rho_after);

W = make_uniform_weight(n);
if isfield(sim, 'source_graph') && isfield(sim.source_graph, 'L'), L_raw = sim.source_graph.L;
else, L_raw = laplacian_placeholder(n); end
[Lsp, Linfo] = normalize_graph_laplacian(L_raw, 'spectral');
fprintf('Spatial L normalized: eig[min,max]=[%.2e, %.2e] -> [%.2e, %.2e]\n', ...
    Linfo.min_eig_before, Linfo.max_eig_before, Linfo.min_eig_after, Linfo.max_eig_after);

%% === [ADD-2] VARETA 子空间（预条件/降维） ===
[U,Sv,V] = svd(L, 'econ');
sing2 = diag(Sv).^2; cum = cumsum(sing2)/sum(sing2);
r = find(cum>=0.99, 1, 'first'); if isempty(r), r = size(Sv,1); end
Ur = U(:,1:r); Sr = Sv(1:r,1:r); Vr = V(:,1:r);
fprintf('[SUBSPACE] rank r=%d / n=%d (%.1f%% energy)\n', r, size(L,2), 100*cum(r));

%% 3) EM loop settings
em = struct();
em.max_em_iter      = 10;
em.tol_Omega        = 5e-3;
em.tol_S            = 1e-3;
em.min_hold_iters   = 2;
em.update_act_every = 2;
em.recfg_hp_every   = 1;
em.lambda2_anneal   = struct('enable', false, 'grow', 1.02, 'max_scale', 1.3);
em.lambda3_ratio    = 0.3;
em.inner = struct( ...
    'alpha0',  0.1, 'alpha_max', 5.0, 'alpha_up', 1.2, 'alpha_down', 0.6, ...
    'alpha_grow_patience', 1, 'max_iter', 30, 'obj_improve_tol', 5e-6, ...
    'weight_mode','hadamard', 'use_graph_laplacian', true, ...
    'diag', struct('enable', true, 'update_every', 1, 'metrics_every', 1, ...
                   'print_every', 1, 'f_view', 1, 'log_csv', 'prox_trace.csv', ...
                   'keep_fig_open', true));

%% 4) State for EM
A_masks = []; if ~exist('Omega_prev','var')||isempty(Omega_prev), Omega_prev=[]; end
if ~exist('Sjj_prev','var')||isempty(Sjj_prev), Sjj_prev=[]; end
hold_counter = 0; hp_last = [];
live_plot_cfg = struct('enable', true, 'f_view', 1, 'plot_every', 5, ...
    'value_mode','abs', 'ground_truth_domain','source');
live_plot_cfg.ground_truth_precision = Omega_true;

fprintf('\n===== EM starts (max %d iters) =====\n', em.max_em_iter);
for t = 1:em.max_em_iter
    fprintf('\n--- EM iter %d ---\n', t);

    %% E-step (Module 2)
    estep_out = module2_estep(estep_in, struct( ...
        'ensure_hermitian', true, 'ensure_real_diag', true, ...
        'ensure_psd', true, 'psd_tol', 1e-10, 'diag_loading', 1e-10 ));
    Sjj_hat = estep_out.source_second_moments;

    %% ===== Noise M-step（带子空间）+ 记录节点方差 =====
    I_p = eye(p); S_xixi_accum = zeros(p,p,'like',estep_in.noise_covariance);
    node_var_cell = cell(F,1);
    for f = 1:F
        Sigma_prior_f = estep_in.source_prior_covariances{f};
        Sigma_prior_f = 0.5*(Sigma_prior_f + Sigma_prior_f');
        Omega_prior_f = inv_psd_robust_(Sigma_prior_f, 1e-8, 1e-12);

        % 子空间版
        SigInvUr = (estep_in.noise_covariance \ Ur);
        M        = Ur' * SigInvUr;
        A_f      = Omega_prior_f + Vr * (Sr * (M * Sr)) * Vr';
        A_f      = (A_f + A_f')/2;
        Sigma_post_f = inv_psd_robust_(A_f, 1e-8, 1e-12);
        node_var_cell{f} = real(diag(Sigma_post_f));

        LinvX = L' / estep_in.noise_covariance;
        T_jv_f  = Sigma_post_f * LinvX;
        T_xi_v_f = I_p - L * T_jv_f;
        S_res_f = T_xi_v_f * emp_cov_cell{f} * T_xi_v_f' + L * Sigma_post_f * L';
        S_res_f = 0.5*(S_res_f + S_res_f');
        S_xixi_accum = S_xixi_accum + S_res_f;
    end
    Sigma_xixi_new = S_xixi_accum / F;
    Sigma_xixi_new = 0.5*(Sigma_xixi_new + Sigma_xixi_new'); 
    ridge = 1e-10 * trace(Sigma_xixi_new)/p;
    Sigma_xixi_new = Sigma_xixi_new + ridge * eye(p,'like',Sigma_xixi_new);
    eta = 0.3; Sigma_xixi = (1-eta)*estep_in.noise_covariance + eta*Sigma_xixi_new;
    [~,chol_flag] = chol(0.5*(Sigma_xixi + Sigma_xixi'),'lower');
    if chol_flag ~= 0, Sigma_xixi = Sigma_xixi + 1e-6*eye(p,'like',Sigma_xixi); end
    estep_in.noise_covariance = Sigma_xixi;
    fprintf('[NOISE] trace(Sigma_xixi)=%.3g\n', trace(Sigma_xixi));

    %% Whitening
    pre = module1_preproc_from_covset(Sjj_hat, struct( ...
        'smoothing_method','moving_average', 'loading_factor',1e-6, 'min_power',1e-10, 'verbose',false ));
    D_src     = pre.D;
    Sjj_tilde = pre.Sigma_tilde;
    diag_whitening_sanity(Sjj_hat, Sjj_tilde, D_src, t);

    %% Active set
    if t <= 2 || isempty(A_masks) || mod(t, em.update_act_every) == 0
        input_data_m3 = struct();
        input_data_m3.whitened_covariances = Sjj_tilde;
        input_data_m3.frequencies          = 1:F;
        act = module3_active_set(input_data_m3, struct( ...
            'proxy_method','correlation', 'quantile_level',0.10, ...
            'force_diagonal_active', true, 'verbose', false));
        A_masks = arrayfun(@(f) logical(act.combined_active_mask(:,:,f)), 1:F, 'uni', 0);
        act_mark = '*updated*';
    else, act_mark = 'kept';
    end
    fprintf('[EM] Active set %s.\n', act_mark);

    %% Hyperparameters
    if mod(t, em.recfg_hp_every) == 0 || isempty(hp_last)
        input_data_m6 = struct();
        input_data_m6.whitened_covariances = Sjj_tilde;
        input_data_m6.kernel_matrix        = K;
        input_data_m6.weight_matrix        = W;
        input_data_m6.active_set_mask      = {A_masks{:}};
        hp = module6_hyperparameter_config(input_data_m6, struct('use_gershgorin', true));
        hp_last = hp;
    else, hp = hp_last; end
    lambda1 = hp.lambda1;
    lambda2_effective = hp.lambda2_suggested;
    lambda3 = em.lambda3_ratio * lambda1; alpha0  = hp.alpha;
    fprintf('[HP ] λ1=%.3g, λ2(sugg)=%.3g, λ3=%.3g, α=%.3g\n', lambda1, lambda2_effective, lambda3, alpha0);

    %% M-step（Module 5，白化域）
    Gamma_init = transport_init(Omega_prev, D_src, Sjj_tilde);

    input_data_m5 = struct();
    input_data_m5.whitened_covariances = Sjj_tilde;
    input_data_m5.initial_precision    = Gamma_init;
    input_data_m5.smoothing_kernel     = K;
    input_data_m5.weight_matrix        = W;
    input_data_m5.active_set_mask      = {A_masks{:}};
    input_data_m5.whitening_matrices   = D_src;

    % Lipschitz 估计
    [L_byS, L_byG] = estimate_L_candidates(Gamma_init, Sjj_tilde);
    L_est     = min([L_byS, L_byG]); alpha_rec = 0.9 / max(L_est, 1e-12);
    alpha0_used = min([hp.alpha, 0.5*alpha_rec, em.inner.alpha_max]); 
    alpha0_used = max(alpha0_used, 5e-4);
    em.inner.alpha0 = alpha0_used;
    fprintf('[STEP] alpha0(rec)=%.3g (L≈%.3g); alpha0(used)=%.3g, alpha_max=%.2f\n', ...
            alpha_rec, L_est, em.inner.alpha0, em.inner.alpha_max);

    % === [ADD-3] 自适应 L1 权重 ===
    l1w_cell = cell(F,1);
    for f=1:F
        wi = 1 ./ (node_var_cell{f} + 1e-8);  wi = wi / median(wi);
        if isempty(Omega_prev), Om_for_deg = warm.Omega_init{f}; else, Om_for_deg = Omega_prev{f}; end
        sup = triu(abs(Om_for_deg) > 1e-8, 1);
        deg = sum(sup | sup', 2);
        Dg  = 1 ./ sqrt((deg + 1) * (deg + 1).');
        Wnode = sqrt(wi * wi.');
        Wpen  = Wnode .* Dg; Wpen(1:n+1:end)=0;
        l1w_cell{f} = Wpen;
    end

    params5 = struct( ...
        'lambda1', lambda1, 'lambda2', lambda2_effective, ...
        'lambda2_suggested', hp.lambda2_suggested, ...
        'alpha0', em.inner.alpha0, 'max_iter', em.inner.max_iter, 'verbose', true, ...
        'active_set_update_freq', 10, 'alpha_max', em.inner.alpha_max, ...
        'alpha_up', em.inner.alpha_up, 'alpha_down', em.inner.alpha_down, ...
        'alpha_grow_patience', em.inner.alpha_grow_patience, ...
        'obj_improve_tol', em.inner.obj_improve_tol, ...
        'weight_mode', em.inner.weight_mode, ...
        'use_graph_laplacian', em.inner.use_graph_laplacian, ...
        'diag', em.inner.diag, ...
        'alpha_min', 1e-5, 'armijo_c1', 1e-5, 'backtrack_beta', 0.5, ...
        'max_backtrack_per_iter', 25, ...
        'backtrack_patience', Inf, 'lambda2_decay_factor', 1.0, ...
        'lambda2_min', lambda2_effective, ...
        'penalize_diagonal', false, ...
        'l1_weights', {l1w_cell}, ...
        'use_single_step', true );   % <<< 新增参数

    params5.diag.live_plot = live_plot_cfg;

    if lambda3 > 0
        params5.lambda3                    = lambda3;
        params5.spatial_graph_matrix       = Lsp;
        params5.spatial_graph_is_laplacian = true;
        params5.spatial_weight_mode        = 'node';
    end

    diag_alpha_sanity(Gamma_init, Sjj_tilde, hp, params5, t);

    [Gamma_tilde_star, prox_res] = module5_proximal(input_data_m5, params5);

    %% Recolor
    input_data_m8 = [];
    input_data_m8.whitened_precision_matrices = Gamma_tilde_star;
    input_data_m8.whitening_matrices = D_src;
    recol = module8_recoloring(input_data_m8, struct());
    Omega_src = recol.recolored_precision_matrices;

    %% 外层收敛统计
    [delta_Omega, delta_S] = compute_deltas(Omega_src, Omega_prev, Sjj_hat, Sjj_prev);
    gnorm = read_grad_norm_(prox_res);
    fprintf('[EM] ΔΩ=%.3e, ΔS=%.3e, ||grad||≈%.3g\n', delta_Omega, delta_S, gnorm);

    if delta_Omega < em.tol_Omega && delta_S < em.tol_S, hold_counter = hold_counter + 1; 
    else, hold_counter = 0; end
    if hold_counter >= em.min_hold_iters
        fprintf('>>> EM converged by (ΔΩ,ΔS) for %d consecutive iters at t=%d.\n', hold_counter, t);
        break;
    end

    estep_in.source_prior_covariances = invert_and_fix(Omega_src, 1e-10);
    Omega_prev = Omega_src; Sjj_prev = Sjj_hat;
end
fprintf('===== EM finished at iter %d =====\n', t);

%% === [ADD-4] Support-refit 去偏 ===
support_cell = cell(F,1);
for f=1:F
    Gf = Gamma_tilde_star{f};
    M  = abs(Gf) > 1e-8; M(1:n+1:end) = true;
    support_cell{f} = M;
end
Gamma_init_refit = cell(F,1);
for f=1:F
    G0 = Gamma_tilde_star{f};
    G0(~support_cell{f}) = 0;
    Gamma_init_refit{f} = (G0+G0')/2;
end
input_data_m5_refit = input_data_m5;
input_data_m5_refit.initial_precision = Gamma_init_refit;
params_refit = params5; 
params_refit.max_iter = 10; 
params_refit.lambda2 = 1e-8;
params_refit.l1_weights = []; 
params_refit.active_set_update_freq = Inf;
params_refit.use_single_step = true;  % <<< 保持一致
input_data_m5_refit.active_set_mask = support_cell;
[Gamma_tilde_refit, ~] = module5_proximal(input_data_m5_refit, params_refit);
input_data_m8 = []; 
input_data_m8.whitened_precision_matrices = Gamma_tilde_refit;
input_data_m8.whitening_matrices = D_src;
recol2 = module8_recoloring(input_data_m8, struct());
Omega_src = recol2.recolored_precision_matrices;

%% 5) 评估/可视化
f_view = 1;
Om     = Omega_src{f_view};
pcorr  = abs(-Om) ./ sqrt((abs(diag(Om))+eps) * (abs(diag(Om))+eps)');
pcorr(1:n+1:end) = 0;
fprintf('Done. Example partial coherence at f=%d: max=%g, median=%g\n', ...
    f_view, max(pcorr(:)), median(pcorr(pcorr>0)));

viz_pipeline_summary(prox_res, Gamma_tilde_star, Sjj_tilde, K, A_masks, Omega_src, Omega_true);

opt_compare = struct('mode','match_sparsity', 'f_view', 1, 'show_pr', true, 'show_roc', true);
opt_compare.Gamma = Omega_src;
metrics = gt_compare_and_plot(Omega_true, Omega_src, opt_compare);

opts_sc = struct('zero_tol',1e-12, 'use_abs',true, 'per_freq',true);
sc1 = sanity_edge_ranking(Omega_src, Omega_true, opts_sc, 'raw_abs');
sc2 = sanity_edge_ranking(Omega_src, Omega_true, opts_sc, 'pcorr_abs');
fprintf('\n[SELF-CHECK] ===== Summary =====\n');
fprintf('RAW |Ω|     : AUROC=%.3f (vs -score: %.3f) | AUPR=%.3f | ρ(|Ω_est|,|Ω_true|)=%.3f\n', ...
    sc1.auroc_pos, sc1.auroc_neg, sc1.aupr, sc1.corr_pos);
fprintf('PCorr |Ω|   : AUROC=%.3f (vs -score: %.3f) | AUPR=%.3f | ρ(|Ω_est|,|Ω_true|)=%.3f\n', ...
    sc2.auroc_pos, sc2.auroc_neg, sc2.aupr, sc2.corr_pos);
if sc1.auroc_neg > sc1.auroc_pos + 0.03 || sc2.auroc_neg > sc2.auroc_pos + 0.03
    fprintf('>>> [ALERT] 排序方向疑似反了。\n');
end
if isfield(sc1,'perF') && ~isempty(sc1.perF)
    fprintf('[Per-freq AUROC] raw |Ω|:'); fprintf(' %.3f', sc1.perF); fprintf('\n');
    fprintf('[Per-freq AUROC] pcorr  :'); fprintf(' %.3f', sc2.perF); fprintf('\n');
end
[pe,se] = edge_corr_abs(Omega_src, Omega_true);
fprintf('[Density] est>0 ratio=%.3f | CorrΩ : Pearson=%.3f | Spearman=%.3f\n', ...
    sc1.est_density, pe, se);

%% ===== Local helpers (保持不变) =====
function C = coerce_cov_cell(X, F_hint)
    if isa(X,'cell'), C = X(:); return; end
    if isnumeric(X) && ndims(X)==3 && size(X,1)==size(X,2)
        F = size(X,3); C = cell(F,1); for f = 1:F, C{f} = X(:,:,f); end, return;
    end
    if isnumeric(X) && ismatrix(X) && size(X,1)==size(X,2)
        if nargin>=2 && ~isempty(F_hint), C = repmat({X}, F_hint, 1); else, C = {X}; end, return;
    end
    error('coerce_cov_cell:unsupported','Expect cell{F,1} | p×p×F | single p×p.');
end
function K = make_frequency_kernel(F, sigma)
    if nargin < 2, sigma = 3.0; end
    [I,J] = ndgrid(1:F,1:F);
    K = exp(-((I-J).^2)/(2*sigma^2)); K = (K + K')/2;
end
function W = make_uniform_weight(n), W = ones(n); W(1:n+1:end) = 0; end
function L = laplacian_placeholder(n), A = ones(n) - eye(n); d = sum(A,2); L = diag(d) - A; end
function [Lnorm, info] = normalize_graph_laplacian(L, mode)
    if nargin<2, mode='spectral'; end
    L = (L+L')/2; ev = eig(full(L));
    info.min_eig_before = min(real(ev)); info.max_eig_before = max(real(ev));
    switch lower(mode), case 'spectral', s = max(1, info.max_eig_before); Lnorm = L / s; otherwise, Lnorm = L; end
    ev2 = eig(full(Lnorm)); info.min_eig_after = min(real(ev2)); info.max_eig_after = max(real(ev2));
end
function Gamma_init = transport_init(Omega_prev, D_src, Sjj_tilde)
    F = numel(D_src); Gamma_init = cell(F,1);
    if isempty(Omega_prev)
        for f=1:F
            St = (Sjj_tilde{f} + Sjj_tilde{f}')/2;
            [U, D] = eig(full(St), 'vector'); d = real(D); d = max(d, 1e-10);
            G = U * diag(1./d) * U'; G = (G + G')/2;
            Gamma_init{f} = G;
        end, return;
    end
    for f=1:F
        D = D_src{f}; Gamma_init{f} = (D \ Omega_prev{f}) / D; Gamma_init{f} = (Gamma_init{f} + Gamma_init{f}')/2;
    end
end
function Sigma_prior = invert_and_fix(Omega_cell, eps_ld)
    if nargin < 2, eps_ld = 1e-10; end
    F = numel(Omega_cell); n = size(Omega_cell{1},1); Sigma_prior = cell(F,1);
    for f=1:F
        Om  = (Omega_cell{f} + Omega_cell{f}')/2; Om(~isfinite(Om)) = 0;
        d   = real(diag(Om)); d = max(d, eps_ld); Om(1:n+1:end) = d;
        [U, S] = eig(full(Om), 'vector'); S = real(S); S = max(S, 2*eps_ld);
        Sigma =U * diag(1./S) * U'; Sigma = (Sigma + Sigma')/2; Sigma(1:n+1:end) = real(diag(Sigma));
        Sigma_prior{f} = Sigma;
    end
end
function [dOmega, dS] = compute_deltas(Omega, Omega_prev, Sjj, Sjj_prev)
    if isempty(Omega_prev), dOmega = inf; else
        num=0; den=0; for f=1:numel(Omega), num=num+norm(Omega{f}-Omega_prev{f},'fro'); den=den+norm(Omega_prev{f},'fro'); end
        dOmega = num / max(1, den);
    end
    if isempty(Sjj_prev), dS = inf; else
        num=0; den=0; for f=1:numel(Sjj), num=num+norm(Sjj{f}-Sjj_prev{f},'fro'); den=den+norm(Sjj_prev{f},'fro'); end
        dS = num / max(1, den);
    end
end
function g = read_grad_norm_(prox_res)
    g = NaN;
    try
        if isfield(prox_res,'grad_norm'), g = prox_res.grad_norm; end
        if isnan(g) && isfield(prox_res,'gradient_norm_history'), h = prox_res.gradient_norm_history; if ~isempty(h), g = h(end); end, end
        if isnan(g) && isfield(prox_res,'grad_norm_mean'), g = prox_res.grad_norm_mean; end
        if isnan(g), g = -1; end
    catch, g = -1; end
end
function [pe, se] = edge_corr_abs(Acell, Bcell)
    F = numel(Acell); ps = zeros(F,1); ss = zeros(F,1);
    for f=1:F
        A = abs(Acell{f}); B = abs(Bcell{f}); idx = triu(true(size(A)),1);
        a = A(idx); b = B(idx);
        if isempty(a) || isempty(b) || all(a==a(1)) || all(b==b(1)), ps(f)=0; ss(f)=0;
        else, ps(f) = corr(a,b,'type','Pearson','rows','complete'); ss(f) = corr(a,b,'type','Spearman','rows','complete'); end
    end
    pe = mean(ps); se = mean(ss);
end

%% ========= Self-check helpers =========
function diag_whitening_sanity(Sjj_hat, Sjj_tilde, D_src, iter_id)
    F = numel(Sjj_hat);
    md_hat = zeros(F,1); sd_hat = zeros(F,1); md_til = zeros(F,1); sd_til = zeros(F,1);
    rel_err = zeros(F,1); is_spd  = true;
    for f = 1:F
        Sh = (Sjj_hat{f} + Sjj_hat{f}')/2;
        St = (Sjj_tilde{f} + Sjj_tilde{f}')/2;
        D  = D_src{f};
        md_hat(f)=mean(real(diag(Sh))); sd_hat(f)=std(real(diag(Sh)));
        md_til(f)=mean(real(diag(St))); sd_til(f)=std(real(diag(St)));
        R = D*Sh*D - St; rel_err(f) = norm(R,'fro') / max(1, norm(St,'fro'));
        try, ev = eig(St); is_spd = is_spd && all(real(ev) > -1e-10); catch, is_spd=false; end
    end
    fprintf(['[DIAG][t=%d] <Whitening>\n    mean(diag S_hat)=%.3g±%.3g | mean(diag S_tilde)=%.3g±%.3g (≈1)\n' ...
             '    max relErr ||D*S*D - S_tilde||_F / ||S_tilde||_F = %.2e | S_tilde SPD? %s\n'], ...
            iter_id, mean(md_hat), mean(sd_hat), mean(md_til), mean(sd_til), max(rel_err), ternary(is_spd,'YES','NO'));
    if mean(md_til) > 2 || mean(md_til) < 0.5
        warning('[DIAG] whitened diagonals far from 1 (mean=%.3g).', mean(md_til));
    end
    if max(rel_err) > 1e-6
        warning('[DIAG] D*S*D and S_tilde mismatch (%.2e).', max(rel_err));
    end
end
function diag_alpha_sanity(Gamma_init, Sjj_tilde, hp, params5, iter_id)
    [L_byS, L_byG] = estimate_L_candidates(Gamma_init, Sjj_tilde);
    L_est = min([L_byS, L_byG]); alpha_rec = 0.9 / max(L_est, 1e-12);
    amax = []; if isfield(params5,'alpha_max') && ~isempty(params5.alpha_max), amax = params5.alpha_max; end
    fprintf('[DIAG][t=%d] <Step-size> L_byS≈%.3g, L_byG≈%.3g -> α_rec≈%.3g | hp.α=%.3g, α_max=%s\n', ...
            iter_id, L_byS, L_byG, alpha_rec, hp.alpha, iif(isempty(amax),'[]',num2str(amax)));
    if ~isempty(amax) && amax < 0.1*alpha_rec
        warning('[DIAG] alpha_max (%.3g) << α_rec (%.3g).', amax, alpha_rec);
    end
end
function [L_byS, L_byG] = estimate_L_candidates(Gamma_init, Sjj_tilde)
    F = numel(Sjj_tilde); Ls = 0; Lg = 0;
    for f = 1:F
        St = (Sjj_tilde{f} + Sjj_tilde{f}')/2; s  = svds(St, 1); Ls = max(Ls, s^2);
        G  = (Gamma_init{f} + Gamma_init{f}')/2;
        try, lam_min = min(real(eig(G))); if ~isfinite(lam_min) || lam_min <= 0, invnorm = svds(pinv(G), 1); else, invnorm = 1/lam_min; end
        catch, invnorm = svds(pinv(G), 1); end
        Lg = max(Lg, invnorm^2);
    end, L_byS = Ls; L_byG = Lg;
end
function out = ternary(cond, a, b), if cond, out=a; else, out=b; end, end
function out = iif(cond, a, b), if cond, out=a; else, out=b; end, end

%% ========= sanity_edge_ranking helper =========
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
        y = abs(t) > ztol;
        sgt = abs(t);

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

    [~,~,~,auroc_pos] = perfcurve(y_all,  s_all,  true);
    [~,~,~,auroc_neg] = perfcurve(y_all, -s_all, true);
    [rec,prec,~,aupr] = perfcurve(y_all, s_all, true, 'xCrit','reca','yCrit','prec');

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

    est_density = mean(s_all > max(ztol, 1e-9));
    pos_ratio   = mean(y_all==1);
    sc.est_density = est_density; sc.pos_ratio = pos_ratio;
end
function v = getfield_with_default(s, name, def)
    if isfield(s,name) && ~isempty(s.(name)), v = s.(name); else, v = def; end
end

%% ========= eLORETA helpers =========
function warm = eloreta_warmstart_from_covs(Svv_cell, L, gamma_grid, opts)
    if nargin < 4, opts = struct(); end
    maxit = getf(opts,'maxit',50); tol = getf(opts,'tol',1e-6); verb = getf(opts,'verbose',false);
    F = numel(Svv_cell);
    Tjv_cell = cell(F,1); Sjj_cell = cell(F,1);
    Om0_cell = cell(F,1); gcv_cell = cell(F,1); gopt_cell = cell(F,1);
    for f = 1:F
        Svv = Svv_cell{f};
        [Tjv, Sjj, ~, gamma_opt, gcv] = eloreta_simple_(Svv, L, gamma_grid, maxit, tol, verb);
        Sjj = psd_project_(Sjj); Om0 = inv_psd_robust_(Sjj, 1e-8, 1e-12);
        Tjv_cell{f}  = Tjv; Sjj_cell{f}  = Sjj; Om0_cell{f}  = Om0;
        gcv_cell{f}  = gcv; gopt_cell{f} = gamma_opt;
    end
    warm = struct('Omega_init',{Om0_cell}, 'Sjj_e',{Sjj_cell}, 'Tjv',{Tjv_cell}, ...
                  'gamma_opt',{gopt_cell}, 'gcv_curve',{gcv_cell});
end
function [Tjv, Sjj, W, gamma_opt, gcv] = eloreta_simple_(Svv, L, gamma_grid, maxit, tol, verb)
    [p,n] = size(L); gcv = zeros(numel(gamma_grid),1);
    best.T = []; best.W = []; best.gamma = NaN; best.score = Inf;
    for k = 1:numel(gamma_grid)
        gamma = gamma_grid(k);
        w = ones(n,1);
        for it=1:maxit
            Winv = diag(1./w); A = hermi_( L*Winv*L' ); alpha= gamma * trace(A)/p;
            M    = inv_psd_( A + alpha*eye(p) );
            w_old = w; for i=1:n, li=L(:,i); mii = real(li' * M * li); w(i)= sqrt(max(mii, eps)); end
            if norm(w-w_old)/max(1,norm(w_old)) < tol, break; end
        end
        Winv = diag(1./w); T = Winv * L' * M;
        Txiv = eye(p) - L*T; num  = real(trace( hermi_(Txiv*Svv*Txiv') ))/p;
        den  = ( real(trace(Txiv))/p )^2 + eps; gcv(k) = num / den;
        if gcv(k) < best.score, best.T=T; best.W=diag(w); best.gamma=gamma; best.score=gcv(k); end
        if verb && (mod(k,10)==1), fprintf('[eLORETA] gamma=%.3g, GCV=%.3g\n', gamma, gcv(k)); end
    end
    Tjv = best.T; gamma_opt = best.gamma; W = best.W; Sjj = Tjv * Svv * Tjv';
end
function A = hermi_(A), A = (A + A')/2; end
function X = inv_psd_(X), X = hermi_(X); [U,S] = eig(X,'vector'); S = max(real(S), eps); X = U*diag(1./S)*U'; X = hermi_(X); end
function S = psd_project_(S), S = hermi_(S); [U,d] = eig(S,'vector'); d = max(real(d), 0); S = U*diag(d)*U'; S = hermi_(S); end
function v = getf(s, f, d), if isfield(s,f) && ~isempty(s.(f)), v=s.(f); else, v=d; end, end
function Om = inv_psd_robust_(A, eps_reg, min_ratio)
    A = (A + A')/2; [V,D] = eig(A); d = real(diag(D)); dmax = max(d);
    floor_val = max(min_ratio * max(dmax, eps), 0);
    d(d<floor_val) = floor_val; if eps_reg > 0, d = (d + eps_reg * dmax) / (1 + eps_reg); end
    Om = V * diag(1./d) * V'; Om = (Om + Om')/2;
end