function [Omega_est, Sigma_src_est, outs] = solver_jspace_3d_opt(Svv_cell, L, GraphLaplacian, cfg)
% SOLVER_JSPACE_3D_OPT - Scheme-1 (eta-search) 3D hyperparameter optimization + stable EM
%
% Revisions:
%   - [OPTIMIZER] Added Smart Bounds Shrinking (eta_shrink) & Warm-Start (InitialPoints).
%   - [FEATURE] Added cfg.postprocess_enable to skip Rayleigh/Debias steps.
%   - [DIAGNOSTICS] Added strict scale/condition checks (dbg_jspace_cells).
%   - [SAFETY] Added automatic fallback if Debias produces invalid matrices.
%
if nargin < 4, cfg = struct(); end

% ============================================================
% 1) Input normalization
% ============================================================
if ~iscell(Svv_cell)
    if ndims(Svv_cell) == 3
        F_in = size(Svv_cell, 3);
        tmp = cell(F_in, 1);
        for f = 1:F_in, tmp{f} = Svv_cell(:, :, f); end
        Svv_cell = tmp;
    else
        Svv_cell = {Svv_cell};
    end
end
F = numel(Svv_cell);
[Ns, Nr] = size(L);

% ============================================================
% 2) Config
% ============================================================
MAX_EM_ITER        = get_cfg(cfg, 'max_em_iter', 10);
VERBOSE            = get_cfg(cfg, 'verbose', true);
USE_GPU            = get_cfg(cfg, 'use_gpu', false);
DEBUG_PRINT        = get_cfg(cfg, 'debug_print', false);
DEBUG_RAY          = get_cfg(cfg, 'debug_ray', true);

% Switch: Enable/Disable Robust Post-Processing
DO_POST            = get_cfg(cfg, 'postprocess_enable', true); 

% Surrogate settings
MAX_OPT_EVALS      = get_cfg(cfg, 'opt_max_evals', 50);
MIN_OPT_POINTS     = get_cfg(cfg, 'opt_min_points', 10);
USE_PARALLEL       = get_cfg(cfg, 'opt_use_parallel', true);

% Inner-loop objective speed
OBJ_MAX_ITER       = get_cfg(cfg, 'obj_max_iter', 30);
OBJ_TOL            = get_cfg(cfg, 'obj_tol', 3e-3);
REFIT_MAX_ITER     = get_cfg(cfg, 'obj_refit_max_iter', 12);
REFIT_TOL          = get_cfg(cfg, 'obj_refit_tol', 2e-2);

% Density definition (partial-corr)
PCOR_EPS           = get_cfg(cfg, 'dens_pcor_eps', 5e-3);
DENS_PEN_W         = get_cfg(cfg, 'density_penalty_weight', 1e7);

% Mask stabilization
MASK_DENS_FLOOR    = get_cfg(cfg, 'mask_density_floor', max(0.03, log(Nr)/Nr));
MASK_UNION_MODE    = get_cfg(cfg, 'mask_union', true);

% EM inertia update
UPDATE_RATE        = get_cfg(cfg, 'update_rate', 0.35);

% Rayleigh Range (Will be overridden by Adaptive Logic)
RAY_RANGE          = get_cfg(cfg, 'rayleigh_range', 2.0:0.2:5.0);

% Samples for information criterion
M_SAMPLES          = get_cfg(cfg, 'm_samples', 100 * Nr);

% Optimizer choice
optimizer = lower(get_cfg(cfg, 'optimizer', 'fista'));
if isfield(cfg, 'use_fista') && cfg.use_fista, optimizer = 'fista'; end
if strcmp(optimizer, 'fista')
    mstep_solver = @module5_fista_main;
else
    mstep_solver = @module5_proximal_main;
end

% Graph Laplacian
if isempty(GraphLaplacian)
    GraphLaplacian = zeros(Nr, Nr, 'like', L);
end

% Frequency kernel
if isfield(cfg, 'freq_kernel') && ~isempty(cfg.freq_kernel)
    K_freq = cfg.freq_kernel;
else
    K_freq = build_default_freq_kernel_(F);
end

% Weight matrix W_gamma (DWI soft prior)
W_gamma = get_cfg(cfg, 'weight_matrix', ones(Nr));
dwi_mask = [];
dwi_connectivity_mask = [];
if isfield(cfg, 'dwi_C') && ~isempty(cfg.dwi_C)
    if VERBOSE
        fprintf('[J-SPACE-3D] Using DWI soft prior: building W_gamma from cfg.dwi_C.\n');
    end
    [W_gamma, dwi_mask] = build_dwi_soft_prior_(cfg.dwi_C, Nr, cfg);
    C_conn = gather(cfg.dwi_C);
    C_conn = real(C_conn);
    C_conn = (C_conn ~= 0);
    C_conn = C_conn | C_conn.';
    C_conn(1:Nr+1:end) = true;
    dwi_connectivity_mask = C_conn;
end

% ============================================================
% 3) Noise covariance
% ============================================================
tr_S = 0;
for f = 1:F, tr_S = tr_S + trace(Svv_cell{f}); end
noise_cov = (tr_S / (F * Ns)) * 0.05 * eye(Ns);

% ============================================================
% 4) GPU init
% ============================================================
if USE_GPU
    try
        L = gpuArray(L);
        noise_cov = gpuArray(noise_cov);
        for f = 1:F, Svv_cell{f} = gpuArray(Svv_cell{f}); end
        GraphLaplacian = gpuArray(GraphLaplacian);
        if VERBOSE, fprintf('[J-SPACE-3D] GPU Acceleration ENABLED (E-step).\n'); end
    catch
        USE_GPU = false;
        if VERBOSE, fprintf('[J-SPACE-3D] GPU init failed. Falling back to CPU.\n'); end
    end
end

% ============================================================
% 5) Initialization (eLORETA)
% ============================================================
Sigma_source_curr = [];
if isfield(cfg, 'init_sigma_source') && ~isempty(cfg.init_sigma_source)
    Sigma_source_curr = cfg.init_sigma_source;
    if ~iscell(Sigma_source_curr), Sigma_source_curr = {Sigma_source_curr}; end
else
    if VERBOSE, fprintf('[Init] Running eLORETA...\n'); end
    S_avg = zeros(Ns, Ns, 'like', Svv_cell{1});
    for f = 1:F, S_avg = S_avg + Svv_cell{f}; end
    S_avg = S_avg / F;
    [T_eloreta, ~] = run_eloreta_core(L, S_avg, 0.05);
    Sigma_source_curr = cell(F, 1);
    for f = 1:F
        Sjj_eloreta = utils_math.make_hermitian(T_eloreta * Svv_cell{f} * T_eloreta');
        [Sjj_spd, ~] = utils_math.project_spd(Sjj_eloreta, 1e-8);
        Sigma_source_curr{f} = Sjj_spd;
    end
end
Gamma_warm_start = cell(F, 1);
for f = 1:F, Gamma_warm_start{f} = eye(Nr, 'like', L); end

% ============================================================
% 6) One E-step + Whitening
% ============================================================
if VERBOSE, fprintf('[J-SPACE-3D] Preparing data for optimization...\n'); end
[Psijj_cell, ~] = module2_estep(Svv_cell, L, Sigma_source_curr, noise_cov);
[Sjj_tilde0, D_cell, ~] = module1_data_whitening(Psijj_cell, 'smoothing_window', 1);

% [DIAGNOSIS POINT A] Check Scale of Input Covariances
dbg_jspace_cells(Sjj_tilde0, [], 'Sjj_tilde BEFORE any optimization');

% Base CPU input
m5_in_base = struct();
m5_in_base.whitened_covariances = cell(F, 1);
for f = 1:F, m5_in_base.whitened_covariances{f} = gather(Sjj_tilde0{f}); end
m5_in_base.smoothing_kernel = gather(K_freq);
m5_in_base.weight_matrix    = gather(W_gamma);

% Base params
m5_p_base = struct();
m5_p_base.alpha0 = 0.1;
m5_p_base.min_eig = get_cfg(cfg, 'min_eig', 1e-6);
m5_p_base.spatial_graph_matrix = gather(GraphLaplacian);
m5_p_base.spatial_graph_is_laplacian = true;
m5_p_base.max_iter = OBJ_MAX_ITER;
m5_p_base.tol = OBJ_TOL;
m5_p_base.verbose = false;
m5_p_base.weight_mode = 'hadamard';
m5_p_base.auto_tune = false;
if strcmp(optimizer, 'fista')
    m5_p_base.backtracking_beta = 0.5;
    m5_p_base.max_backtracking  = 20;
    m5_p_base.monotone          = true;
    m5_p_base.use_restart       = true;
end

% ============================================================
% 7) Compute scales/thresholds
% ============================================================
if VERBOSE, fprintf('[J-SPACE-3D] Computing scales/thresholds for eta-search...\n'); end
[sc0, thr0] = compute_scales_thresholds_(m5_in_base.whitened_covariances, m5_in_base.smoothing_kernel, m5_in_base.weight_matrix, Nr, cfg);

min_density_limit = log(Nr) / Nr;
max_density_limit = thr0.data_density;
max_density_cap   = get_cfg(cfg, 'max_density_cap', 0.40);
max_density_limit = min(max_density_limit, max_density_cap);
max_density_limit = max(max_density_limit, min_density_limit * 1.5);

if VERBOSE
    fprintf('  [Auto] Nr=%d => log(N)/N=%.2f%%\n', Nr, min_density_limit*100);
    fprintf('  [Auto] t_active=%.4f | t_rescue=%.4f | data_density=%.2f%% (clamped max=%.2f%%)\n', ...
        thr0.t_active, thr0.t_rescue, thr0.data_density*100, max_density_cap*100);
end

% Initial active mask
act_params = struct('quantile_level', -thr0.t_active, 'strategy', 'intersection', 'force_diagonal', true);
[mask0, ~] = module3_active_set(m5_in_base.whitened_covariances, act_params);
if ~isempty(dwi_connectivity_mask)
    for f = 1:F
        mask0{f} = mask0{f} & dwi_connectivity_mask;
        mask0{f}(1:Nr+1:end) = true;
    end
end
m5_in_base.active_mask = mask0;

% ============================================================
% 8) Surrogate optimization on eta (WITH shrink + warm-start)
% ============================================================

% ---- 8.1 Read global bounds (keep your original defaults) ----
eta1_lb = get_cfg(cfg, 'eta1_lb', 1e-3);
eta1_ub = get_cfg(cfg, 'eta1_ub', 0.5);
eta2_lb = get_cfg(cfg, 'eta2_lb', 1e-3);
eta2_ub = get_cfg(cfg, 'eta2_ub', 1.0);
eta3_lb = get_cfg(cfg, 'eta3_lb', 1e-3);
eta3_ub = get_cfg(cfg, 'eta3_ub', 0.5);

% (Support vector input if provided)
if isfield(cfg,'eta_lb') && numel(cfg.eta_lb)==3
    eta1_lb = cfg.eta_lb(1); eta2_lb = cfg.eta_lb(2); eta3_lb = cfg.eta_lb(3);
end
if isfield(cfg,'eta_ub') && numel(cfg.eta_ub)==3
    eta1_ub = cfg.eta_ub(1); eta2_ub = cfg.eta_ub(2); eta3_ub = cfg.eta_ub(3);
end

eta_global_lb = [eta1_lb, eta2_lb, eta3_lb];
eta_global_ub = [eta1_ub, eta2_ub, eta3_ub];

% ---- 8.2 Default init eta0 ----
eta0 = [get_cfg(cfg,'eta1_init',0.1), get_cfg(cfg,'eta2_init',0.2), get_cfg(cfg,'eta3_init',0.1)];
if isfield(cfg,'eta0') && numel(cfg.eta0)==3
    eta0 = cfg.eta0(:).';
end
eta0 = max(eta0, eta_global_lb);
eta0 = min(eta0, eta_global_ub);
x0_log = log10(eta0);

% ---- 8.3 Shrink bounds using prev/bank ----
DO_SHRINK  = get_cfg(cfg, 'eta_shrink_enable', true);
SH_DEC     = get_cfg(cfg, 'eta_shrink_decades', 0.7);        % ± in log10
SH_MODE    = lower(get_cfg(cfg, 'eta_shrink_mode', 'prev')); % 'prev'|'bank_range'

eta_prev = [];
if isfield(cfg,'eta_prev') && ~isempty(cfg.eta_prev)
    eta_prev = cfg.eta_prev(:).';
end
eta_bank = [];
if isfield(cfg,'eta_bank') && ~isempty(cfg.eta_bank) && size(cfg.eta_bank,2)==3
    eta_bank = cfg.eta_bank;
end

lb = eta_global_lb;
ub = eta_global_ub;

if DO_SHRINK
    if strcmp(SH_MODE,'bank_range') && ~isempty(eta_bank)
        % Cover the range of historical best points
        lb = min(eta_bank,[],1) ./ (10.^SH_DEC);
        ub = max(eta_bank,[],1) .* (10.^SH_DEC);
    else
        % Shrink around previous point only
        if ~isempty(eta_prev), c = eta_prev; else, c = eta0; end
        lb = c ./ (10.^SH_DEC);
        ub = c .* (10.^SH_DEC);
    end

    % Clamp to global bounds
    lb = max(lb, eta_global_lb);
    ub = min(ub, eta_global_ub);

    % Prevent inversion if shrinkage is too aggressive
    if any(lb >= ub)
        lb = eta_global_lb;
        ub = eta_global_ub;
    end
end

lb_log = log10(lb);
ub_log = log10(ub);

if VERBOSE
    fprintf('[J-SPACE-3D] eta bounds (log10): lb=[%.2f %.2f %.2f], ub=[%.2f %.2f %.2f]\n', ...
        lb_log(1),lb_log(2),lb_log(3), ub_log(1),ub_log(2),ub_log(3));
end

% ---- 8.4 Build InitialPoints (warm-start + a few LHS points) ----
N_INIT    = get_cfg(cfg, 'opt_init_points', 8);   
BANK_K    = get_cfg(cfg, 'opt_init_bank_k', 5);
SEED      = get_cfg(cfg, 'opt_rng_seed', 0);

init_pts = x0_log;  % always include eta0

% Include eta_prev
if ~isempty(eta_prev)
    init_pts = [init_pts; log10(max(min(eta_prev, ub), lb))];
end

% Include last BANK_K from eta_bank
if ~isempty(eta_bank)
    k = min(BANK_K, size(eta_bank,1));
    pick = eta_bank(end-k+1:end, :);
    pick = max(min(pick, ub), lb);
    init_pts = [init_pts; log10(pick)];
end

% Remove out-of-bound & duplicates
inb = all(init_pts >= lb_log & init_pts <= ub_log, 2);
init_pts = init_pts(inb,:);
init_pts = unique(round(init_pts, 6), 'rows', 'stable');

% Fill with LHS until N_INIT
n_missing = max(0, N_INIT - size(init_pts,1));
if n_missing > 0
    rng(SEED, 'twister');
    U = lhsdesign(n_missing, 3, 'criterion', 'maximin', 'iterations', 20);
    lhs_pts = lb_log + U .* (ub_log - lb_log);
    init_pts = [init_pts; lhs_pts];
end

% Cap initial points to avoid consuming entire budget
max_init_safe = min(10, max(3, floor(MAX_OPT_EVALS/5)));
if size(init_pts,1) > max_init_safe
    init_pts = init_pts(1:max_init_safe, :);
end

if VERBOSE
    fprintf('[J-SPACE-3D] surrogateopt InitialPoints = %d (budget=%d)\n', size(init_pts,1), MAX_OPT_EVALS);
end

% ---- 8.5 Run surrogateopt ----
if VERBOSE
    fprintf('[J-SPACE-3D] Running surrogateopt on eta (log-space)...\n');
end

obj_fun = @(x_log) wrapper_objective_eta_aic_density_( ...
    x_log, m5_in_base, m5_p_base, sc0, ...
    M_SAMPLES, min_density_limit, max_density_limit, ...
    PCOR_EPS, DENS_PEN_W, ...
    Nr, F, mstep_solver, ...
    REFIT_MAX_ITER, REFIT_TOL);

% ---- [MOD] allow multi initial points in log-space ----
initPts = x0_log;
if isfield(cfg,'opt_initial_points_log') && ~isempty(cfg.opt_initial_points_log)
    initPts = cfg.opt_initial_points_log;   % K x 3 in log10-space
end

opts = optimoptions('surrogateopt', ...
    'MaxFunctionEvaluations', MAX_OPT_EVALS, ...
    'MinSurrogatePoints', MIN_OPT_POINTS, ...
    'UseParallel', USE_PARALLEL, ...
    'PlotFcn', [], ...
    'InitialPoints', initPts);


tic_opt = tic;
[x_best_log, fval, exitflag, output_opt, trials] = surrogateopt(obj_fun, lb_log, ub_log, opts);
time_opt = toc(tic_opt);
eta_best = 10.^x_best_log;

outs = struct();
outs.opt_results.fval = fval;
outs.opt_results.exitflag = exitflag;
outs.opt_results.output = output_opt;
outs.opt_results.trials = trials;
outs.opt_results.trials.X = 10.^trials.X;
outs.best_eta.eta1 = eta_best(1);
outs.best_eta.eta2 = eta_best(2);
outs.best_eta.eta3 = eta_best(3);

if VERBOSE
    fprintf('  > Opt finished in %.2fs. Best AIC+Penalty: %.3e\n', time_opt, fval);
    fprintf('  > BEST eta: eta1=%.3e | eta2=%.3e | eta3=%.3e\n', eta_best(1), eta_best(2), eta_best(3));
end

outs.thresholds = thr0;
outs.density_limits.min = min_density_limit;
outs.density_limits.max = max_density_limit;

% ============================================================
% 9) Full EM with DYNAMIC lambdas
% ============================================================
if VERBOSE, fprintf('\n[J-SPACE-3D] Running Full EM with DYNAMIC lambdas...\n'); end
m5_p_em = m5_p_base;
m5_p_em.max_iter = get_cfg(cfg, 'em_mstep_max_iter', 80);
m5_p_em.tol      = get_cfg(cfg, 'em_mstep_tol', 1e-3);
outs.loglik = zeros(MAX_EM_ITER, 1);
outs.em_density_trace = zeros(MAX_EM_ITER, 1);
outs.em_lambdas = zeros(MAX_EM_ITER, 3);
active_mask_prev = m5_in_base.active_mask;

for em_iter = 1:MAX_EM_ITER
    iter_tic = tic;
    
    % ---- E-step ----
    [Psijj_cell, e_stats] = module2_estep(Svv_cell, L, Sigma_source_curr, noise_cov);
    outs.loglik(em_iter) = e_stats.log_likelihood;
    
    % ---- Whitening ----
    [Sjj_tilde, D_cell, ~] = module1_data_whitening(Psijj_cell, 'smoothing_window', 1);
    for f = 1:F
        m5_in_base.whitened_covariances{f} = gather(Sjj_tilde{f});
    end
    
    % ---- Recompute scales ----
    [sc_it, ~] = compute_scales_thresholds_(m5_in_base.whitened_covariances, m5_in_base.smoothing_kernel, m5_in_base.weight_matrix, Nr, cfg);
    lam1 = eta_best(1) * sc_it.lambda1_scale;
    lam2 = eta_best(2) * sc_it.lambda2_max;
    lam3 = 0;
    if sc_it.lk_norm > eps
        lam3 = eta_best(3) * sc_it.lambda3_scale / (sc_it.lk_norm + eps);
    end
    m5_p_em.lambda1 = lam1;
    m5_p_em.lambda2 = lam2;
    m5_p_em.lambda3 = lam3;
    outs.em_lambdas(em_iter, :) = [lam1, lam2, lam3];
    
    % ---- Build new mask ----
    [mask_new, ~] = module3_active_set(m5_in_base.whitened_covariances, act_params);
    if ~isempty(dwi_connectivity_mask)
        for f = 1:F
            mask_new{f} = mask_new{f} & dwi_connectivity_mask;
            mask_new{f}(1:Nr+1:end) = true;
        end
    end
    if MASK_UNION_MODE
        for f = 1:F
            active_mask_prev{f} = active_mask_prev{f} | mask_new{f};
            active_mask_prev{f}(1:Nr+1:end) = true;
        end
        m5_in_base.active_mask = active_mask_prev;
    else
        m5_in_base.active_mask = mask_new;
        active_mask_prev = mask_new;
    end
    
    % ---- Warm start ----
    local_in = m5_in_base;
    local_in.precision_matrices = cell(F, 1);
    for f = 1:F, local_in.precision_matrices{f} = gather(Gamma_warm_start{f}); end
    
    % ---- M-step ----
    [Gamma_new, ~] = mstep_solver(local_in, m5_p_em);
    
    % [DIAGNOSIS POINT B] Check Gamma Evolution
    dbg_jspace_cells([], Gamma_new, sprintf('After M-step EM iter %d', em_iter));
    
    Gamma_warm_start = Gamma_new;
    
    % ---- Diagnostics ----
    den = density_from_G_pcor_(Gamma_new, PCOR_EPS);
    outs.em_density_trace(em_iter) = den;
    
    % ---- Recolor + inertia ----
    m8_in = struct();
    m8_in.whitened_precision_matrices = Gamma_new;
    m8_in.whitening_matrices          = D_cell;
    m8_in.original_covariances        = Sigma_source_curr;
    recol = module8_recoloring(m8_in, struct('verbose', false));
    Omega_temp = recol.recolored_precision_matrices;
    
    for f = 1:F
        Om = (Omega_temp{f} + Omega_temp{f}')/2;
        [Om_spd, ~] = utils_math.project_spd(Om, 1e-8);
        S_next = inv(Om_spd);
        S_next = (S_next + S_next')/2;
        Sigma_source_curr{f} = (1-UPDATE_RATE)*Sigma_source_curr{f} + UPDATE_RATE*S_next;
    end
    
    if VERBOSE
        fprintf('  > EM %d/%d: LogLik=%.3e | dens=%.2f%% | lambdas=[%.2e %.2e %.2e] | time=%.2fs\n', ...
            em_iter, MAX_EM_ITER, outs.loglik(em_iter), den*100, lam1, lam2, lam3, toc(iter_tic));
    end
end

% ============================================================
% 10) Robust post-processing (Optional)
% ============================================================
if USE_GPU
    target_G = cell(F, 1);
    for f=1:F, target_G{f} = gpuArray(Gamma_new{f}); end
    target_S = Sjj_tilde;
    K_dev = gpuArray(K_freq);
    W_dev = gpuArray(W_gamma);
else
    target_G = Gamma_new;
    target_S = m5_in_base.whitened_covariances;
    K_dev = gather(K_freq);
    W_dev = gather(W_gamma);
end

if DO_POST
    if VERBOSE, fprintf('[J-SPACE-3D] Robust Post-Processing...\n'); end

    % [DIAGNOSIS POINT C1] Before Debias
    dbg_jspace_cells([], target_G, 'Gamma_hat BEFORE debias');

    % 10.1 Debias
    [Gamma_debiased, ~, ~, ~] = module_debias(target_G, target_S, M_SAMPLES);

    % [DIAGNOSIS POINT C2] After Debias
    dbg_jspace_cells([], Gamma_debiased, 'Gamma_debiased AFTER debias');

    % --- Debias sanity guard (Safety Fallback) ---
    bad = false;
    for f = 1:numel(Gamma_debiased)
        Gd = (Gamma_debiased{f}+Gamma_debiased{f}')/2;
        dd = real(diag(Gd));
        if any(~isfinite(dd)) || median(dd) < 1e-6 || median(dd) > 1e6 || rcond(Gd) < 1e-12
            bad = true; break;
        end
        P = -Gd ./ sqrt(max(dd,eps) * max(dd,eps)');
        P(1:size(Gd,1)+1:end) = 0;
        if max(abs(P(:))) > 1.05   % pcor safety check
            bad = true; break;
        end
    end

    if bad
        warning('[JSPACE] Debias produced invalid precision/pcor. Fallback to Gamma_hat.');
        Gamma_debiased = target_G;
    end

    % 10.2 Rayleigh search setup
    ray_params = struct();
    ray_params.lambda1 = outs.em_lambdas(end,1);
    ray_params.lambda3 = outs.em_lambdas(end,3);
    ray_params.weight_mode = 'hadamard';
    ray_params.variance_source = 'hat';
    ray_params.Gamma_hat = target_G;
    ray_params.threshold_domain = 'pcor'; 
    ray_params.r_range = RAY_RANGE;       
    ray_params.density_min = min_density_limit;
    ray_params.density_max = max_density_limit;
    ray_params.density_penalty_weight = get_cfg(cfg, 'density_penalty_weight', 1e5);

    % Adaptive r_grid for pcor-domain Rayleigh
    if strcmpi(ray_params.threshold_domain, 'pcor')
        pvals = [];
        maskL = tril(true(Nr), -1);

        for ff = 1:F
            Gd = Gamma_debiased{ff};
            d  = real(diag(Gd));
            d  = max(d, 1e-12);
            denom = sqrt(d * d.');
            P = -Gd ./ denom;
            P(1:Nr+1:end) = 0;

            tmp = abs(P(maskL));
            tmp = tmp(isfinite(tmp));
            pvals = [pvals; tmp]; %#ok<AGROW>
        end

        if isempty(pvals)
            ray_params.r_range = 0.05:0.05:1.0;
        else
            p90 = quantile(pvals, 0.90);
            p99 = quantile(pvals, 0.99);
            T_samp = M_SAMPLES;
            r_lo = max(0.01, 0.30 * p90 * sqrt(T_samp));
            r_hi = max(r_lo*1.5, 1.50 * p99 * sqrt(T_samp));
            ray_params.r_range = linspace(r_lo, r_hi, 40);
        end

        if VERBOSE
            thr_min = ray_params.r_range(1)/sqrt(M_SAMPLES);
            thr_max = ray_params.r_range(end)/sqrt(M_SAMPLES);
            fprintf('[Dbg-Ray] Auto r_grid(pcor): r=[%.3f, %.3f], thr=[%.3e, %.3e]\n', ...
                ray_params.r_range(1), ray_params.r_range(end), thr_min, thr_max);
        end
    end

    % Call local module_rayleigh_search
    [best_r, mask_cell, ~, ~] = module_rayleigh_search(Gamma_debiased, target_S, M_SAMPLES, K_dev, W_dev, ray_params);

    [ray_den_vec, ray_den_med, ray_den_min, ray_den_max] = mask_density_stats_(mask_cell);
    outs.post.best_r      = best_r;
    outs.post.ray_den     = ray_den_med;
    outs.post.ray_den_vec = ray_den_vec;

    if VERBOSE
        fprintf('  [Rayleigh] best_r=%.2f | ray_den(med)=%.2f%% [min=%.2f%%, max=%.2f%%] | target=[%.2f%%, %.2f%%]\n', ...
            best_r, ray_den_med*100, ray_den_min*100, ray_den_max*100, min_density_limit*100, max_density_limit*100);
    end

    % 10.3 Rescue decision
    need_rescue = (ray_den_med < min_density_limit) || (ray_den_med > max_density_limit);
    refit_mask = cell(F, 1);

    if need_rescue
        if VERBOSE
            fprintf('  [Rescue] Outlier density (%.2f%%). Using correlation rescue at t_rescue=%.4f\n', ...
                ray_den_med*100, thr0.t_rescue);
        end
        thr_rescue = thr0.t_rescue;
        for f = 1:F
            refit_mask{f} = abs(target_S{f}) > thr_rescue;
            refit_mask{f}(1:Nr+1:end) = true;
        end
    else
        for f = 1:F
            refit_mask{f} = mask_cell{f};
            refit_mask{f}(1:Nr+1:end) = true;
        end
    end

    % 10.4 Refit
    refit_in = struct();
    refit_in.whitened_covariances = target_S;
    refit_in.smoothing_kernel     = K_dev;
    refit_in.weight_matrix        = W_dev;
    refit_in.precision_matrices   = target_G;
    refit_in.active_mask          = refit_mask;
    refit_params = m5_p_em;
    refit_params.lambda2  = 0;
    refit_params.lambda3  = max(refit_params.lambda3, 1e-2);
    refit_params.max_iter = get_cfg(cfg, 'post_refit_max_iter', 100);
    refit_params.tol      = get_cfg(cfg, 'post_refit_tol', 1e-4);
    refit_params.verbose  = false;
    try
        [Gamma_refit, ~] = mstep_solver(refit_in, refit_params);
    catch
        if VERBOSE, fprintf('  [Refit] Failed. Fallback to target_G.\n'); end
        Gamma_refit = target_G;
    end

else
    % =========================
    % Post-processing DISABLED
    % =========================
    if VERBOSE, fprintf('[J-SPACE-3D] Post-processing DISABLED. Returning EM Gamma_hat directly.\n'); end
    
    % Pass EM result as the "Refitted" result (for recoloring)
    Gamma_refit = target_G;
    
    % Mask is effectively full (or undefined), needed for recoloring if active_set_masks used
    refit_mask = cell(F,1);
    for f=1:F, refit_mask{f} = true(Nr); end
    
    % Fill outs with NaNs
    outs.post.best_r = NaN;
    outs.post.ray_den = NaN;
    outs.post.ray_den_vec = nan(F,1);
end

% 10.5 Final recolor (Must happen regardless, to restore original scale)
m8_in = struct();
m8_in.whitened_precision_matrices = Gamma_refit;
m8_in.whitening_matrices          = D_cell;
m8_in.active_set_masks            = refit_mask;
m8_in.original_covariances        = Sigma_source_curr;
recol = module8_recoloring(m8_in, struct('verbose', false));
Omega_final = recol.recolored_precision_matrices;

% ============================================================
% 11) Return
% ============================================================
if USE_GPU
    Omega_est = cell(F,1);
    Sigma_src_est = cell(F,1);
    for f = 1:F
        Omega_est{f} = gather(Omega_final{f});
        Sigma_src_est{f} = gather(Sigma_source_curr{f});
    end
else
    Omega_est = Omega_final;
    Sigma_src_est = Sigma_source_curr;
end

outs.global_hyperparams = struct();
outs.global_hyperparams.eta1 = eta_best(1);
outs.global_hyperparams.eta2 = eta_best(2);
outs.global_hyperparams.eta3 = eta_best(3);
outs.global_hyperparams.lambda1_final = outs.em_lambdas(end,1);
outs.global_hyperparams.lambda2_final = outs.em_lambdas(end,2);
outs.global_hyperparams.lambda3_final = outs.em_lambdas(end,3);

if get_cfg(cfg, 'plot', false)
    plot_eta_trials_(outs.opt_results);
end
if DEBUG_PRINT
    fprintf('[Dbg] EM density: start=%.2f%% end=%.2f%%\n', outs.em_density_trace(1)*100, outs.em_density_trace(end)*100);
end
end

% ============================================================
% Helper Functions
% ============================================================
function score = wrapper_objective_eta_aic_density_( ...
    x_log, m5_in, m5_p, sc0, ...
    M, min_den, max_den, ...
    pcor_eps, penW, ...
    Nr, F, solver, ...
    refit_max_iter, refit_tol)

eta = 10.^x_log(:)';
lambda1 = eta(1) * sc0.lambda1_scale;
lambda2 = eta(2) * sc0.lambda2_max;
lambda3 = 0;
if sc0.lk_norm > eps
    lambda3 = eta(3) * sc0.lambda3_scale / (sc0.lk_norm + eps);
end

local_p = m5_p;
local_p.lambda1 = lambda1;
local_p.lambda2 = lambda2;
local_p.lambda3 = lambda3;

try
    [G_biased, ~] = solver(m5_in, local_p);
    [den_med, total_edges] = density_from_G_pcor_(G_biased, pcor_eps);
    
    penalty = 0;
    if den_med < min_den
        penalty = penW * ((min_den - den_med) / max(min_den, eps))^2;
    elseif den_med > max_den
        penalty = penW * ((den_med - max_den) / max(max_den, eps))^2;
    end
    
    do_refit = (penalty < 1e5) && (den_med > 0);
    
    if do_refit
        refit_mask = cell(F,1);
        for f = 1:F
            Gf = G_biased{f};
            Mf = mask_from_G_pcor_(Gf, pcor_eps);
            Mf(1:Nr+1:end) = true;
            refit_mask{f} = Mf;
        end
        refit_in = m5_in;
        refit_in.active_mask = refit_mask;
        refit_in.precision_matrices = G_biased;
        refit_p = local_p;
        refit_p.lambda2 = 0;
        refit_p.lambda3 = 1e-4;
        refit_p.max_iter = refit_max_iter;
        refit_p.tol      = refit_tol;
        refit_p.verbose  = false;
        [G_refit, ~] = solver(refit_in, refit_p);
    else
        G_refit = G_biased;
    end
    
    ld_sum = 0;
    tr_sum = 0;
    for ff = 1:F
        Gf = G_refit{ff};
        [ld_f, valid] = utils_math.safe_log_det(Gf);
        if ~valid
            score = 1e15;
            return;
        end
        ld_sum = ld_sum + ld_f;
        tr_sum = tr_sum + real(trace(m5_in.whitened_covariances{ff} * Gf));
    end
    minus_2_ll = M * (tr_sum - ld_sum);
    aic_score  = minus_2_ll + 2 * total_edges;
    score = aic_score + penalty;
    if isnan(score) || isinf(score), score = 1e15; end
catch
    score = 1e15;
end
end

function [den_med, total_edges] = density_from_G_pcor_(G_cell, pcor_eps)
F = numel(G_cell);
p = size(G_cell{1}, 1);
dens = zeros(F,1);
total_edges = 0;
for f = 1:F
    G = G_cell{f};
    M = mask_from_G_pcor_(G, pcor_eps);
    M(1:p+1:end) = false;
    edges = nnz(M) / 2;
    total_edges = total_edges + edges;
    dens(f) = edges / (p*(p-1)/2);
end
den_med = median(dens);
end

function M = mask_from_G_pcor_(G, pcor_eps)
p = size(G,1);
d = real(diag(G));
d = max(d, 1e-12);
denom = sqrt(d * d.');
P = -G ./ denom;
P(1:p+1:end) = 0;
M = abs(P) > pcor_eps;
end

function [active_mask, stats] = threshold_active_mask_with_floor_(Sjj_cell, thresh, dens_floor)
if ~iscell(Sjj_cell), Sjj_cell = {Sjj_cell}; end
F = numel(Sjj_cell);
p = size(Sjj_cell{1}, 1);
vals = [];
maskL = tril(true(p), -1);
for f = 1:F
    Sf = Sjj_cell{f};
    vals = [vals; abs(Sf(maskL))]; %#ok<AGROW>
end
vals = vals(isfinite(vals));
if isempty(vals)
    t_use = thresh;
else
    dens_now = mean(vals > thresh);
    if dens_now < dens_floor
        q = max(0.0, min(1.0, 1 - dens_floor));
        t_q = quantile(vals, q);
        t_use = min(thresh, t_q);
    else
        t_use = thresh;
    end
end
active_mask = cell(F,1);
stats.num_active_edges = zeros(F,1);
for f = 1:F
    M = abs(Sjj_cell{f}) > t_use;
    M(1:p+1:end) = true;
    active_mask{f} = M;
    stats.num_active_edges(f) = (nnz(M) - p) / 2;
end
stats.density = mean(stats.num_active_edges) / (p*(p-1)/2);
stats.t_used = t_use;
end

function [active_mask, stats] = threshold_active_mask_(Sjj_cell, thresh)
if ~iscell(Sjj_cell), Sjj_cell = {Sjj_cell}; end
F = numel(Sjj_cell);
p = size(Sjj_cell{1}, 1);
active_mask = cell(F, 1);
stats.num_active_edges = zeros(F, 1);
for f = 1:F
    M = abs(Sjj_cell{f}) > thresh;
    M(1:p+1:end) = true;
    active_mask{f} = M;
    stats.num_active_edges(f) = (nnz(M) - p) / 2;
end
stats.density = mean(stats.num_active_edges) / (p*(p-1)/2);
end

function [sc, thr] = compute_scales_thresholds_(Sjj_cpu_cell, K_freq_cpu, W_gamma_cpu, Nr, cfg)
if ~iscell(Sjj_cpu_cell), Sjj_cpu_cell = {Sjj_cpu_cell}; end
F = numel(Sjj_cpu_cell);
p = size(Sjj_cpu_cell{1}, 1);
Sbar = zeros(p, p);
for f = 1:F, Sbar = Sbar + Sjj_cpu_cell{f}; end
Sbar = utils_math.make_hermitian(Sbar / F);
I = eye(p);
if p < 2
    lambda2_max = 0;
else
    maskL = tril(true(p), -1);
    S_off = abs(Sbar(maskL));
    W_off = abs(W_gamma_cpu(maskL));
    ratios = S_off ./ max(W_off, eps);
    ratios = ratios(isfinite(ratios));
    if isempty(ratios), lambda2_max = 0; else, lambda2_max = max(ratios); end
end
lambda1_scale = norm(Sbar - I, 'fro') * F;
if ~isfinite(lambda1_scale), lambda1_scale = 0; end
if F <= 1
    lambda3_scale = 0;
else
    Gbar = zeros(p, p);
    for f = 1:F, Gbar = Gbar + (Sjj_cpu_cell{f} - I); end
    Gbar = Gbar / F;
    gdiff = zeros(F, 1);
    for f = 1:F
        gdiff(f) = norm((Sjj_cpu_cell{f} - I) - Gbar, 'fro');
    end
    lambda3_scale = median(gdiff);
end
if ~isfinite(lambda3_scale), lambda3_scale = 0; end
lk_norm = 0;
if F > 1
    Ksym = (K_freq_cpu + K_freq_cpu') / 2;
    Lk = diag(sum(Ksym, 2)) - Ksym;
    lk_norm = norm(Lk, 2);
end
vals = [];
maskL = tril(true(p), -1);
for f = 1:F
    Sf = Sjj_cpu_cell{f};
    vals = [vals; abs(Sf(maskL))]; %#ok<AGROW>
end
[t_active, ~] = utils_stats_gmm_threshold_1d(vals);
if isempty(vals)
    q95 = 0;
    data_density = 0.05;
else
    q95 = quantile(vals, 0.95);
    data_density = mean(vals > t_active);
end
sc = struct();
sc.lambda1_scale = lambda1_scale;
sc.lambda2_max   = lambda2_max;
sc.lambda3_scale = lambda3_scale;
sc.lk_norm       = lk_norm;
thr = struct();
thr.t_active     = t_active;
thr.t_rescue     = max(t_active, q95);
thr.data_density = data_density;
if isfield(cfg,'t_active_override') && ~isempty(cfg.t_active_override)
    thr.t_active = cfg.t_active_override;
end
end

% DWI Prior Builder with Debug Print
function [W_gamma, dwi_mask] = build_dwi_soft_prior_(C_in, Nr, cfg)
C = gather(C_in);
C = real(C);
C = max(C, 0);
C(1:Nr+1:end) = 0;
mx = max(C(:));
if mx <= 0
    W_gamma = ones(Nr);
    dwi_mask = [];
    return;
end
C_norm = C / mx;
if isfield(cfg,'dwi_weight_mode') && ~isempty(cfg.dwi_weight_mode)
    mode = lower(cfg.dwi_weight_mode);
else
    mode = 'power';
end
if isfield(cfg,'dwi_weight_alpha') && ~isempty(cfg.dwi_weight_alpha)
    alpha = cfg.dwi_weight_alpha;
else
    alpha = get_cfg(cfg, 'dwi_alpha', 2.0);
end
w_floor = get_cfg(cfg, 'dwi_w_floor', 0.2);
if isfield(cfg,'dwi_weight_clip') && ~isempty(cfg.dwi_weight_clip)
    w_clip = cfg.dwi_weight_clip;
else
    w_clip = [];
end
if isfield(cfg,'dwi_weight_normalize_median') && ~isempty(cfg.dwi_weight_normalize_median)
    do_med_norm = logical(cfg.dwi_weight_normalize_median);
else
    do_med_norm = false;
end
q = get_cfg(cfg, 'dwi_mask_quantile', 0.90);
% Debug Print for DWI
if get_cfg(cfg,'verbose',true)
    fprintf('[DWI] alpha=%.3f | w_floor=%.3f | mask_quantile=%.2f\n', alpha, w_floor, q);
end
switch mode
    case 'power'
        base = (1 - C_norm).^alpha;
    case 'exp'
        base = exp(-alpha * C_norm);
    otherwise
        base = (1 - C_norm).^alpha;
end
W_gamma = w_floor + (1 - w_floor) * base;
W_gamma(1:Nr+1:end) = 1;
if do_med_norm
    off = W_gamma(tril(true(Nr), -1));
    off = off(isfinite(off) & off>0);
    if ~isempty(off)
        med = median(off);
        if med > 0
            W_gamma = W_gamma / med;
            W_gamma(1:Nr+1:end) = 1;
        end
    end
end
if ~isempty(w_clip) && numel(w_clip)==2
    W_gamma = min(max(W_gamma, w_clip(1)), w_clip(2));
    W_gamma(1:Nr+1:end) = 1;
end
vals = C_norm(tril(true(Nr), -1));
vals = vals(vals > 0 & isfinite(vals));
if isempty(vals)
    dwi_mask = [];
else
    thr = quantile(vals, q);
    dwi_mask = (C_norm >= thr);
    dwi_mask = dwi_mask | dwi_mask.';
    dwi_mask(1:Nr+1:end) = true;
end
end

function K = build_default_freq_kernel_(F)
if F <= 1
    K = 1;
    return;
end
K = eye(F);
for i = 1:(F-1)
    K(i, i+1) = 1;
    K(i+1, i) = 1;
end
end

function [T, W] = run_eloreta_core(L, Svv, regu)
[nchan, ndum] = size(L);
if nargin < 3, regu = 0.05; end
W = eye(ndum, 'like', L);
for k = 1:15
    K = (L * W) * L.';
    I_n = eye(nchan, 'like', L);
    alpha = regu * trace(K) / nchan;
    M = (K + alpha * I_n) \ I_n;
    W_old = W;
    for i = 1:ndum
        li = L(:, i);
        val = real(li' * M * li);
        W(i, i) = sqrt(complex(max(val, 1e-12)));
    end
    if norm(diag(W)-diag(W_old))/(norm(diag(W_old))+1e-12) < 1e-3
        break;
    end
end
I_n = eye(nchan, 'like', L);
K_final = (L * W) * L.';
alpha = regu * trace(K_final) / nchan;
T = W * (L.' * ((K_final + alpha * I_n) \ I_n));
end

function plot_eta_trials_(opt_results)
if ~isfield(opt_results,'trials') || ~isfield(opt_results.trials,'X'), return; end
X = opt_results.trials.X;
Fv = opt_results.trials.Fval;
if isempty(X) || isempty(Fv), return; end
figure('Name','Eta Trials','Color','w');
scatter3(log10(X(:,1)), log10(X(:,2)), log10(X(:,3)), 40, Fv, 'filled');
grid on; colorbar;
xlabel('log10(eta1)'); ylabel('log10(eta2)'); zlabel('log10(eta3)');
title('surrogateopt trials (eta space)');
end

function val = get_cfg(s, f, d)
if isfield(s, f), val = s.(f); else, val = d; end
end

function [dens_vec, dens_med, dens_min, dens_max] = mask_density_stats_(mask_cell)
F = numel(mask_cell);
p = size(mask_cell{1},1);
dens_vec = zeros(F,1);
for f = 1:F
    M = mask_cell{f};
    M(1:p+1:end) = false; % remove diag
    edges = nnz(M)/2;
    dens_vec(f) = edges / (p*(p-1)/2);
end
dens_med = median(dens_vec);
dens_min = min(dens_vec);
dens_max = max(dens_vec);
end

function dbg_print_G_stats_(tag, G_cell)
F = numel(G_cell);
p = size(G_cell{1},1);
diag_min = zeros(F,1); diag_med = zeros(F,1); diag_max = zeros(F,1);
pcor_q = zeros(F,3); 
off_q  = zeros(F,3);
for f=1:F
    G = G_cell{f};
    d = real(diag(G)); d = max(d,1e-12);
    diag_min(f)=min(d); diag_med(f)=median(d); diag_max(f)=max(d);
    maskL = tril(true(p),-1);
    off = abs(G(maskL)); off = off(isfinite(off));
    if isempty(off), off_q(f,:) = [NaN NaN NaN];
    else, off_q(f,:) = quantile(off,[0.5 0.9 0.99]);
    end
    denom = sqrt(d*d');
    P = -G ./ denom; P(1:p+1:end)=0;
    poff = abs(P(maskL)); poff = poff(isfinite(poff));
    if isempty(poff), pcor_q(f,:) = [NaN NaN NaN];
    else, pcor_q(f,:) = quantile(poff,[0.5 0.9 0.99]);
    end
end
fprintf('  [%s] diag(min/med/max): med=%.2e | [min=%.2e, max=%.2e]\n', ...
    tag, median(diag_med), min(diag_min), max(diag_max));
fprintf('  [%s] |G_off| q50/q90/q99 (median over f): %.2e / %.2e / %.2e\n', ...
    tag, median(off_q(:,1),'omitnan'), median(off_q(:,2),'omitnan'), median(off_q(:,3),'omitnan'));
fprintf('  [%s] |P_off| q50/q90/q99 (median over f): %.2e / %.2e / %.2e\n', ...
    tag, median(pcor_q(:,1),'omitnan'), median(pcor_q(:,2),'omitnan'), median(pcor_q(:,3),'omitnan'));
end

function dbg_jspace_cells(Sigma_cells, Gamma_cells, tag)
% dbg_jspace_cells(Sigma_cells, Gamma_cells, tag)
% Sigma_cells / Gamma_cells 可以传 []，函数会自动跳过

    fprintf('\n[DBG-JSPACE] %s\n', tag);

    if ~isempty(Sigma_cells)
        F = numel(Sigma_cells);
        dmed = zeros(F,1);
        dmin = zeros(F,1);
        dmax = zeros(F,1);
        rc   = zeros(F,1);
        trv  = zeros(F,1);

        for f = 1:F
            S = (Sigma_cells{f}+Sigma_cells{f}')/2;
            dd = real(diag(S));
            dmed(f) = median(dd);
            dmin(f) = min(dd);
            dmax(f) = max(dd);
            trv(f)  = real(trace(S));
            rc(f)   = rcond(S);
        end

        fprintf('  [Sigma] diag med(min/med/max) over f: min=%.2e | med=%.2e | max=%.2e\n', ...
            min(dmed), median(dmed), max(dmed));
        fprintf('  [Sigma] diag min(min/med/max) over f: min=%.2e | med=%.2e | max=%.2e\n', ...
            min(dmin), median(dmin), max(dmin));
        fprintf('  [Sigma] rcond(min/med/max) over f   : min=%.1e | med=%.1e | max=%.1e\n', ...
            min(rc), median(rc), max(rc));
        fprintf('  [Sigma] trace(min/med/max) over f   : min=%.2e | med=%.2e | max=%.2e\n', ...
            min(trv), median(trv), max(trv));
    end

    if ~isempty(Gamma_cells)
        F = numel(Gamma_cells);
        dmed = zeros(F,1);
        dmin = zeros(F,1);
        dmax = zeros(F,1);
        rc   = zeros(F,1);
        pmax = zeros(F,1);

        for f = 1:F
            G = (Gamma_cells{f}+Gamma_cells{f}')/2;
            dd = real(diag(G));
            dmed(f) = median(dd);
            dmin(f) = min(dd);
            dmax(f) = max(dd);
            rc(f)   = rcond(G);

            % partial correlation max(|P_ij|) offdiag
            P = -G ./ sqrt(max(dd,eps) * max(dd,eps)');
            P(1:size(G,1)+1:end) = 0;
            pmax(f) = max(abs(P(:)));
        end

        fprintf('  [Gamma] diag med(min/med/max) over f: min=%.2e | med=%.2e | max=%.2e\n', ...
            min(dmed), median(dmed), max(dmed));
        fprintf('  [Gamma] rcond(min/med/max) over f   : min=%.1e | med=%.1e | max=%.1e\n', ...
            min(rc), median(rc), max(rc));
        fprintf('  [Gamma] max|P_off|(min/med/max) over f: min=%.2e | med=%.2e | max=%.2e\n', ...
            min(pmax), median(pmax), max(pmax));
    end
end
