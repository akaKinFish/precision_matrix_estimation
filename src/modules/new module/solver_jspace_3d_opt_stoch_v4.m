function [Omega_est, Sigma_src_est, outs] = solver_jspace_3d_opt_stoch_v4(Svv_cell, L, GraphLaplacian, cfg)
% SOLVER_JSPACE_3D_OPT_STOCH
%   J-SPACE Scheme-1: eta-search (surrogateopt) + stable EM + Stochastic M-step.
%
%   Key Features:
%   - Uses "Common Random Numbers" (CRN) for surrogate objective stability.
%   - M-step delegated to STOCH-FISTA (module5_stoch_fista_main) on mini-batches.
%   - Configurable Mode A (independent freq) / Mode B (post-hoc smoothing).
%
% Dependencies (External Files):
%   - module2_estep
%   - module1_data_whitening
%   - module5_stoch_fista_main
%   - module_debias
%   - module_rayleigh_search
%   - module8_recoloring
%   - utils_math.* (project_spd, make_hermitian, safe_log_det)
%   - utils_stats_gmm_threshold_1d
%
if nargin < 4, cfg = struct(); end

% ============================================================
% 1) Input Normalization & Config
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

% Global Config
MAX_EM_ITER        = get_cfg(cfg, 'max_em_iter', 10);
VERBOSE            = get_cfg(cfg, 'verbose', true);
USE_GPU            = get_cfg(cfg, 'use_gpu', false);
DO_POST            = get_cfg(cfg, 'postprocess_enable', true);

% Surrogate Config
MAX_OPT_EVALS      = get_cfg(cfg, 'opt_max_evals', 50);
MIN_OPT_POINTS     = get_cfg(cfg, 'opt_min_points', 10);
USE_PARALLEL       = get_cfg(cfg, 'opt_use_parallel', true);

% Solver & Penalty Config
OBJ_STOCH_MAX_ITER = get_cfg(cfg, 'obj_stoch_max_iter', 10);
EM_STOCH_MAX_ITER  = get_cfg(cfg, 'em_stoch_max_iter', 30);
PCOR_EPS           = get_cfg(cfg, 'dens_pcor_eps', 3e-3);
DENS_PEN_W         = get_cfg(cfg, 'density_penalty_weight', 1e7);
MASK_DENS_FLOOR    = get_cfg(cfg, 'mask_density_floor', max(0.03, log(Nr)/Nr));
MASK_UNION_MODE = get_cfg(cfg, 'mask_union', true);

RAY_RANGE          = get_cfg(cfg, 'rayleigh_range', 0.7:0.1:3.16);
M_SAMPLES          = get_cfg(cfg, 'm_samples', 100 * Nr);
UPDATE_RATE        = get_cfg(cfg, 'update_rate', 0.35);
UPDATE_RATE = 1;

% Post-EM post-processing config
POST_RUN_SEARCH    = get_cfg(cfg, 'post_run_search', true);
POST_RESCUE_MODE   = lower(get_cfg(cfg, 'post_rescue_mode', 'auto')); % off | auto | always
POST_RESCUE_UNION  = get_cfg(cfg, 'post_rescue_union', false);
POST_REFIT_MODE    = lower(get_cfg(cfg, 'post_refit_mode', 'fused_ridge_surrogate')); % debias_only | single_ridge | fused_ridge_surrogate | stoch_mstep
POST_SEARCH_SCORE  = lower(get_cfg(cfg, 'post_search_score', 'hybrid')); % higgs | jspace | hybrid
POST_RIDGE_PENALTY = get_cfg(cfg, 'post_ridge_penalty', []);
POST_MIN_EIG       = get_cfg(cfg, 'post_min_eig', get_cfg(cfg, 'min_eig', 1e-8));
POST_STORE_AUX     = get_cfg(cfg, 'post_store_aux', false);

% Stoch-Specific Config (helper at bottom)
stoch_cfg = get_stoch_cfg_(cfg);

% ============================================================
% 2) Pre-computation (Laplacian, Kernel, Weights)
% ============================================================
if isempty(GraphLaplacian)
    GraphLaplacian = zeros(Nr, Nr, 'like', L);
end

if isfield(cfg, 'freq_kernel') && ~isempty(cfg.freq_kernel)
    K_freq = cfg.freq_kernel;
else
    K_freq = build_default_freq_kernel_(F);
end

W_gamma = get_cfg(cfg, 'weight_matrix', ones(Nr));
W_gamma(1:Nr+1:end) = 0;
dwi_aux = struct();
dwi_mask = [];
if isfield(cfg, 'dwi_C') && ~isempty(cfg.dwi_C)
    if VERBOSE, fprintf('[J-SPACE-3D-STOCH] Using DWI soft prior.\n'); end
    [W_gamma, dwi_mask, dwi_aux] = build_dwi_soft_prior_(cfg.dwi_C, Nr, cfg);
end

% Noise Covariance
tr_S = 0; for f = 1:F, tr_S = tr_S + trace(Svv_cell{f}); end
noise_cov = (tr_S / (F * Ns)) * 0.05 * eye(Ns);

% GPU Setup
if USE_GPU
    try
        L = gpuArray(L);
        noise_cov = gpuArray(noise_cov);
        for f = 1:F, Svv_cell{f} = gpuArray(Svv_cell{f}); end
        GraphLaplacian = gpuArray(GraphLaplacian);
        if VERBOSE, fprintf('[J-SPACE-3D-STOCH] GPU Acceleration ENABLED.\n'); end
    catch
        USE_GPU = false;
        if VERBOSE, fprintf('[J-SPACE-3D-STOCH] GPU init failed. Falling back to CPU.\n'); end
    end
end

% ============================================================
% 3) Initialization (eLORETA)
% ============================================================
Sigma_source_curr = [];
if isfield(cfg, 'init_sigma_source') && ~isempty(cfg.init_sigma_source)
    Sigma_source_curr = cfg.init_sigma_source;
    if ~iscell(Sigma_source_curr), Sigma_source_curr = {Sigma_source_curr}; end
else
    if VERBOSE, fprintf('[Init] Running eLORETA...\n'); end
    % [EXTERNAL DEPENDENCY] We use a simple local version or external 'run_eloreta_core' if available
    S_avg = zeros(Ns, Ns, 'like', Svv_cell{1});
    for f = 1:F, S_avg = S_avg + Svv_cell{f}; end
    S_avg = S_avg / F;
    [T_eloreta, ~] = run_eloreta_core_(L, S_avg, 0.05); % Local helper below

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
% 4) One E-step + Whitening (Base Data for Opt)
% ============================================================
if VERBOSE, fprintf('[J-SPACE-3D-STOCH] Preparing data for optimization...\n'); end

% [EXTERNAL DEPENDENCY] module2_estep
[Psijj_cell, ~] = module2_estep(Svv_cell, L, Sigma_source_curr, noise_cov);

% [EXTERNAL DEPENDENCY] module1_data_whitening
[Sjj_tilde0, D_cell, ~] = module1_data_whitening(Psijj_cell, 'smoothing_window', 1);

% Prepare M-step inputs
m5_in_base = struct();
m5_in_base.whitened_covariances = cell(F, 1);
for f = 1:F, m5_in_base.whitened_covariances{f} = gather(Sjj_tilde0{f}); end
m5_in_base.smoothing_kernel = gather(K_freq);
m5_in_base.weight_matrix    = gather(W_gamma);

m5_p_base = struct();
m5_p_base.alpha0 = 0.1;
m5_p_base.min_eig = get_cfg(cfg, 'min_eig', 1e-6);
m5_p_base.spatial_graph_matrix = gather(GraphLaplacian);
m5_p_base.spatial_graph_is_laplacian = true;
m5_p_base.verbose = false;
m5_p_base.weight_mode = 'hadamard';
m5_p_base.auto_tune = false;
m5_p_base.stoch = stoch_cfg; % Pass stoch config to module5
m5_p_base.stoch.max_iter = OBJ_STOCH_MAX_ITER;

% [EXTERNAL DEPENDENCY] Module 5 handle
% Ensure 'module5_stoch_fista_main.m' is on your path
mstep_solver = @module5_stoch_fista_main;

% ============================================================
% 5) Compute Thresholds & Eta Bounds
% ============================================================
if VERBOSE, fprintf('[J-SPACE-3D-STOCH] Computing scales/thresholds...\n'); end

[sc0, thr0] = compute_scales_thresholds_(m5_in_base.whitened_covariances, m5_in_base.smoothing_kernel, m5_in_base.weight_matrix, Nr, cfg);

min_den_limit = log(Nr) / Nr;
max_den_limit = min(thr0.data_density, get_cfg(cfg, 'max_density_cap', 0.40));
max_den_limit = max(max_den_limit, min_den_limit * 1.5);

% Initial Mask
[mask0, ~] = threshold_active_mask_(m5_in_base.whitened_covariances, thr0.t_active);
active_initial_packed = [];
if get_cfg(cfg, 'store_active_masks', false)
    active_initial_packed = pack_mask_cell_lower_(mask0);
end
if ~isempty(dwi_mask)
    for f = 1:F, mask0{f} = mask0{f} | dwi_mask; mask0{f}(1:Nr+1:end) = true; end
end
m5_in_base.active_mask = mask0;

% Eta Bounds Setup
eta_bounds = prepare_eta_bounds_(cfg);
initPts = eta_bounds.x0_log;
if isfield(cfg,'opt_initial_points_log') && ~isempty(cfg.opt_initial_points_log)
    initPts = cfg.opt_initial_points_log;
end

%%% [CHECK-1] print REAL bounds used by surrogateopt (after shrink)
if VERBOSE
    fprintf('[eta bounds] lb=[%.3e %.3e %.3e] ub=[%.3e %.3e %.3e]\n', ...
        eta_bounds.lb(1), eta_bounds.lb(2), eta_bounds.lb(3), ...
        eta_bounds.ub(1), eta_bounds.ub(2), eta_bounds.ub(3));
    fprintf('[eta bounds log10] lb=[%.3f %.3f %.3f] ub=[%.3f %.3f %.3f]\n', ...
        eta_bounds.lb_log(1), eta_bounds.lb_log(2), eta_bounds.lb_log(3), ...
        eta_bounds.ub_log(1), eta_bounds.ub_log(2), eta_bounds.ub_log(3));
    fprintf('[eta init log10] x0=[%.3f %.3f %.3f]  (#initPts=%d)\n', ...
        eta_bounds.x0_log(1), eta_bounds.x0_log(2), eta_bounds.x0_log(3), size(initPts,1));
end

% ============================================================
% 6) Surrogate Optimization (CRN Seeded)
% ============================================================
opt_seed_base = get_cfg(cfg,'opt_rng_seed', 0);

% Objective Wrapper with Stoch Config
obj_fun = @(x_log) wrapper_objective_stoch_( ...
    x_log, m5_in_base, m5_p_base, sc0, ...
    M_SAMPLES, min_den_limit, max_den_limit, ...
    PCOR_EPS, DENS_PEN_W, ...
    Nr, F, mstep_solver, ...
    opt_seed_base);

%%% quick determinism check
xchk = eta_bounds.x0_log(:).';
s1 = obj_fun(xchk);
s2 = obj_fun(xchk);
fprintf('[obj determinism] f(x0) run1=%.6e run2=%.6e | rel_diff=%.3e\n', ...
    s1, s2, abs(s1-s2)/max(1,abs(s1)));

opts = optimoptions('surrogateopt', ...
    'MaxFunctionEvaluations', MAX_OPT_EVALS, ...
    'MinSurrogatePoints', MIN_OPT_POINTS, ...
    'UseParallel', USE_PARALLEL, ...
    'PlotFcn', [], ...
    'InitialPoints', initPts);

if VERBOSE, fprintf('[J-SPACE-3D-STOCH] Running surrogateopt on eta...\n'); end

tic_opt = tic;
[x_best_log, fval, exitflag, output_opt, trials] = surrogateopt(obj_fun, eta_bounds.lb_log, eta_bounds.ub_log, opts);
time_opt = toc(tic_opt);

%%% [CHECK-2] analyze trials: is it still improving near the end?
try
    fv = trials.Fval(:);
    nfv = numel(fv);
    k0 = max(1, nfv-9);
    best_all = min(fv);
    best_last10 = min(fv(k0:end));
    fprintf('[surrogateopt trials] evals=%d | best_all=%.3e | best_last10=%.3e\n', ...
        nfv, best_all, best_last10);

    % optional quick "still improving?" indicator
    if best_last10 < 0.99 * best_all
        fprintf('[surrogateopt trials] WARNING: last10 best is still improving a lot -> MaxFunctionEvaluations may be too small.\n');
    end
catch ME
    fprintf('[surrogateopt trials] (analysis failed) %s\n', ME.message);
end

eta_best = 10.^x_best_log;

%%% [CHECK-3] boundary hit check (relative tolerance)
tol_rel = 1e-3;  % 0.1% relative tolerance
lb = eta_bounds.lb(:).';
ub = eta_bounds.ub(:).';
eta = eta_best(:).';

hit_lb = abs(eta - lb) ./ max(abs(lb), eps) < tol_rel;
hit_ub = abs(eta - ub) ./ max(abs(ub), eps) < tol_rel;

if VERBOSE
    fprintf('[eta best] eta=[%.3e %.3e %.3e]\n', eta(1), eta(2), eta(3));
    fprintf('[eta hit]  hit_lb=[%d %d %d] hit_ub=[%d %d %d] (tol_rel=%.1e)\n', ...
        hit_lb(1), hit_lb(2), hit_lb(3), hit_ub(1), hit_ub(2), hit_ub(3), tol_rel);
end

outs = struct();
outs.meta.Ns = Ns;
outs.meta.Nr = Nr;
outs.meta.F = F;
outs.meta.T_spec = M_SAMPLES;
outs.meta.m_samples = M_SAMPLES;
outs.meta.freq_kernel = gather(K_freq);

outs.dwi_prior = dwi_aux;
outs.dwi_prior.has_dwi_mask = ~isempty(dwi_mask);
if ~isempty(dwi_mask)
    Mtmp = dwi_mask;
    Mtmp(1:Nr+1:end) = false;
    outs.dwi_prior.mask_density = nnz(tril(Mtmp, -1)) / (Nr*(Nr-1)/2);
else
    outs.dwi_prior.mask_density = NaN;
end
if get_cfg(cfg, 'store_weight_matrix', false)
    outs.dwi_prior.W_gamma = gather(W_gamma);
end
outs.eta_bounds = eta_bounds; % Save bounds info from [CHECK-1]
outs.opt_results.fval = fval;
outs.opt_results.exitflag = exitflag;
outs.opt_results.output = output_opt;
outs.opt_results.trials = trials;
outs.opt_results.trials.X_log = trials.X;     % Save log-space trials from [CHECK-2]
outs.opt_results.trials.X = 10.^trials.X;     % Save linear-space trials from [CHECK-2]
outs.best_eta.eta1 = eta_best(1);
outs.best_eta.eta2 = eta_best(2);
outs.best_eta.eta3 = eta_best(3);
outs.best_eta.hit_lb = hit_lb; % Save hit info from [CHECK-3]
outs.best_eta.hit_ub = hit_ub; % Save hit info from [CHECK-3]
outs.thresholds = thr0;
outs.density_limits.min = min_den_limit;
outs.density_limits.max = max_den_limit;

if VERBOSE
    fprintf('  > Opt finished in %.2fs. Best Score: %.3e\n', time_opt, fval);
    fprintf('  > BEST eta: eta1=%.3e | eta2=%.3e | eta3=%.3e\n', eta_best(1), eta_best(2), eta_best(3));
end

% ============================================================
% 7) Full EM with Dynamic Lambdas (M-step = Stoch-FISTA)
% ============================================================
if VERBOSE, fprintf('\n[J-SPACE-3D-STOCH] Running Full EM...\n'); end

m5_p_em = m5_p_base;
m5_p_em.stoch.max_iter = EM_STOCH_MAX_ITER;

outs.loglik = zeros(MAX_EM_ITER, 1);
outs.em_density_trace = zeros(MAX_EM_ITER, 1);
outs.em_lambdas = zeros(MAX_EM_ITER, 3);
% Pre-allocate stats arrays [CHECK-4]
outs.em_active_mask_density = zeros(MAX_EM_ITER,1);
outs.em_active_t_used = zeros(MAX_EM_ITER,1);

active_mask_prev = m5_in_base.active_mask;

for em_iter = 1:MAX_EM_ITER
    iter_tic = tic;

    % ---- E-step ----
    [Psijj_cell, e_stats] = module2_estep(Svv_cell, L, Sigma_source_curr, noise_cov);
    outs.loglik(em_iter) = e_stats.log_likelihood;

    % ---- Whitening ----
    [Sjj_tilde, D_cell, ~] = module1_data_whitening(Psijj_cell, 'smoothing_window', 1);
    for f = 1:F, m5_in_base.whitened_covariances{f} = gather(Sjj_tilde{f}); end

    % ---- Lambdas ----
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

    % ---- Update Active Mask ----
    %%% [CHECK-4] capture mask stats (density + t_used)
    [mask_new, mask_stats] = threshold_active_mask_with_floor_( ...
        m5_in_base.whitened_covariances, thr0.t_active, MASK_DENS_FLOOR);

    outs.em_active_mask_density(em_iter) = mask_stats.density;
    outs.em_active_t_used(em_iter) = mask_stats.t_used;

    if ~isempty(dwi_mask)
        for f = 1:F, mask_new{f} = mask_new{f} | dwi_mask; mask_new{f}(1:Nr+1:end) = true; end
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

    % ---- M-step (STOCH) ----
    local_in = m5_in_base;
    local_in.precision_matrices = cell(F, 1);
    for f = 1:F, local_in.precision_matrices{f} = gather(Gamma_warm_start{f}); end

    % Set deterministic seed based on EM iter for reproducibility
    m5_p_em.stoch.seed = stoch_cfg.seed + em_iter;
    % m5_p_em.stoch.band_groups = {1:F};   % 只有一个 band，包含全部频点
    m5_p_em.stoch.m_per_band = 3;  
    [Gamma_new, res] = mstep_solver(local_in, m5_p_em);
    %%% batch coverage in this M-step call
    covered = false(1,F);
    for k = 1:numel(res.batch_trace)
        covered(res.batch_trace{k}) = true;
    end
    cov_rate = mean(covered);
    if VERBOSE
        fprintf('    [stoch coverage] covered=%.1f%% of freqs in this M-step call\n', 100*cov_rate);
    end

    Gamma_warm_start = Gamma_new;

    % ---- Diagnostics ----
    den = density_from_G_pcor_(Gamma_new, PCOR_EPS);
    outs.em_density_trace(em_iter) = den;

    %%% pcor magnitude diagnostics
    pcor_max = 0;
    pcor_p99 = 0;
    for ff = 1:F
        G = Gamma_new{ff};
        p = size(G,1);
        d = max(real(diag(G)), 1e-12);
        P = -G ./ sqrt(d*d.');
        P(1:p+1:end) = 0;
        v = abs(P(tril(true(p),-1)));
        v = v(isfinite(v));
        if ~isempty(v)
            pcor_max = max(pcor_max, max(v));
            pcor_p99 = max(pcor_p99, quantile(v, 0.99));
        end
    end
    if VERBOSE
        fprintf('    [pcor stats] max=%.2e | p99=%.2e | eps=%.2e\n', pcor_max, pcor_p99, PCOR_EPS);
    end

    % ---- Recolor + Inertia ----
    % [EXTERNAL DEPENDENCY] module8_recoloring
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
        % Updated print with mask stats [CHECK-4]
        fprintf('  > EM %d/%d: LogLik=%.3e | dens=%.2f%% | maskDen=%.2f%% | t_used=%.2e | lambdas=[%.2e %.2e %.2e] | time=%.2fs\n', ...
            em_iter, MAX_EM_ITER, outs.loglik(em_iter), den*100, ...
            outs.em_active_mask_density(em_iter)*100, outs.em_active_t_used(em_iter), ...
            lam1, lam2, lam3, toc(iter_tic));
    end
end

% ============================================================
% 8) Post-processing (Debias -> Rayleigh Search -> Rescue -> HIGGS-like Refit)
% ============================================================
target_G = Gamma_new;
target_S = m5_in_base.whitened_covariances;
K_dev = gather(K_freq);
W_dev = gather(W_gamma);

if get_cfg(cfg, 'store_active_masks', false)
    outs.active.final_mask = pack_mask_cell_lower_(active_mask_prev);
    outs.active.final_mask_format = 'lower-triangular linear indices per frequency';
    outs.active.final_mask_p = Nr;
    outs.active.final_mask_F = F;
    outs.active.initial_mask = active_initial_packed;

end

if DO_POST
    if VERBOSE
        fprintf('[J-SPACE-3D-STOCH] Post-processing (debias -> search -> rescue -> refit)...\n');
    end

    ridge_penalty_post = POST_RIDGE_PENALTY;
    if isempty(ridge_penalty_post)
        ridge_penalty_post = max(outs.em_lambdas(end,2)^2, 1e-3);
    end

    % ---- Debias only: no fixed Rayleigh threshold inside module_debias ----
    debias_params = struct();
    debias_params.variance_source = 'hat';
    [Gamma_debiased, ~, Var_proxies, ~] = module_debias_v2(target_G, target_S, M_SAMPLES, debias_params);

    % ---- HIGGS-like single-frequency ridge/Frobenius template ----
    [Gamma_base, base_aux] = module_higgs_ridge_template(target_S, ...
        struct('ridge_penalty', ridge_penalty_post, 'min_eig', POST_MIN_EIG));

    % ---- Rayleigh search: masks from debiased stats, scoring on ridge template ----
    if POST_RUN_SEARCH
        ray_params = struct();
        ray_params.lambda1 = outs.em_lambdas(end,1);
        ray_params.lambda2 = outs.em_lambdas(end,2);
        ray_params.lambda3 = outs.em_lambdas(end,3);
        ray_params.weight_mode = 'hadamard';
        ray_params.variance_source = 'hat';
        ray_params.Gamma_hat = target_G;
        ray_params.threshold_domain = 'entry';
        ray_params.template_cell = Gamma_base;
        ray_params.r_range = RAY_RANGE;
        ray_params.score_mode = POST_SEARCH_SCORE;
        ray_params.density_min = min_den_limit;
        ray_params.density_max = max_den_limit;
        ray_params.density_penalty_weight = DENS_PEN_W;
        ray_params.min_eig = POST_MIN_EIG;

        [best_r, mask_cell, Var_proxies, ray_stats] = module_rayleigh_search_v2( ...
            Gamma_debiased, target_S, M_SAMPLES, K_dev, W_dev, ray_params);

        [ray_den_vec, ray_den_med, ray_den_min, ray_den_max] = mask_density_stats_(mask_cell);
    else
        best_r = NaN;
        mask_cell = cell(F,1);
        for f = 1:F
            mask_cell{f} = true(Nr);
        end
        ray_den_vec = NaN(F,1);
        ray_den_med = NaN;
        ray_den_min = NaN;
        ray_den_max = NaN;
        ray_stats = struct();
    end

    % ---- Rescue thresholds MUST be recomputed on the FINAL whitened covariances ----
    [~, thr_end] = compute_scales_thresholds_(target_S, K_dev, W_dev, Nr, cfg);

    bad_ray = POST_RUN_SEARCH && ( ...
        (~isfinite(best_r)) || ...
        (ray_den_med < min_den_limit) || ...
        (ray_den_med > max_den_limit));

    switch POST_RESCUE_MODE
        case 'off'
            rescue_used = false;
        case 'always'
            rescue_used = true;
        otherwise % auto
            rescue_used = bad_ray;
    end

    refit_mask = cell(F,1);
    if rescue_used
        thr_rescue = thr_end.t_rescue;
        for f = 1:F
            refit_mask{f} = abs(target_S{f}) > thr_rescue;
            refit_mask{f}(1:Nr+1:end) = true;

            if POST_RESCUE_UNION && POST_RUN_SEARCH
                refit_mask{f} = refit_mask{f} | mask_cell{f};
                refit_mask{f}(1:Nr+1:end) = true;
            end

            if ~isempty(dwi_mask)
                refit_mask{f} = refit_mask{f} | dwi_mask;
                refit_mask{f}(1:Nr+1:end) = true;
            end
        end
    else
        for f = 1:F
            refit_mask{f} = mask_cell{f};
            refit_mask{f}(1:Nr+1:end) = true;
        end
    end

    % ---- Refit ----
    refit_aux = struct();
    switch POST_REFIT_MODE
        case 'debias_only'
            Gamma_refit = Gamma_debiased;

        case 'single_ridge'
            [Gamma_refit, refit_aux] = module_higgs_refit(target_S, refit_mask, K_dev, W_dev, ...
                struct('mode', 'single_ridge', ...
                       'lambda1', outs.em_lambdas(end,1), ...
                       'lambda3', outs.em_lambdas(end,3), ...
                       'ridge_penalty', ridge_penalty_post, ...
                       'min_eig', POST_MIN_EIG));

        case 'fused_ridge_surrogate'
            [Gamma_refit, refit_aux] = module_higgs_refit(target_S, refit_mask, K_dev, W_dev, ...
                struct('mode', 'fused_ridge_surrogate', ...
                       'lambda1', outs.em_lambdas(end,1), ...
                       'lambda3', outs.em_lambdas(end,3), ...
                       'ridge_penalty', ridge_penalty_post, ...
                       'min_eig', POST_MIN_EIG));

        otherwise % 'stoch_mstep'
            refit_in = struct();
            refit_in.whitened_covariances = target_S;
            refit_in.smoothing_kernel     = K_dev;
            refit_in.weight_matrix        = W_dev;
            refit_in.precision_matrices   = target_G;
            refit_in.active_mask          = refit_mask;

            refit_params = m5_p_em;
            refit_params.lambda2 = 0; % no sparsity penalty after support selection
            refit_params.stoch.max_iter = EM_STOCH_MAX_ITER;

            try
                [Gamma_refit, refit_aux] = mstep_solver(refit_in, refit_params);
            catch
                if VERBOSE
                    fprintf('  [Refit] Stochastic refit failed. Falling back to HIGGS-like single_ridge.\n');
                end
                [Gamma_refit, refit_aux] = module_higgs_refit(target_S, refit_mask, K_dev, W_dev, ...
                    struct('mode', 'single_ridge', ...
                           'lambda1', outs.em_lambdas(end,1), ...
                           'lambda3', outs.em_lambdas(end,3), ...
                           'ridge_penalty', ridge_penalty_post, ...
                           'min_eig', POST_MIN_EIG));
            end
    end

    % ---- Diagnostics ----
    outs.post.best_r            = best_r;
    outs.post.ray_den           = ray_den_med;
    outs.post.ray_den_vec       = ray_den_vec;
    outs.post.ray_den_min       = ray_den_min;
    outs.post.ray_den_max       = ray_den_max;
    outs.post.rescue_used       = rescue_used;
    outs.post.rescue_mode       = POST_RESCUE_MODE;
    outs.post.search_score_mode = POST_SEARCH_SCORE;
    outs.post.refit_mode        = POST_REFIT_MODE;
    outs.post.post_ridge_penalty = ridge_penalty_post;
    outs.post.t_active_final    = thr_end.t_active;
    outs.post.t_rescue_final    = thr_end.t_rescue;
    outs.post.variance_proxy_source = debias_params.variance_source;
    if get_cfg(cfg, 'store_masks', true)
        outs.post.rayleigh_mask = pack_mask_cell_lower_(mask_cell);
        outs.post.final_mask = pack_mask_cell_lower_(refit_mask);
        outs.post.mask_format = 'lower-triangular linear indices per frequency';
        outs.post.mask_p = Nr;
        outs.post.mask_F = F;

        [final_den_vec, final_den_med, final_den_min, final_den_max] = mask_density_stats_(refit_mask);
        outs.post.final_den_vec = final_den_vec;
        outs.post.final_den = final_den_med;
        outs.post.final_den_min = final_den_min;
        outs.post.final_den_max = final_den_max;
    end
    if POST_STORE_AUX
        outs.post.base_aux    = base_aux;
        outs.post.ray_stats   = ray_stats;
        outs.post.refit_aux   = refit_aux;
        outs.post.Var_proxies = Var_proxies;
    end
else
    % Post-process disabled
    Gamma_refit = target_G;
    refit_mask = cell(F,1);
    for f = 1:F
        refit_mask{f} = true(Nr);
    end
    outs.post.best_r             = NaN;
    outs.post.ray_den            = NaN;
    outs.post.ray_den_vec        = NaN(F,1);
    outs.post.ray_den_min        = NaN;
    outs.post.ray_den_max        = NaN;
    outs.post.rescue_used        = false;
    outs.post.rescue_mode        = 'off';
    outs.post.search_score_mode  = POST_SEARCH_SCORE;
    outs.post.refit_mode         = 'none';
    outs.post.post_ridge_penalty = NaN;
    outs.post.t_active_final     = NaN;
    outs.post.t_rescue_final     = NaN;
    if get_cfg(cfg, 'store_masks', true)
        outs.post.rayleigh_mask = pack_mask_cell_lower_(refit_mask);
        outs.post.final_mask = pack_mask_cell_lower_(refit_mask);
        outs.post.mask_format = 'lower-triangular linear indices per frequency';
        outs.post.mask_p = Nr;
        outs.post.mask_F = F;
    end
end

% Final Recolor
m8_in = struct();
m8_in.whitened_precision_matrices = Gamma_refit;
m8_in.whitening_matrices          = D_cell;
m8_in.active_set_masks            = refit_mask;
m8_in.original_covariances        = Sigma_source_curr;
recol = module8_recoloring(m8_in, struct('verbose', false));
Omega_final = recol.recolored_precision_matrices;

% ============================================================
% 9) Return
% ============================================================
% ============================================================
% 9) Return
% ============================================================
% ============================================================
% 9) Return
% ============================================================
Omega_est = Omega_final;

% ---- 终极修复：基于自适应条件数 (Adaptive Condition Number) 的安全求逆 ----
Sigma_src_est = cell(F, 1);
for f = 1:F
    % 1. 严格对称化最终的精度矩阵
    Om = (Omega_final{f} + Omega_final{f}') / 2;
    
    % 2. 转换到对角线为 1 的标准相关空间 (隔离物理尺度的畸变)
    d_Om = abs(diag(Om)); 
    d_Om(d_Om < 1e-15) = 1e-15; 
    D_scale = diag(1 ./ sqrt(d_Om));
    Om_scaled = D_scale * Om * D_scale;
    Om_scaled = (Om_scaled + Om_scaled') / 2;
    
    % ==========================================================
    % 3. 【自适应正则化核心】动态计算最优 epsilon
    % ==========================================================
    % 提取所有特征值 (由于矩阵是 360x360，eig 计算在毫秒级，毫无压力)
    eig_vals = real(eig(Om_scaled));
    max_eig = max(eig_vals);
    min_eig = min(eig_vals);
    
    % 设定目标最大条件数 (1e4 是复杂网络求逆的黄金安全水位)
    target_cond = 1e4; 
    
    % 计算要达到该条件数，矩阵需要的"最小特征值底线"
    safe_min_eig = max_eig / target_cond;
    
    % 如果当前最小特征值已经跌破底线 (甚至为负)，则自适应补足差值
    if min_eig < safe_min_eig
        epsilon = safe_min_eig - min_eig;
    else
        epsilon = 1e-12; % 如果原本就很健康，只加一个机器极小值防浮点误差
    end
    
    % 施加自适应扰动
    Om_safe = Om_scaled + epsilon * eye(Nr, 'like', Om_scaled);
    % ==========================================================
    
    % 4. 此时求逆极度安全，且信号无损
    S_scaled = inv(Om_safe);
    
    % 5. 还原到真实的物理量级
    S_paired = D_scale * S_scaled * D_scale;
    
    % 6. 保证输出协方差矩阵的绝对对称性
    Sigma_src_est{f} = (S_paired + S_paired') / 2;
end

outs.global_hyperparams.eta1 = eta_best(1);
outs.global_hyperparams.eta2 = eta_best(2);
outs.global_hyperparams.eta3 = eta_best(3);
outs.global_hyperparams.lambda1_final = outs.em_lambdas(end,1);
outs.global_hyperparams.lambda2_final = outs.em_lambdas(end,2);
outs.global_hyperparams.lambda3_final = outs.em_lambdas(end,3);

end

% ============================================================
% LOCAL HELPERS
% ============================================================

function score = wrapper_objective_stoch_(x_log, m5_in, m5_p, sc0, M, min_den, max_den, pcor_eps, penW, Nr, F, solver, seed_base)
% Surrogate Objective Function with CRN Seed
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

% TRUE CRN: same random batches for all x during surrogateopt
local_p.stoch.seed = seed_base;

try
    [G_biased, ~] = solver(m5_in, local_p);

    [den_med, total_edges] = density_from_G_pcor_(G_biased, pcor_eps);

    penalty = 0;
    if den_med < min_den
        penalty = penW * ((min_den - den_med) / max(min_den, eps))^2;
    elseif den_med > max_den
        penalty = penW * ((den_med - max_den) / max(max_den, eps))^2;
    end

    ld_sum = 0;
    tr_sum = 0;
    for ff = 1:F
        Gf = G_biased{ff};
        % [EXTERNAL DEPENDENCY] utils_math.safe_log_det
        [ld_f, valid] = utils_math.safe_log_det(Gf);
        if ~valid
            score = 1e15; return;
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

function h = hash_seed_(x_log)
x = x_log(:);
v = sum(abs(x) .* (1:numel(x))') + 13*sum(x.^2);
h = mod(floor(1e6 * abs(v)), 2^31-1);
end

function cfg_stoch = get_stoch_cfg_(cfg)
stoch = struct();
if isfield(cfg,'stoch') && ~isempty(cfg.stoch)
    stoch = cfg.stoch;
end
stoch.mode = get_cfg(stoch,'mode','B');
stoch.seed = get_cfg(stoch,'seed', get_cfg(cfg,'opt_rng_seed',0));
stoch.neighbor_closure = get_cfg(stoch,'neighbor_closure', 1);
stoch.use_backtracking = get_cfg(stoch,'use_backtracking', true);
stoch.backtracking_factor = get_cfg(stoch,'backtracking_factor', 2.0);
stoch.max_backtracking = get_cfg(stoch,'max_backtracking', 25);

% [D1] Explicit Stoch Controls (New Fields)
stoch.use_nesterov = get_cfg(stoch, 'use_nesterov', false);      % Default: Stable
stoch.L_min = get_cfg(stoch, 'L_min', 1e-3);
stoch.L_max = get_cfg(stoch, 'L_max', 1e6);
stoch.L_shrink_on_success = get_cfg(stoch, 'L_shrink_on_success', true);

if isfield(cfg,'freq'), stoch.freq = cfg.freq; end
cfg_stoch = stoch;
end

function [den_med, total_edges] = density_from_G_pcor_(G_cell, pcor_eps)
F = numel(G_cell);
dens = zeros(F,1);
total_edges = 0;
p = size(G_cell{1},1);
for f = 1:F
    G = G_cell{f};
    M = mask_from_G_pcor_(G, pcor_eps);
    M(1:p+1:end) = false;
    edges = nnz(M)/2;
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
vals = vals(isfinite(vals));

% [EXTERNAL DEPENDENCY] utils_stats_gmm_threshold_1d
[t_active, ~] = utils_stats_gmm_threshold_1d(vals);

if isempty(vals)
    q95 = 0; data_density = 0.05;
else
    q95 = quantile(vals, 0.95);
    data_density = mean(vals > t_active);
end

sc.lambda1_scale = lambda1_scale;
sc.lambda2_max   = lambda2_max;
sc.lambda3_scale = lambda3_scale;
sc.lk_norm       = lk_norm;
thr.t_active     = t_active;
thr.t_rescue     = max(t_active, q95);
thr.data_density = data_density;

if isfield(cfg,'t_active_override') && ~isempty(cfg.t_active_override)
    thr.t_active = cfg.t_active_override;
end
end

function [W_gamma, dwi_mask, aux] = build_dwi_soft_prior_(C_in, Nr, cfg)
C = gather(C_in);
C = real(C);
C = max(C, 0);

% Undirected precision graph: symmetrize DWI affinity.
C = 0.5 * (C + C.');
C(1:Nr+1:end) = 0;

maskL = tril(true(Nr), -1);
mx = max(C(maskL));

aux = struct();
aux.mode = get_cfg(cfg, 'dwi_weight_mode', 'power');
aux.alpha = get_cfg(cfg, 'dwi_weight_alpha', 1.0);
aux.w_floor = get_cfg(cfg, 'dwi_w_floor', 0.2);
aux.mask_quantile = get_cfg(cfg, 'dwi_mask_quantile', 0.90);

if mx <= 0
    W_gamma = ones(Nr);
    W_gamma(1:Nr+1:end) = 0;
    dwi_mask = [];
    aux.available = false;
    aux.W_min = NaN;
    aux.W_median = NaN;
    aux.W_max = NaN;
    aux.dwi_mask_density = NaN;
    return;
end

aux.available = true;
C_norm = C / mx;

mode = aux.mode;
alpha = aux.alpha;
w_floor = aux.w_floor;

switch lower(mode)
    case 'power'
        base = (1 - C_norm).^alpha;
    case 'exp'
        base = exp(-alpha * C_norm);
    case 'linear'
        base = max(0, 1 - alpha * C_norm);
    otherwise
        error('Unknown dwi_weight_mode: %s', mode);
end

W_gamma = w_floor + (1 - w_floor) * base;
W_gamma = 0.5 * (W_gamma + W_gamma.');
W_gamma(1:Nr+1:end) = 0;

wvals = W_gamma(maskL);
aux.W_min = min(wvals);
aux.W_median = median(wvals);
aux.W_max = max(wvals);

q = aux.mask_quantile;
cvals = C_norm(maskL);
tq = quantile(cvals, q);
M = (C_norm >= tq);
M = M | M.';
M(1:Nr+1:end) = true;
dwi_mask = logical(M);

Moff = dwi_mask;
Moff(1:Nr+1:end) = false;
aux.dwi_mask_density = nnz(tril(Moff, -1)) / (Nr*(Nr-1)/2);
end

function K = build_default_freq_kernel_(F)
if F <= 1, K = 1; return; end
K = zeros(F);
for i = 1:(F-1), K(i, i+1) = 1; K(i+1, i) = 1; end
end

function [dens_vec, dens_med, dens_min, dens_max] = mask_density_stats_(mask_cell)
F = numel(mask_cell);
p = size(mask_cell{1},1);
dens_vec = zeros(F,1);
for f = 1:F
    M = mask_cell{f};
    M(1:p+1:end) = false;
    edges = nnz(M)/2;
    dens_vec(f) = edges / (p*(p-1)/2);
end
dens_med = median(dens_vec);
dens_min = min(dens_vec);
dens_max = max(dens_vec);
end

% Local eLORETA Core (Minimal Implementation)
function [T, W] = run_eloreta_core_(L, Svv, regu)
[nchan, ndum] = size(L);
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
    if norm(diag(W)-diag(W_old))/(norm(diag(W_old))+1e-12) < 1e-3, break; end
end
I_n = eye(nchan, 'like', L);
K_final = (L * W) * L.';
alpha = regu * trace(K_final) / nchan;
T = W * (L.' * ((K_final + alpha * I_n) \ I_n));
end

function val = get_cfg(s, f, d)
if isfield(s, f), val = s.(f); else, val = d; end
end

function varargout = prepare_eta_bounds_(cfg, eta_global_lb, eta_global_ub)
% PREPARE_ETA_BOUNDS_  Build eta bounds & eta0 in log space with optional shrink.
%
% Usage examples (all supported):
%   info = prepare_eta_bounds_(cfg);
%   info = prepare_eta_bounds_(cfg, eta_global_lb, eta_global_ub);
%   [lb, ub, lb_log, ub_log, x0_log, eta0] = prepare_eta_bounds_(cfg, eta_global_lb, eta_global_ub);
%
% Reads from cfg (if exists):
%   cfg.eta1_lb/ub cfg.eta2_lb/ub cfg.eta3_lb/ub   (global bounds override)
%   cfg.eta_lb / cfg.eta_ub   (3-vector bounds override)
%   cfg.eta1_init cfg.eta2_init cfg.eta3_init or cfg.eta0 (3-vector)
%   cfg.eta_prev (3-vector)   cfg.eta_bank (Kx3)
%   cfg.eta_shrink_enable (default true)
%   cfg.eta_shrink_decades (default 0.7)
%   cfg.eta_shrink_mode: 'prev' | 'bank_range' (default 'prev')

if nargin < 1 || isempty(cfg), cfg = struct(); end

% -------- default global bounds (if not provided) --------
if nargin < 2 || isempty(eta_global_lb)
    eta_global_lb = [1e-3, 1e-3, 1e-3];
end
if nargin < 3 || isempty(eta_global_ub)
    eta_global_ub = [0.5, 1.0, 0.5];
end

% Support vector bounds if provided in cfg
if isfield(cfg,'eta_lb') && numel(cfg.eta_lb)==3
    eta_global_lb = cfg.eta_lb(:).';
end
if isfield(cfg,'eta_ub') && numel(cfg.eta_ub)==3
    eta_global_ub = cfg.eta_ub(:).';
end

% Support scalar-per-dimension bounds if provided in cfg
if isfield(cfg,'eta1_lb'), eta_global_lb(1) = cfg.eta1_lb; end
if isfield(cfg,'eta2_lb'), eta_global_lb(2) = cfg.eta2_lb; end
if isfield(cfg,'eta3_lb'), eta_global_lb(3) = cfg.eta3_lb; end
if isfield(cfg,'eta1_ub'), eta_global_ub(1) = cfg.eta1_ub; end
if isfield(cfg,'eta2_ub'), eta_global_ub(2) = cfg.eta2_ub; end
if isfield(cfg,'eta3_ub'), eta_global_ub(3) = cfg.eta3_ub; end

% -------- eta0 --------
eta0 = [get_cfg_(cfg,'eta1_init',0.1), get_cfg_(cfg,'eta2_init',0.2), get_cfg_(cfg,'eta3_init',0.1)];
if isfield(cfg,'eta0') && numel(cfg.eta0)==3
    eta0 = cfg.eta0(:).';
end
eta0 = max(eta0, eta_global_lb);
eta0 = min(eta0, eta_global_ub);

% -------- prev/bank --------
eta_prev = [];
if isfield(cfg,'eta_prev') && ~isempty(cfg.eta_prev) && numel(cfg.eta_prev)==3
    eta_prev = cfg.eta_prev(:).';
end
eta_bank = [];
if isfield(cfg,'eta_bank') && ~isempty(cfg.eta_bank) && size(cfg.eta_bank,2)==3
    eta_bank = cfg.eta_bank;
end

% -------- shrink options --------
DO_SHRINK = get_cfg_(cfg,'eta_shrink_enable', true);
SH_DEC    = get_cfg_(cfg,'eta_shrink_decades', 0.7);          % +/- decades in log10
SH_MODE   = lower(get_cfg_(cfg,'eta_shrink_mode','prev'));    % 'prev'|'bank_range'

% Start from global bounds
lb = eta_global_lb;
ub = eta_global_ub;

if DO_SHRINK
    if strcmp(SH_MODE,'bank_range') && ~isempty(eta_bank)
        lb = min(eta_bank,[],1) ./ (10.^SH_DEC);
        ub = max(eta_bank,[],1) .* (10.^SH_DEC);
    else
        if ~isempty(eta_prev), c = eta_prev; else, c = eta0; end
        lb = c ./ (10.^SH_DEC);
        ub = c .* (10.^SH_DEC);
    end
    % clamp to global bounds
    lb = max(lb, eta_global_lb);
    ub = min(ub, eta_global_ub);

    % safety: avoid inversion
    if any(lb >= ub)
        lb = eta_global_lb;
        ub = eta_global_ub;
    end
end

% clamp eta0 again into final bounds
eta0 = max(eta0, lb);
eta0 = min(eta0, ub);

lb_log = log10(lb);
ub_log = log10(ub);
x0_log = log10(eta0);

info = struct();
info.eta0 = eta0;
info.lb = lb;
info.ub = ub;
info.lb_log = lb_log;
info.ub_log = ub_log;
info.x0_log = x0_log;
info.eta_global_lb = eta_global_lb;
info.eta_global_ub = eta_global_ub;
info.eta_prev = eta_prev;
info.eta_bank = eta_bank;
info.do_shrink = DO_SHRINK;
info.shrink_decades = SH_DEC;
info.shrink_mode = SH_MODE;

% -------- flexible outputs --------
if nargout <= 1
    varargout{1} = info;
    return;
end

% Common multi-output orders:
% [lb, ub, lb_log, ub_log, x0_log, eta0]
varargout{1} = lb;
if nargout >= 2, varargout{2} = ub; end
if nargout >= 3, varargout{3} = lb_log; end
if nargout >= 4, varargout{4} = ub_log; end
if nargout >= 5, varargout{5} = x0_log; end
if nargout >= 6, varargout{6} = eta0; end
if nargout >= 7, varargout{7} = info; end
end

function val = get_cfg_(s, f, d)
if isfield(s, f), val = s.(f); else, val = d; end
end

function packed = pack_mask_cell_lower_(mask_cell)
%PACK_MASK_CELL_LOWER Pack logical masks as lower-triangular linear indices.
%   packed{f} contains uint32 indices into the p-by-p matrix for edges
%   in tril(mask,-1). Diagonal entries are not stored.

F = numel(mask_cell);
packed = cell(F, 1);

for f = 1:F
    M = logical(mask_cell{f});
    p = size(M, 1);
    M(1:p+1:end) = false;
    idx = find(tril(M, -1));
    packed{f} = uint32(idx);
end
end