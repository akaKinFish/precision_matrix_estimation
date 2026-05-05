
function [Omega_est, Sigma_src_est, outs] = solver_jspace_3d_opt_stoch_v2(Svv_cell, L, GraphLaplacian, cfg)
% SOLVER_JSPACE_3D_OPT_STOCH_V2
% Experimental J-SPACE variant:
%   1) Search eta in log-space by surrogateopt;
%   2) Full EM recomputes scale-adaptive lambdas and active masks EACH iteration;
%   3) Debiasing is applied INSIDE the EM loop after every M-step;
%   4) Final output is the last internally-debiased, recolored precision.
%
% Key cfg fields added/used:
%   cfg.internal_debias               (default true)
%   cfg.debias_blend                  (default 1.0, in [0,1])
%   cfg.surrogate_mode                'mini_em' or 'mstep1' (default 'mini_em')
%   cfg.surrogate_em_iter             (default 2; ignored when surrogate_mode='mstep1')
%   cfg.surrogate_internal_debias     (default = internal_debias)
%   cfg.warm_start_from_search        (default true)
%   cfg.warm_start_sigma_from_search  (default false)
%
% External dependencies:
%   - module2_estep
%   - module1_data_whitening
%   - module5_stoch_fista_main
%   - module_debias
%   - module8_recoloring
%   - utils_math.* (project_spd, make_hermitian, safe_log_det)
%   - utils_stats_gmm_threshold_1d
%
% Notes:
%   - This version intentionally removes the old post-EM Rayleigh/refit branch.
%   - If internal debias becomes unstable, try cfg.debias_blend = 0.3 ~ 0.7.
%
if nargin < 4, cfg = struct(); end

% ============================================================
% 1) Input normalization & config
% ============================================================
if ~iscell(Svv_cell)
    if ndims(Svv_cell) == 3
        F_in = size(Svv_cell, 3);
        tmp = cell(F_in, 1);
        for f = 1:F_in
            tmp{f} = Svv_cell(:, :, f);
        end
        Svv_cell = tmp;
    else
        Svv_cell = {Svv_cell};
    end
end

F = numel(Svv_cell);
[Ns, Nr] = size(L);

for f = 1:F
    Svv_cell{f} = sanitize_hermitian_cov_(Svv_cell{f});
end

% Global config
MAX_EM_ITER        = get_cfg(cfg, 'max_em_iter', 10);
VERBOSE            = get_cfg(cfg, 'verbose', true);
USE_GPU            = get_cfg(cfg, 'use_gpu', false);

% Surrogate config
MAX_OPT_EVALS      = get_cfg(cfg, 'opt_max_evals', 50);
MIN_OPT_POINTS     = get_cfg(cfg, 'opt_min_points', 10);
USE_PARALLEL       = get_cfg(cfg, 'opt_use_parallel', true);
SURR_MODE          = lower(get_cfg(cfg, 'surrogate_mode', 'mini_em'));
SURR_EM_ITER       = round(get_cfg(cfg, 'surrogate_em_iter', 2));
if any(strcmp(SURR_MODE, {'mstep1','single_mstep','truncated_mstep','warm_mstep'}))
    SURR_MODE = 'mstep1';
    SURR_EM_ITER = 1;
else
    SURR_MODE = 'mini_em';
    SURR_EM_ITER = max(1, SURR_EM_ITER);
end

% Solver & penalty config
OBJ_STOCH_MAX_ITER = get_cfg(cfg, 'obj_stoch_max_iter', 10);
EM_STOCH_MAX_ITER  = get_cfg(cfg, 'em_stoch_max_iter', 30);
PCOR_EPS           = get_cfg(cfg, 'dens_pcor_eps', 3e-3);
DENS_PEN_W         = get_cfg(cfg, 'density_penalty_weight', 1e7);
MASK_DENS_FLOOR    = get_cfg(cfg, 'mask_density_floor', max(0.03, log(max(Nr,2))/max(Nr,2)));
MASK_UNION_MODE    = get_cfg(cfg, 'mask_union', true);
M_SAMPLES          = get_cfg(cfg, 'm_samples', 100 * Nr);
UPDATE_RATE        = get_cfg(cfg, 'update_rate', 0.35);
UPDATE_RATE        = min(max(UPDATE_RATE, 0), 1);

% Experimental switches
INTERNAL_DEBIAS            = get_cfg(cfg, 'internal_debias', true);
DEBIAS_BLEND               = get_cfg(cfg, 'debias_blend', 1.0);
DEBIAS_BLEND               = min(max(DEBIAS_BLEND, 0), 1);
SURR_INTERNAL_DEBIAS       = get_cfg(cfg, 'surrogate_internal_debias', INTERNAL_DEBIAS);
WARM_START_FROM_SEARCH     = get_cfg(cfg, 'warm_start_from_search', true);
WARM_START_SIGMA_FROM_SEARCH = get_cfg(cfg, 'warm_start_sigma_from_search', false);

if get_cfg(cfg, 'postprocess_enable', false) && VERBOSE
    fprintf('[Post] postprocess_enable is ignored in this internal-debias EM solver.\n');
end

% Stochastic solver config
stoch_cfg = get_stoch_cfg_(cfg);

% ============================================================
% 2) Pre-computation (Laplacian, kernel, weights)
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
dwi_mask = [];
if isfield(cfg, 'dwi_C') && ~isempty(cfg.dwi_C)
    if VERBOSE, fprintf('[J-SPACE-3D-STOCH] Using DWI soft prior.\n'); end
    [W_gamma, dwi_mask] = build_dwi_soft_prior_(cfg.dwi_C, Nr, cfg);
end

% Noise covariance
tr_S = 0;
for f = 1:F
    tr_S = tr_S + real(trace(Svv_cell{f}));
end
noise_scale = (tr_S / max(F * Ns, 1)) * 0.05;
noise_cov = noise_scale * eye(Ns);

% GPU setup (used only for the E-step/whitening path of the full EM loop)
if USE_GPU
    try
        L = gpuArray(L);
        noise_cov = gpuArray(noise_cov);
        GraphLaplacian = gpuArray(GraphLaplacian);
        for f = 1:F
            Svv_cell{f} = gpuArray(Svv_cell{f});
        end
        if VERBOSE, fprintf('[J-SPACE-3D-STOCH] GPU acceleration enabled for E-step/whitening.\n'); end
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
    S_avg = zeros(Ns, Ns, 'like', Svv_cell{1});
    for f = 1:F
        S_avg = S_avg + Svv_cell{f};
    end
    S_avg = S_avg / F;
    [T_eloreta, ~] = run_eloreta_core_(L, S_avg, 0.05);

    Sigma_source_curr = cell(F, 1);
    for f = 1:F
        Sjj_eloreta = utils_math.make_hermitian(T_eloreta * Svv_cell{f} * T_eloreta');
        [Sjj_spd, ~] = utils_math.project_spd(Sjj_eloreta, 1e-8);
        Sigma_source_curr{f} = Sjj_spd;
    end
end

if numel(Sigma_source_curr) == 1 && F > 1
    Sigma_source_curr = repmat(Sigma_source_curr(:), F, 1);
elseif numel(Sigma_source_curr) ~= F
    error('init_sigma_source must be either a single matrix or a cell array with one entry per frequency.');
end

Sigma_source_curr = sanitize_covariance_list_(Sigma_source_curr, 1e-8);
if USE_GPU
    Sigma_source_curr = maybe_to_device_cell_(gather_cell_(Sigma_source_curr), true);
else
    Sigma_source_curr = gather_cell_(Sigma_source_curr);
end

Sigma_source_init_cpu = gather_cell_(Sigma_source_curr);
Gamma_warm_start = make_identity_cell_(F, Nr);

% ============================================================
% 4) Warm-start E-step + whitening (base data for search)
% ============================================================
if VERBOSE, fprintf('[J-SPACE-3D-STOCH] Preparing warm-start data for eta search...\n'); end

[Psijj_cell0, ~] = module2_estep(Svv_cell, L, Sigma_source_curr, noise_cov);
[Sjj_tilde0, D0_cell, ~] = module1_data_whitening(Psijj_cell0, 'smoothing_window', 1);

Sjj_tilde0_cpu = gather_cell_(Sjj_tilde0);
D0_cell_cpu    = gather_cell_(D0_cell);
K_cpu          = gather(K_freq);
W_cpu          = gather(W_gamma);
GraphL_cpu     = gather(GraphLaplacian);
Svv_cpu        = gather_cell_(Svv_cell);
L_cpu          = gather(L);
noise_cov_cpu  = gather(noise_cov);

% M-step static inputs (CPU side)
m5_in_base = struct();
m5_in_base.whitened_covariances = Sjj_tilde0_cpu;
m5_in_base.smoothing_kernel     = K_cpu;
m5_in_base.weight_matrix        = W_cpu;

m5_p_base = struct();
m5_p_base.alpha0 = 0.1;
m5_p_base.min_eig = get_cfg(cfg, 'min_eig', 1e-6);
m5_p_base.spatial_graph_matrix = GraphL_cpu;
m5_p_base.spatial_graph_is_laplacian = true;
m5_p_base.verbose = false;
m5_p_base.weight_mode = 'hadamard';
m5_p_base.auto_tune = false;
m5_p_base.stoch = stoch_cfg;
m5_p_base.stoch.max_iter = OBJ_STOCH_MAX_ITER;

mstep_solver = @module5_stoch_fista_main;

% ============================================================
% 5) Warm-start thresholds/scales and eta bounds
% ============================================================
if VERBOSE, fprintf('[J-SPACE-3D-STOCH] Computing warm-start scales/thresholds...\n'); end

[sc0, thr0] = compute_scales_thresholds_(m5_in_base.whitened_covariances, ...
    m5_in_base.smoothing_kernel, m5_in_base.weight_matrix, Nr, cfg);

min_den_limit = log(max(Nr,2)) / max(Nr,2);
max_den_limit = min(thr0.data_density, get_cfg(cfg, 'max_density_cap', 0.40));
max_den_limit = max(max_den_limit, min_den_limit * 1.5);

[mask0, mask0_stats] = build_em_active_mask_(m5_in_base.whitened_covariances, ...
    thr0.t_active, MASK_DENS_FLOOR, dwi_mask, [], false);
m5_in_base.active_mask = mask0;

eta_bounds = prepare_eta_bounds_(cfg);
initPts = eta_bounds.x0_log;
if isfield(cfg, 'opt_initial_points_log') && ~isempty(cfg.opt_initial_points_log)
    initPts = cfg.opt_initial_points_log;
end

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
% 6) Surrogate optimization
% ============================================================
opt_seed_base = get_cfg(cfg, 'opt_rng_seed', 0);

sur_ctx = struct();
sur_ctx.Svv_cell              = Svv_cpu;
sur_ctx.L                     = L_cpu;
sur_ctx.noise_cov             = noise_cov_cpu;
sur_ctx.Sigma_source_init     = Sigma_source_init_cpu;
sur_ctx.whitened_covariances0 = Sjj_tilde0_cpu;
sur_ctx.D0                    = D0_cell_cpu;
sur_ctx.mask0                 = mask0;
sur_ctx.sc0                   = sc0;
sur_ctx.thr0                  = thr0;
sur_ctx.smoothing_kernel      = K_cpu;
sur_ctx.weight_matrix         = W_cpu;
sur_ctx.dwi_mask              = dwi_mask;
sur_ctx.Nr                    = Nr;
sur_ctx.F                     = F;
sur_ctx.cfg                   = cfg;
sur_ctx.update_rate           = UPDATE_RATE;
sur_ctx.mask_density_floor    = MASK_DENS_FLOOR;
sur_ctx.mask_union_mode       = MASK_UNION_MODE;
sur_ctx.pcor_eps              = PCOR_EPS;
sur_ctx.min_eig               = m5_p_base.min_eig;
sur_ctx.internal_debias       = SURR_INTERNAL_DEBIAS;
sur_ctx.debias_blend          = DEBIAS_BLEND;
sur_ctx.surrogate_em_iter     = SURR_EM_ITER;
sur_ctx.obj_stoch_max_iter    = OBJ_STOCH_MAX_ITER;
sur_ctx.M_samples             = M_SAMPLES;
sur_ctx.Gamma_warm_start0     = make_identity_cell_(F, Nr);

obj_fun = @(x_log) wrapper_objective_stoch_(x_log, sur_ctx, m5_p_base, ...
    M_SAMPLES, min_den_limit, max_den_limit, DENS_PEN_W, mstep_solver, opt_seed_base);

if VERBOSE
    xchk = eta_bounds.x0_log(:).';
    s1 = obj_fun(xchk);
    s2 = obj_fun(xchk);
    fprintf('[obj determinism] f(x0) run1=%.6e run2=%.6e | rel_diff=%.3e\n', ...
        s1, s2, abs(s1-s2)/max(1, abs(s1)));
end

if VERBOSE
    if strcmp(SURR_MODE, 'mstep1')
        fprintf('[J-SPACE-3D-STOCH] Search mode: one warm-start truncated M-step.\n');
    else
        fprintf('[J-SPACE-3D-STOCH] Search mode: mini-EM surrogate (%d iterations).\n', SURR_EM_ITER);
    end
end

tic_opt = tic;
if exist('surrogateopt', 'file') == 2
    opts = optimoptions('surrogateopt', ...
        'MaxFunctionEvaluations', MAX_OPT_EVALS, ...
        'MinSurrogatePoints', MIN_OPT_POINTS, ...
        'UseParallel', USE_PARALLEL, ...
        'PlotFcn', [], ...
        'InitialPoints', initPts);

    [x_best_log, fval, exitflag, output_opt, trials] = surrogateopt( ...
        obj_fun, eta_bounds.lb_log, eta_bounds.ub_log, opts);
else
    if VERBOSE
        fprintf('[surrogateopt] Not found. Falling back to evaluating initial points only.\n');
    end
    X = initPts;
    if isempty(X), X = eta_bounds.x0_log; end
    fv = inf(size(X,1),1);
    for ii = 1:size(X,1)
        fv(ii) = obj_fun(X(ii,:));
    end
    [fval, idx_best] = min(fv);
    x_best_log = X(idx_best, :);
    exitflag = 0;
    output_opt = struct('message', 'surrogateopt not found; evaluated initial points only.');
    trials = struct('X', X, 'Fval', fv);
end
time_opt = toc(tic_opt);

try
    fv = trials.Fval(:);
    nfv = numel(fv);
    k0 = max(1, nfv - 9);
    best_all = min(fv);
    best_last10 = min(fv(k0:end));
    if VERBOSE
        fprintf('[surrogateopt trials] evals=%d | best_all=%.3e | best_last10=%.3e\n', ...
            nfv, best_all, best_last10);
        if best_last10 < 0.99 * best_all
            fprintf('[surrogateopt trials] WARNING: last10 best is still improving -> consider larger opt_max_evals.\n');
        end
    end
catch ME
    if VERBOSE
        fprintf('[surrogateopt trials] analysis failed: %s\n', ME.message);
    end
end

eta_best = 10.^x_best_log;
lb = eta_bounds.lb(:).';
ub = eta_bounds.ub(:).';
eta = eta_best(:).';
tol_rel = 1e-3;
hit_lb = abs(eta - lb) ./ max(abs(lb), eps) < tol_rel;
hit_ub = abs(eta - ub) ./ max(abs(ub), eps) < tol_rel;

if VERBOSE
    fprintf('[eta best] eta=[%.3e %.3e %.3e]\n', eta(1), eta(2), eta(3));
    fprintf('[eta hit]  hit_lb=[%d %d %d] hit_ub=[%d %d %d]\n', ...
        hit_lb(1), hit_lb(2), hit_lb(3), hit_ub(1), hit_ub(2), hit_ub(3));
    fprintf('  > Search finished in %.2fs. Best score: %.3e\n', time_opt, fval);
end

outs = struct();
outs.eta_bounds = eta_bounds;
outs.opt_results.fval = fval;
outs.opt_results.exitflag = exitflag;
outs.opt_results.output = output_opt;
outs.opt_results.trials = trials;
if isfield(trials, 'X')
    outs.opt_results.trials.X_log = trials.X;
    outs.opt_results.trials.X = 10.^trials.X;
end
outs.best_eta.eta1 = eta_best(1);
outs.best_eta.eta2 = eta_best(2);
outs.best_eta.eta3 = eta_best(3);
outs.best_eta.hit_lb = hit_lb;
outs.best_eta.hit_ub = hit_ub;
outs.thresholds.warm = thr0;
outs.thresholds.warm_mask_density = mask0_stats.density;
outs.density_limits.min = min_den_limit;
outs.density_limits.max = max_den_limit;
outs.surrogate.mode = SURR_MODE;
outs.surrogate.em_iter = SURR_EM_ITER;
outs.surrogate.internal_debias = SURR_INTERNAL_DEBIAS;

% Optional replay of the best surrogate path to build a better warm start
prev_mask_em = mask0;
if WARM_START_FROM_SEARCH
    try
        if VERBOSE, fprintf('[Warm Start] Replaying best surrogate path...\n'); end
        best_path = run_surrogate_path_(eta_best, sur_ctx, m5_p_base, mstep_solver, opt_seed_base);
        Gamma_warm_start = best_path.Gamma_final;
        prev_mask_em = best_path.mask_last;
        if WARM_START_SIGMA_FROM_SEARCH
            Sigma_source_curr = maybe_to_device_cell_(best_path.Sigma_final, USE_GPU);
        end
        outs.search_warm_start.used = true;
        outs.search_warm_start.lambda_last = best_path.lambdas_last;
        outs.search_warm_start.mask_density_last = best_path.mask_density_last;
        outs.search_warm_start.graph_density_last = best_path.graph_density_last;
        outs.search_warm_start.threshold_last = best_path.thr_last.t_active;
    catch ME
        if VERBOSE
            fprintf('[Warm Start] Replay failed. Using identity warm start. Reason: %s\n', ME.message);
        end
        Gamma_warm_start = make_identity_cell_(F, Nr);
        prev_mask_em = mask0;
        outs.search_warm_start.used = false;
        outs.search_warm_start.error = ME.message;
    end
else
    outs.search_warm_start.used = false;
end

% ============================================================
% 7) Full EM with INTERNAL debias
% ============================================================
if VERBOSE
    fprintf('\n[J-SPACE-3D-STOCH] Running full EM with INTERNAL debias...\n');
end

outs.loglik           = nan(MAX_EM_ITER, 1);
outs.em_lambdas       = nan(MAX_EM_ITER, 3);
outs.em_scales        = nan(MAX_EM_ITER, 3);
outs.em_thresholds    = nan(MAX_EM_ITER, 1);
outs.em_mask_density  = nan(MAX_EM_ITER, 1);
outs.em_graph_density = nan(MAX_EM_ITER, 1);
outs.em_time          = nan(MAX_EM_ITER, 1);

Gamma_biased_last = [];
Gamma_used_last = [];
Omega_last = [];
Sjj_last = [];
D_last = [];
last_iter_done = 0;

m5_p_em = m5_p_base;
m5_p_em.stoch.max_iter = EM_STOCH_MAX_ITER;

for em_iter = 1:MAX_EM_ITER
    iter_tic = tic;

    % ---- A. E-step ----
    [Psijj_cell, e_stats] = module2_estep(Svv_cell, L, Sigma_source_curr, noise_cov);
    if isstruct(e_stats) && isfield(e_stats, 'log_likelihood')
        outs.loglik(em_iter) = e_stats.log_likelihood;
    else
        outs.loglik(em_iter) = NaN;
    end

    % ---- B. Whitening ----
    [Sjj_tilde, D_cell, ~] = module1_data_whitening(Psijj_cell, 'smoothing_window', 1);
    Sjj_tilde_cpu = gather_cell_(Sjj_tilde);
    D_cell_cpu    = gather_cell_(D_cell);
    Sigma_curr_cpu = gather_cell_(Sigma_source_curr);

    % ---- C. Dynamic scales + dynamic active mask ----
    [sc_k, thr_k] = compute_scales_thresholds_(Sjj_tilde_cpu, K_cpu, W_cpu, Nr, cfg);
    [mask_k, mask_stats] = build_em_active_mask_(Sjj_tilde_cpu, thr_k.t_active, ...
        MASK_DENS_FLOOR, dwi_mask, prev_mask_em, MASK_UNION_MODE);
    [lambda1_k, lambda2_k, lambda3_k] = eta_to_lambdas_(eta_best, sc_k);

    outs.em_scales(em_iter, :) = [sc_k.lambda1_scale, sc_k.lambda2_max, sc_k.lambda3_scale];
    outs.em_lambdas(em_iter, :) = [lambda1_k, lambda2_k, lambda3_k];
    outs.em_thresholds(em_iter) = thr_k.t_active;
    outs.em_mask_density(em_iter) = mask_stats.density;

    % ---- D. M-step ----
    local_in = struct();
    local_in.whitened_covariances = Sjj_tilde_cpu;
    local_in.smoothing_kernel     = K_cpu;
    local_in.weight_matrix        = W_cpu;
    local_in.active_mask          = mask_k;
    local_in.precision_matrices   = gather_cell_(Gamma_warm_start);

    m5_p_em.lambda1 = lambda1_k;
    m5_p_em.lambda2 = lambda2_k;
    m5_p_em.lambda3 = lambda3_k;
    m5_p_em.stoch.seed = stoch_cfg.seed + em_iter;
    m5_p_em.stoch.max_iter = EM_STOCH_MAX_ITER;

    try
        [Gamma_new, ~] = mstep_solver(local_in, m5_p_em);
    catch ME
        if VERBOSE
            fprintf('[EM %d] M-step failed, fallback to warm start. Reason: %s\n', em_iter, ME.message);
        end
        Gamma_new = local_in.precision_matrices;
    end
    Gamma_new = sanitize_precision_list_(Gamma_new, m5_p_em.min_eig);

    % ---- E. Internal debias ----
    if INTERNAL_DEBIAS
        Gamma_db = debias_and_project_(Gamma_new, Sjj_tilde_cpu, M_SAMPLES, m5_p_em.min_eig);
        if DEBIAS_BLEND < 1
            Gamma_used = blend_precision_lists_(Gamma_new, Gamma_db, DEBIAS_BLEND, m5_p_em.min_eig);
        else
            Gamma_used = Gamma_db;
        end
    else
        Gamma_used = Gamma_new;
    end
    Gamma_warm_start = Gamma_used;

    % ---- F. Recolor + inertia update ----
    m8_in = struct();
    m8_in.whitened_precision_matrices = Gamma_used;
    m8_in.whitening_matrices          = D_cell_cpu;
    m8_in.original_covariances        = Sigma_curr_cpu;

    recol = module8_recoloring(m8_in, struct('verbose', false));
    Omega_temp = sanitize_precision_list_(recol.recolored_precision_matrices, 1e-8);

    Sigma_next_cpu = cell(F, 1);
    for f = 1:F
        Om = utils_math.make_hermitian(Omega_temp{f});
        [Om_spd, ~] = utils_math.project_spd(Om, 1e-8);
        S_cand = pinv(Om_spd);
        S_new = (1 - UPDATE_RATE) * Sigma_curr_cpu{f} + UPDATE_RATE * S_cand;
        S_new = utils_math.make_hermitian(S_new);
        [S_new, ~] = utils_math.project_spd(S_new, 1e-8);
        Sigma_next_cpu{f} = S_new;
    end
    Sigma_source_curr = maybe_to_device_cell_(Sigma_next_cpu, USE_GPU);

    % ---- G. Diagnostics ----
    [den_graph, ~] = density_from_G_pcor_(Gamma_used, PCOR_EPS);
    outs.em_graph_density(em_iter) = den_graph;
    outs.em_time(em_iter) = toc(iter_tic);

    if VERBOSE
        fprintf('  > EM %d/%d: LogLik=%.3e | lambda=[%.3e %.3e %.3e] | mask=%.2f%% | graph=%.2f%% | time=%.2fs\n', ...
            em_iter, MAX_EM_ITER, outs.loglik(em_iter), ...
            lambda1_k, lambda2_k, lambda3_k, ...
            100 * outs.em_mask_density(em_iter), 100 * outs.em_graph_density(em_iter), outs.em_time(em_iter));
    end

    prev_mask_em = mask_k;
    Gamma_biased_last = Gamma_new;
    Gamma_used_last = Gamma_used;
    Omega_last = Omega_temp;
    Sjj_last = Sjj_tilde_cpu;
    D_last = D_cell_cpu;
    last_iter_done = em_iter;
end

% ============================================================
% 8) Return final internal-debias EM output
% ============================================================
if isempty(Omega_last)
    % Defensive fallback
    m8_in = struct();
    m8_in.whitened_precision_matrices = gather_cell_(Gamma_warm_start);
    m8_in.whitening_matrices          = D0_cell_cpu;
    m8_in.original_covariances        = Sigma_source_init_cpu;
    recol = module8_recoloring(m8_in, struct('verbose', false));
    Omega_last = sanitize_precision_list_(recol.recolored_precision_matrices, 1e-8);
    D_last = D0_cell_cpu;
    Sjj_last = Sjj_tilde0_cpu;
    Gamma_biased_last = gather_cell_(Gamma_warm_start);
    Gamma_used_last = gather_cell_(Gamma_warm_start);
end

Omega_est = Omega_last;
Sigma_src_est = gather_cell_(Sigma_source_curr);

outs.final_iter = last_iter_done;
outs.final_whitened_covariances = Sjj_last;
outs.final_whitening_matrices   = D_last;
outs.final_Gamma_biased         = Gamma_biased_last;
outs.final_Gamma_used           = Gamma_used_last;
outs.final_Omega                = Omega_est;
outs.internal_debias = INTERNAL_DEBIAS;
outs.debias_blend = DEBIAS_BLEND;
outs.update_rate = UPDATE_RATE;

if last_iter_done >= 1
    outs.global_hyperparams.eta1 = eta_best(1);
    outs.global_hyperparams.eta2 = eta_best(2);
    outs.global_hyperparams.eta3 = eta_best(3);
    outs.global_hyperparams.lambda1_final = outs.em_lambdas(last_iter_done, 1);
    outs.global_hyperparams.lambda2_final = outs.em_lambdas(last_iter_done, 2);
    outs.global_hyperparams.lambda3_final = outs.em_lambdas(last_iter_done, 3);
else
    outs.global_hyperparams.eta1 = eta_best(1);
    outs.global_hyperparams.eta2 = eta_best(2);
    outs.global_hyperparams.eta3 = eta_best(3);
    outs.global_hyperparams.lambda1_final = NaN;
    outs.global_hyperparams.lambda2_final = NaN;
    outs.global_hyperparams.lambda3_final = NaN;
end

end

% ============================================================
% LOCAL HELPERS
% ============================================================

function score = wrapper_objective_stoch_(x_log, sur_ctx, m5_p_base, M, min_den, max_den, penW, solver, seed_base)
eta = 10.^x_log(:).';
try
    aux = run_surrogate_path_(eta, sur_ctx, m5_p_base, solver, seed_base);

    G_eval = aux.Gamma_final;
    S_eval = aux.Sjj_final;

    [den_med, total_edges] = density_from_G_pcor_(G_eval, sur_ctx.pcor_eps);

    penalty = 0;
    if den_med < min_den
        penalty = penW * ((min_den - den_med) / max(min_den, eps))^2;
    elseif den_med > max_den
        penalty = penW * ((den_med - max_den) / max(max_den, eps))^2;
    end

    ld_sum = 0;
    tr_sum = 0;
    F = numel(G_eval);
    for ff = 1:F
        Gf = G_eval{ff};
        [ld_f, valid] = utils_math.safe_log_det(Gf);
        if ~valid || ~isfinite(ld_f)
            score = 1e15;
            return;
        end
        ld_sum = ld_sum + ld_f;
        tr_sum = tr_sum + real(trace(S_eval{ff} * Gf));
    end

    minus_2_ll = M * (tr_sum - ld_sum);
    aic_like = minus_2_ll + 2 * total_edges;
    score = aic_like + penalty;

    if isnan(score) || isinf(score)
        score = 1e15;
    end
catch
    score = 1e15;
end
end

function aux = run_surrogate_path_(eta, sur_ctx, m5_p_base, solver, seed_base)
% Runs either:
%   - one warm-start truncated M-step (surrogate_em_iter = 1), or
%   - a small multi-EM surrogate path (surrogate_em_iter > 1)
%
Sigma_local = sur_ctx.Sigma_source_init;
Gamma_local = sur_ctx.Gamma_warm_start0;
prev_mask = sur_ctx.mask0;

Sjj_cpu = sur_ctx.whitened_covariances0;
D_cpu   = sur_ctx.D0;
sc_k    = sur_ctx.sc0;
thr_k   = sur_ctx.thr0;

for kk = 1:sur_ctx.surrogate_em_iter
    if kk > 1
        [Psijj_k, ~] = module2_estep(sur_ctx.Svv_cell, sur_ctx.L, Sigma_local, sur_ctx.noise_cov);
        [Sjj_tilde_k, D_k, ~] = module1_data_whitening(Psijj_k, 'smoothing_window', 1);
        Sjj_cpu = gather_cell_(Sjj_tilde_k);
        D_cpu   = gather_cell_(D_k);
        [sc_k, thr_k] = compute_scales_thresholds_(Sjj_cpu, ...
            sur_ctx.smoothing_kernel, sur_ctx.weight_matrix, sur_ctx.Nr, sur_ctx.cfg);
    end

    [mask_k, mask_stats] = build_em_active_mask_(Sjj_cpu, thr_k.t_active, ...
        sur_ctx.mask_density_floor, sur_ctx.dwi_mask, prev_mask, sur_ctx.mask_union_mode);
    [lambda1_k, lambda2_k, lambda3_k] = eta_to_lambdas_(eta, sc_k);

    local_in = struct();
    local_in.whitened_covariances = Sjj_cpu;
    local_in.smoothing_kernel     = sur_ctx.smoothing_kernel;
    local_in.weight_matrix        = sur_ctx.weight_matrix;
    local_in.active_mask          = mask_k;
    local_in.precision_matrices   = Gamma_local;

    local_p = m5_p_base;
    local_p.lambda1 = lambda1_k;
    local_p.lambda2 = lambda2_k;
    local_p.lambda3 = lambda3_k;
    local_p.stoch.seed = seed_base + kk - 1;
    local_p.stoch.max_iter = sur_ctx.obj_stoch_max_iter;

    [Gamma_biased, ~] = solver(local_in, local_p);
    Gamma_biased = sanitize_precision_list_(Gamma_biased, local_p.min_eig);

    if sur_ctx.internal_debias
        Gamma_db = debias_and_project_(Gamma_biased, Sjj_cpu, sur_ctx.M_samples, local_p.min_eig);
        if sur_ctx.debias_blend < 1
            Gamma_used = blend_precision_lists_(Gamma_biased, Gamma_db, sur_ctx.debias_blend, local_p.min_eig);
        else
            Gamma_used = Gamma_db;
        end
    else
        Gamma_used = Gamma_biased;
    end

    Gamma_local = Gamma_used;

    m8_in = struct();
    m8_in.whitened_precision_matrices = Gamma_used;
    m8_in.whitening_matrices          = D_cpu;
    m8_in.original_covariances        = Sigma_local;
    recol = module8_recoloring(m8_in, struct('verbose', false));
    Omega_k = sanitize_precision_list_(recol.recolored_precision_matrices, 1e-8);

    Sigma_next = cell(numel(Sigma_local), 1);
    for f = 1:numel(Sigma_local)
        Om = utils_math.make_hermitian(Omega_k{f});
        [Om_spd, ~] = utils_math.project_spd(Om, 1e-8);
        S_cand = pinv(Om_spd);
        S_new = (1 - sur_ctx.update_rate) * Sigma_local{f} + sur_ctx.update_rate * S_cand;
        S_new = utils_math.make_hermitian(S_new);
        [S_new, ~] = utils_math.project_spd(S_new, 1e-8);
        Sigma_next{f} = S_new;
    end
    Sigma_local = Sigma_next;
    prev_mask = mask_k;
end

[graph_den, ~] = density_from_G_pcor_(Gamma_local, sur_ctx.pcor_eps);

aux = struct();
aux.Gamma_final = Gamma_local;
aux.Sigma_final = Sigma_local;
aux.Sjj_final = Sjj_cpu;
aux.D_final = D_cpu;
aux.mask_last = prev_mask;
aux.sc_last = sc_k;
aux.thr_last = thr_k;
aux.lambdas_last = [lambda1_k, lambda2_k, lambda3_k];
aux.mask_density_last = mask_stats.density;
aux.graph_density_last = graph_den;
end

function Gamma_out = debias_and_project_(Gamma_in, Sjj_cell, M_samples, min_eig)
try
    [Gamma_db, ~, ~, ~] = module_debias(Gamma_in, Sjj_cell, M_samples);
    Gamma_out = sanitize_precision_list_(Gamma_db, min_eig);
catch
    Gamma_out = sanitize_precision_list_(Gamma_in, min_eig);
end
end

function Gamma_out = blend_precision_lists_(Gamma_a, Gamma_b, alpha, min_eig)
alpha = min(max(alpha, 0), 1);
F = numel(Gamma_a);
Gamma_out = cell(F,1);
for f = 1:F
    G = (1 - alpha) * Gamma_a{f} + alpha * Gamma_b{f};
    G = utils_math.make_hermitian(G);
    [G_spd, ~] = utils_math.project_spd(G, min_eig);
    Gamma_out{f} = G_spd;
end
end

function Sigma_out = sanitize_covariance_list_(Sigma_in, min_eig)
F = numel(Sigma_in);
Sigma_out = cell(F,1);
for f = 1:F
    S = utils_math.make_hermitian(Sigma_in{f});
    [S_spd, ~] = utils_math.project_spd(S, min_eig);
    Sigma_out{f} = S_spd;
end
end

function Gamma_out = sanitize_precision_list_(Gamma_in, min_eig)
F = numel(Gamma_in);
Gamma_out = cell(F,1);
for f = 1:F
    G = utils_math.make_hermitian(Gamma_in{f});
    [G_spd, ~] = utils_math.project_spd(G, min_eig);
    Gamma_out{f} = G_spd;
end
end

function [mask_k, stats] = build_em_active_mask_(Sjj_cell, thresh, dens_floor, dwi_mask, prev_mask, union_mode)
[mask_k, stats0] = threshold_active_mask_with_floor_(Sjj_cell, thresh, dens_floor);
F = numel(mask_k);
p = size(mask_k{1}, 1);
stats = stats0;
stats.num_active_edges_final = zeros(F, 1);
for f = 1:F
    if ~isempty(dwi_mask)
        mask_k{f} = mask_k{f} | dwi_mask;
    end
    if union_mode && ~isempty(prev_mask)
        mask_k{f} = mask_k{f} | prev_mask{f};
    end
    mask_k{f}(1:p+1:end) = true;
    stats.num_active_edges_final(f) = (nnz(mask_k{f}) - p) / 2;
end
stats.density = mean(stats.num_active_edges_final) / max(p*(p-1)/2, 1);
end

function [lambda1, lambda2, lambda3] = eta_to_lambdas_(eta, sc)
lambda1 = eta(1) * sc.lambda1_scale;
lambda2 = eta(2) * sc.lambda2_max;
lambda3 = 0;
if sc.lk_norm > eps
    lambda3 = eta(3) * sc.lambda3_scale / (sc.lk_norm + eps);
end

if ~isfinite(lambda1), lambda1 = 0; end
if ~isfinite(lambda2), lambda2 = 0; end
if ~isfinite(lambda3), lambda3 = 0; end
end

function cell_out = gather_cell_(cell_in)
if ~iscell(cell_in)
    cell_in = {cell_in};
end
cell_out = cell(size(cell_in));
for i = 1:numel(cell_in)
    try
        cell_out{i} = gather(cell_in{i});
    catch
        cell_out{i} = cell_in{i};
    end
end
end

function cell_out = maybe_to_device_cell_(cell_in, use_gpu)
cell_out = cell(size(cell_in));
if use_gpu
    for i = 1:numel(cell_in)
        cell_out{i} = gpuArray(cell_in{i});
    end
else
    for i = 1:numel(cell_in)
        cell_out{i} = cell_in{i};
    end
end
end

function G0 = make_identity_cell_(F, p)
G0 = cell(F, 1);
for f = 1:F
    G0{f} = eye(p);
end
end

function S = sanitize_hermitian_cov_(S)
S = utils_math.make_hermitian(S);
d = real(diag(S));
d(~isfinite(d)) = 0;
d = max(d, 0);
n = size(S, 1);
S(1:n+1:end) = d;
end

function cfg_stoch = get_stoch_cfg_(cfg)
stoch = struct();
if isfield(cfg,'stoch') && ~isempty(cfg.stoch)
    stoch = cfg.stoch;
end
stoch.mode = get_cfg(stoch, 'mode', 'B');
stoch.seed = get_cfg(stoch, 'seed', get_cfg(cfg, 'opt_rng_seed', 0));
stoch.neighbor_closure = get_cfg(stoch, 'neighbor_closure', 1);
stoch.use_backtracking = get_cfg(stoch, 'use_backtracking', true);
stoch.backtracking_factor = get_cfg(stoch, 'backtracking_factor', 2.0);
stoch.max_backtracking = get_cfg(stoch, 'max_backtracking', 25);

stoch.use_nesterov = get_cfg(stoch, 'use_nesterov', false);
stoch.L_min = get_cfg(stoch, 'L_min', 1e-3);
stoch.L_max = get_cfg(stoch, 'L_max', 1e6);
stoch.L_shrink_on_success = get_cfg(stoch, 'L_shrink_on_success', true);

if isfield(cfg, 'freq')
    stoch.freq = cfg.freq;
end
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
    dens(f) = edges / max(p*(p-1)/2, 1);
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
stats.density = mean(stats.num_active_edges) / max(p*(p-1)/2, 1);
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
stats.density = mean(stats.num_active_edges) / max(p*(p-1)/2, 1);
stats.t_used = t_use;
end

function [sc, thr] = compute_scales_thresholds_(Sjj_cpu_cell, K_freq_cpu, W_gamma_cpu, Nr, cfg)
if ~iscell(Sjj_cpu_cell), Sjj_cpu_cell = {Sjj_cpu_cell}; end
F = numel(Sjj_cpu_cell);
p = size(Sjj_cpu_cell{1}, 1);

Sbar = zeros(p, p);
for f = 1:F
    Sbar = Sbar + Sjj_cpu_cell{f};
end
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
    if isempty(ratios)
        lambda2_max = 0;
    else
        lambda2_max = max(ratios);
    end
end

lambda1_scale = norm(Sbar - I, 'fro') * F;
if ~isfinite(lambda1_scale), lambda1_scale = 0; end

if F <= 1
    lambda3_scale = 0;
else
    Gbar = zeros(p, p);
    for f = 1:F
        Gbar = Gbar + (Sjj_cpu_cell{f} - I);
    end
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

fallback_q = get_cfg(cfg, 't_active_fallback_quantile', 0.85);
if isempty(vals)
    t_active = 0;
    q95 = 0;
    data_density = 0.05;
else
    try
        [t_active, ~] = utils_stats_gmm_threshold_1d(vals);
        if ~isfinite(t_active) || t_active < 0
            error('Invalid GMM threshold.');
        end
    catch
        t_active = quantile(vals, fallback_q);
    end
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

if isfield(cfg, 't_active_override') && ~isempty(cfg.t_active_override)
    thr.t_active = cfg.t_active_override;
end
end

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
mode = get_cfg(cfg, 'dwi_weight_mode', 'power');
alpha = get_cfg(cfg, 'dwi_weight_alpha', 2.0);
w_floor = get_cfg(cfg, 'dwi_w_floor', 0.2);
q = get_cfg(cfg, 'dwi_mask_quantile', 0.90);

switch mode
    case 'power'
        base = (1 - C_norm).^alpha;
    case 'exp'
        base = exp(-alpha * C_norm);
    case 'linear'
        base = max(0, 1 - alpha * C_norm);
    otherwise
        error('Unknown dwi_weight_mode');
end

W_gamma = w_floor + (1 - w_floor) * base;
W_gamma(1:Nr+1:end) = 1;

dwi_mask = [];
maskL = tril(true(Nr), -1);
cvals = C_norm(maskL);
if ~isempty(cvals)
    tq = quantile(cvals, q);
    M = (C_norm >= tq);
    M = M | M';
    M(1:Nr+1:end) = true;
    dwi_mask = logical(M);
end
end

function K = build_default_freq_kernel_(F)
if F <= 1
    K = 1;
    return;
end
K = eye(F);
for i = 1:(F - 1)
    K(i, i+1) = 1;
    K(i+1, i) = 1;
end
end

function [dens_vec, dens_med, dens_min, dens_max] = mask_density_stats_(mask_cell)
F = numel(mask_cell);
p = size(mask_cell{1},1);
dens_vec = zeros(F,1);
for f = 1:F
    M = mask_cell{f};
    M(1:p+1:end) = false;
    edges = nnz(M)/2;
    dens_vec(f) = edges / max(p*(p-1)/2, 1);
end
dens_med = median(dens_vec);
dens_min = min(dens_vec);
dens_max = max(dens_vec);
end

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
    if norm(diag(W) - diag(W_old)) / (norm(diag(W_old)) + 1e-12) < 1e-3
        break;
    end
end
I_n = eye(nchan, 'like', L);
K_final = (L * W) * L.';
alpha = regu * trace(K_final) / nchan;
T = W * (L.' * ((K_final + alpha * I_n) \ I_n));
end

function val = get_cfg(s, f, d)
if isfield(s, f)
    val = s.(f);
else
    val = d;
end
end

function varargout = prepare_eta_bounds_(cfg, eta_global_lb, eta_global_ub)
if nargin < 1 || isempty(cfg), cfg = struct(); end

if nargin < 2 || isempty(eta_global_lb)
    eta_global_lb = [1e-3, 1e-3, 1e-3];
end
if nargin < 3 || isempty(eta_global_ub)
    eta_global_ub = [0.5, 1.0, 0.5];
end

if isfield(cfg, 'eta_lb') && numel(cfg.eta_lb) == 3
    eta_global_lb = cfg.eta_lb(:).';
end
if isfield(cfg, 'eta_ub') && numel(cfg.eta_ub) == 3
    eta_global_ub = cfg.eta_ub(:).';
end

if isfield(cfg, 'eta1_lb'), eta_global_lb(1) = cfg.eta1_lb; end
if isfield(cfg, 'eta2_lb'), eta_global_lb(2) = cfg.eta2_lb; end
if isfield(cfg, 'eta3_lb'), eta_global_lb(3) = cfg.eta3_lb; end
if isfield(cfg, 'eta1_ub'), eta_global_ub(1) = cfg.eta1_ub; end
if isfield(cfg, 'eta2_ub'), eta_global_ub(2) = cfg.eta2_ub; end
if isfield(cfg, 'eta3_ub'), eta_global_ub(3) = cfg.eta3_ub; end

eta0 = [get_cfg_(cfg, 'eta1_init', 0.1), get_cfg_(cfg, 'eta2_init', 0.2), get_cfg_(cfg, 'eta3_init', 0.1)];
if isfield(cfg, 'eta0') && numel(cfg.eta0) == 3
    eta0 = cfg.eta0(:).';
end
eta0 = max(eta0, eta_global_lb);
eta0 = min(eta0, eta_global_ub);

eta_prev = [];
if isfield(cfg, 'eta_prev') && ~isempty(cfg.eta_prev) && numel(cfg.eta_prev) == 3
    eta_prev = cfg.eta_prev(:).';
end
eta_bank = [];
if isfield(cfg, 'eta_bank') && ~isempty(cfg.eta_bank) && size(cfg.eta_bank, 2) == 3
    eta_bank = cfg.eta_bank;
end

DO_SHRINK = get_cfg_(cfg, 'eta_shrink_enable', true);
SH_DEC    = get_cfg_(cfg, 'eta_shrink_decades', 0.7);
SH_MODE   = lower(get_cfg_(cfg, 'eta_shrink_mode', 'prev'));

lb = eta_global_lb;
ub = eta_global_ub;

if DO_SHRINK
    if strcmp(SH_MODE, 'bank_range') && ~isempty(eta_bank)
        lb = min(eta_bank, [], 1) ./ (10.^SH_DEC);
        ub = max(eta_bank, [], 1) .* (10.^SH_DEC);
    else
        if ~isempty(eta_prev)
            c = eta_prev;
        else
            c = eta0;
        end
        lb = c ./ (10.^SH_DEC);
        ub = c .* (10.^SH_DEC);
    end

    lb = max(lb, eta_global_lb);
    ub = min(ub, eta_global_ub);

    if any(lb >= ub)
        lb = eta_global_lb;
        ub = eta_global_ub;
    end
end

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

if nargout <= 1
    varargout{1} = info;
    return;
end

varargout{1} = lb;
if nargout >= 2, varargout{2} = ub; end
if nargout >= 3, varargout{3} = lb_log; end
if nargout >= 4, varargout{4} = ub_log; end
if nargout >= 5, varargout{5} = x0_log; end
if nargout >= 6, varargout{6} = eta0; end
if nargout >= 7, varargout{7} = info; end
end

function val = get_cfg_(s, f, d)
if isfield(s, f)
    val = s.(f);
else
    val = d;
end
end

