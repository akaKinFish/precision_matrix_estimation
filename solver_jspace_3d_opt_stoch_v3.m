function [Omega_est, Sigma_src_est, outs] = solver_jspace_3d_opt_stoch_v3(Svv_cell, L, GraphLaplacian, cfg)
% SOLVER_JSPACE_3D_OPT_STOCH
% Fast factor-augmented JSPACE (rank-1 low-rank common-mode deflation)
%
% Core idea:
%   1) Keep the original outer EM: scalp -> posterior source second moment.
%   2) Keep the original whitening step.
%   3) Estimate ONE shared low-rank alpha-mode direction u in whitened space.
%   4) At each outer EM iteration, estimate per-frequency amplitudes alpha_f
%      cheaply and subtract alpha_f * u*u' from the whitened covariance.
%   5) Feed the residual covariances into the existing stochastic FISTA M-step.
%   6) Reconstruct:
%        - Sigma_src_est : full covariance = low-rank + residual covariance
%        - Omega_est     : direct precision from the residual process
%
% Important semantic change vs. the original solver:
%   - Sigma_src_est is the FULL covariance estimate used for spectrum/coherence.
%   - Omega_est is the DIRECT precision estimate used for partial coherence.
%   - In the factor-augmented model, they are generally NOT inverse pairs.
%
% External dependencies required on path:
%   - module2_estep
%   - module1_data_whitening
%   - module5_stoch_fista_main
%   - utils_math.project_spd, utils_math.make_hermitian, utils_math.safe_log_det
%   - utils_stats_gmm_threshold_1d
%   - module_debias    (optional, if postprocess_enable=true)
%
% Recommended usage:
%   cfg = struct();
%   cfg.freq = freq;                     % required for alpha-band detection
%   cfg.lowrank_enable = true;          % default true
%   cfg.lowrank_alpha_band = [8 12];    % default [8 12]
%   cfg.lowrank_alpha_smooth = 0.25;    % default 0.25
%   cfg.postprocess_enable = true;      % optional
%   cfg.debias_only = true;             % optional
%   [Omega_est, Sjj_est, outs] = solver_jspace_3d_opt_stoch(Svv_cross, L, [], cfg);
%
% Notes:
%   - This solver is designed to be a drop-in replacement for the current
%     stochastic solver in run_jspace_real / run_jspace_stoch_age_cohort.
%   - For compatibility with the wrappers, the second output remains named
%     Sigma_src_est, and can be saved as Sjj_est by the wrapper.
%
if nargin < 4, cfg = struct(); end

% ============================================================
% 1) Input normalization
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

% ============================================================
% 2) Config
% ============================================================
MAX_EM_ITER        = get_cfg(cfg, 'max_em_iter', 10);
VERBOSE            = get_cfg(cfg, 'verbose', true);
USE_GPU            = get_cfg(cfg, 'use_gpu', false);
DO_POST            = get_cfg(cfg, 'postprocess_enable', true);
DEBIAS_ONLY        = get_cfg(cfg, 'debias_only', true);

% Surrogate settings
MAX_OPT_EVALS      = get_cfg(cfg, 'opt_max_evals', 50);
MIN_OPT_POINTS     = get_cfg(cfg, 'opt_min_points', 10);
USE_PARALLEL       = get_cfg(cfg, 'opt_use_parallel', true);

% Stochastic M-step settings
OBJ_STOCH_MAX_ITER = get_cfg(cfg, 'obj_stoch_max_iter', 10);
EM_STOCH_MAX_ITER  = get_cfg(cfg, 'em_stoch_max_iter', 30);
PCOR_EPS           = get_cfg(cfg, 'dens_pcor_eps', 3e-3);
DENS_PEN_W         = get_cfg(cfg, 'density_penalty_weight', 1e7);
MASK_DENS_FLOOR    = get_cfg(cfg, 'mask_density_floor', max(0.03, log(Nr)/Nr));
MASK_UNION_MODE    = get_cfg(cfg, 'mask_union', true);
M_SAMPLES          = get_cfg(cfg, 'm_samples', 100 * Nr);
UPDATE_RATE        = get_cfg(cfg, 'update_rate', 0.35);
MIN_EIG            = get_cfg(cfg, 'min_eig', 1e-6);

% Factor-augmented low-rank part
LOWRANK_ENABLE       = get_cfg(cfg, 'lowrank_enable', true);
LOWRANK_ALPHA_BAND   = get_cfg(cfg, 'lowrank_alpha_band', [8 12]);   % mode-u estimation band
LOWRANK_ALPHA_SMOOTH = get_cfg(cfg, 'lowrank_alpha_smooth', 0.25);   % alpha_f smoothing across freq
LOWRANK_MIN_EIG      = get_cfg(cfg, 'lowrank_resid_min_eig', MIN_EIG);

% ===== NEW: adaptive alpha gate =====
LOWRANK_GATE_ENABLE        = get_cfg(cfg, 'lowrank_gate_enable', true);
LOWRANK_GATE_SEARCH_BAND   = get_cfg(cfg, 'lowrank_gate_search_band', [6 14]); % broad candidate band
LOWRANK_GATE_SMOOTH        = get_cfg(cfg, 'lowrank_gate_smooth', 0.15);        % smooth b_f^(0) before peak finding
LOWRANK_GATE_FRACTION      = get_cfg(cfg, 'lowrank_gate_fraction', 0.45);      % width determined at 45% of peak
LOWRANK_GATE_FLOOR         = get_cfg(cfg, 'lowrank_gate_floor', 0.02);         % keep tiny residual subtraction outside alpha
LOWRANK_GATE_MIN_WIDTH_HZ  = get_cfg(cfg, 'lowrank_gate_min_width_hz', 0.8);
LOWRANK_GATE_MAX_WIDTH_HZ  = get_cfg(cfg, 'lowrank_gate_max_width_hz', 3.0);

% Stochastic config helper
stoch_cfg = get_stoch_cfg_(cfg);

% ============================================================
% 3) Pre-computation (kernel, weights, DWI prior)
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
    if VERBOSE
        fprintf('[JSPACE-LR1] Using DWI soft prior.\n');
    end
    [W_gamma, dwi_mask] = build_dwi_soft_prior_(cfg.dwi_C, Nr, cfg);
end

% Noise covariance
tr_S = 0;
for f = 1:F
    tr_S = tr_S + trace(Svv_cell{f});
end
noise_cov = (tr_S / (F * Ns)) * 0.05 * eye(Ns);

% ============================================================
% 4) GPU init (E-step only; all optimization kept on CPU)
% ============================================================
if USE_GPU
    try
        L_gpu = gpuArray(L);
        noise_cov_gpu = gpuArray(noise_cov);
        Svv_gpu = cell(F, 1);
        for f = 1:F
            Svv_gpu{f} = gpuArray(Svv_cell{f});
        end
        if VERBOSE
            fprintf('[JSPACE-LR1] GPU enabled for E-step tensors.\n');
        end
    catch ME
        warning('[JSPACE-LR1] GPU initialization failed: %s. Falling back to CPU.', ME.message);
        USE_GPU = false;
    end
end

% ============================================================
% 5) Initialization (eLORETA)
% ============================================================
Sigma_source_curr = [];
if isfield(cfg, 'init_sigma_source') && ~isempty(cfg.init_sigma_source)
    Sigma_source_curr = cfg.init_sigma_source;
    if ~iscell(Sigma_source_curr)
        Sigma_source_curr = {Sigma_source_curr};
    end
else
    if VERBOSE
        fprintf('[Init] Running eLORETA...\n');
    end
    if USE_GPU
        S_avg = zeros(Ns, Ns, 'like', Svv_gpu{1});
        for f = 1:F
            S_avg = S_avg + Svv_gpu{f};
        end
        S_avg = S_avg / F;
        [T_eloreta, ~] = run_eloreta_core_(L_gpu, S_avg, 0.05);
    else
        S_avg = zeros(Ns, Ns, 'like', Svv_cell{1});
        for f = 1:F
            S_avg = S_avg + Svv_cell{f};
        end
        S_avg = S_avg / F;
        [T_eloreta, ~] = run_eloreta_core_(L, S_avg, 0.05);
    end

    Sigma_source_curr = cell(F, 1);
    for f = 1:F
        if USE_GPU
            Sjj_eloreta = utils_math.make_hermitian(T_eloreta * Svv_gpu{f} * T_eloreta');
        else
            Sjj_eloreta = utils_math.make_hermitian(T_eloreta * Svv_cell{f} * T_eloreta');
        end
        [Sjj_spd, ~] = utils_math.project_spd(Sjj_eloreta, 1e-8);
        Sigma_source_curr{f} = gather_if_needed_(Sjj_spd);
    end
end

Gamma_warm_start = cell(F, 1);
for f = 1:F
    Gamma_warm_start{f} = eye(Nr);
end

% ============================================================
% 6) Base E-step + whitening
% ============================================================
if VERBOSE
    fprintf('[JSPACE-LR1] Preparing base whitened data...\n');
end
if USE_GPU, L_estep = L_gpu; noise_cov_estep = noise_cov_gpu; else, L_estep = L; noise_cov_estep = noise_cov; end
[Psijj_cell0, ~] = run_estep_cpu_or_gpu_(Svv_cell, Sigma_source_curr, L_estep, noise_cov_estep, USE_GPU);
[Sjj_tilde0, D_cell0, ~] = module1_data_whitening(Psijj_cell0, 'smoothing_window', 3);
Sjj_tilde0 = gather_cell_(Sjj_tilde0);
D_cell0    = gather_cell_(D_cell0);

% ============================================================
% 7) Estimate shared alpha mode u ONCE from base whitened covariance
% ============================================================
freq_axis = get_cfg(cfg, 'freq', []);
if LOWRANK_ENABLE
    [u_shared, u_info] = estimate_shared_mode_once_(Sjj_tilde0, freq_axis, LOWRANK_ALPHA_BAND, VERBOSE);
    if VERBOSE
        fprintf('[JSPACE-LR1] Shared alpha mode estimated once. top_eig=%.3e, alpha_bins=%d\n', ...
            u_info.top_eig, numel(u_info.alpha_idx));
    end

    [gate_vec, gate_info] = build_adaptive_alpha_gate_once_( ...
        Sjj_tilde0, u_shared, freq_axis, K_freq, ...
        LOWRANK_GATE_ENABLE, LOWRANK_GATE_SEARCH_BAND, ...
        LOWRANK_GATE_SMOOTH, LOWRANK_GATE_FRACTION, ...
        LOWRANK_GATE_FLOOR, LOWRANK_GATE_MIN_WIDTH_HZ, ...
        LOWRANK_GATE_MAX_WIDTH_HZ, VERBOSE);

else
    u_shared = [];
    u_info = struct('enabled', false, 'top_eig', NaN, 'alpha_idx', []);
    gate_vec = ones(F,1);
    gate_info = struct('enabled', false, 'peak_freq', NaN, 'peak_idx', NaN);
end

% Base low-rank residuals for eta-search
[alpha0, lowrank0, resid0] = compute_rank1_residuals_( ...
    Sjj_tilde0, u_shared, K_freq, LOWRANK_ALPHA_SMOOTH, ...
    LOWRANK_MIN_EIG, LOWRANK_ENABLE, gate_vec);
% ============================================================
% 8) Build base M-step inputs (on residual covariances)
% ============================================================
m5_in_base = struct();
m5_in_base.whitened_covariances = resid0;
m5_in_base.smoothing_kernel = gather(K_freq);
m5_in_base.weight_matrix    = gather(W_gamma);

m5_p_base = struct();
m5_p_base.alpha0 = 0.1;
m5_p_base.min_eig = MIN_EIG;
m5_p_base.spatial_graph_matrix = gather(GraphLaplacian);
m5_p_base.spatial_graph_is_laplacian = true;
m5_p_base.verbose = false;
m5_p_base.weight_mode = 'hadamard';
m5_p_base.auto_tune = false;
m5_p_base.stoch = stoch_cfg;
m5_p_base.stoch.max_iter = OBJ_STOCH_MAX_ITER;

mstep_solver = @module5_stoch_fista_main;

% ============================================================
% 9) Compute thresholds & eta bounds on residual covariances
% ============================================================
if VERBOSE
    fprintf('[JSPACE-LR1] Computing scales/thresholds on residual covariances...\n');
end
[sc0, thr0] = compute_scales_thresholds_(m5_in_base.whitened_covariances, m5_in_base.smoothing_kernel, m5_in_base.weight_matrix, Nr, cfg);

min_den_limit = log(Nr) / Nr;
max_den_limit = min(thr0.data_density, get_cfg(cfg, 'max_density_cap', 0.40));
max_den_limit = max(max_den_limit, min_den_limit * 1.5);

[mask0, ~] = threshold_active_mask_(m5_in_base.whitened_covariances, thr0.t_active);
if ~isempty(dwi_mask)
    for f = 1:F
        mask0{f} = mask0{f} | dwi_mask;
        mask0{f}(1:Nr+1:end) = true;
    end
end
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
end

% ============================================================
% 10) Surrogate optimization on eta (using residual covariances)
% ============================================================
opt_seed_base = get_cfg(cfg, 'opt_rng_seed', 0);
obj_fun = @(x_log) wrapper_objective_stoch_( ...
    x_log, m5_in_base, m5_p_base, sc0, ...
    M_SAMPLES, min_den_limit, max_den_limit, ...
    PCOR_EPS, DENS_PEN_W, ...
    Nr, F, mstep_solver, ...
    opt_seed_base);

opts = optimoptions('surrogateopt', ...
    'MaxFunctionEvaluations', MAX_OPT_EVALS, ...
    'MinSurrogatePoints', MIN_OPT_POINTS, ...
    'UseParallel', USE_PARALLEL, ...
    'PlotFcn', [], ...
    'InitialPoints', initPts);

if VERBOSE
    fprintf('[JSPACE-LR1] Running surrogateopt on eta...\n');
end

tic_opt = tic;
[x_best_log, fval, exitflag, output_opt, trials] = surrogateopt(obj_fun, eta_bounds.lb_log, eta_bounds.ub_log, opts);
time_opt = toc(tic_opt);
eta_best = 10.^x_best_log;

outs = struct();
outs.eta_bounds = eta_bounds;
outs.opt_results.fval = fval;
outs.opt_results.exitflag = exitflag;
outs.opt_results.output = output_opt;
outs.opt_results.trials = trials;
outs.best_eta.eta1 = eta_best(1);
outs.best_eta.eta2 = eta_best(2);
outs.best_eta.eta3 = eta_best(3);
outs.thresholds = thr0;
outs.density_limits.min = min_den_limit;
outs.density_limits.max = max_den_limit;
outs.lowrank.mode_u = u_shared;
outs.lowrank.mode_info = u_info;
outs.lowrank.alpha0 = alpha0;

if VERBOSE
    fprintf('  > Opt finished in %.2fs. Best score: %.3e\n', time_opt, fval);
    fprintf('  > BEST eta: eta1=%.3e | eta2=%.3e | eta3=%.3e\n', eta_best(1), eta_best(2), eta_best(3));
end

% ============================================================
% 11) Full EM (residual graphical M-step)
% ============================================================
if VERBOSE
    fprintf('\n[JSPACE-LR1] Running Full EM...\n');
end

m5_p_em = m5_p_base;
m5_p_em.stoch.max_iter = EM_STOCH_MAX_ITER;

outs.loglik = zeros(MAX_EM_ITER, 1);
outs.em_density_trace = zeros(MAX_EM_ITER, 1);
outs.em_lambdas = zeros(MAX_EM_ITER, 3);
outs.em_active_mask_density = zeros(MAX_EM_ITER, 1);
outs.em_active_t_used = zeros(MAX_EM_ITER, 1);
outs.lowrank.alpha_trace = zeros(MAX_EM_ITER, F);

active_mask_prev = m5_in_base.active_mask;

for em_iter = 1:MAX_EM_ITER
    iter_tic = tic;

    % ---- E-step ----
    [Psijj_cell, e_stats] = run_estep_cpu_or_gpu_(Svv_cell, Sigma_source_curr, L_estep, noise_cov_estep, USE_GPU);
    outs.loglik(em_iter) = e_stats.log_likelihood;

    % ---- Whitening ----
    [Sjj_tilde, D_cell, ~] = module1_data_whitening(Psijj_cell, 'smoothing_window', 1);
    Sjj_tilde = gather_cell_(Sjj_tilde);
    D_cell    = gather_cell_(D_cell);

    % ---- Rank-1 low-rank subtraction ----
    [alpha_it, lowrank_it, resid_it] = compute_rank1_residuals_( ...
        Sjj_tilde, u_shared, K_freq, LOWRANK_ALPHA_SMOOTH, ...
        LOWRANK_MIN_EIG, LOWRANK_ENABLE, gate_vec);
    outs.lowrank.alpha_trace(em_iter, :) = alpha_it(:).';

    for f = 1:F
        m5_in_base.whitened_covariances{f} = resid_it{f};
    end

    % ---- Dynamic lambdas on residual covariances ----
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

    % ---- Active mask on residual covariances ----
    [mask_new, mask_stats] = threshold_active_mask_with_floor_(m5_in_base.whitened_covariances, thr0.t_active, MASK_DENS_FLOOR);
    outs.em_active_mask_density(em_iter) = mask_stats.density;
    outs.em_active_t_used(em_iter) = mask_stats.t_used;

    if ~isempty(dwi_mask)
        for f = 1:F
            mask_new{f} = mask_new{f} | dwi_mask;
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

    % ---- M-step (residual precision) ----
    local_in = m5_in_base;
    local_in.precision_matrices = cell(F, 1);
    for f = 1:F
        local_in.precision_matrices{f} = Gamma_warm_start{f};
    end
    m5_p_em.stoch.seed = stoch_cfg.seed + em_iter;

    [Gamma_new, res] = mstep_solver(local_in, m5_p_em);
    Gamma_warm_start = Gamma_new;

    % ---- Diagnostics on residual precision ----
    den = density_from_G_pcor_(Gamma_new, PCOR_EPS);
    outs.em_density_trace(em_iter) = den;

    if VERBOSE
        [pcor_max, pcor_p99] = summarize_pcor_cell_(Gamma_new);
        fprintf('    [pcor residual-Γ] max=%.2e | p99=%.2e\n', pcor_max, pcor_p99);
    end

    % ---- Build full covariance candidate ----
    Sigma_full_temp = build_full_covariance_from_rank1_(Gamma_new, lowrank_it, D_cell, LOWRANK_MIN_EIG);

    % ---- Inertia update of FULL covariance state ----
    for f = 1:F
        S_next = Sigma_full_temp{f};
        Sigma_source_curr{f} = (1-UPDATE_RATE) * Sigma_source_curr{f} + UPDATE_RATE * S_next;
        Sigma_source_curr{f} = utils_math.make_hermitian(Sigma_source_curr{f});
        [Sigma_source_curr{f}, ~] = utils_math.project_spd(Sigma_source_curr{f}, LOWRANK_MIN_EIG);
    end

    if VERBOSE
        fprintf('  > EM %d/%d: LogLik=%.3e | dens=%.2f%% | maskDen=%.2f%% | alphaMed=%.2e | lambdas=[%.2e %.2e %.2e] | time=%.2fs\n', ...
            em_iter, MAX_EM_ITER, outs.loglik(em_iter), den*100, ...
            outs.em_active_mask_density(em_iter)*100, median(alpha_it), lam1, lam2, lam3, toc(iter_tic));
        if isfield(res, 'batch_trace')
            covered = false(1, F);
            for kk = 1:numel(res.batch_trace)
                covered(res.batch_trace{kk}) = true;
            end
            fprintf('    [stoch coverage] covered=%.1f%% of freqs in this M-step call\n', 100*mean(covered));
        end
    end
end

% Final point estimates from last EM iteration
Gamma_point = Gamma_new;
D_final     = D_cell;
lowrank_final = lowrank_it;
resid_final = resid_it;
Sigma_full_point = build_full_covariance_from_rank1_(Gamma_point, lowrank_final, D_final, LOWRANK_MIN_EIG);
Omega_point      = recolor_precision_from_whitened_(Gamma_point, D_final, MIN_EIG);

% ============================================================
% 12) Optional debias (direct graph only)
% ============================================================
Gamma_debiased = [];
Omega_debiased = [];
if DO_POST
    if VERBOSE
        fprintf('[JSPACE-LR1] Debiasing direct precision on residual covariances...\n');
    end
    try
        [Gamma_debiased, ~, ~, ~] = module_debias(Gamma_point, resid_final, M_SAMPLES);
        Omega_debiased = recolor_precision_from_whitened_(Gamma_debiased, D_final, MIN_EIG);
    catch ME
        warning('[JSPACE-LR1] Debias failed: %s. Falling back to point estimate.', ME.message);
        Gamma_debiased = [];
        Omega_debiased = [];
    end
end

if DO_POST && DEBIAS_ONLY && ~isempty(Omega_debiased)
    Omega_est = Omega_debiased;
else
    if DO_POST && ~DEBIAS_ONLY && VERBOSE
        warning('[JSPACE-LR1] Only debias-only postprocess is implemented in this fast factor-augmented solver. Returning point-estimate Omega.');
    end
    Omega_est = Omega_point;
end

% Full covariance output is ALWAYS the full covariance point estimate
Sigma_src_est = Sigma_full_point;

% Gather to CPU for saving/visualization safety
Omega_est      = gather_cell_(Omega_est);
Sigma_src_est  = gather_cell_(Sigma_src_est);
Omega_point    = gather_cell_(Omega_point);
Sigma_source_curr = gather_cell_(Sigma_source_curr);
if ~isempty(Omega_debiased), Omega_debiased = gather_cell_(Omega_debiased); end
if ~isempty(Gamma_debiased), Gamma_debiased = gather_cell_(Gamma_debiased); end
Gamma_point = gather_cell_(Gamma_point);
resid_final = gather_cell_(resid_final);
lowrank_final = gather_cell_(lowrank_final);
D_final = gather_cell_(D_final);

% ============================================================
% 13) Outputs
% ============================================================
outs.model_variant = 'factor_augmented_rank1_fast';
outs.output_semantics = [ ...
    'Sigma_src_est is FULL covariance (low-rank common mode + residual covariance). ' ...
    'Omega_est is DIRECT precision from the residual process. They are generally not inverse pairs.' ];

outs.lowrank.enabled = LOWRANK_ENABLE;
outs.lowrank.mode_u = gather_if_needed_(u_shared);
outs.lowrank.mode_info = u_info;
outs.lowrank.gate_vec = gate_vec;
outs.lowrank.gate_info = gate_info;
outs.lowrank.alpha_final = alpha_it;
outs.lowrank.alpha_trace = outs.lowrank.alpha_trace;
outs.lowrank.lowrank_final = lowrank_final;
outs.lowrank.residual_cov_final = resid_final;
outs.lowrank.smooth_xi = LOWRANK_ALPHA_SMOOTH;

outs.sigma_state_em = Sigma_source_curr;      % internal EM state (full covariance state)
outs.Sigma_full_point = Sigma_src_est;        % final full covariance point estimate
outs.Omega_point = Omega_point;               % final direct precision point estimate
outs.Gamma_point = Gamma_point;               % whitened residual precision point estimate
outs.D_final = D_final;

outs.post.debias_only = DEBIAS_ONLY;
outs.post.Omega_debiased = Omega_debiased;
outs.post.Gamma_debiased = Gamma_debiased;
outs.post.best_r = NaN;
outs.post.ray_den = NaN;
outs.post.ray_den_vec = NaN;

outs.global_hyperparams.eta1 = eta_best(1);
outs.global_hyperparams.eta2 = eta_best(2);
outs.global_hyperparams.eta3 = eta_best(3);
outs.global_hyperparams.lambda1_final = outs.em_lambdas(end, 1);
outs.global_hyperparams.lambda2_final = outs.em_lambdas(end, 2);
outs.global_hyperparams.lambda3_final = outs.em_lambdas(end, 3);
end

% ============================================================
% LOCAL HELPERS
% ============================================================
function [Psijj_cell, e_stats] = run_estep_cpu_or_gpu_(Svv_cell, Sigma_source_curr, L, noise_cov, USE_GPU)
    if USE_GPU
        F = numel(Svv_cell);
        Sig_gpu = cell(F,1);
        for f = 1:F
            Sig_gpu{f} = gpuArray(Sigma_source_curr{f});
        end
        [Psijj_cell, e_stats] = module2_estep(cellfun(@gpuArray, Svv_cell, 'UniformOutput', false), L, Sig_gpu, noise_cov);
    else
        [Psijj_cell, e_stats] = module2_estep(Svv_cell, L, Sigma_source_curr, noise_cov);
    end
end

function [u, info] = estimate_shared_mode_once_(Sjj_tilde_cell, freq_axis, alpha_band, VERBOSE)
    F = numel(Sjj_tilde_cell);
    p = size(Sjj_tilde_cell{1},1);
    if isempty(freq_axis)
        alpha_idx = 1:F;
    else
        freq_axis = freq_axis(:).';
        freq_axis = freq_axis(1:min(numel(freq_axis),F));
        alpha_idx = find(freq_axis >= alpha_band(1) & freq_axis <= alpha_band(2));
        if isempty(alpha_idx)
            alpha_idx = 1:F;
        end
    end

    M = zeros(p,p);
    I = eye(p);
    for k = 1:numel(alpha_idx)
        f = alpha_idx(k);
        A = utils_math.make_hermitian(Sjj_tilde_cell{f} - I);
        M = M + A;
    end
    M = utils_math.make_hermitian(M / max(numel(alpha_idx),1));

    [V,D] = eig(M);
    d = real(diag(D));
    [top_eig, idx] = max(d);
    u = V(:,idx);
    if norm(u) < 1e-12 || ~isfinite(norm(u))
        u = ones(p,1) / sqrt(p);
        top_eig = 0;
    else
        u = u / norm(u);
    end

    info = struct();
    info.alpha_idx = alpha_idx;
    info.top_eig = top_eig;
    info.matrix = M;

    if VERBOSE
        fprintf('[JSPACE-LR1] Shared mode estimated on %d alpha-band bins.\n', numel(alpha_idx));
    end
end

function [alpha_vec, lowrank_cell, resid_cell] = compute_rank1_residuals_( ...
    Sjj_tilde_cell, u, K_freq, xi, min_eig, enabled, gate_vec)

    F = numel(Sjj_tilde_cell);
    p = size(Sjj_tilde_cell{1},1);
    I = eye(p);
    lowrank_cell = cell(F,1);
    resid_cell = cell(F,1);

    if nargin < 7 || isempty(gate_vec)
        gate_vec = ones(F,1);
    else
        gate_vec = gate_vec(:);
        if numel(gate_vec) ~= F
            gate_vec = ones(F,1);
        end
    end

    if ~enabled || isempty(u)
        alpha_vec = zeros(F,1);
        for f = 1:F
            lowrank_cell{f} = zeros(p,p);
            Rf = utils_math.make_hermitian(Sjj_tilde_cell{f});
            [Rf, ~] = utils_math.project_spd(Rf, min_eig);
            resid_cell{f} = Rf;
        end
        return;
    end

    b = zeros(F,1);
    uuH = u*u';
    for f = 1:F
        Af = utils_math.make_hermitian(Sjj_tilde_cell{f} - I);
        b(f) = real(u' * Af * u);
    end

    % smooth amplitudes across frequency
    if F > 1 && xi > 0
        Ksym = (K_freq + K_freq') / 2;
        Lk = diag(sum(Ksym, 2)) - Ksym;
        alpha_s = (eye(F) + xi * Lk) \ b;
    else
        alpha_s = b;
    end

    alpha_s = max(0, real(alpha_s));

    % ===== NEW: adaptive gate =====
    alpha_vec = gate_vec .* alpha_s;

    for f = 1:F
        lowrank_cell{f} = alpha_vec(f) * uuH;
        Rf = utils_math.make_hermitian(Sjj_tilde_cell{f} - lowrank_cell{f});
        [Rf, ~] = utils_math.project_spd(Rf, min_eig);
        resid_cell{f} = Rf;
    end
end


function Omega_cell = recolor_precision_from_whitened_(Gamma_cell, D_cell, min_eig)
    F = numel(Gamma_cell);
    Omega_cell = cell(F,1);
    for f = 1:F
        Gf = utils_math.make_hermitian(Gamma_cell{f});
        [Gf, ~] = utils_math.project_spd(Gf, min_eig);
        Df = D_cell{f};
        Of = Df * Gf * Df;
        Of = utils_math.make_hermitian(Of);
        [Of, ~] = utils_math.project_spd(Of, min_eig);
        Omega_cell{f} = Of;
    end
end

function Sigma_cell = build_full_covariance_from_rank1_(Gamma_cell, lowrank_cell, D_cell, min_eig)
    F = numel(Gamma_cell);
    Sigma_cell = cell(F,1);
    for f = 1:F
        Gf = utils_math.make_hermitian(Gamma_cell{f});
        [Gf, ~] = utils_math.project_spd(Gf, min_eig);
        Rinv = Gf \ eye(size(Gf,1));
        Rinv = utils_math.make_hermitian(Rinv);

        Sfull_tilde = utils_math.make_hermitian(lowrank_cell{f} + Rinv);
        [Sfull_tilde, ~] = utils_math.project_spd(Sfull_tilde, min_eig);

        Df = D_cell{f};
        d = real(diag(Df));
        d = max(d, sqrt(min_eig));
        dinv = 1 ./ d;
        Scaling = dinv * dinv';
        Sf = Sfull_tilde .* Scaling;
        Sf = utils_math.make_hermitian(Sf);
        [Sf, ~] = utils_math.project_spd(Sf, min_eig);
        Sigma_cell{f} = Sf;
    end
end

function [pcor_max, pcor_p99] = summarize_pcor_cell_(G_cell)
    F = numel(G_cell);
    pcor_max = 0;
    pcor_p99 = 0;
    for f = 1:F
        G = G_cell{f};
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
end

function C = gather_cell_(C)
    if isempty(C), return; end
    if iscell(C)
        for i = 1:numel(C)
            C{i} = gather_if_needed_(C{i});
        end
    else
        C = gather_if_needed_(C);
    end
end

function X = gather_if_needed_(X)
    try
        if isa(X, 'gpuArray')
            X = gather(X);
        end
    catch
    end
end

function score = wrapper_objective_stoch_(x_log, m5_in, m5_p, sc0, M, min_den, max_den, pcor_eps, penW, Nr, F, solver, seed_base)
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
    local_p.stoch.seed = seed_base; % common random numbers

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
            [ld_f, valid] = utils_math.safe_log_det(Gf);
            if ~valid
                score = 1e15;
                return;
            end
            ld_sum = ld_sum + ld_f;
            tr_sum = tr_sum + real(trace(m5_in.whitened_covariances{ff} * Gf));
        end

        minus_2_ll = M * (tr_sum - ld_sum);
        aic_score = minus_2_ll + 2 * total_edges;
        score = aic_score + penalty;

        if isnan(score) || isinf(score)
            score = 1e15;
        end
    catch
        score = 1e15;
    end
end

function cfg_stoch = get_stoch_cfg_(cfg)
    stoch = struct();
    if isfield(cfg, 'stoch') && ~isempty(cfg.stoch)
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
    if isfield(cfg, 'freq'), stoch.freq = cfg.freq; end
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

function [active_mask, stats] = threshold_active_mask_(Sjj_cell, thresh)
    if ~iscell(Sjj_cell), Sjj_cell = {Sjj_cell}; end
    F = numel(Sjj_cell);
    p = size(Sjj_cell{1},1);
    active_mask = cell(F,1);
    stats.num_active_edges = zeros(F,1);
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
    p = size(Sjj_cell{1},1);
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
    p = size(Sjj_cpu_cell{1},1);

    Sbar = zeros(p,p);
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
        lambda2_mode = lower(get_cfg(cfg, 'lambda2_scale_mode', 'weighted'));
        switch lambda2_mode
            case 'unweighted'
                ratios = S_off;
            otherwise
                W_off = abs(W_gamma_cpu(maskL));
                ratios = S_off ./ max(W_off, eps);
        end
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
        Gbar = zeros(p,p);
        for f = 1:F
            Gbar = Gbar + (Sjj_cpu_cell{f} - I);
        end
        Gbar = Gbar / F;
        gdiff = zeros(F,1);
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

    [t_active, ~] = utils_stats_gmm_threshold_1d(vals);

    if isempty(vals)
        q95 = 0;
        data_density = 0.05;
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

function [W_gamma, dwi_mask] = build_dwi_soft_prior_(C_in, Nr, cfg)
    C = gather_if_needed_(C_in);
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
    for i = 1:(F-1)
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
        edges = nnz(M) / 2;
        dens_vec(f) = edges / (p*(p-1)/2);
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
        if norm(diag(W)-diag(W_old)) / (norm(diag(W_old))+1e-12) < 1e-3
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

    if isfield(cfg,'eta_lb') && numel(cfg.eta_lb)==3
        eta_global_lb = cfg.eta_lb(:).';
    end
    if isfield(cfg,'eta_ub') && numel(cfg.eta_ub)==3
        eta_global_ub = cfg.eta_ub(:).';
    end

    if isfield(cfg,'eta1_lb'), eta_global_lb(1) = cfg.eta1_lb; end
    if isfield(cfg,'eta2_lb'), eta_global_lb(2) = cfg.eta2_lb; end
    if isfield(cfg,'eta3_lb'), eta_global_lb(3) = cfg.eta3_lb; end
    if isfield(cfg,'eta1_ub'), eta_global_ub(1) = cfg.eta1_ub; end
    if isfield(cfg,'eta2_ub'), eta_global_ub(2) = cfg.eta2_ub; end
    if isfield(cfg,'eta3_ub'), eta_global_ub(3) = cfg.eta3_ub; end

    eta0 = [get_cfg_(cfg,'eta1_init',0.1), get_cfg_(cfg,'eta2_init',0.2), get_cfg_(cfg,'eta3_init',0.1)];
    if isfield(cfg,'eta0') && numel(cfg.eta0)==3
        eta0 = cfg.eta0(:).';
    end
    eta0 = max(eta0, eta_global_lb);
    eta0 = min(eta0, eta_global_ub);

    eta_prev = [];
    if isfield(cfg,'eta_prev') && ~isempty(cfg.eta_prev) && numel(cfg.eta_prev)==3
        eta_prev = cfg.eta_prev(:).';
    end
    eta_bank = [];
    if isfield(cfg,'eta_bank') && ~isempty(cfg.eta_bank) && size(cfg.eta_bank,2)==3
        eta_bank = cfg.eta_bank;
    end

    DO_SHRINK = get_cfg_(cfg,'eta_shrink_enable', true);
    SH_DEC    = get_cfg_(cfg,'eta_shrink_decades', 0.7);
    SH_MODE   = lower(get_cfg_(cfg,'eta_shrink_mode','prev'));

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

function [gate_vec, info] = build_adaptive_alpha_gate_once_( ...
    Sjj_tilde_cell, u, freq_axis, K_freq, ...
    gate_enable, search_band, smooth_xi, frac_level, gate_floor, ...
    min_width_hz, max_width_hz, VERBOSE)

    F = numel(Sjj_tilde_cell);
    gate_vec = ones(F,1);

    info = struct();
    info.enabled = false;
    info.peak_idx = NaN;
    info.peak_freq = NaN;
    info.left_width_hz = NaN;
    info.right_width_hz = NaN;
    info.env_raw = [];
    info.env_smooth = [];

    if ~gate_enable || isempty(u)
        return;
    end

    if isempty(freq_axis)
        freq_axis = (1:F).';
    else
        freq_axis = freq_axis(:);
        if numel(freq_axis) < F
            freq_axis = [freq_axis; (numel(freq_axis)+1:F).'];
        else
            freq_axis = freq_axis(1:F);
        end
    end

    p = size(Sjj_tilde_cell{1},1);
    I = eye(p);

    % raw envelope from the shared mode projection
    b = zeros(F,1);
    for f = 1:F
        Af = utils_math.make_hermitian(Sjj_tilde_cell{f} - I);
        b(f) = max(0, real(u' * Af * u));
    end

    % small smoothing before peak detection
    if F > 1 && smooth_xi > 0
        Ksym = (K_freq + K_freq') / 2;
        Lk = diag(sum(Ksym,2)) - Ksym;
        env = (eye(F) + smooth_xi * Lk) \ b;
    else
        env = b;
    end
    env = max(0, real(env));

    info.env_raw = b;
    info.env_smooth = env;

    % broad candidate alpha band
    idx_band = find(freq_axis >= search_band(1) & freq_axis <= search_band(2));
    if isempty(idx_band)
        idx_band = 1:F;
    end

    env_band = env(idx_band);
    [peak_val, loc] = max(env_band);
    if isempty(loc) || peak_val <= 1e-10
        if VERBOSE
            fprintf('[JSPACE-LR1] Adaptive gate not activated: no clear alpha peak in candidate band.\n');
        end
        return;
    end

    peak_idx = idx_band(loc);
    peak_freq = freq_axis(peak_idx);

    % level for width estimation
    frac_level = min(max(frac_level, 0.10), 0.90);
    level = frac_level * peak_val;

    f_left = find_left_crossing_(freq_axis, env, peak_idx, level);
    f_right = find_right_crossing_(freq_axis, env, peak_idx, level);

    if isnan(f_left)
        left_width = 1.5;
    else
        left_width = peak_freq - f_left;
    end
    if isnan(f_right)
        right_width = 1.5;
    else
        right_width = f_right - peak_freq;
    end

    left_width = min(max(left_width, min_width_hz), max_width_hz);
    right_width = min(max(right_width, min_width_hz), max_width_hz);

    % convert width at level=frac_level into sigma
    c = sqrt(2 * log(1 / frac_level));
    sigma_left = left_width / c;
    sigma_right = right_width / c;

    gate = zeros(F,1);
    for f = 1:F
        df = freq_axis(f) - peak_freq;
        if df <= 0
            gate(f) = exp(-0.5 * (df / sigma_left)^2);
        else
            gate(f) = exp(-0.5 * (df / sigma_right)^2);
        end
    end

    % soft floor, not hard truncation
    gate_vec = gate_floor + (1 - gate_floor) * gate;
    gate_vec = min(max(gate_vec, gate_floor), 1);
    gate_vec(peak_idx) = 1;

    info.enabled = true;
    info.peak_idx = peak_idx;
    info.peak_freq = peak_freq;
    info.left_width_hz = left_width;
    info.right_width_hz = right_width;

    if VERBOSE
        fprintf('[JSPACE-LR1] Adaptive alpha gate enabled. peak=%.2f Hz | left=%.2f Hz | right=%.2f Hz\n', ...
            peak_freq, left_width, right_width);
    end
end

function f_left = find_left_crossing_(freq_axis, env, peak_idx, level)
    f_left = NaN;
    if peak_idx <= 1, return; end

    for i = peak_idx:-1:2
        y1 = env(i-1);
        y2 = env(i);
        if (y1 <= level && y2 > level) || (y1 >= level && y2 < level)
            if abs(y2 - y1) < 1e-12
                f_left = freq_axis(i-1);
            else
                t = (level - y1) / (y2 - y1);
                f_left = freq_axis(i-1) + t * (freq_axis(i) - freq_axis(i-1));
            end
            return;
        end
    end
end

function f_right = find_right_crossing_(freq_axis, env, peak_idx, level)
    f_right = NaN;
    F = numel(env);
    if peak_idx >= F, return; end

    for i = peak_idx:(F-1)
        y1 = env(i);
        y2 = env(i+1);
        if (y1 >= level && y2 < level) || (y1 <= level && y2 > level)
            if abs(y2 - y1) < 1e-12
                f_right = freq_axis(i+1);
            else
                t = (level - y1) / (y2 - y1);
                f_right = freq_axis(i) + t * (freq_axis(i+1) - freq_axis(i));
            end
            return;
        end
    end
end