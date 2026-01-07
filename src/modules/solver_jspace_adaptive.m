function [Omega_est, Sigma_src_est, outs] = solver_jspace_adaptive(Svv_cell, L, GraphLaplacian, cfg)
% SOLVER_JSPACE_ADAPTIVE - Robust J-SPACE with Correlation Rescue
%
% Debug/diagnostic prints are gated by cfg.debug_print = true.
if nargin < 4, cfg = struct(); end

% ============================================================
% 1. Setup & Input Normalization
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

% --- Configuration ---
MAX_EM_ITER = get_cfg(cfg, 'max_em_iter', 5);
GRID_SIZE   = get_cfg(cfg, 'grid_size', 20);
USE_GPU     = get_cfg(cfg, 'use_gpu', false);
VERBOSE     = get_cfg(cfg, 'verbose', true);
DEBUG       = get_cfg(cfg, 'debug_print', true);
DBG_THR     = get_cfg(cfg, 'debug_thr', 1e-5);

optimizer = lower(get_cfg(cfg, 'optimizer', 'fista'));
if isfield(cfg, 'use_fista') && cfg.use_fista, optimizer = 'fista'; end

if strcmp(optimizer, 'fista')
    mstep_solver = @module5_fista_main;
else
    mstep_solver = @module5_proximal_main;
end

% Annealing Parameters
SWITCH_ITER = 4;
UPDATE_RATE = 0.4;
RESCUE_DENSITY = 0.15; % Target 15% density if model collapses
L1_RATIO    = get_cfg(cfg, 'lambda1_ratio', 1.0);
L3_RATIO    = get_cfg(cfg, 'lambda3_ratio', 0.1);
ACT_Q       = get_cfg(cfg, 'active_set_q', 0.50);

if isfield(cfg, 'm_samples'), M_SAMPLES = cfg.m_samples; else, M_SAMPLES = 100 * Nr; end

HYPER_MODE  = lower(get_cfg(cfg, 'hyperparam_mode', 'external_fixed'));
if strcmpi(HYPER_MODE, 'global_doc7')
    HYPER_MODE = 'external_fixed';
end
IS_EXTERNAL = strcmpi(HYPER_MODE, 'external_fixed');
IS_LEGACY   = strcmpi(HYPER_MODE, 'legacy_local');

K_freq  = get_cfg(cfg, 'freq_kernel', eye(F));
W_gamma = get_cfg(cfg,'weight_matrix', ones(Nr));

if isempty(GraphLaplacian)
    if VERBOSE, fprintf('[J-SPACE] No Laplacian provided. Spatial smoothing disabled.\n'); end
    GraphLaplacian = zeros(Nr, Nr, 'like', L);
end

% --- GPU Initialization ---
tr_S = 0; for f=1:F, tr_S = tr_S + trace(Svv_cell{f}); end
noise_cov = (tr_S / (F * Ns)) * 0.05 * eye(Ns);

if USE_GPU
    try
        try
            parallel.gpu.enableCUDAForwardCompatibility(true);
        catch
        end
        L = gpuArray(L);
        noise_cov = gpuArray(noise_cov);
        for f=1:F, Svv_cell{f} = gpuArray(Svv_cell{f}); end
        if ~isempty(GraphLaplacian), GraphLaplacian = gpuArray(GraphLaplacian); end
        if VERBOSE, fprintf('[J-SPACE] GPU Acceleration ENABLED (Client-side).\n'); end
    catch
        USE_GPU = false;
        warning('[J-SPACE] GPU initialization failed. Falling back to CPU.');
    end
end

% ============================================================
% 2. Initialization (eLORETA or provided warm-start)
% ============================================================
Sigma_source_curr = [];
if isfield(cfg, 'init_sigma_source') && ~isempty(cfg.init_sigma_source)
    Sigma_source_curr = cfg.init_sigma_source;
    if ~iscell(Sigma_source_curr)
        if ndims(Sigma_source_curr) == 3
            tmp = cell(F, 1);
            for f = 1:F, tmp{f} = Sigma_source_curr(:, :, f); end
            Sigma_source_curr = tmp;
        else
            Sigma_source_curr = {Sigma_source_curr};
        end
    end
end

if isempty(Sigma_source_curr)
    if VERBOSE, fprintf('\n[Init] Running eLORETA for initialization...\n'); end
    S_avg = zeros(Ns, Ns, 'like', Svv_cell{1});
    for f = 1:F, S_avg = S_avg + Svv_cell{f}; end
    S_avg = S_avg / F;
    
    [T_eloreta, ~] = run_eloreta_core(L, S_avg, 0.05);
    
    Sigma_source_curr = cell(F, 1);
    for f = 1:F
        Svv_f = Svv_cell{f};
        Sjj_eloreta = T_eloreta * Svv_f * T_eloreta';
        Sjj_eloreta = utils_math.make_hermitian(Sjj_eloreta);
        [Sjj_spd, ~] = utils_math.project_spd(Sjj_eloreta, 1e-8);
        Sigma_source_curr{f} = Sjj_spd;
    end
end

if isfield(cfg, 'init_gamma') && ~isempty(cfg.init_gamma)
    Gamma_warm_start = cfg.init_gamma;
    if ~iscell(Gamma_warm_start)
        if ndims(Gamma_warm_start) == 3
            tmp = cell(F, 1);
            for f = 1:F, tmp{f} = Gamma_warm_start(:, :, f); end
            Gamma_warm_start = tmp;
        else
            Gamma_warm_start = {Gamma_warm_start};
        end
    end
else
    Gamma_warm_start = cell(F, 1);
    for f=1:F, Gamma_warm_start{f} = eye(Nr, 'like', L); end
end

outs = struct();
outs.loglik = [];
outs.metrics = {};
outs.grid_history = struct();
outs.global_hyperparams = struct();
outs.thresholds = struct();
outs.hyperparam_mode = HYPER_MODE;
thresholds_ready = false;

t_active_fixed = get_cfg(cfg, 't_active', []);
t_rescue_fixed = get_cfg(cfg, 't_rescue', []);
lambda1_fixed = get_cfg(cfg, 'lambda1', []);
lambda2_fixed = get_cfg(cfg, 'lambda2', []);
lambda3_fixed = get_cfg(cfg, 'lambda3', []);

% ============================================================
% 3. EM Algorithm Loop
% ============================================================
for em_iter = 1:MAX_EM_ITER
    iter_tic = tic;
    
    % --- Determine Metric Strategy ---
    if em_iter < SWITCH_ITER
        current_metric = 'aic';
        metric_gamma = 0;
    else
        current_metric = 'bic'; % Standard BIC (here using EBIC formula with gamma=0)
        metric_gamma = 0.0;
    end
    
    if VERBOSE
        fprintf('\n=== EM Iteration %d/%d [Metric: %s] ===\n', ...
            em_iter, MAX_EM_ITER, upper(current_metric));
    end
    
    % -------------------------
    % E-Step
    % -------------------------
    [Psijj_cell, e_stats] = module2_estep(Svv_cell, L, Sigma_source_curr, noise_cov);
    
    % -------------------------
    % M-Step Prep: Whitening
    % -------------------------
    [Sjj_tilde, D_cell, ~] = module1_data_whitening(Psijj_cell, 'smoothing_window', 1);
    
    if DEBUG
        S1 = Sjj_tilde{1};
        d  = real(diag(S1)); d = max(d, 1e-12);
        C  = abs(S1) ./ sqrt(d*d');
        C(1:Nr+1:end) = 0;
        vals = C(tril(true(Nr),-1));
        fprintf('  [Dbg] Sjj_tilde{1}: maxCorr=%.4f  medCorr=%.4f  q99=%.4f  offFroRatio=%.3e\n', ...
            max(vals), median(vals), prctile(vals,99), norm(S1-diag(diag(S1)),'fro')/(norm(S1,'fro')+eps));
    end
    
    % [DIAGNOSTIC] Check Data Correlations
    if VERBOSE
        S_check = Sjj_tilde{1};
        S_check(1:Nr+1:end) = 0;
        max_corr_val = max(abs(S_check(:)));
        fprintf('  [Diag] Max Off-Diag Correlation: %.4f (If <0.01, signal is too weak)\n', max_corr_val);
    end

    % --- Prepare CPU Data for M-step ---
    m5_in_base = struct();
    m5_in_base.whitened_covariances = cell(F, 1);
    for f=1:F, m5_in_base.whitened_covariances{f} = gather(Sjj_tilde{f}); end
    m5_in_base.smoothing_kernel = gather(K_freq);
    m5_in_base.weight_matrix = gather(W_gamma);

    % M-step params (base)
    m5_p_base = struct();
    m5_p_base.alpha0 = 0.1;
    m5_p_base.min_eig = get_cfg(cfg, 'min_eig', 1e-6);
    m5_p_base.spatial_graph_matrix = gather(GraphLaplacian);
    m5_p_base.spatial_graph_is_laplacian = true;
    m5_p_base.max_iter = 40;
    m5_p_base.tol = 1e-3;
    m5_p_base.verbose = false;
    m5_p_base.weight_mode = 'hadamard';
    m5_p_base.auto_tune = false;
    
    if strcmp(optimizer, 'fista')
        m5_p_base.backtracking_beta = 0.5;
        m5_p_base.max_backtracking = 25;
        m5_p_base.monotone = true;
        m5_p_base.use_restart = true;
    end
    
    % Warm-start for all lambdas
    common_start_G_cpu = cell(F, 1);
    for f=1:F, common_start_G_cpu{f} = gather(Gamma_warm_start{f}); end

    % ============================================================
    % Hyperparameters & Optimization (External Fixed OR Auto-Grid)
    % ============================================================
    dbg_alpha_final = NaN;
    dbg_tau_prox = NaN;
    dbg_diagLikeCnt = NaN;
    dbg_dens1 = NaN;
    dbg_densMed = NaN;
    dbg_validLD = false;
    dbg_ld = NaN;
    dbg_tr = NaN;

    if IS_EXTERNAL
        % -------------------------
        % MODE 1: EXTERNAL / FIXED
        % -------------------------
        missing_lambdas = isempty(lambda1_fixed) || isempty(lambda2_fixed) || isempty(lambda3_fixed);
        missing_thresh  = isempty(t_active_fixed) || isempty(t_rescue_fixed);
        
        if missing_lambdas || missing_thresh
            if VERBOSE
                fprintf('  [Hyper] Auto-calculating parameters based on Data Gradients (Doc Sec.7)...\n');
            end
            
            [auto_hp, auto_thr] = compute_global_hyperparams_(Sjj_tilde, K_freq, W_gamma, cfg);
            
            if isempty(lambda1_fixed), lambda1_fixed = auto_hp.lambda1; end
            if isempty(lambda2_fixed), lambda2_fixed = auto_hp.lambda2; end
            if isempty(lambda3_fixed), lambda3_fixed = auto_hp.lambda3; end
            
            if isempty(t_active_fixed), t_active_fixed = auto_thr.t_active; end
            if isempty(t_rescue_fixed), t_rescue_fixed = auto_thr.t_rescue; end
            
            if VERBOSE
                 fprintf('          Calculated: L1=%.2e | L2=%.2e | L3=%.2e | t_act=%.2e\n', ...
                    lambda1_fixed, lambda2_fixed, lambda3_fixed, t_active_fixed);
            end
        end
        
        if isempty(lambda1_fixed) || isempty(lambda2_fixed) || isempty(lambda3_fixed)
            error('solver_jspace_adaptive:missing_lambdas', ...
                'Auto-calculation failed and no external lambdas provided.');
        end
        
        % Set Base Params
        m5_p_base.lambda1 = lambda1_fixed;
        m5_p_base.lambda2 = lambda2_fixed;
        m5_p_base.lambda3 = lambda3_fixed;
        
        % Active set
        [m5_in_base.active_mask, act_stats] = threshold_active_mask_( ...
            m5_in_base.whitened_covariances, t_active_fixed);
            
        if ~thresholds_ready
            outs.thresholds.t_active = t_active_fixed;
            outs.thresholds.t_rescue = t_rescue_fixed;
            thresholds_ready = true;
            if VERBOSE
                fprintf('  [GlobalTH] Used: t_active=%.2e t_rescue=%.2e\n', t_active_fixed, t_rescue_fixed);
            end
        end
        if em_iter == 1
            outs.global_hyperparams = struct('lambda1', lambda1_fixed, 'lambda2', lambda2_fixed, 'lambda3', lambda3_fixed);
        end

        % Single Solve
        local_in = m5_in_base;
        local_in.precision_matrices = common_start_G_cpu;
        [best_G_cpu, st5] = mstep_solver(local_in, m5_p_base);
        
        % Stats Calculation (Single)
        ld_sum = 0; tr_sum = 0; n_edges_total = 0;
        dens_all = zeros(F, 1); valid_all = true;
        if DEBUG, maxoff = zeros(F, 1); offratio = zeros(F, 1); end
        
        for ff = 1:F
            Gf = best_G_cpu{ff};
            if any(isnan(Gf(:))) || any(isinf(Gf(:))), ld_f = -Inf; valid_f = false;
            else, [ld_f, valid_f] = utils_math.safe_log_det(Gf); end
            
            if ~valid_f, ld_f = -1e10; valid_all = false; end
            ld_sum = ld_sum + ld_f;
            tr_sum = tr_sum + real(trace(local_in.whitened_covariances{ff} * Gf));
            
            Off = Gf; Off(1:Nr+1:end) = 0;
            edges_f = sum(abs(Off(:)) > 1e-5) / 2;
            n_edges_total = n_edges_total + edges_f;
            dens_all(ff) = edges_f / (Nr * (Nr - 1) / 2);
            
            if DEBUG
                maxoff(ff) = max(abs(Off(:)));
                offratio(ff) = norm(Off, 'fro') / (norm(Gf, 'fro') + eps);
            end
        end
        
        minus_2_ll = M_SAMPLES * (tr_sum - ld_sum);
        best_den = mean(dens_all);
        score_aic  = minus_2_ll + 2 * n_edges_total;
        score_ebic = minus_2_ll + n_edges_total * log(M_SAMPLES) + 4 * n_edges_total * metric_gamma * log(Nr);
        [full_obj, ~] = module_objective.compute(best_G_cpu, local_in.whitened_covariances, K_freq, W_gamma, m5_p_base);
        
        hist_lam = m5_p_base.lambda2; hist_den = best_den; hist_edges = n_edges_total;
        hist_aic = score_aic; hist_ebic = score_ebic; hist_nll = minus_2_ll; hist_full = full_obj;
        grid_Gamma_results = {best_G_cpu};
        
        if strcmp(current_metric, 'aic'), best_score = score_aic; else, best_score = score_ebic; end
        best_idx = 1; best_lam = m5_p_base.lambda2;
        
        if DEBUG
            alpha_final = NaN;
            if isstruct(st5) && isfield(st5,'final_alpha'), alpha_final = st5.final_alpha; end
            dbg_alpha_final = alpha_final; dbg_tau_prox = alpha_final * best_lam;
            dbg_diagLikeCnt = sum(maxoff < DBG_THR); dbg_dens1 = dens_all(1); dbg_densMed = median(dens_all);
            dbg_validLD = valid_all; dbg_ld = ld_sum; dbg_tr = tr_sum;
        end

    else
        % ------------------------------------------------------------
        % MODE 2: DYNAMIC GRID SEARCH (Based on Document Sec.7)
        % ------------------------------------------------------------
        if VERBOSE, fprintf('  [AutoGrid] Calculating critical thresholds from data...\n'); end
        
        % 1. Calculate limits
        [auto_hp, auto_thr] = compute_global_hyperparams_(Sjj_tilde, K_freq, W_gamma, cfg);
        
        % 2. Set Lambda 2 bounds (The Search Parameter)
        upper_bound = auto_hp.lambda2_max;
        if upper_bound <= 1e-9, upper_bound = 0.5; end 
        lower_bound = upper_bound * 1e-3;
        lambda_grid = logspace(log10(upper_bound), log10(lower_bound), GRID_SIZE);
        
        if VERBOSE
            fprintf('  [AutoGrid] Lambda2 Range: %.2e (Max) -> %.2e\n', upper_bound, lower_bound);
        end
        
        % 3. Fix Lambda 1 and 3 based on data characteristics
        m5_p_base.lambda1 = auto_hp.lambda1;
        m5_p_base.lambda3 = auto_hp.lambda3;
        if VERBOSE
            fprintf('  [AutoGrid] Fixed Base: Lambda1=%.2e | Lambda3=%.2e\n', ...
                m5_p_base.lambda1, m5_p_base.lambda3);
        end
        
        % 4. Set Thresholds
        outs.thresholds.t_active = auto_thr.t_active;
        outs.thresholds.t_rescue = auto_thr.t_rescue;
        
        % Set Active Mask for this iteration (using auto threshold)
        [m5_in_base.active_mask, act_stats] = threshold_active_mask_( ...
            m5_in_base.whitened_covariances, auto_thr.t_active);
        max_density = get_cfg(cfg, 'max_density', 0.35); % for selection filter
        
        % -------------------------
        % Grid Search Loop
        % -------------------------
        if VERBOSE, fprintf('  [Grid] Starting PARFOR (CPU Mode)...\n'); end
        
        hist_lam = zeros(GRID_SIZE, 1); hist_den = zeros(GRID_SIZE, 1); hist_edges = zeros(GRID_SIZE, 1);
        hist_aic = zeros(GRID_SIZE, 1); hist_ebic = zeros(GRID_SIZE, 1); hist_nll = zeros(GRID_SIZE, 1);
        hist_full = zeros(GRID_SIZE, 1); grid_Gamma_results = cell(GRID_SIZE, 1);
        
        % Debug arrays
        dbg_alpha_final = NaN(GRID_SIZE,1); dbg_tau_prox = NaN(GRID_SIZE,1); dbg_diagLikeCnt = NaN(GRID_SIZE,1);
        dbg_dens1 = NaN(GRID_SIZE,1); dbg_densMed = NaN(GRID_SIZE,1); dbg_validLD = false(GRID_SIZE,1);
        dbg_ld = NaN(GRID_SIZE,1); dbg_tr = NaN(GRID_SIZE,1);

        parfor k = 1:GRID_SIZE
            lam = lambda_grid(k);
            local_p = m5_p_base;
            local_p.lambda2 = lam;
            
            local_in = m5_in_base;
            local_in.precision_matrices = common_start_G_cpu;
            
            [G_temp, st5] = mstep_solver(local_in, local_p);
            
            % Stats
            ld_sum = 0; tr_sum = 0; n_edges_total = 0; dens_all = zeros(F, 1); valid_all = true;
            if DEBUG, maxoff = zeros(F, 1); end
            
            for ff = 1:F
                Gf = G_temp{ff};
                if any(isnan(Gf(:))) || any(isinf(Gf(:))), ld_f = -Inf; valid_f = false;
                else, [ld_f, valid_f] = utils_math.safe_log_det(Gf); end
                
                if ~valid_f, ld_f = -1e10; valid_all = false; end
                ld_sum = ld_sum + ld_f;
                tr_sum = tr_sum + real(trace(local_in.whitened_covariances{ff} * Gf));
                
                Off = Gf; Off(1:Nr+1:end) = 0;
                edges_f = sum(abs(Off(:)) > 1e-5) / 2;
                n_edges_total = n_edges_total + edges_f;
                dens_all(ff) = edges_f / (Nr * (Nr - 1) / 2);
                
                if DEBUG, maxoff(ff) = max(abs(Off(:))); end
            end
            
            minus_2_ll = M_SAMPLES * (tr_sum - ld_sum);
            den = mean(dens_all);
            score_aic  = minus_2_ll + 2 * n_edges_total;
            score_ebic = minus_2_ll + n_edges_total * log(M_SAMPLES) + 4 * n_edges_total * metric_gamma * log(Nr);
            [full_obj_val, ~] = module_objective.compute(G_temp, local_in.whitened_covariances, K_freq, W_gamma, local_p);
            
            hist_lam(k) = lam; hist_den(k) = den; hist_edges(k) = n_edges_total;
            hist_aic(k) = score_aic; hist_ebic(k) = score_ebic; hist_nll(k) = minus_2_ll;
            hist_full(k) = full_obj_val;
            grid_Gamma_results{k} = G_temp;
            
            if DEBUG
                if isstruct(st5) && isfield(st5,'final_alpha')
                    dbg_alpha_final(k) = st5.final_alpha;
                end
                dbg_tau_prox(k) = dbg_alpha_final(k) * lam;
                dbg_diagLikeCnt(k) = sum(maxoff < DBG_THR);
                dbg_dens1(k) = dens_all(1); dbg_densMed(k) = median(dens_all);
                dbg_validLD(k) = valid_all; dbg_ld(k) = ld_sum; dbg_tr(k) = tr_sum;
            end
        end
        
        % Selection
        if strcmp(current_metric, 'aic'), score_sel = hist_aic; metric_name = 'AIC';
        else, score_sel = hist_ebic; metric_name = 'EBIC'; end
        
        [best_idx, best_score, select_note] = select_lambda_( ...
            score_sel, hist_den, hist_full, hist_nll, max_density, metric_name);
            
        if VERBOSE && ~isempty(select_note), fprintf('  [Select] %s\n', select_note); end
        best_lam   = hist_lam(best_idx);
        best_den   = hist_den(best_idx);
        best_G_cpu = grid_Gamma_results{best_idx};
        
        outs.global_hyperparams = struct('lambda1', m5_p_base.lambda1, ...
                                         'lambda2', best_lam, ...
                                         'lambda3', m5_p_base.lambda3);
    end

    % --- Record History ---
    outs.grid_history(em_iter).lambda = hist_lam;
    outs.grid_history(em_iter).density = hist_den;
    outs.grid_history(em_iter).edges_total = hist_edges;
    outs.grid_history(em_iter).aic = hist_aic;
    outs.grid_history(em_iter).ebic = hist_ebic;
    outs.grid_history(em_iter).nll = hist_nll;
    outs.grid_history(em_iter).full_obj = hist_full;
    outs.grid_history(em_iter).selected_lambda = best_lam;
    
    if VERBOSE
        fprintf('  [Select] Lam=%.2e | Den=%.2f%% | %s=%.2e\n', ...
            best_lam, best_den*100, upper(current_metric), best_score);
    end
    if DEBUG && numel(dbg_validLD) == 1
        % Single run debug print
         fprintf('  [DbgGrid] selIdx=%d | alpha=%.2e | tau=%.2e | validLD=%d | ld=%.3e | tr=%.3e | diagLikeCnt=%d/%d | dens1=%.2f%% densMed=%.2f%%\n', ...
            best_idx, dbg_alpha_final, dbg_tau_prox, dbg_validLD, dbg_ld, dbg_tr, ...
            dbg_diagLikeCnt, F, dbg_dens1*100, dbg_densMed*100);
    elseif DEBUG
        % Grid run debug print
        fprintf('  [DbgGrid] selIdx=%d | alpha=%.2e | tau=%.2e | validLD=%d | ld=%.3e | tr=%.3e | diagLikeCnt=%d/%d | dens1=%.2f%% densMed=%.2f%%\n', ...
            best_idx, dbg_alpha_final(best_idx), dbg_tau_prox(best_idx), dbg_validLD(best_idx), dbg_ld(best_idx), dbg_tr(best_idx), ...
            dbg_diagLikeCnt(best_idx), F, dbg_dens1(best_idx)*100, dbg_densMed(best_idx)*100);
    end

    % ============================================================
    % Post-Processing with "Correlation Rescue"
    % ============================================================
    if USE_GPU
        best_G_gpu = cell(F, 1);
        for f=1:F, best_G_gpu{f} = gpuArray(best_G_cpu{f}); end
        target_G = best_G_gpu;
        target_S = Sjj_tilde;
    else
        target_G = best_G_cpu;
        target_S = m5_in_base.whitened_covariances;
    end
    
    % 1) Debias
    [Gamma_debiased, ~, Var_proxies, ~] = module_debias(target_G, target_S, M_SAMPLES);
    
    % 2) Rayleigh
    ray_params = struct();
    if IS_EXTERNAL
        ray_params.lambda1 = lambda1_fixed;
        ray_params.lambda3 = lambda3_fixed;
    else
        % For grid search mode, reuse the auto-calculated base
        ray_params.lambda1 = m5_p_base.lambda1;
        ray_params.lambda3 = m5_p_base.lambda3;
    end
    
    ray_params.weight_mode = 'hadamard';
    ray_params.variance_source = 'hat';
    ray_params.Gamma_hat = target_G;
    if isfield(cfg, 'rayleigh_range')
        ray_params.r_range = cfg.rayleigh_range;
    else
        ray_params.r_range = 2.0:0.2:5.0;
    end
    
    if USE_GPU
        K_dev = gpuArray(K_freq);
        W_dev = gpuArray(W_gamma);
    else
        K_dev = gather(K_freq);
        W_dev = gather(W_gamma);
    end
    
    [best_r, mask_cell, Var_proxies, ~] = module_rayleigh_search( ...
        Gamma_debiased, target_S, M_SAMPLES, K_dev, W_dev, ray_params);
        
    G_mask = mask_cell{1}; G_mask(1:Nr+1:end) = 0;
    ray_den = (nnz(G_mask)/2) / (Nr*(Nr-1)/2);
    proceed_to_refit = true;
    refit_mask = cell(F, 1);
    
    % [CRITICAL FIX] CORRELATION-BASED RESCUE
    if best_den < 0.01 || ray_den > 0.20
        if IS_EXTERNAL
            threshold_rescue = t_rescue_fixed;
            if VERBOSE
                fprintf('  [Rescue] Model failed (Den=%.1f%%, Ray=%.1f%%). ACTIVATING CORRELATION RESCUE.\n', ...
                    best_den*100, ray_den*100);
                fprintf('           Using fixed t_rescue=%.2e from warm-start.\n', threshold_rescue);
            end
        else
            if VERBOSE
                fprintf('  [Rescue] Model failed (Den=%.1f%%, Ray=%.1f%%). ACTIVATING CORRELATION RESCUE.\n', ...
                    best_den*100, ray_den*100);
                fprintf('           Using t_rescue=%.2e from Auto-Calc.\n', outs.thresholds.t_rescue);
            end
            threshold_rescue = outs.thresholds.t_rescue;
        end
        for f=1:F
            refit_mask{f} = abs(target_S{f}) > threshold_rescue;
            refit_mask{f}(1:Nr+1:end) = true;
        end
    else
        for f=1:F
            refit_mask{f} = mask_cell{f};
            refit_mask{f}(1:Nr+1:end) = true;
        end
    end
    
    % 3) Refit
    if proceed_to_refit
        refit_in = struct();
        refit_in.whitened_covariances = target_S;
        refit_in.smoothing_kernel = K_dev;
        refit_in.weight_matrix = W_dev;
        refit_in.precision_matrices = target_G;
        refit_in.active_mask = refit_mask;
        refit_params = struct();
        refit_params.lambda1 = ray_params.lambda1;
        refit_params.lambda2 = 0;                    % No L1 in refit
        refit_params.lambda3 = max(ray_params.lambda3, 1e-2);
        refit_params.max_iter = 100;
        refit_params.tol = 1e-4;
        refit_params.verbose = false;
        refit_params.spatial_graph_matrix = if_gpu(GraphLaplacian, USE_GPU);
        refit_params.spatial_graph_is_laplacian = true;
        refit_params.weight_mode = 'hadamard';
        try
            [Gamma_new_refit, ~] = mstep_solver(refit_in, refit_params);
            Gamma_refit = Gamma_new_refit;
        catch
            if VERBOSE, fprintf('  [Refit] Failed. Fallback to Target.\n'); end
            Gamma_refit = target_G;
        end
    else
        Gamma_refit = target_G;
    end
    
    % warm start for next EM
    Gamma_warm_start = Gamma_refit;
    
    % 4) Recolor
    m8_in = struct();
    m8_in.whitened_precision_matrices = Gamma_refit;
    m8_in.whitening_matrices = D_cell;
    m8_in.active_set_masks     = refit_mask;
    m8_in.original_covariances = Sigma_source_curr;
    
    recol = module8_recoloring(m8_in, struct('verbose', false));
    Omega_new = recol.recolored_precision_matrices;
    
    if DEBUG && isfield(recol,'transformation_stats') && isstruct(recol.transformation_stats)
        ts = recol.transformation_stats;
        invMed  = safe_nested_stat(ts, 'inv_error_stats', 'median');
        invP95  = safe_nested_stat(ts, 'inv_error_stats', 'percentile_95');
        condGm  = safe_nested_stat(ts, 'condition_stats', 'geometric_mean');
        fprintf('  [DbgRecolor] successRate=%.2f spdRate=%.2f invMed=%.2e invP95=%.2e condGMean=%.2e quality=%.2f\n', ...
            safe_field(ts,'success_rate',NaN), safe_field(ts,'spd_success_rate',NaN), invMed, invP95, condGm, safe_field(ts,'overall_quality_score',NaN));
    end
    
    % 5) Inertia Update
    for f=1:F
        Om = (Omega_new{f} + Omega_new{f}')/2;
        [Om_spd, ~] = utils_math.project_spd(Om, 1e-8);
        S_next = inv(Om_spd);
        S_next = (S_next + S_next')/2;
        Sigma_source_curr{f} = (1-UPDATE_RATE)*Sigma_source_curr{f} + UPDATE_RATE*S_next;
    end
    
    if DEBUG
        maxoff_Om = zeros(F,1); maxoff_S = zeros(F,1);
        for ff=1:F
            Om = Omega_new{ff}; Om=(Om+Om')/2; Off=Om; Off(1:Nr+1:end)=0;
            maxoff_Om(ff)=max(abs(Off(:)));
            Sf = Sigma_source_curr{ff}; Sf=(Sf+Sf')/2; OffS=Sf; OffS(1:Nr+1:end)=0;
            maxoff_S(ff)=max(abs(OffS(:)));
        end
        fprintf('  [DbgUpdate] Omega diagLike=%d/%d | Sigma diagLike=%d/%d  (thr=%.1e)\n', ...
            sum(maxoff_Om<DBG_THR), F, sum(maxoff_S<DBG_THR), F, DBG_THR);
    end
    
    outs.loglik(em_iter) = e_stats.log_likelihood;
    outs.density(em_iter) = best_den;
    outs.metrics{em_iter} = current_metric;
    if isfield(recol, 'transformation_stats')
        outs.recolor_stats{em_iter} = recol.transformation_stats;
    end
    if VERBOSE, fprintf('  Time: %.2fs\n', toc(iter_tic)); end
end

% ============================================================
% 4. Return
% ============================================================
if USE_GPU
    Omega_est = cell(F,1);
    Sigma_src_est = cell(F,1);
    for f=1:F
        Omega_est{f} = gather(Omega_new{f});
        Sigma_src_est{f} = gather(Sigma_source_curr{f});
    end
else
    Omega_est = Omega_new; 
    Sigma_src_est = Sigma_source_curr;
end

if IS_LEGACY && get_cfg(cfg, 'plot', true)
    plot_all_em_iterations(outs);
end
end

% ============================================================
% Helpers
% ============================================================
function out = if_gpu(in, use_gpu)
if use_gpu, out = gpuArray(in); else, out = gather(in); end
end
function val = get_cfg(s, f, d)
if isfield(s, f), val = s.(f); else, val = d; end
end
function v = safe_field(s, f, d)
if isstruct(s) && isfield(s,f)
    v = s.(f);
else
    v = d;
end
end
function v = safe_nested_stat(s, f1, f2)
v = NaN;
if ~isstruct(s) || ~isfield(s,f1), return; end
t = s.(f1);
if ~isstruct(t) || ~isfield(t,f2), return; end
v = t.(f2);
end

function [active_mask, stats] = threshold_active_mask_(Sjj_cell, thresh)
% Build active masks using a fixed absolute threshold on Sjj_tilde.
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

function [hp, thr] = compute_global_hyperparams_(Sjj_tilde, K_freq, W_gamma, cfg)
% Compute global lambdas and fixed thresholds from the first Sjj_tilde.
if ~iscell(Sjj_tilde), Sjj_tilde = {Sjj_tilde}; end
F = numel(Sjj_tilde);
p = size(Sjj_tilde{1}, 1);

S_cpu = cell(F, 1);
for f=1:F, S_cpu{f} = gather(Sjj_tilde{f}); end
K_cpu = gather(K_freq);
W_cpu = gather(W_gamma);

% Mean covariance across frequencies
Sbar = zeros(p, p);
for f=1:F, Sbar = Sbar + S_cpu{f}; end
Sbar = utils_math.make_hermitian(Sbar / F);
I = eye(p);

% Lambda2 upper bound from off-diagonal entries
if p < 2
    lambda2_max = 0;
else
    mask = tril(true(p), -1);
    S_off = abs(Sbar(mask));
    W_off = abs(W_cpu(mask));
    denom = max(W_off, eps);
    ratios = S_off ./ denom;
    if isempty(ratios) || ~any(isfinite(ratios))
        lambda2_max = 0;
    else
        lambda2_max = max(ratios(isfinite(ratios)));
    end
end

% Lambda1 scale from gradient at identity
lambda1_scale = norm(Sbar - I, 'fro') * F;
if ~isfinite(lambda1_scale), lambda1_scale = 0; end

% Lambda3 scale from frequency variability
if F <= 1
    lambda3_scale = 0;
else
    Gbar = zeros(p, p);
    for f=1:F, Gbar = Gbar + (S_cpu{f} - I); end
    Gbar = Gbar / F;
    
    gdiff = zeros(F, 1);
    for f=1:F
        gdiff(f) = norm((S_cpu{f} - I) - Gbar, 'fro');
    end
    lambda3_scale = median(gdiff);
end
if ~isfinite(lambda3_scale), lambda3_scale = 0; end

% Frequency Laplacian norm
if F <= 1
    lk_norm = 0;
else
    Ksym = (K_cpu + K_cpu') / 2;
    Lk = diag(sum(Ksym, 2)) - Ksym;
    lk_norm = norm(Lk, 2);
end

% Global coefficients
c1 = get_cfg(cfg, 'global_c1', 0.1);
c2 = get_cfg(cfg, 'global_c2', 0.2);
c3 = get_cfg(cfg, 'global_c3', 0.1);

lambda1 = c1 * lambda1_scale;
lambda2 = c2 * lambda2_max;
if F <= 1 || lk_norm <= eps
    lambda3 = 0;
else
    lambda3 = c3 * lambda3_scale / (lk_norm + eps);
end

% Optional overrides
if isfield(cfg, 'lambda1_global'), lambda1 = cfg.lambda1_global; end
if isfield(cfg, 'lambda2_global'), lambda2 = cfg.lambda2_global; end
if isfield(cfg, 'lambda3_global'), lambda3 = cfg.lambda3_global; end

hp = struct();
hp.lambda1 = lambda1;
hp.lambda2 = lambda2;
hp.lambda3 = lambda3;
hp.lambda2_max = lambda2_max;
hp.lambda1_scale = lambda1_scale;
hp.lambda3_scale = lambda3_scale;
hp.mode = 'global_doc7';

vals = collect_offdiag_abs_(S_cpu);
[t_active, meta_active] = utils_stats_gmm_threshold_1d(vals);
if isempty(vals)
    q95 = 0;
else
    q95 = quantile(vals, 0.95);
end
t_rescue = max(t_active, q95);
meta_rescue = struct('fallback', false, 'shared_with_active', false, 'quantile_95', q95);

thr = struct();
thr.t_active = t_active;
thr.t_rescue = t_rescue;
thr.meta_active = meta_active;
thr.meta_rescue = meta_rescue;
end

function vals = collect_offdiag_abs_(S_cell)
% Collect absolute off-diagonal values across frequencies.
if ~iscell(S_cell), S_cell = {S_cell}; end
F = numel(S_cell);
p = size(S_cell{1}, 1);
vals = [];
mask = tril(true(p), -1);
for f=1:F
    Sf = S_cell{f};
    v = abs(Sf(mask));
    vals = [vals; v(:)]; %#ok<AGROW>
end
end

function [active_mask, stats, q_used] = tune_active_set_quantile_(InputMatrices, q_init, min_density, max_density)
q_used = q_init;
active_mask = {};
stats = struct('density', NaN);
if q_used <= 0 || q_used >= 1
    [active_mask, stats] = module3_active_set(InputMatrices, struct('quantile_level', q_used));
    return;
end
for iter = 1:3
    [active_mask, stats] = module3_active_set(InputMatrices, struct('quantile_level', q_used));
    den = stats.density;
    if ~isfinite(den)
        break;
    end
    if den > max_density
        q_next = q_used * max_density / max(den, eps);
        q_used = max(q_next, 1e-3);
    elseif den < min_density
        q_next = q_used * min_density / max(den, 1e-6);
        q_used = min(q_next, 0.95);
    else
        break;
    end
end
end

function [best_idx, best_score, note] = select_lambda_(scores, density, full_obj, nll, max_density, metric_name)
note = '';
best_idx = 1;
best_score = Inf;
valid_score = isfinite(scores);
valid = valid_score & density <= max_density;
if any(valid)
    candidates = find(valid);
    [best_score, rel] = min(scores(candidates));
    best_idx = candidates(rel);
    return;
end
if any(valid_score)
    candidates = find(valid_score);
    [best_score, rel] = min(scores(candidates));
    best_idx = candidates(rel);
    note = sprintf('Density filter rejected all candidates; selecting by %s without density constraint.', metric_name);
    return;
end
valid_full = isfinite(full_obj);
if any(valid_full)
    candidates = find(valid_full);
    [best_score, rel] = min(full_obj(candidates));
    best_idx = candidates(rel);
    note = 'Primary metric invalid; selecting by finite full objective.';
    return;
end
valid_nll = isfinite(nll);
if any(valid_nll)
    candidates = find(valid_nll);
    [best_score, rel] = min(nll(candidates));
    best_idx = candidates(rel);
    note = 'Primary metric invalid; selecting by finite minus_2_ll.';
    return;
end
[~, best_idx] = min(density);
best_score = scores(best_idx);
note = 'All metrics invalid; selecting lowest-density candidate.';
end

function plot_all_em_iterations(outs)
if ~isfield(outs, 'grid_history'), return; end
n_iters = length(outs.grid_history);
if n_iters == 0, return; end
figure('Name', 'J-SPACE EM Convergence Analysis', 'Color', 'w', 'Position', [50, 50, 1200, 800]);
cols = 3; rows = ceil(n_iters / cols);
for i = 1:n_iters
    subplot(rows, cols, i);
    hist = outs.grid_history(i);
    lam = hist.lambda; [lam, idx] = sort(lam, 'ascend');
    aic = hist.aic(idx); ebic = hist.ebic(idx); full = hist.full_obj(idx); den = hist.density(idx) * 100;
    
    aic_n = aic - min(aic(isfinite(aic)));
    ebic_n = ebic - min(ebic(isfinite(ebic)));
    full_n = full - min(full(isfinite(full)));
    
    yyaxis left
    p1 = semilogx(lam, aic_n, 'b-', 'LineWidth', 1.5); hold on;
    p2 = semilogx(lam, ebic_n, 'g--', 'LineWidth', 1.5);
    p3 = semilogx(lam, full_n, 'r:', 'LineWidth', 1.2);
    ylabel('Norm Score'); set(gca, 'YColor', 'k');
    
    yyaxis right
    p4 = semilogx(lam, den, 'm-', 'LineWidth', 1.0);
    ylabel('Density (%)'); ylim([0, max(max(den), 35)*1.2]); set(gca, 'YColor', 'm');
    
    xline(hist.selected_lambda, 'k-', 'LineWidth', 1.5);
    xlabel('Lambda 2'); title(sprintf('EM Iter %d', i)); grid on;
    if i == 1, legend([p1, p2, p3, p4], {'AIC', 'EBIC', 'FullObj', 'Density'}, 'Location', 'best'); end
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
    if norm(diag(W)-diag(W_old))/(norm(diag(W_old))+1e-12) < 1e-3, break; end
end
I_n = eye(nchan, 'like', L);
K_final = (L * W) * L.';
alpha = regu * trace(K_final) / nchan;
T = W * (L.' * ((K_final + alpha * I_n) \ I_n));
end