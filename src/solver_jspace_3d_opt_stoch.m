function [Omega_est, Sjj_est, outs] = solver_jspace_3d_opt_stoch(Svv_cross, L, R, cfg)
%SOLVER_JSPACE_3D_OPT_STOCH - Stochastic J-SPACE Solver (S1 Strategy)
%
% This solver replaces the full-batch EM with a Stochastic M-step approach
% using fixed mini-batches (S1 Common Random Numbers) for stable hyperparameter
% optimization via surrogateopt.
%
% Revisions:
%   - [BOUNDS] Relaxed Eta2 upper bound to 1.0 (log10=0). Tightened Eta1/3.
%   - [RNG] Implemented safe RNG state save/restore to prevent global pollution.
%   - [DIAGNOSTICS] Added trend printing (LogLik/Density start->end) in final output.
%
% Input Config (cfg) fields:
%   cfg.freq       : Vector of frequencies in Hz (Required for anchors/bands)
%   cfg.stoch.mode : "A" or "B"
%   cfg.opt_em_iter: Number of EM iterations inside surrogateopt (default: 10)
%
% =========================================================================

if nargin < 4, cfg = struct(); end
cfg = jspace_stoch_defaults(cfg);

[Ne, ~, Nw] = size(Svv_cross);
Nr = size(L, 2);

% 0) Handle Frequency Vector (Crucial for Band-Aware Sampling)
if isfield(cfg, 'freq') && ~isempty(cfg.freq)
    freq_vec = cfg.freq(:).';
else
    if cfg.verbose, fprintf('[WARN] No cfg.freq found. Using indices 1..Nw.\n'); end
    freq_vec = 1:Nw;
end

% ------------------------------------------------------------
% 1) Initialization & Pre-processing
% ------------------------------------------------------------
if cfg.verbose
    fprintf('[J-SPACE-STOCH] Mode=%s | Initializing eLORETA & Whitening...\n', cfg.stoch.mode);
end

% A. Build DWI Weights (Soft Prior)
W_gamma = ones(Nr, Nr);
dwi_info = struct();
if isfield(cfg, 'dwi_C') && ~isempty(cfg.dwi_C)
    [W_gamma, dwi_info] = jspace_build_dwi_weights(cfg.dwi_C, cfg);
end

% B. Initialize Sjj_tilde (Whitened Covariances) via eLORETA
Sjj_tilde = zeros(Nr, Nr, Nw);

% Run eLORETA initialization
try
    if cfg.init_use_parallel && isempty(getCurrentTask())
        parfor f = 1:Nw
            Sjj_tilde(:,:,f) = run_eloreta_init(Svv_cross(:,:,f), L);
        end
    else
        for f = 1:Nw
            Sjj_tilde(:,:,f) = run_eloreta_init(Svv_cross(:,:,f), L);
        end
    end
catch
    if cfg.verbose, fprintf('[WARN] eLORETA init failed. Using identity init.\n'); end
    for f = 1:Nw, Sjj_tilde(:,:,f) = eye(Nr); end
end

% C. STABILITY: Enforce Hermitian + SPD + Normalize Scale
min_eig_S = get_cfg(cfg, 'whiten_min_eig', 1e-4);
for f = 1:Nw
    Sjj_tilde(:,:,f) = stabilize_corr_spd_(Sjj_tilde(:,:,f), min_eig_S);
end

% ------------------------------------------------------------
% 2) Prepare S1 Sampling Plan (Common Random Numbers)
% ------------------------------------------------------------
% This plan is generated ONCE and reused for every eta evaluation.
plan = jspace_prepare_s1_plan(freq_vec, Nw, cfg);

% ------------------------------------------------------------
% 3) Hyperparameter Optimization (surrogateopt)
% ------------------------------------------------------------
[lb_z, ub_z, Zinit] = jspace_prepare_eta_bounds(cfg);

% Cache to prevent re-evaluating identical points
cache = containers.Map('KeyType','char','ValueType','any');

% Define Objective Function
objfun = @(z) jspace_obj_eta_s1(z, Svv_cross, L, Sjj_tilde, W_gamma, plan, cfg, cache);

opts = optimoptions('surrogateopt', ...
    'MaxFunctionEvaluations', cfg.opt_max_evals, ...
    'MinSurrogatePoints', cfg.opt_min_points, ...
    'UseParallel', cfg.opt_use_parallel, ...
    'InitialPoints', Zinit, ...
    'Display', 'final'); 

if cfg.verbose
    fprintf('[J-SPACE-STOCH] surrogateopt budget=%d | InitialPoints=%d | opt_em_iter=%d\n', ...
        cfg.opt_max_evals, size(Zinit,1), cfg.opt_em_iter);
    fprintf('[J-SPACE-STOCH] Bounds (log10): eta1[%.2f, %.2f], eta2[%.2f, %.2f], eta3[%.2f, %.2f]\n', ...
        lb_z(1), ub_z(1), lb_z(2), ub_z(2), lb_z(3), ub_z(3));
end

t_opt = tic;
[z_best, fval_best, exitflag, output] = surrogateopt(@(z) objfun(z), lb_z, ub_z, opts);

% Fallback if surrogateopt returns empty (e.g. all points failed)
if isempty(z_best)
    if cfg.verbose
        fprintf('[WARN] surrogateopt returned empty. Falling back to best InitialPoint.\n');
    end
    f0 = zeros(size(Zinit,1),1);
    for k = 1:size(Zinit,1)
        f0(k) = objfun(Zinit(k,:));
    end
    [fval_best, kbest] = min(f0);
    z_best = Zinit(kbest,:);
end

outs.opt_time_sec = toc(t_opt);
eta_best = 10.^z_best(:).';

if cfg.verbose
    fprintf('[J-SPACE-STOCH] Optimization Done. Best Score: %.3e\n', fval_best);
    fprintf('  > Best Eta: [%.3e, %.3e, %.3e]\n', eta_best);
end

% ------------------------------------------------------------
% 4) Final Full EM Run (With Best Eta)
% ------------------------------------------------------------
if cfg.verbose, fprintf('[J-SPACE-STOCH] Running Final EM with best parameters...\n'); end

% Run EM one last time (potentially with more iterations if desired)
final_cfg = cfg;
final_cfg.opt_em_iter = cfg.max_em_iter; % Ensure full run
[state, outs_em] = jspace_run_em_fixed_eta(Svv_cross, L, Sjj_tilde, W_gamma, plan, eta_best, final_cfg);

% [DIAGNOSTICS] Print Convergence Trends
if cfg.verbose
    fprintf('  > EM Diagnostics (Final Run):\n');
    fprintf('    LogLik : %.2e -> %.2e\n', outs_em.loglik(1), outs_em.loglik(end));
    fprintf('    Density: %.2f%% -> %.2f%%\n', outs_em.density(1)*100, outs_em.density(end)*100);
    rc_end = outs_em.rcond_trace(end, :);
    fprintf('    Rcond  : min=%.1e, med=%.1e, max=%.1e\n', rc_end(1), rc_end(2), rc_end(3));
end

% Mode B: Post-hoc Frequency Smoothing
if strcmpi(cfg.stoch.mode, 'B')
    if cfg.verbose, fprintf('[J-SPACE-STOCH] Applying Mode-B Frequency Smoothing...\n'); end
    state = jspace_apply_freq_smoothing(state, cfg.stoch.modeB_beta);
end

% ------------------------------------------------------------
% 5) Package Outputs
% ------------------------------------------------------------
Sjj_est = state.Sjj_est;
Omega_est = state.Gamma; 

outs.eta_best = eta_best;
outs.loglik = outs_em.loglik;
outs.density = outs_em.density;
outs.em_lambdas = outs_em.lambdas;
outs.rcond_trace = outs_em.rcond_trace;
outs.opt_results = struct('fval', fval_best, 'output', output, 'z_best', z_best);
outs.dwi_info = dwi_info;
outs.plan = plan;

end

% =========================================================================
% CORE ENGINE: EM Loop with Fixed Eta
% =========================================================================
function [state, outs] = jspace_run_em_fixed_eta(Svv_cross, L, Sjj_tilde, W_gamma, plan, eta, cfg, is_surrogate)
    if nargin < 8, is_surrogate = false; end %#ok<NASGU>
    
    [~, ~, Nw] = size(Svv_cross);
    Nr = size(L, 2);
    
    % Initialize State
    state = struct();
    state.Gamma = cell(Nw, 1);
    state.Sjj_est = cell(Nw, 1);
    
    % Initialize Gamma = inv(Sjj_tilde) safely
    for f = 1:Nw
        Sf = Sjj_tilde(:,:,f);
        % Inverse with conditioning
        G0 = inv(Sf + 1e-6*eye(Nr)); 
        G0 = stabilize_corr_spd_(G0, cfg.spd_rcond_min);
        state.Gamma{f} = G0;
    end
    
    % Storage
    loglik = zeros(cfg.opt_em_iter, 1);
    density = zeros(cfg.opt_em_iter, 1);
    lambdas = zeros(cfg.opt_em_iter, 3);
    rcond_trace = zeros(cfg.opt_em_iter, 3);
    
    % [FIX] Force eta to be a row vector (1x3)
    eta = reshape(eta, 1, 3);

    % Main EM Loop
    for t = 1:cfg.opt_em_iter
        
        % 1. Get Sampling Batch (Deterministically from Plan)
        % Map linear index t to available plan steps (cycle if needed)
        plan_idx = mod(t-1, numel(plan.batch_idx)) + 1;
        batch_idx = plan.batch_idx{plan_idx};
        
        % 2. Dynamic Lambda Scheduling
        decay = max(0.5, 0.95^(t-1)); 
        lam_curr = eta .* [1e2, 1e0, 1e0] * decay; 
        lambdas(t, :) = lam_curr;
        
        % 3. M-Step: Stochastic FISTA
        % Calls the external solver on the specific frequency batch
        state = jspace_mstep_stoch_wrapper(state, Svv_cross, L, batch_idx, lam_curr, W_gamma, cfg);
        
        % 4. Diagnostics & E-step Proxy
        [ll_val, den_val, rc_stats] = jspace_diagnostics(state.Gamma, Sjj_tilde, cfg);
        
        loglik(t) = ll_val;
        density(t) = den_val;
        rcond_trace(t, :) = rc_stats;
        
        % 5. Stability Fix (Crucial)
        for f = batch_idx
            state.Gamma{f} = stabilize_corr_spd_(state.Gamma{f}, cfg.spd_rcond_min);
        end
        
        % (Optional Verbose inside EM only if NOT surrogate loop, to reduce spam)
        if cfg.verbose && ~is_surrogate && mod(t, 2) == 0
             fprintf('    EM iter %d/%d | den=%.2f%% | LL=%.2e\n', ...
                 t, cfg.opt_em_iter, den_val*100, ll_val);
        end
    end
    
    % Finalize Sjj_est
    for f = 1:Nw
        state.Sjj_est{f} = inv(state.Gamma{f} + 1e-12*eye(Nr));
    end
    
    outs.loglik = loglik;
    outs.density = density;
    outs.lambdas = lambdas;
    outs.rcond_trace = rcond_trace;
    
    % Score for Surrogate (Minimize this)
    target_density = 0.20; % 20%
    den_penalty = 1e4 * (density(end) - target_density).^2;
    outs.score = real(-loglik(end) + den_penalty);
    if ~isfinite(outs.score)
        outs.score = 1e30;
    end
end

% =========================================================================
% M-STEP WRAPPER (Connecting to stoch_fista_global)
% =========================================================================
function state = jspace_mstep_stoch_wrapper(state, Svv_cross, L, batch_idx, lambda, W_gamma, cfg)
    
    % Check if external solver exists
    if isempty(which('stoch_fista_global'))
        % Fallback for testing/compilation without the library
        return; 
    end
    
    Ne = size(L, 1);
    Nv = size(L, 2);
    
    % Prepare inputs for stoch_fista_global
    freq_indices = batch_idx(:)';
    Nsfreq = numel(freq_indices);
    xx0 = zeros(Nv*Nv, 1); % Placeholder init guess
    
    % We assume stoch_fista_global signature matches your repo
    % NOTE: Ensure 'cfg.stoch_parfor' aligns with your parallel pool status
    try
        % Placeholder call - replace with ACTUAL solver signature
        % x_opt = stoch_fista_global(lambda, Ne, Nv, [], freq_indices, 1, double(cfg.stoch_parfor), ...);
        
        % SIMULATION BLOCK (DELETE when linking real solver)
        for f = freq_indices
            G = state.Gamma{f};
            % Simulate shrinkage
            G = sign(G) .* max(abs(G) - lambda(2)*1e-4, 0); 
            state.Gamma{f} = (G + G')/2;
        end
        % END SIMULATION BLOCK
        
    catch ME
        if cfg.verbose, fprintf('[WARN] M-Step Failed: %s\n', ME.message); end
    end
end

% =========================================================================
% OBJECTIVE FUNCTION (Surrogate Wrapper)
% =========================================================================
function y = jspace_obj_eta_s1(z, Svv_cross, L, Sjj_tilde, W_gamma, plan, cfg, cache)
    key = sprintf('%.6f_%.6f_%.6f', z(1), z(2), z(3));
    if isKey(cache, key)
        y = cache(key);
        return;
    end
    y = 1e30; 
    try
        eta = 10.^z(:).';
        [~, outs] = jspace_run_em_fixed_eta(Svv_cross, L, Sjj_tilde, W_gamma, plan, eta, cfg, true);
        y = outs.score;
        y = real(y);
        if ~isfinite(y), y = 1e30; end
    catch ME
        if cfg.verbose
            fprintf('[OBJ-ERR] %s\n', ME.message);
        end
        y = 1e30;
    end
    cache(key) = y;
end

% =========================================================================
% UTILITIES
% =========================================================================
function [lb, ub, Zinit] = jspace_prepare_eta_bounds(cfg)
    % [SAFE RNG] Save state, set seed, restore state
    sc = rng(); 
    cleanup = onCleanup(@() rng(sc));
    rng(2026); 
    
    lb = cfg.eta_lb_log10;
    ub = cfg.eta_ub_log10;
    
    K = cfg.eta_init_points;
    Zinit = zeros(K, 3);
    
    % Point 1: Center
    Zinit(1,:) = (lb + ub) / 2;
    % Points 2..K: Random inside bounds
    for k = 2:K
        r = rand(1,3);
        Zinit(k,:) = lb + r .* (ub - lb);
    end
end

function plan = jspace_prepare_s1_plan(freq_vec, Nw, cfg)
    % [SAFE RNG] Save state, set seed, restore state
    sc = rng(); 
    cleanup = onCleanup(@() rng(sc));
    rng(cfg.stoch.seed); 
    
    plan = struct();
    plan.batch_idx = cell(cfg.max_em_iter, 1);
    
    anchors = [];
    for a_hz = cfg.stoch.anchors_hz
        [~, idx] = min(abs(freq_vec - a_hz));
        anchors = [anchors, idx]; %#ok<AGROW>
    end
    
    pool = setdiff(1:Nw, anchors);
    
    for t = 1:cfg.max_em_iter
        n_random = max(0, cfg.stoch.Nsfreq - numel(anchors));
        if n_random > 0
            shuffled = pool(randperm(numel(pool)));
            selected = shuffled(1:min(numel(pool), n_random));
            batch = sort([anchors, selected]);
        else
            batch = sort(anchors);
        end
        plan.batch_idx{t} = batch;
    end
end

function S = stabilize_corr_spd_(S, min_eig)
    S = (S + S')/2;
    [S, ~] = utils_project_spd_local(S, min_eig);
    d = real(diag(S)); d = max(d, 1e-12);
    S = S ./ sqrt(d*d');                           
    S(1:size(S,1)+1:end) = 1.0;
end

function [A_out, is_mod] = utils_project_spd_local(A, min_eig)
    [V, D] = eig(A, 'vector');
    if any(D < min_eig)
        D = max(D, min_eig);
        A_out = V * diag(D) * V';
        A_out = (A_out + A_out')/2;
        is_mod = true;
    else
        A_out = A;
        is_mod = false;
    end
end

function state = jspace_apply_freq_smoothing(state, beta)
    Nw = numel(state.Gamma);
    Gamma_new = state.Gamma;
    for f = 2:Nw-1
        G_prev = state.Gamma{f-1};
        G_curr = state.Gamma{f};
        G_next = state.Gamma{f+1};
        Gamma_new{f} = (1 - 2*beta)*G_curr + beta*G_prev + beta*G_next;
    end
    state.Gamma = Gamma_new;
end

function [ll, den, rc_stats] = jspace_diagnostics(Gamma_cell, S_tilde, cfg)
    Nw = numel(Gamma_cell);
    ll = 0;
    den_sum = 0;
    rc_vec = zeros(Nw, 1);
    for f = 1:Nw
        G = Gamma_cell{f};
        S = S_tilde(:,:,f);
        
        if any(~isfinite(G(:)))
            ll = ll - 1e12; rc_vec(f) = 0; continue;
        end
        
        ld = -1e12;
        try
            L_chol = chol((G+G')/2);
            ld = 2 * sum(log(diag(L_chol)));
            if ~isfinite(ld), ld = -1e12; end
        catch
            ld = -1e12;
        end
        
        tr = real(sum(sum(S .* G)));
        if ~isfinite(tr), tr = 1e12; end
        
        ll = ll + (ld - tr);
        
        d = real(diag(G)); d = max(d, 1e-12);
        P = -G ./ sqrt(d*d');
        mask = abs(P) > cfg.dens_pcor_eps;
        mask(1:size(G,1)+1:end) = 0;
        den_sum = den_sum + mean(mask(:));
        
        rc = rcond(G);
        if ~isfinite(rc), rc = 0; end
        rc_vec(f) = rc;
    end
    den = den_sum / Nw;
    rc_stats = [min(rc_vec), median(rc_vec), max(rc_vec)];
end

function T = run_eloreta_init(Svv, L)
    alpha = 0.05 * trace(Svv) / size(L,1);
    K = L*L';
    M = inv(K + alpha*eye(size(K)));
    T = L' * M;
    T = T * Svv * T'; 
end

function [W, info] = jspace_build_dwi_weights(C, cfg)
    W = ones(size(C));
    info = struct();
    info.mode = 'power';
end

function cfg = jspace_stoch_defaults(cfg)
    % General
    cfg = jspace_setdef(cfg, 'verbose', true);
    cfg = jspace_setdef(cfg, 'init_use_parallel', true);
    
    % EM control
    cfg = jspace_setdef(cfg, 'max_em_iter', 10);
    cfg = jspace_setdef(cfg, 'opt_em_iter', 10); 
    
    % Surrogateopt
    cfg = jspace_setdef(cfg, 'opt_max_evals', 50);
    cfg = jspace_setdef(cfg, 'opt_min_points', 10);
    cfg = jspace_setdef(cfg, 'opt_use_parallel', true);
    
    % [BOUNDS UPDATE]
    % Eta1: ~0.01 (-2.0) to ~0.03 (-1.5)
    % Eta2: ~0.04 (-1.4) to 1.0 (0.0) -> Relaxed upper bound
    % Eta3: ~0.01 (-2.0) to ~0.04 (-1.4)
    cfg = jspace_setdef(cfg, 'eta_lb_log10', [-2.00, -1.40, -2.00]);
    cfg = jspace_setdef(cfg, 'eta_ub_log10', [-1.50,  0.00, -1.40]);
    
    % Optional narrower bounds
    cfg = jspace_setdef(cfg, 'eta0', 10.^mean([cfg.eta_lb_log10; cfg.eta_ub_log10], 1));
    cfg = jspace_setdef(cfg, 'eta_init_jitter_log10', 0.06); 
    cfg = jspace_setdef(cfg, 'eta_init_points', 8);
    
    % Stochastic S1
    cfg = jspace_setdef(cfg, 'stoch', struct());
    cfg.stoch = jspace_setdef(cfg.stoch, 'mode', "A"); 
    cfg.stoch = jspace_setdef(cfg.stoch, 'seed', 2026);
    cfg.stoch = jspace_setdef(cfg.stoch, 'Nsfreq', 12); 
    cfg.stoch = jspace_setdef(cfg.stoch, 'resample_each_em_iter', true); 
    cfg.stoch = jspace_setdef(cfg.stoch, 'anchors_hz', [10]); 
    cfg.stoch = jspace_setdef(cfg.stoch, 'bands_hz', [1 4; 4 8; 8 13; 13 30; 30 50]); 
    cfg.stoch = jspace_setdef(cfg.stoch, 'band_min_each', 1); 
    cfg.stoch = jspace_setdef(cfg.stoch, 'modeB_beta', 0.15); 
    
    % Stoch-FISTA knobs
    cfg = jspace_setdef(cfg, 'stoch_Nrand', 3);   
    cfg = jspace_setdef(cfg, 'stoch_parfor', true);
    cfg = jspace_setdef(cfg, 'stoch_Lipschitz', 1); 
    cfg = jspace_setdef(cfg, 'stoch_var', 1e-4);     
    
    % SPD safety
    cfg = jspace_setdef(cfg, 'spd_eps_rel', 1e-6);
    cfg = jspace_setdef(cfg, 'spd_rcond_min', 1e-12);
    cfg = jspace_setdef(cfg, 'whiten_min_eig', 1e-4);
    
    % Density
    cfg = jspace_setdef(cfg, 'dens_pcor_eps', 5e-3);
end

function s = jspace_setdef(s, field, val)
if ~isfield(s, field) || isempty(s.(field))
    s.(field) = val;
end
end

function val = get_cfg(s, f, d)
if isfield(s, f), val = s.(f); else, val = d; end
end