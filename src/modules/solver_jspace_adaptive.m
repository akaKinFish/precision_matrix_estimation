function [Omega_est, Sigma_src_est, outs] = solver_jspace_adaptive(Svv_cell, L, GraphLaplacian, cfg)
% SOLVER_JSPACE_ADAPTIVE - J-SPACE with Adaptive Annealing Strategy
%
% Strategy: "Annealing EM"
%   1. Iter 1-3: Use AIC (Loose metric) to encourage finding structure.
%   2. Iter 4-N: Switch to EBIC (Strict metric) to prune false positives.
%   3. Use Inertia Update: Sigma_new = (1-rate)*Sigma_old + rate*Sigma_est.
%
% This helps survive the "Death Spiral" where early regularization kills weak signals.

    if nargin < 4, cfg = struct(); end

    % --- 1. Setup & Defaults ---
    if ~iscell(Svv_cell), Svv_cell = {Svv_cell}; end
    F = numel(Svv_cell);
    [Ns, Nr] = size(L);
    
    MAX_EM_ITER = get_cfg(cfg, 'max_em_iter', 10);
    GRID_SIZE   = get_cfg(cfg, 'grid_size', 15);
    USE_GPU     = get_cfg(cfg, 'use_gpu', false);
    VERBOSE     = get_cfg(cfg, 'verbose', true);
    
    % Annealing Params
    SWITCH_ITER = 4;   % Switch from AIC to EBIC at iter 4
    UPDATE_RATE = 0.4; % Keep 60% history, add 40% new info
    
    L1_RATIO    = get_cfg(cfg, 'lambda1_ratio', 1.0);
    L3_RATIO    = get_cfg(cfg, 'lambda3_ratio', 0.1);
    ACT_Q       = get_cfg(cfg, 'active_set_q', 0.50); 
    
    if isfield(cfg, 'm_samples'), M_SAMPLES = cfg.m_samples; else, M_SAMPLES = 100 * Nr; end

    % Handle Empty Laplacian
    if isempty(GraphLaplacian)
        L3_RATIO = 0; GraphLaplacian = zeros(Nr, Nr, 'like', L);
    end

    % GPU Setup
    tr_S = 0; for f=1:F, tr_S = tr_S + trace(Svv_cell{f}); end
    noise_cov = (tr_S / (F * Ns)) * 0.05 * eye(Ns);

    if USE_GPU
        try
            L = gpuArray(L); noise_cov = gpuArray(noise_cov);
            for f=1:F, Svv_cell{f} = gpuArray(Svv_cell{f}); end
            if ~isempty(GraphLaplacian), GraphLaplacian = gpuArray(GraphLaplacian); end
            if VERBOSE, fprintf('[J-SPACE-Adaptive] GPU Acceleration ENABLED.\n'); end
        catch
            USE_GPU = false;
        end
    end

    % --- 2. Initialization ---
    if VERBOSE, fprintf('[Init] Running eLORETA...\n'); end
    S_avg = zeros(Ns, Ns, 'like', Svv_cell{1});
    for f=1:F, S_avg = S_avg + Svv_cell{f}; end
    S_avg = S_avg / F;
    [~, W_eloreta] = run_eloreta_core(L, S_avg, 0.05);
    
    diag_power = diag(W_eloreta).^2;
    L_norm = sum(L.^2, 1)';
    scale_factor = trace(S_avg) / (sum(L_norm .* diag_power) + 1e-10);
    
    Sigma_source_curr = cell(F, 1);
    for f=1:F
        Sigma_source_curr{f} = diag(diag_power * scale_factor);
        if USE_GPU, Sigma_source_curr{f} = gpuArray(Sigma_source_curr{f}); end
    end
    
    Gamma_warm_start = cell(F, 1);
    for f=1:F, Gamma_warm_start{f} = eye(Nr, 'like', L); end
    
    outs = struct();
    outs.loglik = [];
    outs.metrics = {};
    % --- 3. EM Loop ---
    for em_iter = 1:MAX_EM_ITER
        iter_tic = tic;
        
        % Determine Metric Strategy
        if em_iter < SWITCH_ITER
            current_metric = 'aic';  % Loose
            metric_gamma = 0;
        else
            current_metric = 'ebic'; % Strict
            metric_gamma = 0.5;
        end
        
        if VERBOSE
            fprintf('\n=== EM Iteration %d/%d [Metric: %s] ===\n', ...
                    em_iter, MAX_EM_ITER, upper(current_metric)); 
        end
        
        % E-Step
        [Psijj_cell, e_stats] = module2_estep(Svv_cell, L, Sigma_source_curr, noise_cov);
        
        % M-Step Prep
        [Sjj_tilde, D_cell, ~] = module1_data_whitening(Psijj_cell, 'smoothing_window', 1);
        
        % Hyperparams
        m6_in.whitened_covariances = Sjj_tilde;
        m6_in.kernel_matrix = eye(F); m6_in.weight_matrix = eye(Nr);
        hp = module6_hyperparameters(m6_in, 'verbose', false);
        
        % M-Step Grid Search (Internal)
        S_ref = Sjj_tilde{1};
        max_val = max(abs(S_ref(tril(true(Nr), -1))));
        min_val = max_val * 0.01;
        lambda_grid = logspace(log10(max_val), log10(min_val), GRID_SIZE);
        
        m5_in.whitened_covariances = Sjj_tilde;
        m5_in.smoothing_kernel = eye(F);
        m5_in.weight_matrix = eye(Nr);
        m5_in.active_mask = module3_active_set(Sjj_tilde, struct('quantile_level', ACT_Q));
        
        m5_p.lambda1 = hp.lambda1 * L1_RATIO;
        m5_p.lambda3 = m5_p.lambda1 * L3_RATIO;
        m5_p.alpha0 = hp.alpha;
        m5_p.spatial_graph_matrix = GraphLaplacian;
        m5_p.spatial_graph_is_laplacian = true;
        m5_p.max_iter = 100; m5_p.tol = 1e-4; m5_p.verbose = false;
        m5_p.weight_mode = 'hadamard'; m5_p.auto_tune = false;
        
        best_score = Inf; best_G = Gamma_warm_start; best_lam = 0; best_den = 0;
        curr_G = Gamma_warm_start;
        
        for k = 1:GRID_SIZE
            lam = lambda_grid(k);
            m5_p.lambda2 = lam;
            m5_in.precision_matrices = curr_G;
            
            [G_temp, ~] = module5_proximal_main(m5_in, m5_p);
            
            % Compute Score (AIC/EBIC)
            G1 = G_temp{1};
            [ld, valid] = utils_math.safe_log_det(G1);
            if ~valid, ld = -1e10; end
            tr_val = real(trace(Sjj_tilde{1} * G1));
            
            G_off = G1; G_off(1:Nr+1:end) = 0;
            n_edges = sum(abs(G_off(:)) > 1e-5)/2;
            den = n_edges/(Nr*(Nr-1)/2);
            if den > 0.25, break; end
            
            minus_2_ll = M_SAMPLES * (tr_val - ld);
            if strcmp(current_metric, 'aic')
                score = minus_2_ll + 2 * n_edges;
            else
                score = minus_2_ll + n_edges * log(M_SAMPLES) + 4 * n_edges * metric_gamma * log(Nr);
            end
            
            if score < best_score
                best_score = score; best_G = G_temp; best_lam = lam; best_den = den;
            end
            curr_G = G_temp;
        end
        
        Gamma_warm_start = best_G;
        
        if VERBOSE
            fprintf('  Selected lambda2=%.2e | Density=%.2f%% | Score=%.2e\n', best_lam, best_den*100, best_score);
        end
        
        % Debiasing
        [Gamma_debiased, Gamma_ray] = module_debias(best_G, Sjj_tilde, M_SAMPLES);
        
        % Recoloring
        m8_in.whitened_precision_matrices = Gamma_ray;
        m8_in.whitening_matrices = D_cell;
        recol = module8_recoloring(m8_in, struct('verbose', false));
        Omega_new = recol.recolored_precision_matrices;
        
        % Inertia Update
        for f=1:F
            Om = (Omega_new{f} + Omega_new{f}')/2;
            try S_next = inv(Om + 1e-12*eye(Nr)); catch, S_next = pinv(Om); end
            S_next = (S_next + S_next')/2;
            
            % Momentum Update
            Sigma_source_curr{f} = (1-UPDATE_RATE)*Sigma_source_curr{f} + UPDATE_RATE*S_next;
        end
        
        outs.loglik(em_iter) = e_stats.log_likelihood;
        outs.metrics{em_iter} = current_metric;
        if VERBOSE, fprintf('  Time: %.2fs\n', toc(iter_tic)); end
    end
    
    if USE_GPU
        for f=1:F, Omega_new{f} = gather(Omega_new{f}); Sigma_source_curr{f} = gather(Sigma_source_curr{f}); end
    end
    Omega_est = Omega_new;
    Sigma_src_est = Sigma_source_curr;
end

% (Helper functions same as above)
function [T, W] = run_eloreta_core(L, Svv, regu)
    [nchan, ndum] = size(L);
    if nargin < 3, regu = 0.05; end
    W = eye(ndum, 'like', L);
    for k = 1:15
        K = (L * W) * L';
        alpha = regu * trace(K) / nchan;
        M = inv(K + alpha * eye(nchan));
        W_old = W;
        for i = 1:ndum
            li = L(:, i);
            val = real(li' * M * li);
            W(i, i) = sqrt(complex(max(val, 1e-12)));
        end
        if norm(diag(W) - diag(W_old)) / norm(diag(W_old)) < 1e-3, break; end
    end
    T = W * L' * inv((L * W * L') + alpha * eye(nchan));
end
function val = get_cfg(s, f, d), if isfield(s, f), val = s.(f); else, val = d; end, end