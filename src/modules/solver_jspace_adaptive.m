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
    % Normalize input to cell-of-matrices
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

    % Frequency kernel and spatial weights for smoothing
    K_freq = get_cfg(cfg, 'freq_kernel', eye(F));
    W_gamma = get_cfg(cfg, 'weight_matrix', eye(Nr));

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
    % --- 2. Initialization ---
    if VERBOSE, fprintf('[Init] Running eLORETA...\n'); end
    % 频率平均的传感器协方差，类型跟 Svv_cell{1} 保持一致（CPU/GPU）
    S_avg = zeros(Ns, Ns, 'like', Svv_cell{1});
    for f = 1:F
        S_avg = S_avg + Svv_cell{f};
    end
    S_avg = S_avg / F;

    % eLORETA 反演算子（run_eloreta_core 会自己根据 L 的类型返回 CPU/GPU）
    [T_eloreta, ~] = run_eloreta_core(L, S_avg, 0.05);

    % 用 eLORETA 把每个频段的 Svv 投到源域，得到 full covariance 初始化
    Sigma_source_curr = cell(F, 1);
    for f = 1:F
        Svv_f = Svv_cell{f};                              % 已经可能在 GPU 上
        Sjj_eloreta = T_eloreta * Svv_f * T_eloreta';     % 源域协方差估计
        Sjj_eloreta = utils_math.make_hermitian(Sjj_eloreta); % 数值对称

        % 轻微 SPD 投影（project_spd 本身已经在别处用于 GPU，所以直接用）
        [Sjj_spd, ~] = utils_math.project_spd(Sjj_eloreta, 1e-8);

        Sigma_source_curr{f} = Sjj_spd;                   % 类型随 T_eloreta（CPU/GPU）
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
        if VERBOSE
            diag_stats = cellfun(@(S) [mean(real(diag(S))), min(real(diag(S))), max(real(diag(S)))], Sjj_tilde, 'UniformOutput', false);
            dmat = vertcat(diag_stats{:});
            fprintf('  [Diag] Sjj_tilde diag mean/min/max (freq1): %.3e / %.3e / %.3e\n', dmat(1,1), dmat(1,2), dmat(1,3));
        end
        
        % Hyperparams
        m6_in.whitened_covariances = Sjj_tilde;
        m6_in.kernel_matrix = K_freq; 
        m6_in.weight_matrix = W_gamma;
        hp = module6_hyperparameters(m6_in, 'verbose', false);
        
        % M-Step Grid Search (Internal) - data-driven lambda2 range (deeper)
        S_ref = Sjj_tilde{1};
        mask_off = ~eye(Nr);
        abs_S_off = abs(S_ref(mask_off));
        if isempty(abs_S_off)
            max_corr = 1e-3;
        else
            max_corr = max(abs_S_off);
        end
        if max_corr < 1e-4, max_corr = 1e-3; end
        % Search range based on true off-diagonal strength
        min_val = max_corr * 0.01;
        lambda_grid = logspace(log10(max_corr), log10(min_val), GRID_SIZE);
        
        m5_in.whitened_covariances = Sjj_tilde;
        m5_in.smoothing_kernel = K_freq;
        m5_in.weight_matrix = W_gamma;
        m5_in.active_mask = module3_active_set(Sjj_tilde, struct('quantile_level', ACT_Q));
        
    m5_p.lambda1 = hp.lambda1 * L1_RATIO;
    m5_p.lambda3 = m5_p.lambda1 * L3_RATIO;
    m5_p.alpha0 = hp.alpha;
        m5_p.min_eig = get_cfg(cfg, 'min_eig', 1e-6); % floor for SPD projections
        m5_p.spatial_graph_matrix = GraphLaplacian;
        m5_p.spatial_graph_is_laplacian = true;
        m5_p.max_iter = 100; m5_p.tol = 1e-4; m5_p.verbose = true;
        m5_p.weight_mode = 'hadamard'; m5_p.auto_tune = false;

        if VERBOSE
            fprintf('  [Diag] lambda1=%.2e lambda3=%.2e min_eig=%.1e maxK=%.2f maxW=%.2f\n', ...
                m5_p.lambda1, m5_p.lambda3, m5_p.min_eig, max(K_freq(:)), max(W_gamma(:)));
            fprintf('  [Grid] Range: %.2e to %.2e (logspace %d, Max99.9%%=%.2e)\n', lambda_grid(1), lambda_grid(end), GRID_SIZE, max_corr);
            % Active mask size
            nnz_mask = 0; for f=1:F, nnz_mask = nnz_mask + nnz(m5_in.active_mask{f}); end
            fprintf('  [Diag] active mask nnz=%d (per freq avg %.1f)\n', nnz_mask, nnz_mask/F);
            fprintf('  [Diag] lambda2 first three = %.2e, %.2e, %.2e\n', lambda_grid(1), lambda_grid(min(2,GRID_SIZE)), lambda_grid(min(3,GRID_SIZE)));
        end
        
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

            if VERBOSE && (k == 1 || k == GRID_SIZE || n_edges == 0 || ~isfinite(score))
                fprintf('          [Diag] k=%d lam=%.2e dens=%.2f%% edges=%d score=%.2e\n', ...
                    k, lam, den*100, n_edges, score);
            end
        end
        
        Gamma_warm_start = best_G;
        
        if VERBOSE
            fprintf('  Selected lambda2=%.2e | Density=%.2f%% | Score=%.2e\n', best_lam, best_den*100, best_score);
        end
        
        % Debiasing (dense)
        [Gamma_debiased, ~, Var_proxies, masks0] = module_debias(best_G, Sjj_tilde, M_SAMPLES);

        % Rayleigh threshold search (uses debiased matrices, variance from Gamma_hat)
        ray_params.lambda1 = m5_p.lambda1;
        ray_params.lambda3 = m5_p.lambda3;
        ray_params.weight_mode = m5_p.weight_mode;
        ray_params.variance_source = 'hat';
        ray_params.Gamma_hat = best_G;
        if isfield(cfg, 'rayleigh_range')
            ray_params.r_range = cfg.rayleigh_range;
        else
            ray_params.r_range = 2.0:0.2:5.0;
        end

        if VERBOSE
            fprintf('  [Rayleigh] Searching best threshold (Range: %.1f-%.1f)...\n', ...
                    min(ray_params.r_range), max(ray_params.r_range));
        end

        [best_r, mask_cell, Var_proxies, ~] = module_rayleigh_search( ...
            Gamma_debiased, Sjj_tilde, M_SAMPLES, K_freq, W_gamma, ray_params);

        if VERBOSE
            total_support = 0;
            for f=1:F, total_support = total_support + nnz(mask_cell{f}) - Nr; end
            fprintf('  Rayleigh best r=%.2f | support edges=%d | Quadratic refit on support\n', best_r, total_support/2);
        end

        % Rayleigh support density check and robust refit
        % Use mask_cell to count support edges
        G_mask = mask_cell{1};
        G_mask(1:Nr+1:end) = 0;
        ray_edges = nnz(G_mask)/2;
        ray_den = ray_edges / (Nr*(Nr-1)/2);
        if VERBOSE
            fprintf('  [Rayleigh] r=%.1f | Edges=%d (%.1f%%)\n', best_r, ray_edges, ray_den*100);
        end

        proceed_to_refit = true;
        if ray_den > 0.10
            if VERBOSE
                fprintf('  [Warning] Rayleigh density too high (>10%%). Thinning masks before refit.\n');
            end
            % Thin masks to top 10%% edges by magnitude per frequency
            for f=1:F
                msk = mask_cell{f};
                msk(1:Nr+1:end) = false;
                idx = find(msk);
                if isempty(idx), continue; end
                Gf = Gamma_debiased{f};
                vals = abs(Gf(idx));
                [~, order] = sort(vals, 'descend');
                keep_k = max(1, round(0.10 * numel(order)));
                keep_idx = idx(order(1:keep_k));
                new_mask = false(Nr, Nr);
                new_mask(keep_idx) = true;
                new_mask = new_mask | new_mask'; % symmetrize
                new_mask(1:Nr+1:end) = true;
                mask_cell{f} = new_mask;
            end
            % Recompute density after thinning
            G_mask = mask_cell{1}; G_mask(1:Nr+1:end)=0;
            ray_edges = nnz(G_mask)/2;
            ray_den = ray_edges / (Nr*(Nr-1)/2);
            if VERBOSE
                fprintf('  [Rayleigh] After thinning: Edges=%d (%.1f%%)\n', ray_edges, ray_den*100);
            end
            if ray_den > 0.15
                if VERBOSE
                    fprintf('  [Warning] Still too dense (>15%%). Skipping Refit; keep L1 result.\n');
                end
                proceed_to_refit = false;
            end
        end

        if proceed_to_refit
            refit_mask = mask_cell;

            refit_in = m5_in;
            refit_in.precision_matrices = Gamma_debiased;
            refit_in.active_mask = refit_mask;

            refit_params = m5_p;
            refit_params.lambda2 = 0;
            refit_params.lambda3 = max(m5_p.lambda3, 1e-2); % enforce ridge
            refit_params.max_iter = min(50, m5_p.max_iter);
            refit_params.tol = 1e-3;

            try
                [Gamma_refit, ~] = module5_proximal_main(refit_in, refit_params);
                if rcond(Gamma_refit{1}) < 1e-12
                    error('Refit produced singular matrix');
                end
                best_G = Gamma_refit;
                if VERBOSE, fprintf('  [Refit] Completed successfully.\n'); end
            catch ME
                if VERBOSE
                    fprintf('  [Warning] Refit failed (%s). Reverting to Rayleigh estimate.\n', ME.message);
                end
                for f=1:F
                    best_G{f} = utils_math.project_spd(Gamma_debiased{f}, 1e-5);
                end
            end
        end

        Gamma_warm_start = best_G;

        % Recoloring
        m8_in.whitened_precision_matrices = best_G;
        m8_in.whitening_matrices = D_cell;
        recol = module8_recoloring(m8_in, struct('verbose', false));
        Omega_new = recol.recolored_precision_matrices;

        % Inertia Update with SPD projection before inversion
        for f=1:F
            Om = (Omega_new{f} + Omega_new{f}')/2;
            if exist('regularize_spd','file') == 2
                Om_spd = regularize_spd(Om);
            else
                [Om_spd, ~] = utils_math.project_spd(Om, 1e-8);
            end
            S_next = inv(Om_spd);
            S_next = (S_next + S_next')/2;

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
    % RUN_ELORETA_CORE - eLORETA depth weighting & inverse operator
    % L   : [nchan x ndip] leadfield (CPU or GPU)
    % Svv : [nchan x nchan] sensor covariance (同 L 类型)
    % regu: regularization scalar

    [nchan, ndum] = size(L);
    if nargin < 3, regu = 0.05; end

    % W 与 L 同类型（CPU/GPU）
    W = eye(ndum, 'like', L);

    for k = 1:15
        K = (L * W) * L.';   % 注意用 .' 保持类型一致（对实数来说等价于转置）
        % 单位矩阵也要跟 L 同类型
        I_n = eye(nchan, 'like', L);

        alpha = regu * trace(K) / nchan;

        % 解线性系统代替 inv，兼容 GPU
        M = (K + alpha * I_n) \ I_n;

        W_old = W;
        for i = 1:ndum
            li = L(:, i);
            val = real(li' * M * li);
            W(i, i) = sqrt(complex(max(val, 1e-12)));
        end

        % 简单收敛判据
        diff_diag = diag(W) - diag(W_old);
        if norm(diff_diag) / (norm(diag(W_old)) + 1e-12) < 1e-3
            break;
        end
    end

    % 最终 inverse operator: T = W L' (L W L' + alpha I)^{-1}
    I_n = eye(nchan, 'like', L);
    K_final = (L * W) * L.';
    alpha = regu * trace(K_final) / nchan;
    A = K_final + alpha * I_n;

    % 用线性方程求逆算子：X = A \ I_n，然后 T = W L' X
    X = A \ I_n;
    T = W * (L.' * X);
end

function val = get_cfg(s, f, d), if isfield(s, f), val = s.(f); else, val = d; end, end
