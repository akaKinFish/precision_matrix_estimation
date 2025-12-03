function [Omega_est, Sigma_src_est, outs] = solver_jspace(Svv_cell, L, GraphLaplacian, cfg)
% SOLVER_JSPACE - Joint Spatial-Spectral Precision And Connectivity Estimation
%
% Description:
%   A hierarchical EM framework for estimating source connectivity from EEG/MEG data.
%   It combines source imaging (E-step) with structured sparse precision estimation (M-step).
%
%   Key features:
%   1. Initialization via eLORETA (Depth-compensated power estimation).
%   2. Dynamic Active Set selection (Node + Edge sparsity).
%   3. Gershgorin-based hyperparameter bounds.
%   4. Grid Search + EBIC for optimal sparsity selection.
%
% Inputs:
%   Svv_cell       : {F x 1} Sensor Sample Covariance Matrices
%   L              : (Ns x Nr) Leadfield Matrix
%   GraphLaplacian : (Nr x Nr) Spatial Laplacian matrix (optional, for lambda3)
%   cfg            : Configuration struct
%       .m_samples     : (int) Number of time samples (Crucial for EBIC)
%       .max_em_iter   : (int) Max EM iterations (default 5)
%       .grid_size     : (int) Lambda2 search points (default 10)
%       .lambda1_ratio : (double) Scaling for freq smoothing (default 1.0)
%       .lambda3_ratio : (double) Scaling for spatial smoothing (default 0.1)
%       .active_set_q  : (double) Active set quantile (default 0.15)
%       .use_gpu       : (bool) Enable GPU acceleration
%       .verbose       : (bool)
%
% Outputs:
%   Omega_est      : {F x 1} Estimated Source Precision Matrices
%   Sigma_src_est  : {F x 1} Estimated Source Covariance Matrices
%   outs           : Struct with history and statistics

if nargin < 4, cfg = struct(); end

% ============================================================
% 1. Configuration & Input Standardization
% ============================================================

% Ensure inputs are Cell Arrays
if ~iscell(Svv_cell), Svv_cell = {Svv_cell}; end
F = numel(Svv_cell);
[Ns, Nr] = size(L);

% Defaults
MAX_EM_ITER = get_cfg(cfg, 'max_em_iter', 5);
GRID_SIZE   = get_cfg(cfg, 'grid_size', 20);
USE_GPU     = get_cfg(cfg, 'use_gpu', false);
VERBOSE     = get_cfg(cfg, 'verbose', true);

% Hyperparameter Ratios
L1_RATIO    = get_cfg(cfg, 'lambda1_ratio', 1.0); % Frequency smoothing scaling
L3_RATIO    = get_cfg(cfg, 'lambda3_ratio', 0.1); % Spatial smoothing scaling
ACT_Q       = get_cfg(cfg, 'active_set_q', 0.5); % Keep top 15% initially

% Sample size (Critical for EBIC)
if isfield(cfg, 'm_samples')
    M_SAMPLES = cfg.m_samples;
else
    warning('J-SPACE:SampleCountMissing', 'm_samples not provided. Using heuristic 100*p.');
    M_SAMPLES = 100 * Nr;
end

% Noise Regularization for E-step (5% of sensor trace)
tr_S = 0; for f=1:F, tr_S = tr_S + trace(Svv_cell{f}); end
noise_cov = (tr_S / (F * Ns)) * 0.05 * eye(Ns);

% GPU Transfer
if USE_GPU
    try
        L = gpuArray(L);
        noise_cov = gpuArray(noise_cov);
        for f=1:F, Svv_cell{f} = gpuArray(Svv_cell{f}); end
        if ~isempty(GraphLaplacian), GraphLaplacian = gpuArray(GraphLaplacian); end
        if VERBOSE, fprintf('[J-SPACE] GPU Acceleration ENABLED.\n'); end
    catch
        USE_GPU = false;
        warning('GPU failed. Falling back to CPU.');
    end
end

% ============================================================
% [FIX] Handle Empty GraphLaplacian
% ============================================================
% If user passed [], we must ensure lambda3 is effectively disabled
% or the matrix is valid for multiplication.

if isempty(GraphLaplacian)
    if VERBOSE
        fprintf('[J-SPACE] No Spatial Laplacian provided. Disabling Spatial Smoothing (Lambda3).\n');
    end
    % Option A: Disable parameter
    L3_RATIO = 0;

    % Option B: Create dummy zero matrix (safe for multiplication)
    GraphLaplacian = zeros(Nr, Nr, 'like', L);
end

% History Storage
outs = struct('loglik', [], 'best_lambdas', []);

% ============================================================
% 2. Initialization: eLORETA
% ============================================================
% We use eLORETA to get a depth-compensated estimate of source power (Diagonal Sigma).
% This provides a high-quality "Cold Start" for the EM algorithm.

if VERBOSE, fprintf('\n[Init] Running eLORETA for source power initialization...\n'); end

% Average sensor covariance across frequencies for robust initialization
S_avg = zeros(Ns, Ns, 'like', Svv_cell{1});
for f=1:F, S_avg = S_avg + Svv_cell{f}; end
S_avg = S_avg / F;

% Run eLORETA Core (Local Function)
[~, W_eloreta] = run_eloreta_core(L, S_avg, 0.05);

% Construct Initial Source Covariance
% W_eloreta represents weights ~ 1/std(source).
% Thus Covariance ~ W_eloreta^(-2). Or we can use the MNE solution weighted by W.
% For simplicity and robustness, we take the diagonal of the eLORETA solution.

% Reconstruct source power estimate
% InvOp = W^-1 * L' * inv(L*W^-2*L' + alpha*I) ... simplified:
% Let's use the diagonal weights directly as relative power.
diag_power = diag(W_eloreta).^2;

% Normalize global scale to match data trace
% trace(L * Sigma * L') approx trace(Svv)
% trace(L * diag(p) * L') = sum(diag(L'*L) .* p)
L_norm = sum(L.^2, 1)'; % diag(L'*L)
scale_factor = trace(S_avg) / (sum(L_norm .* diag_power) + 1e-10);

Sigma_source_curr = cell(F, 1);
for f=1:F
    % Initialize as diagonal matrix
    Sigma_source_curr{f} = diag(diag_power * scale_factor);
    if USE_GPU, Sigma_source_curr{f} = gpuArray(Sigma_source_curr{f}); end
end

% Initial Gamma for PGD (Warm Start variable)
Gamma_warm_start = cell(F, 1);
for f=1:F, Gamma_warm_start{f} = eye(Nr, 'like', L); end

% ============================================================
% 3. EM Loop
% ============================================================
for em_iter = 1:MAX_EM_ITER
    iter_tic = tic;
    if VERBOSE, fprintf('\n=== EM Iteration %d/%d ===\n', em_iter, MAX_EM_ITER); end

    % --------------------------------------------------------
    % E-Step: Source Second Moments
    % --------------------------------------------------------
    % Calculate E[s*s' | v] using Module 2
    [Psijj_cell, e_stats] = module2_estep(Svv_cell, L, Sigma_source_curr, noise_cov);
    if VERBOSE
        fprintf('  [E-Step] Log-Likelihood: %.4e\n', e_stats.log_likelihood);
    end

    % --------------------------------------------------------
    % M-Step Part 1: Whitening (Module 1)
    % --------------------------------------------------------
    % Transform Psijj to Whitened Space: Sigma_tilde = D * Psijj * D
    [Sjj_tilde, D_cell, ~] = module1_data_whitening(Psijj_cell, 'smoothing_window', 1);

    % --------------------------------------------------------
    % M-Step Part 2: Active Set Selection (Module 3)
    % --------------------------------------------------------
    % Apply "Safe Screening" based on correlation to reduce problem size.
    % We use the "Intersection" strategy (Node Active Set) as discussed.

    act_params.quantile_level = ACT_Q; % e.g., Keep top 15%
    act_params.strategy = 'intersection';
    act_params.force_diagonal = true;

    [active_mask, ~] = module3_active_set(Sjj_tilde, act_params);

    % --------------------------------------------------------
    % M-Step Part 3: Hyperparameter Theory (Module 6)
    % --------------------------------------------------------
    % Calculate theoretical bounds for Alpha and Lambda1
    m6_input.whitened_covariances = Sjj_tilde;
    m6_input.kernel_matrix = eye(F); % Default freq kernel (update if needed)
    m6_input.weight_matrix = eye(Nr);

    hp_theory = module6_hyperparameters(m6_input, 'verbose', false);

    alpha_used = hp_theory.alpha;
    lambda1_used = hp_theory.lambda1 * L1_RATIO;

    % Heuristic for Lambda3 (Spatial): usually smaller than Lambda1
    lambda3_used = lambda1_used * L3_RATIO;

    % --------------------------------------------------------
    % M-Step Part 4: PGD with Grid Search + EBIC (Module 5)
    % --------------------------------------------------------

    % A. Define Grid (保持不变，基于数据)
    S_ref = Sjj_tilde{1};
    mask_off = tril(true(Nr), -1);
    max_val = max(abs(S_ref(mask_off)));

    % 确保 max_val 有意义
    if max_val < 1e-3, max_val = 1.0; end

    % 搜索范围：从最大相关性开始，下探到 0.1%
    min_val = max_val * 0.001;
    lambda_grid = logspace(log10(max_val), log10(min_val), GRID_SIZE);

    % B. Prepare Fixed Params for PGD
    m5_input.whitened_covariances = Sjj_tilde;
    m5_input.smoothing_kernel     = m6_input.kernel_matrix;
    m5_input.weight_matrix        = eye(Nr);
    m5_input.active_mask          = active_mask;

    m5_params.lambda1 = lambda1_used;
    m5_params.lambda3 = lambda3_used;
    m5_params.spatial_graph_matrix = GraphLaplacian;
    m5_params.spatial_graph_is_laplacian = true;

    % [FIX 1] 激进的步长：不要用保守的 alpha_used，给个大初值，靠回溯去缩减
    m5_params.alpha0   = 0.5;

    % [FIX 2] 更严格的容差，防止早停
    m5_params.tol      = 1e-7;

    m5_params.max_iter = 100;
    m5_params.verbose  = false; % 关掉内部打印，避免刷屏
    m5_params.auto_tune = false; % 关掉内部 Gershgorin，我们手动控制了 alpha
    m5_params.weight_mode = 'hadamard';

    % C. Run Grid Search
    best_score = Inf;
    best_Gamma = Gamma_warm_start;
    best_lambda = lambda_grid(1);
    best_density = 0;

    current_G = Gamma_warm_start;

    if VERBOSE
        fprintf('  [M-Step] Grid Search (MaxVal=%.2e) ...\n', max_val);
        fprintf('          %-10s | %-10s | %-10s | %-12s\n', 'Lambda', 'Density', 'Alpha', 'Score');
    end

    % Get Selection Metric (Default EBIC gamma=0 to encourage edges)
    metric_type = get_cfg(cfg, 'selection_metric', 'ebic');
    ebic_gamma  = get_cfg(cfg, 'ebic_gamma', 0.0);

    for k = 1:GRID_SIZE
        lam = lambda_grid(k);

        % Update Params
        m5_params.lambda2 = lam;
        m5_input.precision_matrices = current_G;

        % [FIX 3] 如果上一次结果是全零（Identity），不要用它做 Warm Start
        % 因为在全零点梯度的变化极小，容易再次陷进去。
        % 重新初始化为 Inv(Sigma) 会更有活力。
        if k > 1
            G_prev = current_G{1}; G_prev(1:Nr+1:end)=0;
            if max(abs(G_prev(:))) < 1e-8
                % Reset Warm Start if previous was dead
                % current_G = Gamma_warm_start; % 或者保持 Identity
                % 更好的策略：增大一点 alpha，刺激它跳出来
                m5_params.alpha0 = 1.0;
            end
        end

        % Run PGD
        [G_temp, stats_temp] = module5_proximal_main(m5_input, m5_params);

        % --- Calculate Score ---
        G1 = G_temp{1};
        [ld, valid] = utils_math.safe_log_det(G1);
        if ~valid, ld = -1e10; end
        tr_val = real(trace(Sjj_tilde{1} * G1));

        G_off = G1; G_off(1:Nr+1:end) = 0;
        num_edges = sum(abs(G_off(:)) > 1e-5) / 2;
        density = num_edges / (Nr*(Nr-1)/2);

        % 熔断机制
        if density > 0.25
            if VERBOSE, fprintf('          %.4e | >25%% (Stop)\n', lam); end
            break;
        end

        % Compute Metric
        minus_2_ll = M_SAMPLES * (tr_val - ld);

        switch lower(metric_type)
            case 'aic',  current_score = minus_2_ll + 2 * num_edges;
            case 'bic',  current_score = minus_2_ll + num_edges * log(M_SAMPLES);
            otherwise,   current_score = minus_2_ll + num_edges * log(M_SAMPLES) + ...
                    4 * num_edges * ebic_gamma * log(Nr);
        end

        % 打印调试信息 (显示最终的 alpha，确认它有没有变大)
        final_alpha = stats_temp.final_alpha;
        is_best = '';
        if current_score < best_score
            best_score = current_score;
            best_Gamma = G_temp;
            best_lambda = lam;
            best_density = density;
            is_best = '(*)';
        end

        if VERBOSE
            fprintf('          %.4e | %5.2f%%     | %.2e    | %.4e %s\n', ...
                lam, density*100, final_alpha, current_score, is_best);
        end

        current_G = G_temp;
    end

    Gamma_warm_start = best_Gamma;

    if VERBOSE
        fprintf('           Selected lambda2=%.2e | Density=%.2f%% | %s=%.2e\n', ...
            best_lambda, best_density*100, upper(metric_type), best_score);
    end
    % [NEW] M-Step Part 4.5: Debiasing & Rayleigh Selection

    enable_debias = isfield(cfg, 'enable_debias') && cfg.enable_debias;
    enable_r_search = isfield(cfg, 'enable_rayleigh_search') && cfg.enable_rayleigh_search;

    if enable_debias
        % 1. Debiasing (Dense)
        Gamma_debiased_cell = cell(F, 1);
        for f=1:F
            G = best_Gamma{f};
            S = Sjj_tilde{f};
            G_tilde = 2*G - G*S*G;
            Gamma_debiased_cell{f} = utils_math.make_hermitian(G_tilde);
        end

        % Prepare Params for Search
        ray_params.lambda1 = lambda1_used;
        ray_params.lambda3 = lambda3_used;

        % Define Range
        if isfield(cfg, 'rayleigh_range')
            ray_params.r_range = cfg.rayleigh_range;
        else
            ray_params.r_range = 2.0:0.2:4.5; % Robust Default
        end

        if enable_r_search
            if VERBOSE, fprintf('  [Rayleigh] Searching best threshold (Range: %.1f-%.1f)...\n', ...
                    min(ray_params.r_range), max(ray_params.r_range)); end

            % Pass GraphLaplacian as 'W' argument for spatial smoothing context
            [Gamma_ray, best_r, ~] = module_rayleigh_search(...
                Gamma_debiased_cell, Sjj_tilde, M_SAMPLES, ...
                m6_input.kernel_matrix, GraphLaplacian, ray_params);

            best_Gamma = Gamma_ray;

            if VERBOSE
                fprintf('           Selected r_th=%.1f\n', best_r);
            end
        else
            % Default Fixed
            [Gamma_ray, ~, ~] = module_rayleigh_search(...
                Gamma_debiased_cell, Sjj_tilde, M_SAMPLES, ...
                m6_input.kernel_matrix, GraphLaplacian, struct('r_range', 3.5));
            best_Gamma = Gamma_ray;
        end
    end

    % --------------------------------------------------------
    % M-Step Part 5: Recoloring (Module 8)
    % --------------------------------------------------------
    m8_input.whitened_precision_matrices = best_Gamma;
    m8_input.whitening_matrices = D_cell;

    recol = module8_recoloring(m8_input, struct('verbose', false));

    % Update Source Covariance for next loop
    Omega_new = recol.recolored_precision_matrices;

    for f=1:F
        % 1. Symmetrize Precision Matrix
        Om = (Omega_new{f} + Omega_new{f}') / 2;

        % 2. Robust Inversion via Eigendecomposition
        %    Direct inv() is unstable for sparse/ill-conditioned matrices.
        %    We use eig() to floor tiny eigenvalues before inversion.

        % Ensure we are working with full matrices for eig()
        if issparse(Om), Om = full(Om); end

        [V, D_vec] = eig(Om, 'vector');

        % Floor eigenvalues: Precision eigenvalues correspond to 1/Variance.
        % Extremely small precision eigenvalues (< 1e-9) lead to exploding variance.
        % We clamp them to a safe minimum (e.g., 1e-8).
        min_prec_tol = 1e-8;
        D_safe = max(real(D_vec), min_prec_tol);

        % Reconstruct Covariance: Sigma = V * D^{-1} * V'
        % Optimized multiplication: V * ( (1./D) .* V' )
        S_next = V * ( (1 ./ D_safe) .* V' );

        % 3. Numerical Sanitization
        % Remove imaginary dust and force symmetry
        S_next = real((S_next + S_next') / 2);

        % 4. Physics Constraint: Diagonal Positivity
        % Variance (Power) implies diagonal elements MUST be positive.
        d_diag = diag(S_next);
        if any(d_diag <= 0)
            % Fix invalid diagonals caused by numerical undershoot
            d_diag(d_diag <= 0) = 1e-12;
            S_next(1:Nr+1:end) = d_diag;

            % Optional: strict SPD projection if needed
            % [S_next, ~] = utils_math.project_spd(S_next, 1e-12);
        end

        % 5. Momentum / Inertia Update
        % Sigma_new = (1 - rate) * Sigma_old + rate * Sigma_estimated
        Sigma_source_curr{f} = (1 - UPDATE_RATE) * Sigma_source_curr{f} + UPDATE_RATE * S_next;

    end

    % Store stats
    outs.loglik(em_iter) = e_stats.log_likelihood;
    outs.best_lambdas(em_iter) = best_lambda;

    if VERBOSE
        fprintf('  Time: %.2fs\n', toc(iter_tic));
    end

end % End EM Loop

% ============================================================
% 4. Finalization
% ============================================================
if USE_GPU
    for f=1:F
        Omega_new{f} = gather(Omega_new{f});
        Sigma_source_curr{f} = gather(Sigma_source_curr{f});
    end
end

Omega_est = Omega_new;
Sigma_src_est = Sigma_source_curr;

end

% ============================================================
% Local Helper: eLORETA Core (Initialization)
% ============================================================
function [T, W] = run_eloreta_core(L, Svv, regu)
% Simplified eLORETA implementation for initialization.
% Code adapted from standard Pascual-Marqui implementation logic.
%
% Returns:
%   W: (Nr x Nr) Diagonal weight matrix (Power estimate)

[nchan, ndum] = size(L);
if nargin < 3, regu = 0.05; end

% Initialize weights (Identity)
W = eye(ndum, 'like', L);

% Iterative update (usually converges in < 15 steps)
for k = 1:15
    % M = inv( L*W^-1*L' + alpha*I )
    % Note: standard eLORETA defines W as weights, here W is ~ Covariance
    % Let's use the standard definition: K = L*W*L'.

    % To match the provided eloreta code logic:
    % W are the spatial weights.

    % Calculate K = L * W * L'
    K = (L * W) * L';

    % Regularization
    alpha = regu * trace(K) / nchan;
    M = inv(K + alpha * eye(nchan));

    % Update Weights
    W_old = W;

    % eLORETA update rule: w_i = sqrt( l_i' * M * l_i )
    % This ensures zero localization error.
    for i = 1:ndum
        li = L(:, i);
        val = real(li' * M * li);
        W(i, i) = sqrt(max(val, 1e-12));
    end

    % Check convergence
    diff = norm(diag(W) - diag(W_old)) / norm(diag(W_old));
    if diff < 1e-3, break; end
end

% Compute T (Transfer Function)
% T = W * L' * M
T = W * L' * inv((L * W * L') + alpha * eye(nchan));
end

function val = get_cfg(s, f, d)
if isfield(s, f), val = s.(f); else, val = d; end
end