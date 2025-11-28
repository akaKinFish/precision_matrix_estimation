function [Gamma_best, best_r, stats] = module_rayleigh_search(Gamma_debiased, S_whitened, n_samples, K, W, params)
% MODULE_RAYLEIGH_SEARCH - Grid search for optimal Rayleigh threshold
%
% Theory:
%   Select threshold 'r' that maximizes the Structural Posterior Score:
%   Score(r) = LogLikelihood - (lambda1 * Smooth_Freq + lambda3 * Smooth_Space)
%
%   Note: We intentionally EXCLUDE the L1 term (lambda2) here because the 
%   Rayleigh threshold 'r' itself acts as the sparsity control. Including L1
%   would double-count the penalty and lead to over-shrinkage.
%
% Inputs:
%   Gamma_debiased : {F x 1} Dense matrices
%   S_whitened     : {F x 1} Whitened covariance
%   n_samples      : (int) Sample size T
%   K, W           : Smoothing kernels (Freq, Space)
%   params         : Struct with .lambda1, .lambda3, .r_range
%
% Outputs:
%   Gamma_best     : Optimized sparse precision matrices
%   best_r         : Selected threshold

    % 1. Setup
    F = numel(Gamma_debiased);
    if isfield(params, 'r_range') && ~isempty(params.r_range)
        r_grid = params.r_range;
    else
        r_grid = 2.0:0.2:5.0; % Default robust range
    end
    
    if isfield(params, 'lambda1')
        lambda1 = params.lambda1;
    else
        lambda1 = 0;
    end

    if isfield(params, 'lambda3')
        lambda3 = params.lambda3;
    else
        lambda3 = 0;
    end

    % Pre-calculate Variance Proxy (Eq 7.10)
    Var_proxies = cell(F, 1);
    for f = 1:F
        G = Gamma_debiased{f};
        d = diag(G);
        Var_proxies{f} = real(d * d') + abs(G).^2;
    end

    % Pre-calc Laplacian for Smoothing
    % Freq Laplacian
    K_sym = (K + K')/2;
    L_freq = diag(sum(K_sym, 2)) - K_sym;
    
    % Spatial Laplacian (W is usually the spatial graph or weights)
    % If W is just a weight matrix for L1, we might not have a spatial Laplacian here.
    % Assuming W input is relevant for spatial smoothing logic if lambda3 > 0.
    
    scores = -inf(length(r_grid), 1);
    
    % 2. Grid Search
    for i = 1:length(r_grid)
        r = r_grid(i);
        
        % A. Apply Threshold & Project
        G_candidate = cell(F, 1);
        is_valid_run = true;
        
        for f = 1:F
            G_dense = Gamma_debiased{f};
            V_proxy = Var_proxies{f};
            Threshold = (r / sqrt(complex(n_samples))) * sqrt(complex(V_proxy));
            
            mask = abs(G_dense) >= Threshold;
            mask(1:size(G_dense,1)+1:end) = true; % Keep diagonal
            
            G_sparse = G_dense;
            G_sparse(~mask) = 0;
            
            % SPD Projection
            [G_final, ~] = utils_math.project_spd(G_sparse, 1e-8);
            G_candidate{f} = G_final;
        end
        
        % B. Calculate Score (LogLik - Smooth)
        current_score = calculate_posterior_score(G_candidate, S_whitened, L_freq, W, lambda1, lambda3);
        scores(i) = current_score;
    end
    
    % 3. Select Best
    [max_score, best_idx] = max(scores);
    best_r = r_grid(best_idx);
    
    % Re-construct best result
    Gamma_best = cell(F, 1);
    for f = 1:F
        G_dense = Gamma_debiased{f};
        V_proxy = Var_proxies{f};
        Threshold = (best_r / sqrt(complex(n_samples))) * sqrt(complex(V_proxy));
        mask = abs(G_dense) >= Threshold;
        mask(1:size(G_dense,1)+1:end) = true;
        
        G_sparse = G_dense;
        G_sparse(~mask) = 0;
        [Gamma_best{f}, ~] = utils_math.project_spd(G_sparse, 1e-8);
    end
    
    stats.r_grid = r_grid;
    stats.scores = scores;
end

function score = calculate_posterior_score(Gamma_cell, S_cell, L_freq, W_spatial, lam1, lam3)
    % Internal helper to calculate: LLH - Smooth_Penalty
    F = numel(Gamma_cell);
    p = size(Gamma_cell{1}, 1);
    
    log_lik = 0;
    smooth_penalty = 0;
    
    % 1. Log-Likelihood
    for f = 1:F
        G = Gamma_cell{f};
        S = S_cell{f};
        [ld, valid] = utils_math.safe_log_det(G);
        if ~valid
            score = -inf; return;
        end
        % Term: log|G| - tr(S*G)
        log_lik = log_lik + (ld - real(trace(S * G)));
    end
    
    % 2. Frequency Smoothing (Trace Form)
    if lam1 > 0
        for f = 1:F
            neighbor_sum = zeros(p, p);
            for fp = 1:F
                if L_freq(f, fp) ~= 0
                    neighbor_sum = neighbor_sum + L_freq(f, fp) * Gamma_cell{fp};
                end
            end
            % tr(G_f' * G_neighbor) since weights are usually uniform here
            % or tr(G_f' * W * G_neighbor) if W applies to freq smoothing
            % Simplified: Assuming Uniform Spatial Weight for Freq Smoothness
            smooth_penalty = smooth_penalty + real(trace(Gamma_cell{f}' * neighbor_sum));
        end
    end
    
    % 3. Spatial Smoothing
    if lam3 > 0
        % Assuming W_spatial IS the spatial Laplacian or weights
        for f = 1:F
            G = Gamma_cell{f};
            % tr(G * L_spatial * G)
            smooth_penalty = smooth_penalty + (lam3/lam1) * real(trace(G * W_spatial * G));
        end
    end
    
    % Final Score (Maximize this)
    % Penalty is subtracted
    score = log_lik - lam1 * smooth_penalty;
end