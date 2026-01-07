% ------------------------------------------------------------
% LOCAL Module: Rayleigh Search (Updated for PCOR + Penalty)
% ------------------------------------------------------------
function [best_r, mask_cell, Var_proxies, stats] = module_rayleigh_search(Gamma_debiased, S_whitened, n_samples, K, W, params)
    F = numel(Gamma_debiased);
    if nargin < 6, params = struct(); end
    if isfield(params, 'r_range') && ~isempty(params.r_range)
        r_grid = params.r_range;
    else
        r_grid = 2.0:0.2:5.0;
    end
    if isfield(params, 'lambda1'), lambda1 = params.lambda1; else, lambda1 = 0; end
    if isfield(params, 'lambda3'), lambda3 = params.lambda3; else, lambda3 = 0; end
    if isfield(params, 'weight_mode'), mode = params.weight_mode; else, mode = 'matrix'; end
    if isfield(params, 'variance_source'), var_src = params.variance_source; else, var_src = 'debiased'; end
    if isfield(params, 'Gamma_hat'), Gamma_hat_cells = params.Gamma_hat; else, Gamma_hat_cells = []; end
    if isfield(params, 'threshold_domain'), thresh_domain = params.threshold_domain; else, thresh_domain = 'entry'; end

    Ksym = (K + K') / 2;
    Var_proxies = cell(F, 1);
    for f = 1:F
        if strcmp(var_src, 'hat') && ~isempty(Gamma_hat_cells)
            G_src = Gamma_hat_cells{f};
        else
            G_src = Gamma_debiased{f};
        end
        d = diag(G_src);
        Var_proxies{f} = real(d * d') + abs(G_src).^2;
    end

    scores = -inf(length(r_grid), 1);
    
    for i = 1:length(r_grid)
        r = r_grid(i);
        G_candidate = cell(F, 1);
        for f = 1:F
            G_dense = Gamma_debiased{f};
            p = size(G_dense,1);
            if strcmpi(thresh_domain, 'pcor')
                d = real(diag(G_dense)); d = max(d, 1e-12);
                denom = sqrt(d * d.');
                P = -G_dense ./ denom; P(1:p+1:end) = 0;
                ThresholdP = (r / sqrt(n_samples));
                mask = abs(P) >= ThresholdP;
            else
                V_proxy = Var_proxies{f};
                Threshold = (r / sqrt(complex(n_samples))) * sqrt(complex(V_proxy));
                mask = abs(G_dense) >= Threshold;
            end
            mask(1:p+1:end) = true;
            G_sparse = G_dense;
            G_sparse(~mask) = 0;
            G_sparse(1:size(G_sparse,1)+1:end) = real(diag(G_sparse));
            [G_final, ~] = utils_math.project_spd(G_sparse, 1e-8);
            G_candidate{f} = utils_math.make_hermitian(G_final);
        end
        
        % Calculate Score
        base_score = compute_score(G_candidate, S_whitened, Ksym, W, lambda1, lambda3, mode);
        
        % Add Density Penalty
        penalty = 0;
        if isfield(params,'density_min') && isfield(params,'density_max') && isfield(params,'density_penalty_weight')
            den_list = zeros(F,1);
            for ff = 1:F
                Gtmp = G_candidate{ff};
                Mtmp = abs(Gtmp) > 0; 
                Mtmp(1:size(Gtmp,1)+1:end) = false;
                den_list(ff) = (nnz(Mtmp)/2) / (size(Gtmp,1)*(size(Gtmp,1)-1)/2);
            end
            den_med = median(den_list);
            if den_med < params.density_min
                penalty = params.density_penalty_weight * ((params.density_min - den_med)/max(params.density_min,eps))^2;
            elseif den_med > params.density_max
                penalty = params.density_penalty_weight * ((den_med - params.density_max)/max(params.density_max,eps))^2;
            end
        end
        scores(i) = base_score - penalty;
    end

    [max_score, best_idx] = max(scores);
    best_r = r_grid(best_idx);
    mask_cell = cell(F, 1);
    for f = 1:F
        G_dense = Gamma_debiased{f};
        p = size(G_dense,1);
        if strcmpi(thresh_domain, 'pcor')
             d = real(diag(G_dense)); d = max(d, 1e-12);
             denom = sqrt(d * d.');
             P = -G_dense ./ denom; P(1:p+1:end) = 0;
             ThresholdP = (best_r / sqrt(n_samples));
             mask = abs(P) >= ThresholdP;
        else
             V_proxy = Var_proxies{f};
             Threshold = (best_r / sqrt(complex(n_samples))) * sqrt(complex(V_proxy));
             mask = abs(G_dense) >= Threshold;
        end
        mask(1:p+1:end) = true;
        mask_cell{f} = mask;
    end
    stats.r_grid = r_grid;
    stats.scores = scores;
    stats.best_score = max_score;
end

function score = compute_score(Gamma_cell, S_cell, Ksym, W, lambda1, lambda3, mode)
    F = numel(Gamma_cell);
    log_lik = 0;
    for f = 1:F
        G = Gamma_cell{f};
        S = S_cell{f};
        [ld, valid] = utils_math.safe_log_det(G);
        if ~valid
            score = -inf;
            return;
        end
        log_lik = log_lik + (ld - real(trace(S * G)));
    end
    freq_pen = 0;
    if lambda1 > 0
        for f1 = 1:F
            for f2 = f1+1:F
                if Ksym(f1, f2) == 0, continue; end
                Gdiff = Gamma_cell{f1} - Gamma_cell{f2};
                if strcmp(mode, 'matrix')
                    term = real(trace(Gdiff' * (W * Gdiff)));
                elseif strcmp(mode, 'hadamard')
                    term = real(sum(sum(conj(Gdiff) .* (W .* Gdiff))));
                else
                    error('module_rayleigh_search:UnknownWeightMode', 'Unknown weight mode');
                end
                freq_pen = freq_pen + Ksym(f1, f2) * term;
            end
        end
    end
    space_pen = 0;
    if lambda3 > 0
        for f = 1:F
            G = Gamma_cell{f};
            if strcmp(mode, 'matrix')
                term = real(trace(G' * (W * G)));
            elseif strcmp(mode, 'hadamard')
                term = real(sum(sum(conj(G) .* (W .* G))));
            else
                error('module_rayleigh_search:UnknownWeightMode', 'Unknown weight mode');
            end
            space_pen = space_pen + term;
        end
    end
    score = log_lik - lambda1 * freq_pen - lambda3 * space_pen;
end