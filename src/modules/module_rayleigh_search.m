function [Gamma_best, best_r, stats] = module_rayleigh_search(Gamma_debiased, S_whitened, n_samples, K, W, params)
% MODULE_RAYLEIGH_SEARCH Grid search for optimal Rayleigh threshold.
%
% Score(r) = sum_w [logdet(G_w) - tr(S_w G_w)]
%            - lambda1 * sum_{w<w'} k_{w,w'} ||G_w - G_w'||_W^2
%            - lambda3 * sum_w ||G_w||_W^2
%
% L1 is excluded because the Rayleigh threshold itself sparsifies the graph.

    % ----------------------- setup -----------------------
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
    if isfield(params, 'variance_source')
        var_src = params.variance_source;
    else
        var_src = 'debiased'; % options: 'debiased', 'hat'
    end
    if isfield(params, 'Gamma_hat')
        Gamma_hat_cells = params.Gamma_hat;
    else
        Gamma_hat_cells = [];
    end

    Ksym = (K + K') / 2;

    % ----------------------- variance proxy -----------------------
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

    % ----------------------- grid search -----------------------
    for i = 1:length(r_grid)
        r = r_grid(i);
        G_candidate = cell(F, 1);

        for f = 1:F
            G_dense = Gamma_debiased{f};
            V_proxy = Var_proxies{f};
            Threshold = (r / sqrt(complex(n_samples))) * sqrt(complex(V_proxy));

            mask = abs(G_dense) >= Threshold;
            mask(1:size(G_dense,1)+1:end) = true; % keep diagonal

            G_sparse = G_dense;
            G_sparse(~mask) = 0;

            [G_final, ~] = utils_math.project_spd(G_sparse, 1e-8);
            G_candidate{f} = G_final;
        end

        scores(i) = compute_score(G_candidate, S_whitened, Ksym, W, lambda1, lambda3, mode);
    end

    % ----------------------- select best -----------------------
    [max_score, best_idx] = max(scores);
    best_r = r_grid(best_idx);

    % Reconstruct with best r
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
    stats.best_score = max_score;
end

% ============================================================
function score = compute_score(Gamma_cell, S_cell, Ksym, W, lambda1, lambda3, mode)
    % Helper: structural posterior score
    F = numel(Gamma_cell);
    log_lik = 0;

    % Log-likelihood part
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

    % Frequency smoothing penalty
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

    % Spatial smoothing penalty
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
