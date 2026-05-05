function [best_r, mask_cell, Var_proxies, stats] = module_rayleigh_search_v2(Gamma_stat, S_whitened, n_samples, K, W, params)
% MODULE_RAYLEIGH_SEARCH
%   Search the Rayleigh threshold r using:
%       1) the debiased statistic Gamma_stat to define masks,
%       2) a separate SPD template (typically a HIGGS-like ridge template)
%          to evaluate the masked candidates.
%
% This is the key change that makes the post-EM search closer to HIGGS:
%   - debiased matrix -> support statistic
%   - ridge/Frobenius template -> masked SPD estimator to score
%
% Inputs:
%   Gamma_stat : {F x 1} debiased statistics (typically Gamma_db)
%   S_whitened : {F x 1} final whitened covariances
%   n_samples  : effective sample-size proxy used in the Rayleigh scaling
%   K, W       : JSPACE frequency kernel and weight matrix
%   params     : optional struct
%       .r_range
%       .lambda1, .lambda2, .lambda3
%       .weight_mode        = 'hadamard' | 'matrix'
%       .variance_source    = 'hat' | 'debiased'
%       .Gamma_hat          = {F x 1} biased precision estimates (needed when variance_source='hat')
%       .threshold_domain   = 'entry' | 'pcor'
%       .template_cell      = {F x 1} SPD base template; default = Gamma_stat
%       .score_mode         = 'hybrid' | 'higgs' | 'jspace'
%       .density_min, .density_max, .density_penalty_weight
%       .min_eig
%
% Outputs:
%   best_r     : selected Rayleigh threshold
%   mask_cell  : selected support masks
%   Var_proxies: variance proxies used in thresholding
%   stats      : score / density traces

    if nargin < 6 || isempty(params), params = struct(); end

    r_grid      = get_opt_(params, 'r_range', 0.7:0.1:3.16);
    lambda1     = get_opt_(params, 'lambda1', 0);
    lambda2     = get_opt_(params, 'lambda2', 0);
    lambda3     = get_opt_(params, 'lambda3', 0);
    mode        = get_opt_(params, 'weight_mode', 'hadamard');
    var_src     = get_opt_(params, 'variance_source', 'debiased');
    Gamma_hat   = get_opt_(params, 'Gamma_hat', []);
    thresh_dom  = get_opt_(params, 'threshold_domain', 'entry');
    template_cell = get_opt_(params, 'template_cell', []);
    score_mode  = lower(get_opt_(params, 'score_mode', 'hybrid'));
    min_eig     = get_opt_(params, 'min_eig', 1e-8);

    F = numel(Gamma_stat);
    Ksym = (K + K') / 2;

    if isempty(template_cell)
        template_cell = Gamma_stat;
    end

    % -------- Variance proxy --------
    Var_proxies = cell(F, 1);
    for f = 1:F
        if strcmpi(var_src, 'hat') && ~isempty(Gamma_hat)
            G_src = Gamma_hat{f};
        else
            G_src = Gamma_stat{f};
        end
        d = abs(diag(G_src));
        V = d * d' + abs(G_src).^2;
        Var_proxies{f} = real((V + V') / 2);
    end

    n_r = numel(r_grid);
    scores = -inf(n_r, 1);
    dens_curve = nan(n_r, 1);

    for i = 1:n_r
        r = r_grid(i);
        G_candidate = cell(F, 1);
        mask_current = cell(F, 1);

        for f = 1:F
            Gs = Gamma_stat{f};
            G0 = template_cell{f};
            p = size(Gs, 1);

            if strcmpi(thresh_dom, 'pcor')
                mask = build_mask_pcor_(Gs, r, n_samples);
            else
                V_proxy = Var_proxies{f};
                mask = build_mask_entry_(Gs, V_proxy, r, n_samples);
            end
            mask(1:p+1:end) = true;
            mask_current{f} = mask;

            Gm = G0;
            Gm(~mask) = 0;
            Gm(1:p+1:end) = real(diag(G0));

            Gm = utils_math.make_hermitian(Gm);
            [Gm, ~] = utils_math.project_spd(Gm, min_eig);
            G_candidate{f} = Gm;
        end

        dens_curve(i) = median_density_(mask_current);

        base_score = compute_score_(G_candidate, S_whitened, Ksym, W, ...
            lambda1, lambda2, lambda3, mode, score_mode);

        penalty = density_penalty_(dens_curve(i), params);
        scores(i) = base_score - penalty;
    end

    [best_score, best_idx] = max(scores);
    best_r = r_grid(best_idx);

    % Rebuild the best mask cleanly (no complex() / sqrt(complex()) terms)
    mask_cell = cell(F, 1);
    for f = 1:F
        Gs = Gamma_stat{f};
        p = size(Gs, 1);

        if strcmpi(thresh_dom, 'pcor')
            mask = build_mask_pcor_(Gs, best_r, n_samples);
        else
            V_proxy = Var_proxies{f};
            mask = build_mask_entry_(Gs, V_proxy, best_r, n_samples);
        end
        mask(1:p+1:end) = true;
        mask_cell{f} = mask;
    end

    stats = struct();
    stats.r_grid       = r_grid;
    stats.scores       = scores;
    stats.best_score   = best_score;
    stats.best_idx     = best_idx;
    stats.density_curve = dens_curve;
end

% ============================================================
% helpers
% ============================================================

function mask = build_mask_entry_(G, V_proxy, r, n_samples)
    n_eff = max(real(n_samples), eps);
    Th = (r / sqrt(n_eff)) * sqrt(max(V_proxy, 0));
    mask = abs(G) >= Th;
end

function mask = build_mask_pcor_(G, r, n_samples)
    p = size(G,1);
    d = real(diag(G));
    d = max(d, 1e-12);
    P = -G ./ sqrt(d * d.');
    P(1:p+1:end) = 0;
    Th = r / sqrt(max(real(n_samples), eps));
    mask = abs(P) >= Th;
end

function pen = density_penalty_(den_med, params)
    pen = 0;
    has_fields = isfield(params, 'density_min') && isfield(params, 'density_max') && isfield(params, 'density_penalty_weight');
    if ~has_fields || ~isfinite(den_med), return; end

    dmin = params.density_min;
    dmax = params.density_max;
    wpen = params.density_penalty_weight;

    if den_med < dmin
        pen = wpen * ((dmin - den_med) / max(dmin, eps))^2;
    elseif den_med > dmax
        pen = wpen * ((den_med - dmax) / max(dmax, eps))^2;
    end
end

function den_med = median_density_(mask_cell)
    F = numel(mask_cell);
    p = size(mask_cell{1},1);
    den = zeros(F,1);
    for f = 1:F
        M = mask_cell{f};
        M(1:p+1:end) = false;
        den(f) = (nnz(M) / 2) / (p * (p - 1) / 2);
    end
    den_med = median(den);
end

function score = compute_score_(Gamma_cell, S_cell, Ksym, W, lambda1, lambda2, lambda3, mode, score_mode)
    F = numel(Gamma_cell);

    % --- likelihood term ---
    loglik = 0;
    for f = 1:F
        G = Gamma_cell{f};
        S = S_cell{f};
        [ld, valid] = utils_math.safe_log_det(G);
        if ~valid
            score = -inf;
            return;
        end
        loglik = loglik + (ld - real(trace(S * G)));
    end

    % --- optional l1 term (HIGGS-like search score) ---
    l1_pen = 0;
    if lambda2 > 0
        for f = 1:F
            G = Gamma_cell{f};
            p = size(G,1);
            off = true(p);
            off(1:p+1:end) = false;
            l1_pen = l1_pen + sum(abs(W(off) .* G(off)));
        end
    end

    % --- cross-frequency smoothness ---
    freq_pen = 0;
    if lambda1 > 0
        for f1 = 1:F
            for f2 = (f1+1):F
                if Ksym(f1, f2) == 0, continue; end
                D = Gamma_cell{f1} - Gamma_cell{f2};
                if strcmpi(mode, 'matrix')
                    term = real(trace(D' * (W * D)));
                else
                    term = real(sum(sum(conj(D) .* (W .* D))));
                end
                freq_pen = freq_pen + Ksym(f1, f2) * term;
            end
        end
    end

    % --- within-frequency quadratic term ---
    quad_pen = 0;
    if lambda3 > 0
        for f = 1:F
            G = Gamma_cell{f};
            if strcmpi(mode, 'matrix')
                term = real(trace(G' * (W * G)));
            else
                term = real(sum(sum(conj(G) .* (W .* G))));
            end
            quad_pen = quad_pen + term;
        end
    end

    switch score_mode
        case 'higgs'
            score = loglik - lambda2 * l1_pen;
        case 'jspace'
            score = loglik - lambda1 * freq_pen - lambda3 * quad_pen;
        otherwise % hybrid
            score = loglik - lambda1 * freq_pen - lambda2 * l1_pen - lambda3 * quad_pen;
    end
end

function val = get_opt_(s, f, d)
    if isfield(s, f), val = s.(f); else, val = d; end
end
