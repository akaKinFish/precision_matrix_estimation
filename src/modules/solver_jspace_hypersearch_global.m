function [Omega_best, Sigma_best, outs_best] = solver_jspace_hypersearch_global(Svv_in, L, GraphLaplacian, cfg)
% SOLVER_JSPACE_HYPERSEARCH_GLOBAL Outer-loop hyperparameter search for J-SPACE.
%
% This wrapper performs one warm-start E-step to build a global threshold
% and hyperparameter grid, then runs inner-loop EM with fixed lambdas.

if nargin < 4, cfg = struct(); end
if nargin < 3, GraphLaplacian = []; end

% ============================================================
% 1. Normalize inputs
% ============================================================
Svv_cell = coerce_cell_(Svv_in);
F = numel(Svv_cell);
[Ns, Nr] = size(L);

VERBOSE = get_cfg(cfg, 'verbose', true);
GRID_SIZE = get_cfg(cfg, 'grid_size', 20);
USE_PARFOR = get_cfg(cfg, 'use_parfor', false);
SEARCH_METRIC = lower(get_cfg(cfg, 'search_metric', 'ebic'));
HYPER_MODE = lower(get_cfg(cfg, 'hyper_grid_mode', 'scale_from_sbar'));

alpha2_grid = get_cfg(cfg, 'alpha2_grid', logspace(-3, 0, GRID_SIZE));
alpha1_grid = get_cfg(cfg, 'alpha1_grid', [0.1 0.3 1]);
alpha3_grid = get_cfg(cfg, 'alpha3_grid', [0.1 0.3 1]);

if isempty(alpha1_grid), alpha1_grid = 0.1; end
if isempty(alpha2_grid), alpha2_grid = logspace(-3, 0, GRID_SIZE); end
if isempty(alpha3_grid), alpha3_grid = 0.1; end

alpha1_grid = unique(alpha1_grid(:).');
alpha2_grid = unique(alpha2_grid(:).');
alpha3_grid = unique(alpha3_grid(:).');

W_gamma = get_cfg(cfg, 'weight_matrix', eye(Nr));
if ~ismatrix(W_gamma) || any(size(W_gamma) ~= [Nr Nr])
    W_gamma = eye(Nr);
end

if isempty(GraphLaplacian)
    if VERBOSE, fprintf('[J-SPACE] No Laplacian provided. Spatial smoothing disabled.\n'); end
    GraphLaplacian = zeros(Nr, Nr, 'like', L);
end

% ============================================================
% 2. Warm-start E-step + whitening (single pass)
% ============================================================
Sigma_source_init = init_sigma_source_(Svv_cell, L, cfg, VERBOSE);

tr_S = 0;
for f = 1:F, tr_S = tr_S + trace(Svv_cell{f}); end
noise_cov = (tr_S / (F * Ns)) * 0.05 * eye(Ns, 'like', Svv_cell{1});

[Psijj_cell, ~] = module2_estep(Svv_cell, L, Sigma_source_init, noise_cov);
[Sjj_tilde_init, ~, ~] = module1_data_whitening(Psijj_cell, 'smoothing_window', 1);

% ============================================================
% 3. Global thresholds from warm-start Sjj_tilde
% ============================================================
vals = collect_offdiag_abs_(Sjj_tilde_init);
[t_active, meta_active] = utils_stats_gmm_threshold_1d(vals);
if isempty(vals)
    q95 = 0;
else
    q95 = quantile(vals, 0.95);
end
t_rescue = max(t_active, q95);
meta_rescue = struct('fallback', false, 'shared_with_active', false, 'quantile_95', q95);

% ============================================================
% 4. Build hyperparameter grid from warm-start scale
% ============================================================
Sbar = zeros(Nr, Nr);
for f = 1:F
    Sbar = Sbar + gather(Sjj_tilde_init{f});
end
Sbar = utils_math.make_hermitian(Sbar / F);

lambda2_max = compute_lambda2_max_(Sbar, W_gamma);
[lambda1_scale, lambda3_scale] = compute_scales_(Sjj_tilde_init, Sbar);

lambda1_fixed = get_cfg(cfg, 'lambda1_fixed', []);
lambda3_fixed = get_cfg(cfg, 'lambda3_fixed', []);
if ~isempty(lambda1_fixed), lambda1_fixed = lambda1_fixed(1); end
if ~isempty(lambda3_fixed), lambda3_fixed = lambda3_fixed(1); end

n1 = numel(alpha1_grid);
n2 = numel(alpha2_grid);
n3 = numel(alpha3_grid);
n_combo = n1 * n2 * n3;

alpha1_all = zeros(n_combo, 1);
alpha2_all = zeros(n_combo, 1);
alpha3_all = zeros(n_combo, 1);
lambda1_all = zeros(n_combo, 1);
lambda2_all = zeros(n_combo, 1);
lambda3_all = zeros(n_combo, 1);

for idx = 1:n_combo
    [i1, i2, i3] = ind2sub([n1, n2, n3], idx);
    a1 = alpha1_grid(i1);
    a2 = alpha2_grid(i2);
    a3 = alpha3_grid(i3);
    lam2 = a2 * lambda2_max;

    if strcmpi(HYPER_MODE, 'ratio_to_lambda2')
        if ~isempty(lambda1_fixed), lam1 = lambda1_fixed; else, lam1 = a1 * lam2; end
        if ~isempty(lambda3_fixed), lam3 = lambda3_fixed; else, lam3 = a3 * lam2; end
    else
        if ~isempty(lambda1_fixed), lam1 = lambda1_fixed; else, lam1 = a1 * lambda1_scale; end
        if ~isempty(lambda3_fixed), lam3 = lambda3_fixed; else, lam3 = a3 * lambda3_scale; end
    end

    if F <= 1
        lam3 = 0;
    end

    alpha1_all(idx) = a1;
    alpha2_all(idx) = a2;
    alpha3_all(idx) = a3;
    lambda1_all(idx) = lam1;
    lambda2_all(idx) = lam2;
    lambda3_all(idx) = lam3;
end

if VERBOSE
    fprintf('[Search] t_active=%.2e  t_rescue=%.2e\n', t_active, t_rescue);
    fprintf('[Search] lambda2_max=%.2e  alpha2=[%.1e, %.1e]\n', ...
        lambda2_max, min(alpha2_grid), max(alpha2_grid));
    fprintf('[Search] lambda2 range=[%.2e, %.2e]  combos=%d\n', ...
        min(lambda2_all), max(lambda2_all), n_combo);
    fprintf('[Search] alpha1=[%.1e, %.1e]  alpha3=[%.1e, %.1e]  mode=%s\n', ...
        min(alpha1_grid), max(alpha1_grid), min(alpha3_grid), max(alpha3_grid), HYPER_MODE);
    fprintf('[Search] lambda1 range=[%.2e, %.2e]  lambda3 range=[%.2e, %.2e]\n', ...
        min(lambda1_all), max(lambda1_all), min(lambda3_all), max(lambda3_all));
end

% ============================================================
% 5. Outer-loop search
% ============================================================
score = NaN(n_combo, 1);
hist_aic = NaN(n_combo, 1);
hist_ebic = NaN(n_combo, 1);
hist_nll = NaN(n_combo, 1);
hist_full = NaN(n_combo, 1);
hist_den = NaN(n_combo, 1);
hist_edges = NaN(n_combo, 1);
loglik_curve = cell(n_combo, 1);

cfg_inner = cfg;
cfg_inner.hyperparam_mode = 'external_fixed';
cfg_inner.t_active = t_active;
cfg_inner.t_rescue = t_rescue;
cfg_inner.init_sigma_source = Sigma_source_init;
cfg_inner.init_gamma = get_cfg(cfg, 'init_gamma', []);
if isempty(cfg_inner.init_gamma)
    gamma_init = cell(F, 1);
    for f = 1:F, gamma_init{f} = eye(Nr, 'like', L); end
    cfg_inner.init_gamma = gamma_init;
end
cfg_inner.verbose = get_cfg(cfg, 'inner_verbose', false);
cfg_inner.debug_print = get_cfg(cfg, 'inner_debug_print', false);
cfg_inner.plot = false;

if USE_PARFOR
    parfor idx = 1:n_combo
        [score(idx), hist_aic(idx), hist_ebic(idx), hist_nll(idx), hist_full(idx), ...
            hist_den(idx), hist_edges(idx), loglik_curve{idx}] = run_combo_( ...
            Svv_cell, L, GraphLaplacian, cfg_inner, SEARCH_METRIC, ...
            lambda1_all(idx), lambda2_all(idx), lambda3_all(idx), Nr);
    end
else
    for idx = 1:n_combo
        [score(idx), hist_aic(idx), hist_ebic(idx), hist_nll(idx), hist_full(idx), ...
            hist_den(idx), hist_edges(idx), loglik_curve{idx}] = run_combo_( ...
            Svv_cell, L, GraphLaplacian, cfg_inner, SEARCH_METRIC, ...
            lambda1_all(idx), lambda2_all(idx), lambda3_all(idx), Nr);
        if VERBOSE && mod(idx, max(1, floor(n_combo / 10))) == 0
            fprintf('  [Search] Progress %d/%d\n', idx, n_combo);
        end
    end
end

[best_idx, best_score, select_note] = select_best_combo_(score, hist_full, hist_nll, hist_den, SEARCH_METRIC);
if VERBOSE && ~isempty(select_note)
    fprintf('  [Search] %s\n', select_note);
end

if VERBOSE
    fprintf('[Search] Best score (%s)=%.2e at idx=%d\n', upper(SEARCH_METRIC), best_score, best_idx);
    fprintf('[Search] Selected lambdas: lambda1=%.2e lambda2=%.2e lambda3=%.2e\n', ...
        lambda1_all(best_idx), lambda2_all(best_idx), lambda3_all(best_idx));
    [~, order] = sort(score, 'ascend');
    topk = order(1:min(5, numel(order)));
    fprintf('[Search] Top-%d combos by %s:\n', numel(topk), upper(SEARCH_METRIC));
    for ii = 1:numel(topk)
        k = topk(ii);
        fprintf('  #%d: lambda1=%.2e lambda2=%.2e lambda3=%.2e | score=%.2e | den=%.2f%%\n', ...
            ii, lambda1_all(k), lambda2_all(k), lambda3_all(k), score(k), hist_den(k)*100);
    end
end

% ============================================================
% 6. Final run with the selected hyperparameters
% ============================================================
cfg_final = cfg_inner;
cfg_final.verbose = VERBOSE;
cfg_final.debug_print = get_cfg(cfg, 'debug_print', false);
cfg_final.lambda1 = lambda1_all(best_idx);
cfg_final.lambda2 = lambda2_all(best_idx);
cfg_final.lambda3 = lambda3_all(best_idx);

[Omega_best, Sigma_best, outs_best] = solver_jspace_adaptive(Svv_cell, L, GraphLaplacian, cfg_final);

% ============================================================
% 7. Attach search diagnostics
% ============================================================
search_history = repmat(struct( ...
    'alpha1', [], 'alpha2', [], 'alpha3', [], ...
    'lambda1', [], 'lambda2', [], 'lambda3', [], ...
    'score', [], 'score_metric', '', ...
    'density', [], 'edges', [], ...
    'aic', [], 'ebic', [], 'nll', [], 'full_obj', [], ...
    'loglik_curve', []), n_combo, 1);

for idx = 1:n_combo
    search_history(idx).alpha1 = alpha1_all(idx);
    search_history(idx).alpha2 = alpha2_all(idx);
    search_history(idx).alpha3 = alpha3_all(idx);
    search_history(idx).lambda1 = lambda1_all(idx);
    search_history(idx).lambda2 = lambda2_all(idx);
    search_history(idx).lambda3 = lambda3_all(idx);
    search_history(idx).score = score(idx);
    search_history(idx).score_metric = SEARCH_METRIC;
    search_history(idx).density = hist_den(idx);
    search_history(idx).edges = hist_edges(idx);
    search_history(idx).aic = hist_aic(idx);
    search_history(idx).ebic = hist_ebic(idx);
    search_history(idx).nll = hist_nll(idx);
    search_history(idx).full_obj = hist_full(idx);
    search_history(idx).loglik_curve = loglik_curve{idx};
end

outs_best.search_history = search_history;
outs_best.thresholds = struct('t_active', t_active, 't_rescue', t_rescue, ...
    'meta_active', meta_active, 'meta_rescue', meta_rescue);
outs_best.search_meta = struct('alpha1_grid', alpha1_grid, 'alpha2_grid', alpha2_grid, ...
    'alpha3_grid', alpha3_grid, 'lambda2_max', lambda2_max, ...
    'lambda1_scale', lambda1_scale, 'lambda3_scale', lambda3_scale, ...
    'mode', HYPER_MODE, 'metric', SEARCH_METRIC);
outs_best.best_combo = struct('index', best_idx, 'score', best_score, ...
    'lambda1', lambda1_all(best_idx), 'lambda2', lambda2_all(best_idx), 'lambda3', lambda3_all(best_idx), ...
    'alpha1', alpha1_all(best_idx), 'alpha2', alpha2_all(best_idx), 'alpha3', alpha3_all(best_idx));

if get_cfg(cfg, 'plot', false)
    plot_search_results_(search_history, best_idx, SEARCH_METRIC);
end
end

% ============================================================
% Helpers
% ============================================================
function Svv_cell = coerce_cell_(Svv_in)
    if iscell(Svv_in)
        Svv_cell = Svv_in(:);
        return;
    end
    if ndims(Svv_in) == 3
        F = size(Svv_in, 3);
        Svv_cell = cell(F, 1);
        for f = 1:F, Svv_cell{f} = Svv_in(:, :, f); end
    else
        Svv_cell = {Svv_in};
    end
end

function val = get_cfg(s, f, d)
    if isfield(s, f), val = s.(f); else, val = d; end
end

function Sigma_source = init_sigma_source_(Svv_cell, L, cfg, verbose)
    Sigma_source = [];
    if isfield(cfg, 'init_sigma_source') && ~isempty(cfg.init_sigma_source)
        Sigma_source = cfg.init_sigma_source;
        if ~iscell(Sigma_source)
            if ndims(Sigma_source) == 3
                F = size(Sigma_source, 3);
                tmp = cell(F, 1);
                for f = 1:F, tmp{f} = Sigma_source(:, :, f); end
                Sigma_source = tmp;
            else
                Sigma_source = {Sigma_source};
            end
        end
        return;
    end

    if verbose, fprintf('[Init] Running eLORETA for warm-start...\n'); end
    F = numel(Svv_cell);
    [Ns, ~] = size(Svv_cell{1});

    S_avg = zeros(Ns, Ns, 'like', Svv_cell{1});
    for f = 1:F, S_avg = S_avg + Svv_cell{f}; end
    S_avg = S_avg / F;

    [T_eloreta, ~] = run_eloreta_core_(L, S_avg, 0.05);

    Sigma_source = cell(F, 1);
    for f = 1:F
        Svv_f = Svv_cell{f};
        Sjj_eloreta = T_eloreta * Svv_f * T_eloreta';
        Sjj_eloreta = utils_math.make_hermitian(Sjj_eloreta);
        [Sjj_spd, ~] = utils_math.project_spd(Sjj_eloreta, 1e-8);
        Sigma_source{f} = Sjj_spd;
    end
end

function [score, aic, ebic, nll, full_obj, den, edges, loglik_curve] = run_combo_( ...
    Svv_cell, L, GraphLaplacian, cfg_inner, metric, lambda1, lambda2, lambda3, p)

    cfg_local = cfg_inner;
    cfg_local.lambda1 = lambda1;
    cfg_local.lambda2 = lambda2;
    cfg_local.lambda3 = lambda3;

    [~, ~, outs] = solver_jspace_adaptive(Svv_cell, L, GraphLaplacian, cfg_local);

    if isempty(outs.grid_history)
        aic = Inf; ebic = Inf; nll = Inf; full_obj = Inf; den = NaN; edges = NaN;
    else
        hist = outs.grid_history(end);
        aic = hist.aic;
        ebic = hist.ebic;
        nll = hist.nll;
        full_obj = hist.full_obj;
        den = hist.density;
        if isfield(hist, 'edges_total') && ~isempty(hist.edges_total)
            edges = hist.edges_total;
            if numel(edges) > 1
                edges = edges(end);
            end
        else
            edges = den * (p * (p - 1) / 2);
        end
    end
    loglik_curve = outs.loglik;

    switch metric
        case 'aic'
            score = aic;
        case 'ebic'
            score = ebic;
        case 'bic'
            score = ebic;
        case 'full_obj'
            score = full_obj;
        case 'nll'
            score = nll;
        otherwise
            score = ebic;
    end
end

function [best_idx, best_score, note] = select_best_combo_(score, full_obj, nll, den, metric_name)
    note = '';
    best_idx = 1;
    best_score = Inf;

    valid_score = isfinite(score);
    if any(valid_score)
        candidates = find(valid_score);
        [best_score, rel] = min(score(candidates));
        best_idx = candidates(rel);
        return;
    end

    valid_full = isfinite(full_obj);
    if any(valid_full)
        candidates = find(valid_full);
        [best_score, rel] = min(full_obj(candidates));
        best_idx = candidates(rel);
        note = sprintf('%s invalid; selecting by finite full objective.', upper(metric_name));
        return;
    end

    valid_nll = isfinite(nll);
    if any(valid_nll)
        candidates = find(valid_nll);
        [best_score, rel] = min(nll(candidates));
        best_idx = candidates(rel);
        note = sprintf('%s invalid; selecting by finite minus_2_ll.', upper(metric_name));
        return;
    end

    [~, best_idx] = min(den);
    best_score = score(best_idx);
    note = sprintf('%s invalid; selecting lowest-density candidate.', upper(metric_name));
end

function vals = collect_offdiag_abs_(S_cell)
    if ~iscell(S_cell), S_cell = {S_cell}; end
    F = numel(S_cell);
    p = size(S_cell{1}, 1);
    mask = tril(true(p), -1);
    vals = [];
    for f = 1:F
        Sf = gather(S_cell{f});
        v = abs(Sf(mask));
        vals = [vals; v(:)]; %#ok<AGROW>
    end
end

function lambda2_max = compute_lambda2_max_(Sbar, W)
    p = size(Sbar, 1);
    if p < 2
        lambda2_max = 0;
        return;
    end
    mask = tril(true(p), -1);
    S_off = abs(Sbar(mask));
    W_off = abs(W(mask));
    denom = max(W_off, eps);
    ratios = S_off ./ denom;
    if isempty(ratios) || ~any(isfinite(ratios))
        lambda2_max = 0;
    else
        lambda2_max = max(ratios(isfinite(ratios)));
    end
end

function [scale1, scale3] = compute_scales_(Sjj_tilde, Sbar)
    if ~iscell(Sjj_tilde), Sjj_tilde = {Sjj_tilde}; end
    F = numel(Sjj_tilde);
    p = size(Sjj_tilde{1}, 1);
    I = eye(p);

    scale1 = norm(Sbar - I, 'fro');
    if ~isfinite(scale1), scale1 = 0; end

    if F <= 1
        scale3 = 0;
    else
        diffs = zeros(F, 1);
        for f = 1:F
            Sf = utils_math.make_hermitian(gather(Sjj_tilde{f}));
            diffs(f) = norm(Sf - Sbar, 'fro');
        end
        scale3 = median(diffs);
    end
    if ~isfinite(scale3), scale3 = 0; end
end

function plot_search_results_(history, best_idx, metric)
    if isempty(history), return; end
    lam2 = [history.lambda2];
    score = [history.score];
    loglik = {history.loglik_curve};

    figure('Name', 'J-SPACE Hyperparameter Search', 'Color', 'w', 'Position', [100, 100, 1200, 500]);
    subplot(1, 2, 1);
    if all(lam2 == 0)
        plot(score, 'ko-'); xlabel('Combo Index'); ylabel(upper(metric));
    else
        semilogx(lam2, score, 'ko'); xlabel('Lambda2'); ylabel(upper(metric));
        grid on;
    end
    title('Score vs Lambda2');

    subplot(1, 2, 2);
    [~, order] = sort(score, 'ascend');
    topk = order(1:min(5, numel(order)));
    hold on;
    for k = 1:numel(topk)
        idx = topk(k);
        if isempty(loglik{idx}), continue; end
        plot(loglik{idx}, 'LineWidth', 1.2);
    end
    if best_idx >= 1 && best_idx <= numel(loglik) && ~isempty(loglik{best_idx})
        plot(loglik{best_idx}, 'k-', 'LineWidth', 2.0);
    end
    xlabel('EM Iteration');
    ylabel('Log-likelihood');
    title('Top-5 EM Convergence');
    grid on;
end

function [T, W] = run_eloreta_core_(L, Svv, regu)
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
