function test_jspace_adaptive_smoke()
%TEST_JSPACE_ADAPTIVE_SMOKE Fast smoke test for solver_jspace_hypersearch_global.

    if exist('solver_jspace_hypersearch_global', 'file') ~= 2
        here = fileparts(mfilename('fullpath'));
        addpath(genpath(fullfile(here, '..', 'src')));
        addpath(genpath(fullfile(here, '..', 'utils')));
    end

    rng(7);
    Ns = 6; Nr = 5; F = 2;

    L = randn(Ns, Nr);

    Omega_true = eye(Nr);
    Omega_true(1,2) = 0.20; Omega_true(2,1) = 0.20;
    Omega_true(3,4) = -0.15; Omega_true(4,3) = -0.15;
    Omega_true(2,5) = 0.10; Omega_true(5,2) = 0.10;
    for i = 1:Nr
        Omega_true(i,i) = max(Omega_true(i,i), sum(abs(Omega_true(i,[1:i-1 i+1:end]))) + 0.5);
    end
    Sigma_true = inv(Omega_true);

    Svv_cell = cell(F, 1);
    for f = 1:F
        jitter = 0.02 * randn(Nr);
        jitter = (jitter + jitter') / 2;
        Sigma_f = Sigma_true + jitter;
        [Sigma_f, ~] = utils_math.project_spd(Sigma_f, 1e-6);
        Svv = L * Sigma_f * L';
        noise = 0.01 * trace(Svv) / Ns;
        Svv = Svv + noise * eye(Ns);
        Svv_cell{f} = (Svv + Svv') / 2;
    end

    GraphLaplacian = zeros(Nr);

    cfg = struct();
    cfg.max_em_iter = 2;
    cfg.grid_size = 4;
    cfg.hyper_grid_mode = 'ratio_to_lambda2';
    cfg.alpha1_grid = [0.1 0.3];
    cfg.alpha3_grid = [0 0.1];
    cfg.use_parfor = false;
    cfg.use_gpu = false;
    cfg.verbose = false;
    cfg.inner_verbose = false;
    cfg.debug_print = false;
    cfg.m_samples = 200;
    cfg.plot = false;

    [Omega_est, Sigma_src_est, outs] = solver_jspace_hypersearch_global(Svv_cell, L, GraphLaplacian, cfg);

    assert(iscell(Omega_est) && numel(Omega_est) == F, 'Omega_est size mismatch');
    assert(iscell(Sigma_src_est) && numel(Sigma_src_est) == F, 'Sigma_src_est size mismatch');
    for f = 1:F
        assert(isequal(size(Omega_est{f}), [Nr Nr]), 'Omega_est matrix size mismatch');
        assert(isequal(size(Sigma_src_est{f}), [Nr Nr]), 'Sigma_src_est matrix size mismatch');
    end

    assert(isfield(outs, 'search_history') && ~isempty(outs.search_history), 'search_history missing');
    assert(isfield(outs, 'best_combo') && isfield(outs.best_combo, 'lambda2'), 'best_combo missing');
    assert(isfield(outs, 'thresholds'), 'thresholds missing');
    assert(isfield(outs.thresholds, 't_active') && isfinite(outs.thresholds.t_active), 't_active missing');
    assert(isfield(outs.thresholds, 't_rescue') && isfinite(outs.thresholds.t_rescue), 't_rescue missing');
    scores = [outs.search_history.score];
    assert(any(isfinite(scores)), 'Non-finite search scores');

    max_off = zeros(F, 1);
    max_diag = zeros(F, 1);
    for f = 1:F
        Om = Omega_est{f};
        max_diag(f) = max(abs(diag(Om)));
        Off = Om; Off(1:Nr+1:end) = 0;
        max_off(f) = max(abs(Off(:)));
    end
    rel_thresh = 1e-3;
    assert(any(max_off > rel_thresh * max_diag), 'Off-diagonal collapse detected');

    fprintf('test_jspace_adaptive_smoke: PASS\n');
end
