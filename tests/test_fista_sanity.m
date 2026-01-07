function test_fista_sanity()
%TEST_FISTA_SANITY Quick checks for module5_fista_main on small cell data.

    if exist('module5_fista_main', 'file') ~= 2
        here = fileparts(mfilename('fullpath'));
        addpath(genpath(fullfile(here, '..', 'src')));
        addpath(genpath(fullfile(here, '..', 'utils')));
    end

    rng(0);
    F = 2; p = 6;
    Sigmas = cell(F, 1);
    for f = 1:F
        A = randn(p) + 1i * randn(p);
        Sigmas{f} = utils_math.make_hermitian(A * A' + 0.3 * eye(p));
    end
    Kernel = eye(F);
    W = eye(p);

    % Active mask: forbid the last row/col (except its diagonal)
    base_mask = true(p);
    base_mask(p, :) = false; base_mask(:, p) = false;
    base_mask(1:p+1:end) = true;
    active_mask = cell(F, 1);
    for f = 1:F, active_mask{f} = base_mask; end

    input.whitened_covariances = Sigmas;
    input.smoothing_kernel = Kernel;
    input.weight_matrix = W;
    input.active_mask = active_mask;

    params.lambda1 = 0.05;
    params.lambda2 = 0.1;
    params.lambda3 = 0.02;
    params.min_eig = 1e-5;
    params.max_iter = 60;
    params.tol = 1e-6;
    params.weight_mode = 'hadamard';
    params.verbose = false;

    [Gamma, res] = module5_fista_main(input, params);

    assert(~isempty(res.objective_history), 'Objective history empty');
    diffs = diff(res.objective_history);
    tol_up = max(1e-8, abs(res.objective_history(1)) * 1e-9);
    assert(all(diffs <= tol_up), 'Objective not monotone decreasing');

    for f = 1:F
        G = utils_math.make_hermitian(Gamma{f});
        eigvals = eig(G);
        assert(min(real(eigvals)) > params.min_eig * 0.5, 'SPD floor violated');
        assert(norm(G - G', 'fro') < 1e-10 + 1e-6 * norm(G, 'fro'), 'Hermitian drift detected');
        offmask = ~active_mask{f};
        assert(all(abs(G(offmask)) < 1e-6), 'Active mask not enforced');
    end

    params_unit = params;
    params_unit.unit_diagonal = true;
    params_unit.max_iter = 15;
    [Gamma_unit, ~] = module5_fista_main(input, params_unit);
    for f = 1:F
        d = real(diag(Gamma_unit{f}));
        assert(max(abs(d - 1)) < 1e-3, 'Unit diagonal rescaling failed');
    end

    fprintf('FISTA sanity checks passed. Obj %.3e -> %.3e\n', res.objective_history(1), res.objective_history(end));
end
