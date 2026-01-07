function demo_fista_cells()
%DEMO_FISTA_CELLS Minimal sanity run of fista_with_backtracking_cells on cell matrices.

    addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src')));
    addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils')));

    rng(7);
    F = 2; p = 5;
    Sigmas = cell(F,1); Gamma0 = cell(F,1);
    for f = 1:F
        A = randn(p) + 1i*randn(p);
        Sigmas{f} = utils_math.make_hermitian(A*A' + 0.5*eye(p));
        Gamma0{f}  = eye(p);
    end

    Kernel = eye(F);
    W = eye(p);
    params.lambda1 = 0.01; params.lambda3 = 0.01; params.weight_mode = 'hadamard';
    params.min_eig = 1e-6; params.penalize_diagonal = false; params.enforce_unit_diagonal = false;

    f = @(G) module_objective.compute(G, Sigmas, Kernel, W, setfield(params,'lambda2',0)); %#ok<SFLD>
    grad_f = @(G) module_gradient.compute(G, Sigmas, Kernel, W, setfield(params,'lambda2',0));
    g = @(G) l1_off(G, params);
    prox = @(X, tau) cell_prox_l1_off(X, tau, params);

    lambda = 0.1;
    L0 = 1; eta = 1.5; max_iter = 10; tol = 1e-6; max_bt = 20;
    [G_hat, hist_vals, smoothf] = fista_with_backtracking_cells(f, grad_f, g, prox, Gamma0, lambda, L0, eta, max_iter, tol, max_bt);

    fprintf('History: ');
    fprintf('%.3e ', hist_vals);
    fprintf('\nFinal smooth f = %.3e\n', smoothf);
    assert(isscalar(hist_vals(1)), 'History elements must be scalar');
    assert(hist_vals(end) <= hist_vals(1) + 1e-9, 'Objective did not decrease');
    for fidx = 1:F
        Gh = utils_math.make_hermitian(G_hat{fidx});
        eigvals = eig(Gh);
        assert(min(real(eigvals)) > params.min_eig * 0.5, 'SPD floor violated');
    end
    fprintf('demo_fista_cells passed.\n');
end

function val = l1_off(Gamma_cells, params)
    penalize_diag = isfield(params, 'penalize_diagonal') && params.penalize_diagonal;
    F = numel(Gamma_cells);
    p = size(Gamma_cells{1}, 1);
    val = 0;
    for f = 1:F
        G = Gamma_cells{f};
        if ~penalize_diag, G(1:p+1:end) = 0; end
        val = val + sum(abs(G(:)));
    end
end

function Xp = cell_prox_l1_off(X, tau, params)
    F = numel(X);
    Xp = cell(F,1);
    for f = 1:F
        opts.penalize_diagonal = isfield(params,'penalize_diagonal') && params.penalize_diagonal;
        opts.min_eig = params.min_eig;
        opts.enforce_unit_diagonal = isfield(params,'enforce_unit_diagonal') && params.enforce_unit_diagonal;
        Xp{f} = module_proximal_operator_fista.compute(X{f}, tau, [], opts);
    end
end
