function test_prox_fista()
%TEST_PROX_FISTA Basic sanity check for module_proximal_operator_fista.

    if exist('module_proximal_operator_fista', 'file') ~= 2
        here = fileparts(mfilename('fullpath'));
        addpath(genpath(fullfile(here, '..', 'src')));
        addpath(genpath(fullfile(here, '..', 'utils')));
    end

    rng(42);
    p = 8;
    % Random Hermitian PSD-ish matrix
    A = randn(p) + 1i * randn(p);
    G0 = utils_math.make_hermitian(A * A' + 0.5 * eye(p));

    % Active mask: keep upper-left block + diagonal
    active_mask = false(p);
    active_mask(1:5, 1:5) = true;
    active_mask = active_mask | eye(p);

    tau = 0.2;
    opts.min_eig = 1e-6;
    opts.enforce_unit_diagonal = true;

    [G_new, info] = module_proximal_operator_fista.compute(G0, tau, active_mask, opts);

    % Checks
    assert(norm(G_new - G_new', 'fro') < 1e-10 + 1e-6 * norm(G_new, 'fro'), 'Not Hermitian');
    eigvals = eig((G_new + G_new')/2);
    assert(min(real(eigvals)) > opts.min_eig * 0.5, 'SPD floor violated');
    offmask = ~active_mask;
    assert(all(abs(G_new(offmask)) < 1e-8), 'Active mask not enforced');
    d = real(diag(G_new));
    assert(max(abs(d - 1)) < 1e-3, 'Unit diagonal enforcement failed');

    fprintf('test_prox_fista passed. min eig %.2e (clipped=%d)\n', info.min_eig, info.clipped);
end
