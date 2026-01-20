%TEST_JSPACE_PSPLINE_SMOKE  Fast smoke test for J-SPACE P-spline pipeline.

rootDir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(rootDir, 'JSPACE_Spline_v1'));

rng(0);

Ne = 6;
N = 8;
F = 4;
freq = linspace(1, 40, F).';

Svv = zeros(Ne, Ne, F);
for f = 1:F
    X = randn(Ne) + 1i * randn(Ne);
    S = X * X';
    S = 0.5 * (S + S');
    S(1:Ne+1:end) = real(diag(S));
    Svv(:, :, f) = S;
end

L = randn(Ne, N);
dwi_C = rand(N);
dwi_C = 0.5 * (dwi_C + dwi_C');
dwi_C(1:N+1:end) = 0;

cfg = struct();
cfg.em.max_iter = 2;
cfg.em.verbose = false;
cfg.mstep.fista.max_iter = 8;
cfg.mstep.fista.verbose = false;

[Omega, Sjj, outs] = run_jspace_pspline_real(Svv, L, freq, dwi_C, cfg);

assert(numel(Omega) == F, 'Omega length mismatch.');
assert(numel(Sjj) == F, 'Sjj length mismatch.');
assert(isfield(outs, 'obj_mstep'), 'Missing output field obj_mstep.');

for f = 1:numel(Omega)
    G = 0.5 * (Omega{f} + Omega{f}');
    [~, p] = chol(G);
    assert(p == 0, 'Omega is not SPD at frequency %d.', f);
end
