function test_results = test_module7_leadfield()
% TEST_MODULE7_LEADFIELD - Validation suite for Module 7 leadfield-enabled simulator
%
% This test suite verifies both backward compatibility (no leadfield) and
% the forward-model paths (simple and spherical 3-layer). It focuses on
% dimensions, Hermitian/PD properties, SNR shaping, p>>n rank behavior,
% parameter validation, and determinism with seeds.
%
% Usage:
%   results = test_module7_leadfield();
%   fprintf('Success rate: %.1f%%\n', results.summary.success_rate*100);
%
% Author: <your name>
% Date  : 2025-09-04

fprintf('========================================\n');
fprintf('Testing Module 7 Leadfield (New API)\n');
fprintf('========================================\n\n');

tr = struct();
tr.timestamp = datestr(now);
tr.matlab_version = version();

%% A) Backward compatibility (leadfield OFF)
fprintf('A) Backward compatibility (leadfield OFF) ...\n');
tr.legacy = subtest_legacy_path();

%% B) Simple leadfield, default 19 sensors (10-20)
fprintf('\nB) Simple leadfield (19-ch 10-20) ...\n');
tr.simple19 = subtest_simple_19();

%% C) Spherical 3-layer leadfield
fprintf('\nC) Spherical 3-layer leadfield ...\n');
tr.spherical = subtest_spherical3layer();

%% D) p >> n behavior
fprintf('\nD) p >> n behavior ...\n');
tr.pggn = subtest_p_gg_n();

%% E) Complex Wishart path and Hermitian checks
fprintf('\nE) Complex Wishart & Hermitian checks ...\n');
tr.complexpath = subtest_complex_path();

%% F) SNR shaping (10 dB vs 20 dB)
fprintf('\nF) SNR shaping (trace-based) ...\n');
tr.snr = subtest_snr_shaping();

%% G) Parameter validation
fprintf('\nG) Parameter validation ...\n');
tr.paramval = subtest_param_validation();

%% H) Noise=0 fallback regularization path
fprintf('\nH) Zero-variance noise (regularization path) ...\n');
tr.noisezero = subtest_zero_variance_noise();

%% I) Determinism with random_seed
fprintf('\nI) Determinism with random_seed ...\n');
tr.seed = subtest_seed_determinism();

%% Summary
tr.summary = summarize_results(tr);
print_summary(tr.summary);

test_results = tr;
end

%% -------------------------------------------------------------------------
function r = subtest_legacy_path()
r = init_result('Legacy (no leadfield)');
try
    [Omega, Sigma, Sigma_emp, P] = module7_simulation_improved_complex( ...
        'n_nodes', 8, 'n_freq', 6, 'n_samples', 40, ...
        'generate_leadfield', false, 'random_seed', 42);

    assert(iscell(Omega) && numel(Omega)==6);
    assert(iscell(Sigma) && numel(Sigma)==6);
    assert(iscell(Sigma_emp) && numel(Sigma_emp)==6);

    for f=1:6
        assert(all(size(Omega{f})==[8,8]));
        assert(all(size(Sigma{f})==[8,8]));
        assert(all(size(Sigma_emp{f})==[8,8]));
        assert_hermitian(Omega{f}); assert_spd(Omega{f});
        assert_hermitian(Sigma{f}); assert_spd(Sigma{f});
        assert_hermitian(Sigma_emp{f});
    end
    assert(~isfield(P,'leadfield_matrix') || isempty(P.leadfield_matrix));
    assert(~isfield(P,'electrode_positions') || isempty(P.electrode_positions));

    r.success = true; r.details = 'Legacy path OK (dims, SPD, Hermitian).';
catch ME
    r = fail_result(r, ME);
end
end

function r = subtest_simple_19()
r = init_result('Simple leadfield (19 channels)');
try
    n=8; p=19; F=5; T=80;
    [Omega, Sigma, Sigma_emp, P] = module7_simulation_improved_complex( ...
        'n_nodes', n, 'n_sensors', p, 'n_freq', F, 'n_samples', T, ...
        'generate_leadfield', true, 'leadfield_type','simple', ...
        'sensor_snr_db', 15, 'random_seed', 123);

    % Dimensions
    assert(size(P.leadfield_matrix,1)==p && size(P.leadfield_matrix,2)==n);
    for f=1:F
        assert(all(size(Omega{f})==[n,n]));
        assert(all(size(Sigma{f})==[n,n]));
        assert(all(size(Sigma_emp{f})==[p,p]));
    end

    % L column norms ~ 1 (simple model列归一)
    cn = sqrt(sum(P.leadfield_matrix.^2,1));
    assert(all(abs(cn-1) < 1e-6));

    % Hermitian/SPD checks on sensor-domain covariances
    assert(isfield(P,'Sigma_vv_true') && isfield(P,'Sigma_vv_observed'));
    for f=1:F
        assert_hermitian(P.Sigma_vv_true{f});
        assert_hermitian(P.Sigma_vv_observed{f});
        assert_spd(P.Sigma_vv_observed{f});
    end

    % SNR trace check (~15 dB) —— 用观测减真值得到当频点噪声
    Sxi = P.Sigma_vv_observed{1} - P.Sigma_vv_true{1};
    snr_lin_target = 10^(15/10);
    [ok,msg] = check_trace_snr(P.Sigma_vv_true{1}, Sxi, snr_lin_target, 0.20);
    assert(ok, msg);

    r.success=true; r.details='Simple 19-ch leadfield OK (dims, norms, SNR).';
catch ME
    r = fail_result(r, ME);
end
end

function r = subtest_spherical3layer()
r = init_result('Spherical 3-layer leadfield');
try
    n=12; p=19; F=3; T=60;
    [~, ~, ~, P1] = module7_simulation_improved_complex( ...
        'n_nodes', n, 'n_sensors', p, 'n_freq', F, 'n_samples', T, ...
        'generate_leadfield', true, 'leadfield_type','simple', 'random_seed', 7);

    [~, ~, ~, P2] = module7_simulation_improved_complex( ...
        'n_nodes', n, 'n_sensors', p, 'n_freq', F, 'n_samples', T, ...
        'generate_leadfield', true, 'leadfield_type','spherical3layer', ...
        'head_radius', 0.10, 'conductivities', [0.33, 0.0042, 0.33], ...
        'random_seed', 7);

    Ls = P1.leadfield_matrix; L3 = P2.leadfield_matrix;
    assert(norm(Ls-L3,'fro')/max(norm(Ls,'fro'),1e-9) > 1e-2, '3-layer should differ from simple');

    % Electrode radii ~ head_radius
    r_e = vecnorm(P2.electrode_positions,2,2);
    assert(all(abs(r_e - 0.10) < 5e-3));

    % 当颅骨电导率升高，总幅值应上升（阈值稍放宽到 1.03 以抗数值误差）
    [~, ~, ~, P3] = module7_simulation_improved_complex( ...
        'n_nodes', n, 'n_sensors', p, 'n_freq', 2, 'n_samples', T, ...
        'generate_leadfield', true, 'leadfield_type','spherical3layer', ...
        'conductivities', [0.33, 0.33, 0.33], 'random_seed', 7);
    ratio = norm(P3.leadfield_matrix,'fro') / max(norm(P2.leadfield_matrix,'fro'),1e-12);
    assert(ratio > 1.03, 'Higher skull conductivity should increase overall leadfield magnitude');

    % 频点改变不应改变 L（几何不变）
    [~, ~, ~, P2b] = module7_simulation_improved_complex( ...
        'n_nodes', n, 'n_sensors', p, 'n_freq', 5, 'n_samples', T, ...
        'generate_leadfield', true, 'leadfield_type','spherical3layer', 'random_seed', 7);
    assert(norm(P2b.leadfield_matrix - L3, 'fro') < 1e-12);

    r.success=true; r.details='3-layer leadfield behaves as expected (geometry, conductivity).';
catch ME
    r = fail_result(r, ME);
end
end

function r = subtest_p_gg_n()
r = init_result('p >> n rank & SPD');
try
    n=16; p=64; F=4; T=120;
    [~, ~, ~, P] = module7_simulation_improved_complex( ...
        'n_nodes', n, 'n_sensors', p, 'n_freq', F, 'n_samples', T, ...
        'generate_leadfield', true, 'leadfield_type', 'spherical3layer', 'random_seed', 11);

    for f=1:F
        r_true = rank(P.Sigma_vv_true{f});
        r_obs  = rank(P.Sigma_vv_observed{f});
        assert(r_true <= n + 1e-6);
        assert(r_obs  >= p - 1);  % allow 1 rank drop due to tol
        assert_spd(P.Sigma_vv_observed{f});
    end

    r.success=true; r.details='p>>n: ranks behave; observed SPD.';
catch ME
    r = fail_result(r, ME);
end
end

function r = subtest_complex_path()
r = init_result('Complex path & Hermitian');
try
    n=10; p=19; F=6; T=80;
    [Omega, Sigma, Sigma_emp, P] = module7_simulation_improved_complex( ...
        'n_nodes', n, 'n_sensors', p, 'n_freq', F, 'n_samples', T, ...
        'generate_leadfield', true, 'leadfield_type','simple', ...
        'coefficient_complex_fraction', 1.0, 'random_seed', 5);

    % Some source covariances should have non-zero imaginary parts
    has_imag = false; for f=1:F, has_imag = has_imag || any(abs(imag(Sigma{f}(:)))>1e-10); end
    assert(has_imag, 'Expected complex source covariances');

    % Sensor-domain true/obs/emp covariances Hermitian
    for f=1:F
        assert_hermitian(P.Sigma_vv_true{f});
        assert_hermitian(P.Sigma_vv_observed{f});
        assert_hermitian(Sigma_emp{f});
    end

    r.success = true; r.details = 'Complex path produces Hermitian matrices as expected.';
catch ME
    r = fail_result(r, ME);
end
end

function r = subtest_snr_shaping()
r = init_result('SNR shaping');
try
    n=8; p=19; F=3; T=60;
    % 10 dB
    [~, ~, ~, P10] = module7_simulation_improved_complex( ...
        'n_nodes', n, 'n_sensors', p, 'n_freq', F, 'n_samples', T, ...
        'generate_leadfield', true, 'sensor_snr_db', 10, 'random_seed', 9);
    % 20 dB
    [~, ~, ~, P20] = module7_simulation_improved_complex( ...
        'n_nodes', n, 'n_sensors', p, 'n_freq', F, 'n_samples', T, ...
        'generate_leadfield', true, 'sensor_snr_db', 20, 'random_seed', 9);

    % 直接用观测-真值得到当频点噪声
    Sxi10 = P10.Sigma_vv_observed{1} - P10.Sigma_vv_true{1};
    Sxi20 = P20.Sigma_vv_observed{1} - P20.Sigma_vv_true{1};

    [ok10,msg10] = check_trace_snr(P10.Sigma_vv_true{1}, Sxi10, 10^(10/10), 0.25);
    [ok20,msg20] = check_trace_snr(P20.Sigma_vv_true{1}, Sxi20, 10^(20/10), 0.25);
    assert(ok10, msg10); assert(ok20, msg20);

    r.success=true; r.details='SNR shaping around 10/20 dB OK (trace-based).';
catch ME
    r = fail_result(r, ME);
end
end

function r = subtest_param_validation()
r = init_result('Parameter validation');
try
    nerr = 0; ntotal = 0;

    ntotal=ntotal+1; try
        module7_simulation_improved_complex('generate_leadfield', true, 'leadfield_type','nope'); %#ok<NASGU>
    catch, nerr=nerr+1; end

    ntotal=ntotal+1; try
        module7_simulation_improved_complex('generate_leadfield', true, 'head_radius', -0.1); %#ok<NASGU>
    catch, nerr=nerr+1; end

    ntotal=ntotal+1; try
        module7_simulation_improved_complex('generate_leadfield', true, 'conductivities', [0.33 -0.1 0.33]); %#ok<NASGU>
    catch, nerr=nerr+1; end

    ntotal=ntotal+1; try
        module7_simulation_improved_complex('generate_leadfield', true, 'layer_radii', [0.05 0.10]); %#ok<NASGU>
    catch, nerr=nerr+1; end

    ntotal=ntotal+1; try
        module7_simulation_improved_complex('generate_leadfield', true, 'n_sensors', 0); %#ok<NASGU>
    catch, nerr=nerr+1; end

    ntotal=ntotal+1; try
        module7_simulation_improved_complex('n_nodes', 2); %#ok<NASGU>
    catch, nerr=nerr+1; end

    ok = (nerr==ntotal);
    assert(ok);
    r.success=true; r.details=sprintf('Caught %d/%d invalid parameter attempts', nerr, ntotal);
catch ME
    r = fail_result(r, ME);
end
end

function r = subtest_zero_variance_noise()
r = init_result('Zero-variance noise');
try
    n=8; p=19; F=3; T=40;
    [~, ~, Sigma_emp, P] = module7_simulation_improved_complex( ...
        'n_nodes', n, 'n_sensors', p, 'n_freq', F, 'n_samples', T, ...
        'generate_leadfield', true, 'sensor_noise_mode','variance', 'sensor_noise_variance', 0, ...
        'random_seed', 12);

    for f=1:F
        rtrue = rank(P.Sigma_vv_true{f}); robs = rank(P.Sigma_vv_observed{f});
        assert(robs == rtrue);
        assert_hermitian(Sigma_emp{f});
    end

    r.success=true; r.details='Zero-variance noise path executed; Hermitian empirical covariances produced.';
catch ME
    r = fail_result(r, ME);
end
end

function r = subtest_seed_determinism()
r = init_result('Random seed determinism');
try
    args = {'n_nodes', 8, 'n_sensors', 19, 'n_freq', 4, 'n_samples', 50, 'generate_leadfield', true, 'random_seed', 2025};
    [~, ~, Sigma_emp1, P1] = module7_simulation_improved_complex(args{:});
    [~, ~, Sigma_emp2, P2] = module7_simulation_improved_complex(args{:});

    assert(norm(P1.leadfield_matrix - P2.leadfield_matrix,'fro')==0);
    d=0; for f=1:4, d = d + norm(Sigma_emp1{f} - Sigma_emp2{f}, 'fro'); end
    assert(d==0);

    r.success=true; r.details='Deterministic with fixed seed.';
catch ME
    r = fail_result(r, ME);
end
end

%% -------------------------------------------------------------------------
function r = init_result(name)
r = struct('test_name', name, 'success', false, 'details','');
end

function r = fail_result(r, ME)
r.success = false; r.details = ME.message; fprintf('❌ %s: %s\n', r.test_name, ME.message);
end

function assert_hermitian(A, tol)
if nargin<2, tol=1e-8; end
if ~isequal(size(A,1), size(A,2)), error('Matrix is not square.'); end
if norm(A - A', 'fro') > tol*max(1,norm(A,'fro'))
    error('Matrix is not Hermitian within tolerance.');
end
end

function assert_spd(A, tol)
if nargin<2, tol=1e-10; end
E = eig((A+A')/2);
if min(real(E)) < tol
    error('Matrix is not SPD (min eig = %g).', min(real(E)));
end
end

function [ok,msg] = check_trace_snr(S_true, S_noise, snr_lin_target, rel_tol)
if nargin<4, rel_tol = 0.20; end
rt = trace(S_true); rn = trace(S_noise);
ratio = rt/max(rn, eps);
ok = abs(ratio - snr_lin_target)/snr_lin_target <= rel_tol;
msg = sprintf('Trace-SNR check: got %.3f, target %.3f (tol %.0f%%)', ratio, snr_lin_target, rel_tol*100);
end

function summary = summarize_results(tr)
fields = {'legacy','simple19','spherical','pggn','complexpath','snr','paramval','noisezero','seed'};
passed = 0; total = numel(fields);
for i=1:total
    f = fields{i};
    if isfield(tr, f) && isfield(tr.(f),'success') && tr.(f).success
        passed = passed + 1;
    end
end
rate = passed/total;
if rate==1.0
    status='EXCELLENT'; ready='Ready for production use';
elseif rate>=0.8
    status='GOOD'; ready='Ready with minor issues noted';
elseif rate>=0.6
    status='PARTIAL'; ready='Needs attention before production';
else
    status='POOR'; ready='Requires significant fixes';
end
summary = struct('total_tests', total, 'passed_tests', passed, 'success_rate', rate, ...
                 'overall_status', status, 'readiness_assessment', ready, 'timestamp', datestr(now));
end

function print_summary(s)
fprintf('\n=== Test Summary ===\n');
fprintf('  Total tests   : %d\n', s.total_tests);
fprintf('  Passed tests  : %d\n', s.passed_tests);
fprintf('  Success rate  : %.1f%%%%\n', s.success_rate*100);
fprintf('  Overall status: %s\n', s.overall_status);
fprintf('  Readiness     : %s\n', s.readiness_assessment);
end
