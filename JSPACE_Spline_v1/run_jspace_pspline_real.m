function [Omega_est, Sjj_est, outs] = run_jspace_pspline_real( ...
    Svv_CrossM, L, freq, dwi_C, cfg)
%RUN_JSPACE_PSPLINE_REAL  Top-level entry for J-SPACE (P-spline, Route A).
%
% Responsibilities:
%   1) Add subfolders to MATLAB path
%   2) Validate inputs
%   3) Merge default config with user config
%   4) Sanitize input cross-spectrum
%   5) Call EM solver
%
% Usage:
%   [Omega, Sjj, outs] = run_jspace_pspline_real(Svv, L, freq, dwi_C, cfg);
%
% Inputs:
%   Svv_CrossM : (Ne x Ne x F) tensor or {F x 1} cell
%   L          : (Ne x N) leadfield
%   freq       : (F x 1) frequency vector
%   dwi_C      : (N x N) anatomical connectivity
%   cfg        : optional config struct
%
% Outputs:
%   Omega_est  : {F x 1} precision matrices
%   Sjj_est    : {F x 1} covariance matrices
%   outs       : diagnostics and model info

    % --- 0. Path management ---
    baseDir = fileparts(mfilename('fullpath'));
    if exist(fullfile(baseDir, 'utils'), 'dir')
        addpath(fullfile(baseDir, 'utils'));
    end
    if exist(fullfile(baseDir, 'core'), 'dir')
        addpath(fullfile(baseDir, 'core'));
    end
    if exist(fullfile(baseDir, 'helpers'), 'dir')
        addpath(fullfile(baseDir, 'helpers'));
    end

    if nargin < 5
        cfg = [];
    end

    % --- 1. Input checks ---
    if isempty(Svv_CrossM) || isempty(L) || isempty(freq)
        error('run_jspace_pspline_real:InvalidInput', ...
            'Svv, L, and freq must be provided.');
    end

    [Ne_L, N] = size(L);

    if iscell(Svv_CrossM)
        F = numel(Svv_CrossM);
        if F < 1
            error('run_jspace_pspline_real:InvalidInput', 'Svv cell is empty.');
        end
        for t = 1:F
            Svt = Svv_CrossM{t};
            if ~isequal(size(Svt), [Ne_L Ne_L])
                error('run_jspace_pspline_real:SizeMismatch', ...
                    'Svv{%d} size does not match leadfield.', t);
            end
        end
    else
        if ndims(Svv_CrossM) == 2
            [Ne1, Ne2] = size(Svv_CrossM);
            if Ne1 ~= Ne2
                error('run_jspace_pspline_real:SizeMismatch', ...
                    'Svv must be square.');
            end
            if Ne1 ~= Ne_L
                error('run_jspace_pspline_real:SizeMismatch', ...
                    'Svv channels do not match leadfield.');
            end
            F = 1;
        else
            [Ne1, Ne2, F] = size(Svv_CrossM);
            if Ne1 ~= Ne2
                error('run_jspace_pspline_real:SizeMismatch', ...
                    'Svv must be square per frequency.');
            end
            if Ne1 ~= Ne_L
                error('run_jspace_pspline_real:SizeMismatch', ...
                    'Svv channels do not match leadfield.');
            end
        end
    end

    if numel(freq) ~= F
        error('run_jspace_pspline_real:SizeMismatch', ...
            'freq length must match number of frequencies.');
    end

    if nargin < 4 || isempty(dwi_C)
        error('run_jspace_pspline_real:InvalidInput', ...
            'dwi_C must be provided as N x N.');
    end
    if ~isequal(size(dwi_C), [N N])
        error('run_jspace_pspline_real:SizeMismatch', ...
            'dwi_C must be N x N.');
    end

    % --- 2. Config merge ---
    default_cfg = get_default_cfg_routeA();
    cfg = merge_cfg(default_cfg, cfg);

    % --- 3. Sanitize data ---
    fprintf('=== J-SPACE (P-Spline Route A) Started ===\n');
    fprintf('  Data: Ne=%d, N=%d, F=%d\n', Ne_L, N, F);
    fprintf('  Spline: K=%d, Lambda1=%.2e, Lambda_PS=%.2e\n', ...
        cfg.spline.K, cfg.mstep.lambda1, cfg.mstep.lambda_ps);

    if iscell(Svv_CrossM)
        Svv = cell(F, 1);
        for t = 1:F
            Svv{t} = sanitize_cross_spectrum(Svv_CrossM{t}, cfg.svv);
        end
    else
        Svv = sanitize_cross_spectrum(Svv_CrossM, cfg.svv);
    end

    % --- 4. Call solver ---
    t_total = tic;
    [Omega_est, Sjj_est, outs] = solver_jspace_pspline_em(Svv, L, freq, dwi_C, cfg);
    outs.total_time = toc(t_total);

    fprintf('=== Finished in %.1f seconds. ===\n', outs.total_time);
end

function cfg = get_default_cfg_routeA()
    cfg = struct();

    % Input data cleaning
    cfg.svv.jitter = 1e-10;

    % Noise model
    cfg.noise.sigma2 = 1e-3;

    % EM control
    cfg.em.max_iter = 15;
    cfg.em.update_rate = 0.3;
    cfg.em.m_samples = 1000;
    cfg.em.verbose = true;
    cfg.em.tol = 1e-4;

    % Initialization
    cfg.init.Sjj0 = [];

    % Spline settings
    cfg.spline.K = 8;
    cfg.spline.degree = 3;
    cfg.spline.diff_order = 2;

    % DWI -> weights
    cfg.dwi.alpha = 1.0;
    cfg.dwi.clip = [0.2, 5.0];
    cfg.dwi.eps0 = 1e-5;
    cfg.weights.normalize_weight_median = true;

    % Spatial prior
    cfg.spatial.enable = false;
    cfg.spatial.clip = [0.1, 10];

    % M-step (FISTA)
    cfg.mstep.lambda1 = 0.05;
    cfg.mstep.lambda_ps = 0.1;
    cfg.mstep.lambda2 = 0.0;

    cfg.mstep.fista.max_iter = 50;
    cfg.mstep.fista.tol = 1e-4;
    cfg.mstep.fista.alpha0 = 1e-2;
    cfg.mstep.fista.backtracking_beta = 0.5;
    cfg.mstep.fista.max_backtracking = 20;
    cfg.mstep.fista.verbose = false;

    cfg.mstep.spd.eps_pd = 1e-6;
    cfg.mstep.diag_policy.mode = 'theta0_only';

    % Whitening
    cfg.whiten.eps_diag = 1e-12;
end

function final_cfg = merge_cfg(default_cfg, user_cfg)
    final_cfg = default_cfg;
    if isempty(user_cfg)
        return;
    end
    if ~isstruct(user_cfg)
        error('run_jspace_pspline_real:InvalidInput', 'cfg must be a struct.');
    end
    final_cfg = merge_structs_(final_cfg, user_cfg);
end

function out = merge_structs_(base, override)
    out = base;
    fns = fieldnames(override);
    for i = 1:numel(fns)
        key = fns{i};
        val = override.(key);
        if isstruct(val) && isfield(base, key) && isstruct(base.(key))
            out.(key) = merge_structs_(base.(key), val);
        else
            if ~isempty(val)
                out.(key) = val;
            end
        end
    end
end
