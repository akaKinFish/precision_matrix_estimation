function [Omega_cell, Sjj_cell, outs] = solver_jspace_pspline_em(Svv, L, freq, dwi_C, cfg)
%SOLVER_JSPACE_PSPLINE_EM  J-SPACE EM solver (Route A, P-spline).
%
% Flow:
%   1) Build spline basis and penalty matrices
%   2) E-step: posterior second moments
%   3) Whitening
%   4) M-step: FISTA on spline coefficients
%   5) Recolor + inertia update of Sigma_jj
%
% Inputs:
%   Svv   : (Ne x Ne x F) tensor or {F x 1} cell
%   L     : (Ne x N) leadfield
%   freq  : (F x 1) frequency vector
%   dwi_C : (N x N) anatomical connectivity
%   cfg   : config struct
%
% Outputs:
%   Omega_cell : {F x 1} precision matrices
%   Sjj_cell   : {F x 1} source covariance matrices
%   outs       : diagnostics and model info

    if nargin < 5 || isempty(cfg)
        cfg = struct();
    end

    if isempty(Svv) || isempty(L)
        error('solver_jspace_pspline_em:InvalidInput', 'Svv and L are required.');
    end

    % --- Basic sizes ---
    [Ne, N] = size(L);
    if isnumeric(Svv)
        if ndims(Svv) == 2
            if size(Svv, 1) ~= size(Svv, 2)
                error('solver_jspace_pspline_em:SizeMismatch', 'Svv must be square.');
            end
            if size(Svv, 1) ~= Ne
                error('solver_jspace_pspline_em:SizeMismatch', ...
                    'Svv channels do not match leadfield.');
            end
            F = 1;
        elseif ndims(Svv) == 3
            [Ne_s, Ne_s2, F] = size(Svv);
            if Ne_s ~= Ne_s2
                error('solver_jspace_pspline_em:SizeMismatch', 'Svv must be square per frequency.');
            end
            if Ne_s ~= Ne
                error('solver_jspace_pspline_em:SizeMismatch', ...
                    'Svv channels do not match leadfield.');
            end
        else
            error('solver_jspace_pspline_em:InvalidInput', ...
                'Svv must be 2D or 3D numeric tensor.');
        end
    elseif iscell(Svv)
        F = numel(Svv);
        if F < 1
            error('solver_jspace_pspline_em:InvalidInput', 'Svv cell is empty.');
        end
        for t = 1:F
            Svt = Svv{t};
            if ~isequal(size(Svt), [Ne Ne])
                error('solver_jspace_pspline_em:SizeMismatch', ...
                    'Svv{%d} size does not match leadfield.', t);
            end
        end
    else
        error('solver_jspace_pspline_em:InvalidInput', 'Svv must be numeric or cell.');
    end

    if numel(freq) ~= F
        error('solver_jspace_pspline_em:SizeMismatch', ...
            'freq length must match number of frequencies.');
    end

    % --- Build spline basis and penalties ---
    spline_cfg = get_cfg_value(cfg, {'spline'}, struct());
    if ~isfield(spline_cfg, 'K') || isempty(spline_cfg.K)
        spline_cfg.K = 6;
    end
    if ~isfield(spline_cfg, 'degree') || isempty(spline_cfg.degree)
        spline_cfg.degree = 3;
    end
    if ~isfield(spline_cfg, 'diff_order') || isempty(spline_cfg.diff_order)
        spline_cfg.diff_order = 2;
    end

    if get_cfg_value(cfg, {'em', 'verbose'}, false)
        fprintf('  [Init] Building spline basis (K=%d)...\n', spline_cfg.K);
    end

    [B, basis_info] = build_bspline_basis(freq, spline_cfg);
    K = size(B, 2);
    Ddiff = build_diffmat(K, spline_cfg.diff_order);
    G_spline = Ddiff' * Ddiff;

    % --- Prior weights ---
    if nargin < 4 || isempty(dwi_C)
        w_edge = ones(N);
        w_edge(1:N+1:end) = 0;
        W_spatial = zeros(N);
    else
        if ~isequal(size(dwi_C), [N N])
            error('solver_jspace_pspline_em:SizeMismatch', ...
                'dwi_C must be N x N.');
        end
        dwi_cfg = get_cfg_value(cfg, {'dwi'}, struct());
        w_cfg = get_cfg_value(cfg, {'weights'}, struct());
        w_edge = build_edge_weights_from_dwi(dwi_C, dwi_cfg, w_cfg);

        spatial_cfg = get_cfg_value(cfg, {'spatial'}, struct());
        if ~isfield(spatial_cfg, 'enable') || isempty(spatial_cfg.enable)
            spatial_cfg.enable = false;
        end
        W_spatial = build_spatial_weight_W(dwi_C, spatial_cfg);
    end

    % --- Noise covariance ---
    noise_cfg = get_cfg_value(cfg, {'noise'}, struct());
    Sigma_ee_cell = cell(F, 1);
    if isfield(noise_cfg, 'Sigma_ee_cell') && ~isempty(noise_cfg.Sigma_ee_cell)
        Sigma_ee_cell = normalize_cell_matrix(noise_cfg.Sigma_ee_cell, F, 'Sigma_ee_cell');
    else
        sigma2 = get_cfg_value(cfg, {'noise', 'sigma2'}, 1e-6);
        if numel(sigma2) > 1 && numel(sigma2) ~= F
            error('solver_jspace_pspline_em:SizeMismatch', ...
                'noise.sigma2 must be scalar or length F.');
        end
        for t = 1:F
            s2 = sigma2;
            if numel(sigma2) > 1
                s2 = sigma2(t);
            end
            Sigma_ee_cell{t} = s2 * eye(Ne, 'like', L);
        end
    end

    for t = 1:F
        if ~isequal(size(Sigma_ee_cell{t}), [Ne Ne])
            error('solver_jspace_pspline_em:SizeMismatch', ...
                'Sigma_ee_cell{%d} must be Ne x Ne.', t);
        end
    end

    % --- Initialize Sigma_jj ---
    init_cfg = get_cfg_value(cfg, {'init'}, struct());
    if isfield(init_cfg, 'Sjj0') && ~isempty(init_cfg.Sjj0)
        Sigma_jj_cell = normalize_cell_matrix(init_cfg.Sjj0, F, 'Sjj0');
    else
        Sigma_jj_cell = cell(F, 1);
        for t = 1:F
            Sigma_jj_cell{t} = eye(N, 'like', L);
        end
    end
    for t = 1:F
        if ~isequal(size(Sigma_jj_cell{t}), [N N])
            error('solver_jspace_pspline_em:SizeMismatch', ...
                'Sigma_jj_cell{%d} must be N x N.', t);
        end
    end

    % --- Initialize spline coefficients ---
    Theta_cell = cell(K, 1);
    Theta_cell{1} = eye(N, 'like', L);
    for k = 2:K
        Theta_cell{k} = zeros(N, N, 'like', L);
    end

    % --- EM config ---
    em_max_iter = get_cfg_value(cfg, {'em', 'max_iter'}, 20);
    em_m_samples = get_cfg_value(cfg, {'em', 'm_samples'}, 1);
    em_update_rate = get_cfg_value(cfg, {'em', 'update_rate'}, 0.3);
    em_verbose = get_cfg_value(cfg, {'em', 'verbose'}, false);
    em_tol = get_cfg_value(cfg, {'em', 'tol'}, 1e-4);

    if ~isfinite(em_update_rate)
        em_update_rate = 0.3;
    end
    em_update_rate = min(max(em_update_rate, 0), 1);

    % --- M-step config ---
    mstep_cfg = get_cfg_value(cfg, {'mstep'}, struct());
    if ~isfield(mstep_cfg, 'm_samples') || isempty(mstep_cfg.m_samples)
        mstep_cfg.m_samples = em_m_samples;
    end

    whiten_cfg = get_cfg_value(cfg, {'whiten'}, struct());

    outs.loglik = zeros(em_max_iter, 1);
    outs.obj_mstep = zeros(em_max_iter, 1);
    outs.minEig = zeros(em_max_iter, 1);
    outs.sparsity = zeros(em_max_iter, 1);

    if em_verbose
        fprintf('  [EM] Starting loop (MaxIter=%d)...\n', em_max_iter);
    end

    Omega_cell = cell(F, 1);

    for iter = 1:em_max_iter
        t_start = tic;

        % --- E-step ---
        [Psi_cell, stats_e] = estep_posterior_psi( ...
            Svv, L, Sigma_jj_cell, Sigma_ee_cell, em_m_samples);
        outs.loglik(iter) = stats_e.loglik;

        % --- Whitening ---
        [D_cell, Psi_tilde_cell] = whiten_psi(Psi_cell, whiten_cfg);

        % --- M-step ---
        Theta_init = Theta_cell;
        [Theta_new, Gamma_cell, out_m] = mstep_pspline_precision_fista( ...
            Psi_tilde_cell, B, G_spline, w_edge, W_spatial, Theta_init, mstep_cfg);
        outs.obj_mstep(iter) = out_m.obj;
        outs.minEig(iter) = out_m.minEig_all;

        % --- Recolor ---
        Omega_cell = recolor_precision(Gamma_cell, D_cell);

        % --- Update Sigma_jj (inertia) ---
        Sjj_new_cell = cell(F, 1);
        for t = 1:F
            [S_new, ~, is_spd] = spd_inv_and_logdet(Omega_cell{t});
            if ~is_spd || any(~isfinite(S_new(:)))
                if em_verbose
                    fprintf('  [EM] Warning: Omega not SPD at f=%d. Keeping previous Sjj.\n', t);
                end
                S_new = Sigma_jj_cell{t};
            end
            S_new = 0.5 * (S_new + S_new');
            S_new(1:N+1:end) = real(diag(S_new));
            Sjj_new_cell{t} = S_new;
        end

        rho = em_update_rate;
        for t = 1:F
            S_old = Sigma_jj_cell{t};
            S_new = Sjj_new_cell{t};
            S_updated = (1 - rho) * S_old + rho * S_new;
            S_updated = 0.5 * (S_updated + S_updated');
            S_updated(1:N+1:end) = real(diag(S_updated));
            Sigma_jj_cell{t} = S_updated;
        end

        Theta_cell = Theta_new;

        % --- Diagnostics ---
        [~, ~, current_density] = simple_pcor_density(Omega_cell, 0.01);
        outs.sparsity(iter) = current_density;

        t_iter = toc(t_start);
        if em_verbose
            fprintf('  Iter %2d | LL: %.2e | Obj: %.2e | Dens: %.1f%% | minEig: %.1e | T: %.1fs\n', ...
                iter, stats_e.loglik, out_m.obj, current_density * 100, out_m.minEig_all, t_iter);
        end

        % Convergence check on log-likelihood
        if iter > 1
            denom = max(abs(outs.loglik(iter-1)), eps);
            diff_ll = abs(outs.loglik(iter) - outs.loglik(iter-1)) / denom;
            if diff_ll < em_tol && iter > 5
                if em_verbose
                    fprintf('  [EM] Converged (LL change < %.1e).\n', em_tol);
                end
                break;
            end
        end
    end

    outs.loglik = outs.loglik(1:iter);
    outs.obj_mstep = outs.obj_mstep(1:iter);
    outs.minEig = outs.minEig(1:iter);
    outs.sparsity = outs.sparsity(1:iter);

    outs.model = struct();
    outs.model.B = B;
    outs.model.Theta = Theta_cell;
    outs.model.w_edge = w_edge;
    outs.model.W_spatial = W_spatial;
    outs.basis_info = basis_info;
    outs.final_iter = iter;

    Sjj_cell = Sigma_jj_cell;
end

function val = get_cfg_value(cfg, fields, default_val)
    val = default_val;
    if ~isstruct(cfg)
        return;
    end
    if ischar(fields) || (isstring(fields) && isscalar(fields))
        fields = {char(fields)};
    end
    cur = cfg;
    for i = 1:numel(fields)
        if ~isstruct(cur) || ~isfield(cur, fields{i})
            return;
        end
        cur = cur.(fields{i});
    end
    if ~isempty(cur)
        val = cur;
    end
end

function C = normalize_cell_matrix(x, F, name)
    if nargin < 3
        name = 'input';
    end
    if iscell(x)
        if numel(x) == 1 && F > 1
            C = repmat(x, F, 1);
        elseif numel(x) == F
            C = x(:);
        else
            error('solver_jspace_pspline_em:SizeMismatch', ...
                '%s must have length 1 or F.', name);
        end
    elseif ndims(x) == 3
        if size(x, 3) ~= F
            error('solver_jspace_pspline_em:SizeMismatch', ...
                '%s must have F slices.', name);
        end
        C = cell(F, 1);
        for t = 1:F
            C{t} = x(:, :, t);
        end
    else
        C = repmat({x}, F, 1);
    end
end

function [pmax, p99, dens] = simple_pcor_density(Omega_cell, thr)
    % Estimate average connection density from partial correlations.
    F = numel(Omega_cell);
    N = size(Omega_cell{1}, 1);
    cnt = 0;
    tot = 0;
    vals_all = [];

    check_idx = unique(round(linspace(1, F, min(3, F))));
    for t = check_idx
        Om = Omega_cell{t};
        d = sqrt(max(real(diag(Om)), eps));
        P = -Om ./ (d * d');
        P(1:N+1:end) = 0;

        mask = triu(true(N), 1);
        vals = abs(P(mask));
        cnt = cnt + sum(vals > thr);
        tot = tot + numel(vals);
        vals_all = [vals_all; vals(:)];
    end

    dens = cnt / (tot + eps);
    if isempty(vals_all)
        pmax = 0;
        p99 = 0;
    else
        pmax = max(vals_all);
        vals_sorted = sort(vals_all);
        idx = max(1, ceil(0.99 * numel(vals_sorted)));
        p99 = vals_sorted(idx);
    end
end
