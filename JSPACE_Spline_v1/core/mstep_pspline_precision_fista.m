function [Theta_cell, Gamma_cell, out] = mstep_pspline_precision_fista( ...
    Psi_tilde_cell, B, G_spline, w_edge, W_spatial, Theta_init, mcfg)
%MSTEP_PSPLINE_PRECISION_FISTA  FISTA optimizer for spline precision model.
%
% Objective:
%   Data term + group-lasso (edges) + P-spline + spatial L2
%
% Inputs:
%   Psi_tilde_cell : {F x 1} whitened posterior second moments
%   B              : (F x K) spline basis matrix
%   G_spline       : (K x K) quadratic form matrix (D' * D)
%   w_edge         : (N x N) edge weights for group-lasso
%   W_spatial      : (N x N) spatial weights
%   Theta_init     : {K x 1} initial spline coefficients
%   mcfg           : config struct (lambda1, lambda_ps, lambda2, m_samples, fista, spd, diag_policy)
%
% Outputs:
%   Theta_cell : {K x 1} optimized spline coefficients
%   Gamma_cell : {F x 1} precision matrices at frequencies
%   out        : struct with objective history and diagnostics

    if nargin < 7 || isempty(mcfg)
        mcfg = struct();
    end

    if ~iscell(Psi_tilde_cell)
        if ndims(Psi_tilde_cell) == 3
            F_tmp = size(Psi_tilde_cell, 3);
            tmp_cell = cell(F_tmp, 1);
            for t = 1:F_tmp
                tmp_cell{t} = Psi_tilde_cell(:, :, t);
            end
            Psi_tilde_cell = tmp_cell;
        else
            Psi_tilde_cell = {Psi_tilde_cell};
        end
    end

    F = numel(Psi_tilde_cell);
    N = size(Psi_tilde_cell{1}, 1);
    K = size(B, 2);

    if size(B, 1) ~= F
        error('mstep_pspline_precision_fista:SizeMismatch', ...
            'B must have F rows matching Psi_tilde_cell.');
    end
    if size(G_spline, 1) ~= K || size(G_spline, 2) ~= K
        error('mstep_pspline_precision_fista:SizeMismatch', ...
            'G_spline must be K x K.');
    end

    % Normalize Theta_init to cell format
    if ~iscell(Theta_init)
        if ndims(Theta_init) == 3
            Theta_cell = cell(K, 1);
            for k = 1:K
                Theta_cell{k} = Theta_init(:, :, k);
            end
            Theta_init = Theta_cell;
        else
            Theta_init = repmat({Theta_init}, K, 1);
        end
    end
    if numel(Theta_init) ~= K
        error('mstep_pspline_precision_fista:SizeMismatch', ...
            'Theta_init must have K elements.');
    end

    if isempty(w_edge)
        w_edge = ones(N);
        w_edge(1:N+1:end) = 0;
    end
    if isempty(W_spatial)
        W_spatial = zeros(N);
    end

    % Config defaults
    lam1 = get_field_default(mcfg, 'lambda1', 0);
    lam_ps = get_field_default(mcfg, 'lambda_ps', 0);
    lam2 = get_field_default(mcfg, 'lambda2', 0);
    m_samps = get_field_default(mcfg, 'm_samples', 1);
    if numel(m_samps) > 1 && numel(m_samps) ~= F
        error('mstep_pspline_precision_fista:SizeMismatch', ...
            'm_samples must be scalar or length F.');
    end

    fista = get_field_default(mcfg, 'fista', struct());
    fista.max_iter = get_field_default(fista, 'max_iter', 100);
    fista.alpha0 = get_field_default(fista, 'alpha0', 1e-2);
    fista.backtracking_beta = get_field_default(fista, 'backtracking_beta', 0.5);
    fista.max_backtracking = get_field_default(fista, 'max_backtracking', 20);
    fista.tol = get_field_default(fista, 'tol', 1e-4);
    fista.verbose = get_field_default(fista, 'verbose', false);

    spd_cfg = get_field_default(mcfg, 'spd', struct());
    spd_cfg.eps_pd = get_field_default(spd_cfg, 'eps_pd', 1e-6);
    diag_policy = get_field_default(mcfg, 'diag_policy', 'theta0_only');
    if isstruct(diag_policy)
        if isfield(diag_policy, 'mode') && ~isempty(diag_policy.mode)
            diag_policy = diag_policy.mode;
        else
            diag_policy = 'theta0_only';
        end
    end
    if isstring(diag_policy) && isscalar(diag_policy)
        diag_policy = char(diag_policy);
    end

    % Keep config updated for helper calls
    mcfg.lambda1 = lam1;
    mcfg.lambda_ps = lam_ps;
    mcfg.lambda2 = lam2;
    mcfg.m_samples = m_samps;
    mcfg.fista = fista;
    mcfg.spd = spd_cfg;
    mcfg.diag_policy = diag_policy;

    % Initialize states
    Theta_old = project_to_hermitian(Theta_init);
    Theta_old = enforce_diag_policy(Theta_old, diag_policy);
    Theta_old = spd_shift_on_Theta0(Theta_old, B, spd_cfg.eps_pd);
    Y = Theta_old;
    t_fista = 1;
    alpha = fista.alpha0;

    [obj_old, ~, ok_init] = compute_objective_routeA( ...
        Theta_old, Psi_tilde_cell, B, G_spline, w_edge, W_spatial, mcfg);
    if ~ok_init
        if fista.verbose
            fprintf('    Initial Theta not SPD. Applying shift.\n');
        end
        Theta_old = spd_shift_on_Theta0(Theta_old, B, spd_cfg.eps_pd);
        Y = Theta_old;
        obj_old = compute_objective_routeA( ...
            Theta_old, Psi_tilde_cell, B, G_spline, w_edge, W_spatial, mcfg);
    end

    best_obj = obj_old;
    history_obj = zeros(fista.max_iter, 1);

    for iter = 1:fista.max_iter
        % --- A: Gradient at Y ---
        [Gamma_Y, invGamma_Y, ~, ~, ok_spd_Y] = assemble_Gamma_from_Theta( ...
            Y, B, spd_cfg.eps_pd);

        if ~ok_spd_Y
            Y = spd_shift_on_Theta0(Y, B, spd_cfg.eps_pd);
            [Gamma_Y, invGamma_Y, ~, ~, ok_spd_Y] = assemble_Gamma_from_Theta( ...
                Y, B, spd_cfg.eps_pd);
            t_fista = 1;
        end

        if ~ok_spd_Y
            Y = Theta_old;
            [Gamma_Y, invGamma_Y, ~, ~, ok_spd_Y] = assemble_Gamma_from_Theta( ...
                Y, B, spd_cfg.eps_pd);
            t_fista = 1;
        end
        if ~ok_spd_Y
            warning('mstep_pspline_precision_fista:NonSPD', ...
                'Failed to obtain SPD Gamma at iter %d. Stopping.', iter);
            break;
        end

        Grad = cell(K, 1);
        for k = 1:K
            Grad{k} = zeros(N, N, 'like', Psi_tilde_cell{1});
        end

        % Data term gradient
        for t = 1:F
            m_t = m_samps;
            if numel(m_samps) > 1
                m_t = m_samps(t);
            end
            Diff = Psi_tilde_cell{t} - invGamma_Y{t};
            for k = 1:K
                w_tk = B(t, k);
                if abs(w_tk) > 1e-12
                    Grad{k} = Grad{k} + (m_t * w_tk) * Diff;
                end
            end
        end

        % P-spline gradient
        if lam_ps > 0
            for k = 1:K
                sum_G_Theta = zeros(N, N, 'like', Psi_tilde_cell{1});
                for l = 1:K
                    gkl = G_spline(k, l);
                    if gkl ~= 0
                        sum_G_Theta = sum_G_Theta + gkl * Y{l};
                    end
                end
                Grad{k} = Grad{k} + (2 * lam_ps) * sum_G_Theta;
            end
        end

        % Spatial gradient
        if lam2 > 0
            W2 = W_spatial .^ 2;
            for k = 1:K
                Grad{k} = Grad{k} + (2 * lam2) * (W2 .* Y{k});
            end
        end

        % Enforce diagonal policy on gradients
        Grad = enforce_diag_policy_grad(Grad, diag_policy);

        % Symmetrize gradients for stability
        for k = 1:K
            Gk = Grad{k};
            Gk = 0.5 * (Gk + Gk');
            Gk(1:N+1:end) = real(diag(Gk));
            Grad{k} = Gk;
        end

        % --- B: Backtracking ---
        alpha_curr = alpha;
        Theta_next = Y;
        has_descended = false;

        for bt = 1:fista.max_backtracking
            Theta_try = cell(K, 1);
            for k = 1:K
                Theta_try{k} = Y{k} - alpha_curr * Grad{k};
            end

            Theta_try = prox_group_lasso_edges(Theta_try, lam1, alpha_curr, w_edge);
            Theta_try = project_to_hermitian(Theta_try);
            Theta_try = enforce_diag_policy(Theta_try, diag_policy);
            Theta_try = spd_shift_on_Theta0(Theta_try, B, spd_cfg.eps_pd);

            [obj_try, ~, ok_try] = compute_objective_routeA( ...
                Theta_try, Psi_tilde_cell, B, G_spline, w_edge, W_spatial, mcfg);

            if ok_try && (obj_try <= obj_old + 1e-10)
                Theta_next = Theta_try;
                best_obj = obj_try;
                has_descended = true;
                break;
            end

            alpha_curr = alpha_curr * fista.backtracking_beta;
        end

        if ~has_descended
            if fista.verbose
                fprintf('    Backtracking failed at iter %d. Stopping.\n', iter);
            end
            Theta_next = Theta_old;
            break;
        end

        alpha = alpha_curr;

        % --- C: FISTA momentum update ---
        t_fista_next = (1 + sqrt(1 + 4 * t_fista^2)) / 2;
        momentum_factor = (t_fista - 1) / t_fista_next;

        Y_next = cell(K, 1);
        for k = 1:K
            Diff_Theta = Theta_next{k} - Theta_old{k};
            Y_next{k} = Theta_next{k} + momentum_factor * Diff_Theta;
        end

        Theta_old = Theta_next;
        Y = Y_next;
        t_fista = t_fista_next;
        obj_old = best_obj;
        history_obj(iter) = best_obj;

        % --- D: Convergence check ---
        if iter > 1
            prev = history_obj(iter - 1);
            if isfinite(prev) && prev ~= 0
                rel_change = abs(history_obj(iter) - prev) / (abs(prev) + eps);
                if rel_change < fista.tol
                    if fista.verbose
                        fprintf('    FISTA converged at iter %d (tol=%.1e).\n', ...
                            iter, rel_change);
                    end
                    break;
                end
            end
        end
    end

    Theta_cell = Theta_old;
    [Gamma_cell, ~, ~, minEig_all, ~] = assemble_Gamma_from_Theta( ...
        Theta_cell, B, spd_cfg.eps_pd);

    out = struct();
    out.obj = best_obj;
    out.minEig_all = minEig_all;
    out.history_obj = history_obj(1:iter);
    out.iterations = iter;
    out.alpha = alpha;
end

function Grad = enforce_diag_policy_grad(Grad, policy)
    mode = 'theta0_only';
    if isstruct(policy) && isfield(policy, 'mode') && ~isempty(policy.mode)
        mode = policy.mode;
    elseif ischar(policy) || (isstring(policy) && isscalar(policy))
        mode = char(policy);
    end

    if strcmpi(mode, 'theta0_only')
        N = size(Grad{1}, 1);
        K = numel(Grad);
        for k = 2:K
            Grad{k}(1:N+1:end) = 0;
        end
    end
end

function val = get_field_default(s, field_name, default_val)
    if isstruct(s) && isfield(s, field_name) && ~isempty(s.(field_name))
        val = s.(field_name);
    else
        val = default_val;
    end
end
