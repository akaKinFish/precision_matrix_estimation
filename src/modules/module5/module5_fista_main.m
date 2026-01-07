function [Gamma_cells, results] = module5_fista_main(input_data, params)
% MODULE5_FISTA_MAIN - FISTA optimizer for M-step (cell-valued Hermitian matrices).
%
% Signature is kept compatible with module5_proximal_main:
%   [Gamma_cells, results] = module5_fista_main(input_data, params)
% where input_data contains:
%   .whitened_covariances (Sigmas) : {F x 1}
%   .smoothing_kernel (Kernel)     : (F x F)
%   .weight_matrix (W)             : (p x p)
%   .precision_matrices            : optional warm start {F x 1}
%   .active_mask                   : optional {F x 1} logical masks
%
% This solver separates the smooth part (lambda2=0) and uses
% module_proximal_operator_fista for the nonsmooth L1-off prox.
%
% Typical call from solver_jspace_adaptive:
%   m5_in.whitened_covariances = Sjj_tilde;
%   m5_in.smoothing_kernel     = K_freq;
%   m5_in.weight_matrix        = W_gamma;
%   m5_in.active_mask          = module3_active_set(Sjj_tilde, ...);
%   m5_in.precision_matrices   = Gamma_warm_start;
%   [Gamma_hat, res] = module5_fista_main(m5_in, m5_p);
%
% Outputs (results) mirror PGD fields where possible:
%   .objective_history : f + lambda2 * g per iteration
%   .smooth_history    : f values
%   .l1_history        : lambda2 * g values
%   .final_iter, .final_obj, .final_alpha

    % ----------- Inputs -----------
    Sigmas = input_data.whitened_covariances;
    Kernel = input_data.smoothing_kernel;
    W      = input_data.weight_matrix;
    F = numel(Sigmas);
    p = size(Sigmas{1}, 1);

    if isfield(input_data, 'precision_matrices') && ~isempty(input_data.precision_matrices)
        Gamma_curr = input_data.precision_matrices;
    else
        Gamma_curr = cell(F, 1);
        for f = 1:F, Gamma_curr{f} = eye(p, 'like', Sigmas{1}); end
    end

    if isfield(input_data, 'active_mask')
        active_mask = input_data.active_mask;
    else
        active_mask = [];
    end

    if nargin < 2, params = struct(); end
    params = set_default(params, 'max_iter', 150);
    params = set_default(params, 'tol', 1e-5);
    params = set_default(params, 'verbose', false);
    params = set_default(params, 'lambda1', 0);
    params = set_default(params, 'lambda2', 0);
    params = set_default(params, 'lambda3', 0);
    params = set_default(params, 'weight_mode', 'matrix');
    params = set_default(params, 'min_eig', 1e-8);
    params = set_default(params, 'alpha0', 0.1); % spectral guess: L0 = 1/alpha0
    params = set_default(params, 'backtracking_factor', 2.0);
    params = set_default(params, 'penalize_diagonal', false);
    params = set_default(params, 'enforce_unit_diagonal', false);

    % Smooth-only params (lambda2 stripped)
    params_s = params;
    params_s.lambda2 = 0;

    % Ensure warm start is Hermitian SPD
    for f = 1:F
        Gamma_curr{f} = utils_math.make_hermitian(Gamma_curr{f});
        [Gamma_curr{f}, ~] = utils_math.project_spd(Gamma_curr{f}, params.min_eig);
    end
    Gamma_prev = Gamma_curr;
    Y_curr = Gamma_curr;
    t_curr = 1;

    % Initial step (Lipschitz estimate)
    if isfield(params, 'alpha0') && params.alpha0 > 0
        L_curr = 1 / params.alpha0;
    else
        L_curr = 10; % conservative default
    end
    alpha_curr = 1 / L_curr;

    % Initial objective
    f_curr = module_objective.compute(Gamma_curr, Sigmas, Kernel, W, params_s);
    g_curr = local_l1_off(Gamma_curr, params);
    obj_curr = f_curr + params.lambda2 * g_curr;

    obj_history = obj_curr;
    smooth_history = f_curr;
    l1_history = params.lambda2 * g_curr;
    alpha_history = alpha_curr;

    if params.verbose
        fprintf('%4s | %12s | %12s | %8s\n', 'iter', 'obj', 'rel_diff', 'alpha');
    end

    % ----------- Main loop -----------
    for iter = 1:params.max_iter
        % Gradient at extrapolated point
        Grad_Y = module_gradient.compute(Y_curr, Sigmas, Kernel, W, params_s);
        f_Y = module_objective.compute(Y_curr, Sigmas, Kernel, W, params_s);

        % Backtracking on f
        accepted = false;
        bt_count = 0;
        while ~accepted
            alpha_bt = 1 / L_curr;
            Gamma_candidate = cell(F, 1);
            for f = 1:F
                mask_f = [];
                if ~isempty(active_mask), mask_f = active_mask{f}; end
                Z = Y_curr{f} - alpha_bt * Grad_Y{f};
                [Gamma_candidate{f}, ~] = module_proximal_operator_fista.compute( ...
                    Z, params.lambda2 * alpha_bt, mask_f, params);
            end

            Delta = cell_diff(Gamma_candidate, Y_curr);
            diff_norm_sq = cell_norm_sq(Delta);
            lin_term = cell_inner(Grad_Y, Delta);

            f_new = module_objective.compute(Gamma_candidate, Sigmas, Kernel, W, params_s);
            majorant = f_Y + lin_term + (L_curr/2) * diff_norm_sq;

            if f_new <= majorant * (1 + 1e-12)
                accepted = true;
                Gamma_next = Gamma_candidate;
                f_curr = f_new;
                alpha_curr = alpha_bt;
            else
                L_curr = L_curr * params.backtracking_factor;
                bt_count = bt_count + 1;
                if bt_count > 30
                    % fallback to tiny step
                    L_curr = L_curr * 10;
                end
            end
        end

        g_curr = local_l1_off(Gamma_next, params);
        obj_new = f_curr + params.lambda2 * g_curr;

        % Nesterov update
        t_next = (1 + sqrt(1 + 4 * t_curr^2)) / 2;
        Y_next = cell_axpy(Gamma_next, cell_diff(Gamma_next, Gamma_curr), (t_curr - 1) / t_next);

        % Metrics
        rel_diff = sqrt(cell_norm_sq(cell_diff(Gamma_next, Gamma_curr)) / (cell_norm_sq(Gamma_curr) + 1e-12));
        obj_history(end+1) = obj_new; %#ok<*AGROW>
        smooth_history(end+1) = f_curr;
        l1_history(end+1) = params.lambda2 * g_curr;
        alpha_history(end+1) = alpha_curr;

        if params.verbose && (iter == 1 || mod(iter, 10) == 0)
            fprintf('%4d | %12.4e | %12.3e | %8.1e\n', iter, obj_new, rel_diff, alpha_curr);
        end

        % Convergence
        if rel_diff < params.tol
            Gamma_curr = Gamma_next;
            obj_curr = obj_new;
            break;
        end

        % Prepare next iter
        Gamma_prev = Gamma_curr;
        Gamma_curr = Gamma_next;
        obj_curr = obj_new;
        t_curr = t_next;
        Y_curr = Y_next;
        % keep L_curr as-is (using successful step)
    end

    Gamma_cells = Gamma_curr;

    results = struct();
    results.objective_history = obj_history;
    results.smooth_history = smooth_history;
    results.l1_history = l1_history;
    results.step_history = alpha_history;
    results.final_iter = numel(obj_history);
    results.final_obj = obj_history(end);
    results.final_alpha = alpha_curr;
end

% ----------- Local helpers -----------
function params = set_default(params, field, val)
    if ~isfield(params, field), params.(field) = val; end
end

function g = local_l1_off(Gamma_cells, params)
    % Sum of abs off-diagonal entries across frequencies (no lambda2 scale)
    F = numel(Gamma_cells);
    p = size(Gamma_cells{1}, 1);
    penalize_diag = isfield(params, 'penalize_diagonal') && params.penalize_diagonal;
    g = 0;
    for f = 1:F
        G = Gamma_cells{f};
        if ~penalize_diag
            G(1:p+1:end) = 0;
        end
        g = g + sum(abs(G(:)));
    end
end

function Delta = cell_diff(A, B)
    F = numel(A);
    Delta = cell(F,1);
    for f = 1:F
        Delta{f} = A{f} - B{f};
    end
end

function val = cell_norm_sq(CellMats)
    val = 0;
    F = numel(CellMats);
    for f = 1:F
        val = val + norm(CellMats{f}, 'fro')^2;
    end
end

function val = cell_inner(A, B)
    val = 0;
    F = numel(A);
    for f = 1:F
        val = val + real(sum(sum(conj(A{f}) .* B{f})));
    end
end

function out = cell_axpy(Base, Direction, scale)
    F = numel(Base);
    out = cell(F, 1);
    for f = 1:F
        out{f} = Base{f} + scale * Direction{f};
    end
end
