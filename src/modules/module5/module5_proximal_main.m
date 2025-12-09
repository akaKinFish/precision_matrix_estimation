function [Gamma_cells, results] = module5_proximal_main(input_data, params)
% MODULE5_PROXIMAL_MAIN_ROBUST
% Proximal Gradient Descent with backtracking line search for SPD stability.

    % --- 1. Init ---
    Sigmas = input_data.whitened_covariances;
    Kernel = input_data.smoothing_kernel;
    W      = input_data.weight_matrix;
    F = numel(Sigmas);
    p = size(Sigmas{1}, 1);

    if isfield(input_data, 'precision_matrices') && ~isempty(input_data.precision_matrices)
        Gamma_curr = input_data.precision_matrices;
    else
        Gamma_curr = cell(F, 1);
        for f = 1:F, Gamma_curr{f} = eye(p); end
    end

    if isfield(input_data, 'active_mask'), active_mask = input_data.active_mask; else, active_mask = []; end

    if nargin < 2, params = struct(); end
    % Defaults
    if ~isfield(params, 'lambda1'), params.lambda1 = 0; end
    if ~isfield(params, 'lambda2'), params.lambda2 = 0; end
    if ~isfield(params, 'lambda3'), params.lambda3 = 0; end
    if ~isfield(params, 'max_iter'), params.max_iter = 100; end
    if ~isfield(params, 'tol'), params.tol = 1e-5; end
    if ~isfield(params, 'alpha0'), params.alpha0 = 1.0; end
    if ~isfield(params, 'verbose'), params.verbose = false; end

    obj_history = [];
    alpha = params.alpha0;

    % Initial objective
    [f_curr, ~] = module_objective.compute(Gamma_curr, Sigmas, Kernel, W, params);

    for iter = 1:params.max_iter

        % A. Gradient
        try
            Grad_curr = module_gradient.compute(Gamma_curr, Sigmas, Kernel, W, params);
        catch
            if params.verbose, fprintf('[PGD] Gradient failed. Stop.\n'); end
            break;
        end

        % B. Backtracking line search
        beta = 0.5;
        max_backtrack = 20;
        success = false;

        for bt = 0:max_backtrack
            Gamma_new = cell(F,1);
            is_valid_step = true;

            for f = 1:F
                mask_f = [];
                if ~isempty(active_mask), mask_f = active_mask{f}; end
                try
                    Gamma_new{f} = module_proximal_operator.compute(...
                        Gamma_curr{f}, Grad_curr{f}, alpha, params.lambda2, mask_f, params);
                    if rcond(Gamma_new{f}) < 1e-15
                        is_valid_step = false; break;
                    end
                catch
                    is_valid_step = false; break;
                end
            end

            if ~is_valid_step
                alpha = alpha * beta;
                continue;
            end

            [f_new, ~] = module_objective.compute(Gamma_new, Sigmas, Kernel, W, params);

            if f_new < f_curr || abs(f_new - f_curr)/max(abs(f_curr),1) < 1e-6
                success = true;
                break;
            else
                alpha = alpha * beta;
            end
        end

        if ~success
            if params.verbose
                fprintf('[PGD] Line search failed at iter %d. Step too small.\n', iter);
            end
            break;
        end

        % C. Update
        diff_norm = 0; norm_old = 0;
        for f=1:F
            diff_norm = diff_norm + norm(Gamma_new{f} - Gamma_curr{f}, 'fro')^2;
            norm_old = norm(Gamma_curr{f}, 'fro')^2;
        end
        rel_diff = sqrt(diff_norm) / sqrt(norm_old + 1e-10);

        Gamma_curr = Gamma_new;
        f_curr = f_new;
        obj_history(end+1) = f_curr;

        if params.verbose && (iter==1 || mod(iter,10)==0)
            fprintf('%4d | Obj: %.4e | Diff: %.2e | Alpha: %.1e\n', ...
                iter, f_curr, rel_diff, alpha);
        end

        if rel_diff < params.tol
            break;
        end

        if bt == 0
            alpha = alpha * 1.5;
        end
    end

    Gamma_cells = Gamma_curr;
    results.objective_history = obj_history;
    results.final_iter = iter;
    results.final_alpha = alpha;
end
