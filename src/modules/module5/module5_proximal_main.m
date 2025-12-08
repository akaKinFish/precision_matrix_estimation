function [Gamma_cells, results] = module5_proximal_main(input_data, params)
% MODULE5_PROXIMAL_MAIN Main Orchestrator for Precision Matrix Estimation.
%
% This function implements the Proximal Gradient Descent algorithm.
% It coordinates the Gradient, StepSize, and ProximalOperator modules.
%
% Features:
%   - Gershgorin-based Auto-Tuning (Doc Page 8)
%   - Barzilai-Borwein (BB) Step Size
%   - Active Set Strategy integration
%   - Convergence monitoring
%
% Inputs:
%   input_data : Struct containing:
%       .precision_matrices     (Initial guess, optional)
%       .whitened_covariances   (Required, {F x 1})
%       .smoothing_kernel       (Required, F x F)
%       .weight_matrix          (Required, p x p)
%       .active_mask            (Optional, {F x 1} logical)
%
%   params : Struct containing hyperparameters:
%       .lambda1           : Smoothing penalty (if not auto-tuned)
%       .lambda2           : L1 sparsity penalty
%       .auto_tune         : (bool) Use Gershgorin theory for lambda1/alpha?
%       .max_iter          : Max iterations (default 100)
%       .tol               : Convergence tolerance (default 1e-5)
%       .active_set_update : (bool) Enable dynamic active set?
%
% Outputs:
%   Gamma_cells : Optimized Precision Matrices {F x 1}
%   results     : Struct with history and stats

    % ============================================================
    % 1. Initialization & Validation
    % ============================================================
    
    % Basic input unpacking
    Sigmas = input_data.whitened_covariances;
    Kernel = input_data.smoothing_kernel;
    W      = input_data.weight_matrix;
    
    F = numel(Sigmas);
    p = size(Sigmas{1}, 1);
    
    % Initialize Gamma if not provided (Identity is a safe start for Whitened data)
    if isfield(input_data, 'precision_matrices') && ~isempty(input_data.precision_matrices)
        Gamma_curr = input_data.precision_matrices;
    else
        Gamma_curr = cell(F, 1);
        for f = 1:F, Gamma_curr{f} = eye(p); end
    end
    
    % Initialize Active Mask
    if isfield(input_data, 'active_mask') && ~isempty(input_data.active_mask)
        active_mask = input_data.active_mask;
    else
        active_mask = []; % Empty means "all active"
    end
    
    % Default Parameters
    if nargin < 2, params = struct(); end
    defaults = struct('lambda1', 0.01, ...      % freq smooth
                  'lambda2', 0.01, ...      % L1 sparsity
                  'lambda3', 0.00, ...      % NEW: spatial smooth (ridge-style)
                  'max_iter', 200, ...
                  'tol', 1e-5, ...
                  'verbose', true, ...
                  'auto_tune', false, ...
                  'active_set_freq', 10);
    
    % Merge defaults
    fnames = fieldnames(defaults);
    for k=1:numel(fnames)
        if ~isfield(params, fnames{k}), params.(fnames{k}) = defaults.(fnames{k}); end
    end

    % ============================================================
    % 2. Auto-Tuning via Gershgorin Circle Theorem (Doc Page 8)
    % ============================================================
    if params.auto_tune
        if params.verbose
            fprintf('\n[Init] Performing Gershgorin Auto-Tuning...\n');
        end
        
        % Call module_step_size to calculate theoretical bounds
        % Safety margin delta=0.9 (as per docs)
        [auto_params, gersh_stats] = module_step_size.compute_gershgorin_params(Gamma_curr, Kernel, W, 0.9);
        
        % Override params
        params.lambda1 = auto_params.lambda1;
        params.alpha0  = auto_params.alpha0;
        
        if params.verbose
            fprintf('       Selected lambda1 = %.4e based on K_max=%.2f, R_max=%.2f\n', ...
                    params.lambda1, gersh_stats.K_max, gersh_stats.R_max);
            fprintf('       Selected alpha0  = %.4e based on L_logdet=%.2e\n', ...
                    params.alpha0, gersh_stats.L_logdet);
        end
    else
        if ~isfield(params, 'alpha0'), params.alpha0 = 1e-3; end
    end

    % ============================================================
    % 3. Main Optimization Loop
    % ============================================================
    
    Gamma_prev = [];
    Grad_prev  = [];
    alpha      = params.alpha0;
    
    obj_history = zeros(params.max_iter, 1);
    diff_history = zeros(params.max_iter, 1);
    
    if params.verbose
        fprintf('\n[Main] Starting Proximal Gradient Descent (MaxIter=%d, Tol=%.1e)\n', ...
                params.max_iter, params.tol);
        fprintf('%5s | %12s | %10s | %10s\n', 'Iter', 'Objective', 'Diff', 'Alpha');
        fprintf('-----------------------------------------------------\n');
    end
    
    for iter = 1:params.max_iter
        
        % ---------------------------------------------------
        % A. Gradient Computation (Smooth Part)
        % ---------------------------------------------------
        % module_gradient handles LogDet, Trace, and Smoothing terms
        Grad_curr = module_gradient.compute(Gamma_curr, Sigmas, Kernel, W, params);
        if params.verbose && iter == 1
            % Diagnostic: gradient and iterate norms at start
            grad_norm = 0; gamma_norm = 0;
            for f=1:F
                grad_norm = grad_norm + norm(Grad_curr{f}, 'fro')^2;
                gamma_norm = gamma_norm + norm(Gamma_curr{f}, 'fro')^2;
            end
            grad_norm = sqrt(grad_norm); gamma_norm = sqrt(gamma_norm);
            fprintf('    [Diag] grad_norm=%.3e, gamma_norm=%.3e\n', grad_norm, gamma_norm);
        end
        
        % ---------------------------------------------------
        % B. Step Size Update (Barzilai-Borwein)
        % ---------------------------------------------------
        if iter > 1
            % Calculate adaptive step size based on secant equation
            alpha = module_step_size.compute_bb(Gamma_curr, Gamma_prev, ...
                                                Grad_curr, Grad_prev, params);
        end
        
        % Store history for next BB step
        Gamma_prev = Gamma_curr;
        Grad_prev  = Grad_curr;
        
        % ---------------------------------------------------
        % C. Proximal Update (The "Heavy Lifting")
        % ---------------------------------------------------
        % Gamma_new = Proj_SPD( SoftThresh( Gamma - alpha*Grad ) )
        Gamma_new = cell(F, 1);
        total_diff_sq = 0;
        total_norm_sq = 0;
        
        for f = 1:F
            % Select mask for this frequency
            mask_f = [];
            if ~isempty(active_mask), mask_f = active_mask{f}; end
            
            Gamma_new{f} = module_proximal_operator.compute(...
                Gamma_curr{f}, ...
                Grad_curr{f}, ...
                alpha, ...
                params.lambda2, ...
                mask_f, ...
                params ...
            );
            
            % Accumulate difference for convergence check
            diff_mtx = Gamma_new{f} - Gamma_curr{f};
            total_diff_sq = total_diff_sq + norm(diff_mtx, 'fro')^2;
            total_norm_sq = total_norm_sq + norm(Gamma_curr{f}, 'fro')^2;
        end
        
        rel_diff = sqrt(total_diff_sq) / sqrt(total_norm_sq + 1e-10);
        Gamma_curr = Gamma_new; % Accept step (Simplistic PGD without line search)
        
        % ---------------------------------------------------
        % D. Active Set Update (Optional)
        % ---------------------------------------------------
        % Periodically update the active set to include/exclude edges
        % Call external module5_update_active_set if needed
        if isfield(params, 'active_set_update') && params.active_set_update ...
           && mod(iter, params.active_set_freq) == 0
            
            % Assuming module5_update_active_set exists and implements the KKT logic
            % [active_mask, ~] = module5_update_active_set(active_mask, Gamma_curr, Grad_curr, params);
            % if params.verbose, fprintf('       [ActiveSet] Updated masks.\n'); end
        end
        
        % ---------------------------------------------------
        % E. Logging & Convergence
        % ---------------------------------------------------
        [obj_val, ~] = module_objective.compute(Gamma_curr, Sigmas, Kernel, W, params);
        obj_history(iter) = obj_val;
        diff_history(iter) = rel_diff;
        
        if params.verbose && (mod(iter, 10) == 0 || iter == 1)
            fprintf('%5d | %12.5e | %10.2e | %10.2e\n', ...
                    iter, obj_val, rel_diff, alpha);
        end
        
        if rel_diff < params.tol
            if params.verbose
                fprintf('\n[Converged] Relative difference %.2e < tolerance %.2e at iter %d.\n', ...
                        rel_diff, params.tol, iter);
            end
            break;
        end
        
    end
    
    % ============================================================
    % 4. Final Results
    % ============================================================
    Gamma_cells = Gamma_curr;
    
    results.objective_history = obj_history(1:iter);
    results.diff_history = diff_history(1:iter);
    results.final_iter = iter;
    results.final_alpha = alpha;
    if params.auto_tune
        results.gershgorin_stats = gersh_stats;
    end

end
