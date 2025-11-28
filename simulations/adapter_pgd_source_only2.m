function out = adapter_pgd_source_only2(args, opts)
% ADAPTER_PGD_SOURCE_ONLY2 Bridge for PGD Core with Grid Search + EBIC
%
% Strategy:
%   Performs a "Grid Search" over a reasonable range of lambdas and selects 
%   the best model using EBIC (Extended Bayesian Information Criterion).
%
% Inputs:
%   args.Sjj_tilde   : (p x p) or Cell{F} - Whitened Covariance
%   args.K           : (F x F) - Smoothing Kernel
%   args.W           : (p x p) - Weight Matrix
%   args.D_src       : (p x p) or Cell{F} - Whitening Matrix (for recolor)
%   args.active_mask : (optional) Active set mask
%
%   opts.m_samples   : (int) Sample size (Required for EBIC calculation)
%   opts.recolor     : bool (default true)
%   opts.weight_mode : 'matrix' or 'hadamard' (default 'matrix')

    if nargin < 2, opts = struct(); end

    % ============================================================
    % 1. Data Preparation & Validation
    % ============================================================
    [Sjj, F, p] = ensure_cell_covariances(args);
    
    K = args.K;
    W = args.W;
    
    D_src = [];
    if isfield(args, 'D_src') && ~isempty(args.D_src)
        D_src = ensure_cell_general(args.D_src);
    end
    
    Gamma_init = cell(F, 1);
    for f=1:F, Gamma_init{f} = eye(p); end
    
    if isfield(opts, 'm_samples')
        m_val = opts.m_samples;
    else
        warning('Sample size m not provided. Using heuristic m=100*p for EBIC.');
        m_val = 100 * p;
    end

    % ============================================================
    % 2. Grid Search Setup
    % ============================================================
    % Define search range [lambda_min, lambda_max].
    
    S_temp = Sjj{1};
    mask_tril = tril(true(p), -1);
    off_diag_abs = abs(S_temp(mask_tril));
    lambda_max_data = max(off_diag_abs);
    
    % Lower Bound (Min): Jiang's Limit (Theoretical Noise Floor).
    lambda_min_theory = sqrt(4 * log(p) / m_val);
    
    if lambda_min_theory >= lambda_max_data
        lambda_min_theory = lambda_max_data * 0.1; 
    end
    
    % Define Grid: 5 points logarithmically spaced, Descending
    grid_num = 5;
    lambda_grid = logspace(log10(lambda_max_data * 0.9), log10(lambda_min_theory * 1.1), grid_num);
    lambda_grid = sort(lambda_grid, 'descend'); 
    
    fprintf('[Adapter] Starting Grid Search (Size=%d) for EBIC Selection...\n', grid_num);
    fprintf('          Range: [%.4e, %.4e]\n', min(lambda_grid), max(lambda_grid));
    fprintf('          %-10s | %-10s | %-10s | %-12s\n', 'Lambda', 'Density(%)', 'Time(s)', 'EBIC');
    fprintf('          %s\n', repmat('-',1,55));
    
    best_ebic = Inf;
    best_lambda = lambda_grid(1);
    best_Gamma = Gamma_init; 
    best_stats = [];
    
    current_Gamma_init = Gamma_init;
    
    % ============================================================
    % 3. Grid Loop
    % ============================================================
    for k = 1:grid_num
        lam_curr = lambda_grid(k);
        
        % Prepare Input
        m5_input.whitened_covariances = Sjj;
        m5_input.precision_matrices   = current_Gamma_init; 
        m5_input.smoothing_kernel     = K;
        m5_input.weight_matrix        = W;
        if isfield(args, 'active_mask')
            m5_input.active_mask = ensure_cell_general(args.active_mask); 
        end

        % Prepare Params
        m5_params.lambda1   = 0.01; 
        m5_params.lambda2   = lam_curr;
        m5_params.alpha0    = 1e-3; 
        m5_params.auto_tune = true; 
        m5_params.max_iter  = 150;  
        m5_params.tol       = 1e-4; 
        m5_params.verbose   = false;
        
        % %%% <<< CHANGE HERE: Pass weight_mode from opts to m5_params >>>
        if isfield(opts, 'weight_mode')
            m5_params.weight_mode = opts.weight_mode;
        end
        
        % RUN PGD
        t_start = tic;
        [Gamma_temp, stats_temp] = module5_proximal_main(m5_input, m5_params);
        t_cost = toc(t_start);
        
        % --- EBIC Calculation ---
        G_final = Gamma_temp{1};
        
        try
            L_chol = chol((G_final+G_final')/2);
            ld_val = 2 * sum(log(real(diag(L_chol))));
        catch
            ld_val = -1e10; 
        end
        tr_val = real(trace(Sjj{1} * G_final));
        
        G_off = G_final; 
        G_off(1:p+1:end) = 0; 
        num_edges = sum(abs(G_off(:)) > 1e-5) / 2;
        density = (num_edges / (p*(p-1)/2)) * 100;
        
        % EBIC (High-dimensional Gamma = 0.5)
        gamma_ebic = 0.5; 
        log_lik = (m_val/2) * (ld_val - tr_val);
        minus_2_ll = -2 * log_lik;
        penalty = num_edges * log(m_val) + 4 * num_edges * gamma_ebic * log(p);
        
        current_ebic = minus_2_ll + penalty;
        
        % Select Best
        is_best = '';
        if current_ebic < best_ebic
            best_ebic = current_ebic;
            best_lambda = lam_curr;
            best_Gamma = Gamma_temp;
            best_stats = stats_temp;
            is_best = '(*)';
        end
        
        fprintf('          %.4e | %5.2f%%     | %5.2f      | %.4e %s\n', ...
                lam_curr, density, t_cost, current_ebic, is_best);
        
        current_Gamma_init = Gamma_temp;
    end
    
    fprintf('          %s\n', repmat('-',1,55));
    fprintf('[Adapter] Best Lambda Selected: %.4e\n', best_lambda);

    % ============================================================
    % 4. Recolor
    % ============================================================
    Gamma_final = best_Gamma;
    Omega_src = [];
    
    do_recolor = ~isfield(opts, 'recolor') || opts.recolor;
    
    if do_recolor && ~isempty(D_src)
        Omega_src = cell(F, 1);
        for f = 1:F
            if numel(D_src) == 1 && F > 1
                Df = D_src{1};
            else
                Df = D_src{f};
            end
            Gf = Gamma_final{f};
            Omega_src{f} = Df * Gf * Df;
        end
    end

    out.Gamma_tilde_star = Gamma_final;
    out.Omega_src        = Omega_src;
    out.proximal_stats   = best_stats;
    out.used_params.lambda2 = best_lambda;
    out.used_params.grid_search_best_ebic = best_ebic;
end

% ================= Helpers =================
function [C, F, p] = ensure_cell_covariances(args)
    if ~isfield(args,'Sjj_tilde') || isempty(args.Sjj_tilde)
        error('adapter:missing_Sjj','args.Sjj_tilde required');
    end
    X = args.Sjj_tilde;
    
    if iscell(X)
        C = X; F = numel(C); p = size(C{1},1);
    elseif isnumeric(X) && ndims(X)==3
        [p,~,F] = size(X); C = cell(F,1);
        for f=1:F, Cf = X(:,:,f); C{f} = 0.5*(Cf+Cf'); end
    elseif isnumeric(X) && ismatrix(X) && size(X,1) == size(X,2)
        p = size(X,1);
        F = 1;
        C = {0.5*(X+X')};
    else
        error('adapter:bad_Sjj','Sjj_tilde format error');
    end
end

function C = ensure_cell_general(Input)
    if iscell(Input)
        C = Input;
    else
        C = {Input};
    end
end