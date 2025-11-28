classdef module_step_size
    % MODULE_STEP_SIZE Computes optimization step sizes and hyper-parameters.
    % (FIXED: Removed invalid ternary operators)

    methods (Static)
        
        function [params_out, stats] = compute_gershgorin_params(Gamma_init, Kernel, W, delta)
            % COMPUTE_GERSHGORIN_PARAMS Automatically selects lambda1 and alpha0.
            
            if nargin < 4, delta = 0.9; end
            
            F = numel(Gamma_init);
            
            % 1. Compute L_logdet
            L_logdet = 0;
            for f = 1:F
                G = Gamma_init{f};
                G = (G + G') / 2;
                e = eig(G);
                min_e = min(e);
                if min_e < 1e-12
                    min_e = 1e-12;
                end
                L_f = (1 / min_e)^2;
                L_logdet = max(L_logdet, L_f);
            end
            
            % 2. Compute Smoothing Bounds
            K = (Kernel + Kernel') / 2;
            K_max = max(sum(abs(K), 2));
            
            W_sym = (W + W') / 2;
            R_max = max(sum(abs(W_sym), 2));
            
            % 3. Calculate Parameters
            if K_max == 0 || R_max == 0
                lambda1 = 0; 
                alpha0 = 1 / (L_logdet + delta); 
            else
                lambda1 = delta / (2 * K_max * R_max);
                alpha0 = 1 / (L_logdet + delta);
            end
            
            params_out.lambda1 = lambda1;
            params_out.alpha0 = alpha0;
            
            stats.L_logdet = L_logdet;
            stats.K_max = K_max;
            stats.R_max = R_max;
            stats.delta = delta;
        end
        
        function alpha = compute_bb(Gamma_curr, Gamma_prev, Grad_curr, Grad_prev, opts)
            % COMPUTE_BB Calculates the Barzilai-Borwein step size.
            
            % Defaults (FIXED: Using standard if-else)
            if isfield(opts, 'alpha_min'), alpha_min = opts.alpha_min; else, alpha_min = 1e-12; end
            if isfield(opts, 'alpha_max'), alpha_max = opts.alpha_max; else, alpha_max = 1.0; end
            
            % Fallback for first iteration
            if isempty(Gamma_prev) || isempty(Grad_prev)
                if isfield(opts, 'alpha0'), alpha = opts.alpha0; else, alpha = 1e-3; end
                return;
            end
            
            dot_ss = 0; 
            dot_sy = 0; 
            
            F = numel(Gamma_curr);
            
            for f = 1:F
                S = Gamma_curr{f} - Gamma_prev{f};
                Y = Grad_curr{f} - Grad_prev{f};
                
                s_vec = S(:);
                y_vec = Y(:);
                
                dot_ss = dot_ss + real(dot(s_vec, s_vec)); 
                dot_sy = dot_sy + real(dot(s_vec, y_vec)); 
            end
            
            if abs(dot_sy) < 1e-15
                if isfield(opts, 'alpha_prev')
                    alpha = opts.alpha_prev;
                else
                    alpha = alpha_min;
                end
            else
                alpha = dot_ss / dot_sy;
            end
            
            if alpha <= 0
                alpha = alpha_max; 
            end
            
            alpha = max(alpha_min, min(alpha_max, alpha));
        end
        
        function L_est = estimate_lipschitz(x0, gradient_func, num_samples, var_scale)
            % ESTIMATE_LIPSCHITZ Estimates L via random sampling
            max_ratio = 0;
            F = numel(x0);
            p = size(x0{1}, 1);
            
            for i = 1:num_samples
                dX = cell(F, 1);
                x_perturbed = cell(F, 1);
                norm_dX_sq = 0;
                
                for f = 1:F
                    noise = (randn(p) + 1i*randn(p)) * var_scale;
                    noise = utils_math.make_hermitian(noise);
                    dX{f} = noise;
                    x_perturbed{f} = x0{f} + noise;
                    norm_dX_sq = norm_dX_sq + norm(noise, 'fro')^2;
                end
                norm_dX = sqrt(norm_dX_sq);
                
                grad_0 = gradient_func(x0);
                grad_p = gradient_func(x_perturbed);
                
                norm_dG_sq = 0;
                for f = 1:F
                    dG = grad_p{f} - grad_0{f};
                    norm_dG_sq = norm_dG_sq + norm(dG, 'fro')^2;
                end
                norm_dG = sqrt(norm_dG_sq);
                
                if norm_dX > 1e-15
                    ratio = norm_dG / norm_dX;
                    max_ratio = max(max_ratio, ratio);
                end
            end
            
            L_est = max_ratio;
            if L_est == 0, L_est = 1; end
        end
    end
end