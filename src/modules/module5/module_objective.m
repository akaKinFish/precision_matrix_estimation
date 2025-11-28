classdef module_objective
    % MODULE_OBJECTIVE Computes the scalar value of the objective function.
    %
    % This module is used for monitoring convergence. It implements the full
    % objective function F(Gamma) described in the documentation (Page 7).
    %
    % F(Gamma) = LogLikelihood + Smoothing + Sparsity
    %
    % Dependencies: utils_math.m

    methods (Static)

        function [total_obj, stats] = compute(Gamma_cells, Sigma_cells, Kernel, W, params)
            % COMPUTE Calculates the objective value and its breakdown.
            %
            % Inputs:
            %   Gamma_cells : {F x 1} Current Precision Matrices
            %   Sigma_cells : {F x 1} Whitened Covariances
            %   Kernel      : (F x F) Smoothing kernel
            %   W           : (p x p) Weight matrix
            %   params      : Struct with lambda1, lambda2, weight_mode, etc.
            %
            % Outputs:
            %   total_obj   : Scalar objective value
            %   stats       : Struct containing individual term values

            % 1. Setup
            F = numel(Gamma_cells);
            p = size(Gamma_cells{1}, 1);

            if ~isfield(params, 'lambda1'),  lambda1 = 0;
            else 
                lambda1 = params.lambda1;
            end
            if ~isfield(params, 'lambda2'),  lambda2 = 0;else 
                lambda2 = params.lambda2;
            end

            if ~isfield(params, 'weight_mode'),  mode = 'matrix';
            else 
                mode=params.weight_mode;
            end

            penalize_diag = isfield(params, 'penalize_diagonal') && params.penalize_diagonal;

            % Accumulators
            obj_logdet = 0;
            obj_trace  = 0;
            obj_smooth = 0;
            obj_l1     = 0;

            % -----------------------------------------------------------
            % 2. Log-Likelihood (LogDet + Trace)
            % -----------------------------------------------------------
            for f = 1:F
                G = Gamma_cells{f};
                S = Sigma_cells{f};

                % A. Log-Determinant
                % Use safe computation via Cholesky
                [ld_val, is_valid] = utils_math.safe_log_det(G);

                if ~is_valid
                    % If matrix is not PD, the objective is technically Infinity.
                    % We return Inf to signal the optimizer to backtrack.
                    total_obj = Inf;
                    stats = struct('logdet', Inf, 'trace', 0, 'smooth', 0, 'l1', 0);
                    return;
                end

                % Objective has MINUS log det
                obj_logdet = obj_logdet - ld_val;

                % B. Trace Term
                % real(trace(S * G))
                obj_trace = obj_trace + real(trace(S * G));
            end

            % -----------------------------------------------------------
            % 3. Smoothing Term (Laplacian Form)
            % -----------------------------------------------------------
            if lambda1 > 0
                % Prepare Laplacian L = D - K
                K = (Kernel + Kernel') / 2;
                d = sum(K, 2);
                L = diag(d) - K;

                % Calculate: sum_w tr(G_w' * W * sum_wp(L_w,wp * G_wp))
                % This is the quadratic form x'Lx implemented efficiently.

                for f = 1:F
                    % Compute the "Laplacian Neighbor Sum" for frequency f
                    % neighbor_sum = sum_{w'} L(f, w') * G_{w'}
                    neighbor_sum = zeros(p, p);
                    for fp = 1:F
                        if L(f, fp) ~= 0
                            neighbor_sum = neighbor_sum + L(f, fp) * Gamma_cells{fp};
                        end
                    end

                    % Compute inner product based on weight mode
                    if strcmp(mode, 'matrix')
                        % Term: tr( G_f' * W * neighbor_sum )
                        term = real(trace(Gamma_cells{f}' * (W * neighbor_sum)));
                    elseif strcmp(mode, 'hadamard')
                        % Term: sum( conj(G_f) .* (W.^2) .* neighbor_sum )
                        term = real(sum(sum(conj(Gamma_cells{f}) .* (W.^2 .* neighbor_sum))));
                    else
                        error('ModuleObjective:UnknownMode', 'Unknown weight mode');
                    end

                    obj_smooth = obj_smooth + term;
                end

                % Scaling: The gradient was 2*lambda*..., the objective is lambda*...
                % However, the Laplacian form sum_{i,j} L_{ij} <G_i, G_j> inherently
                % includes the factor of 2 relative to the pairwise sum if not careful.
                % Let's stick to the definition:
                % F_smooth = lambda1 * sum_{w,w'} k_{w,w'} ||G_w - G_w'||^2
                %          = lambda1 * 2 * sum_w tr(G_w L_w G_w) (conceptually)
                % The loop above calculates sum_w G_w (L G)_w.
                % So we just multiply by lambda1.
                obj_smooth = lambda1 * obj_smooth;
            end

            % -----------------------------------------------------------
            % 4. L1 Sparsity Term
            % -----------------------------------------------------------
            if lambda2 > 0
                for f = 1:F
                    G = Gamma_cells{f};

                    if ~penalize_diag
                        % Temporarily set diagonal to 0 for calculation
                        G(1:p+1:end) = 0;
                    end

                    % Sum of absolute values (L1 norm)
                    obj_l1 = obj_l1 + sum(abs(G(:)));
                end
                obj_l1 = lambda2 * obj_l1;
            end

            % -----------------------------------------------------------
            % 5. Final Summation
            % -----------------------------------------------------------
            total_obj = obj_logdet + obj_trace + obj_smooth + obj_l1;

            stats.logdet = obj_logdet;
            stats.trace  = obj_trace;
            stats.smooth = obj_smooth;
            stats.l1     = obj_l1;
        end
    end
end