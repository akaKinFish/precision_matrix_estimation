classdef module_objective
    % MODULE_OBJECTIVE Computes the scalar value of the objective function.
    %
    % Implements:
    %   F(Gamma) = sum_w [-logdet(G_w) + tr(S_w * G_w)]
    %              + lambda1 * sum_{w<w'} k_{w,w'} * ||G_w - G_w'||_W^2
    %              + lambda3 * sum_w ||G_w||_W^2
    %              + lambda2 * L1_offdiag(G_w)

    methods (Static)

        function [total_obj, stats] = compute(Gamma_cells, Sigma_cells, Kernel, W, params)
            % COMPUTE Calculates the objective value and its breakdown.
            %
            % Inputs:
            %   Gamma_cells : {F x 1} current precision matrices
            %   Sigma_cells : {F x 1} whitened covariances
            %   Kernel      : (F x F) smoothing kernel
            %   W           : (p x p) weight matrix
            %   params      : struct with lambda1, lambda2, lambda3, weight_mode, penalize_diagonal
            %
            % Outputs:
            %   total_obj   : scalar objective value
            %   stats       : struct containing individual term values

            % 1. Setup
            F = numel(Gamma_cells);
            p = size(Gamma_cells{1}, 1);

            if ~isfield(params, 'lambda1'), lambda1 = 0; else, lambda1 = params.lambda1; end
            if ~isfield(params, 'lambda2'), lambda2 = 0; else, lambda2 = params.lambda2; end
            if ~isfield(params, 'lambda3'), lambda3 = 0; else, lambda3 = params.lambda3; end

            if ~isfield(params, 'weight_mode')
                mode = 'matrix';
            else
                mode = params.weight_mode;
            end

            penalize_diag = isfield(params, 'penalize_diagonal') && params.penalize_diagonal;

            % Accumulators
            obj_logdet = 0;
            obj_trace  = 0;
            obj_smooth = 0;  % frequency smoothing
            obj_space  = 0;  % spatial smoothing
            obj_l1     = 0;

            % -----------------------------------------------------------
            % 2. Log-likelihood terms
            % -----------------------------------------------------------
            for f = 1:F
                G = Gamma_cells{f};
                S = Sigma_cells{f};

                [ld_val, is_valid] = utils_math.safe_log_det(G);
                if ~is_valid
                    total_obj = Inf;
                    stats = struct('logdet', Inf, 'trace', 0, 'smooth', 0, 'space', 0, 'l1', 0);
                    return;
                end

                obj_logdet = obj_logdet - ld_val;
                obj_trace  = obj_trace + real(trace(S * G));
            end

            % -----------------------------------------------------------
            % 3. Frequency smoothing term (pairwise)
            % -----------------------------------------------------------
            if lambda1 > 0
                Ksym = (Kernel + Kernel') / 2;
                for f1 = 1:F
                    for f2 = f1+1:F
                        if Ksym(f1, f2) == 0, continue; end
                        Gdiff = Gamma_cells{f1} - Gamma_cells{f2};
                        if strcmp(mode, 'matrix')
                            term = real(trace(Gdiff' * (W * Gdiff)));
                        elseif strcmp(mode, 'hadamard')
                            term = real(sum(sum(conj(Gdiff) .* (W .* Gdiff))));
                        else
                            error('ModuleObjective:UnknownMode', 'Unknown weight mode');
                        end
                        obj_smooth = obj_smooth + Ksym(f1, f2) * term; % each pair once, gradient expects 2*lambda1 from squared norm
                    end
                end
                obj_smooth = lambda1 * obj_smooth;
            end

            % -----------------------------------------------------------
            % 4. Spatial smoothing term
            % -----------------------------------------------------------
            if lambda3 > 0
                for f = 1:F
                    G = Gamma_cells{f};
                    if strcmp(mode, 'matrix')
                        term = real(trace(G' * (W * G)));
                    elseif strcmp(mode, 'hadamard')
                        term = real(sum(sum(conj(G) .* (W .* G))));
                    else
                        error('ModuleObjective:UnknownMode', 'Unknown weight mode');
                    end
                    obj_space = obj_space + term;
                end
                obj_space = lambda3 * obj_space;
            end

            % -----------------------------------------------------------
            % 5. L1 sparsity term
            % -----------------------------------------------------------
            if lambda2 > 0
                for f = 1:F
                    G = Gamma_cells{f};
                    if ~penalize_diag
                        G(1:p+1:end) = 0;
                    end
                    obj_l1 = obj_l1 + sum(abs(G(:)));
                end
                obj_l1 = lambda2 * obj_l1;
            end

            % -----------------------------------------------------------
            % 6. Final summation
            % -----------------------------------------------------------
            total_obj = obj_logdet + obj_trace + obj_smooth + obj_space + obj_l1;

            stats.logdet = obj_logdet;
            stats.trace  = obj_trace;
            stats.smooth = obj_smooth;
            stats.space  = obj_space;
            stats.l1     = obj_l1;
        end
    end
end
