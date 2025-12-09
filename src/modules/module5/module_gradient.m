classdef module_gradient
    % MODULE_GRADIENT Computes the gradient of the smooth objective function.
    %
    % Implements the gradient from the documentation:
    %   grad = -Gamma^{-1} + Sigma
    %          + 2 * lambda1 * sum_{w'} k_{w,w'} * W * (Gamma_w - Gamma_w')
    %          + 2 * lambda3 * W * Gamma_w
    %
    % Dependencies:
    %   - utils_math.m

    methods (Static)

        function Grads = compute(Gamma_cells, Sigma_cells, Kernel, W, params)
            % COMPUTE Calculates the gradient matrices for all frequencies.
            %
            % Inputs:
            %   Gamma_cells : {F x 1} Cell array of precision matrices
            %   Sigma_cells : {F x 1} Cell array of whitened covariances
            %   Kernel      : (F x F) Smoothing kernel matrix (k_{w,w'})
            %   W           : (p x p) Weight matrix (W^Gamma)
            %   params      : Struct with fields:
            %                   .lambda1     : Frequency smoothing strength
            %                   .lambda3     : Spatial smoothing strength
            %                   .weight_mode : 'matrix' or 'hadamard'
            %                   .penalize_diagonal : logical, for L1 handling (not used here)
            %
            % Output:
            %   Grads       : {F x 1} Cell array of gradient matrices

            % 1. Validation and setup
            F = numel(Gamma_cells);
            p = size(Gamma_cells{1}, 1);
            Grads = cell(F, 1);

            if ~isfield(params, 'lambda1'), lambda1 = 0; else, lambda1 = params.lambda1; end
            if ~isfield(params, 'lambda3'), lambda3 = 0; else, lambda3 = params.lambda3; end
            penalize_diag = isfield(params, 'penalize_diagonal') && params.penalize_diagonal;

            if isfield(params, 'weight_mode')
                mode = params.weight_mode;
            else
                mode = 'matrix';
            end

            % 2. Pre-calculate Laplacian components for frequency smoothing
            K = (Kernel + Kernel') / 2;      % ensure symmetry
            degrees = sum(K, 2);             % d_w = sum_{w'} k_{w,w'}

            NeighborSum = cell(F, 1);
            for f = 1:F
                NeighborSum{f} = zeros(p, p);
            end

            for row = 1:F
                for col = 1:F
                    if K(row, col) ~= 0
                        NeighborSum{row} = NeighborSum{row} + K(row, col) * Gamma_cells{col};
                    end
                end
            end

            % 3. Gradient per frequency
            for f = 1:F
                G = Gamma_cells{f};
                S = Sigma_cells{f};

                % A. Log-det and trace parts (use SPD projection + stable solve)
                if isfield(params, 'min_eig'), min_eig = params.min_eig; else, min_eig = 1e-8; end
                [G_pd, ~] = utils_math.project_spd(G, min_eig);
                % Use linear solve instead of explicit inv to improve conditioning
                invG = G_pd \ eye(p);
                grad_data_fitting = -invG + S;

                % B. Frequency smoothing gradient
                grad_smooth = zeros(p, p);
                if lambda1 > 0
                    DiffTerm = degrees(f) * G - NeighborSum{f};  % d_w * G_w - sum k G_w'

                    if strcmp(mode, 'matrix')
                        grad_smooth = 2 * lambda1 * (W * DiffTerm);
                    elseif strcmp(mode, 'hadamard')
                        grad_smooth = 2 * lambda1 * (W .* DiffTerm);
                    else
                        error('ModuleGradient:UnknownWeightMode', ...
                            'weight_mode must be "matrix" or "hadamard"');
                    end
                end

                % C. Spatial smoothing gradient
                grad_space = zeros(p, p);
                if lambda3 > 0
                    if strcmp(mode, 'matrix')
                        grad_space = 2 * lambda3 * (W * G);
                    elseif strcmp(mode, 'hadamard')
                        grad_space = 2 * lambda3 * (W .* G);
                    else
                        error('ModuleGradient:UnknownWeightMode', ...
                            'weight_mode must be "matrix" or "hadamard"');
                    end
                end

                % D. Combine and enforce symmetry
                TotalGrad = grad_data_fitting + grad_smooth + grad_space;

                % Keep the diagonal gradient; diagonal L1 handling is done in the proximal step.
                if ~penalize_diag
                    % No diagonal zeroing here.
                end

                Grads{f} = utils_math.make_hermitian(TotalGrad);
            end
        end
    end
end
