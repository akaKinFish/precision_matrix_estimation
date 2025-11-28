classdef module_gradient
    % MODULE_GRADIENT Computes the gradient of the smooth objective function.
    %
    % This module strictly implements the gradient derivation found on Page 7
    % of the provided documentation ("Simplified Objective and Updates").
    %
    % The gradient formula is:
    %   grad = -Gamma^{-1} + Sigma + 2 * lambda1 * sum(k * W * (Gamma_w - Gamma_wp))
    %
    % Dependencies: 
    %   - utils_math.m

    methods (Static)
        
        function Grads = compute(Gamma_cells, Sigma_cells, Kernel, W, params)
            % COMPUTE Calculates the gradient matrices for all frequencies.
            %
            % Inputs:
            %   Gamma_cells : {F x 1} Cell array of Precision matrices (Gamma)
            %   Sigma_cells : {F x 1} Cell array of Whitened Covariances (Sigma)
            %   Kernel      : (F x F) Smoothing kernel matrix (k_{w,w'})
            %   W           : (p x p) Weight matrix (W^Gamma)
            %   params      : Struct containing hyper-parameters:
            %       .lambda1     : Smoothing strength (scalar)
            %       .weight_mode : 'matrix' (standard) or 'hadamard' (element-wise)
            %
            % Outputs:
            %   Grads       : {F x 1} Cell array of gradient matrices
            
            % 1. Validation and Setup
            F = numel(Gamma_cells);
            p = size(Gamma_cells{1}, 1);
            Grads = cell(F, 1);
            
            % Parameters (Fixed: Standard if-else)
            if isfield(params, 'lambda1')
                lambda1 = params.lambda1;
            else
                lambda1 = 0;
            end
            
            if isfield(params, 'weight_mode')
                mode = params.weight_mode;
            else
                mode = 'matrix'; % Default to matrix multiplication as per Page 7 formula
            end
            
            % 2. Pre-calculate Laplacian Matrix components for Efficient Smoothing
            % Ensure Kernel is symmetric
            K = (Kernel + Kernel') / 2;
            
            % Degree vector: d_w = sum_{w'} k_{w,w'}
            degrees = sum(K, 2);
            
            % Pre-calculate the weighted sum of neighbors: sum_{w'} k_{w,w'} * G_{w'}
            % This allows us to compute the sum term efficiently.
            NeighborSum = cell(F, 1);
            
            for f = 1:F
                NeighborSum{f} = zeros(p, p);
            end
            
            % Accumulate neighbor contributions (O(F^2) loop outside, fast matrix add inside)
            for row = 1:F
                for col = 1:F
                    if K(row, col) ~= 0
                        % NeighborSum{row} += k_{row,col} * Gamma{col}
                        NeighborSum{row} = NeighborSum{row} + K(row, col) * Gamma_cells{col};
                    end
                end
            end
            
            % 3. Compute Gradient for Each Frequency
            for f = 1:F
                G = Gamma_cells{f};
                S = Sigma_cells{f};
                
                % -----------------------------------------------------------
                % Part A: Log-Determinant Gradient
                % Formula: -Gamma^{-1} (Page 7)
                % -----------------------------------------------------------
                % Calculate Inverse. G is SPD (guaranteed by proximal step).
                % Using inv(G) is generally acceptable for gradients.
                % For extreme stability, G \ eye(p) is better but slower.
                invG = inv(G);
                
                % -----------------------------------------------------------
                % Part B: Trace Gradient
                % Formula: +Sigma (Page 7)
                % -----------------------------------------------------------
                grad_data_fitting = -invG + S;
                
                % -----------------------------------------------------------
                % Part C: Smoothing Gradient
                % Formula: + 2 * lambda1 * sum_{w'} k_{w,w'} * W * (G_w - G_w')
                % Rewritten via Laplacian: + 2 * lambda1 * W * (d_w * G_w - NeighborSum{w})
                % -----------------------------------------------------------
                grad_smooth = zeros(p, p);
                
                if lambda1 > 0
                    % Calculate the difference term: (d_w * G_w - sum(k*G'))
                    DiffTerm = degrees(f) * G - NeighborSum{f};
                    
                    if strcmp(mode, 'matrix')
                        % Formula from Page 7: W * (...)
                        % Matrix multiplication implies spatial correlation weighting
                        grad_smooth = 2 * lambda1 * (W * DiffTerm);
                        
                    elseif strcmp(mode, 'hadamard')
                        % Element-wise weighting: W .* (...)
                        % Implies edge-specific penalties (often used in GLasso)
                        % Note: If mode is hadamard, W corresponds to W_{ijkl} in docs
                        grad_smooth = 2 * lambda1 * (W .* DiffTerm);
                        
                    else
                        error('ModuleGradient:UnknownWeightMode', ...
                              'weight_mode must be "matrix" or "hadamard"');
                    end
                end
                
                % -----------------------------------------------------------
                % Part D: Combine and Enforce Symmetry
                % -----------------------------------------------------------
                TotalGrad = grad_data_fitting + grad_smooth;
                
                % The gradient of a real-valued function w.r.t a Hermitian matrix 
                % must be Hermitian. We enforce this to prevent numerical drift.
                Grads{f} = utils_math.make_hermitian(TotalGrad);
            end
        end
    end
end