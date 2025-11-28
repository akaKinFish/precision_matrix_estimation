classdef utils_math
    % UTILS_MATH Static utility class for low-level mathematical operations.
    %
    % This class encapsulates stateless mathematical functions required for
    % the Precision Matrix Estimation algorithm, including Hermitian enforcement,
    % SPD projection, and complex proximal operators.
    %
    % It corresponds to the mathematical definitions in the provided documentation,
    % specifically ensuring numerical stability for the whitened optimization problem.

    methods (Static)
        
        function A_out = make_hermitian(A_in)
            % MAKE_HERMITIAN Forces a matrix to be Hermitian.
            %
            % Mathematical formulation:
            %   A_out = (A_in + A_in^H) / 2
            %
            % Usage:
            %   Ensures that numerical noise does not violate the Hermitian
            %   property required for Cholesky decomposition.
            
            A_out = (A_in + A_in') / 2;
        end

        function [val, is_valid] = safe_log_det(A)
            % SAFE_LOG_DET Computes log(det(A)) using Cholesky decomposition.
            %
            % Mathematical formulation:
            %   If A = L * L^H, then det(A) = det(L) * det(L^H) = (prod(diag(L)))^2
            %   log(det(A)) = 2 * sum(log(real(diag(L))))
            %
            % Inputs:
            %   A - Square matrix (intended to be Hermitian PD)
            %
            % Outputs:
            %   val      - The log-determinant value. Returns -Inf if not PD.
            %   is_valid - Boolean flag indicating if A was Positive Definite.
            
            % Ensure input is Hermitian for the check to be valid
            A_sym = (A + A') / 2;
            
            % Attempt Cholesky decomposition
            % 'lower' produces L such that A = L * L'
            [L, fail] = chol(A_sym, 'lower');
            
            if fail
                val = -Inf;
                is_valid = false;
            else
                % Sum of log of diagonal elements avoids underflow/overflow of det()
                val = 2 * sum(log(real(diag(L))));
                is_valid = true;
            end
        end

        function [A_out, stats] = project_spd(A_in, min_eig)
            % PROJECT_SPD Projects a matrix onto the Symmetric Positive Definite (SPD) cone.
            %
            % Logic:
            %   This function ensures the matrix is valid for the log-det term.
            %   It performs an eigen-decomposition and floors small eigenvalues
            %   to 'min_eig'.
            %
            % Mathematical formulation:
            %   A = V * D * V^H
            %   D_new = max(D, min_eig)
            %   A_out = V * D_new * V^H
            %
            % Inputs:
            %   A_in    - Input matrix (pxp)
            %   min_eig - Minimum eigenvalue floor (default: 1e-8)
            %
            % Outputs:
            %   A_out - The projected SPD matrix
            %   stats - Struct containing projection info (e.g., if clipping occurred)
            
            if nargin < 2, min_eig = 1e-8; end
            
            % 1. Force symmetry first
            A_sym = (A_in + A_in') / 2;
            
            % 2. Try fast check with Cholesky
            % If it succeeds, the matrix is already PD (eigenvalues > 0)
            [~, p] = chol(A_sym);
            
            stats.clipped = false;
            
            if p == 0 
                % Matrix is technically PD.
                % Note: Strictly speaking, we might still want to enforce >= min_eig
                % if the condition number is bad, but for speed, we usually accept
                % Cholesky success as "good enough".
                A_out = A_sym;
            else
                % 3. Spectral Projection (The heavy lifting)
                stats.clipped = true;
                
                % Use 'vector' to handle potential hermitian numerical issues
                % making sure output D is real.
                [V, D] = eig(A_sym, 'vector');
                
                % Ensure eigenvalues are real (they should be for Hermitian)
                D = real(D);
                
                % Floor eigenvalues
                D = max(D, min_eig);
                
                % Reconstruct: A = V * diag(D) * V'
                % Using BSXFUN or explicit multiplication
                A_out = V * (D .* V');
                
                % Final symmetry enforcement after reconstruction
                A_out = (A_out + A_out') / 2;
            end
        end

        function Z_out = soft_threshold_complex(Z_in, tau)
            % SOFT_THRESHOLD_COMPLEX Proximal operator for L1 norm on complex numbers.
            %
            % Corresponds to the operator "soft" mentioned in Page 7 of the documentation.
            %
            % Mathematical formulation:
            %   prox_tau(z) = z * max(0, 1 - tau / |z|)
            %               = (z / |z|) * max(0, |z| - tau)
            %
            % Inputs:
            %   Z_in - Input matrix or vector (Complex or Real)
            %   tau  - Threshold value (Scalar or Matrix of same size as Z_in)
            %
            % Outputs:
            %   Z_out - Thresholded matrix
            
            % Compute magnitude
            abs_Z = abs(Z_in);
            
            % Compute the scaling factor: (|z| - tau)_+
            % We compute max(0, |z| - tau) first.
            magnitude_shrunk = max(0, abs_Z - tau);
            
            % Avoid division by zero where Z_in is 0.
            % Create a safe mask for non-zero elements
            idx = abs_Z > 0;
            
            Z_out = zeros(size(Z_in), 'like', Z_in);
            
            % Apply scaling:
            % Z_out = (Z_in ./ abs_Z) .* magnitude_shrunk
            % This preserves the phase (Z_in ./ abs_Z) and scales the magnitude.
            Z_out(idx) = Z_in(idx) .* (magnitude_shrunk(idx) ./ abs_Z(idx));
        end
    end
end