classdef module_proximal_operator
    % MODULE_PROXIMAL_OPERATOR Implements the proximal gradient update step.
    %
    % This module executes the core update logic:
    %   Gamma_{k+1} = Proj_SPD( SoftThresh( Gamma_k - alpha * Grad ) )
    %
    % It handles:
    %   1. Gradient Descent
    %   2. L1 Regularization (via Soft Thresholding)
    %   3. Diagonal handling (usually unpenalized)
    %   4. Active Set constraints (forcing zeros)
    %   5. Hermitian & SPD constraints
    %
    % Dependencies: utils_math.m

    methods (Static)
        
        function Gamma_new = compute(Gamma_curr, Grad, alpha, lambda2, active_mask, opts)
            % COMPUTE Executes one proximal gradient update step.
            %
            % Inputs:
            %   Gamma_curr  : (p x p) Current Precision Matrix
            %   Grad        : (p x p) Gradient of the smooth part
            %   alpha       : (scalar) Step size
            %   lambda2     : (scalar) L1 penalty coefficient
            %   active_mask : (logical p x p) Mask of active edges (optional)
            %   opts        : Struct with options:
            %       .penalize_diagonal : (bool) Apply L1 to diagonal? (Default: false)
            %       .min_eig           : (double) Floor for SPD projection (Default: 1e-8)
            %
            % Output:
            %   Gamma_new   : (p x p) Updated, valid Precision Matrix
            
            if nargin < 6, opts = struct(); end
            
            % -------------------------------------------------------
            % 0. Parse Options
            % -------------------------------------------------------
            % Standard GLasso usually does NOT penalize the diagonal.
            penalize_diag = isfield(opts, 'penalize_diagonal') && opts.penalize_diagonal;
            
            % Minimum eigenvalue floor to ensure log-det validity
            if ~isfield(opts, 'min_eig'),  min_eig = 1e-8;
            else 
                min_eig = opts.min_eig; 
            end
            
            % -------------------------------------------------------
            % 1. Gradient Descent Step
            % -------------------------------------------------------
            % Moving in the direction of negative gradient
            G_step = Gamma_curr - alpha * Grad;
            
            % -------------------------------------------------------
            % 2. Proximal Operator (L1 Soft Thresholding)
            % -------------------------------------------------------
            % Threshold value depends on step size
            tau = alpha * lambda2;
            
            if lambda2 > 0
                % Apply complex soft thresholding to the ENTIRE matrix first
                Gamma_prox = utils_math.soft_threshold_complex(G_step, tau);
                
                % Handle Diagonal Logic
                if ~penalize_diag
                    % If diagonal is NOT penalized, it should not be shrunk.
                    % We restore the diagonal elements from the gradient descent step (G_step)
                    % effectively applying prox_0(z) = z.
                    p = size(Gamma_curr, 1);
                    diag_idx = 1:p+1:p*p;
                    Gamma_prox(diag_idx) = G_step(diag_idx);
                end
            else
                % No L1 penalty, pure gradient descent step
                Gamma_prox = G_step;
            end
            
            % -------------------------------------------------------
            % 3. Active Set Projection
            % -------------------------------------------------------
            % If an active set is provided, we strictly enforce zeros outside it.
            % This is crucial for the efficiency of the Active Set method.
            if ~isempty(active_mask)
                % active_mask is true for elements to keep, false for zeros
                Gamma_prox(~active_mask) = 0;
            end
            
            % -------------------------------------------------------
            % 4. Enforce Symmetry and Real Diagonal
            % -------------------------------------------------------
            % Numerical operations might introduce slight asymmetry or complex diagonals
            Gamma_prox = utils_math.make_hermitian(Gamma_prox);
            
            % Precision matrices must have real diagonals
            p = size(Gamma_prox, 1);
            diag_idx = 1:p+1:p*p;
            Gamma_prox(diag_idx) = real(Gamma_prox(diag_idx));
            
            % -------------------------------------------------------
            % 5. SPD Projection (Constraint Enforcement)
            % -------------------------------------------------------
            % The most expensive but necessary step: ensure matrix is Positive Definite.
            % Uses eigen-decomposition if Cholesky fails.
            [Gamma_new, ~] = utils_math.project_spd(Gamma_prox, min_eig);
            
        end
    end
end