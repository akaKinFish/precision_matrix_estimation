classdef module_gradient
    % MODULE_GRADIENT Computes the gradient of the smooth objective function.
    %
    % Implements the gradient from the documentation:
    %   grad = -Gamma^{-1} + Sigma
    %          + 2 * lambda1 * sum_{w'} k_{w,w'} * W * (Gamma_w - Gamma_w')
    %          + 2 * lambda3 * W * Gamma_w
    %
    % Revisions:
    %   - Added robust SPD projection with condition number clipping.
    %   - Added Cholesky-based inversion to prevent RCOND warnings.
    %   - Added debug diagnostics for numerical instability.
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
            %                   .min_eig     : Minimum eigenvalue floor
            %                   .cond_cap    : Max allowed condition number (optional)
            %                   .debug       : Enable verbose diagnostics
            %
            % Output:
            %   Grads       : {F x 1} Cell array of gradient matrices
            
            % 1. Validation and setup
            F = numel(Gamma_cells);
            p = size(Gamma_cells{1}, 1);
            Grads = cell(F, 1);
            
            if ~isfield(params, 'lambda1'), lambda1 = 0; else, lambda1 = params.lambda1; end
            if ~isfield(params, 'lambda3'), lambda3 = 0; else, lambda3 = params.lambda3; end
            if isfield(params, 'weight_mode'), mode = params.weight_mode; else, mode = 'matrix'; end
            
            % Gradient calc params
            if isfield(params, 'min_eig'), min_eig = params.min_eig; else, min_eig = 1e-6; end
            if isfield(params, 'cond_cap'), cond_cap = params.cond_cap; else, cond_cap = 1e12; end
            if isfield(params, 'rcond_min'), rcond_min = params.rcond_min; else, rcond_min = 1e-12; end
            dbg = isfield(params,'debug') && params.debug;
            iter_id = -1; if isfield(params,'iter'), iter_id = params.iter; end
            
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
                
                % --- A. Log-det gradient term: robust SPD + cond control + chol solve ---
                % 1) Hermitianize
                G0 = utils_math.make_hermitian(G);

                % 2) Basic Min-Eig Projection
                [G_pd, ~] = utils_math.project_spd(G0, min_eig);

                % 3) Condition Number Clipping (if extremely ill-conditioned)
                rc = rcond(G_pd);
                
                if rc < rcond_min
                    rel_proj = norm(G_pd - G0,'fro') / max(norm(G0,'fro'), 1e-12);
                    [G_pd2, info2] = project_spd_cond_clip_(G0, min_eig, cond_cap);
                    
                    if dbg
                        rc2 = rcond(G_pd2);
                        rel2 = norm(G_pd2 - G0,'fro') / max(norm(G0,'fro'), 1e-12);
                        fprintf('[GradSPD] iter=%d f=%d rcond=%.1e -> %.1e | rel_change=%.2e -> %.2e | eig[min,max]=[%.2e, %.2e] cond=%.2e\n', ...
                            iter_id, f, rc, rc2, rel_proj, rel2, info2.eig_min, info2.eig_max, info2.cond_est);
                    end
                    G_pd = G_pd2;
                end

                % 4) Safe Inversion using Cholesky
                invG = inv_spd_chol_(G_pd);
                
                grad_data_fitting = -invG + S;
                
                % --- B. Frequency smoothing gradient ---
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
                
                % --- C. Spatial smoothing gradient ---
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
                
                % --- D. Combine and enforce symmetry ---
                TotalGrad = grad_data_fitting + grad_smooth + grad_space;
                
                % Keep the diagonal gradient; diagonal L1 handling is done in proximal step.
                Grads{f} = utils_math.make_hermitian(TotalGrad);
            end
        end
    end
end

% ============================================================
% Local Helper Functions
% ============================================================

function [A_clip, info] = project_spd_cond_clip_(A, min_eig, cond_cap)
% PROJECT_SPD_COND_CLIP Projects Hermitian matrix to constrained eigenvalues.
%   Constraints:
%     1. lambda >= min_eig
%     2. lambda <= min_eig * cond_cap (implicitly capping condition number)
%   This prevents the matrix from becoming numerically singular (RCOND ~ 0).

    A = (A + A')/2;
    [V,D] = eig(A,'vector');
    d = real(D);

    d = max(d, min_eig);                    % Lower bound
    max_eig = min_eig * cond_cap;           % Upper bound determined by cond_cap
    d = min(d, max_eig);                    % Upper bound clip

    A_clip = V * diag(d) * V';
    A_clip = (A_clip + A_clip')/2;

    info.eig_min = min(d);
    info.eig_max = max(d);
    info.cond_est = info.eig_max / max(info.eig_min, eps);
end

function invA = inv_spd_chol_(A)
% INV_SPD_CHOL Computes inverse of SPD matrix using Cholesky.
%   Falls back to diagonal loading or standard backslash if Cholesky fails.
%   Numerically safer than inv() or \ for ill-conditioned SPD matrices.

    A = (A + A')/2;
    [R,flag] = chol(A);
    
    if flag ~= 0
        % Cholesky failed (not SPD). Try slight diagonal loading.
        tau = 1e-6 * trace(A)/size(A,1);
        A2 = A + tau*eye(size(A),'like',A);
        [R,flag2] = chol(A2);
        
        if flag2 ~= 0
            % Fallback: standard solve (risky but better than crashing)
            invA = A \ eye(size(A,1),'like',A);
            return;
        end
        R = chol(A2);
    end
    
    % inv(A) = R\(R'\I)
    I = eye(size(A,1),'like',A);
    invA = R \ (R' \ I);
end