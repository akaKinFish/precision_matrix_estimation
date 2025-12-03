function [Gamma_debiased, Gamma_rayleigh] = module_debias(Gamma_hat_cell, S_whitened_cell, n_samples)
% MODULE_DEBIAS - One-step Debiasing & Rayleigh Thresholding
%
% Purpose:
%   Corrects the shrinkage bias introduced by L1 regularization (Lasso)
%   and applies a statistically grounded threshold to remove noise.
%
% Theory (Docs Chapter 7):
%   1. Debiasing: Gamma_tilde = 2*Gamma - Gamma * S * Gamma
%   2. Variance:  Var_ij = G_ii*G_jj + |G_ij|^2
%   3. Threshold: |G_ij| > (r_th/sqrt(T)) * sqrt(Var_ij)
%
% Inputs:
%   Gamma_hat_cell  : {F x 1} Sparse Precision Matrices (from PGD)
%   S_whitened_cell : {F x 1} Whitened Covariance Matrices
%   n_samples       : (int) Sample size T
%
% Outputs:
%   Gamma_debiased  : {F x 1} Dense, bias-corrected matrices
%   Gamma_rayleigh  : {F x 1} Sparse, thresholded matrices

    F = numel(Gamma_hat_cell);
    Gamma_debiased = cell(F, 1);
    Gamma_rayleigh = cell(F, 1);
    
    % Robust default threshold (r=3.0 corresponds to p < 0.01 significance)
    r_th = 3.0; 
    
    for f = 1:F
        G = Gamma_hat_cell{f};
        S = S_whitened_cell{f};
        
        % Ensure inputs are on the same device (CPU/GPU)
        if isa(G, 'gpuArray') && ~isa(S, 'gpuArray')
            S = gpuArray(S);
        elseif ~isa(G, 'gpuArray') && isa(S, 'gpuArray')
            S = gather(S);
        end
        
        % --- 1. One-Step Debiasing (Eq 7.1) ---
        % Removes L1 shrinkage bias
        G_tilde = 2 * G - G * S * G;
        
        % Enforce Hermitian symmetry
        G_tilde = (G_tilde + G_tilde') / 2;
        Gamma_debiased{f} = G_tilde;
        
        % --- 2. Variance Estimation (Eq 7.10) ---
        % Asymptotic variance proxy for the debiased estimator
        d = diag(G_tilde);
        Var_proxy = real(d * d') + abs(G_tilde).^2;
        
        % --- 3. Rayleigh Thresholding (Eq 7.12) ---
        % Statistical significance test
        Threshold = (r_th / sqrt(complex(n_samples))) * sqrt(complex(Var_proxy));
        
        % Create Mask
        mask = abs(G_tilde) > Threshold;
        
        % Always keep diagonal (self-loops are essential)
        p = size(G, 1);
        mask(1:p+1:end) = true;
        
        % Apply Mask
        G_ray = G_tilde;
        G_ray(~mask) = 0;
        
        Gamma_rayleigh{f} = G_ray;
    end
end