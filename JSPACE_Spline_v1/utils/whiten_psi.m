function [D_cell, Psi_tilde_cell] = whiten_psi(Psi_cell, whiten_cfg)
%WHITEN_PSI  Diagonal whitening for posterior source covariance.
%
% Inputs:
%   Psi_cell   : {F x 1} cell, each (N x N) Hermitian
%   whiten_cfg : struct with optional field .eps_diag (default 1e-12)
%
% Outputs:
%   D_cell         : {F x 1} diagonal whitening matrices
%   Psi_tilde_cell : {F x 1} whitened matrices D * Psi * D

    if nargin < 2 || isempty(whiten_cfg)
        whiten_cfg = struct();
    end
    if ~isfield(whiten_cfg, 'eps_diag') || isempty(whiten_cfg.eps_diag)
        whiten_cfg.eps_diag = 1e-12;
    end

    F = numel(Psi_cell);
    D_cell = cell(F, 1);
    Psi_tilde_cell = cell(F, 1);

    eps_val = whiten_cfg.eps_diag;

    for t = 1:F
        Psi = Psi_cell{t};
        N = size(Psi, 1);

        diag_p = real(diag(Psi));
        diag_p(~isfinite(diag_p)) = eps_val;
        diag_p = max(diag_p, eps_val);

        inv_sqrt_diag = 1 ./ sqrt(diag_p);
        D = diag(inv_sqrt_diag);

        Psi_tilde = D * Psi * D;
        Psi_tilde = 0.5 * (Psi_tilde + Psi_tilde');
        Psi_tilde(1:N+1:end) = real(diag(Psi_tilde));

        D_cell{t} = D;
        Psi_tilde_cell{t} = Psi_tilde;
    end
end
