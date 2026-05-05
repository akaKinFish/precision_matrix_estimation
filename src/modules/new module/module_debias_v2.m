function [Gamma_debiased, Gamma_rayleigh, Var_proxies, masks] = module_debias_v2(Gamma_hat_cell, S_whitened_cell, n_samples, params)
% MODULE_DEBIAS - One-step debiasing only (fixed Rayleigh threshold removed)
%
% Purpose:
%   Correct the shrinkage bias introduced by the off-diagonal L1 penalty.
%   This module NO LONGER performs any fixed Rayleigh thresholding.
%
% Theory:
%   For each frequency f,
%       Gamma_db,f = 2*Gamma_hat,f - Gamma_hat,f * Sigmae_f * Gamma_hat,f
%
%   A variance proxy is still returned for downstream Rayleigh search:
%       V_f,ij = |Gsrc_ii| |Gsrc_jj| + |Gsrc_ij|^2,
%   where Gsrc is chosen by params.variance_source:
%       'hat'      -> Gamma_hat,f     (default; closest to the current JSPACE docs/code)
%       'debiased' -> Gamma_db,f
%
% Inputs:
%   Gamma_hat_cell   : {F x 1} biased whitened precision estimates
%   S_whitened_cell  : {F x 1} final whitened covariances
%   n_samples        : retained only for API compatibility
%   params           : optional struct
%       .variance_source = 'hat' | 'debiased'
%       .force_real_diag = true | false
%
% Outputs:
%   Gamma_debiased : {F x 1} debiased dense matrices
%   Gamma_rayleigh : [] kept only for backward compatibility
%   Var_proxies    : {F x 1} variance proxies for downstream Rayleigh search
%   masks          : [] kept only for backward compatibility

    if nargin < 4 || isempty(params), params = struct(); end
    variance_source = get_opt_(params, 'variance_source', 'hat');
    force_real_diag = get_opt_(params, 'force_real_diag', true);

    F = numel(Gamma_hat_cell);
    Gamma_debiased = cell(F, 1);
    Var_proxies    = cell(F, 1);

    % Backward-compatibility placeholders:
    Gamma_rayleigh = [];
    masks          = [];

    %#ok<NASGU> % kept for backward-compatibility and possible future use
    n_eff = n_samples;

    for f = 1:F
        G = Gamma_hat_cell{f};
        S = S_whitened_cell{f};

        % Keep CPU/GPU alignment consistent with the original code.
        if isa(G, 'gpuArray') && ~isa(S, 'gpuArray')
            S = gpuArray(S);
        elseif ~isa(G, 'gpuArray') && isa(S, 'gpuArray')
            S = gather(S);
        end

        % --- 1) One-step debiasing ---
        G_db = 2 * G - G * S * G;
        G_db = (G_db + G_db') / 2;

        if force_real_diag
            p = size(G_db, 1);
            G_db(1:p+1:end) = real(diag(G_db));
        end
        Gamma_debiased{f} = G_db;

        % --- 2) Variance proxy for downstream Rayleigh search ---
        switch lower(variance_source)
            case 'debiased'
                G_var = G_db;
            otherwise
                G_var = G;
        end

        d = abs(diag(G_var));
        V = d * d' + abs(G_var).^2;
        V = real((V + V') / 2);
        Var_proxies{f} = V;
    end
end

function val = get_opt_(s, f, d)
    if isfield(s, f), val = s.(f); else, val = d; end
end
