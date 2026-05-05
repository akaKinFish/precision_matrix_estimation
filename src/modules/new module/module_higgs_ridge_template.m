function [Gamma_base, aux] = module_higgs_ridge_template(S_whitened_cell, params)
% MODULE_HIGGS_RIDGE_TEMPLATE
%   Build a HIGGS-like single-frequency ridge / Frobenius template.
%
% Main default path (target-free, spectral):
%   For each frequency f, if Sigmae_f = U diag(s) U^H, define
%       gamma(s) = (sqrt(s.^2 + 4*rho) - s) / (2*rho),
%   and
%       Gamma_base,f = U diag(gamma(s)) U^H.
%
% This is exactly the spectral ridge/Frobenius path used in the HIGGS code
% snippet shared by the user.
%
% Optional target-aware path:
%   If params.target_cell is provided, we also expose the general ridge
%   matrix formula (Bilgrau et al., 2020, Lemma 8):
%       Omega = ( sqrtm(rho*I + 1/4*(S - rho*T)^2) + 1/2*(S - rho*T) )^{-1}.
%
% Inputs:
%   S_whitened_cell : {F x 1} whitened covariance matrices
%   params          : optional struct
%       .ridge_penalty    scalar > 0
%       .min_eig          SPD floor
%       .target_cell      optional {F x 1} targets
%       .use_target_formula logical
%
% Outputs:
%   Gamma_base : {F x 1} SPD ridge templates
%   aux        : diagnostic info

    if nargin < 2 || isempty(params), params = struct(); end

    rho               = get_opt_(params, 'ridge_penalty', 1e-2);
    min_eig           = get_opt_(params, 'min_eig', 1e-8);
    use_target_formula = get_opt_(params, 'use_target_formula', false);
    target_cell       = get_opt_(params, 'target_cell', []);

    rho = max(real(rho), eps);
    F = numel(S_whitened_cell);

    Gamma_base = cell(F, 1);
    aux = struct();
    aux.ridge_penalty = rho;
    aux.used_target_formula = false;
    aux.spectral_eigvals = cell(F,1);

    for f = 1:F
        S = utils_math.make_hermitian(S_whitened_cell{f});
        p = size(S, 1);

        if use_target_formula && ~isempty(target_cell)
            T = utils_math.make_hermitian(target_cell{f});
            A = utils_math.make_hermitian(S - rho * T);

            % Build sqrt term in a numerically stable Hermitian form.
            B = utils_math.make_hermitian(rho * eye(p, 'like', S) + 0.25 * (A * A));
            [UB, dB] = eig(B, 'vector');
            dB = max(real(dB), 0);
            sqrtB = UB * diag(sqrt(dB)) * UB';
            sqrtB = utils_math.make_hermitian(sqrtB);

            Den = utils_math.make_hermitian(sqrtB + 0.5 * A);
            [Den_spd, ~] = utils_math.project_spd(Den, min_eig);

            G = Den_spd \ eye(p, 'like', S);
            G = utils_math.make_hermitian(G);
            [G, ~] = utils_math.project_spd(G, min_eig);

            aux.used_target_formula = true;
            aux.spectral_eigvals{f} = [];
        else
            [U, ds] = eig(S, 'vector');
            ds = real(ds);
            dg = (sqrt(ds.^2 + 4 * rho) - ds) ./ (2 * rho);
            dg = max(real(dg), min_eig);

            G = U * diag(dg) * U';
            G = utils_math.make_hermitian(G);
            [G, ~] = utils_math.project_spd(G, min_eig);

            aux.spectral_eigvals{f} = ds;
        end

        Gamma_base{f} = G;
    end
end

function val = get_opt_(s, f, d)
    if isfield(s, f), val = s.(f); else, val = d; end
end
