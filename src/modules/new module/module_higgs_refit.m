function [Gamma_refit, aux] = module_higgs_refit(S_whitened_cell, mask_cell, K, W, params)
% MODULE_HIGGS_REFIT
%   Analytic / semi-analytic refit paths after Rayleigh or rescue support selection.
%
% Modes:
%   1) 'single_ridge'
%      - HIGGS-like path.
%      - Build a single-frequency ridge template from Sigmae_f.
%      - Apply the final mask.
%      - Hermitianize + SPD-project.
%
%   2) 'fused_ridge_surrogate'
%      - Start from the same single-frequency ridge template.
%      - For each edge (i,j), smooth its frequency trajectory analytically by
%        solving the graph-Tikhonov problem
%
%            min_x  1/2 ||x - b||_2^2
%                  + lambda1 * w_ij * x^H L_K x
%                  + lambda3 * w_ij * ||x||_2^2
%
%        on the active subset of frequencies.
%      - This has the closed-form solution
%
%            x = (I + 2*w_ij*(lambda1*L_K + lambda3*I))^{-1} b
%
%        (restricted to the active subset when masks differ across frequencies).
%
%      This is NOT the exact closed-form solution of the full JSPACE logdet
%      problem. It is an analytic Frobenius-surrogate refit designed to stay
%      close to the HIGGS ridge idea while preserving cross-frequency coupling.
%
% Inputs:
%   S_whitened_cell : {F x 1} whitened covariances
%   mask_cell       : {F x 1} final masks
%   K               : frequency kernel
%   W               : entrywise weight matrix
%   params          : optional struct
%       .mode = 'single_ridge' | 'fused_ridge_surrogate'
%       .ridge_penalty
%       .lambda1
%       .lambda3
%       .min_eig
%
% Outputs:
%   Gamma_refit : {F x 1}
%   aux         : diagnostics

    if nargin < 5 || isempty(params), params = struct(); end

    mode          = lower(get_opt_(params, 'mode', 'fused_ridge_surrogate'));
    ridge_penalty = get_opt_(params, 'ridge_penalty', 1e-2);
    lambda1       = get_opt_(params, 'lambda1', 0);
    lambda3       = get_opt_(params, 'lambda3', 0);
    min_eig       = get_opt_(params, 'min_eig', 1e-8);

    F = numel(S_whitened_cell);
    p = size(S_whitened_cell{1}, 1);

    if isempty(mask_cell)
        mask_cell = cell(F,1);
        for f = 1:F
            mask_cell{f} = true(p);
        end
    end

    [Gamma_base, base_aux] = module_higgs_ridge_template(S_whitened_cell, ...
        struct('ridge_penalty', ridge_penalty, 'min_eig', min_eig));

    switch mode
        case 'single_ridge'
            Gamma_refit = apply_mask_and_project_(Gamma_base, mask_cell, min_eig);
            aux = struct();
            aux.mode = mode;
            aux.base = base_aux;

        case 'fused_ridge_surrogate'
            Ksym = (K + K') / 2;
            if F <= 1
                Lk = 0;
            else
                Lk = diag(sum(Ksym, 2)) - Ksym;
            end

            Gamma_refit = analytic_freq_smoothing_(Gamma_base, mask_cell, Lk, W, lambda1, lambda3, min_eig);

            aux = struct();
            aux.mode = mode;
            aux.base = base_aux;
            aux.Lk = Lk;

        otherwise
            error('module_higgs_refit:UnknownMode', 'Unknown refit mode: %s', mode);
    end
end

% ============================================================
% helpers
% ============================================================

function Gamma_out = apply_mask_and_project_(Gamma_base, mask_cell, min_eig)
    F = numel(Gamma_base);
    Gamma_out = cell(F,1);

    for f = 1:F
        G0 = Gamma_base{f};
        M  = mask_cell{f};
        p  = size(G0,1);

        M(1:p+1:end) = true;

        G = G0;
        G(~M) = 0;
        G(1:p+1:end) = real(diag(G0));
        G = utils_math.make_hermitian(G);
        [G, ~] = utils_math.project_spd(G, min_eig);
        Gamma_out{f} = G;
    end
end

function Gamma_out = analytic_freq_smoothing_(Gamma_base, mask_cell, Lk, W, lambda1, lambda3, min_eig)
    F = numel(Gamma_base);
    p = size(Gamma_base{1}, 1);

    Gamma_raw = cell(F,1);
    for f = 1:F
        Gamma_raw{f} = zeros(p, p, 'like', Gamma_base{f});
        Gamma_raw{f}(1:p+1:end) = real(diag(Gamma_base{f}));
    end

    for i = 1:p
        for j = (i+1):p
            active = false(F,1);
            b = zeros(F,1);

            for f = 1:F
                M = mask_cell{f};
                if M(i,j)
                    active(f) = true;
                    b(f) = Gamma_base{f}(i,j);
                end
            end

            if ~any(active), continue; end

            idx = find(active);
            if isempty(W)
                wij = 1;
            else
                wij = max(real(W(i,j)), 0);
            end

            if numel(idx) == 1 || wij == 0 || (lambda1 <= 0 && lambda3 <= 0)
                x = zeros(F,1);
                x(idx) = b(idx);
            else
                A = eye(numel(idx)) + 2 * wij * (lambda1 * Lk(idx, idx) + lambda3 * eye(numel(idx)));
                A = 0.5 * (A + A');
                x = zeros(F,1);
                x(idx) = A \ b(idx);
            end

            for f = 1:F
                if active(f)
                    Gamma_raw{f}(i,j) = x(f);
                    Gamma_raw{f}(j,i) = conj(x(f));
                end
            end
        end
    end

    Gamma_out = cell(F,1);
    for f = 1:F
        G = utils_math.make_hermitian(Gamma_raw{f});
        [G, ~] = utils_math.project_spd(G, min_eig);
        Gamma_out{f} = G;
    end
end

function val = get_opt_(s, f, d)
    if isfield(s, f), val = s.(f); else, val = d; end
end
