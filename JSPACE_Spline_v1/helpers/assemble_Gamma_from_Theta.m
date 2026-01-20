function [Gamma_cell, invGamma_cell, logdet_cell, minEig_all, ok_spd] = ...
    assemble_Gamma_from_Theta(Theta_cell, B, eps_pd)
%ASSEMBLE_GAMMA_FROM_THETA  Synthesize precision matrices from spline coeffs.
%
% Formula:
%   Gamma(t) = sum_{k=1}^{K} B(t,k) * Theta{k}
%
% Inputs:
%   Theta_cell : {K x 1} cell of (N x N) Hermitian matrices
%   B          : (F x K) Basis matrix
%   eps_pd     : scalar, minimum threshold for SPD check (e.g., 1e-6)
%
% Outputs:
%   Gamma_cell    : {F x 1} cell of Gamma(t)
%   invGamma_cell : {F x 1} cell of inv(Gamma(t))
%   logdet_cell   : {F x 1} cell of log(det(Gamma(t)))
%   minEig_all    : scalar, worst-case minimum eigenvalue across all F (estimate)
%   ok_spd        : boolean, true if ALL Gamma(t) are SPD after stabilization

    if nargin < 3 || isempty(eps_pd)
        eps_pd = 1e-6;
    end

    F = size(B, 1);
    K = size(B, 2);
    if numel(Theta_cell) ~= K
        error('assemble_Gamma_from_Theta:SizeMismatch', ...
            'B has %d columns but Theta_cell has %d elements.', K, numel(Theta_cell));
    end
    N = size(Theta_cell{1}, 1);

    Gamma_cell = cell(F, 1);
    invGamma_cell = cell(F, 1);
    logdet_cell = cell(F, 1);

    ok_spd = true;
    minEig_all = inf;

    I = eye(N, 'like', Theta_cell{1});

    % Loop over frequencies (serial for stability and reproducibility)
    for t = 1:F
        % 1. Assemble Gamma_t
        Gt = zeros(N, N, 'like', Theta_cell{1});
        for k = 1:K
            w = B(t, k);
            if abs(w) > 1e-12
                Gt = Gt + w * Theta_cell{k};
            end
        end

        % Force Hermitian and real diagonal (numerical safety)
        Gt = 0.5 * (Gt + Gt');
        Gt(1:N+1:end) = real(diag(Gt));

        % 2. Check SPD and stabilize per frequency
        [R, p] = chol(Gt); % Gt = R' * R
        if p == 0
            d = real(diag(R));
            minEig_t = min(d)^2;
        else
            minEig_t = estimate_min_eig_(Gt);
        end

        % Apply diagonal shift if needed (per-frequency stabilization)
        shift = max(0, eps_pd - minEig_t);
        if shift > 0
            Gt = Gt + shift * I;
        end

        % Re-factor with jitter escalation if needed
        [R, p] = chol(Gt);
        if p > 0
            jitter = max(shift, eps_pd);
            max_tries = 5;
            for j = 1:max_tries
                Gt_try = Gt + jitter * I;
                [R, p] = chol(Gt_try);
                if p == 0
                    Gt = Gt_try;
                    break;
                end
                jitter = jitter * 10;
            end
        end

        if p > 0
            ok_spd = false;
            invG = [];
            lgDet = -inf;
        else
            d = real(diag(R));
            lgDet = 2 * sum(log(d));
            invG = R \ (R' \ I);
        end

        % Update global stats
        if minEig_t < minEig_all
            minEig_all = minEig_t;
        end

        % Store results
        Gamma_cell{t} = Gt;
        invGamma_cell{t} = invG;
        logdet_cell{t} = lgDet;
    end
end

function minEig = estimate_min_eig_(G)
    % Estimate smallest eigenvalue for Hermitian matrices.
    minEig = NaN;
    try
        minEig = eigs(G, 1, 'smallestreal');
    catch
        try
            minEig = min(real(eig(G)));
        catch
            minEig = NaN;
        end
    end
    if ~isfinite(minEig)
        minEig = -1;
    end
    minEig = real(minEig);
end
