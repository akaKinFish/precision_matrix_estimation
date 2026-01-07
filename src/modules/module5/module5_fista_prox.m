classdef module5_fista_prox
    % MODULE5_FISTA_PROX
    % Pure proximal map for the non-smooth term g(G) = lambda2 * ||G||_{1,off}
    % plus structural constraints (Hermitian, SPD, active mask, optional unit diag).
    %
    % This operator ASSUMES the gradient step has already been taken:
    %   Z = Y - alpha * grad_f(Y)
    % and returns:
    %   prox_{alpha*g}(Z)
    %
    % Differences from module_proximal_operator:
    %   - No embedded gradient step.
    %   - Strictly applies L1 soft-threshold + projections.
    %   - Designed for FISTA / accelerated schemes.
    %
    % Inputs:
    %   Z           : (p x p) matrix after gradient step
    %   alpha       : step size (scalar)
    %   lambda2     : L1 coefficient
    %   active_mask : logical mask (p x p) of allowed entries (optional)
    %   opts        : struct with fields:
    %                   .penalize_diagonal (default: false)
    %                   .min_eig           (default: 1e-8)
    %                   .unit_diagonal     (default: false) re-normalize diag to 1
    %
    % Outputs:
    %   Gamma_new   : prox result (Hermitian SPD)
    %   info        : struct with projection diagnostics
    methods (Static)
        function [Gamma_new, info] = compute(Z, alpha, lambda2, active_mask, opts)
            if nargin < 5, opts = struct(); end
            if nargin < 4, active_mask = []; end

            penalize_diag = isfield(opts, 'penalize_diagonal') && opts.penalize_diagonal;
            if ~isfield(opts, 'min_eig'), min_eig = 1e-8; else, min_eig = opts.min_eig; end
            use_unit_diag = isfield(opts, 'unit_diagonal') && opts.unit_diagonal;

            % 1) Soft-threshold off-diagonals
            if lambda2 > 0
                tau = alpha * lambda2;
                G_soft = utils_math.soft_threshold_complex(Z, tau);
                if ~penalize_diag
                    p = size(Z, 1);
                    diag_idx = 1:p+1:p*p;
                    G_soft(diag_idx) = Z(diag_idx); % keep diagonal untouched
                end
            else
                G_soft = Z;
            end

            % 2) Active set projection (force zeros outside mask)
            if ~isempty(active_mask)
                G_soft(~active_mask) = 0;
            end

            % 3) Hermitian + real diagonal
            G_h = utils_math.make_hermitian(G_soft);
            p = size(G_h, 1);
            diag_idx = 1:p+1:p*p;
            G_h(diag_idx) = real(G_h(diag_idx));

            % 4) SPD projection
            [G_spd, spd_stats] = utils_math.project_spd(G_h, min_eig);

            % 5) Optional unit-diagonal rescaling (correlation-style)
            %    Helps match theoretical diag=1 but changes scale; keep opt-in.
            if use_unit_diag
                floor_val = max(min_eig, eps);
                d = real(diag(G_spd));
                d = max(d, floor_val);
                D_inv_sqrt = diag(1 ./ sqrt(d));
                G_spd = utils_math.make_hermitian(D_inv_sqrt * G_spd * D_inv_sqrt);
                [G_spd, spd_stats] = utils_math.project_spd(G_spd, min_eig);
            end

            Gamma_new = G_spd;
            info = struct();
            info.clipped = spd_stats.clipped;
            lam = eig((Gamma_new + Gamma_new')/2);
            info.min_eig = min(real(lam));
        end
    end
end
