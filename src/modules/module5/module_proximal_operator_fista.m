classdef module_proximal_operator_fista
    % MODULE_PROXIMAL_OPERATOR_FISTA
    % Pure proximal map for the nonsmooth term g(G) = tau * ||G||_{1,off}
    % with structural projections (active mask, Hermitian, SPD, optional unit diag).
    %
    % Usage (single frequency):
    %   Gamma_new = module_proximal_operator_fista.compute(G_point, tau, active_mask, opts)
    %
    % Inputs:
    %   G_point      : (p x p) complex Hermitian candidate after gradient step
    %   tau          : scalar soft-threshold strength (lambda2 / L)
    %   active_mask  : logical p x p mask (true = keep), [] to skip masking
    %   opts         : struct with fields (all optional)
    %       .penalize_diagonal    : default false (do not shrink diag)
    %       .min_eig              : SPD floor (default 1e-8)
    %       .enforce_unit_diagonal: default false; if true, renormalize diag to 1
    %
    % Output:
    %   Gamma_new    : Hermitian SPD matrix after prox and projections
    methods (Static)
        function [Gamma_new, info] = compute(G_point, tau, active_mask, opts)
            if nargin < 4, opts = struct(); end
            if nargin < 3, active_mask = []; end

            penalize_diag = isfield(opts, 'penalize_diagonal') && opts.penalize_diagonal;
            if isfield(opts, 'min_eig'), min_eig = opts.min_eig; else, min_eig = 1e-8; end
            use_unit_diag = isfield(opts, 'enforce_unit_diagonal') && opts.enforce_unit_diagonal;

            % 1) Soft-threshold off-diagonals (complex)
            if tau > 0
                G_soft = utils_math.soft_threshold_complex(G_point, tau);
                if ~penalize_diag
                    p = size(G_point, 1);
                    diag_idx = 1:p+1:p*p;
                    G_soft(diag_idx) = G_point(diag_idx);
                end
            else
                G_soft = G_point;
            end

            % 2) Active mask projection
            if ~isempty(active_mask)
                G_soft(~active_mask) = 0;
            end

            % 3) Hermitian + real diagonal
            G_h = utils_math.make_hermitian(G_soft);
            p = size(G_h, 1);
            diag_idx = 1:p+1:p*p;
            G_h(diag_idx) = real(G_h(diag_idx));

            % 4) SPD projection / eigenvalue floor
            [G_spd, spd_stats] = utils_math.project_spd(G_h, min_eig);

            % 5) Optional unit-diagonal renormalization
            %    Useful when theory assumes correlation-style precision (diag=1).
            %    Risk: rescales overall magnitude and may conflict with likelihood scaling.
            if use_unit_diag
                floor_val = max(min_eig, eps);
                d = real(diag(G_spd));
                d = max(d, floor_val);
                D_inv_sqrt = diag(1 ./ sqrt(d));
                G_spd = utils_math.make_hermitian(D_inv_sqrt * G_spd * D_inv_sqrt);
                [G_spd, spd_stats] = utils_math.project_spd(G_spd, min_eig);
            end

            Gamma_new = G_spd;

            if nargout > 1
                info = struct();
                info.clipped = spd_stats.clipped;
                lam = eig((Gamma_new + Gamma_new')/2);
                info.min_eig = min(real(lam));
            end
        end
    end
end
