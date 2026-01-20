function [Theta_cell, shift_val] = spd_shift_on_Theta0(Theta_cell, B, eps_pd)
%SPD_SHIFT_ON_THETA0  Conservative diagonal shift on spline coefficients.
%
% This routine does NOT assume an explicit intercept basis. Instead, it
% computes per-frequency SPD deficits and distributes a diagonal shift
% across all Theta_k in a least-squares sense using the basis matrix.
%
% Inputs:
%   Theta_cell : {K x 1}
%   B          : Basis matrix
%   eps_pd     : Target minimum eigenvalue
%
% Outputs:
%   Theta_cell : Modified coefficients
%   shift_val  : Amount added (for diagnostics)

    if nargin < 3 || isempty(eps_pd)
        eps_pd = 1e-6;
    end

    F = size(B, 1);
    K = size(B, 2);
    if numel(Theta_cell) ~= K
        error('spd_shift_on_Theta0:SizeMismatch', ...
            'B has %d columns but Theta_cell has %d elements.', K, numel(Theta_cell));
    end
    N = size(Theta_cell{1}, 1);

    shift_val = 0;
    delta_vec = zeros(F, 1);

    for t = 1:F
        Gt = zeros(N, N, 'like', Theta_cell{1});
        for k = 1:K
            Gt = Gt + B(t, k) * Theta_cell{k};
        end
        Gt = 0.5 * (Gt + Gt');
        Gt(1:N+1:end) = real(diag(Gt));

        [R, p] = chol(Gt);
        if p == 0
            d = real(diag(R));
            min_eig = min(d) ^ 2;
        else
            min_eig = estimate_min_eig_(Gt);
        end

        if ~isfinite(min_eig)
            min_eig = -abs(eps_pd);
        end

        if min_eig < eps_pd
            delta_vec(t) = eps_pd - min_eig;
        end
    end

    if all(delta_vec <= 0)
        return;
    end

    % Least-squares distribution of diagonal shifts across basis columns.
    delta_theta = pinv(B) * delta_vec;
    delta_theta(~isfinite(delta_theta)) = 0;
    delta_theta(delta_theta < 0) = 0;

    shift_val = max(delta_theta);
    if shift_val <= 0
        return;
    end

    I = eye(N, 'like', Theta_cell{1});
    for k = 1:K
        if delta_theta(k) ~= 0
            Theta_cell{k} = Theta_cell{k} + delta_theta(k) * I;
            Theta_cell{k} = 0.5 * (Theta_cell{k} + Theta_cell{k}');
            Theta_cell{k}(1:N+1:end) = real(diag(Theta_cell{k}));
        end
    end
end

function minEig = estimate_min_eig_(G)
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
