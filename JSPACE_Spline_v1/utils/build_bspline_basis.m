function [B, basis_info] = build_bspline_basis(freq, spline_cfg)
%BUILD_BSPLINE_BASIS  Build B-spline basis matrix (F x K).
%
% Notes:
%   - Uses Cox-de Boor recursion (no toolbox dependency).
%   - Uses uniform knots on normalized frequency in [0, 1].
%
% Inputs:
%   freq       : (F x 1) frequency vector
%   spline_cfg : struct with fields:
%       .K      : number of basis functions
%       .degree : spline degree (e.g., 3 for cubic)
%
% Outputs:
%   B          : (F x K) basis matrix
%   basis_info : struct with knots and metadata

    if nargin < 2 || isempty(spline_cfg)
        spline_cfg = struct();
    end
    if ~isfield(spline_cfg, 'K') || isempty(spline_cfg.K)
        spline_cfg.K = 6;
    end
    if ~isfield(spline_cfg, 'degree') || isempty(spline_cfg.degree)
        spline_cfg.degree = 3;
    end

    K = spline_cfg.K;
    deg = spline_cfg.degree;
    if K <= deg
        error('build_bspline_basis:InvalidConfig', ...
            'K must be larger than degree.');
    end

    f = freq(:);
    F = numel(f);
    if F == 0
        B = zeros(0, K);
        basis_info = struct('knots', [], 'degree', deg, 'f_range', [NaN, NaN]);
        return;
    end

    f_min = min(f);
    f_max = max(f);
    denom = f_max - f_min;
    if denom <= 0
        denom = eps(class(f));
    end
    x = (f - f_min) / denom;

    n_inner_knots = K - deg;
    dt = 1 / n_inner_knots;
    knots = linspace(-deg * dt, 1 + deg * dt, K + deg + 1);

    B = zeros(F, K, 'like', f);
    for k = 1:K
        B(:, k) = bspline_recur(x, k, deg, knots);
    end

    basis_info = struct();
    basis_info.knots = knots;
    basis_info.degree = deg;
    basis_info.f_range = [f_min, f_max];
end

function y = bspline_recur(x, k, d, t)
    % Cox-de Boor recursion for B-spline basis.
    if d == 0
        y = double(x >= t(k) & x < t(k+1));
        if t(k+1) == max(t)
            y(x == t(k+1)) = 1;
        end
    else
        denom1 = t(k+d) - t(k);
        if denom1 == 0
            w1 = zeros(size(x));
        else
            w1 = (x - t(k)) / denom1;
        end

        denom2 = t(k+d+1) - t(k+1);
        if denom2 == 0
            w2 = zeros(size(x));
        else
            w2 = (t(k+d+1) - x) / denom2;
        end

        y = w1 .* bspline_recur(x, k, d-1, t) + ...
            w2 .* bspline_recur(x, k+1, d-1, t);
    end
end
