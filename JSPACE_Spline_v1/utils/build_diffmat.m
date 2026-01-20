function D = build_diffmat(K, order)
%BUILD_DIFFMAT  Construct finite difference matrix for P-splines.
%
% Inputs:
%   K     : number of spline coefficients
%   order : difference order (typically 1 or 2)
%
% Output:
%   D     : (K-order) x K difference matrix

    if nargin < 2 || isempty(order)
        order = 2;
    end
    if K <= 0 || order <= 0
        error('build_diffmat:InvalidConfig', 'K and order must be positive.');
    end
    if K <= order
        error('build_diffmat:InvalidConfig', ...
            'K must be larger than the difference order.');
    end

    I = eye(K);
    D = diff(I, order);
end
