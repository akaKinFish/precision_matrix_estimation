function Xreg = regularize_spd(X, tol)
% REGULARIZE_SPD Make matrix (numerically) SPD by shifting the spectrum.
%
% Xreg = regularize_spd(X, tol)
%   - symmetrize X
%   - compute min eigenvalue
%   - if min_eig < tol, add (tol - min_eig)*I

    if nargin < 2
        tol = 1e-8;
    end

    % 1) 对称化
    X = (X + X') / 2;

    % 2) 算最小特征值（可以用 eigs 的 'smallestreal' 加速）
    % 对尺寸不大时直接 eig 更简单
    e = eig(X);
    emin = min(real(e));

    % 3) 如已 SPD（或准 SPD），直接返回
    if emin >= tol
        Xreg = X;
        return;
    end

    % 4) 做一个 uniform shift：X + (tol - emin)*I
    shift = tol - emin;      % shift > 0
    Xreg = X + shift * eye(size(X), 'like', X);
    Xreg = (Xreg + Xreg')/2; % 再对称一下防数值误差
end
