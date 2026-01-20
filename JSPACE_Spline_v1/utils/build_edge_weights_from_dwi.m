function w_edge = build_edge_weights_from_dwi(dwi_C, dwi_cfg, w_cfg)
%BUILD_EDGE_WEIGHTS_FROM_DWI  Compute group-lasso edge weights from DWI.
%
% Logic: weight ~ 1 / (connectivity^alpha)
%
% Inputs:
%   dwi_C   : (N x N) connectivity matrix (non-negative)
%   dwi_cfg : struct with optional fields:
%       .alpha : exponent (default 1.0)
%       .clip  : [min, max] clipping range (default [0.1, 10])
%       .eps0  : small constant to avoid divide-by-zero (default 1e-5)
%   w_cfg   : optional struct with .normalize_weight_median (logical)
%
% Output:
%   w_edge  : (N x N) symmetric edge weight matrix with zero diagonal

    if nargin < 2 || isempty(dwi_cfg)
        dwi_cfg = struct();
    end

    if ~isfield(dwi_cfg, 'alpha') || isempty(dwi_cfg.alpha)
        dwi_cfg.alpha = 1.0;
    end
    if ~isfield(dwi_cfg, 'clip') || isempty(dwi_cfg.clip)
        dwi_cfg.clip = [0.1, 10];
    end
    if ~isfield(dwi_cfg, 'eps0') || isempty(dwi_cfg.eps0)
        dwi_cfg.eps0 = 1e-5;
    end

    if ndims(dwi_C) ~= 2 || size(dwi_C, 1) ~= size(dwi_C, 2)
        error('build_edge_weights_from_dwi:SizeMismatch', ...
            'dwi_C must be square.');
    end

    N = size(dwi_C, 1);

    C = abs(dwi_C);
    C = 0.5 * (C + C');
    C(1:N+1:end) = 0;

    vals = C(C > 0);
    if isempty(vals)
        med_val = 1;
    else
        med_val = median(vals);
        if med_val <= 0 || ~isfinite(med_val)
            med_val = 1;
        end
    end
    C_norm = C / med_val;

    w = (C_norm + dwi_cfg.eps0) .^ (-dwi_cfg.alpha);
    w(~isfinite(w)) = dwi_cfg.clip(2);

    w = max(w, dwi_cfg.clip(1));
    w = min(w, dwi_cfg.clip(2));

    if nargin >= 3 && isstruct(w_cfg) && ...
            isfield(w_cfg, 'normalize_weight_median') && w_cfg.normalize_weight_median
        w_vals = w(w < dwi_cfg.clip(2));
        if ~isempty(w_vals)
            w = w / median(w_vals);
        end
    end

    w(1:N+1:end) = 0;
    w_edge = w;
end
