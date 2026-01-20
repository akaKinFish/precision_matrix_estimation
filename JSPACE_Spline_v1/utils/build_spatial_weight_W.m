function W_spatial = build_spatial_weight_W(dwi_C, spatial_cfg)
%BUILD_SPATIAL_WEIGHT_W  Construct spatial L2 penalty weights.
%
% Logic: weaker connectivity -> larger penalty.
%
% Inputs:
%   dwi_C       : (N x N) connectivity matrix
%   spatial_cfg : struct with fields:
%       .enable : true/false
%       .clip   : [min, max] (optional)
%
% Output:
%   W_spatial   : (N x N) weight matrix

    if nargin < 2 || isempty(spatial_cfg) || ...
            ~isfield(spatial_cfg, 'enable') || ~spatial_cfg.enable
        N = size(dwi_C, 1);
        W_spatial = zeros(N);
        return;
    end

    if ndims(dwi_C) ~= 2 || size(dwi_C, 1) ~= size(dwi_C, 2)
        error('build_spatial_weight_W:SizeMismatch', ...
            'dwi_C must be square.');
    end

    clip_range = [0.1, 10];
    if isfield(spatial_cfg, 'clip') && ~isempty(spatial_cfg.clip)
        clip_range = spatial_cfg.clip;
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
    C = C / med_val;

    eps0 = 1e-5;
    W = 1 ./ (C + eps0);
    W(~isfinite(W)) = clip_range(2);

    W = max(W, clip_range(1));
    W = min(W, clip_range(2));

    W(1:N+1:end) = 0;
    W_spatial = W;
end
