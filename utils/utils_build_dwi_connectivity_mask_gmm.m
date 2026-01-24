function [mask, info] = utils_build_dwi_connectivity_mask_gmm(C_in, cfg)
%UTILS_BUILD_DWI_CONNECTIVITY_MASK_GMM Build sparse connectivity mask from DWI weights.

    if nargin < 2, cfg = struct(); end

    C = gather(C_in);
    C = real(C);
    C = max(C, 0);

    p = size(C, 1);
    C(1:p+1:end) = 0;

    use_log = get_cfg_(cfg, 'dwi_gmm_log', false);
    clip_q = get_cfg_(cfg, 'dwi_gmm_clip_quantiles', []);

    vals = C(tril(true(p), -1));
    vals = vals(isfinite(vals));

    if use_log
        vals = log1p(vals);
    end

    if ~isempty(clip_q) && numel(clip_q) == 2
        lo = quantile(vals, clip_q(1));
        hi = quantile(vals, clip_q(2));
        vals = vals(vals >= lo & vals <= hi);
    end

    [t_star, meta] = utils_stats_gmm_threshold_1d(vals);
    if use_log
        t_star = expm1(t_star);
    end

    if ~isfinite(t_star) || t_star < 0
        t_star = 0;
    end

    mask = C > t_star;
    mask = mask | mask.';
    mask(1:p+1:end) = true;

    info = struct();
    info.threshold = t_star;
    info.meta = meta;
    info.use_log = use_log;
    info.clip_quantiles = clip_q;
end

function val = get_cfg_(s, f, d)
    if isfield(s, f), val = s.(f); else, val = d; end
end
