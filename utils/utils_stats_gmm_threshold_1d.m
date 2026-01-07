function [t_star, meta] = utils_stats_gmm_threshold_1d(values)
%UTILS_STATS_GMM_THRESHOLD_1D GMM-EM threshold from 1D samples with fallback.
%
% Implements a 2-component Gaussian mixture, then finds the intersection of
% the two unweighted densities. If EM fails or no root is found, falls back
% to the 90th percentile.

    if nargin < 1
        error('utils_stats_gmm_threshold_1d:missing_input', 'values are required');
    end

    v = values(:);
    v = v(isfinite(v));
    v = v(abs(v) >= 1e-10);

    meta = struct();
    meta.fallback = false;
    meta.reason = '';
    meta.n_samples = numel(v);
    meta.loglik = -Inf;
    meta.iterations = 0;
    meta.params = struct();

    if numel(v) < 5
        [t_star, meta] = fallback_quantile_(v, meta, 'too_few_samples');
        return;
    end

    n_restarts = 20;
    max_iter = 1000;
    tol = 1e-3;
    delta = 1e-3;

    stream = RandStream('mt19937ar', 'Seed', 0);

    best = struct('pi', [], 'mu', [], 'sigma', [], 'loglik', -Inf, 'iter', 0);
    n = numel(v);
    v_std = std(v);
    if ~isfinite(v_std) || v_std <= 0
        v_std = 1;
    end

    for r = 1:n_restarts
        idx = random_pair_(stream, n);
        mu = v(idx);
        if mu(1) == mu(2)
            mu(2) = mu(2) + 0.1 * v_std;
        end
        sigma = [v_std, v_std];
        pi_k = [0.5, 0.5];

        loglik_prev = -Inf;
        loglik = -Inf;
        for iter = 1:max_iter
            p1 = pi_k(1) * normal_pdf_(v, mu(1), sigma(1));
            p2 = pi_k(2) * normal_pdf_(v, mu(2), sigma(2));
            denom = p1 + p2 + realmin;

            r1 = p1 ./ denom;
            r2 = p2 ./ denom;
            n1 = sum(r1);
            n2 = sum(r2);

            pi_k(1) = n1 / n;
            pi_k(2) = n2 / n;
            mu(1) = sum(r1 .* v) / max(n1, eps);
            mu(2) = sum(r2 .* v) / max(n2, eps);
            sigma(1) = sqrt(sum(r1 .* (v - mu(1)).^2) / max(n1, eps) + delta);
            sigma(2) = sqrt(sum(r2 .* (v - mu(2)).^2) / max(n2, eps) + delta);

            p1n = pi_k(1) * normal_pdf_(v, mu(1), sigma(1));
            p2n = pi_k(2) * normal_pdf_(v, mu(2), sigma(2));
            denom_n = p1n + p2n + realmin;
            loglik = sum(log(denom_n));
            if iter > 1
                rel_impr = abs((loglik - loglik_prev) / max(abs(loglik_prev), 1));
                if rel_impr < tol
                    break;
                end
            end
            loglik_prev = loglik;
        end

        if isfinite(loglik) && loglik > best.loglik
            best.pi = pi_k;
            best.mu = mu;
            best.sigma = sigma;
            best.loglik = loglik;
            best.iter = iter;
        end
    end

    if ~isfinite(best.loglik) || any(~isfinite(best.mu)) || any(~isfinite(best.sigma))
        [t_star, meta] = fallback_quantile_(v, meta, 'em_failed');
        return;
    end

    mu1 = best.mu(1); mu2 = best.mu(2);
    s1 = best.sigma(1); s2 = best.sigma(2);
    pi1 = best.pi(1); pi2 = best.pi(2);
    if mu1 > mu2
        [mu1, mu2] = deal(mu2, mu1);
        [s1, s2] = deal(s2, s1);
        [pi1, pi2] = deal(pi2, pi1);
    end

    [t_star, ok] = solve_intersection_(mu1, s1, mu2, s2);
    if ~ok || ~isfinite(t_star)
        [t_star, meta] = fallback_quantile_(v, meta, 'no_intersection');
        return;
    end

    meta.loglik = best.loglik;
    meta.iterations = best.iter;
    meta.params = struct('pi', [pi1, pi2], 'mu', [mu1, mu2], 'sigma', [s1, s2]);
end

% ----------------------- Local helpers -----------------------
function y = normal_pdf_(x, mu, sigma)
    sigma = max(sigma, sqrt(realmin));
    z = (x - mu) ./ sigma;
    y = exp(-0.5 * (z.^2)) ./ (sqrt(2*pi) * sigma);
end

function idx = random_pair_(stream, n)
    r = rand(stream, n, 1);
    [~, order] = sort(r);
    idx = order(1:2);
end

function [t_star, ok] = solve_intersection_(mu1, s1, mu2, s2)
    ok = false;
    if ~isfinite(mu1) || ~isfinite(mu2) || s1 <= 0 || s2 <= 0
        t_star = NaN;
        return;
    end
    if mu1 == mu2
        t_star = mu1;
        ok = true;
        return;
    end

    f = @(t) ((t - mu1).^2) ./ (2 * s1^2) - ((t - mu2).^2) ./ (2 * s2^2) - log(s2 / s1);

    lo = mu1;
    hi = mu2;
    f_lo = f(lo);
    f_hi = f(hi);

    width = hi - lo;
    for k = 1:20
        if sign(f_lo) ~= sign(f_hi)
            ok = true;
            break;
        end
        lo = lo - width;
        hi = hi + width;
        f_lo = f(lo);
        f_hi = f(hi);
    end
    if ~ok
        t_star = NaN;
        return;
    end

    for k = 1:80
        mid = 0.5 * (lo + hi);
        f_mid = f(mid);
        if sign(f_mid) == sign(f_lo)
            lo = mid;
            f_lo = f_mid;
        else
            hi = mid;
            f_hi = f_mid;
        end
        if abs(hi - lo) < 1e-12 * max(1, abs(mid))
            break;
        end
    end
    t_star = 0.5 * (lo + hi);
end

function [t_star, meta] = fallback_quantile_(v, meta, reason)
    meta.fallback = true;
    meta.reason = reason;
    if isempty(v)
        t_star = 0;
    else
        t_star = quantile(v, 0.9);
    end
end
