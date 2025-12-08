function Gamma_refit = module_refit_quadratic(Gamma_debiased, Var_proxies, mask_cell, L_freq, lambda1, lambda3, lambda4, n_samples, W_spatial)
% MODULE_REFIT_QUADRATIC - Edge-wise quadratic refit with cross-frequency smoothing
%
% Inputs:
%   Gamma_debiased : {F x 1} debiased dense matrices (G_tilde for each freq)
%   Var_proxies    : {F x 1} variance proxy matrices
%   mask_cell      : {F x 1} logical support masks for each frequency
%   L_freq         : (F x F) frequency Laplacian
%   lambda1        : scalar, cross-frequency smoothing weight
%   lambda3        : scalar, spatial L2 weight
%   lambda4        : scalar, optional ridge weight
%   n_samples      : scalar, sample size T
%   W_spatial      : (p x p) spatial weight matrix; if empty, weight defaults to 1
%
% Output:
%   Gamma_refit    : {F x 1} refitted SPD precision matrices

    F = numel(Gamma_debiased);
    p = size(Gamma_debiased{1}, 1);

    if nargin < 9, W_spatial = []; end
    if nargin < 8 || isempty(n_samples), n_samples = 1; end
    if nargin < 7 || isempty(lambda4), lambda4 = 0; end

    % === 1) Pre-pack data into 3D arrays，避免深层 cell 索引 ===
    G_all = zeros(p, p, F, 'like', Gamma_debiased{1});
    V_all = zeros(p, p, F, 'like', Var_proxies{1});
    M_all = false(p, p, F);
    for f = 1:F
        G_all(:, :, f) = Gamma_debiased{f};
        V_all(:, :, f) = Var_proxies{f};
        M_all(:, :, f) = mask_cell{f};
    end

    % union 支撑，只在至少一个频率非零的边上做 refit
    union_mask = any(M_all, 3);          % p x p
    union_mask = triu(union_mask, 1);    % 只保留上三角
    [edge_i, edge_j] = find(union_mask); % 所有需要 refit 的边
    n_edges = numel(edge_i);

    % 如果没有边，直接返回原 debiased
    Gamma_refit = Gamma_debiased;
    if n_edges == 0
        return;
    end

    % === 2) Eigendecomposition of frequency Laplacian (do once) ===
    if issparse(L_freq)
        L_freq_full = full(L_freq);
    else
        L_freq_full = L_freq;
    end
    [Ufreq, Dfreq] = eig(L_freq_full);
    lambda_vec = real(diag(Dfreq));   % F x 1 eigenvalues

    % === 3) 对每条边做近似解析 refit ===
    for e = 1:n_edges
        i = edge_i(e);
        j = edge_j(e);

        % 取出这一条边在所有频率上的 debiased 值 & 方差 proxy
        gamma0 = squeeze(G_all(i, j, :));          % F x 1
        v_vec  = squeeze(V_all(i, j, :));          % F x 1
        v_vec  = real(v_vec);
        h_vec  = n_samples ./ max(v_vec, 1e-12);   % F x 1

        % 空间权重
        if ~isempty(W_spatial)
            w_ij = real(W_spatial(i, j));
        else
            w_ij = 1.0;
        end

        % === Approximate H_ij by scalar h_bar * I for speed ===
        h_bar = mean(h_vec);

        % Combine L2 + ridge coefficients
        c_ij = h_bar + 2*lambda3*(w_ij^2) + 2*lambda4;

        % Denominator in frequency eigen space: den_k = c_ij + 2*lambda1*lambda_k
        den = c_ij + 2*lambda1 * lambda_vec;
        den = max(den, 1e-12); % avoid divide-by-zero

        % Spectral filtering of gamma0: gamma_star = h_bar * U * diag(1./den) * U' * gamma0
        tmp = Ufreq' * gamma0;             % F x 1
        tmp = (h_bar ./ den) .* tmp;       % element-wise scale
        gamma_star = Ufreq * tmp;          % F x 1

        % 写回到所有频率上（只更新 mask=1 的频点）
        for f = 1:F
            if M_all(i, j, f)
                G_all(i, j, f) = gamma_star(f);
                G_all(j, i, f) = conj(gamma_star(f));
            end
        end
    end

    % === 4) SPD projection per frequency ===
    for f = 1:F
        Gf = G_all(:, :, f);
        Gf = (Gf + Gf') / 2;
        % 强烈建议直接用 project_spd，不再用旧版 regularize_spd
        [Gf_spd, ~] = utils_math.project_spd(Gf, 1e-8);
        Gamma_refit{f} = Gf_spd;
    end
end
