function [Svv, V, V0, Sjj_sim, Thetajj_sim, J, patch_idx] = ...
    simulate_hggm_realistic(L, seed_idx, vertices, faces, m, options, noise_scale, d0)
% L: p×q; seed_idx: 活动源索引（列向量/下标）；vertices: q×3；faces: F×3
% m: 样本数；options: gen_hggm2 的配置（config/var/extensions/connections等）
% noise_scale: 噪声系数（默认0.1）；d0: 补丁半径（默认5 mm）

if nargin < 7 || isempty(noise_scale), noise_scale = 0.1; end
if nargin < 8 || isempty(d0), d0 = 5e-3; end

[p, q] = size(L);
q_seed = numel(seed_idx);

% 1) 生成源（只在 seed 上有活动）
[Sjj_sim, J_seed, Thetajj_sim] = gen_hggm2(m, q_seed, options);
J = zeros(q, m); J(seed_idx, :) = J_seed;

% 2) 纯净投影
V0 = L * J;

% 3) 汇总补丁索引
patch_idx = false(q,1);
for k = 1:q_seed
    s = seed_idx(k);
    idx_k = get_patch_idx(s, vertices, faces, d0);  % 支持 surfpatch 或 fallback
    patch_idx(idx_k) = true;
end
idx = find(patch_idx);

% 4) 生理噪声（只在补丁处）
bio = randn(numel(idx), m) + 1i*randn(numel(idx), m);
Vbio = L(:, idx) * bio;
Vbio = norm(V0,'fro') * Vbio / max(norm(Vbio,'fro'), 1e-12);

% 5) 传感器噪声
Vsens = randn(p,m) + 1i*randn(p,m);
Vsens = norm(V0,'fro') * Vsens / max(norm(Vsens,'fro'), 1e-12);

% 6) 合成与协方差
V   = V0 + noise_scale * Vbio + noise_scale * Vsens;
Svv = cov(V.');   % 注意转置

end

function idx = get_patch_idx(s, vertices, faces, d0)
% 优先用 surfpatch；若不存在则用欧氏半径阈值近邻作为近似
if exist('surfpatch','file') == 2
    [idx, ~] = surfpatch(s, vertices, faces, d0);
else
    % fallback：与第s个点的欧氏距离<=d0 的顶点
    d = sqrt(sum((vertices - vertices(s,:)).^2, 2));
    idx = find(d <= d0);
end
end
