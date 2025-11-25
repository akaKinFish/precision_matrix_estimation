function out = real_head_model_stimulation(paths, opts)
% REAL_HEAD_MODEL_STIMULATION
% 从真实头模生成仿真传感器数据 Svv，并返回可用于
% hglasso / 你的 PGD-EM 的输入与源域真值。
%
% 输入
%   paths.channel_mat   = '...\channel\channel.mat'      % 可选（目前仅为将来扩展）
%   paths.headmodel_mat = '...\leadfield\headmodel.mat'  % 必需（含 HeadModel.Gain, GridOrient/VertNormals）
%   paths.surf_mat      = '...\surf\surf.mat'            % 必需（含 Sc.Vertices, VertNormals, Atlas）
%
%   opts.mode   = 'roi' | 'fps'         % ROI 聚合（稳）或 FPS 采样（简洁）
%   opts.q      = 32                    % 目标源维度（建议 q ≤ ~2p）
%   opts.roi_ids = []                   % 指定 ROI 下标（仅 'roi' 模式；默认取前 q 个）
%   opts.m      = 2000                  % 样本数（等价频域"段数×峰数"）
%   opts.car    = true                  % 是否做平均参考 CAR（推荐）
%   opts.gamma_diagload = 0.05          % 对角加载系数（稳谱）
%   opts.bio_noise.nsrc = 50            % 生物噪声源数量
%   opts.bio_noise.scale = 0.10         % 生物噪声相对信号的比例（Fro 范数）
%   opts.sens_noise.scale = 0.10        % 传感器噪声相对信号的比例
%   opts.truth.block_num = 4            % 传给 gen_hggm2 的真值网络参数（示例）
%   opts.truth.other_fields ...         % 其余 gen_hggm2 需要的选项
%
% 输出（结构 out）
%   out.Svv        [p x p]   传感器样本协
%   out.v          [p x m]   传感器时序（零均值）
%   out.Lvj        [p x q]   最终用于反演的导联矩阵（ROI 或 FPS 后）
%   out.Kproj      [p x V]   法向投影后的全网格导联（便于复用/显示）
%   out.H          [p x p]   CAR 参考矩阵（若 opts.car=false，则为 I）
%   out.R          [V x q]   ROI 聚合映射（'roi' 模式下有效；'fps' 返回稀疏单位列）
%   out.indms      [1 x q]   ROI 索引或 FPS 选中顶点索引
%   out.Theta_true [q x q]   源域真值精度
%   out.Sjj_true   [q x q]   源域真值协方差
%   out.j_sim      [q x m]   源域样本
%   out.meta       struct    若干元信息（p, V, q, SNR 参数等）
%
% 依赖：gen_hggm2（你仓库已有），若无可改为 gen_hggm1 或自定义。

% ---------- 0) 读取 ----------
if ~isfield(paths, 'headmodel_mat') || ~isfile(paths.headmodel_mat)
    error('headmodel_mat 不存在：%s', paths.headmodel_mat);
end
if ~isfield(paths, 'surf_mat') || ~isfile(paths.surf_mat)
    error('surf_mat 不存在：%s', paths.surf_mat);
end
hm   = load(paths.headmodel_mat);   % -> hm.HeadModel.Gain, .GridOrient (可选)
surf = load(paths.surf_mat);        % -> surf.Sc

% ---------- 1) 三分量→法向一分量（稳条件数） ----------
K3D = hm.HeadModel.Gain;            % [p x (3V)]
[p, n3] = size(K3D);
V = n3/3;
if abs(V - round(V)) > 0
    error('导联矩阵列数不是 3 的倍数。');
end
% 获取法向：优先 Sc.VertNormals，其次 HeadModel.GridOrient
if isfield(surf, 'Sc') && isfield(surf.Sc, 'VertNormals') && ~isempty(surf.Sc.VertNormals)
    nrm = double(surf.Sc.VertNormals);    % [V x 3]
elseif isfield(hm.HeadModel, 'GridOrient') && ~isempty(hm.HeadModel.GridOrient)
    nrm = double(hm.HeadModel.GridOrient); % [V x 3]
else
    error('未找到顶点法向（VertNormals 或 GridOrient）。');
end
Kproj = zeros(p, V);
for ii = 1:p
    Ki = reshape(K3D(ii,:), 3, []).';     % [V x 3]
    Kproj(ii,:) = sum(Ki .* nrm, 2);      % 法向投影
end

% ---------- 2) 平均参考（CAR） ----------
if ~isfield(opts, 'car'); opts.car = true; end
if opts.car
    H = eye(p) - (1/p) * ones(p);
    Kproj = H * Kproj;
else
    H = eye(p);
end

% ---------- 3) ROI 聚合 或 FPS 采样，得到最终 q ----------
if ~isfield(opts, 'mode'); opts.mode = 'roi'; end
if ~isfield(opts, 'q') || isempty(opts.q); opts.q = min(2*p, max(8, floor(p*1.5))); end
q = opts.q;

R = sparse(V, q);
indms = [];
switch lower(opts.mode)
    case 'roi'
        if ~isfield(surf.Sc, 'Atlas') || isempty(surf.Sc.Atlas)
            error('没有 Atlas，无法进行 ROI 聚合。请改用 FPS 模式或提供 Atlas。');
        end
        atlas = surf.Sc.Atlas(surf.Sc.iAtlas);
        scouts = atlas.Scouts;
        roi_all = 1:numel(scouts);
        if isfield(opts, 'roi_ids') && ~isempty(opts.roi_ids)
            roi_sel = opts.roi_ids(:).';
        else
            roi_sel = roi_all(1:min(q, numel(roi_all)));
        end
        q = numel(roi_sel);
        R = sparse(V, q);
        for j = 1:q
            idx = scouts(roi_sel(j)).Vertices(:);
            w = ones(numel(idx),1);     % 可替换为面积权重
            w = w / sum(w + eps);
            R(idx, j) = w;
        end
        Lvj = Kproj * R;                % [p x q]
        indms = roi_sel;

    case 'fps'
        if q > V, error('q(%d) 不能大于 V(%d)。', q, V); end
        if ~isfield(opts, 'fps_seed'); opts.fps_seed = 1; end
        rng(opts.fps_seed);
        P = double(surf.Sc.Vertices);   % [V x 3]
        sel = farthest_point_sampling(P, q);
        indms = sel(:).';
        R = sparse(V, q);
        for j = 1:q, R(indms(j), j) = 1; end
        Lvj = Kproj(:, indms);          % [p x q]

    otherwise
        error('未知 opts.mode：%s（应为 roi 或 fps）', opts.mode);
end

% ---------- 4) 源域真值网络与样本 ----------
if ~isfield(opts, 'm') || isempty(opts.m); opts.m = 2000; end
m = opts.m;
if ~isfield(opts, 'truth'); opts.truth = struct(); end
if ~isfield(opts.truth, 'block_num'); opts.truth.block_num = 4; end

% 统一准备 gen_hggm2 的参数（若调用者没给，就按旧口径设）
gopt = struct('config', 2, 'var', 2, ...
              'extensions', [ceil(q/3); ceil(q/3); q-2*ceil(q/3)], ...
              'connections', [1 2; 2 3]);

% 你仓库里的生成器：gen_hggm2(q,m,options) / gen_hggm1
try
    [~, j_sim, Theta_true] = gen_hggm2(m, q, gopt);
catch
    % 兜底：若无 gen_hggm2，用简单 Toeplitz SPD 生成器
    Theta_true = gallery('tridiag', q, -0.3, 1, -0.3);
    j_sim = (chol(inv(Theta_true), 'lower')) * randn(q, m);  % ~ N(0,Sjj_true)
end
Sjj_true = inv((Theta_true + Theta_true')/2);

% ---------- 5) 前向投影 + 生物噪声 + 传感器噪声 ----------
v_sig = Lvj * j_sim;                       % [p x m] 信号
% 生物噪声
if ~isfield(opts, 'bio_noise'); opts.bio_noise = struct(); end
if ~isfield(opts.bio_noise, 'nsrc');  opts.bio_noise.nsrc  = 50; end
if ~isfield(opts.bio_noise, 'scale'); opts.bio_noise.scale = 0.10; end
pool = setdiff(1:V, find(any(R,2)).');     % 非 ROI/FPS 顶点作为生物噪声候选
nbio = min(opts.bio_noise.nsrc, numel(pool));
bio_id = pool(randperm(numel(pool), nbio));
K_bio = Kproj(:, bio_id);
J_bio = randn(nbio, m);
v_bio = K_bio * J_bio;
% 归一化到与 v_sig 同量级后按比例缩放
nf_sig = norm(v_sig, 'fro') + eps;
v_bio = v_bio / (norm(v_bio, 'fro') + eps) * nf_sig * opts.bio_noise.scale;

% 传感器噪声
if ~isfield(opts, 'sens_noise'); opts.sens_noise = struct(); end
if ~isfield(opts.sens_noise, 'scale'); opts.sens_noise.scale = 0.10; end
v_n = randn(p, m);
v_n = v_n / (norm(v_n, 'fro') + eps) * nf_sig * opts.sens_noise.scale;

% 合成 & CAR 到同一参考
v = v_sig + v_bio + v_n;
if opts.car, v = H * v; end

% 样本协 & 稳谱
Svv = (v * v.') / m;
Svv = (Svv + Svv')/2;
if ~isfield(opts, 'gamma_diagload') || isempty(opts.gamma_diagload)
    opts.gamma_diagload = 0.05;
end
gamma = opts.gamma_diagload;
Svv = (1-gamma)*Svv + gamma*(trace(Svv)/p)*eye(p);

% ---------- 6) 打包输出 ----------
out.Svv      = Svv;
out.v        = v;
out.Lvj      = Lvj;
out.Kproj    = Kproj;
out.H        = H;
out.R        = R;
out.indms    = indms;
out.Theta_true = Theta_true;
out.Sjj_true   = Sjj_true;
out.j_sim      = j_sim;
out.meta = struct('p',p,'V',V,'q',q, ...
    'mode',opts.mode,'m',m, ...
    'bio_scale',opts.bio_noise.scale, ...
    'sens_scale',opts.sens_noise.scale, ...
    'gamma_diagload',gamma);
end

% ====== 辅助函数：最远点采样（FPS）======
function sel = farthest_point_sampling(P, k)
% P: [V x 3] 顶点坐标；k: 选择个数
V = size(P,1);
sel = zeros(1,k);
% 初始选一个几何中心最远点
c = mean(P,1);
d = sum((P - c).^2, 2);
[~, sel(1)] = max(d);
dist = sum((P - P(sel(1),:)).^2, 2);
for i = 2:k
    [~, sel(i)] = max(dist);
    dist = min(dist, sum((P - P(sel(i),:)).^2, 2));
end
end
