function [L_roi, R, roi_info, sens_idx] = build_roi_leadfield_hcp(a, varargin)
% 构造 19×360 的 ROI leadfield（HCP_MMP1）
% - 三向->法向定向
% - ROI 内按面积加权 + 与 seed 列对齐符号后求“有符号均值”
%
% Inputs:
%   a : load('.../surf.mat') 的结构，要求包含：
%       a.Sc.Vertices (8003×3), a.Sc.Faces (15994×3)
%       a.Sc.Atlas(k).Scouts(1×N) —— 选 HCP_MMP1 (N=360)
%       a.HeadModel.Gain (19×24009), a.HeadModel.GridOrient (8003×3)
% Options (name/value):
%   'atlas_idx'     (默认自动寻找 360 scouts 的atlas)
%   'use_area'      (默认 true) 顶点面积加权
%   'sign_align'    (默认 true) 与 seed 列对齐符号
%   'sens_idx'      (默认 1:19) 传感器保留索引（可在这里删一个电极）
%
% Outputs:
%   L_roi   : 19×360 的 leadfield
%   R       : 8003×360 的 ROI 归并矩阵（每列为权重，行和=1）
%   roi_info: struct，含 .scouts (每个 ROI 的顶点索引)、.seed、.labels
%   sens_idx: 实际使用的传感器索引

opts = struct('atlas_idx',[], 'use_area',true, 'sign_align',true, 'sens_idx', []);
if ~isempty(varargin), opts = setOpts(opts, varargin{:}); end

V  = a.Sc.Vertices;      % 8003×3
F  = a.Sc.Faces;         % 15994×3
G  = a.HeadModel.Gain;   % 19×(8003*3)
Nn = a.HeadModel.GridOrient; % 8003×3
p0 = size(G,1); nV = size(Nn,1);

if isempty(opts.sens_idx), sens_idx = 1:p0; else, sens_idx = opts.sens_idx(:).'; end
G = G(sens_idx, :);  % 选传感器
p = size(G,1);

% ---- 1) 选择 HCP_MMP1 atlas（scouts=360）----
atlas_idx = opts.atlas_idx;
if isempty(atlas_idx)
    atlas_idx = find(arrayfun(@(A) numel(A.Scouts)==360, a.Sc.Atlas), 1, 'first');
    if isempty(atlas_idx), error('找不到含 360 scouts 的 Atlas（HCP_MMP1）。'); end
end
Atlas = a.Sc.Atlas(atlas_idx);
if numel(Atlas.Scouts) ~= 360
    error('所选 Atlas 的 ROI 数不是 360。');
end

% ---- 2) 三向 -> 法向定向：得到 19×8003 的 L_norm ----
L_norm = zeros(p, nV);
for v = 1:nV
    Gv = G(:, 3*v-2 : 3*v);       % p×3
    nv = Nn(v, :).';              % 3×1
    L_norm(:, v) = Gv * nv;       % p×1
end

% ---- 3) 顶点面积（用于 ROI 加权）----
if opts.use_area
    Avert = vertex_area(V, F);  % 8003×1
else
    Avert = ones(nV,1);
end

% ---- 4) 构造 ROI 加权矩阵 R（有符号均值）----
R = spalloc(nV, 360, round(nV*1.2));   % 预估稀疏非零个数
roi_info = struct('scouts',{cell(360,1)}, 'seed',zeros(360,1), 'labels',{cell(360,1)});

for r = 1:360
    idx = Atlas.Scouts(r).Vertices(:);
    if isempty(idx), error('ROI %d 为空。', r); end
    roi_info.scouts{r} = idx;
    roi_info.labels{r} = Atlas.Scouts(r).Label;
    seed = Atlas.Scouts(r).Seed;
    if isempty(seed) || seed<=0, seed = idx(1); end
    roi_info.seed(r) = seed;

    % 面积权重 + 归一化
    w = Avert(idx);
    w = w / (sum(w) + eps);

    % 与 seed 列对齐符号，避免抵消
    sgn = ones(numel(idx),1);
    if opts.sign_align
        l0 = L_norm(:, seed);
        for k = 1:numel(idx)
            lv = L_norm(:, idx(k));
            s = real(lv' * l0);
            sgn(k) = 1; if s < 0, sgn(k) = -1; end
        end
    end

    % 写入 R
    R(idx, r) = w .* sgn;
end

% 列归一：确保每列权重和=1（考虑符号后取绝对值总和）
for r = 1:360
    col = full(R(:,r));
    s = sum(abs(col)) + eps;
    R(:,r) = col / s;
end

% ---- 5) 得到 19×360 的 ROI leadfield ----
L_roi = L_norm * R;  % p×360

end

% ===== helpers =====
function A = vertex_area(V, F)
% 给每个三角形的面积 1/3 分摊到其三个顶点
A = zeros(size(V,1),1);
for t = 1:size(F,1)
    v = F(t,:);
    e1 = V(v(2),:) - V(v(1),:);
    e2 = V(v(3),:) - V(v(1),:);
    Atri = 0.5 * norm(cross(e1, e2));
    A(v) = A(v) + Atri/3;
end
end

function o = setOpts(o, varargin)
for k=1:2:numel(varargin)
    o.(varargin{k}) = varargin{k+1};
end
end
