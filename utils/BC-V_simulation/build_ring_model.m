function ring = build_ring_model(p, q, Re, Rs, beta, save_path)
% 构造：p个电极在半径Re圆环上，q个源在半径Rs圆环上，L(i,j)=1/||e_i - s_j||^beta
% 输出并保存：L、伪皮层网格(cortex.Vertices/Faces)、电极坐标electrodes
% 单位建议用米；默认 beta=2（类1/r^2 衰减）

if nargin < 5 || isempty(beta), beta = 2; end

% 1) 圆环坐标
theta_e = linspace(0, 2*pi, p+1); theta_e(end) = [];
theta_s = linspace(0, 2*pi, q+1); theta_s(end) = [];
electrodes = [Re*cos(theta_e(:)), Re*sin(theta_e(:)), zeros(p,1)];
sources   = [Rs*cos(theta_s(:)), Rs*sin(theta_s(:)), zeros(q,1)];

% 2) 解析式 Lead Field（可替换为更物理的偶极模型，这里给教学友好的 1/r^beta）
L = zeros(p,q);
for i = 1:p
    diffs = sources - electrodes(i,:);               % q×3
    r     = sqrt(sum(diffs.^2, 2));                  % q×1
    L(i,:)= (r(:)'.^(-beta));                        % 1×q
end
% 列归一化（数值稳健、与后续噪声归一一致）
L = L ./ max(vecnorm(L,2,1), 1e-12);

% 3) 伪皮层网格（把源点在二维圆上做三角剖分，嵌入 z=0 平面）
DT = delaunayTriangulation(sources(:,1), sources(:,2));
Faces = DT.ConnectivityList;
Vertices = sources;                                  % q×3，z=0

cortex.Vertices = Vertices;
cortex.Faces    = Faces;

% 4) 打包与保存
ring.L = L;
ring.cortex = cortex;
ring.electrodes = electrodes;
ring.meta = struct('p',p,'q',q,'Re',Re,'Rs',Rs,'beta',beta);

if nargin >= 6 && ~isempty(save_path)
    save(save_path, 'L', 'cortex', 'electrodes', 'p', 'q', 'Re', 'Rs', 'beta');
end
end
