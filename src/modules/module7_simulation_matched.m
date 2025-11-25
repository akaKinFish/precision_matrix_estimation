function [Omega_true, Sigma_true, emp_cov_cell, L, T, match_truth, sim] = module7_simulation_matched(cfg)
% MODULE7_SIMULATION_MATCHED (Route B)
% 显式使用 λ1, λ2 生成频域稀疏+平滑的 precision 真值：
%   对每条边 e 的频域轨迹 g_e ∈ R^F，求解（PGD 近似MAP）：
%       min_g  0.5*||g - s*Z||_2^2 + (λ2/2)*g^T Lk g + λ1*||g||_2
%   其中 Z~N(0, I_F)，s 通过二分调节来达到目标稀疏度 edge_density。
%
% 输出：
%   Omega_true : p×p×F   源域 precision 真值（SPD）
%   Sigma_true : p×p×F   源域 covariance 真值
%   emp_cov_cell : {F×1} 传感器域经验协方差
%   L          : n×p 前向矩阵（列单位范数）
%   T          : 每频样本数
%   match_truth: 结构体（lambda1_true, lambda2_true, ... 以及生成用到的G等）
%   sim        : 元信息

%% 0) config 默认
cfg = set_default(cfg, struct( ...
    'p', 10, 'n', 3, 'F', 3, 'T', 4096, 'random_seed', 42, ...
    'edge_density', 0.15, ...      % 目标活跃组比例（||g_e||>0 的边比例）
    'lambda1_star', 0.35, ...      % 组稀疏 λ1
    'lambda2_star', 1.0,  ...      % 频域平滑 λ2（乘在 g^T Lk g）
    'laplacian_type', 'chain', ... % 'chain'|'cycle'|'custom'
    'laplacian_custom', [], ...
    'diag_base', 0.2, ...          % 对角基础值（确保对角占优）
    'diag_spd', 1e-2, ...          % SPD 抬角
    'complex_samples', true, ...
    'noise_type', 'scalar', ...    % 'scalar'|'matrix'|'matrix_per_freq'
    'sigma_xi2', 0.2, ...
    'Sigma_xixi', [], ...
    'Sigma_xixi_per_freq', [] ,...
    'L_fixed', [] ...   % 可选，若提供则使用作为固定 lead field
));
rng(cfg.random_seed);
p = cfg.p; n = cfg.n; F = cfg.F; T = cfg.T;
E = p*(p-1)/2;

%% 1) 频域拉普拉斯 Lk
Lk = make_freq_laplacian(F, cfg.laplacian_type, cfg.laplacian_custom);
% Lipschitz 常数：L = 1 + λ2*λ_max(Lk)
lam_max = max(real(eig((Lk+Lk')/2)));
if ~isfinite(lam_max), lam_max = 0; end

%% 2) 基础模式 Z（每条边一个 F 维向量）
Z = randn(F, E);  % 你可以改为复数，但估计器是实对称，建议先实数

%% 3) 二分搜索尺度 s，使得通过 PGD 得到的 g_e 有 ~edge_density 的活跃比例
target = cfg.edge_density;
lam1 = cfg.lambda1_star;
lam2 = cfg.lambda2_star;

% PGD 参数
Lips = 1 + lam2*lam_max;
eta  = 1 / max(Lips, 1e-12);
maxit = 500;
tol   = 1e-6;

% 二分区间（经验上 0.1~5 基本覆盖；必要可放宽）
s_lo = 0.05; s_hi = 10;
act_tol = 1e-3;  % 判定“非零”的范数阈值
for it=1:20
    s_mid = sqrt(s_lo*s_hi);   % 几何平均更稳定
    G = zeros(F,E);
    for e=1:E
        z = s_mid * Z(:,e);
        G(:,e) = pgd_group(z, Lk, lam1, lam2, eta, maxit, tol);
    end
    dens = nnz(vecnorm(G)>act_tol)/E;
    if dens > target
        % 太稠密 → 降低 s（等价于降低“信号幅度”）
        s_hi = s_mid;
    else
        % 太稀疏 → 提升 s
        s_lo = s_mid;
    end
end
s_final = sqrt(s_lo*s_hi);

% 最终 G
G = zeros(F,E);
for e=1:E
    z = s_final * Z(:,e);
    G(:,e) = pgd_group(z, Lk, lam1, lam2, eta, maxit, tol);
end

%% 4) 由 G 构造 Gamma_true / Omega_true
pairs = find(triu(true(p),1));
Gamma_true = zeros(p,p,F);
for f=1:F
    M = zeros(p);
    vals = -G(f,:);       % 负号，构造 M-型图（便于对角占优）
    M(pairs) = vals;
    M = M + M.';
    diagv = cfg.diag_base + sum(abs(M),2) + cfg.diag_spd;
    M(1:p+1:end) = diagv;
    Gamma_true(:,:,f) = force_spd(M, 1e-10);
end
Omega_true = Gamma_true;

%% 5) 源域协方差 + 传感器观测
Sigma_true = zeros(p,p,F);
for f=1:F
    Ki = (Omega_true(:,:,f)+Omega_true(:,:,f)')/2;
    Sigma_true(:,:,f) = hermitize(inv(Ki));
end

% 前向矩阵 L（列单位范数）
% 前向矩阵 L（列单位范数；若传入 L_fixed 则使用之）
if ~isempty(getd(cfg,'L_fixed',[]))
    L = cfg.L_fixed;
else
    L = randn(n,p);
end
L = L * diag(1./max(sqrt(sum(L.^2,1)),1e-12));  % 列单位范数
% 噪声
switch lower(cfg.noise_type)
    case 'scalar'
        Sig_xi = cfg.sigma_xi2 * eye(n);
        Sig_xi_cell = repmat({Sig_xi}, F,1);
    case 'matrix'
        Sig_xi = hermitize(cfg.Sigma_xixi);
        Sig_xi_cell = repmat({Sig_xi}, F,1);
    case 'matrix_per_freq'
        Sig_xi_cell = cell(F,1);
        for f=1:F, Sig_xi_cell{f} = hermitize(cfg.Sigma_xixi_per_freq{f}); end
    otherwise
        error('noise_type invalid');
end

% 采样 & 经验协方差
emp_cov_cell = cell(F,1);
for f=1:F
    Sf = Sigma_true(:,:,f);
    if cfg.complex_samples
        Ssrc = (randn(p,T)+1i*randn(p,T))/sqrt(2);
    else
        Ssrc = randn(p,T);
    end
    Ssrc = chol_psd(Sf) * Ssrc;
    if cfg.complex_samples
        Noise = (randn(n,T)+1i*randn(n,T))/sqrt(2);
    else
        Noise = randn(n,T);
    end
    Y = L*Ssrc + chol_psd(Sig_xi_cell{f}) * Noise;
    emp_cov_cell{f} = hermitize( (Y*Y')/T );
end

%% 6) 诊断与回传
support = false(p,p,F);
for f=1:F
    A = abs(triu(Gamma_true(:,:,f),1)) > 0;
    support(:,:,f) = A | A.';
end
[edges_per_freq, unique_edges, jacc_min, jacc_med, jacc_max] = support_stats(support);

match_truth = struct();
match_truth.lambda1_true = lam1;
match_truth.lambda2_true = lam2;
match_truth.noise_type   = cfg.noise_type;
match_truth.sigma_xi2_true = cfg.sigma_xi2;
match_truth.Sigma_xixi_true = getd(cfg,'Sigma_xixi',[]);
match_truth.Sigma_xixi_per_freq_true = getd(cfg,'Sigma_xixi_per_freq',[]);
match_truth.support      = support;
match_truth.edges_per_freq = edges_per_freq;
match_truth.unique_edges = unique_edges;
match_truth.jaccard_adj  = [jacc_min, jacc_med, jacc_max];
match_truth.Gamma_true   = Gamma_true;
match_truth.G_per_edge   = G;            % 生成用到的 g_e 轨迹（便于研究）
match_truth.s_final      = s_final;      % 为达密度匹配的尺度
match_truth.Lk           = Lk;

sim = struct('p',p,'m',n,'F',F,'T',T,'L',L);

fprintf('Sim[matched]: p=%d, n=%d, F=%d, T=%d | unique_edges=%d | Jaccard[min/med/max]=[%.2f/%.2f/%.2f]\n',...
    p,n,F,T, unique_edges, jacc_min, jacc_med, jacc_max);

end

%% ====== 内部函数 ======
function g = pgd_group(z, Lk, lam1, lam2, eta, maxit, tol)
% 解决：min_g  0.5*||g - z||^2 + (lam2/2)*g^T Lk g + lam1*||g||_2
% PGD:  v = g - η[(g - z) + lam2*Lk*g];  g <- soft_group(v, η*lam1)
g = z;  % 初始化也可用 zeros
for it=1:maxit
    grad = (g - z) + lam2 * (Lk * g);
    v = g - eta*grad;
    nv = norm(v);
    if nv <= eta*lam1
        g_new = zeros(size(g));
    else
        g_new = (1 - (eta*lam1)/max(nv,eps)) * v;
    end
    if norm(g_new - g) <= tol * max(1, norm(g))
        g = g_new; break;
    end
    g = g_new;
end
end

function Lk = make_freq_laplacian(F, typ, Lcustom)
switch lower(typ)
    case 'chain'
        e = ones(F,1);
        Lk = spdiags([ -e 2*e -e ], -1:1, F,F);
        Lk(1,1)=1; Lk(F,F)=1; Lk = full(Lk);
    case 'cycle'
        e = ones(F,1);
        Lk = spdiags([ -e 2*e -e ], -1:1, F,F);
        Lk(1,F)=-1; Lk(F,1)=-1; Lk = full(Lk);
    case 'custom'
        Lk = Lcustom;
    otherwise
        error('laplacian_type invalid');
end
Lk = (Lk+Lk')/2;
end

function [epf, uniq, jmin, jmed, jmax] = support_stats(supp)
[p,~,F] = size(supp);
mask = triu(true(p),1);
E = nnz(mask);
B = false(E,F);
for f=1:F
    A = supp(:,:,f);
    B(:,f) = A(mask);
end
epf = sum(B,1);
uniq = nnz(any(B,2));
if F>=2
    J = zeros(F-1,1);
    for f=1:F-1
        a = B(:,f); b = B(:,f+1);
        J(f) = nnz(a & b) / max(nnz(a | b), 1);
    end
    jmin=min(J); jmed=median(J); jmax=max(J);
else
    [jmin,jmed,jmax] = deal(NaN);
end
end

function C = chol_psd(S)
S = hermitize(S);
[V,D] = eig(S,'vector'); D = real(D); D(D<0)=0;
C = V*diag(sqrt(D))*V';
end

function M = force_spd(M, eps_floor)
M = hermitize(M);
[V,D] = eig(M,'vector'); D = real(D);
D(D<eps_floor)=eps_floor;
M = V*diag(D)*V'; M = hermitize(M);
end

function A = hermitize(A), A = (A + A')/2; end

function val = getd(s, name, def)
if isstruct(s) && isfield(s,name) && ~isempty(s.(name)), val = s.(name); else, val = def; end
end

function cfg = set_default(cfg, D)
fn = fieldnames(D);
for i=1:numel(fn)
    k = fn{i};
    if ~isfield(cfg,k) || isempty(cfg.(k)), cfg.(k)=D.(k); end
end
end
