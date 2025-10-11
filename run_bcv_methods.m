function out = run_bcv_methods(emp_cov_cell, L, T, opts)
% RUN_BCV_METHODS  在同一 Svv & L 上运行 BC-V 的方法
% 返回:
%   out.higgs / out.eloreta_hglasso / out.lcmv_hglasso / out.vareta
%   每个字段都是 {F×1}，元素为 n×n 源域 "precision" 矩阵（VARETA 我们做了稳健反演）

if nargin < 4, opts = struct(); end

F = numel(emp_cov_cell);
p = size(L,1);
q = size(L,2);

% ---- 公共 param 模版（来自 BC-V demo 的典型设置） ----
base = struct();
base.maxiter_outer = 60;
base.maxiter_inner = 30;

% 尺寸与单位阵
base.p   = p;
base.q   = q;
base.Ip  = eye(p);
base.Iq  = eye(q);

% 关键：传感器噪声初值向量 Op（BC-V 期望 p×1 向量）
base.Op  = ones(p,1);          % *** 修复点 ***

% 采样与正则相关
base.m   = T;
base.aj  = sqrt(log(q)/max(T,1));
base.Ajj = ones(q) - eye(q);
base.axi = 1e-4;
base.Axixi     = eye(p);
base.Axixi_inv = eye(p);

% 其它运行参数（与示例一致）
base.ntry          = 0;
base.prew          = 0;
base.nu            = T;
base.rth1          = 0.7;
base.rth2          = 3.16;
base.eigreg        = 1e-4;     % *** 建议补充 ***
base.run_bash_mode = 1;        % *** 建议补充（不弹waitbar）***
base.use_gpu       = 0;        % 可按需
base.str_band      = "none";   % 可按需
base.method = 'lqa';
base.run_bash_mode = 1;

out = struct();
out.higgs           = cell(F,1);
out.eloreta_hglasso = cell(F,1);
out.lcmv_hglasso    = cell(F,1);
out.vareta          = cell(F,1);

for f = 1:F
    Svv = (emp_cov_cell{f} + emp_cov_cell{f}')/2;
    Lvj = L;
    param = base;

    % === 1) H-HGGM（lasso）===
    param.penalty = 1; % 1=lasso, 2=ridge, 0=naive  :contentReference[oaicite:10]{index=10}
    [Theta_h, ~] = higgs(Svv, Lvj, param);                              % :contentReference[oaicite:11]{index=11}
    out.higgs{f} = symmetrize_(Theta_h);

    % === 2) eLORETA + HG-lasso ===
    param.gamma1      = 0.001;  % 也可暴露为 opts
    param.gamma2      = 0.05;
    param.delta_gamma = 0.001;
    [Theta_el, ~] = eloreta_hg_lasso(Svv, Lvj, param);                  % :contentReference[oaicite:12]{index=12}
    out.eloreta_hglasso{f} = symmetrize_(Theta_el);

    % === 3) LCMV + HG-lasso ===
    param.gamma = sum(abs(diag(Svv)))/(length(Svv)*100);                 % :contentReference[oaicite:13]{index=13}
    [Theta_lc, ~] = lcmv_hg_lasso(Svv, Lvj, param);
    out.lcmv_hglasso{f} = symmetrize_(Theta_lc);

    % === 4) VARETA ===
    [U,S,V] = svd(Lvj,'econ'); S = diag(S);                             % :contentReference[oaicite:14]{index=14}
    [Sjj, ~, ~, ~, ~] = vareta(U, S, V, Svv, 0);  % 返回的是 source covariance
    out.vareta{f} = robust_inv_psd_(Sjj);   % 转 Precision（稳健反演）
end

end

% === helpers ===
function A = symmetrize_(A), A = (A + A')/2; end
function Om = robust_inv_psd_(S)
    S = (S + S')/2;
    [V,D] = eig(S); d = real(diag(D));
    floorv = max(1e-10, 1e-10*max(d));
    d(d<floorv) = floorv;
    Om = V * diag(1./d) * V'; Om = (Om + Om')/2;
end
