function out = compute_spd_and_classification_metrics(args, opts)
% COMPUTE_SPD_AND_CLASSIFICATION_METRICS — 与 Ground Truth 的距离 + 分类指标（可选雷达图）
%
% 该函数独立可复用；输入三套估计矩阵（估计/去偏/雷利矫正）与真值，输出：
%   1) 六种 SPD 距离：Kullback (Stein), Log-Euclid, Riemann (AIRM), JBLD, Bures-Wasserstein (OT2), Alpha-divergence
%   2) ROC/AUC 以及在“最优工作点”的 SENS, SPEC, PREC, F1（与 BC-V/论文 Fig.7 口径一致）
%   3) 可选绘制雷达图（AUC/SENS/SPEC/PREC/F1）
%
% ------------------
% Inputs (struct)
% ------------------
% args.Theta_true : 真值精度矩阵 (p×p) 或含 .X 字段的结构体
% args.Theta_est  : 估计精度矩阵 (p×p) 或含 .X 字段的结构体（可选）
% args.Theta_unb  : 去偏后精度矩阵（可选）
% args.Theta_ray  : Rayleigh 矫正后精度矩阵（可选）
% args.mask       : 逻辑 p×p，评估/阈值/得分所用的边集，
%                   缺省=上三角去对角（triu(true(p),1)）
% args.is_complex : logical，输入是否可能为复 Hermitian（默认 true）
% args.sample_m   : 用于日志/可视化的样本数（可选，统计阈值内部未使用）
%
% opts.normalize_mode : 'maxabs'|'none'，得分归一化（默认 'maxabs'）
% opts.spd.symmetrize : true/false，是否强制 (X+X')/2（默认 true）
% opts.spd.eps        : SPD 投影的最小特征值下界倍数（默认 1e-10）
% opts.spd.project    : true/false，必要时做 SPD 投影（默认 true）
% opts.alpha.value    : α-divergence 的 α（默认 0）
% opts.alpha.use_bcv_bug : true 时按 BC-V 源码中的分母笔误实现（默认 false）
% opts.plot.radar     : logical 是否画雷达图（默认 false）
% opts.plot.title     : 字符串 雷达图标题（可选）
%
% ------------------
% Outputs (struct)
% ------------------
% out.distances.est / .unb / .ray : 每类对真值的六种距离（字段：kullback, logeuclid, riemann, ld, opttransp, alphadiv）
% out.metrics.est / .unb / .ray   : 每类的分类指标（auc, sens, spec, prec, f1, thr_opt, fpr, tpr）
% out.fig.radar                   : 雷达图句柄（若绘制）
%
% 备注：
% - 本函数不会改写输入矩阵；必要时仅在内部做对称化与 SPD 投影（拷贝）。
% - 分类指标依据：
%   真值标签：|Theta_true| 的非零上三角为 1；
%   分数：|Theta_*| 的上三角绝对值（按 normalize_mode 归一化）。

if nargin < 1 || ~isstruct(args)
    error('compute_spd_and_classification_metrics:missing_args','args struct is required');
end
if nargin < 2 || ~isstruct(opts), opts = struct(); end

% ---------- 解包与默认 ----------
A_true = get_matrix(args, 'Theta_true');
A_est  = get_matrix(args, 'Theta_est');
A_unb  = get_matrix(args, 'Theta_unb');
A_ray  = get_matrix(args, 'Theta_ray');

p = size(A_true,1);

if ~isfield(args,'mask') || isempty(args.mask)
    mask = triu(true(p),1);
else
    mask = logical(args.mask);
    if ~isequal(size(mask),[p p])
        error('mask size mismatch: expected %dx%d', p, p);
    end
end

is_complex = true; if isfield(args,'is_complex'), is_complex = logical(args.is_complex); end

% opts defaults
normalize_mode = 'maxabs'; if isfield(opts,'normalize_mode'), normalize_mode = opts.normalize_mode; end

Sdefaults = struct('symmetrize',true,'project',true,'eps',1e-10);
if ~isfield(opts,'spd'), opts.spd = struct(); end
opts.spd = set_defaults(opts.spd, Sdefaults);

Adefaults = struct('value',0,'use_bcv_bug',false);
if ~isfield(opts,'alpha'), opts.alpha = struct(); end
opts.alpha = set_defaults(opts.alpha, Adefaults);

Pdefaults = struct('radar',false,'title','Quality Metrics Radar');
if ~isfield(opts,'plot'), opts.plot = struct(); end
opts.plot = set_defaults(opts.plot, Pdefaults);

% ---------- 预处理：对称化 + SPD 保障（仅用于距离计算） ----------
P = preprocess_spd(A_true, opts.spd);
E = preprocess_spd(A_est,  opts.spd);
U = preprocess_spd(A_unb,  opts.spd);
R = preprocess_spd(A_ray,  opts.spd);

% ---------- 六种 SPD 距离 ----------
alpha = opts.alpha.value; use_bug = opts.alpha.use_bcv_bug;
D_est = compute_all_distances(P,E,alpha,use_bug);
D_unb = compute_all_distances(P,U,alpha,use_bug);
D_ray = compute_all_distances(P,R,alpha,use_bug);

% ---------- 分类指标（AUC/SENS/SPEC/PREC/F1） ----------
Y = abs(A_true); Y = Y - diag(diag(Y)); Y = Y > 0; Yv = Y(mask);

M_est = abs(A_est); M_est = M_est - diag(diag(M_est)); Sv_est = vectorize_by_mask(M_est, mask);
M_unb = abs(A_unb); M_unb = M_unb - diag(diag(M_unb)); Sv_unb = vectorize_by_mask(M_unb, mask);
M_ray = abs(A_ray); M_ray = M_ray - diag(diag(M_ray)); Sv_ray = vectorize_by_mask(M_ray, mask);

Sv_est = normalize_scores(Sv_est, normalize_mode);
Sv_unb = normalize_scores(Sv_unb, normalize_mode);
Sv_ray = normalize_scores(Sv_ray, normalize_mode);

M_est_metrics = compute_classification_metrics(Yv, Sv_est);
M_unb_metrics = compute_classification_metrics(Yv, Sv_unb);
M_ray_metrics = compute_classification_metrics(Yv, Sv_ray);

% ---------- 汇总输出 ----------
out = struct();
out.distances = struct('est',D_est,'unb',D_unb,'ray',D_ray);
out.metrics   = struct('est',M_est_metrics,'unb',M_unb_metrics,'ray',M_ray_metrics);

% ---------- 可选画雷达图 ----------
if opts.plot.radar
    out.fig = struct();
    out.fig.radar = plot_radar( ...
        [M_est_metrics.auc, M_est_metrics.sens, M_est_metrics.spec, M_est_metrics.prec, M_est_metrics.f1;...
         M_unb_metrics.auc, M_unb_metrics.sens, M_unb_metrics.spec, M_unb_metrics.prec, M_unb_metrics.f1; ...
         M_ray_metrics.auc, M_ray_metrics.sens, M_ray_metrics.spec, M_ray_metrics.prec, M_ray_metrics.f1], ...
        {'est','unb','ray'}, {'AUC','SENS','SPEC','PREC','F1'}, opts.plot.title);
end

end  % main

% ===================== helpers =====================
function A = get_matrix(args, field)
    if ~isfield(args, field) || isempty(args.(field))
        A = [];
        return;
    end
    X = args.(field);
    if isstruct(X) && isfield(X,'X'), A = X.X; else, A = X; end
end

function A = preprocess_spd(A, spdopt)
    if isempty(A), A = []; return; end
    if spdopt.symmetrize
        A = 0.5*(A + A');
    end
    if spdopt.project
        % 简单 SPD 投影：特征值下截断
        if issparse(A)
            A = full(A);
        end
        [U,D] = eig(A);
        d = real(diag(D));
        dmin = max(spdopt.eps, 0);
        d = max(d, dmin);
        A = U*diag(d)*U';
        A = 0.5*(A + A');
    end
end

function s = normalize_scores(s, mode)
    if isempty(s), s = []; return; end
    switch lower(mode)
        case 'maxabs'
            m = max(abs(s)); if m>0, s = s./m; end
        case 'none'
            % no-op
        otherwise
            error('unknown normalize_mode: %s', mode);
    end
end

function v = vectorize_by_mask(M, mask)
    if isempty(M), v = []; return; end
    v = M(mask);
end

function D = compute_all_distances(P,A,alpha,use_bug)
    f = @(x) isempty(x);
    if f(A) || f(P)
        D = struct('kullback',NaN,'logeuclid',NaN,'riemann',NaN,'ld',NaN,'opttransp',NaN,'alphadiv',NaN);
        return;
    end
    D = struct();
    D.kullback  = distance_kullback(P,A);
    D.logeuclid = distance_logeuclid(P,A);
    D.riemann   = distance_riemann(P,A);
    D.ld        = distance_ld(P,A);
    D.opttransp = distance_opttransp(P,A);
    D.alphadiv  = distance_alphadiv(P,A,alpha,use_bug);
end

function d = distance_kullback(P,Q)
    % 对称 KL（Stein）：0.5 * tr(Q^{-1}P + P^{-1}Q - 2I)
    p = size(P,1);
    d = 0.5*( trace(Q\P) + trace(P\Q) - 2*p );
    d = real(d);
end

function d = distance_logeuclid(A,B)
    % Log-Euclid: || log(A) - log(B) ||_F
    LA = logm(A); LB = logm(B);
    d = norm(LB - LA, 'fro');
    d = real(d);
end

function d = distance_riemann(A,B)
    % AIRM: || log(A^{-1/2} B A^{-1/2}) ||_F
    As = sqrtm(A);
    As_inv = As\eye(size(A)); % 更稳的 A^{-1/2}
    C = As_inv * B * As_inv';
    d = norm(logm(C), 'fro');
    d = real(d);
end

function d = distance_ld(A,B)
    % JBLD 的平方根：sqrt(logdet((A+B)/2) - 0.5*logdet(A*B))
    C = 0.5*(A+B);
    ldC = logdet_spd(C);
    ldA = logdet_spd(A);
    ldB = logdet_spd(B);
    d = sqrt( max( ldC - 0.5*(ldA+ldB), 0) );
    d = real(d);
end

function d = distance_opttransp(A,B)
    % Bures–Wasserstein
    As = sqrtm(A);
    M  = As * B * As;
    d  = sqrt( real(trace(A) + trace(B) - 2*trace(sqrtm(M))) );
end

function d = distance_alphadiv(A,B,alpha,use_bug)
    % Alpha-divergence（对 SPD）：
    % 标准： (4/(1-α^2)) * log( det( ((1-α)/2)A + ((1+α)/2)B ) / (det(A)^((1-α)/2) det(B)^((1+α)/2)) )
    % 若 use_bug=true，分母按 BC-V 源码笔误：det(A)^((1-α)/2) * det(A)^((1+α)/2)
    if nargin<4, use_bug=false; end
    t1 = (1-alpha)/2; t2 = (1+alpha)/2;
    C = t1*A + t2*B;
    ldC = logdet_spd(C); ldA = logdet_spd(A); ldB = logdet_spd(B);
    if use_bug
        denom = ((1-alpha)/2)*ldA + ((1+alpha)/2)*ldA; % = ldA
    else
        denom = ((1-alpha)/2)*ldA + ((1+alpha)/2)*ldB;
    end
    d = (4/(1-alpha^2)) * ( ldC - denom );
    d = real(d);
end

function v = logdet_spd(A)
    % 通过 chol 更稳（支持复 Hermitian SPD）
    R = chol(A);
    v = 2*sum(log(abs(diag(R))));
    v = real(v);
end

function M = compute_classification_metrics(y, s)
    % y: logical (0/1), s: scores in [0,1]
    if isempty(s)
        M = struct('auc',NaN,'sens',NaN,'spec',NaN,'prec',NaN,'f1',NaN,'thr_opt',NaN,'fpr',[],'tpr',[]);
        return;
    end
    y = logical(y(:)); s = double(s(:));
    % ROC & AUC & optimal point (distance to (0,1))
    [fpr,tpr,thr,auc] = roc_curve_auc(y, s);
    [~,ix] = min( (fpr - 0).^2 + (tpr - 1).^2 );
    thr_opt = thr(ix);
    % confusion at optimal threshold
    yhat = s >= thr_opt;
    TP = sum( yhat & y );
    FP = sum( yhat & ~y );
    TN = sum( ~yhat & ~y );
    FN = sum( ~yhat & y );
    sens = safe_div(TP, TP+FN);
    spec = safe_div(TN, TN+FP);
    prec = safe_div(TP, TP+FP);
    f1   = safe_div(2*prec*sens, (prec+sens));
    M = struct('auc',auc,'sens',sens,'spec',spec,'prec',prec,'f1',f1,'thr_opt',thr_opt,'fpr',fpr,'tpr',tpr);
end

function a = safe_div(x,y)
    if y>0, a = x/y; else, a = NaN; end
end

function [fpr,tpr,thr,auc] = roc_curve_auc(y, s)
    % 无 toolbox 的稳健 ROC/AUC 实现
    % 按分数降序扫阈值
    [s_sorted, idx] = sort(s,'descend');
    y_sorted = y(idx);
    P = sum(y_sorted); N = sum(~y_sorted);
    if P==0 || N==0
        fpr = [0;1]; tpr = [0;1]; thr=[Inf;-Inf]; auc = NaN; return; end
    % 在所有 unique 分数及其下一个值之间取阈值
    [uniq_s, ia] = unique(s_sorted,'stable');
    thr = [Inf; (uniq_s(1:end-1)+uniq_s(2:end))/2; -Inf];
    % 预计算累计 TP/FP
    tp_cum = cumsum(y_sorted);
    fp_cum = cumsum(~y_sorted);
    % 在每个 ia 位置（分数跃迁末尾）记录点
    tp = [0; tp_cum(ia); P];
    fp = [0; fp_cum(ia); N];
    tpr = tp / P; fpr = fp / N;
    % AUC by trapezoid over FPR
    auc = trapz(fpr, tpr);
end

function h = plot_radar(M, legends, axes_labels, ttl)
    % M: K×5 (AUC,SENS,SPEC,PREC,F1)
    K = size(M,1);
    L = numel(axes_labels);
    th = linspace(0, 2*pi, L+1)'; th(end) = th(1);
    M = max(min(M,1),0); % clamp to [0,1]
    h = figure('Name','Radar Metrics');
    polaraxes; hold on;
    for k=1:K
        r = [M(k,:), M(k,1)];
        polarplot(th, r, 'LineWidth', 2);
    end
    thetaticks(rad2deg(th(1:end-1)));
    thetaticklabels(axes_labels);
    rticks(0:0.2:1);
    rlim([0 1]);
    legend(legends,'Location','bestoutside');
    title(ttl);
    hold off;
end

function S = set_defaults(S, D)
    fn = fieldnames(D);
    for i=1:numel(fn)
        k = fn{i};
        if ~isfield(S,k) || isempty(S.(k)), S.(k) = D.(k); end
    end
end
