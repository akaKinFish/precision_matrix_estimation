function summary = run_grid_module7_fixedL()
% 10×10×10×10 网格；首次调用 module7 产出 L 并固定；记录可复现材料
% 产出：
%   results/<timestamp>/summary.csv
%   results/<timestamp>/records.mat（逐 seed 明细，含随机种子与支持）
%   results/<timestamp>/summary.mat, cfg.json

%% ========= 实验/仿真基础配置 =========
base_sim_cfg = struct( ...
    'p', 3, 'n', 10, 'F', 3, 'T', 4096, ...     % 按需改；F≥3
    'random_seed', 1, ...
    'edge_density', 0.15, ...
    'lambda1_star', 0.35, ...
    'lambda2_star', 1.0, ...
    'laplacian_type','chain', ...
    'diag_base', 0.2, 'diag_spd', 1e-2, ...
    'complex_samples', true, ...
    'noise_type','scalar','sigma_xi2',0.2 ...
    );

Nseeds = 2;                     % 跨种子个数（可加大）
seed_list = 1001:(1000+Nseeds);  % 固定记录随机种子

% % 10×10×10×10 超参网格（可按你习惯微调范围）
% grid.lambda1 = [0.02 0.05 0.12];
% grid.lambda2_fac = [0.30 1.10 2.10];
% grid.lambda3_ratio = [0.00 0.25 0.50];
% grid.sigma_f = [1.5 4.0 6.0];
% ===== 精简后的超参网格 =====
% 固定 σ_f，省 3× 计算量
grid.sigma_f = 4.0;   % 标量即可，ndgrid 会处理成 1 维

% λ1：当前最优顶在 0.18 上界 → 向上扩到 0.24–0.26，并在 0.16–0.22 细化
grid.lambda1 = [0.90 1.00 1.05 1.10 1.15 1.20 1.25 1.30 1.40 1.50];

% λ3_ratio：当前最优顶在 0.60 上界 → 向上扩到 0.9，并在 0.5–0.8 细化
grid.lambda3_ratio = [1.10 1.20 1.25 1.30 1.35 1.40 1.50 1.60 1.80 2.00];

% λ2 factor：整体相关性不强，但最优在高端 → 上扩到 3.5–4.0，同时保留 1.x～3.0 的密度
grid.lambda2_fac = [14 16 18 20 22 24 26 30];

% （可选）打印一下当前组合规模，心里有数
fprintf('[GRID cfg] |λ1|=%d, |λ2fac|=%d, |λ3|=%d, |σ_f|=%d\n', ...
    numel(grid.lambda1), numel(grid.lambda2_fac), numel(grid.lambda3_ratio), numel(grid.sigma_f));


q_act = 0.10;                    % active set 分位数
refit_on = true;
warm_method = 'ssblpp';
inner_max_iter = 30;

use_parallel = true;             % 有并行工具箱就开
save_every = 25;                 % 每多少组合落一次盘（断点续跑友好）

%% ========= 输出目录 & 记录 =========
ts = datestr(now,'yyyymmdd_HHMMSS');
outdir = fullfile('results', ['grid_fixedL_',ts]);
if ~exist(outdir,'dir'), mkdir(outdir); end

cfg_dump = struct('base_sim_cfg',base_sim_cfg,'grid',grid, ...
    'Nseeds',Nseeds,'seed_list',{seed_list}, ...
    'q_act',q_act,'refit',refit_on,'warm',warm_method,'inner_max_iter',inner_max_iter);
savejson(fullfile(outdir,'cfg.json'), cfg_dump);

%% ========= 第一次调用：生成并固定 L =========
sim_cfg = base_sim_cfg; sim_cfg.random_seed = seed_list(1);
[Omega_true1, Sigma_true1, emp_cov_cell1, L_fixed, T, ~, simmeta1] = module7_simulation_matched(sim_cfg); %#ok<ASGLU>
% 以后都复用 L_fixed
fprintf('[FIXED L] size(L) = %dx%d 已固定。\n', size(L_fixed,1), size(L_fixed,2));

%% ========= 列举全部组合 =========
[L1, L2F, L3R, SF] = ndgrid(grid.lambda1, grid.lambda2_fac, grid.lambda3_ratio, grid.sigma_f);
combos = [L1(:), L2F(:), L3R(:), SF(:)];
Ncomb = size(combos,1);
fprintf('[GRID] 组合数 = %d\n', Ncomb);

% 预分配记录
Records = struct('lambda1',[],'lambda2_fac',[],'lambda3_ratio',[],'sigma_f',[], ...
    'seed',[],'metrics',[],'support',[]); %#ok<NASGU>
rec_idx = 0;

% 汇总容器
agg = table('Size',[Ncomb 1], ...
    'VariableTypes',"string",'VariableNames',"key"); %#ok<NASGU>
S = struct();  % 聚合器：S.(key).(metric).vals = [Nseeds×1]

%% ========= 并行池（可选）=========
%开并行池
if use_parallel
    p = gcp('nocreate');
    if isempty(p)
        try, parpool; catch, warning('并行池开启失败，改为串行'); use_parallel=false; end
    end
end

%% ========= 进度监控（并行可用）=========
total_tasks = Ncomb * Nseeds;   % 每个组合的每个 seed 都是一条任务
n_done = 0;
t0 = tic;
last_print = 0;
print_interval_sec = 1.0;       % 打印节流：至少隔这么久才打印一次

dq = parallel.pool.DataQueue;          % 工作者 -> 客户端
afterEach(dq, @onTick);                % 每收到一次就更新进度

fprintf('[PROG] 总任务 = %d（组合 %d × seeds %d）。\n', total_tasks, Ncomb, Nseeds);

%% ========= 主循环 =========
for c = 1:Ncomb
    lam1 = combos(c,1); fac2 = combos(c,2); lam3r = combos(c,3); sigf = combos(c,4);
    key = make_key(lam1, fac2, lam3r, sigf);      % -> 合法字段名
    key_raw = sprintf('l1=%.3g|l2f=%.3g|l3r=%.3g|sf=%.3g', lam1, fac2, lam3r, sigf); %#ok<NASGU>

    % 为该组合收集跨 seed 的度量
    dist_LE = zeros(Nseeds,1); dist_AIRM=zeros(Nseeds,1); dist_JBLD=zeros(Nseeds,1);
    dist_BW = zeros(Nseeds,1); dist_KLS = zeros(Nseeds,1);
    AUROC   = zeros(Nseeds,1); AUPR = zeros(Nseeds,1);
    F1K     = zeros(Nseeds,1); Prec = zeros(Nseeds,1); Rec = zeros(Nseeds,1);
    % 用于稳定性：存每个 seed 的支持（对 F 做 union）
    supp_bin = [];

    if use_parallel
        supp_loc = cell(Nseeds,1);
        parfor s = 1:Nseeds
            [dist_metrics, id_metrics, supp_vec] = ...
                run_one_seed(seed_list(s), L_fixed, lam1, fac2, lam3r, sigf, ...
                base_sim_cfg, q_act, refit_on, warm_method, inner_max_iter);
            LE = dist_metrics(1); AIRM=dist_metrics(2); JBLD=dist_metrics(3);
            BW = dist_metrics(4); KLS = dist_metrics(5);
            AU = id_metrics(1); AP  = id_metrics(2); F1 = id_metrics(3);
            PR = id_metrics(4); RC  = id_metrics(5);

            dist_LE(s)=LE; dist_AIRM(s)=AIRM; dist_JBLD(s)=JBLD;
            dist_BW(s)=BW; dist_KLS(s)=KLS;
            AUROC(s)=AU; AUPR(s)=AP; F1K(s)=F1; Prec(s)=PR; Rec(s)=RC;

            supp_loc{s} = supp_vec;

            % === 进度：该 seed 完成一次就 tick 一下（并行安全）===
            send(dq, 1);
        end
        supp_bin = logical(cell2mat(supp_loc(:)));

    else
        supp_list = cell(Nseeds,1);
        for s = 1:Nseeds
            [dist_metrics, id_metrics, supp_vec] = ...
                run_one_seed(seed_list(s), L_fixed, lam1, fac2, lam3r, sigf, ...
                base_sim_cfg, q_act, refit_on, warm_method, inner_max_iter);
            LE = dist_metrics(1); AIRM=dist_metrics(2); JBLD=dist_metrics(3);
            BW = dist_metrics(4); KLS = dist_metrics(5);
            AU = id_metrics(1); AP  = id_metrics(2); F1 = id_metrics(3);
            PR = id_metrics(4); RC  = id_metrics(5);

            dist_LE(s)=LE; dist_AIRM(s)=AIRM; dist_JBLD(s)=JBLD;
            dist_BW(s)=BW; dist_KLS(s)=KLS;
            AUROC(s)=AU; AUPR(s)=AP; F1K(s)=F1; Prec(s)=PR; Rec(s)=RC;

            supp_list{s} = supp_vec;

            % 串行也复用同一计数逻辑
            onTick(1);
        end
        supp_bin = logical(cell2mat(supp_list(:)));
    end

    % 稳定性：跨 seed 两两 Jaccard
    J = jaccard_pairs(supp_bin);
    J = J(:);
    J = J(isfinite(J));
    if isempty(J)
        Jmin = NaN; Jmed = NaN; Jmax = NaN; %#ok<NASGU>
    else
        Jmin = min(J); Jmed = median(J); Jmax = max(J); %#ok<NASGU>
    end

    % 聚合保存
    S.(key).dist_LE = dist_LE;    S.(key).dist_AIRM = dist_AIRM;
    S.(key).dist_JBLD = dist_JBLD;S.(key).dist_BW   = dist_BW;
    S.(key).dist_KLS = dist_KLS;
    S.(key).AUROC = AUROC; S.(key).AUPR=AUPR; S.(key).F1K=F1K; S.(key).Prec=Prec; S.(key).Rec=Rec;
    S.(key).Jaccard = J;

    % 追加到 records（轻量：仅存支持和指标）
    rec_idx = rec_idx + 1;
    Records(rec_idx,1).lambda1 = lam1; %#ok<AGROW>
    Records(rec_idx,1).lambda2_fac = fac2;
    Records(rec_idx,1).lambda3_ratio = lam3r;
    Records(rec_idx,1).sigma_f = sigf;
    Records(rec_idx,1).seed = seed_list;
    Records(rec_idx,1).metrics = struct('LE',{dist_LE},'AIRM',{dist_AIRM},'JBLD',{dist_JBLD}, ...
        'BW',{dist_BW},'KLS',{dist_KLS},'AUROC',{AUROC},'AUPR',{AUPR},'F1K',{F1K},'Prec',{Prec},'Rec',{Rec});
    Records(rec_idx,1).support = supp_bin;

    % 定期落盘
    if mod(c, save_every)==0 || c==Ncomb
        save(fullfile(outdir,'records.mat'),'Records','-v7.3');
        save(fullfile(outdir,'summary.mat'),'S','-v7.3');
        writetable(make_summary_table(S), fullfile(outdir,'summary.csv'));
        fprintf('[SAVE] 已保存到 %s （进度 %d/%d 组合）\n', outdir, c, Ncomb);
    end
end

% 最终 summary
summary = make_summary_table(S);
writetable(summary, fullfile(outdir,'summary.csv'));
save(fullfile(outdir,'summary.mat'),'S','-v7.3');
fprintf('[DONE] 总表写入 %s\n', fullfile(outdir,'summary.csv'));

%% ====== 嵌套：进度更新回调 ======
    function onTick(~)
        % 每次收到一个“seed 完成”的 tick
        n_done = n_done + 1;

        % 节流打印（每秒最多一次；或全部完成时必打）
        t_now = toc(t0);
        if (t_now - last_print) >= print_interval_sec || n_done == total_tasks
            pct = 100 * n_done / max(total_tasks,1);
            rate = n_done / max(t_now, eps);              % 任务/秒
            remain = (total_tasks - n_done) / max(rate,1e-9);
            fprintf('[PROG] %6d/%6d  (%.1f%%)  |  elapsed %s  |  ETA %s\n', ...
                n_done, total_tasks, pct, fmt_time(t_now), fmt_time(remain));
            last_print = t_now;
        end
    end

    function s = fmt_time(t)
        if ~isfinite(t) || t<0, t = 0; end
        hh = floor(t/3600); t = t - 3600*hh;
        mm = floor(t/60);   ss = t - 60*mm;
        s = sprintf('%02d:%02d:%05.2f', hh, mm, ss);
    end

end % ===== 主函数结束 =====



%% ====== 单 seed/单组合的实验 ======
function [dist_metrics, id_metrics, supp_vec] = run_one_seed(seed, L_fixed, lam1, fac2, lam3r, sigf, ...
    base_sim_cfg, q_act, refit_on, warm_method, inner_max_iter)

% 1) 仿真（固定 L）
sim_cfg = base_sim_cfg;
sim_cfg.random_seed = seed;
sim_cfg.L_fixed = L_fixed;   % 关键：固定 lead field
[Omega_true, Sigma_true, emp_cov_cell, ~, T, ~, simmeta] = module7_simulation_matched(sim_cfg); %#ok<ASGLU>

% 2) 估计
cfg = struct();
cfg.kernel_sigma   = sigf;
cfg.lambda3_ratio  = lam3r;
cfg.do_support_refit = refit_on;
cfg.warm_method    = warm_method;
cfg.q_act          = q_act;           % 需要 adapter 补丁支持
cfg.lambda1        = lam1;            % 覆盖 λ1
cfg.lambda2_factor = fac2;            % 相对建议值放缩
cfg.em = struct(); cfg.em.inner = struct();
cfg.em.inner.max_iter = inner_max_iter;

[Omega_est, Dsrc_est, ~, ~] = my_test2_adapter(emp_cov_cell, L_fixed, T, cfg, Omega_true);

% 3) 距离（SPD 距离对协方差，逐频取均值）
F = numel(emp_cov_cell);
dLE=0; dAIRM=0; dJBLD=0; dBW=0; dKLS=0;
for f=1:F
    A = spd_project(Dsrc_est{f});
    B = spd_project(Sigma_true(:,:,f));
    dLE   = dLE   + dist_LE(A,B);
    dAIRM = dAIRM + dist_AIRM(A,B);
    dJBLD = dJBLD + dist_JBLD(A,B);
    dBW   = dBW   + dist_BW(A,B);
    dKLS  = dKLS  + dist_KLS(A,B);
end
dist_metrics = [dLE,dAIRM,dJBLD,dBW,dKLS]/F;

% 4) 识别（基于 Ω 的支持）
[score, label, supp_vec] = edge_scores_labels(Omega_est, Omega_true);
score = real(score(:));
label = logical(label(:));
keep = isfinite(score) & (label==0 | label==1);
score = score(keep);
label = label(keep);
[auroc, aupr] = roc_pr(score, label);
[Kf1, P, R]   = f1_at_K(score, label);

id_metrics = [auroc, aupr, Kf1, P, R];
end


%% ====== —— 指标与工具函数 —— ======

function A = spd_project(A)
A = (A+A')/2;
[V,D] = eig(A,'vector'); D = real(D); D(D<1e-12)=1e-12;
A = V*diag(D)*V';
A = (A+A')/2;
end

% Log-Euclidean
function d = dist_LE(A,B)
LA = logm(A); LB = logm(B);
d = norm(LA-LB,'fro');
end

% AIRM geodesic distance
function d = dist_AIRM(A,B)
X = sqrtm(A);
Y = X \ (B / X);   % A^{-1/2} * B * A^{-1/2}
eigv = eig((Y+Y')/2); eigv = max(real(eigv), 1e-18);
d = norm(log(eigv));
end

% JBLD
function d = dist_JBLD(A,B)
C = 0.5*(A+B);
d = logdet(C) - 0.5*(logdet(A)+logdet(B));
end

% Bures–Wasserstein
function d = dist_BW(A,B)
A12 = sqrtm(A);
X = A12 * sqrtm(spd_project(A12\B/A12)) * A12;
d = sqrt( max(trace(A + B - 2*X), 0) );
end

% Symmetric KL (covariance form)
function d = dist_KLS(A,B)
n = size(A,1);
AiB = trace(A\B);
BiA = trace(B\A);
d = 0.5*(AiB + BiA - 2*n);
end

function v = logdet(M)
% 安全 logdet
[L,p] = chol((M+M')/2,'lower');
if p>0, M = spd_project(M); L = chol(M,'lower'); end
v = 2*sum(log(diag(L)));
end

function [score, label, supp_vec] = edge_scores_labels(Omega_est, Omega_true)
% 兼容 cell 或 n×n×F 形式；返回：
%   score/label：用于 ROC/PR 的拼接分数和标签
%   supp_vec   ：跨频 union 的支持向量（用于稳定性）

Oe = ensure_cell_(Omega_est);
Ot = ensure_cell_(Omega_true);

F = numel(Oe);
n = size(Oe{1},1);
mask = triu(true(n),1);

scores_all = [];
labels_all = [];
supp_union = false(nnz(mask),1);

for f = 1:F
    Eh = abs(real((Oe{f}+Oe{f}')/2));
    Eh(1:n+1:end) = 0;
    scores_all = [scores_all; Eh(mask)]; %#ok<AGROW>

    Gt = abs(real((Ot{min(f,numel(Ot))}+Ot{min(f,numel(Ot))}')/2));
    Gt(1:n+1:end) = 0;
    labels_all = [labels_all; (Gt(mask) > 1e-12)]; %#ok<AGROW>

    % 稳定性：跨频 union
    supp_union = supp_union | (Eh(mask) > 1e-8);
end

score = scores_all;
label = labels_all;
supp_vec = supp_union;
end

function C = ensure_cell_(X)
% 数组 -> cell；cell 直接返回
if iscell(X)
    C = X(:);
elseif isnumeric(X) && ndims(X)==3 && size(X,1)==size(X,2)
    F = size(X,3);
    C = cell(F,1);
    for f = 1:F
        C{f} = X(:,:,f);
    end
else
    error('edge_scores_labels:invalidInput','Expect cell{F} or n×n×F numeric.');
end
end


function [auroc, aupr, curves] = roc_pr(score, label)
% ROC/PR 评估，兼容：
%  - score/label 任意行列形状（内部会矫正为列向量）
%  - label logical/double 均可
%  - 退化情况（全正/全负）直接返回 NaN
%  - 无 perfcurve 时，使用手写 fallback

% 形状与类型
score = real(score(:));
label = logical(label(:));
assert(numel(score)==numel(label), 'roc_pr: length mismatch');

% 去掉 NaN/Inf
keep = isfinite(score);
score = score(keep); label = label(keep);

% 退化情况
hasPos = any(label);
hasNeg = any(~label);
if ~(hasPos && hasNeg)
    auroc = NaN; aupr = NaN;
    if nargout>2, curves = struct('fpr',[],'tpr',[],'rec',[],'prec',[]); end
    warning('roc_pr:degenerate','Labels are all one class; AUROC/AUPR undefined.');
    return;
end

% 有 perfcurve 先用它
if exist('perfcurve','file')==2
    [fpr, tpr,~,auroc] = perfcurve(label, score, true);
    [rec, prec,~,aupr] = perfcurve(label, score, true, 'xCrit','reca','yCrit','prec');
else
    % --- 手写 fallback ---
    [~, idx] = sort(score,'descend');
    y = label(idx);
    P = sum(y); N = numel(y)-P;
    tp = cumsum(y);
    fp = cumsum(~y);
    tpr = tp / P;
    fpr = fp / N;
    % ROC 面积（带端点）
    auroc = trapz([0; fpr; 1], [0; tpr; 1]);
    % PR 曲线 & 面积（常用作近似）
    prec = tp ./ max(tp+fp, 1);
    rec  = tpr;
    aupr = trapz([0; rec; 1], [prec(1); prec; P/(P+N)]);
end

if nargout>2
    curves = struct('fpr',fpr,'tpr',tpr,'rec',rec,'prec',prec);
end
end


function [F1, P, R] = f1_at_K(score, label)
K = sum(label); K = max(K,1);
[ss, idx] = sort(double(score),'descend'); %#ok<ASGLU>
pick = idx(1:min(K,numel(idx)));
tp = sum(label(pick));
P = tp / numel(pick);
R = tp / sum(label);
if P+R==0, F1=0; else, F1 = 2*P*R/(P+R); end
end

function J = jaccard_pairs(S)
%JACCARD_PAIRS  计算列之间（K 个切片）的两两 Jaccard 指数
% 输入 S 可为：
%   - p×p×K 的逻辑/数值（每个切片是一张对称邻接或支持矩阵）
%   - E×K   的逻辑/数值（每列是展平后的上三角边向量）
% 输出：
%   - J : (K*(K-1)/2)×1 的向量，每一项是两列（两个切片）的 Jaccard

if isempty(S)
    J = NaN; return;
end

if ndims(S) == 3
    % S 是 p×p×K：把上三角展开为 E×K
    [p,~,K] = size(S);
    mask = triu(true(p),1);
    E = nnz(mask);
    B = false(E, K);
    for k = 1:K
        A = S(:,:,k) ~= 0;
        A = A | A.';            % 保险起见，强制对称
        B(:,k) = A(mask);
    end
else
    % S 是 E×K：直接使用
    [~,K] = size(S);
    B = S ~= 0;
end

% K<2：没法两两比较，定义为 NaN 返回
if K < 2
    J = NaN;  % 或者 J = []; 看你后面怎么用
    return;
end

% 逐对计算 Jaccard（不用 nchoosek 防止小 K 报错）
npairs = K*(K-1)/2;
J = zeros(npairs,1);
t = 1;
for i = 1:K-1
    ai = B(:,i);
    for j = i+1:K
        aj = B(:,j);
        u = nnz(ai | aj);
        if u == 0
            J(t) = NaN;         % 两个切片都全零：J 未定义
        else
            J(t) = nnz(ai & aj) / u;
        end
        t = t + 1;
    end
end
end


function T = make_summary_table(S)
% S.(key).<metric> = [Nseeds×1]
keys = fieldnames(S);
n = numel(keys);
% 定义输出列
colnames = ["key", ...
    "LE_mean","LE_std","AIRM_mean","AIRM_std","JBLD_mean","JBLD_std", ...
    "BW_mean","BW_std","KLS_mean","KLS_std", ...
    "AUROC_mean","AUROC_std","AUPR_mean","AUPR_std", ...
    "F1_mean","F1_std","Prec_mean","Prec_std","Rec_mean","Rec_std", ...
    "Jaccard_mean","Jaccard_std"];
T = cell(n, numel(colnames));

for i=1:n
    k = keys{i};
    get = @(f) S.(k).(f);
    row = { ...
        string(k), ...
        mean(get('dist_LE')), std(get('dist_LE')), ...
        mean(get('dist_AIRM')), std(get('dist_AIRM')), ...
        mean(get('dist_JBLD')), std(get('dist_JBLD')), ...
        mean(get('dist_BW')),   std(get('dist_BW')), ...
        mean(get('dist_KLS')),  std(get('dist_KLS')), ...
        mean(get('AUROC')),     std(get('AUROC')), ...
        mean(get('AUPR')),      std(get('AUPR')), ...
        mean(get('F1K')),       std(get('F1K')), ...
        mean(get('Prec')),      std(get('Prec')), ...
        mean(get('Rec')),       std(get('Rec')), ...
        mean(get('Jaccard')),   std(get('Jaccard')) ...
        };
    T(i,:) = row;
end

T = cell2table(T, 'VariableNames', cellstr(colnames));
% 排序提示：优先 LE 低、AUPR/F1 高、Jaccard 高
T = sortrows(T, {'LE_mean','AUPR_mean','F1_mean','Jaccard_mean'}, {'ascend','descend','descend','descend'});
end

function savejson(fname, s)
txt = jsonencode(s, 'PrettyPrint',true);
fid=fopen(fname,'w'); fwrite(fid,txt); fclose(fid);
end

function key = make_key(l1, l2f, l3r, sf)
    % 先做一个可读的 key，再转为合法字段名，避免非法字符
    txt = sprintf('l1_%0.5g__l2f_%0.5g__l3r_%0.5g__sf_%0.5g', l1, l2f, l3r, sf);
    key = matlab.lang.makeValidName(txt, 'ReplacementStyle','underscore');
end

