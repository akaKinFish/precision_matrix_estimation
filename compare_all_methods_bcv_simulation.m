function results = compare_all_methods_bcv_simulation(cfg)
% COMPARE_ALL_METHODS  仿真 → 估计 →（BC-V风格）评测 → 控制台输出 → 保存
%
% 评测改为三块：
%   A) 偏相干并排面板图（per-frequency）
%   B) 似然驱动的阈值扫描（对软阈比例做网格，选最大对数似然）
%   C) SPD矩阵距离表（KLS对称 / Log-Euclidean / AIRM / JBLD / Bures-W）

%% ========= 0) 配置 =========
if nargin < 1, cfg = struct(); end

% === 默认基础配置（与你原来一致） ===
default_cfg = struct();
default_cfg.n_nodes = 10;
default_cfg.n_sensors = 3;
default_cfg.n_freq = 3;
default_cfg.n_samples = 4096;
default_cfg.random_seed = 2025;
default_cfg.leadfield_type = 'simple';
default_cfg.simulator = 'matched';
default_cfg.edge_density = 0.15;
default_cfg.support_strategy = 'threshold';
default_cfg.modeN = 'diff';
default_cfg.lambda1_star = 0.35;
default_cfg.lambda2_star = 1.0;
default_cfg.sigma_xi2 = 0.2;
default_cfg.delta_eps = [0.08, 1e-6];
default_cfg.noise_type = 'scalar';
default_cfg.complex_samples = true;

% sim_matched 子结构
sim_matched_default = struct();
sim_matched_default.laplacian_type = 'chain';
sim_matched_default.tau_gp = 1.0;
sim_matched_default.eps_gp = 1e-3;
sim_matched_default.diag_base = 0.2;
sim_matched_default.diag_spd = 1e-2;
default_cfg.sim_matched = sim_matched_default;

% my_test2 baseline
default_cfg.adapter_sigma_xi2 = [];
default_cfg.adapter_lsq_ridge = 5e-3;
default_cfg.adapter_cov_ridge = 1e-6;
default_cfg.adapter_inv_ridge = 5e-2;
default_cfg.adapter_soft_tau = 0.10;
default_cfg.adapter_spd_delta = 0.03;
default_cfg.adapter_spd_eps = 1e-6;
default_cfg.adapter_enable_gls = true;
default_cfg.adapter_enable_refit = true;
default_cfg.adapter_refit_ridge = 1e-3;

% freq_couple 子结构
freq_couple_default = struct();
freq_couple_default.enabled = true;
freq_couple_default.rho = 0.30;
freq_couple_default.n_iter = 3;
default_cfg.adapter_freq_couple = freq_couple_default;

% 横向扫参
default_cfg.grid_mytest2 = {};

% === 新增：BC-V 风格评测配置 ===
eval_default = struct();
eval_default.enable_prroc = false;
eval_default.tau_scan = linspace(0, 0.3, 31);
eval_default.dist_ref = 'omega_true';
eval_default.panel_caxis = [];
eval_default.panel_clip = 0.99;
default_cfg.eval_bcv = eval_default;

% 保留这些旧参数（兼容你的函数）
default_cfg.num_thresh = 128;
default_cfg.tol_true = 1e-10;

cfg = set_default(cfg, default_cfg);
cfg.eval_bcv = set_default(cfg.eval_bcv, eval_default);

rng(cfg.random_seed);

ts = datestr(now,'yyyymmdd_HHMMSS');
outdir = fullfile('results', ['exp_', ts]);
if ~exist(outdir,'dir'), mkdir(outdir); end

fprintf('================= compare_all_methods =================\n');

%% ========= 1) 仿真 =========
fprintf('==[1/5] Simulation (%s generator) ==\n', cfg.simulator);
process_waitbar = waitbar(0,'Please wait...');
%%
%% Create partial correlations (ThetaJJ)
%  Generate source (state) empirical covariance (SJJ) and source (state) activity (J)
m                    = 600;             % Sample number
q                    = 10;              % Number of generators
p                    = 3;              % Number of sensors
nblocks              = 2;               % Number of blocks in simulation
options.config       = 2;               % (2) overlapping blocks (1) nonoverlapping blocks
options.var          = 2;               % (2) complex variable (1) real variable
options.extensions   = [ceil(q/3); ceil(q/3); q - 2*ceil(q/3)]; % patches extensions
options.connections  = [1 2; 2 3];      % patches connections
% [Sjj_sim,j_sim,Thetajj_sim] = gen_hggm1(m,q,nblocks,options);
% j_sim = transpose(j_sim);

[Sjj_sim,j_sim,Thetajj_sim] = gen_hggm2(m,q,options);
%% Creating pseudoLead Field (L)
Lvj             = zeros(p,q);
radj            = 60;
radv            = 85;
angj            = 2*pi/q;
angv            = 2*pi/p;
for contv = 1:p
    for contj = 1:q
        waitbar((contv*contj)/(p*q),process_waitbar,strcat('Creating pseudoLead Field (L)'));
        vectv            = [radv*cos((contv-1)*angv); radv*sin((contv-1)*angv)];
        vectj            = [radj*cos((contj-1)*angj); radj*sin((contj-1)*angj)];
        r                = vectv - vectj;
        r_unit           = r/sqrt(sum(abs(r).^2));
        miu              = vectj/sqrt(sum(abs(vectj).^2));
        Lvj(contv,contj) = (1/(4*pi))*miu'*r_unit/sqrt(sum(abs(r).^2))^2;
    end
end
delete(process_waitbar);


% LeadFields      = {Lvj};
% save('LeadFields_pseudo','LeadFields')
%% pseudo-cortex

process_waitbar = waitbar(0,'Please wait...');
vertices        = zeros(q,2);
for contj = 1:q
    waitbar((contj)/(q),process_waitbar,strcat('pseudo-cortex'));
    vertices(contj,:)    = [radj*cos((contj-1)*angj) radj*sin((contj-1)*angj)];
end
delete(process_waitbar);

process_waitbar = waitbar(0,'Please wait...');
faces           = [[1:q]' [2:q 1]'];
cortex.vertices = vertices;
cortex.faces    = faces;
coor            = zeros(p,2);
for contv = 1:p
    waitbar((contv)/(p),process_waitbar,strcat('HeadModel-pseudo ',' cortex ',' coor'));
    coor(contv,:)    = [radv*cos((contv-1)*angv) radv*sin((contv-1)*angv)];
end
delete(process_waitbar);

%%
process_waitbar = waitbar(0,'Please wait...');

%% Generate data
v0              = Lvj*j_sim; % data (observation)
%% Biological noise
bionoise        = randn(q,m) + 1i*randn(q,m);
bionoise        = Lvj*bionoise;
bionoise        = sum(abs(v0(:)).^2)^(1/2)*bionoise/sum(abs(bionoise(:)).^2)^(1/2);
%% Sensor noise
sensnoise       = randn(p,m) + 1i*randn(p,m);
sensnoise       = sum(abs(v0(:)).^2)^(1/2)*sensnoise/sum(abs(sensnoise(:)).^2)^(1/2);
%% Corrupted data
v               = v0 + 0.1*bionoise + 0.1*sensnoise;
%% Data empirical covariance
Svv             = cov(v');
%% Data empirical covariance
Svv             = cov(v');

% ====== [新增]：把“真值”补齐，统一为 BC-V 语义 ======
T               = m;                 % 样本数 = 你用来造 Svv 的时段数
F               = 1;                 % 这一版仿真是单频
emp_cov_cell    = {Svv};
L               = Lvj;

% 【关键修正】Omega_true 应当是精度真值（Thetajj_sim），不是协方差
Omega_true      = { hermitize(Thetajj_sim) };

% 如果后面评测会用到 Sigma_true，就顺便给出真值协方差（Ω 的稳健逆）
Sigma_true      = { inv_psd( Omega_true{1} ) };

% 【新增】构造 bayes_truth（给 GLS 白化/似然打分用）
% 只把“传感器噪声”当作观测噪声真值；生理噪声属于源域扰动，不并入 Σ_xixi
Sigma_xixi_true = cov( (0.1*sensnoise)' );  % p×p
Sigma_xixi_true = (Sigma_xixi_true + Sigma_xixi_true')/2;    % Hermitian
bayes_truth = struct();
bayes_truth.noise_type        = 'matrix';    % 告诉后续：我提供了完整的 Σ_xixi 真值
bayes_truth.Sigma_xixi_true   = Sigma_xixi_true;

% （可选，仅用于日志展示，不参与计算）
bayes_truth.lambda1_true = getd(cfg,'lambda1_star',NaN);
bayes_truth.lambda2_true = getd(cfg,'lambda2_star',NaN);

% （可选）记录仿真元信息，方便保存/复现
sim = struct('p', q, 'm', m, 'F', F, 'T', T);


%% ========= 2) 你的算法（my_test2 扫参） =========
fprintf('==[2/5] Run YOUR algorithm (my_test2_adapter, hyperparam sweep) ==\n');

% baseline 适配器配置
cfg_est_base = struct();
bayes_truth.noise_type = 'matrix';
switch lower(bayes_truth.noise_type)
    case 'scalar'
        % 修正：使用 getd 替代 ifempty
        sigma_val = getd(cfg, 'adapter_sigma_xi2', []);
        if isempty(sigma_val)
            cfg_est_base.sigma_xi2 = bayes_truth.sigma_xi2_true;
        else
            cfg_est_base.sigma_xi2 = sigma_val;
        end
    case 'matrix'
        cfg_est_base.Sigma_xixi = bayes_truth.Sigma_xixi_true;
    case 'matrix_per_freq'
        cfg_est_base.Sigma_xixi_per_freq = bayes_truth.Sigma_xixi_per_freq_true;
end
cfg_est_base.lsq_ridge   = cfg.adapter_lsq_ridge;
cfg_est_base.cov_ridge   = cfg.adapter_cov_ridge;
cfg_est_base.inv_ridge   = cfg.adapter_inv_ridge;
cfg_est_base.soft_tau    = cfg.adapter_soft_tau;
cfg_est_base.spd_delta   = cfg.adapter_spd_delta;
cfg_est_base.spd_eps     = cfg.adapter_spd_eps;
cfg_est_base.enable_gls  = cfg.adapter_enable_gls;
cfg_est_base.enable_refit= cfg.adapter_enable_refit;
cfg_est_base.refit_ridge = cfg.adapter_refit_ridge;
cfg_est_base.freq_couple = cfg.adapter_freq_couple;

% 组装 runs
runs = struct('name',{},'cfg_est',{});
runs(end+1).name = 'my_test2_base';
runs(end).cfg_est = cfg_est_base;

grid_list = getd(cfg,'grid_mytest2',{});
for i = 1:numel(grid_list)
    delta = grid_list{i};
    nm = getd(delta,'name',sprintf('my2_%02d',i));
    if isfield(delta,'name'), delta = rmfield(delta,'name'); end
    cfg_i = override_struct(cfg_est_base, delta);
    runs(end+1).name = nm; %#ok<AGROW>
    runs(end).cfg_est = cfg_i;
end

% 逐个运行（合法化+去重）
Omega_all = struct();
outs_map  = struct();
for i = 1:numel(runs)
    nm    = runs(i).name;
    cfg_i = runs(i).cfg_est;
    try
        [Om_i, ~, ~, outs_i] = my_test2_adapter(emp_cov_cell, L, T, cfg_i, Omega_true);
    catch
        try
            [Om_i, ~, ~] = my_test2_adapter(emp_cov_cell, L, T, cfg_i);
            outs_i = struct();
        catch ME
            warning('Method %s failed: %s', nm, ME.message);
            continue;
        end
    end
    nm_field = make_valid_name(nm);
    nm_field = uniquify_field(nm_field, fieldnames(Omega_all));
    Omega_all.(nm_field) = Om_i;
    outs_map.(nm_field)  = outs_i;
end

%% ========= 3) （可选）BC-V 方法 =========
if exist('run_bcv_methods','file') == 2
    fprintf('==[3/5] Run BC-V methods ==\n');
    try
        bcv = run_bcv_methods(emp_cov_cell, L, T, struct());
        if exist('bcv','var')
            fbcv = fieldnames(bcv);
            for ii=1:numel(fbcv)
                if isfield(bcv,fbcv{ii})
                    nm_field = make_valid_name(['bcv_' fbcv{ii}]);
                    nm_field = uniquify_field(nm_field, fieldnames(Omega_all));
                    Omega_all.(nm_field) = bcv.(fbcv{ii});
                end
            end
        end
    catch ME
        warning('[BC-V] failed: %s', ME.message);
    end
else
    fprintf('==[3/5] Skip BC-V (run_bcv_methods not found) ==\n');
end

method_names = fieldnames(Omega_all);

%% ========= 4) 评测（BC-V风格） =========
fprintf('==[4/5] Evaluate (BC-V style) ==\n');

% 统一精度格式
OmTrue = coerce_precision_stack(Omega_true);
SiTrue = coerce_precision_stack(Sigma_true);

% ---- A) 偏相干并排面板图 ----
panel_files = plot_pcoh_panels(Omega_all, OmTrue, cfg.eval_bcv, outdir);

% ---- B) 似然驱动阈值扫描 ----
Sii_proxy = make_source_cov_proxy(emp_cov_cell, L, bayes_truth, cfg_est_base);
lik_scan = struct();
for k = 1:numel(method_names)
    name = method_names{k};
    OmEst = coerce_precision_stack(Omega_all.(name));
    [best_tau, tau_grid, ll_curve] = likelihood_scan_tau(OmEst, Sii_proxy, cfg.eval_bcv.tau_scan, cfg_est_base);
    lik_scan.(name) = struct('best_tau',best_tau, 'tau_grid',tau_grid, 'll_curve',ll_curve);
    fprintf('  [lik-scan] %s: best tau=%.4g (grid in [%.3g, %.3g], %d pts)\n', ...
        name, best_tau, min(tau_grid), max(tau_grid), numel(tau_grid));
end
save(fullfile(outdir, 'likelihood_scan.mat'), 'lik_scan');

% ---- C) SPD 距离表 ----
dist_summary = struct();
for k = 1:numel(method_names)
    name = method_names{k};
    OmEst = coerce_precision_stack(Omega_all.(name));
    % 用 Sigma_hat 与 Sigma_true 比更稳（SPD）
    SigHat = zeros(size(OmEst));
    for f=1:size(OmEst,3), SigHat(:,:,f) = inv_psd(OmEst(:,:,f)); end
    if strcmpi(cfg.eval_bcv.dist_ref,'omega_true')
        Astack = OmEst;  Bstack = OmTrue;
    else
        Astack = SigHat; Bstack = SiTrue;
    end
    D = spd_distance_pack(Astack, Bstack);
    dist_summary.(name) = summarize_distances(D);
end

% 导出 CSV
csv_path = fullfile(outdir,'spd_distance_summary.csv');
write_distance_csv(csv_path, dist_summary);
fprintf('  [dist] SPD distance summary saved: %s\n', csv_path);

% ---- （可选）保留 PR/ROC ----
metrics = struct();
method_summaries = struct();
if cfg.eval_bcv.enable_prroc
    fprintf('  [opt] Also computing PR/ROC ...\n');
    % 真值支持与其跨频Jaccard
    gt_stats = support_stats(OmTrue, cfg.tol_true);
    print_support_stats('GroundTruth', gt_stats);
    for k = 1:numel(method_names)
        name = method_names{k};
        OmEst = coerce_precision_stack(Omega_all.(name));
        eval_cfg = struct(); eval_cfg.num_thresh = cfg.num_thresh; eval_cfg.tol_true = cfg.tol_true;
        metrics.(name) = eval_support_curves(OmEst, OmTrue, eval_cfg);
        [thr_match, F1_match, P_match, R_match] = threshold_matching_sparsity_per_freq(OmEst, OmTrue, cfg.tol_true);
        relF_off = rel_frob_offdiag(OmEst, OmTrue);
        Mask_est  = threshold_support_per_freq(OmEst, OmTrue, cfg.tol_true);
        est_stats = support_stats_mask(Mask_est);
        summary_struct = struct();
        summary_struct.relF_off  = relF_off;
        summary_struct.thr_match = thr_match;
        summary_struct.F1_match  = F1_match;
        summary_struct.P_match   = P_match;
        summary_struct.R_match   = R_match;
        summary_struct.est_stats = est_stats;
        method_summaries.(name)  = summary_struct;
        print_method_summary(name, metrics.(name), summary_struct);
    end
    plot_eval_curves_multi(metrics, outdir);
end

%% ========= 5) 保存 =========
results = struct();
results.cfg = cfg;
results.outdir = outdir;
results.Omega_true = OmTrue;
results.Omega_all = Omega_all;
results.L = L; results.T = T;
results.bayes_truth = bayes_truth;
results.sim = sim;
results.adapter_outs = outs_map;

% 新增 BC-V 风格评测产物
results.bcv_eval = struct();
results.bcv_eval.panel_files = panel_files;
results.bcv_eval.lik_scan = lik_scan;
results.bcv_eval.dist_tables = dist_summary;

if cfg.eval_bcv.enable_prroc
    results.metrics = metrics;
    results.method_summaries = method_summaries;
end

save(fullfile(outdir, 'results_all.mat'), '-struct', 'results');
fprintf('\n[OK] Done. Artifacts saved under: %s\n', outdir);
end

%% ==================== Local helpers ====================

function Sii_proxy = make_source_cov_proxy(emp_cov_cell, L, bayes_truth, cfg_est_base)
% 用 GLS 白化 + 最小范数反演，构造"源域样本协方差"的 proxy
F = numel(emp_cov_cell);
[n, p] = size(L);
% 白化
Wcell = cell(F,1);
switch lower(bayes_truth.noise_type)
    case 'matrix_per_freq'
        for f=1:F, Wcell{f} = inv_chol(bayes_truth.Sigma_xixi_per_freq_true{f}); end
    case 'matrix'
        W = inv_chol(bayes_truth.Sigma_xixi_true); for f=1:F, Wcell{f} = W; end
    otherwise
        W = eye(n)/sqrt(bayes_truth.sigma_xi2_true); for f=1:F, Wcell{f} = W; end
end
Sii_proxy = zeros(p,p,F);
for f=1:F
    Svv = hermitize(emp_cov_cell{f});
    W   = Wcell{f};
    Lw  = W * L;
    Svvw= hermitize(W * Svv * W');
    G = (Lw'*Lw + cfg_est_base.lsq_ridge*eye(p)) \ (Lw');
    Sii_proxy(:,:,f) = hermitize(G * Svvw * G') + cfg_est_base.cov_ridge*eye(p);
end
end

function [best_tau, tau_grid, ll_curve] = likelihood_scan_tau(OmEst, Sii_proxy, tau_scan, cfg_est)
% 对每个 tau（相对偏相干中位数）进行软阈收缩，计算源域对数似然并求和取最大
[p,~,F] = size(OmEst);
up = triu(true(p),1);
% 基准尺度：各频的偏相干中位数
meds = zeros(F,1);
for f=1:F
    S = pcoh_scores(OmEst(:,:,f)); v = abs(S(up));
    m = median(v);
    mad = median(abs(v - m)) / 0.6745;
    meds(f) = max(mad, 1e-6);  % 防止零尺度
end
tau_grid = tau_scan;
ll_curve = zeros(numel(tau_grid),1);
for it = 1:numel(tau_grid)
    tau_rel = tau_grid(it);
    ll_sum = 0;
    for f=1:F
        K0 = hermitize(OmEst(:,:,f));
        tau_eff = tau_rel * meds(f);
        K = shrink_offdiag(K0, tau_eff, cfg_est);
        Sf = hermitize(Sii_proxy(:,:,f));
        % 对数似然（忽略常数项）： log|K| - tr(S*K)
        ll_sum = ll_sum + logdet_spd(K) - real(trace(Sf*K));
    end
    ll_curve(it) = ll_sum;
end
[~,idx] = max(ll_curve); best_tau = tau_grid(idx);
end

function panel_files = plot_pcoh_panels(Omega_all, OmTrue, eval_cfg, outdir)
% 画偏相干面板图（每频一图，包含真值）
names = fieldnames(Omega_all);
[~,~,F] = size(OmTrue);
panel_files = cell(F,1);

for f = 1:F
    M = numel(names) + 1;  % +1 for ground truth
    Scells = cell(M,1);
    labels = cell(M,1);
    vmax = 0;
    
    % Ground truth first
    S = pcoh_scores(OmTrue(:,:,f));
    S(1:size(S,1)+1:end) = 0;
    v = abs(S(:)); v = v(~isnan(v));
    if ~isempty(v)
        q = quantile(v, eval_cfg.panel_clip);
        S = max(min(S, q), 0);
        vmax = max(vmax, q);
    end
    Scells{1} = S;
    labels{1} = 'Ground Truth';
    
    % Estimated methods
    for i=1:numel(names)
        K_raw = Omega_all.(names{i});
        % 关键：统一格式（cell 或 3D 数组）
        K_stack = coerce_precision_stack(K_raw);
        
        S = pcoh_scores(K_stack(:,:,f));
        S(1:size(S,1)+1:end) = 0;
        v = abs(S(:)); v = v(~isnan(v));
        if ~isempty(v)
            q = quantile(v, eval_cfg.panel_clip);
            S = max(min(S, q), 0);
            vmax = max(vmax, q);
        end
        Scells{i+1} = S;
        labels{i+1} = names{i};
    end
    
    if ~isempty(eval_cfg.panel_caxis), vmax = eval_cfg.panel_caxis(2); end

    % 绘制
    figure('Color','w','Name',sprintf('PCoh panel (f=%d)', f),'Position',[100 100 260*M 220]);
    tiledlayout(1,M,'Padding','compact','TileSpacing','compact');
    for i=1:M
        nexttile;
        imagesc(Scells{i});
        axis image off; 
        title(strrep(labels{i},'_','\_'),'Interpreter','none','FontSize',10);
        if ~isempty(vmax), caxis([0 vmax]); end
    end
    colormap(parula);
    cb = colorbar; cb.Layout.Tile = 'east'; cb.Label.String = '|partial coherence|';
    file_i = fullfile(outdir, sprintf('pcoh_panel_f%02d.png', f));
    saveas(gcf, file_i); 
    panel_files{f} = file_i;
end
end

function D = spd_distance_pack(Astack, Bstack)
% 计算多种 SPD 距离（逐频）
[~,~,F] = size(Astack);
D = struct();
D.KLS  = zeros(F,1);
D.LE   = zeros(F,1);
D.AIRM = zeros(F,1);
D.JBLD = zeros(F,1);
D.BW   = zeros(F,1);

for f=1:F
    A = hermitize(Astack(:,:,f)); A = spd_project_small(A);
    B = hermitize(Bstack(:,:,f)); B = spd_project_small(B);
    
    % KLS 对称
    D.KLS(f) = 0.5*( trace(B\A) + trace(A\B) - 2*size(A,1) );
    
    % Log-Euclidean
    LA = logm(A); LB = logm(B);
    D.LE(f) = norm(LA - LB, 'fro');
    
    % AIRM
    C = A^(-1/2) * B * A^(-1/2);
    C = spd_project_small(C);
    D.AIRM(f) = norm(logm(C), 'fro');
    
    % JBLD
    D.JBLD(f) = logdet_spd((A+B)/2) - 0.5*logdet_spd(A) - 0.5*logdet_spd(B);
    
    % Bures-Wasserstein
    R = A^(1/2) * B * A^(1/2);
    R = spd_project_small(R);
    D.BW(f) = sqrt( trace(A) + trace(B) - 2*trace(real(sqrtm(R))) );
end
end

function S = summarize_distances(D)
% 对每个距离给出 mean/median/std
fn = fieldnames(D);
S = struct();
for i=1:numel(fn)
    x = D.(fn{i});
    S.(fn{i}) = struct('mean',mean(x), 'median',median(x), 'std',std(x), ...
        'per_freq', x(:).');
end
end

function write_distance_csv(path, dist_summary)
% 导出距离汇总到 CSV
methods = fieldnames(dist_summary);
metrics = {'KLS','LE','AIRM','JBLD','BW'};
fid = fopen(path,'w');
fprintf(fid, 'method,metric,mean,median,std\n');
for i=1:numel(methods)
    m = methods{i};
    S = dist_summary.(m);
    for j=1:numel(metrics)
        mt = metrics{j};
        if isfield(S, mt)
            fprintf(fid, '%s,%s,%.6g,%.6g,%.6g\n', m, mt, ...
                S.(mt).mean, S.(mt).median, S.(mt).std);
        end
    end
end
fclose(fid);
end

%% ===== 补充缺失的函数 =====

function S = support_stats_mask(Mask)
% 从布尔掩码计算支持统计量
[p,~,F] = size(Mask);
up = triu(true(p),1);
S = struct();
edges_per_freq = zeros(F,1);
sup = false(nnz(up), F);
for f=1:F
    af = Mask(repmat(up,[1,1,F]));
    af = reshape(af, nnz(up), F);
    sup(:,f) = af(:,f);
    edges_per_freq(f) = nnz(af(:,f));
end
S.edges_per_freq = edges_per_freq(:).';
S.unique_edges = nnz(any(sup,2));
if F>=2
    J = zeros(F-1,1);
    for f=1:F-1
        a = sup(:,f); b = sup(:,f+1);
        J(f) = nnz(a & b) / max(nnz(a | b), 1);
    end
    S.jaccard_adj_min = min(J);
    S.jaccard_adj_med = median(J);
    S.jaccard_adj_max = max(J);
else
    [S.jaccard_adj_min, S.jaccard_adj_med, S.jaccard_adj_max] = deal(NaN);
end
end

%% ===== 基础工具函数 =====

function print_sim_summary(bayes_truth, sim, Omega_true)
[p,~,F] = size(Omega_true);
fprintf('Sim done: p=%d, m=%d, F=%d, T=%d\n', sim.p, getd(sim,'m',getd(sim,'n',NaN)), sim.F, sim.T);
fprintf('  Hyper-params (truth): lambda1*=%.4g, lambda2*=%.4g\n', ...
    getf(bayes_truth,'lambda1_true',NaN), getf(bayes_truth,'lambda2_true',NaN));
switch lower(bayes_truth.noise_type)
    case 'scalar'
        fprintf('  Noise (scalar): sigma_xi2=%.4g\n', bayes_truth.sigma_xi2_true);
    case 'matrix'
        Sx = hermitize(bayes_truth.Sigma_xixi_true);
        fprintf('  Noise (matrix): tr=%.4g, ||·||_F=%.4g, cond=%.3g\n', ...
            trace(Sx), norm(Sx,'fro'), cond(Sx));
    case 'matrix_per_freq'
        vals = zeros(F,3);
        for f=1:F
            Sx = hermitize(bayes_truth.Sigma_xixi_per_freq_true{f});
            vals(f,:) = [trace(Sx), norm(Sx,'fro'), cond(Sx)];
        end
        fprintf('  Noise (matrix_per_freq): tr[%.3g/%.3g/%.3g], cond[%.3g/%.3g/%.3g]\n', ...
            min(vals(:,1)), median(vals(:,1)), max(vals(:,1)), ...
            min(vals(:,3)), median(vals(:,3)), max(vals(:,3)));
end
gt = support_stats(Omega_true, 1e-12);
print_support_stats('GroundTruth', gt);
end

function S = support_stats(OmStack, tol)
[p,~,F] = size(OmStack);
mask = triu(true(p),1);
S = struct();
edges_per_freq = zeros(F,1);
sup = false(nnz(mask), F);
for f=1:F
    A = score_matrix(OmStack(:,:,f));
    af = A(mask) > tol;
    sup(:,f) = af;
    edges_per_freq(f) = nnz(af);
end
S.edges_per_freq = edges_per_freq(:).';
S.unique_edges = nnz(any(sup,2));
if F>=2
    J = zeros(F-1,1);
    for f=1:F-1
        a = sup(:,f); b = sup(:,f+1);
        J(f) = nnz(a & b) / max(nnz(a | b), 1);
    end
    S.jaccard_adj_min = min(J);
    S.jaccard_adj_med = median(J);
    S.jaccard_adj_max = max(J);
else
    [S.jaccard_adj_min,S.jaccard_adj_med,S.jaccard_adj_max] = deal(NaN);
end
end

function metrics = eval_support_curves(Omega_est_in, Omega_true_in, cfg)
num_thresh = getd(cfg,'num_thresh',128);
tol_true   = getd(cfg,'tol_true',1e-10);
Omega_est  = coerce_precision_stack(Omega_est_in);
Omega_true = coerce_precision_stack(Omega_true_in);
[p,~,F] = size(Omega_true);
mask = triu(true(p),1);
true_scores = []; est_scores = [];
for f=1:F
    Tm = score_matrix(Omega_true(:,:,f));
    Em = score_matrix(Omega_est(:,:,f));
    tvec = Tm(mask); evec = Em(mask);
    true_scores = [true_scores; tvec(:)]; %#ok<AGROW>
    est_scores  = [est_scores;  evec(:)]; %#ok<AGROW>
end
gt = true_scores > tol_true;
q = linspace(0,1,num_thresh);
ths = quantile(est_scores, q);
P=zeros(num_thresh,1); R=P; TPR=P; FPR=P; TP=P; FP=P; TN=P; FN=P;
Npos = sum(gt); Nneg = numel(gt) - Npos;
for k = 1:num_thresh
    pred = est_scores >= ths(k);
    TP(k) = sum( pred &  gt);
    FP(k) = sum( pred & ~gt);
    FN(k) = sum(~pred &  gt);
    TN(k) = sum(~pred & ~gt);
    P(k)  = TP(k) / max(TP(k)+FP(k), 1);
    R(k)  = TP(k) / max(TP(k)+FN(k), 1);
    TPR(k)= TP(k) / max(Npos,1);
    FPR(k)= FP(k) / max(Nneg,1);
end
[Rs, idxR] = sort(R); Ps = P(idxR); AP  = trapz(Rs, Ps);
[FPRs, idxF] = sort(FPR); TPRs = TPR(idxF); AUC = trapz(FPRs, TPRs);
F1 = 2*P.*R ./ max(P+R, eps); [bestF1, ib] = max(F1);
best_struct = struct();
best_struct.F1  = bestF1;
best_struct.P   = P(ib);
best_struct.R   = R(ib);
best_struct.thr = ths(ib);
counts_struct = struct();
counts_struct.TP   = TP; counts_struct.FP   = FP;
counts_struct.FN   = FN; counts_struct.TN   = TN;
counts_struct.Npos = Npos; counts_struct.Nneg = Nneg;
metrics = struct();
metrics.thresholds = ths; metrics.P   = P; metrics.R   = R; metrics.AP  = AP;
metrics.FPR = FPR; metrics.TPR = TPR; metrics.AUC = AUC;
metrics.bestF1 = best_struct; metrics.counts = counts_struct;
end

function S = score_matrix(Om)
Om = hermitize(Om);
d  = real(diag(Om)); d(d<=0) = eps;
s  = sqrt(d);
D12 = s * s.';
S = abs(Om) ./ max(D12, eps);
S(1:size(S,1)+1:end) = 0;
S = hermitize(S);
end

function X3 = coerce_precision_stack(X)
if iscell(X)
    F = numel(X); p = size(X{1},1);
    X3 = zeros(p,p,F,class(X{1}));
    for f = 1:F, Xi = hermitize(X{f}); X3(:,:,f) = Xi; end
elseif isnumeric(X) && ndims(X)==3 && size(X,1)==size(X,2)
    X3 = X;
    for f=1:size(X3,3), X3(:,:,f) = hermitize(X3(:,:,f)); end
else
    error('Unsupported precision format.');
end
end

function W = inv_chol(S)
S = hermitize(S);
[V,D] = eig(S,'vector'); D = real(D); D(D<1e-12) = 1e-12;
S = V*diag(D)*V';
C = chol(S,'lower');
W = C \ eye(size(S));
end

function K = shrink_offdiag(K, tau, cfg)
K = hermitize(K);
p = size(K,1);
d = real(diag(K)); d(d<=0)=eps;
D12 = sqrt(d)*sqrt(d).';
Mag = abs(K);
Phs = complex(real(K), imag(K)) ./ max(Mag, eps);
Kshr = K;
idx = triu(true(p),1);
v = Mag(idx) - tau * D12(idx);
v = max(v, 0);
Kshr(idx) = Phs(idx) .* v;
Kshr = hermitize(Kshr);
Kshr = spd_project_small(Kshr) + cfg.spd_delta*eye(p);
end

function S = inv_psd(K)
K = hermitize(K);
[V,D] = eig(K,'vector'); D = real(D);
D(D < 1e-9*median(abs(D)+eps)) = 1e-9*median(abs(D)+eps);
S = V*diag(1./D)*V'; S = hermitize(S);
end

function A = hermitize(A)
A = (A + A')/2;
end

function K = spd_project_small(K)
K = hermitize(K);
[V,D] = eig(K,'vector');
D = real(D); 
floorv = max(1e-9, 1e-9*median(abs(D)+eps));
D(D<floorv) = floorv;
K = V*diag(D)*V';
K = hermitize(K);
end

function s = logdet_spd(A)
A = spd_project_small(A);
try
    C = chol(hermitize(A), 'lower');
    s = 2*sum(log(diag(C)));
catch
    [V,D] = eig(hermitize(A), 'vector');
    D = real(D); D(D<1e-12)=1e-12;
    s = sum(log(D));
end
end

function S = pcoh_scores(K)
K = hermitize(K);
d = real(diag(K)); d(d<=0) = eps;
sc = sqrt(d);
S = abs(K) ./ max(sc*sc.', eps);
S(1:size(S,1)+1:end) = 0;
S = hermitize(S);
end

function print_support_stats(tag, s)
fprintf('[%s] edges_per_freq=%s | unique=%d | Jaccard [%.2f/%.2f/%.2f]\n', ...
    tag, mat2str(s.edges_per_freq), s.unique_edges, ...
    s.jaccard_adj_min, s.jaccard_adj_med, s.jaccard_adj_max);
end

function A = override_struct(A, B)
if ~isstruct(B), return; end
fb = fieldnames(B);
for k=1:numel(fb), A.(fb{k}) = B.(fb{k}); end
end

function val = set_default(S, D)
val = S; f = fieldnames(D);
for i=1:numel(f)
    if ~isfield(val,f{i}) || isempty(val.(f{i}))
        val.(f{i}) = D.(f{i}); 
    end
end
end

function val = getd(s, name, def)
if isstruct(s) && isfield(s,name) && ~isempty(s.(name))
    val = s.(name); 
else
    val = def; 
end
end

function z = getf(s, f, def)
if isfield(s,f) && ~isempty(s.(f))
    z = s.(f); 
else
    z = def; 
end
end

function nm2 = make_valid_name(nm)
if exist('matlab.lang.makeValidName','file')
    nm2 = matlab.lang.makeValidName(nm);
else
    nm2 = regexprep(nm,'[^a-zA-Z0-9_]','_');
    if isempty(nm2) || (~isletter(nm2(1)) && nm2(1)~='_')
        nm2 = ['x_' nm2];
    end
end
end

function nm_u = uniquify_field(nm, existing)
nm_u = nm; c = 2;
while any(strcmp(existing, nm_u))
    nm_u = sprintf('%s_%d', nm, c); c = c + 1;
end
end

%% ===== PR/ROC 相关函数（仅在启用时使用） =====

function print_method_summary(name, m, s)
fprintf('\n--- %s ---\n', name);
fprintf('PR:   AP=%.3f | ROC: AUC=%.3f\n', m.AP, m.AUC);
fprintf('BestF1=%.3f @thr=%.3g (P=%.3f, R=%.3f)\n', ...
    m.bestF1.F1, m.bestF1.thr, m.bestF1.P, m.bestF1.R);
if isscalar(s.thr_match)
    fprintf('Match-sparsity: thr=%.3g | F1=%.3f (P=%.3f, R=%.3f)\n', ...
        s.thr_match, s.F1_match, s.P_match, s.R_match);
else
    tmn = min(s.thr_match); tmd = median(s.thr_match); tmx = max(s.thr_match);
    fprintf('Match-sparsity: thr[%.3g/%.3g/%.3g] | F1=%.3f (P=%.3f, R=%.3f)\n', ...
        tmn, tmd, tmx, s.F1_match, s.P_match, s.R_match);
end
fprintf('Rel Fro offdiag: %.3f\n', s.relF_off);
end

function Mask = threshold_support_per_freq(OmEst, OmTrue, tol_true)
[p,~,F] = size(OmTrue);
Mask = false(p,p,F);
up = triu(true(p),1);
Npos_f = zeros(F,1);
for f=1:F
    Tm = score_matrix(OmTrue(:,:,f));
    Npos_f(f) = nnz(Tm(up) > tol_true);
end
for f=1:F
    Em = score_matrix(OmEst(:,:,f));
    sc = Em(up);
    sel = false(size(sc));
    if ~isempty(sc) && Npos_f(f)>0
        [~, ord] = sort(sc, 'descend');
        k = min(Npos_f(f), numel(ord));
        sel(ord(1:k)) = true;
    end
    B = false(p); B(up) = sel; B = B | B.'; 
    Mask(:,:,f) = B;
end
end

function [thr_per_f, F1, P, R] = threshold_matching_sparsity_per_freq(OmEst, OmTrue, tol_true)
[p,~,F] = size(OmTrue);
mask = triu(true(p),1);
Npos_f = zeros(F,1);
for f=1:F
    Tm = score_matrix(OmTrue(:,:,f)); 
    Npos_f(f) = nnz(Tm(mask) > tol_true); 
end
thr_per_f = zeros(F,1); TP=0; FP=0; FN=0;
for f=1:F
    Em = score_matrix(OmEst(:,:,f)); 
    Tm = score_matrix(OmTrue(:,:,f));
    gt = Tm(mask) > tol_true; 
    scores = Em(mask);
    if isempty(scores) || Npos_f(f)==0
        thr_per_f(f) = Inf; 
        pr = false(size(scores));
    else
        [scores_sorted, idx_sorted] = sort(scores, 'descend');
        k = min(max(Npos_f(f),1), numel(scores_sorted));
        topk_idx = idx_sorted(1:k); 
        pr = false(size(scores)); 
        pr(topk_idx) = true;
        thr_per_f(f) = scores_sorted(k);
    end
    TP = TP + sum(pr & gt); 
    FP = FP + sum(pr & ~gt); 
    FN = FN + sum(~pr & gt);
end
P = TP / max(TP+FP,1); 
R = TP / max(TP+FN,1);
F1 = 2*P*R / max(P+R, eps);
end

function rel = rel_frob_offdiag(OmEst, OmTrue)
[p,~,F] = size(OmTrue);
mask = triu(true(p),1); 
num=0; den=0;
for f=1:F
    A = score_matrix(OmEst(:,:,f)); 
    B = score_matrix(OmTrue(:,:,f));
    ae = A(mask); 
    be = B(mask);
    num = num + norm(ae - be)^2; 
    den = den + norm(be)^2 + eps;
end
rel = sqrt(num/den);
end

function plot_eval_curves_multi(metrics, outdir)
names = fieldnames(metrics);
% PR
figure('Name','PR (all methods)','Color','w'); 
hold on; grid on;
for i=1:numel(names)
    m = metrics.(names{i}); 
    [Rs, idxR] = sort(m.R); 
    Ps = m.P(idxR);
    plot(Rs, Ps, 'LineWidth', 1.8, ...
        'DisplayName', sprintf('%s (AP=%.3f)',names{i}, m.AP));
end
xlabel('Recall'); ylabel('Precision'); 
xlim([0 1]); ylim([0 1]);
legend('Location','southwest'); 
title('Precision–Recall');
saveas(gcf, fullfile(outdir, 'PR_all.png')); 
close(gcf);

% ROC
figure('Name','ROC (all methods)','Color','w'); 
hold on; grid on; axis square;
plot([0 1],[0 1],'k--','HandleVisibility','off');
for i=1:numel(names)
    m = metrics.(names{i}); 
    [FPRs, idxF] = sort(m.FPR); 
    TPRs = m.TPR(idxF);
    plot(FPRs, TPRs, 'LineWidth', 1.8, ...
        'DisplayName', sprintf('%s (AUC=%.3f)',names{i}, m.AUC));
end
xlabel('FPR'); ylabel('TPR'); 
xlim([0 1]); ylim([0 1]);
legend('Location','southeast'); 
title('ROC');
saveas(gcf, fullfile(outdir, 'ROC_all.png')); 
close(gcf);

% BestF1
figure('Name','Best F1 (all methods)','Color','w');
vals = zeros(numel(names),1); 
for i=1:numel(names)
    vals(i) = metrics.(names{i}).bestF1.F1; 
end
bar(vals); grid on; ylim([0 1]);
set(gca,'XTick',1:numel(names),'XTickLabel',names,'XTickLabelRotation',20);
ylabel('Best F1'); 
title('Best F1 per method');
saveas(gcf, fullfile(outdir, 'BestF1_bar.png')); 
close(gcf);
end