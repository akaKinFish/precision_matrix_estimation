function result = Main_sim2_multi(output_source, cfg)
%% Main_sim2_multi — Sim2 的多频（F>1）仿真 + 多方法对比（含我的跨频PGD）
%
% 生成思路（源域真值）：
%   1) 用 gen_hggm2 先得到一个“基准”精度 Θ0（块结构 + 复数）。
%   2) 在频轴上为非零边生成平滑的幅值曲线（少量边允许局部增/减/翻转），得到 {Θ(f)}。
%   3) 对每个频点做 Hermitian + SPD 投影，Σ(f)=Θ(f)^{-1}，采样 J(f)，投影/加噪得到 Svv(f)。
%
% 方法对比：
%   (1) H-HGGM / EM：penalty = lasso (1), ridge (2), naive (0)【逐频独立】
%   (2) eLORETA + hg-lasso【逐频独立】
%   (3) LCMV + hg-lasso【逐频独立】
%   (4) VARETA【逐频独立】
%   (5) 我的方法 my_test2_adapter_bcvE（E=BC-V，M=PGD，跨频Kernel + 支撑再拟合 + recolor）
%
% 评估：逐频 SPD 距离 + AUC/SEN/SPEC/PREC/F1；提供若干可视化（可关）
%
% 依赖（均为你仓库已有）：
%   gen_hggm2.m, higgs.m, eloreta_hg_lasso.m, lcmv_hg_lasso.m, vareta.m
%   my_test2_adapter_bcvE.m, compute_spd_and_classification_metrics.m
%
% 作者：你我协作版本（2025-10）

if nargin<1 || isempty(output_source), output_source = pwd; end
if nargin<2, cfg = struct(); end

%% -------------------------
%% 全局配置（可按需改）
%% -------------------------
% 频轴与真值平滑
F           = getf(cfg,'F',9);               % 频点数（≥2 才能体现跨频）
sigma_f     = getf(cfg,'sigma_f',3.0);       % 频轴平滑核宽 ~ 2~3
flip_rate   = getf(cfg,'flip_rate',0.02);    % 每个频点可变边翻转率（少量差异）
vary_frac   = getf(cfg,'vary_frac',0.20);    % 在非零边中，允许随频变化的比例
amp_cv      = getf(cfg,'amp_cv',0.25);       % 幅值相对变动（coefficient of variation）
phase_cv    = getf(cfg,'phase_cv',0.10);     % 相位扰动幅度（弧度的比例尺度）
use_phase   = getf(cfg,'use_phase',false);   % 默认只变幅值，若需也可轻度变相位

% 样本/维度/噪声
m           = getf(cfg,'m',600);             % 每频样本数
q           = getf(cfg,'q',22);              % 源数
p           = getf(cfg,'p',30);              % 传感器数
options     = getf(cfg,'options',struct('config',2,'var',2, ...
                 'extensions',[ceil(q/3); ceil(q/3); q-2*ceil(q/3)], ...
                 'connections',[1 2; 2 3]));
bio_noise   = getf(cfg,'bio_noise',0.10);    % 生物噪声比例（投影后）
sens_noise  = getf(cfg,'sens_noise',0.10);   % 传感器噪声比例
noise_mod   = getf(cfg,'noise_freq_mod',0.15); % 噪声随频缓变幅度（0~0.3）

% 可视化
do_panel_maps   = getf(cfg,'plot_panel_maps',true);   % 小 multiples 面板（可关）
do_edge_tracks  = getf(cfg,'plot_edge_tracks',true);  % 顶K边随频曲线（可关）
topK_edges      = getf(cfg,'topK_edges',20);

% 我的适配器配置（与单频版一致）
cfg_my = struct();
cfg_my.output_format    = 'cell';
cfg_my.verbose          = false;
cfg_my.do_support_refit = true;
cfg_my.warm_method      = 'eloreta';
cfg_my.kernel_sigma     = getf(cfg,'kernel_sigma',sigma_f); % 与真值平滑同量级
cfg_my.lambda3_ratio    = getf(cfg,'lambda3_ratio',0.3);
cfg_my.use_subspace     = true;
cfg_my.noise_model      = 'scalar';
cfg_my.do_scale_L       = true;
cfg_my.do_scale_data    = true;

%% 输出目录
out_dir = fullfile(output_source,'Sim2_multi');
if ~isfolder(out_dir), mkdir(out_dir); end

%% -------------------------
%% 伪 leadfield（同圆布点）
%% -------------------------
Lvj = build_pseudo_leadfield(p,q);

%% -------------------------
%% 1) 基准源域真值（单频生成 Θ0）
%% -------------------------
[~, j_sim0, Theta0] = gen_hggm2(m, q, options); %#ok<ASGLU>
Theta0 = 0.5*(Theta0 + Theta0');  % Hermitian

% 非零边集合（去对角）
M0 = abs(Theta0)>0; M0(1:q+1:end)=false;
idx_all = find(M0);
n_all = numel(idx_all); 
n_vary = max(1, round(vary_frac*n_all));
rng('shuffle');
idx_vary = randsample(idx_all, n_vary);
% 频轴平滑核（用于生成平滑扰动）
K = make_frequency_kernel(F, sigma_f); % 行和≈1

%% -------------------------
%% 2) 频轴上生成 {Theta(f), Sigma(f)} 与样本/传感器数据
%% -------------------------
Theta_true  = cell(F,1);
Sigma_true  = cell(F,1);
J_cell      = cell(F,1);
Sjj_true    = cell(F,1);
Svv_cell    = cell(F,1);

% 生成随频平滑的“幅值比例”轨迹（对 idx_vary）
scale_trajs = gen_smooth_scales(length(idx_vary), F, amp_cv, K);
% 生成轻度相位扰动（若启用）
phase_trajs = gen_smooth_scales(length(idx_vary), F, phase_cv, K);

for f = 1:F
    % 基于 Θ0 构造 Θ(f)：对 idx_vary 的边做平滑幅值（与可选相位）扰动
    Tf = Theta0;
    [I_v,J_v] = ind2sub([q,q], idx_vary);
    for k = 1:length(idx_vary)
        i = I_v(k); j = J_v(k);
        base = Theta0(i,j);
        mag  = abs(base);
        ph   = angle(base);
        mag_new = mag * (1 + scale_trajs(k,f));
        if use_phase
            ph = ph + phase_trajs(k,f);
        end
        Tf(i,j) = mag_new * exp(1i*ph);
        Tf(j,i) = conj(Tf(i,j));
    end
    % 随机极少量“翻转/断边/新边”
    Tf = random_edge_flips(Tf, flip_rate, idx_all, q);

    % Hermitian + SPD 投影
    Tf = 0.5*(Tf + Tf');
    Tf = project_spd(Tf, 1e-8, 1e-12);
    Theta_true{f} = Tf;

    % 协方差与样本（复同构）
    Sf = inv_psd_robust(Tf, 1e-8, 1e-12);
    Sigma_true{f} = Sf;
    Jf = sample_complex_gaussian(Sf, m);   % q×m
    Sjj_true{f} = (Jf*Jf')/m; Sjj_true{f} = 0.5*(Sjj_true{f}+Sjj_true{f}');
    J_cell{f} = Jf;

    % 传感器域：v=L*J + 生物噪声(投影后) + 传感器噪声
    v0   = Lvj * Jf; % p×m
    % 频率缓变噪声比例
    modf = 1 + noise_mod*cos(2*pi*(f-1)/max(1,F-1));
    bio  = randn(q,m) + 1i*randn(q,m);
    bio  = Lvj*bio;
    bio  = norm(v0(:))*bio/max(norm(bio(:)),eps);
    sens = randn(p,m) + 1i*randn(p,m);
    sens = norm(v0(:))*sens/max(norm(sens(:)),eps);
    v    = v0 + (bio_noise*modf)*bio + (sens_noise*modf)*sens;

    Svv  = (v*v')/m; Svv = 0.5*(Svv+Svv');
    Svv_cell{f} = Svv;
end

%% -------------------------
%% 3) 基线方法（逐频独立） + 我的跨频方法
%% -------------------------
% 公共参数
param               = struct();
param.use_gpu       = 0;
param.run_bash_mode = 0;
param.str_band      = "none";
param.maxiter_outer = 60;
param.maxiter_inner = 30;
param.p             = p;
param.q             = q;
param.Ip            = eye(p);
param.Op            = ones(p,1);
param.Iq            = eye(q);
param.eigreg        = 1E-4;
param.method        = 'lqa';
aj                  = sqrt(log(q)/m);
param.aj            = aj;
param.Ajj           = (ones(q)-eye(q));  % 与单频 Sim2 口径一致
param.axi           = 1E-4;
param.Axixi         = eye(p);
param.Axixi_inv     = eye(p);
param.ntry          = 0;
param.prew          = 1;
param.rth1          = 0.7;
param.rth2          = 3.16;
param.m = m;
param.nu            = m;


% 容器
method_names = {'higgs_lasso','higgs_ridge','higgs_naive','eloreta_hglasso','lcmv_hglasso','vareta','my_adapter'};
n_methods = numel(method_names);
Theta_est = cell(n_methods,1);
Sigma_est = cell(n_methods,1);
for k=1:n_methods, Theta_est{k} = cell(F,1); Sigma_est{k} = cell(F,1); end

% ---- HIGGS 三种惩罚、eLORETA、LCMV、VARETA（逐频独立） ----
for f = 1:F
    Svv = Svv_cell{f};

    % HIGGS：penalty = lasso(1), ridge(2), naive(0)
    penalty_vec = [1 2 0];
    for pen = 1:numel(penalty_vec)
        param.penalty = penalty_vec(pen);
        [Th, Sj] = higgs(Svv, Lvj, param);  % Sj: 源域协方差, Th: 源域精度
        Theta_est{pen}{f} = 0.5*(Th+Th');
        Sigma_est{pen}{f} = 0.5*(Sj+Sj');
    end

    % eLORETA + hg-lasso
    param.gamma1 = 0.001; param.gamma2 = 0.05; param.delta_gamma = 0.001;
    [Th4, Sj4] = eloreta_hg_lasso(Svv, Lvj, param);
    Theta_est{4}{f} = 0.5*(Th4+Th4');  Sigma_est{4}{f} = 0.5*(Sj4+Sj4');

    % LCMV + hg-lasso
    param.gamma = sum(abs(diag(Svv)))/(length(Svv)*100);
    [Th5, Sj5] = lcmv_hg_lasso(Svv, Lvj, param);
    Theta_est{5}{f} = 0.5*(Th5+Th5');  Sigma_est{5}{f} = 0.5*(Sj5+Sj5');

    % VARETA（只给 Σ，需稳健逆得到 Θ）
    [U, Sv, V] = svd(Lvj, 'econ'); svals = diag(Sv);
    Sj6 = vareta(U, svals, V, Svv, 0);
    Sigma_est{6}{f} = 0.5*(Sj6+Sj6');
    S6 = Sigma_est{6}{f};
    Theta_est{6}{f} = inv_psd_robust(S6, 1e-8, 1e-12);
end


% ---- 我的跨频 PGD 方法（一次性喂全频 cell） ----
[Omega_cell, Dsrc_cell, Gamma_cell, outs_my] = my_test2_adapter_bcvE(Svv_cell, Lvj, m, cfg_my); %#ok<ASGLU>
for f=1:F
    Theta_est{7}{f} = 0.5*(Omega_cell{f}+Omega_cell{f}');
    Sigma_est{7}{f} = 0.5*(Dsrc_cell{f}+Dsrc_cell{f}');
    % 自检：配对是否自洽
    cond_err = norm(Theta_est{7}{f}*Sigma_est{7}{f} - eye(q),'fro')/max(1,q);
    if cond_err>1e-2
        warning('[my-adapter f=%d] ||ΘΣ-I||_F/q = %.3g (>1e-2)', f, cond_err);
    end
end

%% -------------------------
%% 4) 逐频指标（SPD + 分类 + Rayleigh）
%% -------------------------
opts_eval = struct();
opts_eval.normalize_mode = 'maxabs';
opts_eval.spd = struct('symmetrize',true,'project',true,'eps',1e-8);
opts_eval.alpha = struct('value',0,'use_bcv_bug',false);
opts_eval.plot = struct('radar',false,'title','');

eval_results = struct();
for im=1:n_methods
    eval_results.(method_names{im}) = cell(F,1);
end

for f=1:F
    Theta_true_f = Theta_true{f};
    for im=1:n_methods
        T_est = 0.5*(Theta_est{im}{f} + Theta_est{im}{f});
        S_est = 0.5*(Sigma_est{im}{f} + Sigma_est{im}{f});

        % 解析式去偏
        T_unb = 2*T_est - T_est*S_est*T_est;

        % Rayleigh 矫正（rth=3.16）
        dv = abs(diag(T_est));
        Tvar = sqrt(dv*dv.' + abs(T_est).^2);
        rth = 3.16;
        T_ray = T_unb;
        mask_ray = abs(T_ray) < (rth/sqrt(m))*(Tvar - diag(diag(Tvar)));
        T_ray(mask_ray) = 0;

        args_eval = struct();
        args_eval.Theta_true = Theta_true_f;
        args_eval.Theta_est  = T_est;
        args_eval.Theta_unb  = T_unb;
        args_eval.Theta_ray  = T_ray;
        args_eval.is_complex = true;

        out_eval = compute_spd_and_classification_metrics(args_eval, opts_eval);
        eval_results.(method_names{im}){f} = out_eval;
    end
end

%% -------------------------
%% 5) 可视化（可关）
%% -------------------------
if do_panel_maps
    % 画：真值 PCoh & 我的方法 / HIGGS-L1 在所有频点的小面板（示意）
    fig1 = figure('Position',[120,80,1200,180+180*ceil(F/3)]);
    try, load('colormap1'); colormap(cmap); catch, colormap('hot'); end
    nrow = 3; ncol = ceil(F);
    for f=1:F
        % 真值 PCoh（用 Σ_true 计算）
        X = Sjj_true{f}; X = X - diag(diag(X)); X = X/max(abs(X(:))+eps);
        subplot(nrow,F,f); imagesc(abs(X)); title(sprintf('True PCoh f=%d',f));
        % HIGGS-L1
        X2 = Sigma_est{1}{f}; X2 = X2 - diag(diag(X2)); X2 = X2/max(abs(X2(:))+eps);
        subplot(nrow,F,F+f); imagesc(abs(X2)); title(sprintf('HIGGS-L1 PCoh f=%d',f));
        % 我的方法
        X3 = Sigma_est{7}{f}; X3 = X3 - diag(diag(X3)); X3 = X3/max(abs(X3(:))+eps);
        subplot(nrow,F,2*F+f); imagesc(abs(X3)); title(sprintf('My-Adapter PCoh f=%d',f));
    end
    saveas(fig1, fullfile(out_dir,'panel_pcoh_true_higgs_my.fig'));
    close(fig1);
end

if do_edge_tracks
    % 选真值强度 topK 边，画“边权随频”曲线（真值 vs 我的 / HIGGS-L1）
    [edge_list, weights] = topK_edges_from_true(Theta_true, topK_edges);
    fig2 = figure('Position',[120,80,1200,600]); tiledlayout('flow');
    for k=1:size(edge_list,1)
        i=edge_list(k,1); j=edge_list(k,2);
        w_true = zeros(1,F); w_h1=zeros(1,F); w_my=zeros(1,F);
        for f=1:F
            w_true(f) = abs(Theta_true{f}(i,j));
            w_h1(f)   = abs(Theta_est{1}{f}(i,j));
            w_my(f)   = abs(Theta_est{7}{f}(i,j));
        end
        nexttile; plot(1:F, w_true,'-o'); hold on; plot(1:F, w_h1,'--x'); plot(1:F, w_my,'-.*'); hold off;
        xlabel('freq idx'); ylabel('|edge|'); title(sprintf('edge (%d,%d)',i,j));
        legend({'true','higgs-l1','my'},'Location','best');
    end
    saveas(fig2, fullfile(out_dir,'edge_tracks_true_higgs_my.fig'));
    close(fig2);
end

%% -------------------------
%% 5.5  跨频指标汇总并导出（确保包含 my_adapter）
%% -------------------------
% Rayleigh 支撑（与评估口径一致）
rth = 3.16;
Theta_ray = cell(numel(method_names),1);
for im = 1:numel(method_names)
    Theta_ray{im} = cell(F,1);
    for f = 1:F
        T = 0.5*(Theta_est{im}{f}+Theta_est{im}{f}');     % Hermitian
        S = 0.5*(Sigma_est{im}{f}+Sigma_est{im}{f}');
        % 解析式去偏 + Rayleigh 阈值
        Tunb = 2*T - T*S*T;
        dv   = abs(diag(T));
        Tvar = sqrt(dv*dv.' + abs(T).^2);
        Tray = Tunb;
        mask = abs(Tray) < (rth/sqrt(m))*(Tvar - diag(diag(Tvar)));
        Tray(mask) = 0;
        Theta_ray{im}{f} = 0.5*(Tray+Tray');
    end
end

% 真实频轴总变差 TV_true
TV_true = 0;
for f = 2:F
    TV_true = TV_true + norm(Theta_true{f} - Theta_true{f-1}, 'fro');
end
TV_true = max(TV_true, eps);  % 防零保护

% 选取真值 across-f 的 Top-K 边用于轨迹相关
K_for_tracks = min(20, q);    % 你可按需调整
[edge_list, ~] = topK_edges_from_true(Theta_true, K_for_tracks);

% 若逐频评估里带有 cls.AUC / cls.F1，这里自动求均值
has_AUC = false; has_F1 = false;
try
    tmp = eval_results.(method_names{1}){1};
    has_AUC = isfield(tmp,'cls') && isfield(tmp.cls,'AUC');
    has_F1  = isfield(tmp,'cls') && isfield(tmp.cls,'F1');
catch
end

% 汇总各方法
summary_rows = cell(0,1);
for im = 1:numel(method_names)
    % TV_ratio
    TV_est = 0;
    for f = 2:F
        TV_est = TV_est + norm(Theta_est{im}{f} - Theta_est{im}{f-1}, 'fro');
    end
    TV_ratio = TV_est / TV_true;

    % 邻频 Jaccard（基于 Rayleigh 支撑）
    Jacc = nan(F-1,1);
    for f = 2:F
        S1 = abs(Theta_ray{im}{f})   > 0; S1(1:q+1:end) = false;
        S0 = abs(Theta_ray{im}{f-1}) > 0; S0(1:q+1:end) = false;
        inter = nnz(S1 & S0); uni = nnz(S1 | S0);
        if uni>0, Jacc(f-1) = inter/uni; end
    end
    Jacc_mean = mean(Jacc,'omitnan');

    % Top-K 边轨迹相关
    C = nan(size(edge_list,1),1);
    for k = 1:size(edge_list,1)
        i = edge_list(k,1); j = edge_list(k,2);
        t_true = zeros(F,1); t_est = zeros(F,1);
        for f = 1:F
            t_true(f) = abs(Theta_true{f}(i,j));
            t_est(f)  = abs(Theta_est{im}{f}(i,j));
        end
        if std(t_true)>0 && std(t_est)>0
            cc = corr(t_true, t_est, 'type','Pearson', 'rows','complete');
            C(k) = cc;
        end
    end
    EdgeTrackCorr_mean = mean(C,'omitnan');

    % AUC/F1 的均值（若存在）
    AUC_mean = NaN; F1_mean = NaN;
    if has_AUC
        aucv = nan(F,1);
        for f=1:F
            try
                aucv(f) = eval_results.(method_names{im}){f}.cls.AUC;
            catch, aucv(f) = NaN;
            end
        end
        AUC_mean = mean(aucv,'omitnan');
    end
    if has_F1
        f1v = nan(F,1);
        for f=1:F
            try
                f1v(f) = eval_results.(method_names{im}){f}.cls.F1;
            catch, f1v(f) = NaN;
            end
        end
        F1_mean = mean(f1v,'omitnan');
    end

    summary_rows{end+1,1} = table( string(method_names{im}), TV_ratio, Jacc_mean, ...
        EdgeTrackCorr_mean, AUC_mean, F1_mean, ...
        'VariableNames', {'method','TV_ratio','Jaccard_mean','EdgeTrackCorr_mean','AUC_mean','F1_mean'});
end
summary_tbl = vertcat(summary_rows{:});

% 导出 CSV（确保覆盖旧文件，且一定包含 my_adapter 这一行）
out_dir = fullfile(output_source,'Sim2_multi');  % 若变量名不同请对齐
if ~isfolder(out_dir), mkdir(out_dir); end
writetable(summary_tbl, fullfile(out_dir,'summary_metrics.csv'));

%（可选）也塞进返回结构
result_summary = summary_tbl;  %#ok<NASGU>


%% -------------------------
%% 6) 汇总并返回
%% -------------------------
result = struct();
result.out_dir     = out_dir;
result.F           = F;
result.Lvj         = Lvj;
result.m           = m;
result.q           = q;
result.p           = p;
result.K           = K;
result.Theta_true  = Theta_true;
result.Sigma_true  = Sigma_true;
result.J_cell      = J_cell;
result.Svv_cell    = Svv_cell;
result.methods     = method_names;
result.Theta_est   = Theta_est;
result.Sigma_est   = Sigma_est;
result.eval        = eval_results;
result.cfg         = cfg;
result.cfg_my      = cfg_my;

% 也可以在此处添加“跨频平滑指标”（TV、邻频Jaccard等）
% 例如：result.crossfreq = compute_crossfreq_summaries(Theta_true, Theta_est);

end % === 主函数结束 ===


%% ================= helpers =================
function L = build_pseudo_leadfield(p,q)
    L   = zeros(p,q);
    radj  = 60;   radv = 85;
    angj  = 2*pi/q; angv = 2*pi/p;
    for contv = 1:p
        for contj = 1:q
            vectv            = [radv*cos((contv-1)*angv); radv*sin((contv-1)*angv)];
            vectj            = [radj*cos((contj-1)*angj); radj*sin((contj-1)*angj)];
            r                = vectv - vectj;
            r_unit           = r/sqrt(sum(abs(r).^2));
            miu              = vectj/sqrt(sum(abs(vectj).^2));
            L(contv,contj)   = (1/(4*pi))*miu'*r_unit/sqrt(sum(abs(r).^2))^2;
        end
    end
end

function K = make_frequency_kernel(F, sigma)
    [I,J]=ndgrid(1:F,1:F);
    K = exp(-((I-J).^2)/(2*sigma^2));
    K = 0.5*(K+K');   % 对称
    K = K / max(1,max(sum(K,2))); % 行和归一
end

function S = inv_psd_robust(Theta, eps_reg, min_ratio)
    Theta = 0.5*(Theta+Theta');
    [U,D]=eig(full(Theta),'vector'); d=real(D);
    dmax = max(d);
    floor_val = max(min_ratio*max(dmax,eps), 0);
    d(d<floor_val) = floor_val;
    if eps_reg>0, d = (d + eps_reg*dmax)/(1+eps_reg); end
    S = U*diag(1./d)*U'; S = 0.5*(S+S');
end

function T = project_spd(T, eps_reg, min_ratio)
    T = 0.5*(T+T');
    [U,D]=eig(full(T),'vector'); d=real(D);
    dmax = max(d);
    floor_val = max(min_ratio*max(dmax,eps), 0);
    d(d<floor_val) = floor_val;
    if eps_reg>0, d = (d + eps_reg*dmax)/(1+eps_reg); end
    T = U*diag(d)*U'; T = 0.5*(T+T');
end

function J = sample_complex_gaussian(Sigma, m)
    % 复同构采样：Wisomph = [Re Σ, -Im Σ; Im Σ, Re Σ]
    q = size(Sigma,1);
    Wisomph = [real(Sigma) -imag(Sigma); imag(Sigma) real(Sigma)];
    X = mvnrnd(zeros(1,2*q), Wisomph, m);  % m×(2q)
    DataRe = X(:,1:q); DataIm = X(:,q+1:2*q);
    J = (DataRe + 1i*DataIm).';            % q×m
end

function scale = gen_smooth_scales(n_edge, F, cv, K)
    % 生成零均值的平滑尺度扰动，每列对应一个频点
    if cv<=0, scale = zeros(n_edge,F); return; end
    raw = randn(n_edge,F);
    smooth = raw*K; % 在频轴上平滑
    % 标准化到目标 CV（相对幅度）
    for e=1:n_edge
        s = smooth(e,:);
        s = s - mean(s);
        if std(s)>0, s = s/std(s); end
        smooth(e,:) = (cv)*s;
    end
    scale = smooth;
end

function Tnew = random_edge_flips(T, flip_rate, idx_all, q)
    % 在极少量边上做“断开/新连”或小扰动，以制造少量结构差异
    if flip_rate<=0, Tnew=T; return; end
    nflip = binornd(numel(idx_all), flip_rate);
    if nflip==0, Tnew=T; return; end
    sel = randsample(idx_all, nflip);
    Tnew = T;
    for k=1:numel(sel)
        [i,j] = ind2sub([q,q], sel(k));
        if abs(Tnew(i,j))>0
            % 断开（置更小）
            Tnew(i,j)=0; Tnew(j,i)=0;
        else
            % 新连一条小边
            ph = 2*pi*rand;
            Tnew(i,j)= 0.15*exp(1i*ph); % 小幅新边
            Tnew(j,i)= conj(Tnew(i,j));
        end
    end
end

function [edge_list, weights] = topK_edges_from_true(Theta_true, K)
    % 以 across-f 的平均 |Θ| 选 TopK 非零边
    q = size(Theta_true{1},1);
    acc = zeros(q);
    for f=1:numel(Theta_true)
        acc = acc + abs(Theta_true{f});
    end
    acc = acc/numel(Theta_true);
    acc(1:q+1:end)=0;
    [vals,idx] = sort(acc(:),'descend');
    keep = idx(vals>0);
    K = min(K, numel(keep));
    idxK = keep(1:K);
    [I,J] = ind2sub([q,q], idxK);
    edge_list = [I(:),J(:)];
    weights = vals(1:K);
end

function v=getf(s,f,def), if isfield(s,f)&&~isempty(s.(f)), v=s.(f); else, v=def; end, end
