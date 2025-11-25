function result = Main_realistic_em_penalty_test_multi(output_source)
% === Realistic EM Penalty Test (Multiple Frequencies)
%     + My Adapter (BCV-E + PGD)
%     + Per-frequency & Cross-frequency Metrics/Plots
%     + PCoh Panels per Frequency (like Sim2 multiple)
%
% 依赖：
%   higgs.m, eloreta_hg_lasso.m, lcmv_hg_lasso.m, my_test2_adapter_bcvE.m
%   compute_spd_and_classification_metrics.m（若无，则逐频指标只导出Fro误差）

if nargin < 1 || isempty(output_source), output_source = pwd; end
out_dir = fullfile(output_source,'Realistic_em_penalty_test_multi');
if ~isfolder(out_dir), mkdir(out_dir); end

%% -------------------------
%% A) 载入几何/Leadfield/Seeders（与单频版一致）
%% -------------------------
sens_system = 'small';    % 'large' or 'small'
try
    if strcmp(sens_system,'large')
        load('HeadModel_large.mat');
    else
        load('HeadModel_small.mat');
    end
catch
    warning('HeadModel_*.mat not found. Generating assets ...');
    make_sim3_assets();
    if strcmp(sens_system,'large')
        load('HeadModel_large.mat');
    else
        load('HeadModel_small.mat');
    end
end
vertices = cortex.vertices; faces = cortex.faces;

try
    if strcmp(sens_system,'large')
        load('data_tips24_large.mat');
    else
        load('data_tips22_small.mat');
    end
catch
    warning('data_tips*.mat not found. Regenerating assets ...');
    make_sim3_assets();
    if strcmp(sens_system,'large')
        load('data_tips24_large.mat');
    else
        load('data_tips22_small.mat');
    end
end

try
    if strcmp(sens_system,'large')
        load('LeadFields_large.mat'); L_full = LeadFields{1,1};
    else
        load('LeadFields_small.mat'); L_full = LeadFields{1,1};
    end
catch
    warning('LeadFields_*.mat not found. Regenerating assets ...');
    make_sim3_assets();
    if strcmp(sens_system,'large')
        load('LeadFields_large.mat'); L_full = LeadFields{1,1};
    else
        load('LeadFields_small.mat'); L_full = LeadFields{1,1};
    end
end

% 由 tips 选 Seeder（与单频版一致）
Nseed  = numel(data_tips);
Seeders = zeros(Nseed,1);
for k = 1:Nseed
    coord = data_tips(k).Position;
    Seeders(k) = pickpoint(coord(1), coord(2), coord(3), vertices, 1E-3);
end
Seeders = Seeders(randperm(Nseed)); %#ok<RANDP>
L = L_full(:, Seeders);

% 邻域索引（生物噪声用）
d0  = 5E-3;
index_full = [];
for point = 1:Nseed
    Source = Seeders(point);
    [index, ~] = surfpatch(Source, vertices, faces, d0);
    index_full = [index_full; index]; %#ok<AGROW>
end
index_full = unique(index_full(:));
p = size(L,1); q = numel(Seeders);

%% -------------------------
%% B) 多频真值构造 + 采样（等距同构）
%% -------------------------
% 频点设置
F  = 10;                     % 频点数（可改 8~12）
m  = strcmp(sens_system,'large') * 6000 + strcmp(sens_system,'small') * 600;
m  = max(m, 600);            % 样本数（每频）

% 先用 gen_hggm2 得到“结构基准”
options = struct();
options.config      = 2;  % overlapping
options.var         = 2;  % complex
options.extensions  = [ceil(Nseed/3); ceil(Nseed/3); Nseed - 2*ceil(Nseed/3)];
options.connections = [1 2; 2 3];
[~, ~, Theta_base] = gen_hggm2(m, q, options);     % 仅取结构参考

% 多频真值（平滑幅度 + 轻微扰动）
cfg_true = struct();
cfg_true.sigma_time         = 2.5;       % 频轴平滑核宽（越大越平滑）
cfg_true.amp_range          = [0.7,1.3]; % 幅值相对比例范围
cfg_true.phase_jitter_rad   = pi/18;     % 相位轻微扰动（弧度）
cfg_true.support_flip_ratio = 0.05;      % 少量边在少数频点开/关
cfg_true.eps_eig_floor      = 1e-6;      % SPD 地板
cfg_true.rng_seed           = 42;

[Theta_true, Sigma_true] = make_multifreq_true_thetas(Theta_base, F, cfg_true);

% 等距同构采样（每频 j_sim_f ~ CN(0,Sigma_true{f})）
J_cell  = cell(F,1);
V_cell  = cell(F,1);
Svv3D   = zeros(p,p,F);
for f = 1:F
    jf   = sample_complex_gaussian_isomorph(Sigma_true{f}, m);   % q×m
    v0   = L * jf;                                              % p×m

    % 生物噪声（Seeder邻域）与传感器噪声（每频独立）
    bio  = randn(numel(index_full),m) + 1i*randn(numel(index_full),m);
    bio  = L_full(:, index_full) * bio;
    bio  = norm(v0,'fro') * bio / max(norm(bio,'fro'), eps);

    sns  = randn(p,m) + 1i*randn(p,m);
    sns  = norm(v0,'fro') * sns / max(norm(sns,'fro'), eps);

    v    = v0 + 0.1*bio + 0.1*sns;
    Svvf = (v*v')/m; Svvf = 0.5*(Svvf + Svvf');
    Svv3D(:,:,f) = Svvf;

    J_cell{f} = jf; V_cell{f} = v; %#ok<NASGU>
end

emp_cov_cell = squeeze(num2cell(Svv3D, [1 2])).';  % 1×F cell，每格 p×p

%% -------------------------
%% C) 估计（基线逐频 + 我的方法跨频）
%% -------------------------
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
param.m             = m;
param.nu            = m;
param.eigreg        = 1E-4;
aj                  = sqrt(log(q)/m);
Ajj_diag            = 0; 
Ajj_ndiag           = 1;
Ajj                 = Ajj_diag*eye(q)+Ajj_ndiag*(ones(q)-eye(q));
param.aj            = aj;
param.Ajj           = Ajj;
param.axi           = 1E-4;
param.Axixi         = eye(p);
param.Axixi_inv     = eye(p);
param.ntry          = 0;
param.prew          = 0;
param.rth1          = 0.7;
param.rth2          = 3.16;
param.method        = 'lqa';

methods = {'higgs_lasso','higgs_ridge','higgs_naive','eloreta_hglasso','lcmv_hglasso','my_adapter'};
M = numel(methods);

Theta_est = cell(M,1); Sigma_est = cell(M,1);
for im=1:M, Theta_est{im} = cell(F,1); Sigma_est{im} = cell(F,1); end

% 1..3 H-HGGM by freq
for f=1:F
    % LASSO
    param.penalty = 1;
    [T1,S1] = higgs(Svv3D(:,:,f), L, param);
    Theta_est{1}{f} = 0.5*(T1+T1'); Sigma_est{1}{f} = 0.5*(S1+S1');
    % RIDGE
    param.penalty = 2;
    [T2,S2] = higgs(Svv3D(:,:,f), L, param);
    Theta_est{2}{f} = 0.5*(T2+T2'); Sigma_est{2}{f} = 0.5*(S2+S2');
    % NAIVE
    param.penalty = 0;
    [T3,S3] = higgs(Svv3D(:,:,f), L, param);
    Theta_est{3}{f} = 0.5*(T3+T3'); Sigma_est{3}{f} = 0.5*(S3+S3');
end

% 4 eLORETA + hg-lasso (by freq)
param.gamma1      = 0.001;
param.gamma2      = 0.05;
param.delta_gamma = 0.001;
for f=1:F
    [T4,S4] = eloreta_hg_lasso(Svv3D(:,:,f), L, param);
    Theta_est{4}{f} = 0.5*(T4+T4'); Sigma_est{4}{f} = 0.5*(S4+S4');
end

% 5 LCMV + hg-lasso (by freq)
for f=1:F
    param.gamma    = sum(abs(diag(Svv3D(:,:,f))))/(p*100);
    [T5,S5]        = lcmv_hg_lasso(Svv3D(:,:,f), L, param);
    Theta_est{5}{f} = 0.5*(T5+T5'); Sigma_est{5}{f} = 0.5*(S5+S5');
end

% 6 我的方法：一次跨频
cfg_my = struct( ...
  'output_format','cell','verbose',true, ...
  'do_support_refit',true, ...
  'warm_method','ssblpp', ...
  'kernel_sigma',3.0, 'lambda3_ratio',0.0, ...
  'use_subspace',true, 'do_scale_L',true, ...
  'noise_model','scalar', ...
  'sigma2xi', median(arrayfun(@(ff) trace(Svv3D(:,:,ff))/p, 1:F)) * 1e-2 ...
);
[Om_cell, Sig_cell] = my_test2_adapter_bcvE(emp_cov_cell, L, m, cfg_my);
for f=1:F
    Theta_est{6}{f} = 0.5*(Om_cell{f}+Om_cell{f}'); 
    Sigma_est{6}{f} = 0.5*(Sig_cell{f}+Sig_cell{f}');
end

%% -------------------------
%% D) 去偏 + Rayleigh，逐频指标
%% -------------------------
rth = 3.16;
Theta_unb = cell(M,1); Theta_ray = cell(M,1);
for im=1:M
    Theta_unb{im} = cell(F,1); Theta_ray{im} = cell(F,1);
    for f=1:F
        T = 0.5*(Theta_est{im}{f}+Theta_est{im}{f}');
        S = 0.5*(Sigma_est{im}{f}+Sigma_est{im}{f}');
        Tunb = 2*T - T*S*T;
        dv   = abs(diag(T));
        Tvar = sqrt(dv*dv.' + abs(T).^2);
        Tray = Tunb;
        mask = abs(Tray) < (rth/sqrt(m))*(Tvar - diag(diag(Tvar)));
        Tray(mask) = 0; Tray(1:q+1:end) = 0;
        Theta_unb{im}{f} = 0.5*(Tunb+Tunb');
        Theta_ray{im}{f} = 0.5*(Tray+Tray');
    end
end

% 逐频指标表
per_rows = {};
have_cls = (exist('compute_spd_and_classification_metrics.m','file')==2);
for im=1:M
    for f=1:F
        % 基本 Fro 相对误差
        fr_theta = norm(Theta_est{im}{f}-Theta_true{f},'fro')/max(norm(Theta_true{f},'fro'),eps);
        fr_sigma = norm(Sigma_est{im}{f}-Sigma_true{f},'fro')/max(norm(Sigma_true{f},'fro'),eps);
        AUC=NaN; F1=NaN; SEN=NaN; SPE=NaN; PREC=NaN;
        SPD_LE=NaN; SPD_AIRM=NaN; SPD_JBLD=NaN; SPD_Bures=NaN;
        if have_cls
            args = struct('Theta_true',Theta_true{f}, ...
                          'Theta_est', Theta_est{im}{f}, ...
                          'Theta_unb', Theta_unb{im}{f}, ...
                          'Theta_ray', Theta_ray{im}{f}, ...
                          'is_complex', true);
            opts = struct('normalize_mode','maxabs', ...
                          'spd', struct('symmetrize',true,'project',true,'eps',1e-10), ...
                          'alpha', struct('value',0,'use_bcv_bug',false), ...
                          'plot', struct('radar',false,'title',''));
            out = compute_spd_and_classification_metrics(args, opts);
            try
                AUC  = out.metrics.est.auc;  F1   = out.metrics.est.f1;
                SEN  = out.metrics.est.sens; SPE  = out.metrics.est.spec;
                PREC = out.metrics.est.prec;
            catch, end
            try
                SPD_LE   = out.spd.est.logeuclid;
                SPD_AIRM = out.spd.est.airm;
                SPD_JBLD = out.spd.est.jbld;
                SPD_Bures= out.spd.est.bures;
            catch, end
        end
        per_rows(end+1,1) = {table( string(methods{im}), f, fr_theta, fr_sigma, ...
            AUC, F1, SEN, SPE, PREC, SPD_LE, SPD_AIRM, SPD_JBLD, SPD_Bures, ...
            'VariableNames',{'method','freq','FrobThetaRel','FrobSigmaRel', ...
                             'AUC','F1','SEN','SPE','PREC','SPD_LE','SPD_AIRM','SPD_JBLD','SPD_Bures'})}; %#ok<AGROW>
    end
end
per_tbl = vertcat(per_rows{:});
writetable(per_tbl, fullfile(out_dir, sprintf('per_freq_metrics_(%s).csv', sens_system)));

%% -------------------------
%% E) 跨频指标（TV_ratio / 邻频Jaccard / 边轨迹相关）
%% -------------------------
% 真值TV
TV_true = 0;
for f=2:F
    TV_true = TV_true + norm(Theta_true{f}-Theta_true{f-1},'fro');
end
TV_true = max(TV_true, eps);

% 真值Top-K边（频均绝对值排序）
K_for_tracks = min(16, q);
edge_list = topK_edges_from_true(Theta_true, K_for_tracks);

sum_rows = {};
for im=1:M
    % TV_est
    TV_est = 0;
    for f=2:F
        TV_est = TV_est + norm(Theta_est{im}{f}-Theta_est{im}{f-1},'fro');
    end
    TV_ratio = TV_est / TV_true;

    % 邻频 Jaccard（Rayleigh 支撑）
    Jv = nan(F-1,1);
    for f=2:F
        S1 = abs(Theta_ray{im}{f})   > 0; S1(1:q+1:end)=false;
        S0 = abs(Theta_ray{im}{f-1}) > 0; S0(1:q+1:end)=false;
        inter = nnz(S1 & S0); uni = nnz(S1 | S0);
        if uni>0, Jv(f-1)=inter/uni; end
    end
    Jaccard_mean = mean(Jv,'omitnan');

    % Top-K 轨迹相关
    C = nan(size(edge_list,1),1);
    for k=1:size(edge_list,1)
        i=edge_list(k,1); j=edge_list(k,2);
        t_true = zeros(F,1); t_est = zeros(F,1);
        for f=1:F
            t_true(f) = abs(Theta_true{f}(i,j));
            t_est(f)  = abs(Theta_est{im}{f}(i,j));
        end
        if std(t_true)>0 && std(t_est)>0
            C(k) = corr(t_true, t_est, 'type','Pearson','rows','complete');
        end
    end
    EdgeTrackCorr_mean = mean(C,'omitnan');

    % 若逐频评估里带了 AUC/F1（est 态），这里求均值
    AUC_mean = NaN; F1_mean = NaN;
    try
        A = per_tbl.AUC( per_tbl.method==string(methods{im}) );
        Ff= per_tbl.F1(  per_tbl.method==string(methods{im}) );
        AUC_mean = mean(A,'omitnan'); F1_mean = mean(Ff,'omitnan');
    catch, end

    sum_rows{end+1,1} = table( string(methods{im}), TV_ratio, Jaccard_mean, EdgeTrackCorr_mean, ...
        AUC_mean, F1_mean, 'VariableNames', ...
        {'method','TV_ratio','Jaccard_mean','EdgeTrackCorr_mean','AUC_mean','F1_mean'}); %#ok<AGROW>
end
summary_tbl = vertcat(sum_rows{:});
writetable(summary_tbl, fullfile(out_dir, sprintf('summary_metrics_(%s).csv', sens_system)));

%% -------------------------
%% F) 可视化（PCoh panels + Top-K 轨迹 + 逐频误差曲线）
%% -------------------------
% 1) PCoh panels（与 Sim2-multiple 一致风格）
save_pcoh_panels(Theta_true, Theta_est, Theta_ray, methods, out_dir, sens_system);

% 2) Top-K 边轨迹
fig1 = figure('Name','TopK edge tracks','Position',[120 120 1000 560]);
Kshow = min(8, size(edge_list,1));  % 展示前K条
cols = lines(M);  % 方法配色
for kk=1:Kshow
    subplot(2,ceil(Kshow/2),kk); hold on;
    i=edge_list(kk,1); j=edge_list(kk,2);
    t_true = zeros(F,1);
    for f=1:F, t_true(f) = abs(Theta_true{f}(i,j)); end
    plot(1:F, t_true,'k-','LineWidth',2,'DisplayName','true');
    for im=1:M
        t_est = zeros(F,1);
        for f=1:F, t_est(f) = abs(Theta_est{im}{f}(i,j)); end
        plot(1:F, t_est, '-','LineWidth',1.5,'Color',cols(im,:), 'DisplayName',methods{im});
    end
    title(sprintf('edge (%d,%d)',i,j)); xlabel('freq idx'); ylabel('|Theta_{ij}|');
    if kk==1, legend('Location','best'); end
    hold off;
end
saveas(fig1, fullfile(out_dir, sprintf('edge_tracks_topK_(%s).fig', sens_system)));
close(fig1);

% 3) 逐频 Fro 误差曲线（Theta）
fig2 = figure('Name','Per-frequency Fro errors','Position',[100 100 900 520]); hold on;
for im=1:M
    y = zeros(F,1);
    for f=1:F
        y(f) = norm(Theta_est{im}{f}-Theta_true{f},'fro')/max(norm(Theta_true{f},'fro'),eps);
    end
    plot(1:F, y, '-o','LineWidth',1.8,'DisplayName',methods{im});
end
xlabel('freq idx'); ylabel('Rel Fro error (Theta)'); grid on; legend('Location','best');
saveas(fig2, fullfile(out_dir, sprintf('perfreq_fro_theta_(%s).fig', sens_system)));
close(fig2);

%% -------------------------
%% G) 汇总输出
%% -------------------------
result = struct();
result.sens_system = sens_system;
result.F           = F;
result.Svv3D       = Svv3D;
result.L           = L;
result.Seeders     = Seeders;
result.Theta_true  = Theta_true;
result.Sigma_true  = Sigma_true;
result.Theta_est   = Theta_est;
result.Sigma_est   = Sigma_est;
result.Theta_unb   = Theta_unb;
result.Theta_ray   = Theta_ray;
result.paths = struct( ...
    'per_freq_metrics', fullfile(out_dir, sprintf('per_freq_metrics_(%s).csv', sens_system)), ...
    'summary_metrics',  fullfile(out_dir, sprintf('summary_metrics_(%s).csv', sens_system)), ...
    'pcoh_panels_fig',  fullfile(out_dir, sprintf('partial_coherence_maps_(%s).fig', sens_system)) ...
);
save(fullfile(out_dir, sprintf('result_multi_(%s).mat', sens_system)),'-struct','result');

disp(['[OK] Results saved under: ', out_dir]);

end  % ===== main =====


%% ===== helpers =====

function save_pcoh_panels(Theta_true, Theta_est, Theta_ray, methods, out_dir, sens_system)
% 生成与 Sim2-multiple 一致的 PCoh 面板：行=真值+各方法(raw)+my_adapter(ray)；列=频点
F = numel(Theta_true);
M = numel(methods);
nrows = 1 + M + 1;  % true + M methods (raw) + my-adapter(ray)

fig = figure('Name','PCoh panels','Position',[100,80, 160*F, 130*nrows]);
try, load('colormap3'); colormap(cmap); catch, colormap('hot'); end

% 一行的绘制器
    function draw_row(ridx, Xcell, ttlprefix)
        for f=1:F
            X = Xcell{f};
            X = X - diag(diag(X));
            mx = max(abs(X(:))+eps);
            X = X / mx;
            subplot(nrows, F, (ridx-1)*F + f);
            imagesc(abs(X));
            axis square tight;
            if ridx==nrows, xlabel(sprintf('f=%d',f)); end
            if f==1, ylabel('generators'); end
            if f==ceil(F/2)
                title(sprintf('%s', ttlprefix),'FontWeight','normal');
            end
        end
    end

% 1) true
draw_row(1, Theta_true, 'true');

% 2..(1+M) methods raw
for im=1:M
    draw_row(1+im, Theta_est{im}, methods{im});
end

% last) my-adapter (ray)
kmy = find(strcmp(methods,'my_adapter'),1);
if ~isempty(kmy)
    draw_row(nrows, Theta_ray{kmy}, 'my_adapter (ray)');
else
    % 占位（若未找到 my_adapter）
    draw_row(nrows, Theta_est{end}, [methods{end} ' (raw)']);
end

saveas(fig, fullfile(out_dir, sprintf('partial_coherence_maps_(%s).fig', sens_system)));
close(fig);
end

function [Theta_true, Sigma_true] = make_multifreq_true_thetas(Theta_base, F, cfg)
% 基于 Theta_base 的支撑/相位，构造跨频平滑的 Θ*_f
q = size(Theta_base,1);
Theta_true = cell(F,1); Sigma_true = cell(F,1);

rng(getf(cfg,'rng_seed',42));
sigma_time = getf(cfg,'sigma_time',2.5);
arng       = getf(cfg,'amp_range',[0.7,1.3]);
phi_jit    = getf(cfg,'phase_jitter_rad', pi/18);
flip_ratio = getf(cfg,'support_flip_ratio', 0.05);
eps_floor  = getf(cfg,'eps_eig_floor',1e-6);

M = (abs(Theta_base)>0); M(1:q+1:end)=false;       % off-diag 支撑
Phi0 = angle(Theta_base);

% 生成一个 1D 高斯核（时间平滑）
tt = (-3*ceil(sigma_time)):(3*ceil(sigma_time));
gk = exp(-(tt.^2)/(2*sigma_time^2)); gk = gk / sum(gk);

% 支撑翻转位置（少量边、少数频点）
[I,J] = find(triu(M,1));
E = numel(I);
nflip_edges = round(flip_ratio * E);
flip_edges = randsample(E, max(0,nflip_edges));
flip_fmask = false(E,F);
for idx = flip_edges(:).'
    % 随机 1~2 个频点被“关掉”，并用核平滑过渡（此处简化为硬开关）
    nf = randi([1,2]);
    pos = randsample(F, nf);
    flip_fmask(idx,pos) = true;
end

for f=1:F
    Theta_f = zeros(q,q);
    for e=1:E
        i=I(e); j=J(e);
        mag0 = abs(Theta_base(i,j));
        % 平滑随机幅度（映射到 [a_min,a_max]）
        z = randn(1,F); z = conv(z, gk, 'same');
        z = (z - mean(z)) / (std(z)+eps);
        a = arng(1) + (arng(2)-arng(1)) * ( (z - min(z)) / (max(z)-min(z) + eps) );
        aij = a(f) * mag0;

        % 支撑翻转：若该边在 f 被“关掉”，则 aij=0
        if flip_fmask(e,f), aij = 0; end

        % 相位（轻微抖动）
        phi = Phi0(i,j) + phi_jit * randn();

        Theta_f(i,j) = aij * exp(1i*phi);
        Theta_f(j,i) = conj(Theta_f(i,j));
    end
    % 对角项
    dbase = real(diag(Theta_base));
    if all(dbase>0)
        Theta_f(1:q+1:end) = dbase;
    else
        Theta_f(1:q+1:end) = max(1.0, mean(abs(Theta_base(:))));
    end

    % Hermitian + SPD 投影
    Theta_f = 0.5*(Theta_f + Theta_f');
    Theta_f = proj_spd(Theta_f, eps_floor);

    Theta_true{f} = Theta_f;
    Sigma_true{f} = inv_psd_robust_(Theta_f, 1e-8, 1e-12);  % 稳健逆
end
end

function X = proj_spd(X, eps_floor)
X = 0.5*(X+X'); 
[U,D] = eig(full(X),'vector'); d = real(D);
dmax = max(d);
floor_val = max(eps_floor*max(dmax,1), 1e-10);
d(d<floor_val) = floor_val;
X = U*diag(d)*U'; X = 0.5*(X+X');
end

function Data = sample_complex_gaussian_isomorph(Sigma, m)
% 复高斯等距同构采样：R^{2q} ~ N(0, Wisomph) → C^q
q = size(Sigma,1);
Sigma = 0.5*(Sigma + Sigma');   % Hermitian
W = [real(Sigma) -imag(Sigma); imag(Sigma) real(Sigma)];
% Chol 分解，若数值不稳则做特征值下界
[Q,pflag] = chol(W,'lower');
if pflag~=0
    [U,D] = eig((W+W')/2,'vector'); d = real(D);
    dmax = max(d); d = max(d, 1e-12*max(dmax,1));
    Q = U*diag(sqrt(d));
end
Z = randn(2*q, m);
Y = Q * Z;    % 2q×m
DataRe = Y(1:q, :);
DataIm = Y(q+1:2*q, :);
Data   = DataRe + 1i*DataIm;     % q×m
end

function A = inv_psd_robust_(S, eps_reg, min_ratio)
S = 0.5*(S+S'); [U,D]=eig(full(S),'vector'); d=real(D);
dmax = max(d); floor_val = max(min_ratio*max(dmax,eps), 0);
d(d<floor_val) = floor_val;
if eps_reg>0, d = (d + eps_reg*dmax)/(1+eps_reg); end
A = U*diag(1./d)*U'; A = 0.5*(A+A');
end

function E = topK_edges_from_true(Theta_true, K)
% 频均绝对值排序取前K条边（上三角）
q = size(Theta_true{1},1);
F = numel(Theta_true);
Amean = zeros(q,q);
for f=1:F
    Amean = Amean + abs(Theta_true{f});
end
Amean = Amean / F; Amean(1:q+1:end) = 0;
[IU,JU,V] = find(triu(Amean,1));
[~,ord] = sort(V,'descend');
IU = IU(ord); JU = JU(ord);
K = min(K, numel(IU));
E = [IU(1:K), JU(1:K)];
end

function v=getf(s,f,def), if isfield(s,f)&&~isempty(s.(f)), v=s.(f); else, v=def; end, end
