function result = Main_realistic_em_penalty_test_multi_large(output_source)
% === Realistic (High-dim 19×360) EM Penalty Test (Multiple Frequencies)
%     + HCP_MMP1 ROI leadfield
%     + Baselines (per freq) & My Adapter (cross-freq)
%     + Unified colormap/color-scale; one window per method
%     + Per-frequency & Cross-frequency metrics
%     + Per-method wall-clock time logging
%
% 依赖：
%   - higgs.m, eloreta_hg_lasso.m, lcmv_hg_lasso.m
%   - my_test2_adapter_bcvE.m
%   - compute_spd_and_classification_metrics.m（若无，只导出Fro误差）
%   - 本文件内 helpers: build_roi_leadfield_hcp, make_multifreq_true_thetas, etc.

if nargin < 1 || isempty(output_source), output_source = pwd; end
out_dir = fullfile(output_source,'Realistic_em_penalty_test_multi_large');
if ~isfolder(out_dir), mkdir(out_dir); end
sens_tag = 'large';   % 文件名标记

%% -------------------------
%% A) 载入真实头模 & 构造 19×360 ROI leadfield
%% -------------------------
try
    a = load('F:\test_BCV\BC-V_Structure\ICBM152_2023\surf\surf.mat');
    b = load('F:\test_BCV\BC-V_Structure\ICBM152_2023\leadfield\headmodel.mat');
catch
    error('未找到 surf.mat，请确认路径。');
end

% 19 通道（如需 18 通道，可传 sens_idx 指定保留的行）
[L, R, roi_info, sens_idx] = build_roi_leadfield_hcp(a, b,...
    'atlas_idx', [], 'use_area', true, 'sign_align', true, 'sens_idx', 1:19); %#ok<ASGLU>

V  = a.Sc.Vertices; %#ok<NASGU>
Ff = a.Sc.Faces;    %#ok<NASGU>

p = size(L,1);      % 19
q = size(L,2);      % 360

%% -------------------------
%% B) 多频真值 & 采样（等距同构）
%% -------------------------
F = 10;                        % 频点数
m = 600;                       % 每频样本数（可调 600~1200）
fprintf('[config] p=%d, q=%d, F=%d, m=%d\n', p, q, F, m);

% 先用简化块结构得到一个 base（仅作支撑与相位参考）
options = struct('config',2,'var',2, ...
                 'extensions',[ceil(q/3); ceil(q/3); q-2*ceil(q/3)], ...
                 'connections',[1 2; 2 3]);
[~,~,Theta_base] = gen_hggm2(m, q, options);

% 生成跨频真值
cfg_true = struct('sigma_time',2.5, 'amp_range',[0.7,1.3], ...
                  'phase_jitter_rad',pi/18, 'support_flip_ratio',0.05, ...
                  'eps_eig_floor',1e-6, 'rng_seed',42);
[Theta_true, Sigma_true] = make_multifreq_true_thetas(Theta_base, F, cfg_true);

% 采样 + 生成观测
Svv3D = zeros(p,p,F);
for f = 1:F
    jf   = sample_complex_gaussian_isomorph(Sigma_true{f}, m);   % q×m
    v0   = L * jf;                                              % p×m
    % 仅传感器噪声（high-dim 下已足够）
    sns  = randn(p,m) + 1i*randn(p,m);
    sns  = norm(v0,'fro') * sns / max(norm(sns,'fro'), eps);
    v    = v0 + 0.1*sns;
    Svvf = (v*v')/m; Svvf = 0.5*(Svvf+Svvf');
    Svv3D(:,:,f) = Svvf;
end
emp_cov_cell = squeeze(num2cell(Svv3D,[1 2])).';

%% -------------------------
%% C) 估计（基线逐频 + 我的跨频），并记录耗时
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

time_log = struct('method', methods, 'seconds', num2cell(zeros(1, M)));

% 1..3 H-HGGM per freq + timing
tic;
for f=1:F
    param.penalty = 1; [T1,S1] = higgs(Svv3D(:,:,f), L, param);
    Theta_est{1}{f}=0.5*(T1+T1'); Sigma_est{1}{f}=0.5*(S1+S1');
end
time_log(1).seconds = toc;

tic;
for f=1:F
    param.penalty = 2; [T2,S2] = higgs(Svv3D(:,:,f), L, param);
    Theta_est{2}{f}=0.5*(T2+T2'); Sigma_est{2}{f}=0.5*(S2+S2');
end
time_log(2).seconds = toc;

tic;
for f=1:F
    param.penalty = 0; [T3,S3] = higgs(Svv3D(:,:,f), L, param);
    Theta_est{3}{f}=0.5*(T3+T3'); Sigma_est{3}{f}=0.5*(S3+S3');
end
time_log(3).seconds = toc;

% 4 eLORETA + hg-lasso
param.gamma1      = 0.001;
param.gamma2      = 0.05;
param.delta_gamma = 0.001;
tic;
for f=1:F
    [T4,S4] = eloreta_hg_lasso(Svv3D(:,:,f), L, param);
    Theta_est{4}{f}=0.5*(T4+T4'); Sigma_est{4}{f}=0.5*(S4+S4');
end
time_log(4).seconds = toc;

% 5 LCMV + hg-lasso
tic;
for f=1:F
    param.gamma    = sum(abs(diag(Svv3D(:,:,f))))/(p*100);
    [T5,S5] = lcmv_hg_lasso(Svv3D(:,:,f), L, param);
    Theta_est{5}{f}=0.5*(T5+T5'); Sigma_est{5}{f}=0.5*(S5+S5');
end
time_log(5).seconds = toc;

% 6 我的跨频
cfg_my = struct( ...
  'output_format','cell','verbose',true, ...
  'do_support_refit',true, ...
  'warm_method','ssblpp', ...
  'kernel_sigma',3.0, 'lambda3_ratio',0.0, ...
  'use_subspace',true, 'do_scale_L',true, ...
  'noise_model','scalar', ...
  'sigma2xi', median(arrayfun(@(ff) trace(Svv3D(:,:,ff))/p, 1:F)) * 1e-2 ...
);
tic;
[Om_cell, Sig_cell] = my_test2_adapter_bcvE(emp_cov_cell, L, m, cfg_my);
time_log(6).seconds = toc;
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

%% -------------------------
%% E) 逐频 & 跨频数值指标（写CSV）
%% -------------------------
per_rows = {};
have_cls = (exist('compute_spd_and_classification_metrics.m','file')==2);
for im=1:M
    for f=1:F
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
                          'plot', struct('radar',false));
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
        per_rows{end+1,1} = table( string(methods{im}), f, fr_theta, fr_sigma, ...
            AUC, F1, SEN, SPE, PREC, SPD_LE, SPD_AIRM, SPD_JBLD, SPD_Bures, ...
            'VariableNames',{'method','freq','FrobThetaRel','FrobSigmaRel', ...
                             'AUC','F1','SEN','SPE','PREC','SPD_LE','SPD_AIRM','SPD_JBLD','SPD_Bures'}); %#ok<AGROW>
    end
end
per_tbl = vertcat(per_rows{:});
writetable(per_tbl, fullfile(out_dir, sprintf('per_freq_metrics_(%s).csv', sens_tag)));

% 跨频 TV_ratio / 邻频 Jaccard / Top-K 轨迹相关
TV_true = 0; for f=2:F, TV_true = TV_true + norm(Theta_true{f}-Theta_true{f-1},'fro'); end
TV_true = max(TV_true, eps);

edge_list = topK_edges_from_true(Theta_true, min(16, q));
sum_rows = {};
for im=1:M
    TV_est = 0; for f=2:F, TV_est = TV_est + norm(Theta_est{im}{f}-Theta_est{im}{f-1},'fro'); end
    TV_ratio = TV_est / TV_true;

    Jv = nan(F-1,1);
    for f=2:F
        S1 = abs(Theta_ray{im}{f})   > 0; S1(1:q+1:end)=false;
        S0 = abs(Theta_ray{im}{f-1}) > 0; S0(1:q+1:end)=false;
        inter = nnz(S1 & S0); uni = nnz(S1 | S0);
        if uni>0, Jv(f-1)=inter/uni; end
    end
    Jaccard_mean = mean(Jv,'omitnan');

    C = nan(size(edge_list,1),1);
    for k=1:size(edge_list,1)
        i=edge_list(k,1); j=edge_list(k,2);
        t_true = zeros(F,1); t_est = zeros(F,1);
        for f=1:F, t_true(f) = abs(Theta_true{f}(i,j)); t_est(f) = abs(Theta_est{im}{f}(i,j)); end
        if std(t_true)>0 && std(t_est)>0, C(k)=corr(t_true,t_est,'type','Pearson','rows','complete'); end
    end
    EdgeTrackCorr_mean = mean(C,'omitnan');

    AUC_mean = NaN; F1_mean = NaN;
    try
        A = per_tbl.AUC( per_tbl.method==string(methods{im}) );
        Ff= per_tbl.F1(  per_tbl.method==string(methods{im}) );
        AUC_mean = mean(A,'omitnan'); F1_mean = mean(Ff,'omitnan');
    catch, end

    sum_rows{end+1,1} = table(string(methods{im}), TV_ratio, Jaccard_mean, EdgeTrackCorr_mean, ...
        AUC_mean, F1_mean, 'VariableNames', ...
        {'method','TV_ratio','Jaccard_mean','EdgeTrackCorr_mean','AUC_mean','F1_mean'}); %#ok<AGROW>
end
summary_tbl = vertcat(sum_rows{:});
writetable(summary_tbl, fullfile(out_dir, sprintf('summary_metrics_(%s).csv', sens_tag)));

% 写出时间日志
times_tbl = struct2table(time_log);
writetable(times_tbl, fullfile(out_dir, sprintf('method_times_(%s).csv', sens_tag)));

%% -------------------------
%% F) 可视化：每个方法一个窗口 + 统一色标
%% -------------------------
% 先统计全体要画的 |Theta| 的“全局 vmax”（去对角）
global_vmax = 0;
% true
for f=1:F, X=Theta_true{f}; X(1:q+1:end)=0; global_vmax=max(global_vmax, max(abs(X(:)))); end
% methods (raw)
for im=1:M
    for f=1:F, X=Theta_est{im}{f}; X(1:q+1:end)=0; global_vmax=max(global_vmax, max(abs(X(:)))); end
end
% my-adapter (ray) 也画一个
kmy = find(strcmp(methods,'my_adapter'),1);
if ~isempty(kmy)
    for f=1:F, X=Theta_ray{kmy}{f}; X(1:q+1:end)=0; global_vmax=max(global_vmax, max(abs(X(:)))); end
end
if ~isfinite(global_vmax) || global_vmax<=0, global_vmax = 1; end

% 统一 colormap（若有 colormap3 就用）
try, load('colormap3'); cmap_used = cmap; catch, cmap_used = hot; end

% 真值
save_pcoh_panels_onewindow(Theta_true, 'True PCoh', cmap_used, [0 global_vmax], ...
    fullfile(out_dir, sprintf('pcoh_panels_true_(%s).fig', sens_tag)));

% 各方法（raw）
name_for_title = {'HIGGS-L1','HIGGS-Ridge','HIGGS-Naive','eLORETA+HG-LASSO','LCMV+HG-LASSO','My-Adapter (raw)'};
for im=1:M
    save_pcoh_panels_onewindow(Theta_est{im}, [name_for_title{im} ' PCoh'], ...
        cmap_used, [0 global_vmax], ...
        fullfile(out_dir, sprintf('pcoh_panels_%s_(%s).fig', methods{im}, sens_tag)));
end

% my-adapter (ray)
if ~isempty(kmy)
    save_pcoh_panels_onewindow(Theta_ray{kmy}, 'My-Adapter (ray) PCoh', ...
        cmap_used, [0 global_vmax], ...
        fullfile(out_dir, sprintf('pcoh_panels_my_adapter_ray_(%s).fig', sens_tag)));
end

%% -------------------------
%% G) 汇总输出（.mat）
%% -------------------------
result = struct();
result.sens_tag     = sens_tag;
result.F            = F;
result.Svv3D        = Svv3D;
result.L            = L;
result.Theta_true   = Theta_true;
result.Sigma_true   = Sigma_true;
result.Theta_est    = Theta_est;
result.Sigma_est    = Sigma_est;
result.Theta_unb    = Theta_unb;
result.Theta_ray    = Theta_ray;
result.times        = time_log;
result.paths = struct( ...
    'per_freq_metrics', fullfile(out_dir, sprintf('per_freq_metrics_(%s).csv', sens_tag)), ...
    'summary_metrics',  fullfile(out_dir, sprintf('summary_metrics_(%s).csv', sens_tag)), ...
    'method_times',     fullfile(out_dir, sprintf('method_times_(%s).csv', sens_tag)) ...
);
save(fullfile(out_dir, sprintf('result_multi_(%s).mat', sens_tag)), '-struct', 'result');

disp(['[OK] Results saved under: ', out_dir]);

end  % ===== main =====


%% ====== helpers ======

function save_pcoh_panels_onewindow(Theta_cell, title_str, cmap_used, clim, save_path)
% 将某一组（真值或某方法）的 F 个频率点 PCoh 面板画在一个窗口
F = numel(Theta_cell);
q = size(Theta_cell{1},1);
fig = figure('Name',title_str,'Position',[100,80, 160*F, 160]);
colormap(cmap_used);
for f=1:F
    X = 0.5*(Theta_cell{f}+Theta_cell{f}');
    X(1:q+1:end)=0;
    subplot(1,F,f);
    imagesc(abs(X));
    axis square tight;
    % caxis(clim);
    if f==ceil(F/2), title(title_str,'FontWeight','normal'); end
    xlabel(sprintf('f=%d',f)); ylabel('generators');
end
colorbar('Position',[0.93 0.11 0.02 0.815]); % 统一色标
saveas(fig, save_path);
close(fig);
end

function [L_roi, R, roi_info, sens_idx] = build_roi_leadfield_hcp(a, b, varargin)
% 19×360 的 ROI leadfield（HCP_MMP1）
opts = struct('atlas_idx',[], 'use_area',true, 'sign_align',true, 'sens_idx', []);
if ~isempty(varargin), opts = setOpts(opts, varargin{:}); end

V  = a.Sc.Vertices; 
F  = a.Sc.Faces;
G  = b.HeadModel.Gain;        % (nsens × 3*nV)
Nn = b.HeadModel.GridOrient;  % (nV × 3)
p0 = size(G,1); nV = size(Nn,1);

if isempty(opts.sens_idx), sens_idx = 1:p0; else, sens_idx = opts.sens_idx(:).'; end
G = G(sens_idx, :); 
p = size(G,1);

atlas_idx = opts.atlas_idx;
if isempty(atlas_idx)
    atlas_idx = find(arrayfun(@(A) numel(A.Scouts)==360, a.Sc.Atlas), 1, 'first');
    if isempty(atlas_idx), error('找不到含 360 scouts 的 Atlas（HCP_MMP1）。'); end
end
Atlas = a.Sc.Atlas(atlas_idx);
if numel(Atlas.Scouts) ~= 360, error('所选 Atlas 的 ROI 数不是 360。'); end

% 法向定向 -> 19×nV
L_norm = zeros(p, nV);
for v = 1:nV
    Gv = G(:, 3*v-2 : 3*v);
    nv = Nn(v, :).';
    L_norm(:, v) = Gv * nv;
end

% 顶点面积
if opts.use_area, Avert = vertex_area(V, F); else, Avert = ones(nV,1); end

R = spalloc(nV, 360, round(nV*1.2));
roi_info = struct('scouts',{cell(360,1)}, 'seed',zeros(360,1), 'labels',{cell(360,1)});
for r = 1:360
    idx = Atlas.Scouts(r).Vertices(:);
    if isempty(idx), error('ROI %d 为空。', r); end
    roi_info.scouts{r} = idx; 
    roi_info.labels{r} = Atlas.Scouts(r).Label;
    seed = Atlas.Scouts(r).Seed; if isempty(seed)||seed<=0, seed = idx(1); end
    roi_info.seed(r) = seed;

    w = Avert(idx); w = w/(sum(w)+eps);

    sgn = ones(numel(idx),1);
    if opts.sign_align
        l0 = L_norm(:, seed);
        for k=1:numel(idx)
            lv = L_norm(:, idx(k));
            s = real(lv' * l0);
            sgn(k) = (s>=0) - (s<0);
        end
    end
    R(idx, r) = w .* sgn;
end
for r = 1:360
    col = full(R(:,r)); s = sum(abs(col)) + eps; R(:,r) = col / s;
end
L_roi = L_norm * R;
end

function A = vertex_area(V, F)
A = zeros(size(V,1),1);
for t = 1:size(F,1)
    v = F(t,:);
    e1 = V(v(2),:) - V(v(1),:);
    e2 = V(v(3),:) - V(v(1),:);
    Atri = 0.5 * norm(cross(e1, e2));
    A(v) = A(v) + Atri/3;
end
end

function [Theta_true, Sigma_true] = make_multifreq_true_thetas(Theta_base, F, cfg)
q = size(Theta_base,1);
Theta_true = cell(F,1); Sigma_true = cell(F,1);
rng(getf(cfg,'rng_seed',42));
sigma_time = getf(cfg,'sigma_time',2.5);
arng       = getf(cfg,'amp_range',[0.7,1.3]);
phi_jit    = getf(cfg,'phase_jitter_rad', pi/18);
flip_ratio = getf(cfg,'support_flip_ratio', 0.05);
eps_floor  = getf(cfg,'eps_eig_floor',1e-6);

M = (abs(Theta_base)>0); M(1:q+1:end)=false;
Phi0 = angle(Theta_base);
[II,JJ] = find(triu(M,1)); E = numel(II);

tt = (-3*ceil(sigma_time)):(3*ceil(sigma_time));
gk = exp(-(tt.^2)/(2*sigma_time^2)); gk = gk / sum(gk);

nflip_edges = round(flip_ratio * E);
flip_edges = randsample(E, max(0,nflip_edges));
flip_fmask = false(E,F);
for idx = flip_edges(:).'
    nf = randi([1,2]); pos = randsample(F, nf);
    flip_fmask(idx,pos) = true;
end

% 为每条边生成一个平滑幅度序列
Zbank = randn(E, F); 
for e=1:E, Zbank(e,:) = conv(Zbank(e,:), gk, 'same'); end
Zbank = (Zbank - mean(Zbank,2)) ./ (std(Zbank,0,2)+eps);

for f=1:F
    Theta_f = zeros(q,q);
    for e=1:E
        i=II(e); j=JJ(e);
        mag0 = abs(Theta_base(i,j));
        z = Zbank(e,:);
        a = arng(1) + (arng(2)-arng(1)) * ( (z - min(z)) / (max(z)-min(z) + eps) );
        aij = a(f) * mag0;
        if flip_fmask(e,f), aij = 0; end
        phi = Phi0(i,j) + phi_jit * randn();
        Theta_f(i,j) = aij * exp(1i*phi);
        Theta_f(j,i) = conj(Theta_f(i,j));
    end
    dbase = real(diag(Theta_base));
    if all(dbase>0), Theta_f(1:q+1:end) = dbase;
    else, Theta_f(1:q+1:end) = max(1.0, mean(abs(Theta_base(:))));
    end
    Theta_f = 0.5*(Theta_f+Theta_f');
    Theta_f = proj_spd(Theta_f, eps_floor);
    Theta_true{f} = Theta_f;
    Sigma_true{f} = inv_psd_robust_(Theta_f, 1e-8, 1e-12);
end
end

function X = proj_spd(X, eps_floor)
X = 0.5*(X+X');
[U,D] = eig(full(X),'vector'); d = real(D);
dmax = max(d); floor_val = max(eps_floor*max(dmax,1), 1e-10);
d(d<floor_val) = floor_val;
X = U*diag(d)*U'; X = 0.5*(X+X');
end

function Data = sample_complex_gaussian_isomorph(Sigma, m)
q = size(Sigma,1); Sigma = 0.5*(Sigma+Sigma');
W = [real(Sigma) -imag(Sigma); imag(Sigma) real(Sigma)];
[Q,pflag] = chol(W,'lower');
if pflag~=0
    [U,D] = eig((W+W')/2,'vector'); d = real(D);
    dmax = max(d); d = max(d, 1e-12*max(dmax,1));
    Q = U*diag(sqrt(d));
end
Z = randn(2*q, m); Y = Q * Z;
Data = Y(1:q,:) + 1i*Y(q+1:2*q,:);
end

function A = inv_psd_robust_(S, eps_reg, min_ratio)
S = 0.5*(S+S'); [U,D]=eig(full(S),'vector'); d=real(D);
dmax = max(d); floor_val = max(min_ratio*max(dmax,eps), 0);
d(d<floor_val) = floor_val;
if eps_reg>0, d = (d + eps_reg*dmax)/(1+eps_reg); end
A = U*diag(1./d)*U'; A = 0.5*(A+A');
end

function E = topK_edges_from_true(Theta_true, K)
q = size(Theta_true{1},1);
F = numel(Theta_true);
Amean = zeros(q,q);
for f=1:F, Amean = Amean + abs(Theta_true{f}); end
Amean = Amean / F; Amean(1:q+1:end)=0;
[IU,JU,V] = find(triu(Amean,1));
[~,ord] = sort(V,'descend'); IU = IU(ord); JU = JU(ord);
K = min(K, numel(IU)); E = [IU(1:K), JU(1:K)];
end

function v=getf(s,f,def), if isfield(s,f)&&~isempty(s.(f)), v=s.(f); else, v=def; end, end
function o=setOpts(o,varargin), for k=1:2:numel(varargin), o.(varargin{k})=varargin{k+1}; end, end
