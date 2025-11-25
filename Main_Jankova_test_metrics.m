function result = Main_Jankova_test_metrics(output_source)
%% Main_Jankova_test_metrics — HG-LASSO(LQA) vs PGD-only + 距离/分类指标评估
%
% 保留原可视化（10 子图/规模），并为每种方法（hg/pgd）的三类矩阵：
%   - 估计 (est)
%   - 去偏 (unb)
%   - Rayleigh 矫正 (ray)
% 计算与真值的 6 种 SPD 距离 + ROC/AUC + 最优点(SENS,SPEC,PREC,F1)，并绘制雷达图。
%
% 依赖：
%   gen_hggm1.m, hg_lasso_lqa1.m, adapter_pgd_source_only.m,
%   module6_hyperparameter_config.m, module5_proximal.m, module8_recoloring.m,
%   compute_spd_and_classification_metrics.m
%
% --------------------------------------
% 初始化输出目录
% --------------------------------------
output_source = strcat(output_source, filesep,'Jankova_test');
if(~isfolder(output_source)), mkdir(output_source); end

sizes  = [30,100,1000, 3000];
trial  = 1; ntrial = length(sizes);

fig = figure('Name','HG-LASSO vs PGD-only (Jankova)','Position',[50 50 1800 900]);
process_waitbar = waitbar(0,'Please wait...');

% 结果结构（追加 eval_* 字段）
result = struct('q',{},'m',{},'Theta_sim',{}, ...
                'Theta_hg',{},'Theta_unb_hg',{},'Theta_var_hg',{},'Theta_ray_hg',{}, ...
                'Theta_pgd',{},'Theta_unb_pgd',{},'Theta_var_pgd',{},'Theta_ray_pgd',{}, ...
                'hist_bins_hg',{},'hist_pdf_hg',{}, ...
                'hist_bins_pgd',{},'hist_pdf_pgd',{}, ...
                'eval_hg',{},'eval_pgd',{});

for q = sizes
    waitbar(trial/ntrial,process_waitbar,sprintf('size # %d',q));

    % ---------------------------
    % 仿真设置（与原版一致）
    % ---------------------------
    m              = 100*q;                 % Sample number
    nblocks        = 3*floor(log10(q));     % Number of blocks in simulation
    options.config = 2;                     % (2) overlapping blocks (1) nonoverlapping blocks
    options.var    = 2;                     % (2) complex variable (1) real variable
    maxiter        = 30;                    % Maximum number of iterations
    a              = sqrt(log(q)/m);        % Global Regularization parameter
    a_diag         = 0;                     % Selection Mask diagonal
    a_ndiag        = 1;                     % Selection Mask nondiagonal
    A              = a_diag*eye(q)+a_ndiag*(ones(q)-eye(q));  % Selection Mask squared
    nu             = m;

    % ---------------------------
    % 生成数据（新接口：gen_hggm1(param)）
    % ---------------------------
    param = struct('m',m,'q',q,'nblocks',nblocks,'config',options.config,'var',options.var,...
                   'use_gpu',1,'eigreg',1e-4,'run_bash_mode',0);
    [S,Data,Theta_sim] = gen_hggm1(param); %#ok<NASGU>
    Psi       = S.X;                  % 结构体取 .X 作为矩阵
    Theta_sim = Theta_sim.X;

    % ---------------------------
    % 统一缩放（原版口径）
    % ---------------------------
    scale = sqrt(sum(abs(diag(Psi*Psi')))/length(Psi));
    Psi   = Psi / scale;

    % ---------------------------
    % (A) HG-LASSO (LQA)
    % ---------------------------
    [Theta_hg,~] = hg_lasso_lqa1(Psi,m,A,a,nu,maxiter,param);
    Theta_unb_hg   = 2*Theta_hg - Theta_hg*(Psi)*Theta_hg;
    Theta_var_hg   = sqrt(abs(diag(Theta_hg))*abs(diag(Theta_hg))' + abs(Theta_hg).^2);
    % Rayleigh 矫正（保留矩阵，不做归一化）
    rth = 3.16; 
    Theta_ray_hg = Theta_unb_hg;
    mask_hg = abs(Theta_ray_hg) < (rth/sqrt(m))*(Theta_var_hg - diag(diag(Theta_var_hg)));
    Theta_ray_hg(mask_hg) = 0;

    % ---------------------------
    % (B) PGD-only 源域适配器（Module 6 → 5 → (8)）
    % ---------------------------
    % 对角白化：Σ̃ = D * Σ * D, 其中 D = diag(1./sqrt(diag(Σ)))
    d    = real(diag(Psi)); d(d<=0) = eps;      % 数值安全
    D    = diag(1./sqrt(d));                    % whitening matrix (S^{-1/2})
    Sjj_tilde = D * Psi * D;                    % whitened covariance (correlation-like)

    args = struct();
    args.Sjj_tilde = {Sjj_tilde};               % F=1
    args.K         = 1;                         % 无频域平滑
    args.W         = eye(q);                    % 中性权重
    args.D_src     = {D};                       % 供 recolor 使用（Ω = D * Γ̃ * D）

    opts = struct();
    opts.recolor  = true;                       % 需要源域 Ω
    % 若需与 Jankova 的 a 同标尺，可启用：
    % opts.jankova = struct('override',true,'c',1.0,'q',q,'m',m);

    out_pgd = adapter_pgd_source_only(args, opts);

    % 取源域精度矩阵 Θ_pgd（Ω）
    if ~isempty(out_pgd.Omega_src)
        Theta_pgd = out_pgd.Omega_src{1};
    else
        Theta_pgd = D * out_pgd.Gamma_tilde_star{1} * D; % F=1 与 recolor 等价
    end

    Theta_unb_pgd = 2*Theta_pgd - Theta_pgd*(Psi)*Theta_pgd;
    Theta_var_pgd = sqrt(abs(diag(Theta_pgd))*abs(diag(Theta_pgd))' + abs(Theta_pgd).^2);
    Theta_ray_pgd = Theta_unb_pgd;
    mask_pgd = abs(Theta_ray_pgd) < (rth/sqrt(m))*(Theta_var_pgd - diag(diag(Theta_var_pgd)));
    Theta_ray_pgd(mask_pgd) = 0;

 % ---------------------------
% 画图：每个 trial 占 1 行 × 10 列（前 5: HG；后 5: PGD）
% ---------------------------
figure(fig);  % 确保在主图上绘图（防止外部函数改了当前 figure）
r1 = (trial-1)*10 + (1:5);   % 当前行的前五格（HG）
r2 = (trial-1)*10 + (6:10);  % 当前行的后五格（PGD）

% 1) simulated PCoh（各自一张以便左右对照）
X = Theta_sim; X = X - diag(diag(X)); X = X / max(abs(X(:))+eps);
subplot(ntrial,10,r1(1)); imagesc(abs(X)); title('simulated PCoh'); ylabel('generators'); xlabel('generators');
subplot(ntrial,10,r2(1)); imagesc(abs(X)); title('simulated PCoh'); ylabel('generators'); xlabel('generators');

% 2) hg-lasso PCoh
X = Theta_hg; X = X - diag(diag(X)); X = X / max(abs(X(:))+eps);
subplot(ntrial,10,r1(2)); imagesc(abs(X)); title('hg-lasso PCoh'); ylabel('generators'); xlabel('generators');

% 3) hg-lasso unbiased PCoh
X = Theta_unb_hg; X = X - diag(diag(X)); X = X / max(abs(X(:))+eps);
subplot(ntrial,10,r1(3)); imagesc(abs(X)); title('unbiased PCoh'); ylabel('generators'); xlabel('generators');

% 4) hg-lasso Rayleigh z-stat 直方图 + 理论曲线
X_hg     = Theta_unb_hg; V_hg = Theta_var_hg;
ind_true = find(abs(Theta_sim) > 0); X_hg(ind_true) = 0; V_hg(ind_true) = 0;
ind_var  = find(V_hg/max(V_hg(:)+eps) > 1E-2);
z_hg     = sqrt(m)*abs(X_hg(ind_var))./V_hg(ind_var);

subplot(ntrial,10,r1(4));
h = histogram(z_hg(:),'Normalization','pdf'); bins = h.BinEdges; pdf = 2*bins.*exp(-bins.^2);
hold on; plot(bins,pdf,'LineWidth',2); xlim([0 5]);
title('hg Rayleigh pdf'); xlabel('z-stat'); ylabel('pdf(z)'); legend('z-stat','Rayleigh'); hold off;

% 5) hg-lasso Rayleigh 校正 PCoh（可视化）
X = Theta_ray_hg; X = X - diag(diag(X)); X = X / max(abs(X(:))+eps);
subplot(ntrial,10,r1(5)); imagesc(abs(X)); title('Rayleigh corrected'); ylabel('generators'); xlabel('generators');

% 6) pgd simulated（与左侧保持对照）
X = Theta_sim; X = X - diag(diag(X)); X = X / max(abs(X(:))+eps);
subplot(ntrial,10,r2(1)); imagesc(abs(X)); title('simulated PCoh'); ylabel('generators'); xlabel('generators');

% 7) pgd PCoh
X = Theta_pgd; X = X - diag(diag(X)); X = X / max(abs(X(:))+eps);
subplot(ntrial,10,r2(2)); imagesc(abs(X)); title('pgd PCoh'); ylabel('generators'); xlabel('generators');

% 8) pgd unbiased PCoh
X = Theta_unb_pgd; X = X - diag(diag(X)); X = X / max(abs(X(:))+eps);
subplot(ntrial,10,r2(3)); imagesc(abs(X)); title('unbiased PCoh'); ylabel('generators'); xlabel('generators');

% 9) pgd Rayleigh z-stat 直方图 + 理论曲线
X_pg     = Theta_unb_pgd; V_pg = Theta_var_pgd;
X_pg(ind_true) = 0; V_pg(ind_true) = 0;
ind_var2 = find(V_pg/max(V_pg(:)+eps) > 1E-2);
z_pg     = sqrt(m)*abs(X_pg(ind_var2))./V_pg(ind_var2);

subplot(ntrial,10,r2(4));
h2 = histogram(z_pg(:),'Normalization','pdf'); bins2 = h2.BinEdges; pdf2 = 2*bins2.*exp(-bins2.^2);
hold on; plot(bins2,pdf2,'LineWidth',2); xlim([0 5]);
title('pgd Rayleigh pdf'); xlabel('z-stat'); ylabel('pdf(z)'); legend('z-stat','Rayleigh'); hold off;

% 10) pgd Rayleigh 校正 PCoh（可视化）
X = Theta_ray_pgd; X = X - diag(diag(X)); X = X / max(abs(X(:))+eps);
subplot(ntrial,10,r2(5)); imagesc(abs(X)); title('Rayleigh corrected'); ylabel('generators'); xlabel('generators');

% 美化（对当前主图生效）
try
    load('colormap3'); colormap(cmap); %#ok<NODEF>
catch
    colormap('parula');
end


    % ---------------------------
    % 评估：六种 SPD 距离 + AUC/SENS/SPEC/PREC/F1 + 雷达图
    % ---------------------------
    mask_eval = triu(true(q),1);

    args_hg = struct('Theta_true',Theta_sim,'Theta_est',Theta_hg,'Theta_unb',Theta_unb_hg,'Theta_ray',Theta_ray_hg,'mask',mask_eval);
    opts_eval = struct();
    opts_eval.normalize_mode = 'maxabs';
    opts_eval.spd = struct('symmetrize',true,'project',true,'eps',1e-10);
    opts_eval.alpha = struct('value',0,'use_bcv_bug',false);
    opts_eval.plot = struct('radar',true,'title',sprintf('HG-lasso (q=%d, m=%d)',q,m));
    eval_hg = compute_spd_and_classification_metrics(args_hg, opts_eval);
    if isfield(eval_hg,'fig') && isfield(eval_hg.fig,'radar') && ishghandle(eval_hg.fig.radar)
        saveas(eval_hg.fig.radar, fullfile(output_source, sprintf('Radar_HG_q%d.fig',q)));
        close(eval_hg.fig.radar);
    end

    args_pgd = struct('Theta_true',Theta_sim,'Theta_est',Theta_pgd,'Theta_unb',Theta_unb_pgd,'Theta_ray',Theta_ray_pgd,'mask',mask_eval);
    opts_eval.plot.title = sprintf('PGD-only (q=%d, m=%d)',q,m);
    eval_pgd = compute_spd_and_classification_metrics(args_pgd, opts_eval);
    if isfield(eval_pgd,'fig') && isfield(eval_pgd.fig,'radar') && ishghandle(eval_pgd.fig.radar)
        saveas(eval_pgd.fig.radar, fullfile(output_source, sprintf('Radar_PGD_q%d.fig',q)));
        close(eval_pgd.fig.radar);
    end

    % 记录结果
    k = trial;
    result(k).q = q; result(k).m = m; result(k).Theta_sim = Theta_sim;
    result(k).Theta_hg = Theta_hg; result(k).Theta_unb_hg = Theta_unb_hg; result(k).Theta_var_hg = Theta_var_hg; result(k).Theta_ray_hg = Theta_ray_hg;
    result(k).Theta_pgd = Theta_pgd; result(k).Theta_unb_pgd = Theta_unb_pgd; result(k).Theta_var_pgd = Theta_var_pgd; result(k).Theta_ray_pgd = Theta_ray_pgd;
    result(k).hist_bins_hg = bins; result(k).hist_pdf_hg = pdf;
    result(k).hist_bins_pgd = bins2; result(k).hist_pdf_pgd = pdf2;
    result(k).eval_hg = eval_hg; result(k).eval_pgd = eval_pgd;

    trial = trial + 1;
end

% 收尾：保存图
delete(process_waitbar);
saveas(fig, fullfile(output_source,'PCoh_compare.fig'));
disp(['Saving figure ---->  PCoh_compare.fig  to  ---> ', output_source]);
close(fig);

end
