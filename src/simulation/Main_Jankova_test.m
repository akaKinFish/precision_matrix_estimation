function result = Main_Jankova_test(output_source)
%% Main_Jankova_test — HG-LASSO(LQA) vs PGD-only 源域适配器 对比
%
% 说明：在原版 Jankova 仿真基础上，保留 hg_lasso_lqa1 流程，
%       追加调用 adapter_pgd_source_only，从同一 Psi 出发对比两种估计。
%       命名/字段遵循 all_code 习惯；PGD 超参数来自 Module 6（默认），
%       不做手工覆盖（如需同题同标尺，可在 opts.jankova.override=true）。
%
% 依赖：gen_hggm1.m, hg_lasso_lqa1.m, adapter_pgd_source_only.m,
%       module6_hyperparameter_config.m, module5_proximal.m, module8_recoloring.m
%
% --------------------------------------
% 初始化输出目录
% --------------------------------------
output_source = strcat(output_source, filesep,'Jankova_test');
if(~isfolder(output_source))
    mkdir(output_source);
end

sizes  = [30,100,1000];
trial  = 1; ntrial = length(sizes);

fig = figure('Name','HG-LASSO vs PGD-only (Jankova)','Position',[50 50 1800 900]);
process_waitbar = waitbar(0,'Please wait...');

result = struct('q',{},'m',{},'Theta_sim',{}, ...
                'Theta_hg',{},'Theta_unb_hg',{},'Theta_var_hg',{}, ...
                'Theta_pgd',{},'Theta_unb_pgd',{},'Theta_var_pgd',{}, ...
                'hist_bins_hg',{},'hist_pdf_hg',{}, ...
                'hist_bins_pgd',{},'hist_pdf_pgd',{});

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
    param = struct('m',m,'q',q,'nblocks',nblocks,'config',options.config,'var',options.var);
    param.use_gpu = false;
    param.eigreg        = 1E-4;
    param.run_bash_mode = 0;
    [S,Data,Theta_sim] = gen_hggm1(param); %#ok<NASGU>
    % Psi = S;  % 后续变量名保持 Psi 以兼容下游代码
    Psi = S.X;
    Theta_sim = Theta_sim.X;
    % ---------------------------
    % 统一缩放（原版口径）
    % ---------------------------
    scale = sqrt(sum(abs(diag(Psi*Psi')))/length(Psi));
    Psi   = Psi / scale;

    % ---------------------------
    % (A) HG-LASSO (LQA)
    % ---------------------------
    [Theta_hg,llh] = hg_lasso_lqa1(Psi,m,A,a,nu,maxiter,param); %#ok<ASGLU>
    Theta_unb_hg   = 2*Theta_hg - Theta_hg*(Psi)*Theta_hg;
    Theta_var_hg   = sqrt(abs(diag(Theta_hg))*abs(diag(Theta_hg))' + abs(Theta_hg).^2);

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
    % args.A_masks = [] 默认全非对角由内部处理/或不使用
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
        % 若未 recolor，回落到白化域 Γ̃ 并做粗略回色
        Theta_pgd = D * out_pgd.Gamma_tilde_star{1} * D; % 与 recolor 等价（F=1）
    end

    Theta_unb_pgd = 2*Theta_pgd - Theta_pgd*(Psi)*Theta_pgd;
    Theta_var_pgd = sqrt(abs(diag(Theta_pgd))*abs(diag(Theta_pgd))' + abs(Theta_pgd).^2);

    % ---------------------------
    % 画图：每个 q 占两行×5列（上：HG，下：PGD）
    % ---------------------------
    % 行索引
    r1 = (trial-1)*10 + (1:5);   % hg 一行
    r2 = (trial-1)*10 + (6:10);  % pgd 一行

    % 1) simulated PCoh（公用）
    X = Theta_sim; X = X - diag(diag(X)); X = X / max(abs(X(:))+eps);
    subplot(ntrial,10,r1(1)); imagesc(abs(X)); title('simulated PCoh'); ylabel('generators'); xlabel('generators');
    subplot(ntrial,10,r2(1)); imagesc(abs(X)); title('simulated PCoh'); ylabel('generators'); xlabel('generators');

    % 2) hg-lasso PCoh
    X = Theta_hg; X = X - diag(diag(X)); X = X / max(abs(X(:))+eps);
    subplot(ntrial,10,r1(2)); imagesc(abs(X)); title('hg-lasso PCoh'); ylabel('generators'); xlabel('generators');

    % 2') pgd PCoh
    X = Theta_pgd; X = X - diag(diag(X)); X = X / max(abs(X(:))+eps);
    subplot(ntrial,10,r2(2)); imagesc(abs(X)); title('pgd PCoh'); ylabel('generators'); xlabel('generators');

    % 3) hg-lasso unbiased PCoh
    X = Theta_unb_hg; X = X - diag(diag(X)); X = X / max(abs(X(:))+eps);
    subplot(ntrial,10,r1(3)); imagesc(abs(X)); title('unbiased PCoh'); ylabel('generators'); xlabel('generators');

    % 3') pgd unbiased PCoh
    X = Theta_unb_pgd; X = X - diag(diag(X)); X = X / max(abs(X(:))+eps);
    subplot(ntrial,10,r2(3)); imagesc(abs(X)); title('unbiased PCoh'); ylabel('generators'); xlabel('generators');

    % 4) hg-lasso Rayleigh 直方图
    X_hg     = Theta_unb_hg; V_hg = Theta_var_hg;
    ind_true = find(abs(Theta_sim) > 0); X_hg(ind_true) = 0; V_hg(ind_true) = 0;
    ind_var  = find(V_hg/max(V_hg(:)+eps) > 1E-2);
    z_hg     = sqrt(m)*abs(X_hg(ind_var))./V_hg(ind_var);

    subplot(ntrial,10,r1(4));
    h = histogram(z_hg(:),'Normalization','pdf'); bins = h.BinEdges; pdf = 2*bins.*exp(-bins.^2);
    hold on; plot(bins,pdf,'LineWidth',2); xlim([0 5]);
    title('hg Rayleigh pdf'); xlabel('z-stat'); ylabel('pdf(z)'); legend('z-stat','Rayleigh'); hold off;

    % 4') pgd Rayleigh 直方图
    X_pg     = Theta_unb_pgd; V_pg = Theta_var_pgd;
    X_pg(ind_true) = 0; V_pg(ind_true) = 0;
    ind_var2 = find(V_pg/max(V_pg(:)+eps) > 1E-2);
    z_pg     = sqrt(m)*abs(X_pg(ind_var2))./V_pg(ind_var2);

    subplot(ntrial,10,r2(4));
    h2 = histogram(z_pg(:),'Normalization','pdf'); bins2 = h2.BinEdges; pdf2 = 2*bins2.*exp(-bins2.^2);
    hold on; plot(bins2,pdf2,'LineWidth',2); xlim([0 5]);
    title('pgd Rayleigh pdf'); xlabel('z-stat'); ylabel('pdf(z)'); legend('z-stat','Rayleigh'); hold off;

    % 5) hg-lasso Rayleigh 校正 PCoh
    rth = 3.16; X = Theta_unb_hg; X_var = Theta_var_hg;
    rayleigh_mask = find(abs(X) < (rth/sqrt(m))*(X_var - diag(diag(X_var))));
    X(rayleigh_mask) = 0; X = X - diag(diag(X)); X = X / max(abs(X(:))+eps);
    subplot(ntrial,10,r1(5)); imagesc(abs(X)); title('Rayleigh corrected'); ylabel('generators'); xlabel('generators');

    % 5') pgd Rayleigh 校正 PCoh
    X = Theta_unb_pgd; X_var = Theta_var_pgd;
    rayleigh_mask = find(abs(X) < (rth/sqrt(m))*(X_var - diag(diag(X_var))));
    X(rayleigh_mask) = 0; X = X - diag(diag(X)); X = X / max(abs(X(:))+eps);
    subplot(ntrial,10,r2(5)); imagesc(abs(X)); title('Rayleigh corrected'); ylabel('generators'); xlabel('generators');

    % 美化
    try
        load('colormap3'); colormap(cmap); %#ok<NODEF>
    catch
        colormap('parula');
    end

    % 记录结果
    k = trial;
    result(k).q = q; result(k).m = m; result(k).Theta_sim = Theta_sim;
    result(k).Theta_hg = Theta_hg; result(k).Theta_unb_hg = Theta_unb_hg; result(k).Theta_var_hg = Theta_var_hg;
    result(k).Theta_pgd = Theta_pgd; result(k).Theta_unb_pgd = Theta_unb_pgd; result(k).Theta_var_pgd = Theta_var_pgd;
    result(k).hist_bins_hg = bins; result(k).hist_pdf_hg = pdf;
    result(k).hist_bins_pgd = bins2; result(k).hist_pdf_pgd = pdf2;

    trial = trial + 1;
end


% 收尾：保存图
delete(process_waitbar);
saveas(fig, fullfile(output_source,'PCoh_compare.fig'));
disp(['Saving figure ---->  PCoh_compare.fig  to  ---> ', output_source]);
close(fig);

end
