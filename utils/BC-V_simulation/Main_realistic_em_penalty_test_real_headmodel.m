function Main_realistic_em_penalty_test_real_headmodel(output_sourse)
% === Realistic EM Penalty Test + My Adapter + Metrics (with real_head_model_stimulation) ===

if nargin < 1 || isempty(output_sourse)
    output_sourse = pwd;
end
out_dir = fullfile(output_sourse,'Realistic_em_penalty_test');
if ~isfolder(out_dir), mkdir(out_dir); end
rng(20251105)

for ii = 1:1
    %% ---------------------- 场景与头模路径 ----------------------
    sens_system = 'real_head';  % 使用真实头模流程
    % process_waitbar = waitbar(0,'Please wait... prepare simulation ...');
    %
    % % 真实头模（Brainstorm导出）根目录
    % % 例：D:\result\Cuba2003\2pre\BC-V_Structure\F4XW6TBLHF3Q
    % root_dir = 'D:\result\Cuba2003\2pre\BC-V_Structure\F4XW6TBLHF3Q';  % <=== 修改这里
    %
    % paths.channel_mat   = fullfile(root_dir, 'channel', 'channel.mat');
    % paths.headmodel_mat = fullfile(root_dir, 'leadfield', 'headmodel.mat');
    % paths.surf_mat      = fullfile(root_dir, 'surf', 'surf.mat');
    %
    % %% ---------------------- 仿真配置（交给新方法） ----------------------
    % % q 目标维度按场景给"期望值"，最终会 min(期望, 2p) 以保证可辨识
    % % m 为样本数（等价频域段数×峰数）
    % switch sens_system
    %     case 'small',      q_target = 22;  m = 600;
    %     case 'large',      q_target = 36;  m = 6000;
    %     case 'extraLarge', q_target = 64;  m = 6000;
    %     case 'ultraLarge', q_target = 128; m = 6000;
    %     otherwise,         q_target = 32;  m = 2000;
    % end
    %
    % simopts = struct();
    % simopts.mode            = 'fps';   % 'roi' 更稳（需图谱），'fps' 简洁通用
    % simopts.q               = q_target;
    % simopts.m               = m;
    % simopts.car             = true;    % 平均参考
    % simopts.gamma_diagload  = 0.05;    % 对角加载
    % simopts.bio_noise       = struct('nsrc', 50, 'scale', 0.10);
    % simopts.sens_noise      = struct('scale', 0.10);
    % simopts.truth           = struct('block_num', 4, 'var', 2);  % var=2 表示复变量口径（若 gen_hggm2 支持）
    %
    % %% ---------------------- NEW：调用真实头模仿真 ----------------------
    % out_sim = real_head_model_stimulation(paths, simopts);

    % 1. 设置路径
    paths.channel_mat = 'D:\result\Cuba2003\2pre\BC-V_Structure\F4XW6TBLHF3Q\channel\channel.mat';
    paths.headmodel_mat = 'D:\result\Cuba2003\2pre\BC-V_Structure\F4XW6TBLHF3Q\leadfield\headmodel.mat';
    paths.surf_mat = 'D:\result\Cuba2003\2pre\BC-V_Structure\F4XW6TBLHF3Q\surf\surf.mat';

    % 2. 运行仿真（使用推荐配置）
simulation_real_headmodel('real_head', './data', paths, ...
    'Nseed', 32, 'Nsamp', 1000, 'd0', 1e-2, ...
    'seed_method', 'maxdist', 'random_seed', 2024);

    InverseSolvers(output_sourse, paths, sens_system);
    Results(output_sourse);

    % 停在这里：仅按你的要求保留仿真+求逆+结果三步
    return;


    % 统一命名以兼容后续流程
    Svv         = out_sim.Svv;
    Lvj         = out_sim.Lvj;          % 已是降维后的 p×q 导联矩阵
    Seeders     = out_sim.indms(:);     % 记录（仅用于保存/可视化回填）
    Thetajj_sim = out_sim.Theta_true;   % 源域真值精度
    Sjj_sim     = out_sim.Sjj_true;     % 源域真值协
    p           = out_sim.meta.p;
    q           = out_sim.meta.q;
    m           = out_sim.meta.m;
    H_ref       = out_sim.H;            % 参考矩阵（如需要与其他代码对齐）
    waitbar(0.15,process_waitbar,'Simulation ready. Running estimators ...');

    %% ---------------------- 统一参数（按维度自适应） ----------------------
    param               = struct();
    param.use_gpu       = 1;
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

    % 数值稳谱的地板建议稍大些（大 q 时 1e-2~1e-1 更稳）
    param.eigreg        = 1E-2;   % 由 1e-4 调高
    param.axi           = 1E-2;   % 由 1e-4 调高

    aj                  = sqrt(log(q)/m);
    Ajj_diag            = 0;
    Ajj_ndiag           = 1;
    Ajj                 = Ajj_diag*eye(q)+Ajj_ndiag*(ones(q)-eye(q));
    param.aj            = aj;
    param.Ajj           = Ajj;
    param.Axixi         = eye(p);
    param.Axixi_inv     = eye(p);
    param.ntry          = 0;
    param.prew          = 0;
    param.rth1          = 0.7;
    param.rth2          = 3.16;
    param.method        = 'lqa';

    %% ---------------------- 估计：H-HGGM 三惩罚 ----------------------
    methods = {'higgs-lasso','higgs-ridge','higgs-naive', ...
        'eloreta-hglasso','lcmv-hglasso','my-adapter'};
    K = numel(methods);
    Theta_est = zeros(q,q,K);  % raw
    Sigma_est = zeros(q,q,K);  % 对应源协方差（用于去偏）

    waitbar(0.30,process_waitbar,['H-HGGM penalties ... SS:' sens_system]);
    for k = 1:3
        switch k
            case 1, param.penalty = 1; % l1
            case 2, param.penalty = 2; % l2
            case 3, param.penalty = 0; % naive
        end
        % ~~ 直接用 Lvj（已是 p×q），不再取 Lvj(:,Seeders) ~~
        [Theta_est(:,:,k), Sjj_tmp] = higgs(Svv, Lvj, param);   % 返回 Sjj
        Sigma_est(:,:,k) = Sjj_tmp;
    end

    %% ---------------------- eLORETA + HG-LASSO ----------------------
    waitbar(0.50,process_waitbar,['eLORETA+HG-LASSO ... SS:' sens_system]);
    param.gamma1        = 0.001;
    param.gamma2        = 0.05;
    param.delta_gamma   = 0.001;
    [Theta_est(:,:,4), Sjj_tmp] = eloreta_hg_lasso(Svv, Lvj, param);
    Sigma_est(:,:,4) = Sjj_tmp;

    %% ---------------------- LCMV + HG-LASSO ----------------------
    waitbar(0.65,process_waitbar,['LCMV+HG-LASSO ... SS:' sens_system]);
    param.gamma         = sum(abs(diag(Svv)))/(length(Svv)*100);
    [Theta_est(:,:,5), Sjj_tmp] = lcmv_hg_lasso(Svv, Lvj, param);
    Sigma_est(:,:,5) = Sjj_tmp;

    %% ---------------------- 你的方法（my_test2_adapter_bcvE） ----------------------
    waitbar(0.80,process_waitbar,['my-test2-adapter ... SS:' sens_system]);
    emp_cov_cell = {Svv};

    cfg_my = struct( ...
        'output_format','cell','verbose',true, ...
        'do_support_refit',true, ...
        'warm_method','ssblpp', ...      % 仅用于日志
        'kernel_sigma',3.0, 'lambda3_ratio',0.0, ...
        'use_subspace',true, 'do_scale_L',true, ...
        'sigma2xi', (trace(Svv)/size(Svv,1))*1e-2, ...  % E-step σ^2 初值
        'use_higgs_init', true ...
        );
    [Om_cell, Sig_cell] = my_test2_adapter_bcvE_v2(emp_cov_cell, Lvj, m, cfg_my);

    Theta_est(:,:,6) = 0.5*(Om_cell{1}+Om_cell{1}');
    Sigma_est(:,:,6) = Sig_cell{1};

    %% ---------------------- 去偏 & Rayleigh ----------------------
    % === 在获得 Theta_est(:,:,k) 与 Sigma_est(:,:,k) 之后，做 rth 交叉验证 ===
    % 先准备 Psijj（做一次轻量 E 步得到 ESEC）
    [Psijj_for_cv] = get_psijj_for_cv(Svv, Lvj, param);

    rth_grid = param.rth1 : 0.1 : param.rth2;  % 例如 0.7:0.1:3.16
    for k = 1:K
        T_raw  = 0.5*(Theta_est(:,:,k)+Theta_est(:,:,k)');
        S_raw  = 0.5*(Sigma_est(:,:,k)+Sigma_est(:,:,k)');
        T_unb  = 2*T_raw - T_raw*S_raw*T_raw;                 % 去偏
        Tvar   = sqrt(abs(diag(T_raw))*abs(diag(T_raw))' + abs(T_raw).^2);
        Tvar   = Tvar - diag(diag(Tvar));                     % 去掉对角

        % 基线用 Ridge（稳 SPD）
        Theta_ridge_X = ridge_from_psijj(Psijj_for_cv, param);

        % === 交叉验证 rth ===
        [rth_best, mask_best] = cv_rth_by_partial_llh( ...
            rth_grid, T_unb, Tvar, m, Theta_ridge_X, Psijj_for_cv, param);

        % 应用最优掩膜并做一次 SPD 投影
        Theta_mask_X       = Theta_ridge_X;
        Theta_mask_X(mask_best) = 0;  Theta_mask_X(1:q+1:end) = 0;
        Theta_mask         = higgs_eigendecomposition(Theta_mask_X, param); % 你已有此函数

        % 更新 "rayleigh"/"unb" 两个输出（保持你原有变量名）
        Theta_unb(:,:,k)   = T_unb;
        Theta_ray(:,:,k)   = Theta_mask.X;
    end


    %% ---------------------- 可视化（地图） ----------------------
    figure_partial_coherence_maps = figure('Position',[180,110,1200,560]);
    try, load('colormap3'); colormap(cmap); catch, colormap('hot'); end

    % 1) simulated PCoh（真值）
    X = Thetajj_sim; X = X - diag(diag(X)); X = X/max(abs(X(:))+eps);
    subplot(2,4,1); imagesc(abs(X)); title('simulated PCoh'); ylabel('generators'); xlabel('generators');

    % 2) higgs-lasso (raw)
    X = Theta_est(:,:,1); X = X - diag(diag(X)); X = X/max(abs(X(:))+eps);
    subplot(2,4,2); imagesc(abs(X)); title('higgs-lasso (raw)');

    % 3) higgs-ridge (raw)
    X = Theta_est(:,:,2); X = X - diag(diag(X)); X = X/max(abs(X(:))+eps);
    subplot(2,4,3); imagesc(abs(X)); title('higgs-ridge (raw)');

    % 4) higgs-naive (raw)
    X = Theta_est(:,:,3); X = X - diag(diag(X)); X = X/max(abs(X(:))+eps);
    subplot(2,4,4); imagesc(abs(X)); title('higgs-naive (raw)');

    % 5) eloreta-hglasso (raw)
    X = Theta_est(:,:,4); X = X - diag(diag(X)); X = X/max(abs(X(:))+eps);
    subplot(2,4,5); imagesc(abs(X)); title('eloreta-hglasso (raw)'); ylabel('generators'); xlabel('generators');

    % 6) lcmv-hglasso (raw)
    X = Theta_est(:,:,5); X = X - diag(diag(X)); X = X/max(abs(X(:))+eps);
    subplot(2,4,6); imagesc(abs(X)); title('lcmv-hglasso (raw)');

    % 7) my-adapter (raw)
    X = Theta_est(:,:,6); X = X - diag(diag(X)); X = X/max(abs(X(:))+eps);
    subplot(2,4,7); imagesc(abs(X)); title('my-adapter (raw)');

    % 8) my-adapter (ray)
    X = Theta_ray(:,:,6); X = X - diag(diag(X)); X = X/max(abs(X(:))+eps);
    subplot(2,4,8); imagesc(abs(X)); title('my-adapter (ray)');

    saveas(figure_partial_coherence_maps, fullfile(out_dir, sprintf('partial_coherence_maps_(%s).fig',sens_system)));
    close(figure_partial_coherence_maps);

    %% ---------------------- 评估（SPD距离 + 分类指标） ----------------------
    waitbar(0.92,process_waitbar,['metrics ... SS:' sens_system]);

    if exist('compute_spd_and_classification_metrics.m','file') ~= 2
        warning('compute_spd_and_classification_metrics.m not found on path. Metrics skipped.');
        evals = [];
    else
        evals = struct([]);
        mask_eval = triu(true(q),1);
        is_complex = ~isreal(Svv);  % 自动判断复/实
        for k=1:K
            args_k = struct('Theta_true',Thetajj_sim, ...
                'Theta_est', Theta_est(:,:,k), ...
                'Theta_unb', Theta_unb(:,:,k), ...
                'Theta_ray', Theta_ray(:,:,k), ...
                'mask', mask_eval, 'is_complex', is_complex);
            opts_k = struct('normalize_mode','maxabs', ...
                'spd', struct('symmetrize',true,'project',true,'eps',1e-10), ...
                'alpha', struct('value',0,'use_bcv_bug',false), ...
                'plot', struct('radar',false));
            evals(k).name = methods{k};
            evals(k).out  = compute_spd_and_classification_metrics(args_k, opts_k);
        end

        % ---- 整齐的对比表（raw 指标）----
        Tab = cell(K, 6);  % name + 5 指标
        for k=1:K
            mtr = evals(k).out.metrics.est;
            Tab{k,1} = methods{k};
            Tab{k,2} = mtr.auc;
            Tab{k,3} = mtr.sens;
            Tab{k,4} = mtr.spec;
            Tab{k,5} = mtr.prec;
            Tab{k,6} = mtr.f1;
        end
        T = cell2table(Tab, 'VariableNames', {'method','AUC','SENS','SPEC','PREC','F1'});
        writetable(T, fullfile(out_dir, sprintf('metrics_est_(%s).csv', sens_system)));

        % ---- 跨方法雷达图（raw）----
        radarM = zeros(K,5);
        for k=1:K
            mtr = evals(k).out.metrics.est;
            radarM(k,:) = [mtr.auc, mtr.sens, mtr.spec, mtr.prec, mtr.f1];
        end
        fig_radar = figure('Name','Radar across methods','Position',[100 100 560 560]);
        plot_radar_methods(radarM, methods, {'AUC','SENS','SPEC','PREC','F1'}, ...
            sprintf('Quality (raw) — %s', sens_system));
        saveas(fig_radar, fullfile(out_dir, sprintf('radar_metrics_est_(%s).fig', sens_system)));
        close(fig_radar);

        % ---- 你的方法：est/unb/ray 三态雷达（可选）----
        kmy=find(strcmp(methods,'my-adapter'),1);
        if ~isempty(kmy)
            m_est = evals(kmy).out.metrics.est;
            m_unb = evals(kmy).out.metrics.unb;
            m_ray = evals(kmy).out.metrics.ray;
            M3 = [m_est.auc, m_est.sens, m_est.spec, m_est.prec, m_est.f1; ...
                m_unb.auc, m_unb.sens, m_unb.spec, m_unb.prec, m_unb.f1; ...
                m_ray.auc, m_ray.sens, m_ray.spec, m_ray.prec, m_ray.f1];
            fig_radar2 = figure('Name','Radar my-adapter est/unb/ray','Position',[100 100 560 560]);
            plot_radar_methods(M3, {'est','unb','ray'}, {'AUC','SENS','SPEC','PREC','F1'}, ...
                sprintf('my-adapter states — %s', sens_system));
            saveas(fig_radar2, fullfile(out_dir, sprintf('radar_my_states_(%s).fig', sens_system)));
            close(fig_radar2);
        end
    end

    %% ---------------------- 汇总输出 ----------------------
    result = struct();
    result.sens_system  = sens_system;
    result.Svv          = Svv;
    result.Lvj          = Lvj;
    result.Seeders      = Seeders;
    result.Sjj_sim      = Sjj_sim;
    result.Thetajj_sim  = Thetajj_sim;
    result.Theta_est    = Theta_est;
    result.Sigma_est    = Sigma_est;
    result.Theta_unb    = Theta_unb;
    result.Theta_ray    = Theta_ray;
    result.evals        = evals;
    result.paths = struct( ...
        'maps_fig', fullfile(out_dir, sprintf('partial_coherence_maps_(%s).fig',sens_system)), ...
        'radar_raw_fig', fullfile(out_dir, sprintf('radar_metrics_est_(%s).fig', sens_system)) ...
        );

    save(fullfile(out_dir, sprintf('result_(%s).mat',sens_system)), '-struct', 'result');

    if isvalid(process_waitbar), delete(process_waitbar); end
end

end  % main

% ===== 本地雷达图（跨方法）=====
function plot_radar_methods(M, legends, axes_labels, ttl)
K = size(M,1); L = numel(axes_labels);
th = linspace(0, 2*pi, L+1)'; th(end) = th(1);
M = max(min(M,1),0); % clamp to [0,1]
figure(gcf); clf(gcf);
polaraxes; hold on;
for k=1:K
    r = [M(k,:), M(k,1)];
    polarplot(th, r, 'LineWidth', 2);
end
thetaticks(rad2deg(th(1:end-1)));
thetaticklabels(axes_labels);
rticks(0:0.2:1); rlim([0 1]);
legend(legends,'Location','bestoutside');
title(ttl);
hold off;
end

function Psijj_X = get_psijj_for_cv(Svv, Lvj, param)
% 做一次轻量 E 步拿到 Psijj（ESEC）
% 用 BC-V 的初始化套路取 sigma2xi0 / Sigmajj0
[Svv_s, Lvj_s, ~, ~, sigma2xi0, Sigmajj0] = higgs_initial_values(Svv, Lvj, param);
[~,~,~,Psijj,~] = higgs_expectation(Svv_s, Lvj_s, sigma2xi0, Sigmajj0, param);
Psijj_X = Psijj.X;
end

function Theta_ridge_X = ridge_from_psijj(Psijj_X, param)
% 解析 ridge：先对 Psijj 特征分解，再用 ridge 公式回组装
q  = size(Psijj_X,1);
e  = higgs_eigendecomposition(Psijj_X, param);  % 返回 U,d,X
aj = param.aj;
d_ridge = (sqrt(e.d.^2 + 4*aj^2) - e.d) / (2*aj^2);
Theta_ridge_X = e.U * spdiags(d_ridge,0,q,q) * e.U';
Theta_ridge_X = 0.5*(Theta_ridge_X + Theta_ridge_X');
end

function [rth_best, mask_best] = cv_rth_by_partial_llh(rth_grid, T_unb, Tvar, m, Theta_ridge_X, Psijj_X, param)
% 用 partial jj-likelihood 在 rth 网格上选最优
q = size(T_unb,1);
llh_vals = -inf(numel(rth_grid),1);
for ii = 1:numel(rth_grid)
    rth  = rth_grid(ii);
    mask = (abs(T_unb) < (rth/sqrt(m)) .* Tvar);
    mask(1:q+1:end) = false;  % 不动对角

    Theta_tmp_X       = Theta_ridge_X;
    Theta_tmp_X(mask) = 0;  Theta_tmp_X(1:q+1:end) = 0;

    % SPD 投影（避免非正定）
    Theta_tmp = higgs_eigendecomposition(Theta_tmp_X, param);

    % partial jj-likelihood（与 BC-V 一致）
    llh_vals(ii) = partial_jj_llh(Theta_tmp.X, Psijj_X, param);
end
[~, idx]  = max(llh_vals);
rth_best  = rth_grid(idx);
mask_best = (abs(T_unb) < (rth_best/sqrt(m)) .* Tvar);
mask_best(1:q+1:end) = false;
end

function val = partial_jj_llh(Theta_X, Psijj_X, param)
% 只用 jj 块的似然（带正则项），对应你给的 BC-V 公式
[U,D] = eig(0.5*(Theta_X+Theta_X')); d = real(diag(D));
if any(d<=0 | ~isfinite(d)), val = -inf; return; end
core = sum(log(d)) - sum(abs(sum(Theta_X .* Psijj_X.', 2))); % "trace(Theta*Psijj)"的逐行安全版本

aj  = param.aj; Ajj = param.Ajj;
switch param.penalty
    case 1  % lasso
        pen = aj * sum(abs(Ajj(:) .* Theta_X(:)));
    case 2  % ridge
        pen = 0.5*(aj^2) * sum(sum((Ajj .* Theta_X).^2));
    otherwise
        pen = 0;
end
val = core - pen;
end
