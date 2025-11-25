function result = Main_sim2(output_source)
%% Main_sim2 — Simplified EM Penalty Test + My Adapter 对比（含 metrics 版）
%
% 对比方法：
%   1) H-HGGM / EM：penalty = lasso (1), ridge (2), naive (0)
%   2) eLORETA + hg-lasso
%   3) LCMV + hg-lasso
%   4) VARETA
%   5) 我的方法 my_test2_adapter（EM + PGD-only 源域估计）
%
% 备注：我的方法内部含 support-refit（L1去偏的一种），但不含解析式去偏与 Rayleigh 矫正；
%       因此在调用处（本脚本）对所有方法统一做解析式去偏与 Rayleigh 矫正，并计算指标。
%
% 依赖：gen_hggm2.m, higgs.m, eloreta_hg_lasso.m, lcmv_hg_lasso.m, vareta.m,
%       my_test2_adapter_bcvE.m, compute_spd_and_classification_metrics.m

if nargin < 1 || isempty(output_source)
    output_source = pwd;
end

%% 输出目录
out_dir = fullfile(output_source,'Simplified_em_penalty_test');
if ~isfolder(out_dir), mkdir(out_dir); end

%% -------------------------
%% 仿真配置
%% -------------------------
m                    = 600;   % 样本数
q                    = 22;    % 源个数
p                    = 30;    % 传感器个数
nblocks              = 2;     % 块数
options.config       = 2;     % (2) overlapping blocks (1) nonoverlapping blocks
options.var          = 2;     % (2) complex variable (1) real variable
options.extensions   = [ceil(q/3); ceil(q/3); q - 2*ceil(q/3)]; % patches extensions
options.connections  = [1 2; 2 3];      % patches connections

%% 生成源域真值（Sjj_sim, j_sim, Thetajj_sim）
[Sjj_sim, j_sim, Thetajj_sim] = gen_hggm2(m, q, options);
% 确认 j_sim 形状为 q×m
if size(j_sim,1) ~= q && size(j_sim,2) == q
    j_sim = j_sim.'; % 转为 q×m
end

%% 生成伪 leadfield Lvj（几何同圆布点）
Lvj   = zeros(p,q);
radj  = 60;   radv = 85;
angj  = 2*pi/q; angv = 2*pi/p;
wb = waitbar(0,'Creating pseudo Lead Field (L) ...');
for contv = 1:p
    for contj = 1:q
        if ishandle(wb)
            waitbar(((contv-1)*q+contj)/(p*q), wb);
        end
        vectv            = [radv*cos((contv-1)*angv); radv*sin((contv-1)*angv)];
        vectj            = [radj*cos((contj-1)*angj); radj*sin((contj-1)*angj)];
        r                = vectv - vectj;
        r_unit           = r/sqrt(sum(abs(r).^2));
        miu              = vectj/sqrt(sum(abs(vectj).^2));
        Lvj(contv,contj) = (1/(4*pi))*miu'*r_unit/sqrt(sum(abs(r).^2))^2;
    end
end
if ishandle(wb), delete(wb); end

%% 伪 cortex/头模（用于可视化坐标存档）
vertices        = zeros(q,2);
for contj = 1:q
    vertices(contj,:) = [radj*cos((contj-1)*angj) radj*sin((contj-1)*angj)];
end
faces           = [[1:q]' [2:q 1]'];
cortex.vertices = vertices; cortex.faces = faces;
coor            = zeros(p,2);
for contv = 1:p
    coor(contv,:)    = [radv*cos((contv-1)*angv) radv*sin((contv-1)*angv)];
end
hm_dir = fullfile('simulations','Sim2_h_hggm_simplified_head_model');
if ~isfolder(hm_dir), mkdir(hm_dir); end
save(fullfile(hm_dir,'HeadModel_pseudo'),'cortex','coor');

%% 生成观测数据 v，并加噪
v0              = Lvj * j_sim; % p×m
% 生物噪声（源级）
bionoise        = randn(q,m) + 1i*randn(q,m);
bionoise        = Lvj*bionoise;
bionoise        = sum(abs(v0(:)).^2)^(1/2)*bionoise/sum(abs(bionoise(:)).^2)^(1/2);
% 传感器噪声
sensnoise       = randn(p,m) + 1i*randn(p,m);
sensnoise       = sum(abs(v0(:)).^2)^(1/2)*sensnoise/sum(abs(sensnoise(:)).^2)^(1/2);
% 合成
v               = v0 + 0.1*bionoise + 0.1*sensnoise;
% 经验协方差（偏置型）
Svv             = (v*v')/m;          % 与 nu=m 口径一致
Svv             = (Svv + Svv')/2;    % Hermitian

%% 统一参数结构
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
param.method        = 'lqa';
% LASSO 口径
aj                  = sqrt(log(q)/m);
Ajj_diag            = 0;
Ajj_ndiag           = 1;
Ajj                 = Ajj_diag*eye(q)+Ajj_ndiag*(ones(q)-eye(q));
param.aj            = aj;
param.Ajj           = Ajj;
% 传感器噪声先验（若方法用到）
param.axi           = 1E-4;
param.Axixi         = eye(p);
param.Axixi_inv     = eye(p);
param.ntry          = 0;
param.prew          = 1;
param.rth1          = 0.7;
param.rth2          = 3.16;

%% -------------------------
%% 1) H-HGGM（EM）三种惩罚
%% -------------------------
penalty             = [1 2 0]; % 1:l1 2:l2 0:naive
Thetajj_est         = zeros(q,q,length(penalty) + 4); % 3 (higgs) + 1 (eloreta) + 1 (lcmv) + 1 (my)
Sjj_est             = zeros(q,q,length(penalty) + 4); % 3 (higgs) + 1 (eloreta) + 1 (lcmv) + 1 (vareta) + 1 (my)
llh_outer           = cell(1,length(penalty));

wb2 = waitbar(0,'H-HGGM (EM) with penalties ...');
for k = 1:length(penalty)
    if ishandle(wb2), waitbar(k/length(penalty), wb2); end
    param.penalty   = penalty(k);
    [Thetajj_est(:,:,k), Sjj_est(:,:,k), Tjv_est, llh_outer{k}] = higgs(Svv, Lvj, param); %#ok<ASGLU>
end
if ishandle(wb2), delete(wb2); end

%% -------------------------
%% 2) eLORETA + hg-lasso
%% -------------------------
param.gamma1        = 0.001;
param.gamma2        = 0.05;
param.delta_gamma   = 0.001;
[Thetajj_est(:,:,4), Sjj_est(:,:,4), gamma_grid, gamma, gcv] = eloreta_hg_lasso(Svv, Lvj, param); %#ok<ASGLU>

%% -------------------------
%% 3) LCMV + hg-lasso
%% -------------------------
param.gamma         = sum(abs(diag(Svv)))/(length(Svv)*100);
[Thetajj_est(:,:,5), Sjj_est(:,:,5)] = lcmv_hg_lasso(Svv, Lvj, param);

%% -------------------------
%% 4) VARETA
%% -------------------------
[U,Sv,V] = svd(Lvj, 'econ');
Sv = diag(Sv);
[Sjj_est(:,:,6), lambda, Glambda, reg_param, G] = vareta(U, Sv, V, Svv, 0); %#ok<ASGLU>

%% -------------------------
%% 5) 我的方法（my_test2_adapter_bcvE）
%% -------------------------
emp_cov_cell = {Svv};   % F=1 频
cfg_my = struct();
cfg_my.output_format    = 'cell';
cfg_my.verbose          = false;
cfg_my.do_support_refit = true;   % adapter 内部的 L1 支撑再拟合
cfg_my.warm_method      = 'eloreta';
cfg_my.kernel_sigma     = 3.0;
cfg_my.lambda3_ratio    = 0.3;
cfg_my.use_subspace     = true;
cfg_my.noise_model      = 'scalar';
cfg_my.do_scale_L       = true;   % 与适配器默认一致，提升数值稳健
cfg_my.do_scale_data    = true;

[Omega_cell, Dsrc_cell, Gamma_cell, outs_my] = my_test2_adapter_bcvE(emp_cov_cell, Lvj, m, cfg_my); %#ok<ASGLU>

% 适配器已 recolor 完成：第1个输出=Ω，第2个输出=Σ（稳健 Ω^{-1}）
Theta_my = 0.5*(Omega_cell{1} + Omega_cell{1}');   % 源域精度 Ω
Sigma_my = 0.5*(Dsrc_cell{1} + Dsrc_cell{1}');     % 源域协方差 Σ

% 轻量自检：Θ 与 Σ 是否自洽
cond_err = norm(Theta_my*Sigma_my - eye(q), 'fro')/max(1,q);
if cond_err > 1e-2
    warning('[my-adapter] Σ 与 Ω 配对偏差较大：||ΘΣ-I||_F/q = %.3g（>1e-2）', cond_err);
else
    fprintf('[my-adapter] ||Theta*Sigma - I||_F/q = %.3g\n', cond_err);
end

% 解析式去偏（保存到 result 用）
Theta_unb_my = 2*Theta_my - Theta_my * Sigma_my * Theta_my;

% Rayleigh 矫正（同 Sim1 口径；rth=3.16）
Theta_var_my = sqrt(abs(diag(Theta_my))*abs(diag(Theta_my))' + abs(Theta_my).^2);
rth = 3.16;
Theta_ray_my = Theta_unb_my;
mask_ray = abs(Theta_ray_my) < (rth/sqrt(m))*(Theta_var_my - diag(diag(Theta_var_my)));
Theta_ray_my(mask_ray) = 0;

% 可视化占位：第8幅标题是“PCoh”，应画协方差口径；因此这里放 Σ
Sjj_est(:,:,7)     = Sigma_my;    % 用协方差保持与其他方法的“PCoh”口径一致
Thetajj_est(:,:,7) = Theta_my;    % 同时保存精度矩阵供 metrics

%% -------------------------
%% 指标计算（SPD 距离 + AUC/SENS/SPEC/PREC/F1）
%% -------------------------
method_names = {'higgs_lasso','higgs_ridge','higgs_naive','eloreta_hglasso','lcmv_hglasso','vareta','my_adapter'};
Theta_list = cell(1, numel(method_names));
Sigma_list = cell(1, numel(method_names));

% 1..3: H-HGGM
Theta_list{1} = Thetajj_est(:,:,1); Sigma_list{1} = Sjj_est(:,:,1);
Theta_list{2} = Thetajj_est(:,:,2); Sigma_list{2} = Sjj_est(:,:,2);
Theta_list{3} = Thetajj_est(:,:,3); Sigma_list{3} = Sjj_est(:,:,3);
% 4: eLORETA + hg-lasso
Theta_list{4} = Thetajj_est(:,:,4); Sigma_list{4} = Sjj_est(:,:,4);
% 5: LCMV + hg-lasso
Theta_list{5} = Thetajj_est(:,:,5); Sigma_list{5} = Sjj_est(:,:,5);
% 6: VARETA（仅有 Sjj，需要稳健逆得到 precision）
Sigma_list{6} = Sjj_est(:,:,6);
S6 = 0.5*(Sigma_list{6}+Sigma_list{6}');
[U6,d6] = eig(full(S6),'vector'); d6 = real(d6); d6max = max(d6); d6 = max(d6, 1e-8*max(d6max,1));
Theta_list{6} = 0.5*(U6*diag(1./d6)*U6' + U6*diag(1./d6)*U6');
% 7: 我的方法（用 Theta_my / Sigma_my）
Theta_list{7} = Theta_my;      Sigma_list{7} = Sigma_my;

% 评估参数
opts_eval = struct();
opts_eval.normalize_mode = 'maxabs';
opts_eval.spd = struct('symmetrize',true,'project',true,'eps',1e-8);
opts_eval.alpha = struct('value',0,'use_bcv_bug',false);
opts_eval.plot = struct('radar',false,'title','');

eval_results   = struct();
methods_details = struct();
for ii = 1:numel(method_names)
    Theta_i = 0.5*(Theta_list{ii} + Theta_list{ii}');
    Sigma_i = 0.5*(Sigma_list{ii} + Sigma_list{ii}');
    % 解析式去偏
    Theta_unb_i = 2*Theta_i - Theta_i*Sigma_i*Theta_i;
    % 方差与 Rayleigh 阈值
    dv = abs(diag(Theta_i));
    Theta_var_i = sqrt(dv*dv.' + abs(Theta_i).^2);
    rth_local = 3.16;
    Theta_ray_i = Theta_unb_i;
    mask_ray_i  = abs(Theta_ray_i) < (rth_local/sqrt(m))*(Theta_var_i - diag(diag(Theta_var_i)));
    Theta_ray_i(mask_ray_i) = 0;
    % 评估
    args_eval = struct();
    args_eval.Theta_true = Thetajj_sim;
    args_eval.Theta_est  = Theta_i;
    args_eval.Theta_unb  = Theta_unb_i;
    args_eval.Theta_ray  = Theta_ray_i;
    args_eval.is_complex = true;
    out_eval = compute_spd_and_classification_metrics(args_eval, opts_eval);
    % 汇总
    eval_results.(method_names{ii}) = out_eval;
    methods_details.(method_names{ii}) = struct( ...
        'Theta',Theta_i,'Sigma',Sigma_i,'Theta_unb',Theta_unb_i,'Theta_ray',Theta_ray_i,'Theta_var',Theta_var_i);
end

%% -------------------------
%% 可视化
%% -------------------------
figure_partial_coherence_maps = figure('Position',[182,114,1000,560]);
try
    load('colormap1'); colormap(cmap);
catch
    colormap('hot');
end

% 1) simulated PCoh（真值）
X = Sjj_sim; X = X - diag(diag(X)); X = X/max(abs(X(:))+eps);
subplot(2,4,1); imagesc(abs(X)); title('simulated PCoh'); ylabel('generators'); xlabel('generators');

% 2) higgs-lasso
X = Sjj_est(:,:,1); X = X - diag(diag(X)); X = X/max(abs(X(:))+eps);
subplot(2,4,2); imagesc(abs(X)); title('higgs-lasso PCoh'); ylabel('generators'); xlabel('generators');

% 3) higgs-ridge
X = Sjj_est(:,:,2); X = X - diag(diag(X)); X = X/max(abs(X(:))+eps);
subplot(2,4,3); imagesc(abs(X)); title('higgs-ridge PCoh'); ylabel('generators'); xlabel('generators');

% 4) higgs-naive
X = Sjj_est(:,:,3); X = X - diag(diag(X)); X = X/max(abs(X(:))+eps);
subplot(2,4,4); imagesc(abs(X)); title('higgs-naive PCoh'); ylabel('generators'); xlabel('generators');

% 5) eloreta-hglasso
X = Sjj_est(:,:,4); X = X - diag(diag(X)); X = X/max(abs(X(:))+eps);
subplot(2,4,5); imagesc(abs(X)); title('eloreta-hglasso PCoh'); ylabel('generators'); xlabel('generators');

% 6) lcmv-hglasso
X = Sjj_est(:,:,5); X = X - diag(diag(X)); X = X/max(abs(X(:))+eps);
subplot(2,4,6); imagesc(abs(X)); title('lcmv-hglasso PCoh'); ylabel('generators'); xlabel('generators');

% 7) Vareta
X = Sjj_est(:,:,6); X = X - diag(diag(X)); X = X/max(abs(X(:))+eps);
subplot(2,4,7); imagesc(abs(X)); title('Vareta PCoh'); ylabel('generators'); xlabel('generators');

% 8) My Adapter
X = Sjj_est(:,:,7); X = X - diag(diag(X)); X = X/max(abs(X(:))+eps);
subplot(2,4,8); imagesc(abs(X)); title('my-adapter PCoh'); ylabel('generators'); xlabel('generators');

% 保存图
saveas(figure_partial_coherence_maps, fullfile(out_dir,'partial_coherence_maps.fig'));
disp(['Saving figure ---->  partial_coherence_maps.fig  to  ---> ', out_dir]);
close(figure_partial_coherence_maps);

%% -------------------------
%% 汇总结果
%% -------------------------
result = struct();
result.Svv            = Svv;
result.Lvj            = Lvj;
result.Sjj_sim        = Sjj_sim;
result.Thetajj_sim    = Thetajj_sim;
result.Sjj_est        = Sjj_est;         % q×q×7（第7项为 my-adapter 的 Σ：用于PCoh作图）
result.Thetajj_est    = Thetajj_est;     % q×q×7（第7项为 my-adapter 的 Θ：用于metrics）
result.llh_outer      = llh_outer;       % 仅 H-HGGM 三类
result.my             = struct('Theta',Theta_my,'Sigma',Sigma_my,'Theta_unb',Theta_unb_my,'Theta_ray',Theta_ray_my);
result.paths          = struct('fig', fullfile(out_dir,'partial_coherence_maps.fig'));
result.param          = param;
% 新增：指标
result.eval           = eval_results;
result.methods        = methods_details;

end
