function [] = InverseSolvers(output_source, paths, sens_system)

%% Run Inverse Solvers
%% Loading Simulation Substrate or Real Data.

load('./data/Pseudorand_Net.mat')
%%
if strcmp(sens_system,'real_head') == 1
    headmodel_data = load(paths.headmodel_mat);

    % 提取导联场矩阵
    if isfield(headmodel_data, 'HeadModel') && isfield(headmodel_data.HeadModel, 'Gain')
        K = headmodel_data.HeadModel.Gain;
    elseif isfield(headmodel_data, 'Gain')
        K = headmodel_data.Gain;
    else
        error('头模型文件中未找到 Gain 矩阵');
    end

end
%%
subject    = [1]; %User defined, pick numbers from 1-number of Lead Fields
LeadFields = LeadFields(1,subject);
% 'Vsim' is a cell array containing Simulated or Real Data at different conditions,
% every cell contains a tensor of (Number of Sensors, Time points or Frequencies, Subjects X Samples)
%% Run h-hggm
sol_higgs = InverseSolver_higgs(Svv_sim,LeadFields,Seeders_sim,Nsamp,sens_system);

save(strcat('result',filesep,'Solutions_higgs'),'sol_higgs', '-v7.3')

%% Run my-adapter (external debias + Rayleigh) and save final
try
    if exist('my_test2_adapter_bcvE','file') == 2
        Nsim = size(Svv_sim, 2);
        q = size(Seeders_sim, 1);
        p = size(LeadFields{1}, 1);
        Lvj0 = LeadFields{1};

        % parameters for debias + CV-thresholding
        m = Nsamp;
        param = struct();
        param.use_gpu = 1;          % for BC-V E-step helpers
        param.eigreg  = 1e-2;       % SPD safety for eigendecomposition
        param.rth1    = 0.7;
        param.rth2    = 3.16;
        param.penalty = 2;          % ridge baseline in partial llh
        aj            = sqrt(log(q)/max(1,m));
        Ajj           = ones(q) - eye(q);   % only off-diagonal penalized
        param.aj      = aj;
        param.Ajj     = Ajj;
        param.m = m;
        param.q = q;
        param.p = p;
        param.Op = ones(p,1);
        param.Iq = eye(q);
        param.Ip = eye(p);
        param.axi           = 1E-4;
        param.Axixi         = eye(p);
        param.Axixi_inv     = eye(p);
        param.ntry          = 0;
        Theta_raw_all = nan(q,q,Nsim);
        Sigma_raw_all = nan(q,q,Nsim);
        Theta_unb_all = nan(q,q,Nsim);
        Theta_ray_all = nan(q,q,Nsim);
        rth_best_vec  = nan(Nsim,1);


        for sim = 1:Nsim
            try
                Svv_this   = Svv_sim{1, sim}{1};
                seeders    = Seeders_sim(:, sim);
                L_sel      = Lvj0(:, seeders);
                T          = Nsamp;

                % run adapter (no internal Rayleigh)
                cfg_my = struct('verbose', false);
                [Omega_cell, Dsrc_cell, ~, ~] = my_test2_adapter_bcvE({Svv_this}, L_sel, T, cfg_my);

                if iscell(Omega_cell) && ~isempty(Omega_cell)
                    Theta_raw = 0.5*(Omega_cell{1} + Omega_cell{1}');
                    if iscell(Dsrc_cell) && ~isempty(Dsrc_cell)
                        Sigma_raw = 0.5*(Dsrc_cell{1} + Dsrc_cell{1}');
                    else
                        % fallback if Dsrc not returned: use robust inverse
                        Sigma_raw = pinv(0.5*(Theta_raw+Theta_raw'));
                    end

                    % compute Psijj for CV in jj-likelihood
                    Psijj_for_cv = get_psijj_for_cv_local(Svv_this, L_sel, param);

                    % debias in original domain
                    T_unb = 2*Theta_raw - Theta_raw*Sigma_raw*Theta_raw;

                    % Rayleigh variance proxy and grid search
                    Tvar = sqrt(abs(diag(Theta_raw))*abs(diag(Theta_raw))' + abs(Theta_raw).^2);
                    Tvar = Tvar - diag(diag(Tvar));
                    Theta_ridge_X = ridge_from_psijj_local(Psijj_for_cv, param);
                    rth_grid = param.rth1 : 0.1 : param.rth2;
                    [rth_best, mask_best] = cv_rth_by_partial_llh_local( ...
                        rth_grid, T_unb, Tvar, m, Theta_ridge_X, Psijj_for_cv, param);

                    Theta_mask_X       = Theta_ridge_X;
                    Theta_mask_X(mask_best) = 0;  Theta_mask_X(1:q+1:end) = 0;
                    Theta_mask         = higgs_eigendecomposition(Theta_mask_X, param);

                    % collect
                    Theta_raw_all(:,:,sim) = Theta_raw;
                    Sigma_raw_all(:,:,sim) = Sigma_raw;
                    Theta_unb_all(:,:,sim) = T_unb;
                    Theta_ray_all(:,:,sim) = Theta_mask.X;
                    rth_best_vec(sim)      = rth_best;
                end
            catch MEi
                warning('my-adapter failed on sim %d: %s', sim, MEi.message);
            end
        end

        sol_my_adapter = struct();
        sol_my_adapter.name        = 'my-adapter';
        sol_my_adapter.Theta       = Theta_ray_all;   % final (Rayleigh)
        sol_my_adapter.Theta_unb   = Theta_unb_all;   % optional
        sol_my_adapter.rth_best    = rth_best_vec;
        sol_my_adapter.Seeders_sim = Seeders_sim;
        sol_my_adapter.Nsamp       = Nsamp;
        save(strcat('result',filesep,'Solutions_my_adapter'),'sol_my_adapter', '-v7.3');
    else
        warning('my_test2_adapter_bcvE.m not found on path. Skipping my-adapter solutions.');
    end
catch ME
    warning('my-adapter run encountered an error: %s', ME.message);
end

end

function Psijj_X = get_psijj_for_cv_local(Svv, Lvj, param)
% 做一次轻量 E 步拿到 Psijj（ESEC）
% 用 BC-V 的初始化套路取 sigma2xi0 / Sigmajj0
param.prew = 0;
[Svv_s, Lvj_s, ~, ~, sigma2xi0, Sigmajj0] = higgs_initial_values(Svv, Lvj, param);
[~,~,~,Psijj,~] = higgs_expectation(Svv_s, Lvj_s, sigma2xi0, Sigmajj0, param);
Psijj_X = Psijj.X;
end

function Theta_ridge_X = ridge_from_psijj_local(Psijj_X, param)
% 解析 ridge：先对 Psijj 特征分解，再用 ridge 公式回组装
q  = size(Psijj_X,1);
e  = higgs_eigendecomposition(Psijj_X, param);  % 返回 U,d,X
aj = param.aj;
d_ridge = (sqrt(e.d.^2 + 4*aj^2) - e.d) / (2*aj^2);
Theta_ridge_X = e.U * spdiags(d_ridge,0,q,q) * e.U';
Theta_ridge_X = 0.5*(Theta_ridge_X + Theta_ridge_X');
end

function [rth_best, mask_best] = cv_rth_by_partial_llh_local(rth_grid, T_unb, Tvar, m, Theta_ridge_X, Psijj_X, param)
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
    llh_vals(ii) = partial_jj_llh_local(Theta_tmp.X, Psijj_X, param);
end
[~, idx]  = max(llh_vals);
rth_best  = rth_grid(idx);
mask_best = (abs(T_unb) < (rth_best/sqrt(m)) .* Tvar);
mask_best(1:q+1:end) = false;
end

function val = partial_jj_llh_local(Theta_X, Psijj_X, param)
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
