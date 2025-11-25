function [Omega_est, Dsrc_est, Gamma_est, outs] = my_test2_adapter_bcvE_v2(emp_cov_cell, L, T, cfg)
% MY_TEST2_ADAPTER_BCVE
% E-step fully delegated to BC-V's higgs_expectation, then PGD (Module 6/5/8).
%
% Key新增:
%   cfg.rayleigh.enable  (default: true)
%   cfg.rayleigh.rth1    (default: 0.7)
%   cfg.rayleigh.rth2    (default: 3.16)
%   cfg.rayleigh.dr      (default: 0.1)
%   cfg.rayleigh.verbose (default: false)

% ---------- defaults ----------
if nargin<4, cfg = struct(); end
emp_cov_cell = coerce_cov_cell(emp_cov_cell);
F = numel(emp_cov_cell);
[p,n] = size(L);

verbose       = getf(cfg,'verbose',true);
do_scale_L    = getf(cfg,'do_scale_L',true);
do_scale_data = getf(cfg,'do_scale_data',true);
eigreg        = getf(cfg,'eigreg',1e-4);
noise_model   = getf(cfg,'noise_model','scalar'); %#ok<NASGU>
Sigma_xixi0   = getf(cfg,'Sigma_xixi',[]);
Sigmajj0_in   = getf(cfg,'Sigmajj0',[]);
kernel_sigma  = getf(cfg,'kernel_sigma',3.0);
lambda3_ratio = getf(cfg,'lambda3_ratio',0.0);
do_support_refit = getf(cfg,'do_support_refit',true);

% Rayleigh配置
ray = getf(cfg,'rayleigh',struct());
ray.enable  = getf(ray,'enable', true);
ray.rth1    = getf(ray,'rth1', 0.7);
ray.rth2    = getf(ray,'rth2', 3.16);
ray.dr      = getf(ray,'dr',   0.1);
ray.verbose = getf(ray,'verbose', false);

em = getf(cfg,'em',struct());
em.max_estep_iter = getf(em,'max_estep_iter',2);
em.eta_sigma      = getf(em,'eta_sigma',0.3);

% ---------- scale L ----------
if do_scale_L
    sL = sqrt(trace(L*L')/size(L,1));
    if isfinite(sL) && sL>0, L = L/sL; if verbose, fprintf('[scale] L by 1/%.3g\n',sL); end, end
end

% ---------- data pre-scale ----------
if do_scale_data
    gamma_grid = logspace(-4,1,30);
    gopt = zeros(F,1); scl = ones(F,1);
    for f=1:F
        [~, Sjj_e, ~, gamma_opt] = eloreta_simple_(emp_cov_cell{f}, L, gamma_grid, 30, 1e-6, false); %#ok<ASGLU>
        scl(f) = mean(real(diag(Sjj_e))); if ~isfinite(scl(f))||scl(f)<=0, scl(f)=1; end
        emp_cov_cell{f} = emp_cov_cell{f} / scl(f);
        gopt(f) = gamma_opt;
    end
    if verbose, fprintf('[scale] data per-freq by eLORETA mean power, median=%.3g\n', median(scl)); end
end

% ---------- initial sigma2xi ----------
if isempty(Sigma_xixi0)
    trS = 0; for f=1:F, trS = trS + trace(emp_cov_cell{f}); end
    sigma2xi = real(trS/(F*p)) * 1e-2;
else
    if isscalar(Sigma_xixi0)
        sigma2xi = Sigma_xixi0;
    else
        sigma2xi = real(mean(diag(Sigma_xixi0)));
    end
end
sigma2xi = max(sigma2xi, 1e-10);

% ---------- initial Sigmajj prior ----------
Sigmajj0 = struct('X', []);
if ~isempty(Sigmajj0_in)
    if isstruct(Sigmajj0_in) && isfield(Sigmajj0_in,'X')
        Sigmajj0.X = Sigmajj0_in.X;
    else
        Sigmajj0.X = Sigmajj0_in;
    end
else
    Sigmajj0.X = eye(n);
end
Sigmajj0.X = 0.5*(Sigmajj0.X + Sigmajj0.X');

% ---------- BC-V param ----------
param = struct();
param.Ip     = eye(p);
param.eigreg = eigreg;
param.use_gpu = 1;

% ---------- E-step ----------
Psijj  = cell(F,1);
Sjj    = cell(F,1);
Psixixi_cell = cell(F,1);
Sigmajj_post_cell = cell(F,1);

for tE = 1:em.max_estep_iter
    sacc = 0;
    for f=1:F
        Svv_f = emp_cov_cell{f};
        [Sxixi,Psixixi,Sjj_f,Psijj_f,Sigmajj_post,~] = higgs_expectation(Svv_f, L, sigma2xi*ones(p,1), Sigmajj0, param); %#ok<ASGLU>
        Psijj{f}  = Psijj_f.X;
        Sjj{f}    = Sjj_f;
        Psixixi_cell{f} = Psixixi;
        Sigmajj_post_cell{f} = Sigmajj_post.X;
        sacc = sacc + trace(Psixixi);
    end
    sigma2xi_new = real(sacc/(F*p));
    sigma2xi = max( (1-em.eta_sigma)*sigma2xi + em.eta_sigma*sigma2xi_new, 1e-10 );
    Sigmajj0.X = mean_cat_cell_(Sigmajj_post_cell);
    if verbose
        fprintf('[E] iter %d: sigma2xi=%.3g\n', tE, sigma2xi);
    end
end

% ---------- Whitening ----------
pre = module1_preproc_from_covset(Psijj, struct('smoothing_method','moving_average', ...
                                                'loading_factor',1e-6,'min_power',1e-10,'verbose',false));
D_src     = pre.D;
Sjj_tilde = pre.Sigma_tilde;

% ---------- Active set ----------
input_data_m3 = struct();
input_data_m3.whitened_covariances = Sjj_tilde;
input_data_m3.frequencies = 1:F;
act = module3_active_set(input_data_m3, struct('proxy_method','correlation', ...
                         'quantile_level', getf(cfg,'q_act',0.10), 'force_diagonal_active', true, 'verbose', false));
A_masks = arrayfun(@(f) logical(act.combined_active_mask(:,:,f)), 1:F, 'uni', 0);

% ---------- Hyper-parameters ----------
K = make_frequency_kernel(F, kernel_sigma);
K = real(0.5*(K+K')); K = K / max(1,max(sum(K,2)));
W = ones(n); W(1:n+1:end)=0;
input_data_m6 = struct();
        input_data_m6.whitened_covariances = Sjj_tilde;
        input_data_m6.kernel_matrix        = K;
        input_data_m6.weight_matrix        = W;
        input_data_m6.active_set_mask      = {A_masks{:}};
hp = module6_hyperparameter_config(input_data_m6, struct());

lambda1 = getf(cfg,'lambda1', hp.lambda1);
lambda2 = getf(cfg,'lambda2_factor',1.0) * hp.lambda2_suggested;
lambda3 = lambda3_ratio * lambda1;
alpha0  = hp.alpha;

if verbose
    fprintf('[HP] λ1=%.3g, λ2=%.3g, λ3=%.3g, α=%.3g\n', lambda1, lambda2, lambda3, alpha0);
end

% ---------- M-step：PGD ----------
Gamma_init = cell(F,1);
for f=1:F
    St = 0.5*(Sjj_tilde{f}+Sjj_tilde{f}');
    Gamma_init{f} = inv_psd_robust_(St, 1e-8, 1e-12);
end

params5 = struct();
params5.lambda1 = lambda1;
params5.lambda2 = lambda2;
params5.lambda2_suggested = hp.lambda2_suggested;
params5.lambda3 = lambda3;
params5.alpha0  = min(alpha0, 0.5);
params5.max_iter = getf(getf(cfg,'inner',struct()),'max_iter',30);
params5.verbose  = verbose;
params5.active_set_update_freq = 10;
params5.alpha_max = 5.0; params5.alpha_up=1.2; params5.alpha_down=0.6; params5.alpha_grow_patience=1;
params5.obj_improve_tol = 5e-6;
params5.weight_mode = 'hadamard';
params5.use_graph_laplacian = true;
params5.spatial_graph_matrix = eye(n);
params5.spatial_graph_is_laplacian = true;
params5.spatial_weight_mode = 'node';
params5.diag = struct('enable',false);
params5.alpha_min = 1e-5; params5.armijo_c1=1e-5;
params5.backtrack_beta=0.5; params5.max_backtrack_per_iter=25;
params5.backtrack_patience=Inf; params5.lambda2_decay_factor=1.0; params5.lambda2_min=lambda2;
params5.penalize_diagonal = false;
params5.l1_weights = [];
params5.use_single_step = true;

input_data_m5 = struct();
input_data_m5.whitened_covariances = Sjj_tilde;
input_data_m5.initial_precision    = Gamma_init;
input_data_m5.smoothing_kernel     = K;
input_data_m5.weight_matrix        = W;
input_data_m5.active_set_mask      = {A_masks{:}};
input_data_m5.whitening_matrices   = D_src;

[Gamma_tilde, prox_res] = module5_proximal(input_data_m5, params5);

% ---------- NEW: Rayleigh 网格收缩（白化域），PGD之后、refit之前 ----------
ray_info = struct('grid', [], 'best_rth', [], 'obj', []);
if ray.enable
    rgrid = ray.rth1:ray.dr:ray.rth2;
    ray_info.grid = rgrid;
    best_rths = zeros(F,1);
    objs_all  = cell(F,1);

    for f=1:F
        G = 0.5*(Gamma_tilde{f}+Gamma_tilde{f}');
        S = 0.5*(Sjj_tilde{f}+Sjj_tilde{f}');
        [Gbest, best_rth, obj_vec] = rayleigh_grid_shrink_whitened(G, S, T, rgrid);
        Gamma_tilde{f} = Gbest;  % 用最优阈值后的 Γ̃
        best_rths(f) = best_rth;
        objs_all{f}  = obj_vec;
        if ray.verbose
            fprintf('[Rayleigh] f=%d best rth=%.2f, obj=%.4g\n', f, best_rth, max(obj_vec));
        end
    end
    ray_info.best_rth = best_rths;
    ray_info.obj      = objs_all;
end

% ---------- recolor ----------
input_data_m8 = struct();
    input_data_m8.whitened_precision_matrices = Gamma_tilde;
    input_data_m8.whitening_matrices = D_src;
recol = module8_recoloring(input_data_m8, struct());
Omega_src = recol.recolored_precision_matrices;

% ---------- 支撑再拟合（可选） ----------
if do_support_refit
    support_cell = cell(F,1);
    for f=1:F
        G0 = Gamma_tilde{f}; M = abs(G0)>1e-8; M(1:n+1:end)=true; support_cell{f} = M;
    end
    params_refit = params5; params_refit.max_iter = 10; params_refit.lambda2 = 1e-8;
    params_refit.l1_weights = []; params_refit.active_set_update_freq = Inf; params_refit.use_single_step = true;
    Ginit_ref = cell(F,1);
    for f=1:F, Ginit_ref{f} = (Gamma_tilde{f} .* support_cell{f}); Ginit_ref{f} = 0.5*(Ginit_ref{f}+Ginit_ref{f}'); end
    input_data_m5_ref = input_data_m5; input_data_m5_ref.initial_precision = Ginit_ref; input_data_m5_ref.active_set_mask = support_cell;
    [Gamma_tilde_refit, ~] = module5_proximal(input_data_m5_ref, params_refit);
    input_data_m8_refit = struct();
    input_data_m8_refit.whitened_precision_matrices = Gamma_tilde_refit;
    input_data_m8_refit.whitening_matrices = D_src;
    recol2 = module8_recoloring(input_data_m8_refit, struct());
    Omega_src = recol2.recolored_precision_matrices;
end

% ---------- outputs ----------
Gamma_est = Gamma_tilde;
Omega_est = Omega_src;
Dsrc_est  = invert_cell_spd_(Omega_src, 1e-8, 1e-12);

outs = struct();
outs.Psijj = Psijj; outs.Sjj = Sjj; outs.Psixixi = Psixixi_cell;
outs.Sigmajj_post = Sigmajj_post_cell;
outs.hp = hp; outs.prox_res = prox_res;
outs.sigma2xi_final = sigma2xi;
outs.whitening = struct('D', {D_src}, 'Sjj_tilde', {Sjj_tilde});
outs.K = K; outs.W = W;
outs.rayleigh = ray_info;

end % main

% ================= Rayleigh helper =================
function [Gbest, best_rth, obj_vec] = rayleigh_grid_shrink_whitened(Gamma, Swh, T, rgrid)
% 在白化域做 Rayleigh 网格阈值选择，返回最优Γ与rth
% 目标：max logdet(G) - trace(G*Swh)（不加罚）
    q = size(Gamma,1);
    obj_vec   = -inf(numel(rgrid),1);
    Gcand_set = cell(numel(rgrid),1);

    % 去偏与方差
    Tunb = 2*Gamma - Gamma*Swh*Gamma;             % debiased
    dv   = abs(diag(Gamma));
    Tvar = sqrt(dv*dv.' + abs(Gamma).^2);

    for k = 1:numel(rgrid)
        rth = rgrid(k);
        thr = (rth/sqrt(max(T,1))) * (Tvar - diag(diag(Tvar)));  % 只对非对角
        Mzero = (abs(Tunb) < thr);
        % 保留对角，不置零
        Mzero(1:q+1:end) = false;

        Gk = Gamma;
        Gk(Mzero) = 0;
        % 数值保险：Hermitian + SPD 投影
        Gk = 0.5*(Gk+Gk');
        Gk = proj_spd_safe(Gk, 1e-10);

        % 目标：logdet - trace(G*S)
        [ok,ld] = try_chol_logdet(Gk);
        if ~ok
            % 退路：用特征值近似并加微小地板
            [U,D] = eig(0.5*(Gk+Gk')); d = real(diag(D));
            d = max(d, 1e-12*max(d));
            ld = sum(log(d));
        end
        obj_vec(k) = ld - real(trace(Gk*Swh));
        Gcand_set{k} = Gk;
    end

    [~,ix] = max(obj_vec);
    best_rth = rgrid(ix);
    Gbest    = Gcand_set{ix};
end

function [ok,ld] = try_chol_logdet(A)
    ok = true; ld = -inf;
    try
        R = chol(0.5*(A+A'),'lower');
        ld = 2*sum(log(real(diag(R))+eps));
    catch
        ok = false;
    end
end

function G = proj_spd_safe(G, eps_floor)
    G = 0.5*(G+G');
    [U,D] = eig(full(G),'vector'); d = real(D);
    dmax = max(d);
    floor_val = max(eps_floor*max(dmax,1), 1e-12);
    d(d<floor_val) = floor_val;
    G = U*diag(d)*U'; G = 0.5*(G+G');
end

% ================= helpers (原有) =================
function C = coerce_cov_cell(X)
    if iscell(X), C = X(:); return; end
    if isnumeric(X) && ndims(X)==3 && size(X,1)==size(X,2)
        F=size(X,3); C=cell(F,1); for f=1:F, C{f}=X(:,:,f); end; return;
    end
    if isnumeric(X) && ismatrix(X) && size(X,1)==size(X,2)
        C={X}; return;
    end
    error('emp_cov_cell must be cell or p×p or p×p×F');
end
function v=getf(s,f,def), if isfield(s,f)&&~isempty(s.(f)), v=s.(f); else, v=def; end, end
function A = inv_psd_robust_(S, eps_reg, min_ratio)
    S = 0.5*(S+S'); [U,D]=eig(full(S),'vector'); d=real(D);
    dmax = max(d); floor_val = max(min_ratio*max(dmax,eps), 0);
    d(d<floor_val) = floor_val;
    if eps_reg>0, d = (d + eps_reg*dmax)/(1+eps_reg); end
    A = U*diag(1./d)*U'; A = 0.5*(A+A');
end
function M = make_frequency_kernel(F, sigma)
    [I,J]=ndgrid(1:F,1:F); M=exp(-((I-J).^2)/(2*sigma^2)); M=0.5*(M+M); M=M/max(1,max(sum(M,2)));
end
function S = invert_cell_spd_(Omega_cell, eps_reg, min_ratio)
    F=numel(Omega_cell); S=cell(F,1);
    for f=1:F, S{f} = inv_psd_robust_(Omega_cell{f}, eps_reg, min_ratio); end
end
function X = mean_cat_cell_(C)
    F=numel(C); X=zeros(size(C{1})); for f=1:F, X=X+real(0.5*(C{f}+C{f}')); end; X=X/F;
end
function [Tjv,Sjj,W,gamma_opt] = eloreta_simple_(Svv, L, gamma_grid, maxit, tol, verb)
    [p,n]=size(L);
    best.score=Inf; best.T=[]; best.W=[]; best.g=[];
    for k=1:numel(gamma_grid)
        gamma=gamma_grid(k); w=ones(n,1);
        for it=1:maxit
            Winv = diag(1./w);
            A = (L*Winv*L'); A=0.5*(A+A');
            alpha = gamma*trace(A)/p;
            M = chol(A+alpha*eye(p)); M = (M\(M'\eye(p)));
            w_old=w;
            for i=1:n, li=L(:,i); w(i)=sqrt(max(real(li'*M*li), eps)); end
            if norm(w-w_old)/max(1,norm(w_old)) < tol, break; end
        end
        Winv=diag(1./w); A=(L*Winv*L'); A=0.5*(A+A'); alpha = gamma*trace(A)/p;
        M = chol(A+alpha*eye(p)); M=(M\(M'\eye(p)));
        T = Winv*L'*M; Txiv = eye(p) - L*T;
        num = real(trace((Txiv*Svv*Txiv')))/p; den=(real(trace(Txiv))/p)^2 + eps;
        gcv=num/den;
        if gcv<best.score, best.score=gcv; best.T=T; best.W=diag(w); best.g=gamma; end
        if verb && mod(k,10)==1, fprintf('[eLORETA] g=%.3g GCV=%.3g\n',gamma,gcv); end
    end
    Tjv = best.T; W=best.W; gamma_opt=best.g; Sjj = Tjv*Svv*Tjv';
end
