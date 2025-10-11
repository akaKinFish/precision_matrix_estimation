function [Omega_src, Dsrc, Gamma_tilde_star, outs] = my_test2_adapter(emp_covariance, L, T, cfg, Omega_true)
% MY_TEST2_ADAPTER  严格复刻 test2 管线（统一接口版）
%
% 入参
%   emp_covariance : {F×1} / p×p×F / 单张 p×p（传感器域经验协方差）
%   L              : p×n leadfield
%   T              : 标称样本数
%   cfg            : 可选配置（详见内部）
%   Omega_true     : {F×1} GT 源域 precision（可选；不推荐用于初值）
%
% 出参
%   Omega_src         : {F×1} 最终源域 precision（recolor + refit 后）
%   Dsrc              : {F×1} 白化矩阵（用于评估/可视化）
%   Gamma_tilde_star  : {F×1} 白化域 precision（最终）
%   outs              : 结构体，包含所有中间产物

if nargin < 4 || isempty(cfg), cfg = struct(); end
if nargin < 5, Omega_true = []; end

%% ========== 0) 统一数据形态 ==========
emp_cov_cell = coerce_cov_cell_(emp_covariance);
F = numel(emp_cov_cell);
p = size(L,1);  
n = size(L,2);

%% ========== 1) 构造 E 步输入（scalar struct） ==========
prior_cell = repmat({eye(n)}, F, 1);

estep_in = struct();
estep_in.leadfield_matrix         = L;
estep_in.empirical_covariances    = emp_cov_cell;
estep_in.source_prior_covariances = prior_cell;
estep_in.noise_covariance         = pick_(cfg,'Sigma_xixi', eye(p));
estep_in.frequencies              = 1:F;

assert(isstruct(estep_in) && isscalar(estep_in));
assert(iscell(estep_in.empirical_covariances) && numel(estep_in.empirical_covariances)==F);

%% ========== [ADD-1] sSSBL风格：尺度统一 ==========
% 1) Leadfield 尺度归一化
scaleLvj = sqrt(trace(L*L')/p);
if isfinite(scaleLvj) && scaleLvj>0
    L = L / scaleLvj;
    estep_in.leadfield_matrix = L;
    if pick_(cfg,'verbose',false)
        fprintf('[SCALE] Leadfield scaled by 1/%.3g\n', scaleLvj);
    end
end

% 2) 数据尺度归一化（关键补充）
gamma_grid = logspace(-4, 1, 30);
opts_el = struct('maxit', 50, 'tol', 1e-6, 'verbose', false);
warm_tmp = eloreta_warmstart_from_covs(emp_cov_cell, L, gamma_grid, opts_el);

scaleJ = zeros(F,1);
for f=1:F
    Sjj_e = warm_tmp.Sjj_e{f};
    scaleJ(f) = mean(abs(diag(Sjj_e)));
    if ~isfinite(scaleJ(f)) || scaleJ(f)<=0
        scaleJ(f) = 1;
    end
    emp_cov_cell{f} = emp_cov_cell{f} / scaleJ(f);
end
estep_in.empirical_covariances = emp_cov_cell;

if pick_(cfg,'verbose',false)
    fprintf('[SCALE] Data scaled per-freq, median factor=%.3g\n', median(scaleJ));
end

%% ========== 1.5) Warm start（关键补充） ==========
warm_method = pick_(cfg, 'warm_method', 'ssblpp');
switch lower(warm_method)
    case 'ssblpp'
        opts_ws = struct(); 
        opts_ws.m = T;
        warm = warmstart_ssblpp_from_covs(emp_cov_cell, L, opts_ws);
    case 'eloreta'
        warm = warm_tmp;  % 复用前面的
    otherwise
        error('Unknown warm_method: %s', warm_method);
end
Omega_prev = warm.Omega_init; 
Sjj_prev   = warm.Sjj_e;

%% ========== 2) E-step（Module 2） ==========
E = module2_estep(estep_in, struct( ...
    'ensure_hermitian',true, 'ensure_real_diag',true, ...
    'ensure_psd',true, 'psd_tol',1e-10, 'diag_loading',1e-10));
Sjj_hat = E.source_second_moments;

%% ========== [ADD-5] Noise M-step + node variance（关键补充） ==========
% VARETA 子空间
[U,Sv,V] = svd(L, 'econ');
sing2 = diag(Sv).^2; 
cum = cumsum(sing2)/sum(sing2);
r = find(cum>=0.99, 1, 'first'); 
if isempty(r), r = size(Sv,1); end
Ur = U(:,1:r); 
Sr = Sv(1:r,1:r); 
Vr = V(:,1:r);

I_p = eye(p); 
S_xixi_accum = zeros(p,p,'like',estep_in.noise_covariance);
node_var_cell = cell(F,1);

for f = 1:F
    Sigma_prior_f = estep_in.source_prior_covariances{f};
    Sigma_prior_f = 0.5*(Sigma_prior_f + Sigma_prior_f');
    Omega_prior_f = inv_psd_robust_(Sigma_prior_f, 1e-8, 1e-12);
    
    % 子空间版后验协方差
    SigInvUr = (estep_in.noise_covariance \ Ur);
    M = Ur' * SigInvUr;
    A_f = Omega_prior_f + Vr * (Sr * (M * Sr)) * Vr';
    A_f = (A_f + A_f')/2;
    Sigma_post_f = inv_psd_robust_(A_f, 1e-8, 1e-12);
    node_var_cell{f} = real(diag(Sigma_post_f));
    
    % 残差噪声协方差
    LinvX = L' / estep_in.noise_covariance;
    T_jv_f = Sigma_post_f * LinvX;
    T_xi_v_f = I_p - L * T_jv_f;
    S_res_f = T_xi_v_f * emp_cov_cell{f} * T_xi_v_f' + L * Sigma_post_f * L';
    S_res_f = 0.5*(S_res_f + S_res_f');
    S_xixi_accum = S_xixi_accum + S_res_f;
end

Sigma_xixi_new = S_xixi_accum / F;
Sigma_xixi_new = 0.5*(Sigma_xixi_new + Sigma_xixi_new');
ridge = 1e-10 * trace(Sigma_xixi_new)/p;
Sigma_xixi_new = Sigma_xixi_new + ridge * eye(p,'like',Sigma_xixi_new);

% 平滑更新
eta = pick_(cfg, 'noise_eta', 0.3);
Sigma_xixi = (1-eta)*estep_in.noise_covariance + eta*Sigma_xixi_new;
[~,chol_flag] = chol(0.5*(Sigma_xixi + Sigma_xixi'),'lower');
if chol_flag ~= 0
    Sigma_xixi = Sigma_xixi + 1e-6*eye(p,'like',Sigma_xixi);
end
estep_in.noise_covariance = Sigma_xixi;

if pick_(cfg,'verbose',false)
    fprintf('[NOISE] trace(Sigma_xixi)=%.3g\n', trace(Sigma_xixi));
end

%% ========== 3) Whitening（Module 1） ==========
pre = module1_preproc_from_covset(Sjj_hat, struct( ...
    'smoothing_method','moving_average', ...
    'loading_factor', 1e-6, ...
    'min_power', 1e-10, ...
    'verbose', false));
Dsrc     = pre.D;              % ⭐ 输出变量1
Sjj_tilde = pre.Sigma_tilde;

%% ========== 4) Active set（Module 3） ==========
input_m3 = struct();
input_m3.whitened_covariances = Sjj_tilde;
input_m3.frequencies = 1:F;

act = module3_active_set(input_m3, struct( ...
    'proxy_method','correlation', ...
    'quantile_level', pick_(cfg,'active_quantile',0.10), ...
    'force_diagonal_active', true, ...
    'verbose', false));
A_masks = arrayfun(@(f) logical(act.combined_active_mask(:,:,f)), 1:F, 'uni', 0);

%% ========== 5) Hyperparameters（Module 6） ==========
K = make_frequency_kernel_(F, pick_(cfg,'kernel_sigma',3.0));
K = real(0.5*(K + K')); 
K = max(K, 0);
row_sums = sum(K,2); 
maxrow = max(row_sums);
if maxrow > 0, K = K / maxrow; end

W = make_uniform_weight_(n);

input_m6 = struct();
input_m6.whitened_covariances = Sjj_tilde;
input_m6.kernel_matrix = K;
input_m6.weight_matrix = W;
input_m6.active_set_mask = {A_masks{:}};

hp = module6_hyperparameter_config(input_m6, struct('use_gershgorin', true));

lambda1 = pick_(cfg,'lambda1', hp.lambda1);
lambda2 = pick_(cfg,'lambda2', hp.lambda2_suggested);
lambda3_ratio = pick_(cfg,'lambda3_ratio', 0.3);
lambda3 = lambda3_ratio * lambda1;
alpha0  = pick_(cfg,'alpha0', hp.alpha);

%% ========== 6) Initial precision（白化域） ==========
Gamma_init = transport_init_(Omega_prev, Dsrc, Sjj_tilde);

%% ========== [ADD-3] 自适应 L1 权重（关键补充） ==========
l1w_cell = cell(F,1);
for f=1:F
    wi = 1 ./ (node_var_cell{f} + 1e-8);  
    wi = wi / median(wi);
    
    % 基于前一步的度数
    if isempty(Omega_prev)
        Om_for_deg = warm.Omega_init{f};
    else
        Om_for_deg = Omega_prev{f};
    end
    sup = triu(abs(Om_for_deg) > 1e-8, 1);
    deg = sum(sup | sup', 2);
    Dg = 1 ./ sqrt((deg + 1) * (deg + 1).');
    
    Wnode = sqrt(wi * wi.');
    Wpen = Wnode .* Dg; 
    Wpen(1:n+1:end) = 0;
    l1w_cell{f} = Wpen;
end

%% ========== 7) Module 5: Proximal（白化域） ==========
input_m5 = struct();
input_m5.whitened_covariances = Sjj_tilde;
input_m5.initial_precision = Gamma_init;
input_m5.smoothing_kernel = K;
input_m5.weight_matrix = W;
input_m5.active_set_mask = {A_masks{:}};
input_m5.whitening_matrices = Dsrc;

params5 = struct( ...
    'lambda1', lambda1, ...
    'lambda2', lambda2, ...
    'lambda2_suggested', hp.lambda2_suggested, ...
    'alpha0', alpha0, ...
    'max_iter', pick_(cfg,'max_iter',30), ...
    'verbose', false, ...
    'alpha_max', 5.0, ...
    'alpha_up', 1.2, ...
    'alpha_down', 0.6, ...
    'alpha_grow_patience', 1, ...
    'obj_improve_tol', 5e-6, ...
    'weight_mode', 'hadamard', ...
    'use_graph_laplacian', true, ...
    'penalize_diagonal', false, ...
    'l1_weights', {l1w_cell}, ...
    'use_single_step', true, ...
    'diag', struct('enable', false));

% 空间正则化
if lambda3 > 0
    % 如果有空间图，使用；否则用占位符
    if isfield(cfg,'spatial_laplacian') && ~isempty(cfg.spatial_laplacian)
        Lsp = cfg.spatial_laplacian;
    else
        L_raw = laplacian_placeholder_(n);
        [Lsp, ~] = normalize_graph_laplacian_(L_raw, 'spectral');
    end
    
    params5.lambda3 = lambda3;
    params5.spatial_graph_matrix = Lsp;
    params5.spatial_graph_is_laplacian = true;
    params5.spatial_weight_mode = 'node';
end

[Gamma_tilde_star, proximal_results] = module5_proximal(input_m5, params5);
% ⭐ Gamma_tilde_star 是输出变量2

%% ========== [ADD-4] Support-refit 去偏（关键补充） ==========
if pick_(cfg,'enable_refit', true)
    support_cell = cell(F,1);
    for f=1:F
        Gf = Gamma_tilde_star{f};
        M = abs(Gf) > 1e-8; 
        M(1:n+1:end) = true;  % 保持对角线
        support_cell{f} = M;
    end
    
    Gamma_init_refit = cell(F,1);
    for f=1:F
        G0 = Gamma_tilde_star{f};
        G0(~support_cell{f}) = 0;
        Gamma_init_refit{f} = (G0 + G0')/2;
    end
    
    input_m5_refit = input_m5;
    input_m5_refit.initial_precision = Gamma_init_refit;
    input_m5_refit.active_set_mask = support_cell;
    
    params_refit = params5;
    params_refit.max_iter = pick_(cfg,'refit_maxiter',10);
    params_refit.lambda2 = 1e-8;
    params_refit.l1_weights = [];
    params_refit.active_set_update_freq = Inf;
    
    [Gamma_tilde_star, ~] = module5_proximal(input_m5_refit, params_refit);
    % ⭐ 更新 Gamma_tilde_star
end

%% ========== 8) Module 8: Recoloring（反白化） ==========
input8 = struct();
input8.whitened_precision_matrices = Gamma_tilde_star;
input8.whitening_matrices = Dsrc;
recol = module8_recoloring(input8, struct());
Omega_src = recol.recolored_precision_matrices;  % ⭐ 输出变量0（主输出）
Omega_src = cellfun(@(A) (A+A')/2, Omega_src, 'uni', 0);

%% ========== 输出 ==========
if nargout >= 4
    outs = struct();
    outs.estep_results = E;
    outs.Sjj_hat = Sjj_hat;
    outs.preprocessing = pre;
    outs.active_results = act;
    outs.hyperparams = hp;
    outs.proximal_results = proximal_results;
    outs.recoloring_results = recol;
    outs.warm = warm;
    outs.node_var_cell = node_var_cell;
    outs.scaleJ = scaleJ;
    outs.K = K;
    outs.W = W;
else
    outs = struct();  % 空结构体
end

end

%% ==================== Helpers ====================
function C = coerce_cov_cell_(X)
    if iscell(X), C = X(:); return; end
    if isnumeric(X) && ndims(X)==3 && size(X,1)==size(X,2)
        F = size(X,3); C = cell(F,1);
        for f=1:F, C{f} = (X(:,:,f)+X(:,:,f)')/2; end
        return;
    end
    if isnumeric(X) && ismatrix(X) && size(X,1)==size(X,2)
        C = {(X+X')/2}; return;
    end
    error('emp_covariance: unsupported format');
end

function v = pick_(s, key, dflt)
    if isfield(s,key) && ~isempty(s.(key))
        v = s.(key);
    else
        v = dflt;
    end
end

function K = make_frequency_kernel_(F, sigma)
    [I,J] = ndgrid(1:F,1:F);
    K = exp(-((I-J).^2)/(2*sigma^2));
    K = (K + K')/2;
end

function W = make_uniform_weight_(n)
    W = ones(n); 
    W(1:n+1:end) = 0;
end

function L = laplacian_placeholder_(n)
    A = ones(n) - eye(n); 
    d = sum(A,2); 
    L = diag(d) - A;
end

function [Lnorm, info] = normalize_graph_laplacian_(L, mode)
    if nargin<2, mode='spectral'; end
    L = (L+L')/2; 
    ev = eig(full(L));
    info.min_eig_before = min(real(ev)); 
    info.max_eig_before = max(real(ev));
    switch lower(mode)
        case 'spectral'
            s = max(1, info.max_eig_before); 
            Lnorm = L / s;
        otherwise
            Lnorm = L;
    end
    ev2 = eig(full(Lnorm)); 
    info.min_eig_after = min(real(ev2)); 
    info.max_eig_after = max(real(ev2));
end

function Gamma_init = transport_init_(Omega_prev, D_src, Sjj_tilde)
    F = numel(D_src); 
    Gamma_init = cell(F,1);
    if isempty(Omega_prev)
        for f=1:F
            St = (Sjj_tilde{f} + Sjj_tilde{f}')/2;
            [U, D] = eig(full(St), 'vector'); 
            d = real(D); 
            d = max(d, 1e-10);
            G = U * diag(1./d) * U'; 
            G = (G + G')/2;
            Gamma_init{f} = G;
        end
        return;
    end
    for f=1:F
        D = D_src{f}; 
        Gamma_init{f} = (D \ Omega_prev{f}) / D; 
        Gamma_init{f} = (Gamma_init{f} + Gamma_init{f}')/2;
    end
end

function Om = inv_psd_robust_(A, eps_reg, min_ratio)
    A = (A + A')/2; 
    [V,D] = eig(A); 
    d = real(diag(D)); 
    dmax = max(d);
    floor_val = max(min_ratio * max(dmax, eps), 0);
    d(d<floor_val) = floor_val; 
    if eps_reg > 0
        d = (d + eps_reg * dmax) / (1 + eps_reg); 
    end
    Om = V * diag(1./d) * V'; 
    Om = (Om + Om')/2;
end