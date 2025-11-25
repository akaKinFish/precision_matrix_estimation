function estep_results = module2_estep_main(input_data, estep_params)
% MODULE2_ESTEP_MAIN - E-step (BC-V aligned), main 只做“叫用+数据换壳”
% 必要输入字段：
%   input_data.leadfield_matrix           (p×n)   L
%   input_data.empirical_covariances      {F}     Svv(f)
%   input_data.source_prior_covariances   {F}|n×n Σjj_prior (可全相同)
%   input_data.noise_covariance           p×p     Σξξ（标量模型时可给 σ²I）
%   input_data.frequencies                1×F
%
% 关键输出（供 M 步/诊断/对白化等使用）：
%   .posterior_source_covariances         {F} Σ_post(f)
%   .transfer_functions                   {F} Tjv(f) = Σ_post L' Σξ^{-1}
%   .effective_source_second_moments      {F} Psijj(f) = Σ_post + T Svv T'
%   .effective_residual_second_moments    {F} Psixixi(f) = Sξξ + L Σ_post L'
%   .initial_precision_matrices           {F} Ω_init(f) = inv_robust(Σ_post)
%   .whitening.D                          {F} D_src(f)
%   .whitening.Sigma_tilde                {F} Sjj_tilde(f) = D Psijj D
%   .noise.sigma2_scalar_final            1×1 σ²（标量）| []（非标量时）
%   .stats                                诊断

% -------- Input checks --------
if nargin < 1, error('module2_estep_main:insufficient_input','input_data is required'); end
if nargin < 2, estep_params = struct(); end
need = {'leadfield_matrix','empirical_covariances','source_prior_covariances','noise_covariance','frequencies'};
for k = 1:numel(need)
    if ~isfield(input_data, need{k})
        error('module2_estep_main:missing_field','Missing input_data.%s', need{k});
    end
end

% -------- Unpack & normalize --------
L           = input_data.leadfield_matrix;     % p×n
Svv_cell    = input_data.empirical_covariances;% {F}
Sigma_prior = input_data.source_prior_covariances;
Sigma_xixi  = input_data.noise_covariance;     % p×p
freqs       = input_data.frequencies;

if ~iscell(Svv_cell),    Svv_cell = {Svv_cell}; end
if ~iscell(Sigma_prior), Sigma_prior = {Sigma_prior}; end

[p,n] = size(L);
F = numel(Svv_cell);
if numel(Sigma_prior)==1 && F>1, Sigma_prior = repmat(Sigma_prior,F,1); end
if numel(Sigma_prior)~=F
    error('module2_estep_main:prior_count_mismatch','#priors(%d)~=F(%d)', numel(Sigma_prior), F);
end
if ~isequal(size(Sigma_xixi),[p p])
    error('module2_estep_main:noise_size','noise_covariance must be %dx%d',p,p);
end
if ~isvector(freqs) || numel(freqs)~=F, freqs = 1:F; end

% -------- Params (含“标量噪声内循环”) --------
P = struct('eig_floor_ratio',1e-12, 'eps_reg',1e-8, ...
           'verbose', false, ...
           'noise_model','scalar', ...     % 'scalar'|'diag'|'full' （本文件实现了 scalar 的内循环）
           'estep_inner_iters', 2, ...
           'eta_sigma', 0.3, ...
           'do_whitening', true, ...
           'whitening_opts', struct('smoothing_method','moving_average','loading_factor',1e-6,'min_power',1e-10,'verbose',false));
fn = fieldnames(P);
for i=1:numel(fn)
    if isfield(estep_params, fn{i}) && ~isempty(estep_params.(fn{i})), P.(fn{i}) = estep_params.(fn{i}); end
end

% -------- Init outputs --------
estep_results = struct();
estep_results.posterior_source_covariances       = cell(F,1);
estep_results.transfer_functions                 = cell(F,1);
estep_results.effective_source_second_moments    = cell(F,1); % Psijj
estep_results.residual_transfer_functions        = cell(F,1);
estep_results.residual_covariances               = cell(F,1);
estep_results.posterior_noise_covariances        = cell(F,1); % L Σ_post L'
estep_results.effective_residual_second_moments  = cell(F,1); % Psixixi
estep_results.initial_precision_matrices         = cell(F,1); % inv(Σ_post)
estep_results.prior_precision_matrices           = cell(F,1); % inv(Σ_prior) 仅诊断
estep_results.whitening                          = struct('D',[],'Sigma_tilde',[]);
estep_results.noise                              = struct('sigma2_scalar_final',[]);
estep_results.stats                              = struct();

% 标量噪声：提取 σ²
sigma2_scalar = [];
if strcmpi(P.noise_model,'scalar')
    sigma2_scalar = real(mean(diag((Sigma_xixi+Sigma_xixi')/2)));
    if ~isfinite(sigma2_scalar) || sigma2_scalar<=0, sigma2_scalar = 1e-4; end
end

% -------- Inner E-step: loop over frequencies + optional scalar-noise update --------
for itE = 1:max(1,P.estep_inner_iters)
    sacc = 0;  % for scalar noise
    for f = 1:F
        % 1) Σ_post
        Sigma_post_f = module2_posterior_source_covariance(Sigma_prior{f}, L, Sigma_xixi, ...
                                struct('regularization_factor',P.eps_reg, ...
                                       'jitter_max_tries',6, ...
                                       'ensure_positive_definite',true, ...
                                       'min_eigenvalue_ratio',P.eig_floor_ratio, ...
                                       'verbose',false));
        estep_results.posterior_source_covariances{f} = Sigma_post_f;

        % 2) Tjv
        Tjv_f = module2_dstf_from_posterior(L, Sigma_post_f, Sigma_xixi); % n×p
        estep_results.transfer_functions{f} = Tjv_f;

        % 3) 有效二阶矩（Psijj / Psixixi）
        [Psijj_f, Txi_v_f, Sxixi_f, Sig_xi_post_f, Psixixi_f] = ...
            module2_effective_moments(Svv_cell{f}, L, Sigma_post_f, Tjv_f);
        estep_results.effective_source_second_moments{f}   = Psijj_f;
        estep_results.residual_transfer_functions{f}       = Txi_v_f;
        estep_results.residual_covariances{f}              = Sxixi_f;
        estep_results.posterior_noise_covariances{f}       = Sig_xi_post_f;
        estep_results.effective_residual_second_moments{f} = Psixixi_f;

        % 4) 初始化精度（后验）& 先验精度（诊断）
        Omega_init_f = inv_psd_robust(Sigma_post_f, P.eps_reg, P.eig_floor_ratio);
        Omega_prior_f  = inv_psd_robust(Sigma_prior{f}, P.eps_reg, P.eig_floor_ratio);
        estep_results.initial_precision_matrices{f} = Omega_init_f;
        estep_results.prior_precision_matrices{f}   = Omega_prior_f;

        % 5) 标量噪声累计
        if strcmpi(P.noise_model,'scalar')
            sacc = sacc + trace(Psixixi_f)/p;
        end
    end

    % 6) 标量噪声更新（与 BC-V 口径一致）+ 指数平滑
    if strcmpi(P.noise_model,'scalar')
        sigma2_new = sacc / F;
        sigma2_scalar = max((1-P.eta_sigma)*sigma2_scalar + P.eta_sigma*sigma2_new, 1e-10);
        Sigma_xixi = sigma2_scalar * eye(p);  % 覆盖下一轮 E 内循环
    else
        % 非标量模型在 adapter 或上层另行处理
    end
end
estep_results.noise.sigma2_scalar_final = sigma2_scalar;

% -------- Whitening（白化 Psijj，得到 Sjj_tilde 与 D_src） --------
if P.do_whitening
    [D_src_cell, Sjj_tilde_cell] = module2_whitening_from_covset( ...
        estep_results.effective_source_second_moments, P.whitening_opts);
    estep_results.whitening.D = D_src_cell;
    estep_results.whitening.Sigma_tilde = Sjj_tilde_cell;
end

% -------- 诊断/统计（简要） --------
try
    cn = zeros(F,1);
    for f=1:F
        cn(f) = cond(L*Sigma_prior{f}*L' + Sigma_xixi);
    end
    estep_results.stats.cond_forward = cn;
catch
end
end
function A = inv_psd_robust(S, eps_reg, min_ratio)
    S=0.5*(S+S'); [U,D]=eig(full(S),'vector'); d=real(D);
    dmax=max(d); floor_val=max(min_ratio*max(dmax,eps),0); d(d<floor_val)=floor_val;
    if eps_reg>0, d=(d+eps_reg*dmax)/(1+eps_reg); end
    A=U*diag(1./d)*U'; A=0.5*(A+A');
end