function [Omega_est, Dsrc_est, Gamma_est, outs] = my_test2_adapter(emp_cov_cell, L, T, cfg, Omega_true_opt)
% 薄适配器：预处理→调用 module2_estep_main（含 BC-V 风格 E 步与噪声内循环、白化）→ M 步（5/8）→ 可选 refit

if nargin < 3, error('至少需要 emp_cov_cell, L, T'); end
if nargin < 4 || isempty(cfg), cfg = struct(); end
if nargin < 5, Omega_true_opt = []; end

% ---------- 基本准备 ----------
emp_cov_cell = coerce_cov_cell(emp_cov_cell);
F = numel(emp_cov_cell);
[p,n] = size(L);

verbose       = getf(cfg,'verbose',true);
do_scale_L    = getf(cfg,'do_scale_L',true);
do_scale_data = getf(cfg,'do_scale_data',true);
q_act         = getf(cfg,'q_act',0.10);
kernel_sigma  = getf(cfg,'kernel_sigma',3.0);
lambda3_ratio = getf(cfg,'lambda3_ratio',0.0);
lambda1_override = getf(cfg,'lambda1',[]);
lambda2_factor   = getf(cfg,'lambda2_factor',1.0);
do_support_refit = getf(cfg,'do_support_refit',true);
noise_model   = getf(cfg,'noise_model','scalar');

% ---------- 预处理（“数据换壳”，非计算性） ----------
if do_scale_L
    sL = sqrt(trace(L*L')/size(L,1));
    if isfinite(sL)&&sL>0, L = L/sL; if verbose, fprintf('[scale] L by 1/%.3g\n',sL); end, end
end
if do_scale_data
    gamma_grid = logspace(-4,1,30);
    warm_tmp   = eloreta_warmstart_from_covs(emp_cov_cell, L, gamma_grid, struct('maxit',50,'tol',1e-6,'verbose',false));
    scl = ones(F,1);
    for f=1:F
        d = mean(real(diag(warm_tmp.Sjj_e{f}))); if ~isfinite(d)||d<=0, d=1; end
        emp_cov_cell{f} = emp_cov_cell{f}/d; scl(f)=d;
    end
    if verbose, fprintf('[scale] data per-freq by eLORETA mean power, median=%.3g\n', median(scl)); end
end

% ---------- E-step（把重计算交给 module2） ----------
estep_in = struct();
estep_in.leadfield_matrix         = L;
estep_in.empirical_covariances    = emp_cov_cell;
estep_in.source_prior_covariances = repmat({eye(n)}, F, 1);
% 初噪声：标量模型用 trace/p * 1e-2
if strcmpi(noise_model,'scalar')
    trS = 0; for f=1:F, trS = trS + trace(emp_cov_cell{f}); end
    sigma2_init = real(trS/(F*p))*1e-2; sigma2_init = max(sigma2_init,1e-10);
    estep_in.noise_covariance = sigma2_init * eye(p);
else
    estep_in.noise_covariance = eye(p)*1e-4; % 其它模型占位
end
estep_in.frequencies = 1:F;

estep_params = struct('noise_model',noise_model, ...
                      'estep_inner_iters', getf(getf(cfg,'em',struct()), 'max_estep_iter', 2), ...
                      'eta_sigma', getf(getf(cfg,'em',struct()), 'eta_sigma', 0.3), ...
                      'do_whitening', true, ...
                      'whitening_opts', struct('smoothing_method','moving_average','loading_factor',1e-6,'min_power',1e-10,'verbose',false));

E = module2_estep_main(estep_in, estep_params);

Sjj_tilde = E.whitening.Sigma_tilde;
D_src     = E.whitening.D;
Psijj     = E.effective_source_second_moments;
Sigma_post_cell = E.posterior_source_covariances;
sigma2_final = E.noise.sigma2_scalar_final; %#ok<NASGU>

% ---------- Active set ----------
input_data_m3 = struct();
        input_data_m3.whitened_covariances = Sjj_tilde;
        input_data_m3.frequencies = 1:F;
act = module3_active_set(input_data_m3, struct('proxy_method','correlation', ...
                     'quantile_level', q_act, 'force_diagonal_active', true, 'verbose', false));
A_masks = arrayfun(@(f) logical(act.combined_active_mask(:,:,f)), 1:F, 'uni', 0);

% ---------- 超参数（Module 6） ----------
K = make_frequency_kernel(F, kernel_sigma);
K = real(0.5*(K+K')); K = K / max(1,max(sum(K,2)));
W = ones(n); W(1:n+1:end)=0;
input_data_m6 = struct();
        input_data_m6.whitened_covariances = Sjj_tilde;
        input_data_m6.kernel_matrix        = K;
        input_data_m6.weight_matrix        = W;
        input_data_m6.active_set_mask      = {A_masks{:}};
hp = module6_hyperparameter_config(input_data_m6, struct('use_gershgorin', true));
lambda1 = hp.lambda1; if ~isempty(lambda1_override), lambda1 = lambda1_override; end
lambda2 = lambda2_factor * hp.lambda2_suggested;
lambda3 = lambda3_ratio * lambda1;
alpha0  = hp.alpha;

if verbose
    fprintf('[HP] λ1=%.3g, λ2=%.3g, λ3=%.3g, α=%.3g\n', lambda1, lambda2, lambda3, alpha0);
end

% ---------- M-step（Module 5 PGD） ----------
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
params5.spatial_graph_matrix = eye(n); % 若有图替换
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
input_data_m5.active_set_mask      = {A_masks{:}}; %#ok<CCAT>
input_data_m5.whitening_matrices   = D_src;

[Gamma_tilde, prox_res] = module5_proximal(input_data_m5, params5);

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
        G0 = Gamma_tilde{f}; M = abs(G0)>1e-8; M(1:n+1:end)=true; support_cell{f}=M;
    end
    params_refit = params5; params_refit.max_iter=10; params_refit.lambda2=1e-8;
    params_refit.l1_weights=[]; params_refit.active_set_update_freq=Inf; params_refit.use_single_step=true;
    Ginit_ref = cell(F,1); for f=1:F, Ginit_ref{f} = 0.5*( (Gamma_tilde{f}.*support_cell{f}) + (Gamma_tilde{f}.*support_cell{f})'); end
    input_data_m5_ref = input_data_m5; input_data_m5_ref.initial_precision=Ginit_ref; input_data_m5_ref.active_set_mask=support_cell;
    [Gamma_tilde_refit, ~] = module5_proximal(input_data_m5_ref, params_refit);
    input_data_m8_refit = struct();
    input_data_m8_refit.whitened_precision_matrices = Gamma_tilde_refit;
    input_data_m8_refit.whitening_matrices = D_src;
    recol2 = module8_recoloring(input_data_m8_refit, struct());
    Omega_src = recol2.recolored_precision_matrices;
end

% ---------- 输出 ----------
Gamma_est = Gamma_tilde;
Omega_est = Omega_src;
Dsrc_est  = invert_cell_spd_(Omega_src, 1e-8, 1e-12);

outs = struct();
outs.estep = E;
outs.prox_res = prox_res;
outs.A_masks = A_masks;
outs.K = K; outs.W = W;

% 小结
if verbose
    f_view = 1; Om = Omega_src{f_view};
    pcorr = abs(-Om) ./ sqrt((abs(diag(Om))+eps) * (abs(diag(Om))+eps)');
    pcorr(1:n+1:end)=0;
    fprintf('Done. Partial coherence@f=%d: max=%g, median=%g\n', f_view, max(pcorr(:)), median(pcorr(pcorr>0)));
end
end

% ===== helpers（与你之前版本一致）=====
%% ======== 辅助函数（保持你原有版本，仅微小稳健化） ========
function v = get_field(s, name, default_val)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default_val;
end
end

function C = coerce_cov_cell(X, F_hint)
if isa(X,'cell')
    C = X(:);
    return;
end
if isnumeric(X) && ndims(X)==3 && size(X,1)==size(X,2)
    F = size(X,3);
    C = cell(F,1);
    for f = 1:F
        C{f} = X(:,:,f);
    end
    return;
end
if isnumeric(X) && ismatrix(X) && size(X,1)==size(X,2)
    if nargin>=2 && ~isempty(F_hint)
        C = repmat({X}, F_hint, 1);
    else
        C = {X};
    end
    return;
end
error('coerce_cov_cell:unsupported','Expect cell{F,1} | p×p×F | single p×p.');
end

function K = make_frequency_kernel(F, sigma)
if nargin < 2, sigma = 3.0; end
[I,J] = ndgrid(1:F,1:F);
K = exp(-((I-J).^2)/(2*sigma^2));
K = (K + K')/2;
end

function W = make_uniform_weight(n)
W = ones(n);
W(1:n+1:end) = 0;
end

function L = laplacian_placeholder(n)
A = ones(n) - eye(n);
d = sum(A,2);
L = diag(d) - A;
end

function [Lnorm, info] = normalize_graph_laplacian(L, mode)
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

function Gamma_init = transport_init(Omega_prev, D_src, Sjj_tilde)
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
    Gamma_init{f} = (D \ Omega_prev{f}) / D;  % 等价于 inv(D)*Ω*inv(D)'
    Gamma_init{f} = (Gamma_init{f} + Gamma_init{f}')/2;
end
end

function Sigma_prior = invert_and_fix(Omega_cell, eps_ld)
if nargin < 2, eps_ld = 1e-10; end
F = numel(Omega_cell);
n = size(Omega_cell{1},1);
Sigma_prior = cell(F,1);

for f=1:F
    Om = (Omega_cell{f} + Omega_cell{f}')/2;
    Om(~isfinite(Om)) = 0;
    d = real(diag(Om));
    d = max(d, eps_ld);
    Om(1:n+1:end) = d;
    [U, S] = eig(full(Om), 'vector');
    S = real(S);
    S = max(S, 2*eps_ld);
    Sigma = U * diag(1./S) * U';
    Sigma = (Sigma + Sigma')/2;
    Sigma(1:n+1:end) = real(diag(Sigma));
    Sigma_prior{f} = Sigma;
end
end

function [dOmega, dS] = compute_deltas(Omega, Omega_prev, Sjj, Sjj_prev)
if isempty(Omega_prev)
    dOmega = inf;
else
    num=0; den=0;
    for f=1:numel(Omega)
        num = num + norm(Omega{f}-Omega_prev{f},'fro');
        den = den + norm(Omega_prev{f},'fro');
    end
    dOmega = num / max(1, den);
end

if isempty(Sjj_prev)
    dS = inf;
else
    num=0; den=0;
    for f=1:numel(Sjj)
        num = num + norm(Sjj{f}-Sjj_prev{f},'fro');
        den = den + norm(Sjj_prev{f},'fro');
    end
    dS = num / max(1, den);
end
end

function g = read_grad_norm_(prox_res)
g = NaN;
try
    if isfield(prox_res,'grad_norm')
        g = prox_res.grad_norm;
    end
    if isnan(g) && isfield(prox_res,'gradient_norm_history')
        h = prox_res.gradient_norm_history;
        if ~isempty(h)
            g = h(end);
        end
    end
    if isnan(g) && isfield(prox_res,'grad_norm_mean')
        g = prox_res.grad_norm_mean;
    end
    if isnan(g)
        g = -1;
    end
catch
    g = -1;
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

function [L_byS, L_byG] = estimate_L_candidates(Gamma_init, Sjj_tilde)
F = numel(Sjj_tilde);
Ls = 0;
Lg = 0;

for f = 1:F
    St = (Sjj_tilde{f} + Sjj_tilde{f}')/2;
    s = svds(St, 1);
    Ls = max(Ls, s^2);

    G = (Gamma_init{f} + Gamma_init{f}')/2;
    try
        lam_min = min(real(eig(G)));
        if ~isfinite(lam_min) || lam_min <= 0
            invnorm = svds(pinv(G), 1);
        else
            invnorm = 1/lam_min;
        end
    catch
        invnorm = svds(pinv(G), 1);
    end
    Lg = max(Lg, invnorm^2);
end

L_byS = Ls;
L_byG = Lg;
end

function diag_whitening_sanity(Sjj_hat, Sjj_tilde, D_src, iter_id)
F = numel(Sjj_hat);
md_hat = zeros(F,1); sd_hat = zeros(F,1);
md_til = zeros(F,1); sd_til = zeros(F,1);
rel_err = zeros(F,1);
is_spd = true;

for f = 1:F
    Sh = (Sjj_hat{f} + Sjj_hat{f}')/2;
    St = (Sjj_tilde{f} + Sjj_tilde{f}')/2;
    D = D_src{f};
    md_hat(f) = mean(real(diag(Sh)));
    sd_hat(f) = std(real(diag(Sh)));
    md_til(f) = mean(real(diag(St)));
    sd_til(f) = std(real(diag(St)));
    R = D*Sh*D - St;
    rel_err(f) = norm(R,'fro') / max(1, norm(St,'fro'));
    try
        ev = eig(St);
        is_spd = is_spd && all(real(ev) > -1e-10);
    catch
        is_spd = false;
    end
end

fprintf(['[DIAG][t=%d] <Whitening>\n    mean(diag S_hat)=%.3g±%.3g | ' ...
    'mean(diag S_tilde)=%.3g±%.3g (≈1)\n' ...
    '    max relErr ||D*S*D - S_tilde||_F / ||S_tilde||_F = %.2e | ' ...
    'S_tilde SPD? %s\n'], ...
    iter_id, mean(md_hat), mean(sd_hat), mean(md_til), mean(sd_til), ...
    max(rel_err), ternary(is_spd,'YES','NO'));

if mean(md_til) > 2 || mean(md_til) < 0.5
    warning('[DIAG] whitened diagonals far from 1 (mean=%.3g).', mean(md_til));
end
if max(rel_err) > 1e-6
    warning('[DIAG] D*S*D and S_tilde mismatch (%.2e).', max(rel_err));
end
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

%% ========= eLORETA warmstart helpers =========
function warm = eloreta_warmstart_from_covs(Svv_cell, L, gamma_grid, opts)
if nargin < 4, opts = struct(); end
maxit = getf(opts,'maxit',50);
tol = getf(opts,'tol',1e-6);
verb = getf(opts,'verbose',false);

F = numel(Svv_cell);
Tjv_cell = cell(F,1);
Sjj_cell = cell(F,1);
Om0_cell = cell(F,1);
gcv_cell = cell(F,1);
gopt_cell = cell(F,1);

for f = 1:F
    Svv = Svv_cell{f};
    [Tjv, Sjj, ~, gamma_opt, gcv] = eloreta_simple_(Svv, L, gamma_grid, maxit, tol, verb);
    Sjj = psd_project_(Sjj);
    Om0 = inv_psd_robust_(Sjj, 1e-8, 1e-12);
    Tjv_cell{f} = Tjv;
    Sjj_cell{f} = Sjj;
    Om0_cell{f} = Om0;
    gcv_cell{f} = gcv;
    gopt_cell{f} = gamma_opt;
end

warm = struct('Omega_init',{Om0_cell}, 'Sjj_e',{Sjj_cell}, 'Tjv',{Tjv_cell}, ...
    'gamma_opt',{gopt_cell}, 'gcv_curve',{gcv_cell});
end

function [Tjv, Sjj, W, gamma_opt, gcv] = eloreta_simple_(Svv, L, gamma_grid, maxit, tol, verb)
[p,n] = size(L);
gcv = zeros(numel(gamma_grid),1);
best.T = [];
best.W = [];
best.gamma = NaN;
best.score = Inf;

for k = 1:numel(gamma_grid)
    gamma = gamma_grid(k);
    w = ones(n,1);
    for it=1:maxit
        Winv = diag(1./w);
        A = hermi_(L*Winv*L');
        alpha = gamma * trace(A)/p;
        M = inv_psd_(A + alpha*eye(p));
        w_old = w;
        for i=1:n
            li = L(:,i);
            mii = real(li' * M * li);
            w(i) = sqrt(max(mii, eps));
        end
        if norm(w-w_old)/max(1,norm(w_old)) < tol
            break;
        end
    end
    Winv = diag(1./w);
    T = Winv * L' * M;
    Txiv = eye(p) - L*T;
    num = real(trace(hermi_(Txiv*Svv*Txiv')))/p;
    den = (real(trace(Txiv))/p)^2 + eps;
    gcv(k) = num / den;
    if gcv(k) < best.score
        best.T = T;
        best.W = diag(w);
        best.gamma = gamma;
        best.score = gcv(k);
    end
    if verb && (mod(k,10)==1)
        fprintf('[eLORETA] gamma=%.3g, GCV=%.3g\n', gamma, gcv(k));
    end
end

Tjv = best.T;
gamma_opt = best.gamma;
W = best.W;
Sjj = Tjv * Svv * Tjv';
end

function A = hermi_(A)
A = (A + A')/2;
end

function X = inv_psd_(X)
X = hermi_(X);
[U,S] = eig(X,'vector');
S = max(real(S), eps);
X = U*diag(1./S)*U';
X = hermi_(X);
end

function S = psd_project_(S)
S = hermi_(S);
[U,d] = eig(S,'vector');
d = max(real(d), 0);
S = U*diag(d)*U';
S = hermi_(S);
end

function v = getf(s, f, d)
if isfield(s,f) && ~isempty(s.(f))
    v = s.(f);
else
    v = d;
end
end
function S = invert_cell_spd_(Omega_cell, eps_reg, min_ratio)
    F=numel(Omega_cell); S=cell(F,1);
    for f=1:F, S{f}=inv_psd_robust_(Omega_cell{f}, eps_reg, min_ratio); end
end
