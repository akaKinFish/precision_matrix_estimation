%% PIPELINE_SOURCE_EM  (Source-domain EM pipeline, scalar estep_in version)
% Build estep_in as a *scalar struct* with per-frequency matrices stored as {F×1} cells.

%% 0) Simulation (Module 7)
n  = 10;  p  = 3;  F  = 3;  T  = 4096;

[Omega_true, Sigma_true, emp_covariance, sim] = module7_simulation_improved_complex( ...
    'n_nodes', n, 'n_sensors', p, 'n_freq', F, 'n_samples', T, ...
    'generate_leadfield', true, ...
    'leadfield_type', 'simple', ...
    'random_seed', 42);

L          = sim.leadfield_matrix;   % p×n
Sigma_xixi = sim.Sigma_xixi;         % p×p

%% 1) Build E-step input as a *scalar struct* (NOT a struct array)
emp_cov_cell = coerce_cov_cell(emp_covariance);     % supports cell / 3D / single
F            = numel(emp_cov_cell);
prior_cell   = repmat({eye(n)}, F, 1);

estep_in = struct();                                  % scalar struct
estep_in.leadfield_matrix         = L;                % p×n
estep_in.empirical_covariances    = emp_cov_cell;     % {F×1}, p×p
estep_in.source_prior_covariances = prior_cell;       % {F×1}, n×n
estep_in.noise_covariance         = Sigma_xixi;       % p×p
estep_in.frequencies              = (1:F);            % 1×F

assert(isstruct(estep_in) && isscalar(estep_in), 'estep_in must be a scalar struct');
assert(isa(estep_in.empirical_covariances,'cell') && numel(estep_in.empirical_covariances)==F);
assert(all(cellfun(@(A) isequal(size(A),[p p]), estep_in.empirical_covariances)));
assert(isa(estep_in.source_prior_covariances,'cell') && numel(estep_in.source_prior_covariances)==F);
assert(all(cellfun(@(S) isequal(size(S),[n n]), estep_in.source_prior_covariances)));

%% 2) E-step (Module 2 + wrapper)  —— 数值保障：开启 PSD & 统一加载
estep_out = module2_estep(estep_in, struct( ...
    'ensure_hermitian', true, ...
    'ensure_real_diag', true, ...
    'ensure_psd',       true, ...      % <--- enable PSD
    'psd_tol',          1e-10, ...
    'diag_loading',     1e-10 ...      % <--- small loading at one place
));
Sjj_hat = estep_out.source_second_moments;   % {F×1}, n×n

%% 3) Preprocessing / Whitening in source domain (Module 1 wrapper)
pre = module1_preproc_from_covset(Sjj_hat, struct( ...
    'smoothing_method','moving_average', ...
    'loading_factor',  1e-6, ...
    'min_power',       1e-10, ...   % <- fixed name here
    'verbose',         false ...
));
D_src     = pre.D;            % {F×1}, n×n
Sjj_tilde = pre.Sigma_tilde;  % {F×1}, n×n

%% 4) Active set (Module 3)
input_data_m3 = struct();
input_data_m3.whitened_covariances = Sjj_tilde;  % 16x1 cell
input_data_m3.frequencies          = 1:F;        % 1xF double

act = module3_active_set(input_data_m3, struct( ...
        'proxy_method',   'correlation', ...
        'quantile_level', 0.25, ...
        'force_diagonal_active', true, ...
        'verbose', false ...
    ));
A_masks = arrayfun(@(f) logical(act.combined_active_mask(:,:,f)), 1:F, 'uni', 0);

%% 5) Hyperparameters (Module 6)
K = make_frequency_kernel(F, 3.0);    % single source of truth
W = make_uniform_weight(n);

input_data_m6 = struct();
input_data_m6.whitened_covariances = Sjj_tilde;
input_data_m6.kernel_matrix        = K;
input_data_m6.weight_matrix        = W;
input_data_m6.active_set_mask      = {A_masks{:}};

hp = module6_hyperparameter_config(input_data_m6, struct('use_gershgorin', true));

lambda1 = hp.lambda1;  
lambda2_suggested = hp.lambda2_suggested;  
alpha   = hp.alpha;

% 默认采用建议值，可在下方覆盖并记录
lambda2_effective = lambda2_suggested;

%% 6) Proximal solver (Module 5) in whitened domain
Gamma0 = cellfun(@(S) diag(1 ./ max(real(diag(S)), 1e-8)), Sjj_tilde, 'uni', 0);

input_data_m5 = struct();
input_data_m5.whitened_covariances = Sjj_tilde;
input_data_m5.initial_precision    = Omega_true;
input_data_m5.smoothing_kernel     = K;
input_data_m5.weight_matrix        = W;
input_data_m5.active_set_mask      = {A_masks{:}};
input_data_m5.whitening_matrices   = D_src;      % <--- 为 Live 可视化白化 GT

% === 空间平滑配置（λ3 + 谱归一 L） ===
% 若你已有 make_source_graph 的 L，可直接谱归一：
%   [Lsp, Linfo] = normalize_graph_laplacian(L_raw, 'spectral');
% 这里示例直接从 Omega_true 的节点图来源处取得 L_raw，自行替换为你的 L_raw。
if isfield(sim, 'source_graph') && isfield(sim.source_graph, 'L')
    L_raw = sim.source_graph.L;
else
    % 若仿真没暴露 L，这里提供一个简单兜底（不会用于推理，只给示例）
    L_raw = eye(n) - (diag(sum(ones(n)-eye(n),2)) \ (ones(n)-eye(n))); % 非常粗糙，占位
end
[Lsp, Linfo] = normalize_graph_laplacian(L_raw, 'spectral');
fprintf('Spatial L normalized: eig[min,max]=[%.2e, %.2e] -> [%.2e, %.2e]\n', ...
    Linfo.min_eig_before, Linfo.max_eig_before, Linfo.min_eig_after, Linfo.max_eig_after);

% 启用/关闭 λ3（默认建议：0.3~1.0 * λ1）
lambda3 = 0.3 * lambda1;

% Hadamard 语义 + 实时诊断
params5 = struct( ...
    'lambda1', 1e-3, ...
    'lambda2', 0.02, ...            % 生效值
    'lambda2_suggested', lambda2_suggested, ...  % 记录来源
    'alpha0',  5e-3, ...
    'max_iter', 100, 'verbose', true, ...
    'active_set_update_freq', 10, ...
    'alpha_max', 2e-3, 'alpha_up', 1.05, ...
    'alpha_down', 0.7, 'alpha_grow_patience', 2, ...
    'obj_improve_tol', 1e-6, ...
    'weight_mode','hadamard', 'use_graph_laplacian', true, ...
    'diag', struct( 'enable', true, 'update_every', 1, 'metrics_every', 1, ...
                    'print_every', 1, 'f_view', 1, 'log_csv', 'prox_trace.csv', ...
                    'keep_fig_open', true, 'weight_mode','hadamard', ...
                    'use_graph_laplacian', true ) );

% 在线可视化（Live）
params5.diag.live_plot = struct( ...
    'enable', true, ...
    'f_view', 1, ...
    'plot_every', 5, ...
    'value_mode', 'abs', ...
    'ground_truth_domain', 'source' ...
);

    params5.diag.live_plot.ground_truth_precision = Omega_true; % 若 GT 在源域

% 单频空间平滑开关
if lambda3 > 0
    params5.lambda3                     = lambda3;
    params5.spatial_graph_matrix        = Lsp;     % 谱归一后 L
    params5.spatial_graph_is_laplacian  = true;
    if ~isfield(params5,'spatial_weight_mode') || isempty(params5.spatial_weight_mode)
        params5.spatial_weight_mode = 'node';      % 'node' or 'hadamard'
    end
end

[Gamma_tilde_star, prox_res] = module5_proximal(input_data_m5, params5);

%% 7) Recolor to source scale (Module 8)
input_data_m8 = struct();
input_data_m8.whitened_precision_matrices = Gamma_tilde_star;
input_data_m8.whitening_matrices          = D_src;
recol = module8_recoloring(input_data_m8, struct());
Omega_src = recol.recolored_precision_matrices;

%% 8) Simple readout (optional)
f_view = 1;
Om     = Omega_src{f_view};
pcorr  = abs(-Om) ./ sqrt((abs(diag(Om))+eps) * (abs(diag(Om))+eps)');
pcorr(1:n+1:end) = 0;
fprintf('Done. Example partial coherence at f=%d: max=%g, median=%g\n', ...
    f_view, max(pcorr(:)), median(pcorr(pcorr>0)));

viz_pipeline_summary(prox_res, Gamma_tilde_star, Sjj_tilde, K, A_masks, Omega_src, Omega_true);
opt_compare = struct('mode','match_sparsity', 'f_view', 1, 'show_pr', true, 'show_roc', true);
opt_compare.Gamma = Gamma_tilde_star;
metrics = gt_compare_and_plot(Omega_true, Omega_src, opt_compare);
Sigma_hat_test   = ensure_sigma_collection(Sjj_hat);
Sigma_tildeteest = ensure_sigma_collection(Sjj_tilde);
debug_em_convergence(Sigma_hat_test,   'weight_mode','hadamard', 'use_laplacian',true);
debug_em_convergence(Sigma_tildeteest, 'weight_mode','hadamard', 'use_laplacian',true);

%% ===== Local helpers =====
function C = coerce_cov_cell(X, F_hint)
% Turn X into {F×1} cell of square matrices.
    if isa(X,'cell'), C = X(:); return; end
    if isnumeric(X) && ndims(X)==3 && size(X,1)==size(X,2)
        F = size(X,3); C = cell(F,1);
        for f = 1:F, C{f} = X(:,:,f); end, return;
    end
    if isnumeric(X) && ismatrix(X) && size(X,1)==size(X,2)
        if nargin>=2 && ~isempty(F_hint), C = repmat({X}, F_hint, 1);
        else, C = {X}; end, return;
    end
    error('coerce_cov_cell:unsupported','Expect cell{F,1} | p×p×F | single p×p.');
end

function K = make_frequency_kernel(F, sigma)
% Gaussian kernel along frequency index (1..F), row-normalized
    if nargin < 2, sigma = 3.0; end
    [I,J] = ndgrid(1:F,1:F);
    K = exp(-((I-J).^2)/(2*sigma^2));
    K = (K + K')/2;              % 保证对称
    K = K / max(sum(K,2));       % 统一尺度（不会破坏对称）
end

function W = make_uniform_weight(n)
% Ones off-diagonal, zeros on diagonal
    W = ones(n); W(1:n+1:end) = 0;
end
