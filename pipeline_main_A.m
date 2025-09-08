function outs = pipeline_main_A(opts)
% PIPELINE A (SOURCE-DOMAIN end-to-end)
% - n = 128 sources, p = 18 sensors
% - Module7/2/1/3/6/5/8 remain unchanged; this file is glue only.

if nargin<1, opts = struct(); end

%% ================== Module 7: Simulation ==================
sim_par = struct();
sim_par.n_nodes   = 128;         % n (sources)
sim_par.n_freq    = 16;          % F
sim_par.n_samples = 4096;        % T
sim_par.generate_leadfield = true;
sim_par.n_sensors = 18;          % p (sensors)
sim_par.leadfield_type = 'simple';
sim_par.electrode_layout = '1020-19';
sim_par.source_space = 'cortical';
sim_par.random_seed = 42;

fprintf('[Pipeline] Generating simulation via Module 7 (leadfield=true)...\n');
[Omega_true, Sigma_true, Sigma_emp, sim_params] = module7_simulation_improved_complex( ...
   'n_nodes',sim_par.n_nodes, 'n_freq',sim_par.n_freq, 'n_samples',sim_par.n_samples, ...
   'generate_leadfield',true, 'n_sensors',sim_par.n_sensors, ...
   'leadfield_type','simple','electrode_layout',sim_par.electrode_layout, ...
   'source_space',sim_par.source_space,'random_seed',sim_par.random_seed);

F  = sim_params.n_freq;
p  = sim_params.n_sensors;
n  = sim_par.n_nodes;

%% ================== Module 2: E-step（源域统计） ==================
% 重要：此处使用 module2_estep_source_stats 来得到源域 Ŝ_jj 与 Ω_init（n×n）。
fprintf('[Pipeline] E-step (Module 2, source stats)...\n');

estep_in = struct();
estep_in.leadfield_matrix         = sim_params.leadfield_matrix;   % L (p×n)
estep_in.empirical_covariances    = sim_params.Sigma_vv_observed;  % {F} 传感器域 Σ_vv(obs)
estep_in.source_prior_covariances = Sigma_true;                    % {F} 源先验 Σ_jj,true (演示用)
estep_in.noise_covariance         = sim_params.Sigma_xixi;         % Σ_ξξ (p×p)
estep_in.frequencies              = (1:F)';

estep_params = struct('regularization_factor',1e-8, ...
                      'condition_threshold',1e12, ...
                      'min_eigenvalue_ratio',1e-12, ...
                      'verbose', true);

% 该包装器内部调用你的 module2_estep / module2_estep_main
ss = module2_estep_source_stats(estep_in, estep_params);

% 源域二阶矩（M-step 的“样本协方差”）：{F}，每个 n×n
Sjj_hat        = ss.source_second_moments;
% 源域初始精度：{F}，每个 n×n（确保你的 module2_estep_main 已改为源域输出）
Gamma_init_cell = ss.estep_results.initial_precision_matrices;

%% ================== Module 1: Whitening（源域） ==================
fprintf('[Pipeline] Whitening in SOURCE domain via Module 1...\n');

preproc_params = struct('smoothing_method','moving_average', 'window_size',5, ...
                        'diagonal_loading',true, 'loading_factor',1e-3, ...
                        'min_power',1e-10, 'target_diagonal',1.0, ...
                        'diagonal_tolerance',0.1, 'force_hermitian',true, ...
                        'check_psd',true, 'verbose',true);

% 用你提供的兼容包装：把 {F} 的源域协方差直接送进 Module1 主流程
pre = module1_preproc_from_covset(Sjj_hat, preproc_params);

% 源域 whiten 后协方差 {F}（n×n）与对角 whitening 矩阵 {F}（n×n）
Sigma_tilde = pre.Sigma_tilde;
D_list      = pre.D;

%% ================== Module 3: Active Set（源域） ==================
fprintf('[Pipeline] Active set selection (Module 3, source domain)...\n');

input3 = struct();
input3.whitened_covariances       = Sigma_tilde;       % {F} n×n
input3.initial_precision_matrices = Gamma_init_cell;   % {F} n×n（用于 precision 代理）
input3.frequencies                = (1:F)';

params3 = struct('proxy_method','precision', ...       % 或 'correlation'
                 'quantile_level',0.25, ...
                 'force_diagonal_active',true, ...
                 'verbose', true);

as_res = module3_active_set_main(input3, params3);

% 聚合为单一掩码（逐频 OR），并去掉对角
active_mask = any(as_res.combined_active_mask, 3);   % n×n
active_mask = active_mask | active_mask';
active_mask(1:n+1:end) = false;

%% ================== Module 6: Hyperparameters（源域） ==================
fprintf('[Pipeline] Hyperparameters (Module 6, source domain)...\n');

% 构造一个简单 RBF 频核（避免依赖 pdist/squareform）
freq = (1:F)'; sig = max(1, round(F/6));
[I,J] = ndgrid(freq,freq);
K = exp(-((I-J).^2) / (2*sig^2));       % F×F

Wg = ones(n) - eye(n);                  % n×n，仅惩罚非对角

input6 = struct();
input6.whitened_covariances = Sigma_tilde;           % {F} n×n
input6.kernel_matrix        = K;                     % F×F
input6.weight_matrix        = Wg;                    % n×n
input6.active_set_mask      = repmat({active_mask}, F, 1);  % {F} n×n

cfg6 = struct('safety_margin',0.90,'verbose',true,'use_gershgorin',true, ...
              'initialization_method','inverse');

hp = module6_hyperparameters(input6, cfg6);

lambda1 = hp.lambda1;
alpha0  = hp.alpha;
lambda2 = max(hp.lambda2_suggested, 5e-3);   % 给个地板

%% ================== Module 5: Proximal Solver（源域） ==================
fprintf('[Pipeline] Proximal solver (Module 5, source domain)...\n');

params5 = struct();
params5.lambda1 = lambda1;
params5.lambda2 = lambda2;
params5.alpha0  = min(alpha0, 5e-3);
params5.backtracking       = true;
params5.bt_shrink          = 0.5;
params5.bt_max_iter        = 20;
params5.max_new_edges_per_sweep = 64;
params5.max_iters          = 200;
params5.tol                = 1e-5;
params5.verbose            = true;

% 这里的 Sigma_tilde / active_mask 均为 源域（n×n）
input5 = struct('Sigma_tilde',Sigma_tilde,'active_mask',active_mask);

[Gamma_tilde_opt, proximal_results] = module5_proximal_main(input5, params5);
% Gamma_tilde_opt: {F}，每个 n×n（whitened precision）

%% ================== Module 8: Recoloring（源域） ==================
fprintf('[Pipeline] Recoloring (Module 8) back to SOURCE-domain scale...\n');

input8 = struct();
input8.whitened_precision_matrices = Gamma_tilde_opt;    % {F} n×n
input8.whitening_matrices          = D_list;             % {F} n×n
input8.original_covariances        = Sjj_hat;            % {F} n×n（回到 Ŝ 的尺度）
input8.active_set_masks            = repmat({active_mask}, F,1);

rec8 = struct('force_hermitian',true,'validate_spd',false, ...
              'compute_quality_metrics',true,'verbose',true);

recoloring_results = module8_recoloring_main(input8, rec8);

%% ================== Reporting（源域一致性） ==================
inv_err = zeros(F,1);
for f=1:F
    Om = recoloring_results.recolored_precision_matrices{f};  % n×n
    S  = Sjj_hat{f};                                          % n×n
    inv_err(f) = norm(Om*S - eye(n), 'fro') / sqrt(n);
end
avg_inv_err = mean(inv_err);
fprintf('\n[Pipeline] Avg inverse-consistency error (SOURCE domain): %.3e\n', avg_inv_err);
fprintf('[Pipeline] Done.\n');

%% ================== Outputs ==================
outs = struct();
outs.sim_params         = sim_params;
outs.Sigma_true         = Sigma_true;                % {F} n×n
outs.Sigma_emp          = Sigma_emp;                 % {F} (sensor-domain, from sim)
outs.estep_results      = ss.estep_results;          % 包含 T_jv / Σ_jj,post / Ω_init(源域)
outs.Sjj_hat            = Sjj_hat;                   % {F} n×n
outs.preprocessing      = pre;                       % Module1（源域）
outs.active_results     = as_res;                    % Module3（源域）
outs.hyperparams        = hp;                        % Module6
outs.proximal_results   = proximal_results;          % Module5
outs.Gamma_tilde_opt    = Gamma_tilde_opt;           % {F} n×n
outs.recoloring_results = recoloring_results;        % Module8（源域）
outs.avg_inv_err        = avg_inv_err;
end
