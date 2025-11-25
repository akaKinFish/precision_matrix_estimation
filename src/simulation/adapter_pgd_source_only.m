function out = adapter_pgd_source_only(args, opts)
% ADAPTER_PGD_SOURCE_ONLY — 只跑“源域（白化域）PGD”的轻量适配器
% 目的：
%   - 从 *源层* 的白化协方差出发（Sjj_tilde），
%   - 用 Module 6 自动估计超参数 (λ₁, α, λ₂_suggested)，
%   - 调 Module 5（proximal PGD）得到 {Γ̃_ω}，
%   - 可选：用 Module 8 进行 recolor 得到源域 {Ω_ω}。
%
% 命名/字段严格沿用现有 all_code 口径：
%   input_data_m6: whitened_covariances, kernel_matrix, weight_matrix, active_set_mask
%   input_data_m5: whitened_covariances, initial_precision, smoothing_kernel, weight_matrix,
%                  active_set_mask, whitening_matrices
%   params5      : lambda1, lambda2, lambda2_suggested, alpha0, max_iter, verbose, ...
%   调用         : module6_hyperparameter_config → module5_proximal → (optional) module8_recoloring
%
% ------------------
% Inputs
% ------------------
% args.Sjj_tilde      : cell{F,1} or (p×p×F)  — 源域白化协方差 Σ̃_ω
% args.K              : (F×F)                 — 频率平滑核
% args.W              : (p×p)                 — 权重矩阵（Graph-Laplacian等使用）
% args.A_masks        : cell{F,1} (optional)  — 活跃边集（可留空[]）
% args.D_src          : cell{F,1} (optional)  — 白化矩阵 D(ω)（用于可视化/可选recolor）
% args.Gamma_init     : cell{F,1} (optional)  — PGD 初值 Γ̃^0_ω（缺省用 I）
%
% opts.inner          : struct (optional)     — 与 test2 相同风格的 PGD 参数子集
%                        字段：max_iter, alpha_max, alpha_up, alpha_down,
%                              alpha_grow_patience, obj_improve_tol,
%                              weight_mode, use_graph_laplacian, diag
%                      若缺省，将设一组保守默认。
%
% opts.recolor        : logical (default true) — 是否调用 Module 8 做 recolor 到 {Ω_ω}
% opts.hp_config      : struct  (optional)     — 透传给 module6_hyperparameter_config 的配置
%
% %%% Jankova 对齐（可选，同 BC-V 口径的 λ₂ 规则）
% opts.jankova.override : logical (default false)
% opts.jankova.c        : double  (default 1.0)    — 缩放系数 c
% opts.jankova.q        : integer (default p)      — 变量数（或有效特征数）
% opts.jankova.m        : integer (default 1)      — 样本段数 m
%
% ------------------
% Outputs (struct)
% ------------------
% out.Gamma_tilde_star     : cell{F,1}  — PGD 输出的白化精度矩阵 {Γ̃_ω}
% out.Omega_src            : cell{F,1}  — （可选）recolor 后的源域精度矩阵 {Ω_ω}
% out.hp                   : struct      — Module 6 超参数输出（含 lambda1, alpha, lambda2_suggested）
% out.params5              : struct      — 实际喂给 Module 5 的参数（含 override 标记）
% out.proximal_results     : struct      — Module 5 返回的详细日志/统计
% out.recoloring_results   : struct      — Module 8 返回（若 opts.recolor=true 且提供 D_src）
%
% 备注：本适配器**不**做 E-step/传感域到源域的反演；它假设你已经在源域拿到了 Σ̃_ω。

if nargin < 1, error('adapter_pgd_source_only:missing_args','args required'); end
if ~isstruct(args), error('adapter_pgd_source_only:invalid_args','args must be struct'); end
if nargin < 2 || ~isstruct(opts), opts = struct(); end

% ---------- 规范化/基本检查 ----------
[Sjj_tilde, F, p] = ensure_cell_covariances(args);
K = must_field(args,'K','(F×F) kernel matrix');
W = must_field(args,'W','(p×p) weight matrix');

A_masks = [];
if isfield(args,'A_masks') && ~isempty(args.A_masks)
    A_masks = normalize_masks(args.A_masks, F, p);
end

D_src = [];
if isfield(args,'D_src') && ~isempty(args.D_src)
    D_src = ensure_cell_square(args.D_src, F, p, 'D_src');
end

Gamma_init = [];
if isfield(args,'Gamma_init') && ~isempty(args.Gamma_init)
    Gamma_init = ensure_cell_square(args.Gamma_init, F, p, 'Gamma_init');
else
    % 缺省用 I 做 warm start（与现有管线兼容）
    I = eye(p);
    Gamma_init = repmat({I}, F, 1);
end

% ---------- (1) 调 Module 6 拿超参数 ----------
input_data_m6 = struct();
input_data_m6.whitened_covariances = Sjj_tilde;
input_data_m6.kernel_matrix        = K;
input_data_m6.weight_matrix        = W;
input_data_m6.active_set_mask      = A_masks;  % 允许 []

hp_cfg = struct();
if isfield(opts,'hp_config') && isstruct(opts.hp_config)
    hp_cfg = opts.hp_config;
end
hp = module6_hyperparameter_config(input_data_m6, hp_cfg);

lambda1            = hp.lambda1;
lambda2_effective  = hp.lambda2_suggested;   % 默认使用建议值
alpha0             = hp.alpha;

% ---------- (1b) 可选：Jankova 对齐覆盖 λ₂ ----------
jankova = struct('override',false,'c',1.0,'q',p,'m',1);
if isfield(opts,'jankova') && isstruct(opts.jankova)
    fn = fieldnames(opts.jankova);
    for i=1:numel(fn), jankova.(fn{i}) = opts.jankova.(fn{i}); end
end
if jankova.override
    lambda2_effective = jankova.c * sqrt(log(max(jankova.q,1)) / max(jankova.m,1));
end

% ---------- (2) 组装 Module 5 输入/参数 ----------
input_data_m5 = struct();
input_data_m5.whitened_covariances = Sjj_tilde;
input_data_m5.initial_precision    = Gamma_init;
input_data_m5.smoothing_kernel     = K;
input_data_m5.weight_matrix        = W;
input_data_m5.active_set_mask      = A_masks;
input_data_m5.whitening_matrices   = D_src;   % 仅供可视化/记录，可为空

% inner 参数：保持 test2 风格字段名；缺省给一套保守值
inner = default_inner(p);
if isfield(opts,'inner') && isstruct(opts.inner)
    inner = set_defaults(opts.inner, inner);
end

% alpha 的安全夹取（与 main 的 alpha_min/max 配合）
alpha0_used = min(max(alpha0, inner.alpha_min), inner.alpha_max);

params5 = struct( ...
    'lambda1', lambda1, ...
    'lambda2', lambda2_effective, ...
    'lambda2_suggested', hp.lambda2_suggested, ...
    'alpha0',  alpha0_used, ...
    'max_iter', inner.max_iter, 'verbose', inner.verbose, ...
    'active_set_update_freq', 10, ...
    'alpha_max', inner.alpha_max, 'alpha_up', inner.alpha_up, ...
    'alpha_down', inner.alpha_down, 'alpha_grow_patience', inner.alpha_grow_patience, ...
    'obj_improve_tol', inner.obj_improve_tol, ...
    'weight_mode', inner.weight_mode, 'use_graph_laplacian', inner.use_graph_laplacian, ...
    'diag', inner.diag);

% 明确禁用 λ₂ 退火（更贴近“同题同标尺”对比）
params5.backtrack_beta           = 0.5;
params5.armijo_c1                = 1e-4;
params5.max_backtrack_per_iter   = 20;
params5.backtrack_patience       = Inf;   % 不触发降 λ₂
params5.lambda2_decay_factor     = 1.0;   % 不衰减
params5.lambda2_min              = params5.lambda2;

% ---------- (3) 跑 PGD ----------
[Gamma_tilde_star, prox_out] = module5_proximal(input_data_m5, params5);

% ---------- (4) 可选 recolor ----------
recolor = true; if isfield(opts,'recolor') && ~opts.recolor, recolor = false; end
recol = struct(); Omega_src = [];
if recolor && ~isempty(D_src)
    input_data_m8 = struct();
    input_data_m8.whitened_precision_matrices = Gamma_tilde_star;
    input_data_m8.whitening_matrices          = D_src;
    recolor_params = struct();
    recolor_params.force_hermitian = true;
    recolor_params.validate_spd    = true;
    recolor_params.verbose         = false;
    recolor_params.compute_quality_metrics = true;
    recolor_params.g_min_threshold = 1e-12;
    recolor_params.inv_error_tolerance = [];
    recolor_params.hermitian_tolerance  = 1e-12;
    recolor_params = set_defaults(recolor_params, struct()); % no-ops,占位保持接口

    recol = module8_recoloring(input_data_m8, recolor_params);
    if isfield(recol,'recolored_precision_matrices')
        Omega_src = recol.recolored_precision_matrices;
    end
end

% ---------- (5) 汇总输出 ----------
out = struct();
out.Gamma_tilde_star   = Gamma_tilde_star;
out.Omega_src          = Omega_src;           % 若未 recolor，则为 []
out.hp                 = hp;
out.params5            = params5;
out.proximal_results   = prox_out;
out.recoloring_results = recol;

% 附带“可复现实验”所需的最小信息（方便落盘）
out.meta = struct();
out.meta.F = F; out.meta.p = p; out.meta.used_jankova_override = logical(jankova.override);
if jankova.override
    out.meta.jankova = jankova;
end

end  % function adapter_pgd_source_only

% ===================== helpers =====================
function [C, F, p] = ensure_cell_covariances(args)
    if ~isfield(args,'Sjj_tilde') || isempty(args.Sjj_tilde)
        error('adapter_pgd_source_only:missing_Sjj','args.Sjj_tilde required');
    end
    X = args.Sjj_tilde;
    if iscell(X)
        C = X; F = numel(C); p = size(C{1},1);
    elseif isnumeric(X) && ndims(X)==3
        [p,~,F] = size(X); C = cell(F,1);
        for f=1:F, Cf = X(:,:,f); C{f} = 0.5*(Cf+Cf'); end
    else
        error('adapter_pgd_source_only:bad_Sjj','Sjj_tilde must be cell{F,1} or p×p×F numeric');
    end
end

function M = must_field(S, name, human)
    if ~isfield(S,name) || isempty(S.(name))
        error('adapter_pgd_source_only:missing_%s', name);
    end
    M = S.(name);
end

function C = ensure_cell_square(X, F, p, human)
    if iscell(X)
        assert(numel(X)==F, '%s must be cell{F,1}', human);
        C = X;
    elseif isnumeric(X) && ndims(X)==3 && size(X,1)==p && size(X,2)==p && size(X,3)==F
        C = cell(F,1); for f=1:F, C{f} = 0.5*(X(:,:,f)+X(:,:,f)'); end
    else
        error('adapter_pgd_source_only:bad_%s', human);
    end
end

function A = normalize_masks(val, F, p)
    if isempty(val), A = []; return; end
    if iscell(val)
        assert(numel(val)==F, 'A_masks cell length mismatch');
        A = cell(F,1);
        for f=1:F
            B = logical(val{f});
            assert(ismatrix(B) && all(size(B)==[p p]), 'mask size mismatch');
            A{f} = (B | B'); % 强制对称
        end
        return;
    end
    if isnumeric(val) || islogical(val)
        if ndims(val)==3 && size(val,1)==p && size(val,2)==p && size(val,3)==F
            A = cell(F,1);
            for f=1:F
                B = logical(val(:,:,f));
                A{f} = (B | B');
            end
            return;
        end
    end
    error('adapter_pgd_source_only:bad_masks','A_masks must be cell{F,1} or logical p×p×F');
end

function S = default_inner(p)
    S = struct();
    S.max_iter = 300;
    S.verbose  = true;
    S.alpha_max = 1.0;  S.alpha_min = 1e-6;
    S.alpha_up  = 1.5;  S.alpha_down = 0.5;
    S.alpha_grow_patience = 10;
    S.obj_improve_tol = 1e-6;
    S.weight_mode = 'matrix';
    S.use_graph_laplacian = true;
    S.diag = struct('enable', false, 'log_csv','', 'print_every',1,'update_every',1, ...
                    'live_plot', struct('enable',false,'f_view',1,'plot_every',5, ...
                                        'value_mode','abs','ground_truth_precision',[], ...
                                        'ground_truth_domain','source'));
end

function S = set_defaults(S, D)
    fn = fieldnames(D);
    for i=1:numel(fn)
        k = fn{i};
        if ~isfield(S,k) || isempty(S.(k)), S.(k) = D.(k); end
    end
end
