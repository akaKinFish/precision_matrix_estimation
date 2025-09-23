function [Gamma_new, step_info] = module5_single_proximal_step(Gamma_current, gradient, alpha, active_mask, aux_data, params)
% MODULE5_SINGLE_PROXIMAL_STEP - Single proximal gradient update with robust PD safeguards
%
% Steps:
%   (1) Gradient step:   Gtmp = Gamma_k - alpha * grad
%   (2) Complex soft-threshold (off-diagonal by default; diagonal optional)
%   (3) Active-set projection
%   (4) Hermitian symmetrization + real diagonal
%   (5) PD stabilization (ridge first, then min-eig lifting if needed)
%   (6) PD check; if still not PD -> backtrack (shrink alpha) and retry (SPD-only backtracking)
%
% Inputs:
%   Gamma_current : (p x p) current precision matrix (Hermitian PD expected)
%   gradient      : (p x p) smooth gradient at current iterate
%   alpha         : scalar step size (provided by outer loop)
%   active_mask   : (p x p logical) active set mask (true = free variable), or []
%   aux_data      : struct with fields
%       .lambda2          : L1 penalty coefficient (required)
%   params        : struct (all optional)
%       .mode               : 'simplified' | 'joint' (passed to threshold helper)
%       .beta_backtrack     : (default 0.5) step shrink on SPD-backtrack
%       .max_backtrack      : (default 20)  max # of SPD backtracks here
%       .ridge_scale        : (default 1e-10) ridge = ridge_scale * max(trace(|G|)/p,1)
%       .ridge_min          : (default 1e-12) minimal ridge lower bound
%       .min_eig_floor      : (default 1e-8)  floor for eigen lifting
%       .penalize_diagonal  : (default false) whether to shrink diagonal in prox
%       .penalize_mask      : (default [])    elementwise logical mask for L1 (true=penalize)
%       .verbose            : (default false)
%
% Outputs:
%   Gamma_new  : (p x p) next iterate (Hermitian PD)
%   step_info  : struct
%       .backtrack_count       % SPD safeguard backtracks
%       .final_step_size       % alpha actually used here
%       .psd_violations        % #non-PD hits encountered
%       .success               % whether we produced PD output
%       .matrix_change_norm    % ||ΔG||_F
%       .relative_change       % ||ΔG||_F / ||G||_F
%       .objective_decrease    % <grad, ΔG> (first-order proxy)
%       .final_is_hermitian
%       .final_diagonal_real
%       .final_condition_number

% -------------------- input checks --------------------
if nargin < 6
    error('module5_single_proximal_step:insufficient_input', 'All 6 arguments are required');
end
p = size(Gamma_current, 1);
assert(isequal(size(Gamma_current), [p p]) && isequal(size(gradient), [p p]), ...
    'Gamma_current and gradient must be square and same size');
if ~isempty(active_mask)
    assert(isequal(size(active_mask), [p p]), 'active_mask must be %dx%d', p, p);
end
if ~(isscalar(alpha) && isfinite(alpha) && alpha > 0)
    error('module5_single_proximal_step:invalid_step_size', 'alpha must be positive, finite scalar');
end
if ~isstruct(aux_data) || ~isfield(aux_data,'lambda2')
    error('module5_single_proximal_step:invalid_aux', 'aux_data.lambda2 is required');
end
assert(all(isfinite(Gamma_current(:))) && all(isfinite(gradient(:))), ...
    'Gamma_current and gradient must contain finite values');

% -------------------- defaults --------------------
params = set_default_(params, 'mode',              'simplified');
params = set_default_(params, 'beta_backtrack',    0.5);
params = set_default_(params, 'max_backtrack',     20);
params = set_default_(params, 'ridge_scale',       1e-10);
params = set_default_(params, 'ridge_min',         1e-12);
params = set_default_(params, 'min_eig_floor',     1e-8);
params = set_default_(params, 'penalize_diagonal', false);
params = set_default_(params, 'penalize_mask',     []);
params = set_default_(params, 'verbose',           false);

lambda2    = aux_data.lambda2;
Gamma_new  = Gamma_current; % fallback

% -------------------- init step_info --------------------
step_info = struct('backtrack_count',0,'final_step_size',alpha,'psd_violations',0, ...
                   'success',false,'matrix_change_norm',0,'relative_change',0, ...
                   'objective_decrease',0,'final_is_hermitian',false, ...
                   'final_diagonal_real',false,'final_condition_number',Inf);

% -------------------- SPD-oriented backtracking loop --------------------
current_alpha = alpha;

for bt_iter = 0:params.max_backtrack

    % (1) gradient step
    Gtmp = Gamma_current - current_alpha * gradient;

    % (2) complex soft-thresholding（默认不惩罚对角）
    tau  = lambda2 * current_alpha;
    Gprox = module5_soft_threshold_complex(Gtmp, tau, params.mode);      % 需存在于你的工程
    if ~params.penalize_diagonal
        Gprox(1:p+1:end) = Gtmp(1:p+1:end);
    end
    if ~isempty(params.penalize_mask)
        % 仅对 mask==true 的元素保留 shrink，其他回退到未收缩值
        M = logical(params.penalize_mask);
        if ~isequal(size(M), [p p]), error('penalize_mask size mismatch'); end
        keep = ~M;  % 不惩罚的元素回退
        Gprox(keep) = Gtmp(keep);
    end

    % (3) active-set 投影（冻结 inactive 位置到当前值）
    if ~isempty(active_mask)
        Gproj = module5_active_set_projection(Gprox, active_mask);       % 需存在于你的工程
    else
        Gproj = Gprox;
    end

    % (4) 共轭对称化 + 实对角
    Gsym = (Gproj + Gproj')/2;
    Gsym(1:p+1:end) = real(diag(Gsym));

    % (5) PD 稳定化（ridge -> min-eig 提升）
    [Gstab, was_pd] = local_stabilize_pd(Gsym, params.ridge_scale, params.ridge_min);
    if ~was_pd
        [Gstab, was_pd] = local_lift_min_eig(Gstab, params.min_eig_floor);
    end

    % (6) PD 检查：通过则接受；否则仅做“SPD backtrack”（缩 alpha）重试
    if was_pd
        Gamma_new = Gstab;
        step_info.success         = true;
        step_info.final_step_size = current_alpha;
        break;
    else
        step_info.psd_violations  = step_info.psd_violations + 1;
        step_info.backtrack_count = step_info.backtrack_count + 1;

        if bt_iter >= params.max_backtrack
            % 兜底：强制特征值抬升，避免死循环
            if params.verbose
                warning('module5_single_proximal_step:max_spd_backtrack', ...
                        'SPD backtrack reached max=%d, forcing eig-lifted PD.', params.max_backtrack);
            end
            Gamma_new = local_force_pd(Gsym, params.min_eig_floor);
            step_info.success         = true;
            step_info.final_step_size = current_alpha;
            break;
        end

        % 仅用于 SPD 的回退（外层 Armijo 回溯由 main 控制）
        current_alpha = max(realmin, params.beta_backtrack * current_alpha);
        continue;
    end
end

% -------------------- metrics --------------------
if step_info.success
    dG = Gamma_new - Gamma_current;
    step_info.matrix_change_norm  = norm(dG, 'fro');
    step_info.relative_change     = step_info.matrix_change_norm / max(norm(Gamma_current,'fro'),1e-12);
    step_info.objective_decrease  = real(trace(gradient' * dG));          % 一阶下降代理
    step_info.final_is_hermitian  = (norm(Gamma_new - Gamma_new','fro') <= 1e-12);
    step_info.final_diagonal_real = all(abs(imag(diag(Gamma_new))) < 1e-12);
    try
        step_info.final_condition_number = cond((Gamma_new+Gamma_new')/2);
    catch
        step_info.final_condition_number = Inf;
    end
end
end

% ==================== local helpers ====================
function params = set_default_(params, k, v)
    if ~isfield(params,k) || isempty(params.(k)), params.(k) = v; end
end

function [Gpd, is_pd] = local_stabilize_pd(G, ridge_scale, ridge_min)
% scale-aware ridge then chol 检查
    p = size(G,1);
    trAbs = trace(abs(G));
    ridge = max(ridge_min, ridge_scale * max(trAbs/p, 1));
    Gpd = (G + G')/2 + ridge * eye(p);
    [~,flag] = chol(Gpd, 'lower');
    is_pd = (flag == 0);
end

function [Gpd, is_pd] = local_lift_min_eig(G, floor_val)
% 若最小特征值 < floor，则整体平移至 floor
    Gsym = (G + G')/2;
    try
        e  = eig(Gsym);
        lm = min(real(e));
        if lm < floor_val
            Gpd = Gsym + (floor_val - lm + eps) * eye(size(Gsym,1));
        else
            Gpd = Gsym;
        end
    catch
        % 极端情况下 eig 失败，保守平移
        Gpd = Gsym + floor_val * eye(size(Gsym,1));
    end
    [~,flag] = chol(Gpd, 'lower');
    is_pd = (flag == 0);
end

function Gpd = local_force_pd(G, floor_val)
% 最后兜底：按特征分解把所有特征值抬到 floor 以上
    Gsym = (G + G')/2;
    [V,D] = eig(Gsym);
    d = real(diag(D)); d = max(d, floor_val);
    Gpd = V * diag(d) * V'; Gpd = (Gpd + Gpd')/2;
end
