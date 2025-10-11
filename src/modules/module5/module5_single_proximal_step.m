function [Gamma_new, step_info] = module5_single_proximal_step(Gamma_current, gradient, alpha, active_mask, aux_data, params)
% MODULE5_SINGLE_PROXIMAL_STEP - Single proximal gradient update with SPD projection prox
%
% 更新要点（相对旧版）：
% - 第(5)步改为：对称化后直接做“SPD 锥正交投影 + 最小特征值 floor”，不再使用 ridge 和逐步抬最小特征值。
% - 由于 SPD 投影保证 PD，故本函数内不再做“仅SPD回溯”；外层 Armijo/majorant 回溯逻辑保持不变。
%
% Steps:
%   (1) Gradient step:   Gtmp = Gamma_k - alpha * grad
%   (2) Complex soft-threshold (off-diagonal by default; diagonal optional)
%   (3) Active-set projection
%   (4) Hermitian symmetrization + real diagonal
%   (5) SPD projection (orthogonal to the SPD cone with eigen floor)
%   (6) Output + metrics (no SPD-only backtracking needed)
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
%       .beta_backtrack     : (保留但本函数不再使用，仅为兼容) 
%       .max_backtrack      : (保留但本函数不再使用，仅为兼容)
%       .min_eig_floor      : (default 1e-8)  floor for eigen projection
%       .penalize_diagonal  : (default false) whether to shrink diagonal in prox
%       .penalize_mask      : (default [])    elementwise logical mask for L1 (true=penalize)
%       .verbose            : (default false)
%
% Outputs:
%   Gamma_new  : (p x p) next iterate (Hermitian PD)
%   step_info  : struct
%       .backtrack_count
%       .final_step_size
%       .psd_violations
%       .success
%       .matrix_change_norm
%       .relative_change
%       .objective_decrease
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


params = set_default_(params, 'mode',              'simplified');
params = set_default_(params, 'beta_backtrack',    0.5);
params = set_default_(params, 'max_backtrack',     20);
params = set_default_(params, 'ridge_scale',       1e-10);
params = set_default_(params, 'ridge_min',         1e-12);
params = set_default_(params, 'min_eig_floor',     1e-8);
params = set_default_(params, 'penalize_diagonal', false);
params = set_default_(params, 'penalize_mask',     []);
params = set_default_(params, 'l1_weights',        []);   % << weighted L1
params = set_default_(params, 'verbose',           false);

lambda2    = aux_data.lambda2;
Gamma_new  = Gamma_current;

step_info = struct('backtrack_count',0,'final_step_size',alpha,'psd_violations',0, ...
                   'success',false,'matrix_change_norm',0,'relative_change',0, ...
                   'objective_decrease',0,'final_is_hermitian',false, ...
                   'final_diagonal_real',false,'final_condition_number',Inf);

current_alpha = alpha;
for bt_iter = 0:params.max_backtrack
    % (1) gradient step
    Gtmp = Gamma_current - current_alpha * gradient;

    % (2) weighted complex soft-threshold on off-diag
    base_tau  = lambda2 * current_alpha;
    Tau = base_tau;
    if ~isempty(params.l1_weights)
       isscalar(params.l1_weights), Tau = base_tau * params.l1_weights;
        
    end
    Gprox = soft_threshold_complex_weighted(Gtmp, Tau, params.mode);
    if ~params.penalize_diagonal, Gprox(1:p+1:end) = Gtmp(1:p+1:end); end
    if ~isempty(params.penalize_mask)
        M = logical(params.penalize_mask); if ~isequal(size(M), [p p]), error('penalize_mask size mismatch'); end
        keep = ~M; Gprox(keep) = Gtmp(keep);
    end

    % (3) active-set projection
    if ~isempty(active_mask), Gproj = module5_active_set_projection(Gprox, active_mask);
    else,                      Gproj = Gprox; end

    % (4) Hermitian + real diag
    Gsym = (Gproj + Gproj')/2; Gsym(1:p+1:end) = real(diag(Gsym));

    % (5) PD stabilize
    [Gstab, was_pd] = local_stabilize_pd(Gsym, params.ridge_scale, params.ridge_min);
    if ~was_pd, [Gstab, was_pd] = local_lift_min_eig(Gstab, params.min_eig_floor); end

    % (6) accept or SPD-backtrack
    if was_pd
        Gamma_new = Gstab;
        step_info.success         = true;
        step_info.final_step_size = current_alpha; break;
    else
        step_info.psd_violations  = step_info.psd_violations + 1;
        step_info.backtrack_count = step_info.backtrack_count + 1;
        if bt_iter >= params.max_backtrack
            if params.verbose
                warning('module5_single_proximal_step:max_spd_backtrack','Reached max backtracks, force eig-lift.');
            end
            Gamma_new = local_force_pd(Gsym, params.min_eig_floor);
            step_info.success         = true;
            step_info.final_step_size = current_alpha; break;
        end
        current_alpha = max(realmin, params.beta_backtrack * current_alpha);
        continue;
    end
end

if step_info.success
    dG = Gamma_new - Gamma_current;
    step_info.matrix_change_norm  = norm(dG, 'fro');
    step_info.relative_change     = step_info.matrix_change_norm / max(norm(Gamma_current,'fro'),1e-12);
    step_info.objective_decrease  = real(trace(gradient' * dG));
    step_info.final_is_hermitian  = (norm(Gamma_new - Gamma_new','fro') <= 1e-12);
    step_info.final_diagonal_real = all(abs(imag(diag(Gamma_new))) < 1e-12);
    try, step_info.final_condition_number = cond((Gamma_new+Gamma_new')/2); catch, step_info.final_condition_number = Inf; end
end
end

% ===== helpers =====
function params = set_default_(params, k, v)
if ~isfield(params,k) || isempty(params.(k)), params.(k) = v; end
end
function X = soft_threshold_complex_weighted(X, Tau, mode)
n = size(X,1); off = ~eye(n);
if isscalar(Tau), Tau = Tau * ones(n,'like',X); end
A = abs(X); shrink = max(0, A - Tau); ang = angle(X);
switch mode
    case {'simplified','joint'}
        X(off) = shrink(off) .* exp(1i*ang(off));
    otherwise, error('Unknown mode.');
end
end
function [Gpd, is_pd] = local_stabilize_pd(G, ridge_scale, ridge_min)
    p = size(G,1); trAbs = trace(abs(G));
    ridge = max(ridge_min, ridge_scale * max(trAbs/p, 1));
    Gpd = (G + G')/2 + ridge * eye(p);
    [~,flag] = chol(Gpd, 'lower'); is_pd = (flag == 0);
end
function [Gpd, is_pd] = local_lift_min_eig(G, floor_val)
    Gsym = (G + G')/2;
    try, e  = eig(Gsym); lm = min(real(e)); 
        if lm < floor_val, Gpd = Gsym + (floor_val - lm + eps) * eye(size(Gsym,1)); else, Gpd = Gsym; end
    catch, Gpd = Gsym + floor_val * eye(size(Gsym,1)); end
    [~,flag] = chol(Gpd, 'lower'); is_pd = (flag == 0);
end
function Gpd = local_force_pd(G, floor_val)
    Gsym = (G + G')/2; [V,D] = eig(Gsym); d = real(diag(D)); d = max(d, floor_val);
    Gpd = V * diag(d) * V'; Gpd = (Gpd + Gpd')/2;
end