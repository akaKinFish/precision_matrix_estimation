function [Gamma_new, step_info] = module5_single_proximal_step(Gamma_current, gradient, alpha, active_mask, aux_data, params)
% MODULE5_SINGLE_PROXIMAL_STEP - Single proximal gradient update with robust PD safeguards
%
% Steps:
%   (1) Gradient step:   Gtmp = Gamma_k - alpha * grad
%   (2) Complex soft-threshold (off-diagonal)
%   (3) Active-set projection
%   (4) Hermitian symmetrization + real diagonal
%   (5) PD stabilization (tiny ridge + min-eig lifting if needed)
%   (6) PD check; if still not PD -> backtrack (shrink alpha) and retry
%
% Inputs:
%   Gamma_current : (p x p) current precision matrix (Hermitian PD expected)
%   gradient      : (p x p) smooth gradient at current iterate
%   alpha         : scalar step size
%   active_mask   : (p x p logical) active set mask (true = free variable)
%   aux_data      : struct with fields
%       .lambda2  : L1 penalty coefficient
%   params        : struct with fields (all optional)
%       .mode               : 'simplified' | 'joint' (passed to thresholding helpers)
%       .beta_backtrack     : (default 0.5) step shrink on backtrack
%       .max_backtrack      : (default 30)  max # of backtracks
%       .psd_tolerance      : (default 1e-12) (kept for compatibility)
%       .verbose            : (default false)
%       .ridge_scale        : (default 1e-10) ridge = ridge_scale * trace(|G|)/p
%       .ridge_min          : (default 1e-12) minimal ridge lower bound
%       .min_eig_floor      : (default 1e-8)  final floor for eigenvalues in last-resort
%       .penalize_diagonal  : (default false) whether to shrink diagonal in prox
%
% Outputs:
%   Gamma_new  : (p x p) next iterate (Hermitian PD)
%   step_info  : struct
%       .backtrack_count
%       .final_step_size
%       .objective_decrease      (first-order approx)
%       .psd_violations
%       .success
%       .matrix_change_norm
%       .relative_change
%       .final_is_hermitian
%       .final_diagonal_real
%       .final_condition_number

% -------------------- input checks --------------------
if nargin < 6
    error('module5_single_proximal_step:insufficient_input', ...
          'All 6 arguments are required');
end

p = size(Gamma_current, 1);
if ~isequal(size(Gamma_current), [p p]) || ~isequal(size(gradient), [p p])
    error('module5_single_proximal_step:dimension_mismatch', ...
          'Gamma_current and gradient must be square and same size');
end
if ~isequal(size(active_mask), [p p])
    error('module5_single_proximal_step:active_mask_size', ...
          'active_mask must be %dx%d', p, p);
end
if ~all(isfinite(Gamma_current(:))) || ~all(isfinite(gradient(:)))
    error('module5_single_proximal_step:invalid_values', ...
          'Gamma_current and gradient must contain finite values');
end
if ~(isscalar(alpha) && isfinite(alpha) && alpha > 0)
    error('module5_single_proximal_step:invalid_step_size', ...
          'alpha must be positive, finite scalar');
end
if ~isstruct(aux_data) || ~isfield(aux_data,'lambda2')
    error('module5_single_proximal_step:invalid_aux', ...
          'aux_data.lambda2 is required');
end

% -------------------- defaults --------------------
if ~isfield(params,'mode'),               params.mode = 'simplified'; end
if ~isfield(params,'beta_backtrack'),     params.beta_backtrack = 0.5; end
if ~isfield(params,'max_backtrack'),      params.max_backtrack = 30;  end
if ~isfield(params,'psd_tolerance'),      params.psd_tolerance = 1e-12; end
if ~isfield(params,'verbose'),            params.verbose = false; end
if ~isfield(params,'ridge_scale'),        params.ridge_scale = 1e-10; end
if ~isfield(params,'ridge_min'),          params.ridge_min   = 1e-12; end
if ~isfield(params,'min_eig_floor'),      params.min_eig_floor = 1e-8; end
if ~isfield(params,'penalize_diagonal'),  params.penalize_diagonal = false; end

lambda2 = aux_data.lambda2;

% -------------------- init outputs --------------------
step_info = struct();
step_info.backtrack_count        = 0;
step_info.psd_violations         = 0;
step_info.success                = false;
step_info.final_step_size        = alpha;
step_info.objective_decrease     = 0;
step_info.matrix_change_norm     = 0;
step_info.relative_change        = 0;
step_info.final_is_hermitian     = false;
step_info.final_diagonal_real    = false;
step_info.final_condition_number = Inf;

Gamma_new = Gamma_current; % fallback

% -------------------- backtracking loop --------------------
current_alpha = alpha;

for bt_iter = 0:params.max_backtrack

    % (1) gradient step
    Gtmp = Gamma_current - current_alpha * gradient;

    % (2) complex soft-thresholding
    %     NOTE: off-diagonal shrink by magnitude; diagonal optionally untouched
    tau = lambda2 * current_alpha;
    Gprox = module5_soft_threshold_complex(Gtmp, tau, params.mode);
    if ~params.penalize_diagonal
        % restore diagonal from pre-threshold values (no L1 on diag)
        Gprox(1:p+1:end) = Gtmp(1:p+1:end);
    end

    % (3) active-set projection (freeze inactive entries to current value)
    Gproj = module5_active_set_projection(Gprox, active_mask);

    % (4) Hermitian symmetrization + real diagonal
    Gsym = (Gproj + Gproj')/2;
    Gsym(1:p+1:end) = real(diag(Gsym));

    % (5) PD stabilization (tiny ridge first, then min-eig lifting if needed)
    [Gstab, was_pd] = local_stabilize_pd(Gsym, params.ridge_scale, params.ridge_min);
    if ~was_pd
        % try an extra eigen lifting (cheap for small p; safe for large p due to early chol)
        [Gstab, was_pd] = local_lift_min_eig(Gstab, params.min_eig_floor);
    end

    % (6) final PD check -> accept or backtrack
    if was_pd
        Gamma_new = Gstab;
        step_info.success         = true;
        step_info.final_step_size = current_alpha;
        break;
    else
        % shrink alpha and retry
        step_info.psd_violations  = step_info.psd_violations + 1;
        step_info.backtrack_count = step_info.backtrack_count + 1;

        if bt_iter == params.max_backtrack
            % last resort: hard eigen-lift to floor
            if params.verbose
                warning('module5_single_proximal_step:max_backtrack_reached', ...
                        'Maximum backtracking iterations reached, using eigenvalue-lifted solution');
            end
            Gamma_new = local_force_pd(Gsym, params.min_eig_floor);
            step_info.success         = true;
            step_info.final_step_size = current_alpha;
            break;
        end

        % reduce step
        current_alpha = current_alpha * params.beta_backtrack;
        continue;
    end
end

% -------------------- quality metrics --------------------
if step_info.success
    dG = Gamma_new - Gamma_current;
    step_info.matrix_change_norm  = norm(dG, 'fro');
    step_info.relative_change     = step_info.matrix_change_norm / max(norm(Gamma_current,'fro'),1e-12);
    % first-order objective decrease proxy: <grad, dG>
    step_info.objective_decrease  = real(trace(gradient' * dG));

    step_info.final_is_hermitian  = (norm(Gamma_new - Gamma_new','fro') <= 1e-12);
    step_info.final_diagonal_real = all(abs(imag(diag(Gamma_new))) < 1e-12);
    try
        step_info.final_condition_number = cond(Gamma_new);
    catch
        step_info.final_condition_number = Inf;
    end
else
    % keep defaults; Gamma_new already set to Gamma_current
end

end

% ==================== local helpers ====================

function [Gpd, is_pd] = local_stabilize_pd(G, ridge_scale, ridge_min)
% Add a tiny, scale-aware ridge and try chol.
% ridge = max(ridge_min, ridge_scale * trace(|G|)/p)
    p = size(G,1);
    trAbs = trace(abs(G));
    ridge = max(ridge_min, ridge_scale * max(trAbs/p, 1));
    Gpd = (G + G')/2;
    Gpd = Gpd + ridge * eye(p);
    [~,flag] = chol(Gpd);
    is_pd = (flag == 0);
end

function [Gpd, is_pd] = local_lift_min_eig(G, floor_val)
% If smallest eigenvalue < floor, shift up to floor.
    Gsym = (G + G')/2;
    try
        ev = eig(Gsym);
        lam_min = min(real(ev));
        if lam_min < floor_val
            shift = floor_val - lam_min + eps;
            Gpd = Gsym + shift * eye(size(Gsym,1));
        else
            Gpd = Gsym;
        end
    catch
        % If eig fails, apply conservative shift
        Gpd = Gsym + floor_val * eye(size(Gsym,1));
    end
    [~,flag] = chol(Gpd);
    is_pd = (flag == 0);
end

function Gpd = local_force_pd(G, floor_val)
% Unconditional eigen lifting to at least floor_val (used only at last resort).
    Gsym = (G + G')/2;
    [V,D] = eig(Gsym);
    d = real(diag(D));
    d = max(d, floor_val);
    Gpd = V * diag(d) * V';
    Gpd = (Gpd + Gpd')/2;
end
