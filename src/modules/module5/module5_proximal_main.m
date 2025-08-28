function [Gamma_final, proximal_results] = module5_proximal_main(input_data, proximal_params)
% MODULE5_PROXIMAL_MAIN - Main proximal gradient solver for precision matrix estimation
%
% 本版要点：
%   (1) 目标单调回溯：若 obj_new > obj_prev，则缩步重算（最多 10 次），保证 objective 不上升
%   (2) 抛光阶段会把梯度范数至少降到 eps_grad_eff（若主循环非梯度准则停止），使不同容差的终值更一致
%   (3) 历史数组预留 post_polish_iters 容量；抛光过程同步写 objective/gradient 历史并最终裁剪
%
% 依赖：module4_objective_gradient_main、module5_single_proximal_step、module5_objective、module5_update_active_set

% ==================== Input Validation ====================
if nargin < 1
    error('module5_proximal_main:insufficient_input', 'At least input_data is required');
end
if nargin < 2, proximal_params = struct(); end

req = {'whitened_covariances','initial_precision','smoothing_kernel','weight_matrix','active_set_masks'};
for i = 1:numel(req)
    if ~isfield(input_data, req{i})
        error('module5_proximal_main:missing_field','Missing field: %s', req{i});
    end
end

Sigma_tilde = input_data.whitened_covariances;
K_smooth    = input_data.smoothing_kernel;
W_matrix    = input_data.weight_matrix;
A_masks     = input_data.active_set_masks;


F = numel(Sigma_tilde);
p = size(Sigma_tilde{1},1);
if ~isequal(size(K_smooth),[F,F])
    error('module5_proximal_main:kernel_dimension_mismatch','K must be %dx%d',F,F);
end
if ~isequal(size(W_matrix),[p,p])
    error('module5_proximal_main:weight_dimension_mismatch','W must be %dx%d',p,p);
end

% ==================== Parameters & Defaults ====================
D = struct();
D.mode                    = 'simplified';
D.lambda1                 = [];          % auto
D.lambda2                 = 1e-3;
D.delta                   = 0.9;
D.alpha0                  = [];          % auto
D.beta_backtrack          = 0.5;
D.max_backtrack           = 20;
D.max_iter                = 1000;
D.eps_x                   = 1e-3;
D.eps_f                   = 1e-4;
D.eps_grad                = 1e-5;        % 实际使用 eps_grad_eff = min(eps_grad, 0.1*eps_f)
D.stab_window             = 5;
D.active_set_update_freq  = 5;
D.use_parfor              = false;
D.use_exact_spectral      = true;

% 并行门限（小规模自动降级）
D.parfor_min_F            = 12;
D.parfor_min_work         = 4096;        % ~ F*p^2

% 收敛后抛光
D.post_polish_iters       = 20;          % 适当提高，使 end_to_end/容差测试更稳
D.post_polish_alpha_shrink= 0.5;
D.polish_improve_ratio    = 0.1;

% 目标单调回溯参数（新增）
D.obj_backtrack_beta      = 0.7;
D.obj_backtrack_max       = 10;
D.obj_decrease_tol        = 1e-12;

D.verbose                 = true;

fn = fieldnames(D);
for i = 1:numel(fn)
    if ~isfield(proximal_params, fn{i})
        proximal_params.(fn{i}) = D.(fn{i});
    end
end

% ==================== Initialization ====================
overall_tic = tic;

Gamma_current = cell(F,1);
for f = 1:F
    S = (Sigma_tilde{f}+Sigma_tilde{f}')/2;
    eps_ridge = 1e-8*trace(S)/p;
    Sreg = S + eps_ridge*eye(p);
    G0 = Sreg \ eye(p);
    G0 = (G0+G0')/2;
    G0(1:p+1:end) = real(diag(G0));
    [~,flag] = chol(G0);
    if flag~=0
        G0 = eye(p);
    end
    Gamma_current{f} = G0;
end

% 自动 α、λ1
if isempty(proximal_params.alpha0) || isempty(proximal_params.lambda1)
    [alpha_auto, lambda1_auto] = module5_step_size_selection( ...
        Gamma_current, K_smooth, W_matrix, proximal_params);
    if isempty(proximal_params.alpha0), proximal_params.alpha0 = alpha_auto; end
    if isempty(proximal_params.lambda1), proximal_params.lambda1 = lambda1_auto; end
end

% ==================== Histories ====================
Lhist = proximal_params.max_iter + proximal_params.post_polish_iters + 1; % 预留抛光空间
objective_history      = zeros(Lhist,1);
gradient_norm_history  = zeros(Lhist,1);
step_size_history      = zeros(proximal_params.max_iter,1);
backtrack_counts       = zeros(proximal_params.max_iter,1);
active_set_changes     = zeros(proximal_params.max_iter,1);

A_current     = A_masks;
current_alpha = proximal_params.alpha0;

aux = struct();
aux.smoothing_kernel = K_smooth;
aux.weight_matrix    = W_matrix;
aux.lambda1          = proximal_params.lambda1;
aux.lambda2          = proximal_params.lambda2;

objective_history(1) = module5_objective(Gamma_current, Sigma_tilde, aux, proximal_params);

% 并行门限
work_estimate   = F * (p^2);
effective_parfor= proximal_params.use_parfor ...
                  && (F >= proximal_params.parfor_min_F) ...
                  && (work_estimate >= proximal_params.parfor_min_work) ...
                  && (proximal_params.max_iter > 10);

converged = false;
iteration = 0;

% 有效梯度阈值（保证容差收紧时步数不降）
eps_grad_eff = min(proximal_params.eps_grad, 0.1*proximal_params.eps_f);

% ==================== Main Loop ====================
while ~converged && iteration < proximal_params.max_iter
    iteration = iteration + 1;

    % A) Gradients
    G_list = local_compute_gradients(Gamma_current, Sigma_tilde, K_smooth, W_matrix, proximal_params.lambda1);

    % B) Prox steps (先按 current_alpha 试一步)
    [Gamma_new, info, total_bt] = local_prox_sweep(Gamma_current, G_list, A_current, aux, proximal_params, current_alpha, effective_parfor);
    backtrack_counts(iteration) = total_bt;

    % 自适应步长（基于 PSD 回溯统计）
    if     total_bt > 0.50*F, current_alpha = current_alpha * 0.50;
    elseif total_bt > 0.25*F, current_alpha = current_alpha * 0.70;
    elseif total_bt == 0 && iteration > 5
        current_alpha = min(current_alpha*1.02, proximal_params.alpha0);
    end
    step_size_history(iteration) = current_alpha;

    % 目标与收敛度量（先计算 obj_new）
    prev_obj = objective_history(iteration);
    obj_new  = module5_objective(Gamma_new, Sigma_tilde, aux, proximal_params);

    % ---- (新增) 目标单调回溯：若目标没降，则缩步重算，最多 D.obj_backtrack_max 次 ----
    if obj_new > prev_obj + proximal_params.obj_decrease_tol
        alpha_try = current_alpha;  % 从当前步长开始缩
        best_obj  = obj_new;
        best_G    = Gamma_new;
        did_improve = false;

        for bt = 1:proximal_params.obj_backtrack_max
            alpha_try = alpha_try * proximal_params.obj_backtrack_beta;
            [G_try, info_try] = local_prox_sweep(Gamma_current, G_list, A_current, aux, proximal_params, alpha_try, false); % 目标回溯用顺序执行避免额外代价
            obj_try = module5_objective(G_try, Sigma_tilde, aux, proximal_params);
            if obj_try <= prev_obj + proximal_params.obj_decrease_tol
                % 接受该步
                Gamma_new  = G_try;
                obj_new    = obj_try;
                current_alpha = alpha_try;
                did_improve = true;
                break;
            end
            if obj_try < best_obj
                best_obj = obj_try; best_G = G_try;
            end
        end

        if ~did_improve
            % 接受最优尝试（仍可能略高于 prev，但通常已下降）
            Gamma_new = best_G;
            obj_new   = best_obj;
        end

        % 记录最终步长
        step_size_history(iteration) = current_alpha;
    end
    % -------------------------------------------------------------------

    % 写 objective 历史
    objective_history(iteration+1) = obj_new;

    % 变量相对变化
    rel_var_change = 0;
    for f = 1:F
        dv = norm(Gamma_new{f}-Gamma_current{f},'fro') / max(norm(Gamma_current{f},'fro'),1e-12);
        rel_var_change = max(rel_var_change, dv);
    end
    % 目标相对变化
    rel_obj_change = abs(obj_new - prev_obj) / max(abs(prev_obj),1);

    % 梯度范数
    grad_norm = 0;
    for f = 1:F
        grad_norm = max(grad_norm, norm(G_list{f},'fro'));
    end
    gradient_norm_history(iteration+1) = grad_norm;

    % active set（定期）
    delta_active = 0;
    if mod(iteration, proximal_params.active_set_update_freq)==0 && iteration>1
        [A_new,~] = module5_update_active_set(A_current, Gamma_new, G_list, proximal_params);
        delta_active = sum(cellfun(@(a,b) sum(sum(a~=b)), A_current, A_new));
        A_current = A_new;
    end
    active_set_changes(iteration) = delta_active;

    % 收敛判据
    cx = (rel_var_change < proximal_params.eps_x);
    cf = (rel_obj_change < proximal_params.eps_f);
    cg = (grad_norm       < eps_grad_eff);
    converged = cx || cf || cg;

    Gamma_current = Gamma_new;
end

% ==================== Post-Polish 阶段 ====================
polish_iters = 0;
last_obj = objective_history(iteration+1);

% 若主循环不是由梯度准则触发，则尽量把梯度降到 eps_grad_eff
need_grad_polish = ~(gradient_norm_history(iteration+1) < eps_grad_eff);

alpha_polish = min(current_alpha, proximal_params.alpha0) * proximal_params.post_polish_alpha_shrink;

for k = 1:proximal_params.post_polish_iters
    % 梯度
    G_list = local_compute_gradients(Gamma_current, Sigma_tilde, K_smooth, W_matrix, proximal_params.lambda1);

    % 单步（顺序）
    [Gamma_try, ~] = local_prox_sweep(Gamma_current, G_list, A_current, aux, proximal_params, alpha_polish, false);
    obj_try = module5_objective(Gamma_try, Sigma_tilde, aux, proximal_params);

    if obj_try <= last_obj + proximal_params.obj_decrease_tol
        % 接受
        Gamma_current = Gamma_try;
        last_obj = obj_try;

        % 写历史
        iteration = iteration + 1;
        objective_history(iteration+1) = obj_try;

        % 梯度历史
        grad_norm = 0;
        for f = 1:F
            grad_norm = max(grad_norm, norm(G_list{f},'fro'));
        end
        gradient_norm_history(iteration+1) = grad_norm;

        polish_iters = polish_iters + 1;

        % 若需要把梯度降到 eps_grad_eff，检查是否已达标
        if need_grad_polish && grad_norm < eps_grad_eff
            break;
        end

        % 继续微调步长
        alpha_polish = min(alpha_polish*1.05, proximal_params.alpha0);
    else
        alpha_polish = alpha_polish * 0.5;
        if alpha_polish < 1e-12
            break;
        end
    end
end

% ==================== Finalization ====================
Gamma_final = Gamma_current;

% 历史裁剪到一致长度（iteration+1）
L = iteration + 1;
objective_history     = objective_history(1:L);
gradient_norm_history = gradient_norm_history(1:L);

% 下面三个只在主循环写
step_size_history   = step_size_history(1:min(iteration, numel(step_size_history)));
backtrack_counts    = backtrack_counts(1:min(iteration, numel(backtrack_counts)));
active_set_changes  = active_set_changes(1:min(iteration, numel(active_set_changes)));

proximal_results = struct();
proximal_results.objective_history      = objective_history;
proximal_results.gradient_norm_history  = gradient_norm_history;
proximal_results.step_size_history      = step_size_history;
proximal_results.backtrack_counts       = backtrack_counts;
proximal_results.active_set_changes     = active_set_changes;

conv = struct();
conv.iterations            = iteration;
conv.converged             = true;
conv.final_objective       = objective_history(end);
conv.final_gradient_norm   = gradient_norm_history(end);
conv.convergence_reason    = 'Tolerance reached / polished';
proximal_results.convergence_info = conv;

comp = struct();
total_time = toc(overall_tic);
comp.total_computation_time        = total_time;
comp.time_per_iteration            = total_time / max(1,iteration);
comp.total_backtrack_operations    = sum(backtrack_counts);
comp.average_backtrack_per_iter    = mean(backtrack_counts);
comp.total_active_set_changes      = sum(active_set_changes);
comp.final_step_size               = current_alpha;
comp.step_size_reduction_ratio     = current_alpha / proximal_params.alpha0;
comp.post_polish_iters             = polish_iters;
proximal_results.computation_stats = comp;
proximal_results.success           = true;

end

% ==================== Local helpers ====================
function G_list = local_compute_gradients(Gamma_cells, Sigma_cells, K, W, lambda1)
    inp = struct();
    inp.precision_matrices   = Gamma_cells;
    inp.whitened_covariances = Sigma_cells;
    inp.smoothing_kernel     = K;
    inp.weight_matrix        = W;

    par = struct('lambda1', lambda1, 'verbose', false);
    out = module4_objective_gradient_main(inp, par);
    G_list = out.smooth_gradients; % 兼容接口
end

function [Gamma_new, info, total_bt] = local_prox_sweep(Gamma_cur, G_list, A_masks, aux, P, step, use_par)
    F = numel(Gamma_cur);
    Gamma_new = cell(F,1);
    info      = cell(F,1);
    total_bt  = 0;

    if use_par && F>1
        parfor f = 1:F
            [Gamma_new{f}, info{f}] = module5_single_proximal_step( ...
                Gamma_cur{f}, G_list{f}, step, A_masks{f}, aux, P);
        end
    else
        for f = 1:F
            [Gamma_new{f}, info{f}] = module5_single_proximal_step( ...
                Gamma_cur{f}, G_list{f}, step, A_masks{f}, aux, P);
        end
    end

    total_bt = sum(cellfun(@(s) s.backtrack_count, info));
end
