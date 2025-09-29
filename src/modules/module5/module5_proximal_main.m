function [Gamma_cells, results] = module5_proximal_main(input_data, params)
% MODULE5_PROXIMAL_MAIN
% Proximal gradient solver for precision matrices {Gamma_f} in the whitened domain.
% == 改进点（与旧版兼容）==
%   • BB 步长：迭代级初值（基于上一迭代 ΔΓ, Δ∇），减少无谓回溯
%   • Armijo 外层回溯：f(new) <= f(old) - c1 * alpha * ||grad||^2
%   • 步长下限 alpha_min，避免几何缩到机器噪声级
%   • 回溯累计触发“优先降 λ2”（退火），而非继续缩步长
%   • CSV/打印格式与旧版保持一致

% ------------ Parse inputs ------------
Sigma  = input_data.whitened_covariances;
Gamma0 = input_data.initial_precision;
K      = input_data.smoothing_kernel;
W      = input_data.weight_matrix;

F = numel(Sigma); n = size(Sigma{1},1);

% Optional passthrough
A_masks = [];
if isfield(input_data,'active_set_mask')
    A_masks = input_data.active_set_mask;
end
D_whiten = [];
if isfield(input_data,'whitening_matrices') && ~isempty(input_data.whitening_matrices)
    D_whiten = input_data.whitening_matrices;   % {F×1}, for GT -> whitened viz
end

% params (keep legacy names), 新增项给默认值保证后向兼容
lambda1     = getf(params,'lambda1',0);
lambda2_eff = getf(params,'lambda2',0);               % effective value (REQUIRED)
lambda2_sug = getf(params,'lambda2_suggested',[]);    % optional provenance
pen_diag    = getf(params,'penalize_diagonal',false); % L1 on diagonal?

alpha0      = getf(params,'alpha0',1e-3);
alpha_min   = getf(params,'alpha_min',1e-6);
alpha_max   = getf(params,'alpha_max',1.0);
armijo_c1   = getf(params,'armijo_c1',1e-4);
bt_beta     = getf(params,'backtrack_beta',0.5);
bt_max_iter = getf(params,'max_backtrack_per_iter',20);

% 回溯累计触发退火（降 λ2）
anneal_patience = getf(params,'backtrack_patience',10);
anneal_factor   = getf(params,'lambda2_decay_factor',0.9);
lambda2_min     = getf(params,'lambda2_min', max(1e-6, 0.01*lambda2_eff));

max_iter = getf(params,'max_iter',300);
verbose  = getf(params,'verbose',true);

weight_mode = lower(getf(params,'weight_mode','matrix'));
use_L       = getf(params,'use_graph_laplacian',true);

diagopt = struct('enable',false,'log_csv','','print_every',1,'update_every',1);
if isfield(params,'diag') && isstruct(params.diag)
    diagopt = set_defaults(params.diag, diagopt);
end

% live-plot defaults
viz = struct('enable',false,'f_view',1,'plot_every',5,'value_mode','abs', ...
             'ground_truth_precision',[],'ground_truth_domain','source');
if isfield(diagopt,'live_plot') && isstruct(diagopt.live_plot)
    viz = set_defaults(diagopt.live_plot, viz);
end

% ------------ CSV init (header + CONFIG) ------------
csv_fid = -1; csv_header_written = false; csv_config_written = false;
if diagopt.enable && isfield(diagopt,'log_csv') && ~isempty(diagopt.log_csv)
    is_new_file = ~exist(diagopt.log_csv,'file');
    csv_fid = fopen(diagopt.log_csv, 'a');
    if csv_fid == -1
        warning('module5_proximal_main:csv_open_failed','Failed to open CSV at %s', diagopt.log_csv);
    else
        if is_new_file
            fprintf(csv_fid, 'row_type,iter,obj,loglik,smooth,l1,spatial,step_size,alpha,grad_norm,lambda2_effective\n');
            csv_header_written = true;
        end
        fprintf(csv_fid, 'CONFIG,NA,NA,NA,NA,NA,NA,NA,NA,NA,%.8g\n', lambda2_eff);
        if ~isempty(lambda2_sug)
            fprintf(csv_fid, '# lambda2_suggested=%.8g\n', lambda2_sug);
        else
            fprintf(csv_fid, '# lambda2_suggested=NA\n');
        end
        fprintf(csv_fid, '# penalize_diagonal=%d\n', pen_diag);
        csv_config_written = true;
    end
end

% console echo if overridden
if ~isempty(lambda2_sug) && isfinite(lambda2_sug) && abs(lambda2_sug - lambda2_eff) > 0
    fprintf('[module5] lambda2 overridden: suggested=%.6g -> effective=%.6g\n', lambda2_sug, lambda2_eff);
else
    if verbose
        fprintf('[module5] lambda2 effective=%.6g (no override)\n', lambda2_eff);
    end
end
if verbose
    fprintf('[module5] penalize_diagonal = %d\n', pen_diag);
end

% ------------ Initialization ------------
Gamma = Gamma0;
alpha = min(max(alpha0, alpha_min), alpha_max);

% %%% FIX: 对 warm start/任意初值做“对称化 + 对角归一 + SPD 投影”（每个频点）
for f=1:F
    Gamma{f} = project_spd_(Gamma{f}, 1e-8, true);  % eig-floor=1e-8, 规范化 diag->1
end

% objective decomposition helpers
aux = struct('lambda1',lambda1,'lambda2',lambda2_eff,'weight_matrix',W,'smoothing_kernel',K);
params_obj = struct('weight_mode',weight_mode,'use_graph_laplacian',use_L, ...
                    'penalize_diagonal', pen_diag);
if isfield(params,'lambda3'), params_obj.lambda3 = params.lambda3; end
if isfield(params,'spatial_graph_matrix'),        params_obj.spatial_graph_matrix = params.spatial_graph_matrix; end
if isfield(params,'spatial_graph_is_laplacian'),  params_obj.spatial_graph_is_laplacian = params.spatial_graph_is_laplacian; end
if isfield(params,'spatial_weight_mode'),         params_obj.spatial_weight_mode = params.spatial_weight_mode; end

% gradient module params
grad_params = struct('lambda1',lambda1,'weight_mode',weight_mode,'use_graph_laplacian',use_L, ...
                     'force_hermitian',true,'symmetrization_tolerance',1e-12);
if isfield(params,'lambda3'), grad_params.lambda3 = params.lambda3; end
if isfield(params,'spatial_graph_matrix'),        grad_params.spatial_graph_matrix = params.spatial_graph_matrix; end
if isfield(params,'spatial_graph_is_laplacian'),  grad_params.spatial_graph_is_laplacian = params.spatial_graph_is_laplacian; end
if isfield(params,'spatial_weight_mode'),         grad_params.spatial_weight_mode = params.spatial_weight_mode; end

% histories
objective_history       = zeros(max_iter,1);
gradient_norm_history   = zeros(max_iter,1);
step_size_history       = zeros(max_iter,1);
backtrack_counts        = zeros(max_iter,1);
active_set_changes      = [];  % optional

% Live figure state (optional)
GT_cells = [];
if viz.enable && ~isempty(viz.ground_truth_precision)
    GT_cells = coerce_cov_cell_local(viz.ground_truth_precision, F);
    switch lower(string(viz.ground_truth_domain))
        case "source"
            if isempty(D_whiten)
                warning('ground_truth_domain=source but input_data.whitening_matrices not provided; using given GT as-is for visualization.');
            else
                for f = 1:F
                    Df = D_whiten{f};
                    GT_cells{f} = Df * GT_cells{f} * Df;
                end
            end
        case "whitened"
        otherwise
            warning('Unknown ground_truth_domain=%s; using GT as-is.', string(viz.ground_truth_domain));
    end
end
viz_state = [];
if viz.enable
    viz_state = init_live_plot_(F, n, viz, GT_cells);
end

% ------------ 初值目标（带兜底） ------------
try
    [loglik_val, smooth_val, l1_val, aux_terms] = module5_objective_terms(Gamma, Sigma, aux, params_obj);
catch ME
    % %%% NEW: 若仍因非 PD 失败，抬高地板再试一次
    if contains(ME.message,'objective_terms')
        for f=1:F, Gamma{f} = project_spd_(Gamma{f}, 1e-6, true); end
        [loglik_val, smooth_val, l1_val, aux_terms] = module5_objective_terms(Gamma, Sigma, aux, params_obj);
    else
        rethrow(ME);
    end
end
spatial_val = 0;
if isfield(aux_terms,'spatial_term'), spatial_val = aux_terms.spatial_term; end
obj_val = loglik_val + smooth_val + l1_val + spatial_val;

best_obj = obj_val;
if verbose
    fprintf('[module5] it=%4d  obj=%.6e  (loglik=%.6e, smooth=%.6e, l1=%.6e, spatial=%.6e)  alpha=%g  ||grad||=N/A\n', ...
        0, obj_val, loglik_val, smooth_val, l1_val, spatial_val, alpha);
end
if csv_fid~=-1
    fprintf(csv_fid, 'STEP,%d,%.10g,%.10g,%.10g,%.10g,%.10g,NA,%.6g,NA,%.8g\n', ...
        0, obj_val, loglik_val, smooth_val, l1_val, spatial_val, alpha, lambda2_eff);
end

% BB 所需缓存
prev_Gamma = []; prev_grad = [];

% 回溯累计（用于退火）
bt_cum = 0;

% ------------ Main loop ------------
for it = 1:max_iter
    % 1) 计算梯度
    input = struct();
    input.precision_matrices   = Gamma;
    input.whitened_covariances = Sigma;
    input.smoothing_kernel     = K;
    input.weight_matrix        = W;

    gres  = module4_objective_gradient_main(input, grad_params);
    Ggrad = gres.smooth_gradients;
    grad_norm = aggregate_grad_norm_(Ggrad);

    % 2) BB 步长（迭代级初值）
    if isempty(prev_Gamma) || isempty(prev_grad)
        alpha_try = min(alpha_max, max(alpha_min, alpha));
    else
        [alpha_bb1, alpha_bb2] = bb_stepsize_(Gamma, prev_Gamma, Ggrad, prev_grad);
        alpha_try = alpha_bb1;                          % BB1 通常更稳
        if ~isfinite(alpha_try) || alpha_try<=0
            alpha_try = alpha;                          % 退回上次
        end
        alpha_try = min(alpha_max, max(alpha_min, alpha_try));
    end

    % 3) 外层 Armijo 回溯（与内部 SPD 保护相独立）
    accept = false; bt_cnt = 0; obj_new = NaN;
    while ~accept
        % forward 梯度步
        Gamma_tmp = cell(F,1);
        for f=1:F, Gamma_tmp{f} = Gamma{f} - alpha_try * Ggrad{f}; end

        % proximal L1（含可选对角、活动集与 SPD 微加载）
        Gamma_cand = Gamma_tmp;
        for f=1:F
            Gc = Gamma_tmp{f};

            % off-diagonal shrink
            U = triu(Gc, 1);
            U = sign(U) .* max(abs(U) - alpha_try*lambda2_eff, 0);
            Gc = diag(diag(Gc)) + U + U';

            % diagonal shrink（可选）
            if pen_diag
                d = diag(Gc);
                d = sign(d) .* max(abs(d) - alpha_try*lambda2_eff, 0);
                Gc(1:n+1:end) = d;
            end

            % 活动集
            if ~isempty(A_masks)
                mask = A_masks{f};
                Gc(~mask) = 0;
                Gc = 0.5*(Gc+Gc');
            end

            % %%% FIX: SPD 保护必须对“每个频点”做完整投影（对称化+对角归一+特征值地板）
           Gc = project_spd_(Gc, 1e-6, true);

            Gamma_cand{f} = Gc;
        end

        % 目标与 Armijo 判据（带兜底：若仍出错则降步重试）
        try
            [loglik_val, smooth_val, l1_val, aux_terms] = module5_objective_terms(Gamma_cand, Sigma, aux, params_obj);
            spatial_val = 0; if isfield(aux_terms,'spatial_term'), spatial_val = aux_terms.spatial_term; end
            obj_trial = loglik_val + smooth_val + l1_val + spatial_val;
        catch ME
            if contains(ME.message,'objective_terms')
                % %%% NEW: 目标评估失败 → 认为步长过大（或谱地板不足），回溯
                bt_cnt  = bt_cnt + 1; bt_cum = bt_cum + 1;
                alpha_try = max(alpha_min, bt_beta * alpha_try);
                if bt_cnt >= bt_max_iter || alpha_try <= alpha_min
                    obj_trial = obj_val;   % 放弃本轮更新
                    break;
                else
                    continue;              % 继续回溯尝试
                end
            else
                rethrow(ME);
            end
        end

        if obj_trial <= obj_val - armijo_c1 * alpha_try * (grad_norm^2)
            accept  = true;
            obj_new = obj_trial;
            Gamma   = Gamma_cand;
        else
            bt_cnt  = bt_cnt + 1;
            bt_cum  = bt_cum + 1;
            alpha_try = max(alpha_min, bt_beta * alpha_try);
            if bt_cnt >= bt_max_iter || alpha_try <= alpha_min
                % 本轮放弃进一步回溯；接受“不更新”以继续到退火逻辑
                obj_new = obj_val;
                break;
            end
        end
    end

    % 4) 回溯累计过多 -> “优先降 λ2”（退火），并温和放宽 α
    if bt_cum >= anneal_patience && lambda2_eff > lambda2_min
        old_l2 = lambda2_eff;
        lambda2_eff = max(lambda2_min, anneal_factor * lambda2_eff);
        aux.lambda2 = lambda2_eff;
        if verbose
            fprintf('[module5][anneal] too many backtracks: lambda2 %.3g -> %.3g; relax alpha\n', old_l2, lambda2_eff);
        end
        bt_cum   = 0;
        alpha_try = min(alpha_max, max(alpha_try, 1.2*alpha));  % 适度放宽
    end

    % 5) 接受/记录
    if accept
        prev_Gamma = Gamma; prev_grad = Ggrad;
        obj_val = obj_new; best_obj = min(best_obj, obj_val);
        alpha   = alpha_try;
    else
        % 未改善的迭代：保持 Γ/obj，不改变 best_obj，仅记录步长尝试
        alpha = alpha_try;
    end

    % histories
    objective_history(it)     = obj_val;
    gradient_norm_history(it) = grad_norm;
    step_size_history(it)     = alpha;
    backtrack_counts(it)      = bt_cnt;

    % logging to console
    if verbose && mod(it, getf(diagopt,'print_every',1))==0
        fprintf('[module5] it=%4d  obj=%.6e  (loglik=%.6e, smooth=%.6e, l1=%.6e, spatial=%.6e)  alpha=%.3g  ||grad||=%.1f\n', ...
            it, obj_val, loglik_val, smooth_val, l1_val, spatial_val, alpha, gradient_norm_history(it));
    end

    % logging to CSV
    if csv_fid~=-1
        fprintf(csv_fid, 'STEP,%d,%.10g,%.10g,%.10g,%.10g,%.10g,%.6g,%.6g,%.6g,%.8g\n', ...
            it, obj_val, loglik_val, smooth_val, l1_val, spatial_val, ...
            getf(params,'obj_improve_tol',1e-6), alpha, gradient_norm_history(it), lambda2_eff);
    end

    % live plot（可选）
    if viz.enable && mod(it, max(1,round(viz.plot_every)))==0
        f_view = max(1, min(F, round(getf(viz,'f_view',1))));
        GTf = [];
        if ~isempty(GT_cells), GTf = GT_cells{f_view}; end
        viz_state = update_live_plot_(viz_state, it, ...
            objective_history(1:it), gradient_norm_history(1:it), ...
            loglik_val, smooth_val, l1_val, spatial_val, ...
            Gamma{f_view}, GTf, viz.value_mode);
    end

    % 简单早停：步长到下限且目标近似平台
    if it>=5
        recent = objective_history(max(1,it-4):it);
        if alpha<=alpha_min && max(abs(diff(recent))) < getf(params,'obj_improve_tol',1e-6)
            if verbose, fprintf('[module5] early stop: alpha near alpha_min and objective plateau.\n'); end
            break;
        end
    end
end

% ------------ Pack results ------------
Gamma_cells = Gamma;
results = struct('best_objective', best_obj, ...
                 'objective_history',       objective_history, ...
                 'gradient_norm_history',   gradient_norm_history, ...
                 'step_size_history',       step_size_history, ...
                 'backtrack_counts',        backtrack_counts, ...
                 'active_set_changes',      active_set_changes, ...
                 'csv_header_written',      csv_header_written, ...
                 'csv_config_written',      csv_config_written);

if csv_fid~=-1, fclose(csv_fid); end
end

% ---------------- helpers ----------------
function v = getf(S, name, def)
    if isfield(S,name) && ~isempty(S.(name)), v = S.(name); else, v = def; end
end
function S = set_defaults(S, D)
    ff = fieldnames(D);
    for i=1:numel(ff)
        k = ff{i};
        if ~isfield(S,k) || isempty(S.(k)), S.(k) = D.(k); end
    end
end
function gn = aggregate_grad_norm_(Gcells)
    s = 0;
    for i=1:numel(Gcells), Gi = Gcells{i}; s = s + sum(abs(Gi(:)).^2); end
    gn = sqrt(s);
end

function C = coerce_cov_cell_local(X, F_hint)
    if isa(X,'cell'), C = X(:); return; end
    if isnumeric(X) && ndims(X)==3 && size(X,1)==size(X,2)
        F = size(X,3); C = cell(F,1); for f=1:F, C{f}=X(:,:,f); end; return;
    end
    if isnumeric(X) && ismatrix(X) && size(X,1)==size(X,2)
        F = F_hint; C = repmat({X}, F, 1); return;
    end
    error('coerce_cov_cell_local:unsupported','Expect cell{F,1} | p×p×F | single p×p.');
end

function [a1, a2] = bb_stepsize_(G, Gprev, dG, dGprev)
% BB1/BB2（聚合所有频率的 Fro 内积）
    s2 = 0; sy = 0; y2 = 0;
    F_ = numel(G);
    for f_=1:F_
        S = G{f_}  - Gprev{f_};
        Y = dG{f_} - dGprev{f_};
        s2 = s2 + sum(real(S(:).*conj(S(:))));
        sy = sy + sum(real(S(:).*conj(Y(:))));
        y2 = y2 + sum(real(Y(:).*conj(Y(:))));
    end
    a1 = s2 / max(1e-16, sy);     % BB1
    a2 = max(1e-16, sy) / max(1e-16, y2); % BB2（备用）
end

% --------- live plot subroutines（与你原版一致，占位实现）---------
function S = init_live_plot_(F, n, viz, GT_cells)
    S = struct();
    if ~viz.enable, return; end
    S.fig = figure('Name','Module5 Live','NumberTitle','off'); clf(S.fig);
    S.tlo = tiledlayout(S.fig,2,3,'TileSpacing','compact','Padding','compact');
    S.ax_conv   = nexttile(S.tlo,1);
    S.ax_comps  = nexttile(S.tlo,2);
    S.ax_legend = nexttile(S.tlo,3); axis(S.ax_legend,'off');
    S.ax_gt     = nexttile(S.tlo,4);
    S.ax_est    = nexttile(S.tlo,5);
    S.ax_diff   = nexttile(S.tlo,6);
    title(S.ax_conv,  'Convergence');
    title(S.ax_comps, 'Objective components');
    title(S.ax_gt,    'GT (whitened/source-as-is)');
    title(S.ax_est,   'Estimate \Gamma (live)');
    title(S.ax_diff,  '|Est| - |GT|');
    axes(S.ax_legend); cla(S.ax_legend);
    txt = {sprintf('F=%d, f_{view}=%d', F, viz.f_view), ...
           sprintf('value_mode=%s', viz.value_mode), ...
           datestr(now)};
    text(0.0,1.0,strjoin(txt, '\n'),'Units','normalized','VerticalAlignment','top');
    S.im_gt   = imagesc(S.ax_gt, zeros(n)); colorbar(S.ax_gt); axis(S.ax_gt,'square');
    S.im_est  = imagesc(S.ax_est, zeros(n)); colorbar(S.ax_est); axis(S.ax_est,'square');
    S.im_diff = imagesc(S.ax_diff, zeros(n)); colorbar(S.ax_diff); axis(S.ax_diff,'square');
    S.GT_cache = []; if ~isempty(GT_cells), S.GT_cache = GT_cells; end
end

function S = update_live_plot_(S, it, obj_hist, grad_hist, ...
    loglik_val, smooth_val, l1_val, spatial_val, G_view, GT_view, value_mode)
    if isempty(S), return; end
    axes(S.ax_conv); cla(S.ax_conv); hold(S.ax_conv,'on');
    yyaxis(S.ax_conv,'left');
    base = min(obj_hist(1:it));
    semilogy(S.ax_conv, max(obj_hist(1:it)-base, eps), '-o','MarkerSize',3);
    ylabel(S.ax_conv,'Objective (shifted, log)');
    yyaxis(S.ax_conv,'right');
    semilogy(S.ax_conv, max(grad_hist(1:it), eps), '-s','MarkerSize',3);
    ylabel(S.ax_conv,'||grad||_F (log)');
    xlabel(S.ax_conv,'iter'); grid(S.ax_conv,'on');

    axes(S.ax_comps); cla(S.ax_comps); hold(S.ax_comps,'on');
    plot(S.ax_comps, 1:it, repmat(loglik_val,1,it), '-', 'DisplayName','loglik');
    plot(S.ax_comps, 1:it, repmat(smooth_val,1,it), '-', 'DisplayName','smooth');
    plot(S.ax_comps, 1:it, repmat(l1_val,1,it),     '-', 'DisplayName','l1');
    plot(S.ax_comps, 1:it, repmat(spatial_val,1,it),'-', 'DisplayName','spatial');
    legend(S.ax_comps,'Location','northeast'); grid(S.ax_comps,'on');
    xlabel('iter'); ylabel('value');

    switch lower(string(value_mode))
        case "abs",   V = abs(G_view);
        case "real",  V = real(G_view);
        case "imag",  V = imag(G_view);
        otherwise,    V = abs(G_view);
    end
    set(S.im_est,  'CData', V);
    if ~isempty(GT_view)
        switch lower(string(value_mode))
            case "abs",   GTV = abs(GT_view);
            case "real",  GTV = real(GT_view);
            case "imag",  GTV = imag(GT_view);
            otherwise,    GTV = abs(GT_view);
        end
        set(S.im_gt,   'CData', GTV);
        set(S.im_diff, 'CData', abs(V) - abs(GTV));
        title(S.ax_gt, 'GT (whitened/source-as-is)');
    else
        set(S.im_gt,   'CData', zeros(size(V)));
        set(S.im_diff, 'CData', zeros(size(V)));
        title(S.ax_gt, 'GT (N/A)');
    end
    drawnow limitrate;
end

% %%% NEW: 统一的 SPD 投影（对称化 + 可选对角归一 + 特征值地板）
function G = project_spd_(G, eig_floor, norm_diag)
    G = 0.5*(G+G');                     % Hermitian
    if norm_diag
        d = real(diag(G));
        d = max(d, 1e-12);
        D = diag(1./sqrt(d));
        G = D*G*D;                      % 使 diag ≈ 1
        G = 0.5*(G+G');
    end
    [U,S] = eig(0.5*(G+G'));
    s = max(diag(S), eig_floor);
    G = U*diag(s)*U';
    G = 0.5*(G+G');
end
