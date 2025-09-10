function [Gamma_cells, results] = module5_proximal_main(input_data, params)
% MODULE5_PROXIMAL_MAIN
% Proximal gradient solver for precision matrices {Gamma_f} in the whitened domain.
% Backward-compatible; diagnostics/logging improvements:
%   - CSV header includes 'lambda2_effective'
%   - CONFIG row logs lambda2_suggested/effective
%   - (NEW) penalize_diagonal switch for L1, default false
%   - (NEW) return histories: objective_history, gradient_norm_history,
%           step_size_history, backtrack_counts, active_set_changes (if any)
%   - (NEW) optional live plot during optimization (diag.live_plot.*)

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

% params (keep legacy names)
lambda1     = getf(params,'lambda1',0);
lambda2_eff = getf(params,'lambda2',0);               % effective value (REQUIRED)
lambda2_sug = getf(params,'lambda2_suggested',[]);    % optional provenance
pen_diag    = getf(params,'penalize_diagonal',false); % (NEW) L1 on diagonal?

alpha0   = getf(params,'alpha0',1e-3);
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
alpha = alpha0;

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
active_set_changes      = [];  % optional; keep empty unless you implement AS updates here

% Live figure state (optional)
GT_cells = [];
if viz.enable && ~isempty(viz.ground_truth_precision)
    % Convert GT to cells
    GT_cells = coerce_cov_cell_local(viz.ground_truth_precision, F);
    switch lower(string(viz.ground_truth_domain))
        case "source"
            if isempty(D_whiten)
                warning('ground_truth_domain=source but input_data.whitening_matrices not provided; using given GT as-is for visualization.');
            else
                % whiten: Γ̃_GT = D * Γ_GT * D
                for f = 1:F
                    Df = D_whiten{f};
                    GT_cells{f} = Df * GT_cells{f} * Df;
                end
            end
        case "whitened"
            % do nothing
        otherwise
            warning('Unknown ground_truth_domain=%s; using GT as-is.', string(viz.ground_truth_domain));
    end
end
viz_state = [];
if viz.enable
    viz_state = init_live_plot_(F, n, viz, GT_cells);
end

% ------------ Main loop ------------
best_obj = inf;
for it = 1:max_iter
    % gradient of smooth part
    input = struct();
    input.precision_matrices   = Gamma;        % {F×1} cell
    input.whitened_covariances = Sigma;        % {F×1} cell
    input.smoothing_kernel     = K;            % single source
    input.weight_matrix        = W;

    gres = module4_objective_gradient_main(input, grad_params);
    Ggrad = gres.smooth_gradients;

    % gradient step
    Gamma_tmp = cell(F,1);
    for f=1:F, Gamma_tmp{f} = Gamma{f} - alpha * Ggrad{f}; end

    % proximal L1 (off-diagonals; optionally diagonals)
    Gamma_next = Gamma_tmp;
    for f=1:F
        G = Gamma_tmp{f};

        % --- off-diagonal shrink (legacy) ---
        U = triu(G, 1);
        U = sign(U) .* max(abs(U) - alpha*lambda2_eff, 0);
        G = diag(diag(G)) + U + U';

        % --- (NEW) diagonal shrink if requested ---
        if pen_diag
            d = diag(G);
            d = sign(d) .* max(abs(d) - alpha*lambda2_eff, 0);
            G(1:n+1:end) = d;
        end

        % optional: mask by active set
        if ~isempty(A_masks)
            mask = A_masks{f};
            G(~mask) = 0;
            G = 0.5*(G+G');  % keep symmetric
        end
        % SPD safeguard (tiny diagonal loading if necessary)
        [~,flag] = chol(0.5*(G+G'),'lower');
        if flag~=0
            G = 0.5*(G+G') + 1e-12*eye(n);
        end
        Gamma_next{f} = G;
    end

    % objective value (loglik + cross-freq smooth + l1 + spatial(optional))
    [loglik_val, smooth_val, l1_val, aux_terms] = module5_objective_terms(Gamma_next, Sigma, aux, params_obj);
    spatial_val = 0;
    if isfield(aux_terms,'spatial_term'), spatial_val = aux_terms.spatial_term; end
    obj_val = loglik_val + smooth_val + l1_val + spatial_val;

    % acceptance & step control (simple Armijo-like without explicit backtrack)
    accepted = (obj_val < best_obj);
    if accepted
        best_obj = obj_val;
        Gamma = Gamma_next;   % accept
        alpha = min(getf(params,'alpha_max',2e-3), getf(params,'alpha_up',1.05)*alpha);
    else
        alpha = max(1e-8, getf(params,'alpha_down',0.7)*alpha);
    end

    % --- histories ---
    objective_history(it)     = obj_val;
    gradient_norm_history(it) = aggregate_grad_norm_(Ggrad);
    step_size_history(it)     = alpha;
    backtrack_counts(it)      = 0;   % no explicit backtracking in this variant

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

    % live plot
    if viz.enable && mod(it, max(1,round(viz.plot_every)))==0
        f_view = max(1, min(F, round(getf(viz,'f_view',1))));
        GTf = [];
        if ~isempty(GT_cells), GTf = GT_cells{f_view}; end
        viz_state = update_live_plot_(viz_state, it, ...
            objective_history(1:it), gradient_norm_history(1:it), ...
            loglik_val, smooth_val, l1_val, spatial_val, ...
            Gamma{f_view}, GTf, viz.value_mode);
    end
end

% ------------ Pack results ------------
Gamma_cells = Gamma;
% Trim histories to actual length (max_iter)
results = struct('best_objective', best_obj, ...
                 'objective_history',       objective_history, ...
                 'gradient_norm_history',   gradient_norm_history, ...
                 'step_size_history',       step_size_history, ...
                 'backtrack_counts',        backtrack_counts, ...
                 'active_set_changes',      active_set_changes, ...
                 'csv_header_written',      csv_header_written, ...
                 'csv_config_written',      csv_config_written);

% close CSV
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

% --------- live plot subroutines ---------
function S = init_live_plot_(F, n, viz, GT_cells)
    % Create tiled layout figure
    S = struct();
    S.fig = figure('Name','Module5 Live','NumberTitle','off'); clf(S.fig);
    S.tlo = tiledlayout(S.fig,2,3,'TileSpacing','compact','Padding','compact');

    % Axes
    S.ax_conv   = nexttile(S.tlo,1);
    S.ax_comps  = nexttile(S.tlo,2);
    S.ax_legend = nexttile(S.tlo,3); axis(S.ax_legend,'off');

    S.ax_gt     = nexttile(S.tlo,4);
    S.ax_est    = nexttile(S.tlo,5);
    S.ax_diff   = nexttile(S.tlo,6);

    % Static titles
    title(S.ax_conv,  'Convergence');
    title(S.ax_comps, 'Objective components');
    title(S.ax_gt,    'GT (whitened/source-as-is)');
    title(S.ax_est,   'Estimate \Gamma (live)');
    title(S.ax_diff,  '|Est| - |GT|');

    % Legend panel
    axes(S.ax_legend); cla(S.ax_legend);
    txt = {sprintf('F=%d, f_{view}=%d', F, viz.f_view), ...
           sprintf('value_mode=%s', viz.value_mode), ...
           datestr(now)};
    text(0.0,1.0,strjoin(txt, '\n'),'Units','normalized','VerticalAlignment','top');

    % Initialize images
    S.im_gt   = imagesc(S.ax_gt, zeros(n)); colorbar(S.ax_gt); axis(S.ax_gt,'square');
    S.im_est  = imagesc(S.ax_est, zeros(n)); colorbar(S.ax_est); axis(S.ax_est,'square');
    S.im_diff = imagesc(S.ax_diff, zeros(n)); colorbar(S.ax_diff); axis(S.ax_diff,'square');

    % Cache GT for color scaling (if provided)
    S.GT_cache = [];
    if ~isempty(GT_cells), S.GT_cache = GT_cells; end
end

function S = update_live_plot_(S, it, obj_hist, grad_hist, ...
    loglik_val, smooth_val, l1_val, spatial_val, G_view, GT_view, value_mode)

    % --- (1) Convergence ---
    axes(S.ax_conv); cla(S.ax_conv); hold(S.ax_conv,'on');
    yyaxis(S.ax_conv,'left');
    base = min(obj_hist(1:it));
    semilogy(S.ax_conv, max(obj_hist(1:it)-base, eps), '-o','MarkerSize',3);
    ylabel(S.ax_conv,'Objective (shifted, log)');
    yyaxis(S.ax_conv,'right');
    semilogy(S.ax_conv, max(grad_hist(1:it), eps), '-s','MarkerSize',3);
    ylabel(S.ax_conv,'||grad||_F (log)');
    xlabel(S.ax_conv,'iter'); grid(S.ax_conv,'on');

    % --- (2) Components ---
    axes(S.ax_comps); cla(S.ax_comps); hold(S.ax_comps,'on');
    plot(S.ax_comps, 1:it, repmat(loglik_val,1,it), '-', 'DisplayName','loglik');      % last-val trace
    plot(S.ax_comps, 1:it, repmat(smooth_val,1,it), '-', 'DisplayName','smooth');
    plot(S.ax_comps, 1:it, repmat(l1_val,1,it),     '-', 'DisplayName','l1');
    plot(S.ax_comps, 1:it, repmat(spatial_val,1,it),'-', 'DisplayName','spatial');
    legend(S.ax_comps,'Location','northeast'); grid(S.ax_comps,'on');
    xlabel(S.ax_comps,'iter'); ylabel(S.ax_comps,'value');

    % --- (3) Heatmaps ---
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
