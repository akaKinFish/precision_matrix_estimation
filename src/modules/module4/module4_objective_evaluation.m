function [objective_values, computation_stats] = module4_objective_evaluation(input_data, evaluation_params)
% MODULE4_OBJECTIVE_EVALUATION - Evaluate complete objective function (revised)
%
% F(Γ) = sum_w[-logdet Γ_w + tr(Σ_w Γ_w)] 
%        + (λ1/2) * sum_{w,w'} k_{w,w'} ||Γ_w - Γ_{w'}||^2_W
%        + λ2 * sum_w ||Γ_w||_L1   (diag excluded by default)
%
% The smoothing term also supports a Laplacian form which is algebraically
% equivalent:  S = λ1 * sum_{w,w'} L_{w,w'} <Γ_w, W Γ_{w'}>,
% where L = D - K, D = diag(K 1).  (No 1/2 in this form.)
%
% Options (evaluation_params) new defaults included:
%   .lambda1                (double, >=0) default 0.01
%   .lambda2                (double, >=0) default 0.01
%   .penalize_diagonal      (logical)     default false
%   .penalize_mask          (pxp logical or cell{F} of logical) default []
%   .chol_tolerance         (double)      default 1e-12
%   .use_graph_laplacian    (logical)     default true
%   .kernel_zero_tol        (double)      default 1e-12
%   .strict_input_validation(logical)     default true
%   .allow_eig_fallback     (logical)     default false (diagnostic only)
%   .verbose                (logical)     default false
%
% 新增 options:
%   .lambda3                     (double, >=0) default 0
%   .spatial_graph_matrix        (pxp)
%   .spatial_graph_is_laplacian  (logical) default true
%   .spatial_weight_mode         ('node'|'hadamard') default 'node'

% Stats fields unified:
%   computation_time (+ alias total_computation_time),
%   component_times.logdet/trace/smoothing/l1,
%   cholesky_failures, parameters_used (echo of params)

% ==================== Input Validation ====================
if nargin < 1
    error('module4_objective_evaluation:insufficient_input', ...
          'At least input_data is required');
end
if nargin < 2, evaluation_params = struct(); end

required_fields = {'precision_matrices','whitened_covariances','smoothing_kernel','weight_matrix'};
for i = 1:numel(required_fields)
    if ~isfield(input_data, required_fields{i})
        error('module4_objective_evaluation:missing_field', ...
             'Required field "%s" not found in input_data', required_fields{i});
    end
end

Gammas = input_data.precision_matrices;
Sigmas = input_data.whitened_covariances;
K      = input_data.smoothing_kernel;
W      = input_data.weight_matrix;

if ~iscell(Gammas) || ~iscell(Sigmas) || isempty(Gammas) || isempty(Sigmas)
    error('module4_objective_evaluation:invalid_matrices', ...
          'precision_matrices and whitened_covariances must be non-empty cell arrays');
end
F = numel(Gammas);
p = size(Gammas{1},1);
if numel(Sigmas) ~= F
    error('module4_objective_evaluation:length_mismatch', ...
          'precision_matrices and whitened_covariances must have same length');
end
for f = 1:F
    if ~isequal(size(Gammas{f}),[p,p]) || ~isequal(size(Sigmas{f}),[p,p])
        error('module4_objective_evaluation:dimension_mismatch', ...
              'All matrices must be %dx%d', p, p);
    end
end
if ~isequal(size(K),[F,F]) || ~isequal(size(W),[p,p])
    error('module4_objective_evaluation:kernel_weight_mismatch', ...
          'kernel must be %dx%d and weight must be %dx%d', F, F, p, p);
end

% ==================== Parameters & Defaults ====================
defaults = struct();
defaults.lambda1                 = 0.01;
defaults.lambda2                 = 0.01;
defaults.penalize_diagonal       = false;
defaults.penalize_mask           = [];
defaults.chol_tolerance          = 1e-12;
defaults.use_graph_laplacian     = true;
defaults.kernel_zero_tol         = 1e-12;
defaults.strict_input_validation = true;
defaults.allow_eig_fallback      = false;   % diagnostics only
defaults.verbose                 = false;

fn = fieldnames(defaults);
for i = 1:numel(fn)
    if ~isfield(evaluation_params, fn{i})
        evaluation_params.(fn{i}) = defaults.(fn{i});
    end
end
if evaluation_params.lambda1 < 0 || evaluation_params.lambda2 < 0
    error('module4_objective_evaluation:invalid_lambda', ...
          'lambda1 and lambda2 must be non-negative');
end

% ==================== Initialize Results ====================
objective_values = struct();
objective_values.total_objective   = 0;
objective_values.smooth_objective  = 0;
objective_values.logdet_terms      = 0;
objective_values.trace_terms       = 0;
objective_values.smoothing_penalty = 0;
objective_values.l1_penalty        = 0;
objective_values.individual_logdet = zeros(F,1);
objective_values.individual_trace  = zeros(F,1);

computation_stats = struct();
computation_stats.computation_time  = 0;
computation_stats.component_times   = struct();
computation_stats.cholesky_failures = 0;
computation_stats.parameters_used   = evaluation_params;

t_total = tic;

% ==================== Preprocess K and W ====================
% Symmetrize K, zero diagonal, drop tiny entries
K = (K + K')/2;
K(1:F+1:end) = 0;
K(abs(K) < evaluation_params.kernel_zero_tol) = 0;

% Symmetrize W (PSD is recommended but not strictly enforced here)
W = (W + W')/2;

% ==================== 1) Log-determinant terms ====================
if evaluation_params.verbose, fprintf('Computing log-determinant terms... '); end
t_log = tic;

for f = 1:F
    Gf = (Gammas{f} + Gammas{f}')/2;  % enforce exact Hermitian before chol

    % Regularized Cholesky with retries
    tol = evaluation_params.chol_tolerance;
    reg = max(tol, tol*mean(diag(Gf)));
    success = false;
    for attempt = 0:4
        [R, flag] = chol(Gf + (attempt>0)*reg*(10^(attempt-1))*eye(p), 'upper');
        if flag == 0
            if attempt>0 && evaluation_params.verbose
                fprintf(' (chol regularized f=%d, attempt=%d)', f, attempt);
            end
            logdet_val = 2*sum(log(diag(R)));
            objective_values.individual_logdet(f) = -logdet_val;
            success = true;
            break;
        end
    end
    if ~success
        computation_stats.cholesky_failures = computation_stats.cholesky_failures + 1;
        if evaluation_params.allow_eig_fallback
            % diagnostic only
            ev = eig(Gf);
            ev = max(real(ev), evaluation_params.chol_tolerance);
            objective_values.individual_logdet(f) = -sum(log(ev));
        else
            error('module4_objective_evaluation:spd_violation', ...
                 'Cholesky failed for frequency %d; matrix not SPD after regularization.', f);
        end
    end
end
objective_values.logdet_terms = sum(objective_values.individual_logdet);
computation_stats.component_times.logdet = toc(t_log);
if evaluation_params.verbose
    fprintf('completed (%.3fs)\n', computation_stats.component_times.logdet);
end

% ==================== 2) Trace terms ====================
if evaluation_params.verbose, fprintf('Computing trace terms... '); end
t_tr = tic;

for f = 1:F
    Sf = (Sigmas{f} + Sigmas{f}')/2;        % use Hermitian part
    objective_values.individual_trace(f) = real(trace(Sf * Gammas{f}));
end
objective_values.trace_terms = sum(objective_values.individual_trace);
computation_stats.component_times.trace = toc(t_tr);
if evaluation_params.verbose
    fprintf('completed (%.3fs)\n', computation_stats.component_times.trace);
end

% ==================== 3) Smoothing penalty ====================
if evaluation_params.verbose, fprintf('Computing smoothing penalty... '); end
t_sm = tic;

if evaluation_params.lambda1 > 0
    if evaluation_params.use_graph_laplacian
        % Laplacian form: S = λ1 * sum_{w,w'} L_{w,w'} <Γ_w, W Γ_{w'}>
        d = sum(K,2);
        L = diag(d) - K;

        % Precompute WG_w = W * Γ_w
        WG = cell(F,1);
        for w = 1:F
            WG{w} = W * Gammas{w};
        end

        S = 0;
        for w = 1:F
            for wp = 1:F
                if L(w,wp) ~= 0
                    S = S + L(w,wp) * real(trace(Gammas{w}' * WG{wp}));
                end
            end
        end
        objective_values.smoothing_penalty = evaluation_params.lambda1 * S; % no 1/2 here
    else
        % Direct double sum with 1/2 to avoid double counting pairs
        S = 0;
        for w = 1:F
            Gw = Gammas{w};
            for wp = 1:F
                kwp = K(w,wp);
                if kwp ~= 0
                    D = Gw - Gammas{wp};
                    S = S + kwp * real(trace(D' * W * D));
                end
            end
        end
        objective_values.smoothing_penalty = 0.5 * evaluation_params.lambda1 * S; % <-- 1/2 fix
    end
else
    objective_values.smoothing_penalty = 0;
end
computation_stats.component_times.smoothing = toc(t_sm);
if evaluation_params.verbose
    fprintf('completed (%.3fs)\n', computation_stats.component_times.smoothing);
end

% ==================== 4) L1 penalty (vectorized, mask-aware) ====================
if evaluation_params.verbose, fprintf('Computing L1 penalty... '); end
t_l1 = tic;

if evaluation_params.lambda2 > 0
    l1_total = 0;

    % Build per-frequency mask M_f (logical pxp)
    has_mask = ~isempty(evaluation_params.penalize_mask);
    if has_mask
        if iscell(evaluation_params.penalize_mask)
            assert(numel(evaluation_params.penalize_mask)==F, ...
                'penalize_mask cell must have length F');
        else
            % broadcast a single pxp mask to all frequencies
            assert(isequal(size(evaluation_params.penalize_mask),[p,p]), ...
                'penalize_mask must be pxp or cell{F} of pxp');
        end
    end

    for f = 1:F
        Gf = Gammas{f};
        if has_mask
            Mf = iscell(evaluation_params.penalize_mask) ...
               && ~isempty(evaluation_params.penalize_mask{f});
            if Mf
                M = logical(evaluation_params.penalize_mask{f});
            else
                M = logical(evaluation_params.penalize_mask);
            end
        else
            M = true(p); % start with all
        end
        if ~evaluation_params.penalize_diagonal
            M(1:p+1:end) = false; % drop diagonal
        end
        if any(M(:))
            l1_total = l1_total + sum(abs(Gf(M)));
        end
    end

    objective_values.l1_penalty = evaluation_params.lambda2 * l1_total;
else
    objective_values.l1_penalty = 0;
end
computation_stats.component_times.l1 = toc(t_l1);
if evaluation_params.verbose
    fprintf('completed (%.3fs)\n', computation_stats.component_times.l1);
end

% ==================== Spatial (λ3) single-frequency penalty ====================
objective_values.spatial_penalty = 0;
if isfield(evaluation_params,'lambda3') && evaluation_params.lambda3 > 0 && ...
   isfield(evaluation_params,'spatial_graph_matrix') && ~isempty(evaluation_params.spatial_graph_matrix)
    try
        sp_out = module5_spatial_smoothing_singlefreq(Gammas, struct( ...
            'lambda3', evaluation_params.lambda3, ...
            'spatial_graph_matrix', evaluation_params.spatial_graph_matrix, ...
            'spatial_graph_is_laplacian', getfield(evaluation_params,'spatial_graph_is_laplacian',true), ...
            'spatial_weight_mode', getfield(evaluation_params,'spatial_weight_mode','node'), ...
            'return_gradient', false, ...
            'validate_inputs', true, ...
            'enforce_hermitian_grad', false ...
        ));
        objective_values.spatial_penalty = sp_out.term;
    catch ME
        warning('module4_objective_evaluation:spatial_eval_failed','%s',ME.message);
        objective_values.spatial_penalty = 0;
    end
end


% ==================== Combine ====================
% ==================== Combine ====================
objective_values.smooth_objective = ...
    objective_values.logdet_terms + ...
    objective_values.trace_terms + ...
    objective_values.smoothing_penalty + ...
    objective_values.spatial_penalty;

objective_values.total_objective = ...
    objective_values.smooth_objective + objective_values.l1_penalty;


computation_stats.computation_time = toc(t_total);
% alias for backward compatibility
computation_stats.total_computation_time = computation_stats.computation_time;

if evaluation_params.verbose
    fprintf('--------------------------------------------\n');
    fprintf('Objective Function Breakdown:\n');
    fprintf('  Log-determinant:  %12.6f\n', objective_values.logdet_terms);
    fprintf('  Trace terms:      %12.6f\n', objective_values.trace_terms);
    fprintf('  Smoothing (λ1):   %12.6f\n', objective_values.smoothing_penalty);
    fprintf('  Spatial (λ3):     %12.6f\n', objective_values.spatial_penalty);

    fprintf('  L1 penalty (λ2):  %12.6f\n', objective_values.l1_penalty);
    fprintf('  ────────────────────────────────\n');
    fprintf('  Smooth objective: %12.6f\n', objective_values.smooth_objective);
    fprintf('  Total objective:  %12.6f\n', objective_values.total_objective);
    fprintf('--------------------------------------------\n');
    fprintf('Computation Summary:\n');
    fprintf('  Total time: %.3fs | Cholesky failures: %d/%d\n', ...
        computation_stats.computation_time, computation_stats.cholesky_failures, F);
    fprintf('============================================\n');
end

end
