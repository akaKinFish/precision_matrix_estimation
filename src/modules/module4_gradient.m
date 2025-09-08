function gradient_results = module4_gradient(input_data, gradient_params)
% MODULE4_GRADIENT - Main interface for Module 4 gradient computation (wrapper, BC-fixed)
%
% 作用：
%   - 与 module4_objective_gradient_main 对接；
%   - 兼容多套字段命名（别名映射），只在 wrapper 内做胶水，不改 main；
%   - 保证输出含 .smooth_gradients 与 .success 字段（向后兼容）。
%
% 期待 main 的输入字段：
%   .precision_matrices    (cell{F,1})  Γ_ω
%   .whitened_covariances  (cell{F,1})  Σ̃_ω
%   .smoothing_kernel      (F×F)
%   .weight_matrix         (p×p)
%
% Pipeline 常见别名：
%   precision  : {'precision_matrices','Gamma_list','Gamma','Gamma_tilde','Gamma_current'}
%   covariances: {'whitened_covariances','Sigma_tilde','Sigma_whitened','Sigma_whitened_list'}
%   kernel     : {'smoothing_kernel','kernel_matrix','K'}
%   weight     : {'weight_matrix','W_gamma','W'}
%
% 注意：不做算法行为改变，只做字段名与形状规整。

if nargin < 2, gradient_params = struct(); end
if nargin < 1 || ~isstruct(input_data)
    error('module4_gradient:insufficient_input','input_data must be a struct');
end

% -------- 组装 core_input（别名规范化）--------
core_input = struct();

% 1) precision_matrices (cell{F,1}, p×p each)
core_input.precision_matrices = pick_and_normalize_cell(input_data, ...
    {'precision_matrices','Gamma_list','Gamma','Gamma_tilde','Gamma_current'}, ...
    'module4_gradient:missing_precision', 'precision matrices');

% 2) whitened_covariances (cell{F,1}, p×p each)
core_input.whitened_covariances = pick_and_normalize_cell(input_data, ...
    {'whitened_covariances','Sigma_tilde','Sigma_whitened','Sigma_whitened_list'}, ...
    'module4_gradient:missing_covariances', 'whitened covariances');

% 3) smoothing_kernel (F×F)
core_input.smoothing_kernel = pick_and_check_matrix(input_data, ...
    {'smoothing_kernel','kernel_matrix','K'}, ...
    size(core_input.whitened_covariances,1), ...
    'module4_gradient:missing_kernel', 'smoothing kernel');

% 4) weight_matrix (p×p)
p = size(core_input.precision_matrices{1},1);
core_input.weight_matrix = pick_and_check_matrix(input_data, ...
    {'weight_matrix','W_gamma','W'}, ...
    p, ...
    'module4_gradient:missing_weight', 'weight matrix');

% -------- 路由到 main --------
try
    gradient_results = module4_objective_gradient_main(core_input, gradient_params);
catch ME
    error('module4_gradient:computation_failed','Module 4 gradient computation failed: %s', ME.message);
end

% -------- 向后兼容输出字段 --------
if ~isfield(gradient_results,'smooth_gradients')
    % 若 main 用了别名如 .gradients，则映射回来
    if isfield(gradient_results,'gradients')
        gradient_results.smooth_gradients = gradient_results.gradients;
    else
        error('module4_gradient:missing_output','Expected field "smooth_gradients" not found in results.');
    end
end
if ~isfield(gradient_results,'success') || isempty(gradient_results.success)
    gradient_results.success = true;
end

end

% ================= 工具函数（仅限本文件） =================

function cell_out = pick_and_normalize_cell(S, names, err_id, human_name)
% 在 S 中按 names 顺序查找字段，接受：
%   - cell{F,1}  (每项 p×p)
%   - 单个矩阵 p×p          -> 复制成 F×1（若能从其它域推断 F，否则报错）
%   - 3D 数组 p×p×F         -> 切成 cell{F,1}
% 如果无法确定 F，尝试从其它已存在 cell 字段推断。

val = [];
for k = 1:numel(names)
    if isfield(S, names{k}) && ~isempty(S.(names{k}))
        val = S.(names{k});
        break;
    end
end
if isempty(val)
    error(err_id, 'Missing %s (tried: %s)', human_name, strjoin(names, ', '));
end

% cell 形式
if iscell(val)
    cell_out = val(:);
    return;
end

% 3D 数组 -> cell
if isnumeric(val) && ndims(val)==3
    F = size(val,3);
    cell_out = cell(F,1);
    for f=1:F
        cell_out{f} = val(:,:,f);
    end
    return;
end

% 单个 2D 矩阵：尝试推断 F
if isnumeric(val) && ismatrix(val) && size(val,1)==size(val,2)
    F = infer_F_from_struct(S);
    if isempty(F)
        error(err_id, ['A single %s matrix was provided but frequency count F could not be inferred. ', ...
                       'Please pass a cell array or 3D stack.'], human_name);
    end
    cell_out = repmat({val}, F, 1);
    return;
end

error(err_id, 'Unsupported format for %s.', human_name);
end

function F = infer_F_from_struct(S)
% 尝试从常见域推断 F（cell 的长度优先，其次 3D 第三维）
F = [];
candidates = {'precision_matrices','Gamma_list','Gamma','Gamma_tilde','Gamma_current', ...
              'whitened_covariances','Sigma_tilde','Sigma_whitened','Sigma_whitened_list'};
for k=1:numel(candidates)
    if isfield(S, candidates{k}) && ~isempty(S.(candidates{k}))
        v = S.(candidates{k});
        if iscell(v), F = numel(v); return; end
        if isnumeric(v) && ndims(v)==3, F = size(v,3); return; end
    end
end
end

function M = pick_and_check_matrix(S, names, dim, err_id, human_name)
% 取第一个命中的矩阵字段，并检查维度是 dim×dim
M = [];
for k = 1:numel(names)
    if isfield(S, names{k}) && ~isempty(S.(names{k}))
        M = S.(names{k});
        break;
    end
end
if isempty(M)
    error(err_id, 'Missing %s (tried: %s)', human_name, strjoin(names, ', '));
end
if ~isnumeric(M) || ~ismatrix(M) || any(size(M)~=[dim dim])
    error(err_id, '%s must be %dx%d, got %dx%d.', human_name, dim, dim, size(M,1), size(M,2));
end
end
