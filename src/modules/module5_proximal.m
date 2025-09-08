function [Gamma_final, proximal_results] = module5_proximal(input_data, proximal_params)
% MODULE5_PROXIMAL - Proximal solver wrapper (BC-fixed, alias mapping)
%
% 作用：
%   - 与 module5_proximal_main 对接；
%   - 只在 wrapper 里做字段别名映射与形状规范化（不改 main）；
%   - 兼容 Pipeline A 常见命名（Sigma_tilde / kernel_matrix / combined_active_mask 等）。
%
% main 期望输入字段：
%   .whitened_covariances   (cell{F,1})   Σ̃_ω
%   .initial_precision      (cell{F,1})   Γ_ω^0
%   .smoothing_kernel       (F×F)
%   .weight_matrix          (p×p)
%   .active_set_masks       (cell{F,1}) or []  (可选)
%
% Pipeline 常见别名：
%   whitened_covariances : {'whitened_covariances','Sigma_tilde','Sigma_whitened','Sigma_whitened_list'}
%   initial_precision    : {'initial_precision','initial_precision_matrices','Gamma_init','initial_Gamma','initial_precision_list'}
%   smoothing_kernel     : {'smoothing_kernel','kernel_matrix','K'}
%   weight_matrix        : {'weight_matrix','W_gamma','W'}
%   active_set_masks     : {'active_set_masks','combined_active_masks','combined_active_mask','edge_active_masks','edge_active_mask'}
%
% 注意：
%   - 若 active_set 以 3D 布尔 p×p×F 传入，会自动拆为 cell{F,1}；
%   - 若完全缺失 active_set_masks，将传空数组 [] 给 main（由 main 决定默认）。

if nargin < 2, proximal_params = struct(); end
if nargin < 1 || ~isstruct(input_data)
    error('module5_proximal:insufficient_input','input_data must be a struct');
end

core_input = struct();

% 1) whitened_covariances
core_input.whitened_covariances = pick_and_normalize_cell(input_data, ...
    {'whitened_covariances','Sigma_tilde','Sigma_whitened','Sigma_whitened_list'}, ...
    'module5_proximal:missing_covariances','whitened covariances');

F = numel(core_input.whitened_covariances);
p = size(core_input.whitened_covariances{1},1);

% 2) initial_precision
core_input.initial_precision = pick_and_normalize_cell_with_defaultF(input_data, ...
    {'initial_precision','initial_precision_matrices','Gamma_init','initial_Gamma','initial_precision_list'}, ...
    F, p, 'module5_proximal:missing_initial_precision','initial precision');

% 3) smoothing_kernel (F×F)
core_input.smoothing_kernel = pick_and_check_matrix(input_data, ...
    {'smoothing_kernel','kernel_matrix','K'}, ...
    F, 'module5_proximal:missing_kernel','smoothing kernel');

% 4) weight_matrix (p×p)
core_input.weight_matrix = pick_and_check_matrix(input_data, ...
    {'weight_matrix','W_gamma','W'}, ...
    p, 'module5_proximal:missing_weight','weight matrix');

% 5) active_set_masks（可选）
core_input.active_set_masks = pick_and_normalize_masks_optional(input_data, ...
    {'active_set_masks','combined_active_masks','combined_active_mask','edge_active_masks','edge_active_mask'}, ...
    F, p);

% ---- 路由到 main ----
[Gamma_final, proximal_results] = module5_proximal_main(core_input, proximal_params);

% ---- 向后兼容一些输出别名（不覆盖 main 的原字段）----
if isstruct(proximal_results)
    if ~isfield(proximal_results,'final_precision_matrices')
        proximal_results.final_precision_matrices = Gamma_final;
    end
    if ~isfield(proximal_results,'success') || isempty(proximal_results.success)
        proximal_results.success = true;
    end
end

end

% ================= 工具函数（仅限本文件） =================

function cell_out = pick_and_normalize_cell(S, names, err_id, human_name)
% 接受 cell{F,1} / 3D p×p×F / 单个 p×p（若无法推断 F 将报错）

val = [];
for k=1:numel(names)
    if isfield(S,names{k}) && ~isempty(S.(names{k}))
        val = S.(names{k}); break;
    end
end
if isempty(val)
    error(err_id, 'Missing %s (tried: %s)', human_name, strjoin(names, ', '));
end

if iscell(val)
    cell_out = val(:);
    return;
end

if isnumeric(val) && ndims(val)==3
    F = size(val,3);
    cell_out = cell(F,1);
    for f=1:F, cell_out{f} = val(:,:,f); end
    return;
end

if isnumeric(val) && ismatrix(val) && size(val,1)==size(val,2)
    % 单频 p×p：尝试从 S 其它域推断 F
    F = infer_F_from_struct(S);
    if isempty(F)
        error(err_id, ['A single %s matrix was provided but frequency count F could not be inferred. ', ...
                       'Please pass a cell array or 3D stack.'], human_name);
    end
    cell_out = repmat({val}, F, 1);
    return;
end

error(err_id,'Unsupported format for %s.', human_name);
end

function cell_out = pick_and_normalize_cell_with_defaultF(S, names, F, p, err_id, human_name)
% 与上类似，但若为单个 p×p 矩阵，直接复制为 F×1（已知 F）。

val = [];
for k=1:numel(names)
    if isfield(S,names{k}) && ~isempty(S.(names{k}))
        val = S.(names{k}); break;
    end
end
if isempty(val)
    error(err_id, 'Missing %s (tried: %s)', human_name, strjoin(names, ', '));
end

if iscell(val)
    if numel(val) ~= F
        error(err_id,'%s cell length mismatch: got %d, expected %d.', human_name, numel(val), F);
    end
    cell_out = val(:);
    return;
end

if isnumeric(val) && ndims(val)==3
    if size(val,3) ~= F
        error(err_id,'%s 3D stack freq mismatch: got %d, expected %d.', human_name, size(val,3), F);
    end
    cell_out = cell(F,1);
    for f=1:F, cell_out{f} = val(:,:,f); end
    return;
end

if isnumeric(val) && ismatrix(val) && size(val,1)==size(val,2)
    if size(val,1) ~= p
        error(err_id,'%s size mismatch: got %dx%d, expected %dx%d.', human_name, size(val,1), size(val,2), p, p);
    end
    cell_out = repmat({val}, F, 1);
    return;
end

error(err_id,'Unsupported format for %s.', human_name);
end

function M = pick_and_check_matrix(S, names, dim, err_id, human_name)
M = [];
for k=1:numel(names)
    if isfield(S,names{k}) && ~isempty(S.(names{k}))
        M = S.(names{k}); break;
    end
end
if isempty(M)
    error(err_id,'Missing %s (tried: %s)', human_name, strjoin(names, ', '));
end
if ~isnumeric(M) || ~ismatrix(M) || any(size(M)~=[dim dim])
    error(err_id,'%s must be %dx%d, got %dx%d.', human_name, dim, dim, size(M,1), size(M,2));
end
end

function F = infer_F_from_struct(S)
F = [];
candidates = {'whitened_covariances','Sigma_tilde','Sigma_whitened','Sigma_whitened_list', ...
              'initial_precision','initial_precision_matrices','Gamma_init','initial_Gamma','initial_precision_list'};
for k=1:numel(candidates)
    if isfield(S,candidates{k}) && ~isempty(S.(candidates{k}))
        v = S.(candidates{k});
        if iscell(v), F = numel(v); return; end
        if isnumeric(v) && ndims(v)==3, F = size(v,3); return; end
    end
end
end

function masks = pick_and_normalize_masks_optional(S, names, F, p)
% 返回 cell{F,1} 的逻辑矩阵；若完全缺失，返回 []
val = [];
for k=1:numel(names)
    if isfield(S,names{k}) && ~isempty(S.(names{k}))
        val = S.(names{k}); break;
    end
end
if isempty(val)
    masks = [];  % 交给 main 决定默认行为
    return;
end

% 已是 cell{F,1}
if iscell(val)
    if numel(val) ~= F
        error('module5_proximal:mask_length_mismatch', ...
              'active_set_masks cell length mismatch: got %d, expected %d.', numel(val), F);
    end
    masks = cell(F,1);
    for f=1:F
        A = logical(val{f});
        if ~ismatrix(A) || any(size(A)~=[p p])
            error('module5_proximal:mask_size_mismatch','Mask at f=%d must be %dx%d.', f, p, p);
        end
        % 轻微数值防御：强制对称
        masks{f} = (A | A'); 
    end
    return;
end

% 3D 逻辑/数值 -> cell
if isnumeric(val) || islogical(val)
    if ndims(val)==3 && size(val,1)==p && size(val,2)==p && size(val,3)==F
        masks = cell(F,1);
        for f=1:F
            A = logical(val(:,:,f));
            masks{f} = (A | A');
        end
        return;
    end
end

error('module5_proximal:unsupported_masks_format', ...
      'active set masks must be cell{F,1} of %dx%d or logical %dx%dx%d.', p,p,p,p,F);
end
