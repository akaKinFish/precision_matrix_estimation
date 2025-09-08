function estep_results = module2_estep(input_data, estep_params)
% MODULE2_ESTEP - E-step computation wrapper (fixed, backward-compatible)
%
% 该包装器仅做“胶水”和“兼容别名”：
%   - 调用 module2_estep_main(input_data, estep_params)
%   - 将主函数产出的字段映射到历史可能使用的命名：
%       transfer_functions                -> source_transfer_functions
%       residual_covariances              -> residual_empirical_covariances
%       posterior_source_covariances      -> posterior_covariances
%       initial_precision_matrices        -> initial_precision
%   - 保证 .success 字段存在（布尔）
%   - 不改动主函数的真实输出内容
%
% 注意：不对 main 的算法和数值做任何修改；如需理论或算法层面的调整，
% 应在 module2_estep_main 中进行（当前版本无需）。

% ==================== Input Validation ====================
if nargin < 1
    error('module2_estep:insufficient_input', 'At least input_data is required');
end
if nargin < 2
    estep_params = struct();
end

% ==================== Route to Main Implementation ====================
try
    core = module2_estep_main(input_data, estep_params);
catch ME
    % 使用稳定的外层错误标识，并保留原始异常信息
    ME2 = MException('module2_estep:computation_failed', ...
                     'Module 2 E-step computation failed: %s', ME.message);
    try
        ME2 = addCause(ME2, ME);
    catch
        % 某些旧版 MATLAB 无 addCause，则忽略
    end
    throw(ME2);
end

% ==================== Backward-compatible Field Aliases ====================
% 直接把 main 的结果作为基底
estep_results = core;

% 成功标志保证存在
if ~isfield(estep_results, 'success') || isempty(estep_results.success)
    estep_results.success = true;
end

% 1) source_transfer_functions（旧名） <- transfer_functions（主名）
if isfield(core, 'transfer_functions')
    estep_results.source_transfer_functions = core.transfer_functions;
elseif isfield(core, 'source_transfer_functions')
    % 反向也兼容一下：如果主函数将来换回旧名
    estep_results.transfer_functions = core.source_transfer_functions;
else
    % 两者都没有，给出清晰错误信息（保持 wrapper 行为可诊断）
    error('module2_estep:missing_output_field', ...
        'Neither "transfer_functions" nor "source_transfer_functions" found in main output.');
end

% 2) residual_empirical_covariances（旧名） <- residual_covariances（主名）
if isfield(core, 'residual_covariances') && ~isfield(estep_results, 'residual_empirical_covariances')
    estep_results.residual_empirical_covariances = core.residual_covariances;
end

% 3) posterior_covariances（旧名） <- posterior_source_covariances（主名）
if isfield(core, 'posterior_source_covariances') && ~isfield(estep_results, 'posterior_covariances')
    estep_results.posterior_covariances = core.posterior_source_covariances;
end

% 4) initial_precision（旧名） <- initial_precision_matrices（主名）
if isfield(core, 'initial_precision_matrices') && ~isfield(estep_results, 'initial_precision')
    estep_results.initial_precision = core.initial_precision_matrices;
end
if isfield(core,'source_second_moments') && ~isfield(estep_results,'Sjj_hat')
    estep_results.Sjj_hat = core.source_second_moments;
end
if isfield(core,'initial_precision_matrices') && ~isfield(estep_results,'Omega_init')
    estep_results.Omega_init = core.initial_precision_matrices;
end
% ==================== Final Sanity (Do NOT force-rename back) ====================
% 不再做“必须存在旧名”之类的强校验；只要主名或旧名其一存在即可。
% 这样既满足通过 wrapper 统一调用的要求，也对老代码保持兼容。

end
