function config_results = module6_hyperparameter_config(input_data, config_params)
% MODULE6_HYPERPARAMETER_CONFIG - Hyperparameter configuration wrapper function
%
% This thin wrapper calls module6_hyperparameters and preserves a stable API.

% ==================== Input Validation ====================
if nargin < 1
    error('module6_hyperparameter_config:insufficient_input', ...
          'At least input_data is required');
end
if nargin < 2
    config_params = struct();
end

% ==================== Route to Main Implementation ====================
try
    config_results = module6_hyperparameters(input_data, config_params);
catch ME_original
    % Re-throw with a consistent identifier and try to attach the original cause
    try
        ME = MException('module6_hyperparameter_config:computation_failed', ...
                        'Module 6 hyperparameter configuration failed.');
        % addCause is available in newer MATLAB versions
        ME = addCause(ME, ME_original);
        throw(ME);
    catch
        error('module6_hyperparameter_config:computation_failed', ...
              'Module 6 hyperparameter configuration failed: %s', ME_original.message);
    end
end

% ==================== Ensure Backward Compatibility ====================
if ~isfield(config_results, 'success')
    config_results.success = true;  % if main returned, treat as success
end
if ~isfield(config_results, 'lambda1')
    error('module6_hyperparameter_config:missing_output', ...
          'Expected lambda1 field not found in results');
end

end