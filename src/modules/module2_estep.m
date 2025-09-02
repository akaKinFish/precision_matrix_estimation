function estep_results = module2_estep(input_data, estep_params)
% MODULE2_INITIALIZATION - E-step computation wrapper function
%
% This thin wrapper calls module2_estep_main and preserves a stable API.

% ==================== Input Validation ====================
if nargin < 1
    error('module2_initialization:insufficient_input', ...
          'At least input_data is required');
end
if nargin < 2
    estep_params = struct();
end

% ==================== Route to Main Implementation ====================
try
    estep_results = module2_estep_main(input_data, estep_params);
catch ME_original
    % Re-throw with a consistent identifier and try to attach the original cause
    try
        ME = MException('module2_initialization:computation_failed', ...
                        'Module 2 E-step computation failed.');
        % addCause is available in newer MATLAB versions
        ME = addCause(ME, ME_original);
        throw(ME);
    catch
        error('module2_initialization:computation_failed', ...
              'Module 2 E-step computation failed: %s', ME_original.message);
    end
end

% ==================== Ensure Backward Compatibility ====================
if ~isfield(estep_results, 'success')
    estep_results.success = true;  % if main returned, treat as success
end
if ~isfield(estep_results, 'source_transfer_functions')
    error('module2_initialization:missing_output', ...
          'Expected source_transfer_functions field not found in results');
end

end