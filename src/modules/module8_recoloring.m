function recoloring_results = module8_recoloring(input_data, recoloring_params)
% MODULE8_RECOLORING - Recoloring wrapper function for inverse whitening transformation
%
% This thin wrapper calls module8_recoloring_main and preserves a stable API.

% ==================== Input Validation ====================
if nargin < 1
    error('module8_recoloring:insufficient_input', ...
          'At least input_data is required');
end
if nargin < 2
    recoloring_params = struct();
end

% ==================== Route to Main Implementation ====================
try
    recoloring_results = module8_recoloring_main(input_data, recoloring_params);
catch ME_original
    % Re-throw with a consistent identifier and try to attach the original cause
    try
        ME = MException('module8_recoloring:computation_failed', ...
                        'Module 8 recoloring computation failed.');
        % addCause is available in newer MATLAB versions
        ME = addCause(ME, ME_original);
        throw(ME);
    catch
        error('module8_recoloring:computation_failed', ...
              'Module 8 recoloring computation failed: %s', ME_original.message);
    end
end

% ==================== Ensure Backward Compatibility ====================
if ~isfield(recoloring_results, 'success')
    recoloring_results.success = true;  % if main returned, treat as success
end
if ~isfield(recoloring_results, 'recolored_precision_matrices')
    error('module8_recoloring:missing_output', ...
          'Expected recolored_precision_matrices field not found in results');
end
end
