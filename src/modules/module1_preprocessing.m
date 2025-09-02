function preprocessing_results = module1_preprocessing(input_data, preprocessing_params)
% MODULE1_PREPROCESSING - Data preprocessing wrapper function
%
% This thin wrapper calls module1_preprocessing_main and preserves a stable API.

% ==================== Input Validation ====================
if nargin < 1
    error('module1_preprocessing:insufficient_input', ...
          'At least input_data is required');
end
if nargin < 2
    preprocessing_params = struct();
end

% ==================== Route to Main Implementation ====================
try
    % Convert struct params to name-value pairs for main function
    if isempty(preprocessing_params)
        preprocessing_results = module1_preprocessing_main(input_data);
    else
        % Convert struct to name-value pairs
        param_names = fieldnames(preprocessing_params);
        param_values = struct2cell(preprocessing_params);
        name_value_args = [param_names'; param_values'];
        name_value_args = name_value_args(:)';
        
        preprocessing_results = module1_preprocessing_main(input_data, name_value_args{:});
    end
    
catch ME_original
    % Re-throw with a consistent identifier and try to attach the original cause
    try
        ME = MException('module1_preprocessing:computation_failed', ...
                        'Module 1 preprocessing computation failed.');
        % addCause is available in newer MATLAB versions
        ME = addCause(ME, ME_original);
        throw(ME);
    catch
        error('module1_preprocessing:computation_failed', ...
              'Module 1 preprocessing computation failed: %s', ME_original.message);
    end
end

% ==================== Ensure Backward Compatibility ====================
if ~isfield(preprocessing_results, 'success')
    preprocessing_results.success = true;  % if main returned, treat as success
end
if ~isfield(preprocessing_results, 'Sigma_tilde')
    error('module1_preprocessing:missing_output', ...
          'Expected Sigma_tilde field not found in results');
end

end