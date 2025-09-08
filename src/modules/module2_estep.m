function estep_results = module2_estep(input_data, estep_params)
% MODULE2_ESTEP - E-step computation wrapper (backward-compatible)
% - Calls module2_estep_main(input_data, estep_params)
% - Adds legacy aliases without changing main outputs.

if nargin < 1
    error('module2_estep:insufficient_input', 'At least input_data is required');
end
if nargin < 2
    estep_params = struct();
end

try
    core = module2_estep_main(input_data, estep_params);
catch ME
    ME2 = MException('module2_estep:computation_failed', ...
                     'Module 2 E-step computation failed: %s', ME.message);
    try, ME2 = addCause(ME2, ME); catch, end
    throw(ME2);
end

estep_results = core;
if ~isfield(estep_results, 'success') || isempty(estep_results.success)
    estep_results.success = true;
end

% Legacy aliases (keep both names available)
if isfield(core, 'transfer_functions')
    estep_results.source_transfer_functions = core.transfer_functions;
elseif isfield(core, 'source_transfer_functions')
    estep_results.transfer_functions = core.source_transfer_functions;
else
    error('module2_estep:missing_output_field', ...
        'Neither "transfer_functions" nor "source_transfer_functions" found.');
end

if isfield(core, 'residual_covariances') ...
   && ~isfield(estep_results, 'residual_empirical_covariances')
    estep_results.residual_empirical_covariances = core.residual_covariances;
end

if isfield(core, 'posterior_source_covariances') ...
   && ~isfield(estep_results, 'posterior_covariances')
    estep_results.posterior_covariances = core.posterior_source_covariances;
end

if isfield(core, 'initial_precision_matrices') ...
   && ~isfield(estep_results, 'initial_precision')
    estep_results.initial_precision = core.initial_precision_matrices;   % source-domain n×n
end

if isfield(core, 'source_second_moments') ...
   && ~isfield(estep_results, 'Sjj_hat')
    estep_results.Sjj_hat = core.source_second_moments;                  % n×n
end

end
