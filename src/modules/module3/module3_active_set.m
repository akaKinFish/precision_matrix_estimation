function active_set_results = module3_active_set(input_data, active_set_params)
% MODULE3_ACTIVE_SET - Active edge set selection for sparse precision estimation
%
% Syntax:
%   active_set_results = module3_active_set(input_data, active_set_params)
%
% Description:
%   Performs active edge set selection to reduce optimization dimension by 
%   selecting plausible edges based on proxy statistics. This wrapper function
%   provides a simplified interface to the complete Module 3 pipeline.
%
%   The module implements both correlation-based and precision-based edge proxy 
%   computation with node-level filtering and combined active set construction.
%
% Input Arguments:
%   input_data - (struct) Required fields:
%     .whitened_covariances        - (cell array, Fx1) Whitened empirical covariances
%     .initial_precision_matrices  - (cell array, Fx1) Initial precision estimates (optional)
%     .frequencies                 - (double array, Fx1) Frequency values
%
%   active_set_params - (struct) Optional parameters:
%     .proxy_method               - ('correlation'|'precision') Method for proxy computation
%     .quantile_level             - (double) Quantile level for threshold (default: 0.1)
%     .force_diagonal_active      - (logical) Always keep diagonal entries active
%     .min_active_edges          - (integer) Minimum number of edges to keep active
%     .verbose                   - (logical) Display progress information
%
% Output Arguments:
%   active_set_results - (struct) Contains:
%     .edge_proxies              - (cell array, Fx1) Edge proxy values c_ij(f)
%     .threshold_value           - (double) Computed threshold tau
%     .edge_active_mask          - (logical array, pxpxF) Edge activity masks
%     .node_active_mask          - (logical array, pxF) Node activity masks  
%     .combined_active_mask      - (logical array, pxpxF) Final combined active masks
%     .active_edge_statistics    - (struct) Statistics about active edges
%     .computation_stats         - (struct) Computation statistics
%     .success                   - (logical) Overall success indicator
%
% Examples:
%   % Basic usage with whitened covariances
%   input_data.whitened_covariances = whitened_cov;
%   input_data.frequencies = frequencies;
%   results = module3_active_set(input_data, struct());
%   
%   % Advanced usage with precision-based proxy
%   input_data.initial_precision_matrices = init_precision;
%   params.proxy_method = 'precision';
%   params.quantile_level = 0.05;
%   results = module3_active_set(input_data, params);
%
% See also: MODULE3_ACTIVE_SET_MAIN, MODULE3_EDGE_PROXY_COMPUTATION,
%           MODULE3_THRESHOLD_DETERMINATION, MODULE3_COMBINED_ACTIVE_SET
%
% Author: [Your Name]
% Date: [Current Date]  
% Version: 1.0

% Input validation
if nargin < 1
    error('module3_active_set:insufficient_input', ...
          'At least input_data is required');
end

if nargin < 2
    active_set_params = struct();
end

% Call the main implementation
active_set_results = module3_active_set_main(input_data, active_set_params);

end