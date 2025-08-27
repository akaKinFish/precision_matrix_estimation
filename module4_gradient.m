function gradient_results = module4_gradient(input_data, gradient_params)
% MODULE4_GRADIENT - Main interface for Module 4 gradient computation
%
% Syntax:
%   gradient_results = module4_gradient(input_data, gradient_params)
%
% Description:
%   Wrapper function providing the main interface for Module 4 objective gradient
%   computation. Routes to the comprehensive implementation while maintaining
%   backward compatibility and providing simplified access.
%   
%   This function serves as the primary entry point for Module 4, automatically
%   handling input validation, parameter setup, and result formatting.
%
% Input Arguments:
%   input_data - (struct) Input data structure containing:
%     .precision_matrices         - (cell array, Fx1) Current precision estimates Γ_ω
%     .whitened_covariances      - (cell array, Fx1) Whitened covariances Σ_ω
%     .smoothing_kernel          - (double array, FxF) Kernel matrix k_{ω,ω'}
%     .weight_matrix             - (double array, pxp) Weight matrix W^Γ
%
%   gradient_params - (struct) Optional parameters (defaults applied if empty):
%     .lambda1                   - (double) Smoothing parameter (default: 0.01)
%     .lambda2                   - (double) L1 penalty parameter (default: 0.01)
%     .penalize_diagonal         - (logical) Include diagonal in L1 penalty (default: false)
%     .use_graph_laplacian       - (logical) Use Laplacian method (default: true)
%     .verbose                   - (logical) Display progress (default: false)
%     [Additional parameters supported by main implementation]
%
% Output Arguments:
%   gradient_results - (struct) Complete gradient computation results:
%     .smooth_gradients          - (cell array, Fx1) Smooth gradient matrices G_ω
%     .gradient_components       - (struct) Breakdown of gradient components
%     .computation_stats         - (struct) Computation statistics and timing
%     .hermitian_violations      - (double array, Fx1) Symmetry constraint violations
%     .success                   - (logical) Overall computation success
%
% Examples:
%   % Basic usage - compute gradients with default parameters
%   results = module4_gradient(input_data);
%   smooth_grads = results.smooth_gradients;
%   
%   % Advanced usage with custom parameters
%   params = struct();
%   params.lambda1 = 0.05;
%   params.use_graph_laplacian = false;
%   params.verbose = true;
%   results = module4_gradient(input_data, params);
%   
%   % Access detailed statistics
%   stats = results.computation_stats;
%   fprintf('Computation time: %.3fs\n', stats.total_computation_time);
%   fprintf('Hermitian enforcements: %d\n', stats.hermitian_enforcement_count);
%
% Notes:
%   This wrapper function:
%   - Provides backward compatibility with existing code
%   - Applies sensible defaults for all parameters
%   - Routes to the full implementation (module4_objective_gradient_main)
%   - Maintains consistent error handling and validation
%   - Supports all advanced features of the complete implementation
%
% See also: MODULE4_OBJECTIVE_GRADIENT_MAIN, MODULE4_OBJECTIVE_EVALUATION,
%           OBJECTIVEGRADIENTCOMPUTER
%
% Author: [Your Name]
% Date: [Current Date]
% Version: 1.0

% ==================== Input Validation ====================
if nargin < 1
    error('module4_gradient:insufficient_input', ...
          'At least input_data is required');
end

if nargin < 2
    gradient_params = struct();
end

% ==================== Route to Main Implementation ====================
% The main implementation handles all validation, computation, and result formatting
try
    gradient_results = module4_objective_gradient_main(input_data, gradient_params);
    
catch ME
    % Re-throw with consistent error identifier
    error('module4_gradient:computation_failed', ...
          'Module 4 gradient computation failed: %s', ME.message);
end

% ==================== Ensure Backward Compatibility ====================
% Add any legacy fields or formatting if needed for compatibility
if ~isfield(gradient_results, 'success')
    gradient_results.success = true;  % If we got here, computation succeeded
end

% Ensure all expected fields are present
if ~isfield(gradient_results, 'smooth_gradients')
    error('module4_gradient:missing_output', ...
          'Expected smooth_gradients field not found in results');
end

end