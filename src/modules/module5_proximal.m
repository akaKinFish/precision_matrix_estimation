function [Gamma_final, proximal_results] = module5_proximal(input_data, proximal_params)
% MODULE5_PROXIMAL - Main entry point for Module 5 proximal gradient updates
%
% Syntax:
%   [Gamma_final, proximal_results] = module5_proximal(input_data, proximal_params)
%
% Description:
%   Entry point for Module 5 proximal gradient solver. This function provides
%   a unified interface for sparse precision matrix estimation using proximal
%   gradient methods with complex amplitude soft thresholding.
%   
%   This function is a wrapper around module5_proximal_main for consistency
%   with the module naming convention.
%
% Input Arguments:
%   input_data - (struct) Required fields:
%     .whitened_covariances     - (cell array, Fx1) Whitened covariances Σ̃_ω
%     .initial_precision        - (cell array, Fx1) Initial precision estimates
%     .smoothing_kernel         - (double array, FxF) Kernel matrix k_{ω,ω'}
%     .weight_matrix           - (double array, pxp) Weight matrix W^Γ
%     .active_set_masks        - (cell array, Fx1) Active set masks A_ω
%
%   proximal_params - (struct) Optional parameters (see module5_proximal_main)
%
% Output Arguments:
%   Gamma_final - (cell array, Fx1) Final precision matrix estimates
%   proximal_results - (struct) Detailed optimization results
%
% Examples:
%   % Basic usage
%   [Gamma_opt, results] = module5_proximal(input_data, struct());
%   
%   % Advanced usage
%   params.lambda2 = 0.01;
%   params.max_iter = 500;
%   params.verbose = true;
%   [Gamma_opt, results] = module5_proximal(input_data, params);
%
% See also: MODULE5_PROXIMAL_MAIN, MODULE5_SINGLE_PROXIMAL_STEP,
%           MODULE4_OBJECTIVE_GRADIENT_MAIN
%
% Author: [Your Name]
% Date: [Current Date]
% Version: 1.0

% Direct call to the main implementation
[Gamma_final, proximal_results] = module5_proximal_main(input_data, proximal_params);

end