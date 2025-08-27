function [trace_gradients, computation_stats] = module4_trace_gradient(covariance_matrices, params)
% MODULE4_TRACE_GRADIENT - Compute gradients of trace terms
%
% Syntax:
%   [trace_gradients, computation_stats] = module4_trace_gradient(covariance_matrices, params)
%
% Description:
%   Computes the gradient of tr(Σ_ω Γ_ω) with respect to each precision matrix Γ_ω.
%   Ensures proper Hermitian symmetry for complex matrices.
%   
%   Mathematical formula:
%   ∇_{Γ_ω} [tr(Σ_ω Γ_ω)] = Σ_ω
%   
%   For Hermitian constraint, uses: ∇ = (Σ_ω + Σ_ω^H) / 2
%
% Input Arguments:
%   covariance_matrices - (cell array, Fx1) Covariance matrices {Σ_ω}
%   params - (struct) Parameters:
%     .force_hermitian          - (logical) Force Hermitian result (default: true)
%     .symmetrization_tolerance - (double) Tolerance for symmetry check (default: 1e-10)
%     .verbose                  - (logical) Display progress information (default: false)
%
% Output Arguments:
%   trace_gradients - (cell array, Fx1) Gradient matrices {Σ_ω} (Hermitian)
%   computation_stats - (struct) Contains:
%     .computation_time         - (double) Total computation time
%     .hermitian_enforcement_count - (integer) Number of symmetrizations
%     .symmetry_errors          - (double array, Fx1) Symmetry errors per frequency
%     .is_complex               - (logical) Whether matrices are complex-valued
%
% Examples:
%   % Basic usage
%   params = struct('force_hermitian', true, 'verbose', false);
%   [grads, stats] = module4_trace_gradient(sigma_cells, params);
%   
%   % Check symmetry enforcement
%   if stats.hermitian_enforcement_count > 0
%       fprintf('Symmetrized %d matrices\n', stats.hermitian_enforcement_count);
%   end
%
% Mathematical Background:
%   The trace term tr(Σ Γ) has gradient:
%   d/dΓ [tr(Σ Γ)] = Σ
%   
%   For Hermitian optimization, we ensure the gradient is Hermitian:
%   ∇ = (Σ + Σ^H) / 2
%   
%   This is essential for maintaining the Hermitian constraint on Γ.
%
% See also: MODULE4_OBJECTIVE_GRADIENT_MAIN, MODULE4_LOG_DET_GRADIENT
%
% Author: [Your Name]  
% Date: [Current Date]
% Version: 1.0

% ==================== Input Validation ====================
if nargin < 1
    error('module4_trace_gradient:insufficient_input', ...
          'covariance_matrices is required');
end

if nargin < 2
    params = struct();
end

if ~iscell(covariance_matrices) || isempty(covariance_matrices)
    error('module4_trace_gradient:invalid_input', ...
          'covariance_matrices must be a non-empty cell array');
end

% Set default parameters
defaults = struct();
defaults.force_hermitian = true;
defaults.symmetrization_tolerance = 1e-10;
defaults.verbose = false;

field_names = fieldnames(defaults);
for i = 1:numel(field_names)
    fname = field_names{i};
    if ~isfield(params, fname)
        params.(fname) = defaults.(fname);
    end
end

% ==================== Initialize ====================
F = length(covariance_matrices);
p = size(covariance_matrices{1}, 1);

trace_gradients = cell(F, 1);
computation_stats = struct();
computation_stats.computation_time = 0;
computation_stats.hermitian_enforcement_count = 0;
computation_stats.symmetry_errors = zeros(F, 1);
computation_stats.is_complex = false;

main_tic = tic;

% ==================== Validate Matrix Dimensions ====================
for f = 1:F
    if ~isnumeric(covariance_matrices{f}) || ~ismatrix(covariance_matrices{f})
        error('module4_trace_gradient:invalid_matrix', ...
              'covariance_matrices{%d} must be a numeric matrix', f);
    end
    
    if ~isequal(size(covariance_matrices{f}), [p, p])
        error('module4_trace_gradient:size_mismatch', ...
              'covariance_matrices{%d} must be %dx%d, got %dx%d', ...
              f, p, p, size(covariance_matrices{f}, 1), size(covariance_matrices{f}, 2));
    end
end

% Check if matrices are complex-valued
computation_stats.is_complex = any(cellfun(@(x) ~isreal(x), covariance_matrices));

% ==================== Process Each Frequency ====================
for f = 1:F
    Sigma_f = covariance_matrices{f};
    
    % The gradient of tr(Σ_ω Γ_ω) with respect to Γ_ω is simply Σ_ω
    trace_grad = Sigma_f;
    
    % ==================== Force Hermitian Symmetry ====================
    if params.force_hermitian
        % Compute Hermitian part: (Σ + Σ^H) / 2
        trace_grad_hermitian = (trace_grad + trace_grad') / 2;
        
        % Measure symmetry error
        symmetry_error = norm(trace_grad - trace_grad_hermitian, 'fro');
        computation_stats.symmetry_errors(f) = symmetry_error;
        
        if symmetry_error > params.symmetrization_tolerance
            computation_stats.hermitian_enforcement_count = computation_stats.hermitian_enforcement_count + 1;
            
            if params.verbose
                fprintf('  Trace gradient freq %d: Hermitian enforcement (error: %.2e)\n', ...
                        f, symmetry_error);
            end
        end
        
        trace_grad = trace_grad_hermitian;
    end
    
    trace_gradients{f} = trace_grad;
end

% ==================== Final Statistics ====================
computation_stats.computation_time = toc(main_tic);
computation_stats.max_symmetry_error = max(computation_stats.symmetry_errors);
computation_stats.mean_symmetry_error = mean(computation_stats.symmetry_errors);

if params.verbose
    fprintf('  Trace gradients: %d/%d matrices symmetrized\n', ...
            computation_stats.hermitian_enforcement_count, F);
    if computation_stats.max_symmetry_error > 1e-12
        fprintf('  Max symmetry error: %.2e\n', computation_stats.max_symmetry_error);
    end
end

end