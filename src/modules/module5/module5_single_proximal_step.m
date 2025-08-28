function [Gamma_new, step_info] = module5_single_proximal_step(Gamma_current, gradient, alpha, active_mask, aux_data, params)
% MODULE5_SINGLE_PROXIMAL_STEP - Single proximal gradient update with backtracking
%
% Syntax:
%   [Gamma_new, step_info] = module5_single_proximal_step(Gamma_current, gradient, alpha, active_mask, aux_data, params)
%
% Description:
%   Performs a single proximal gradient step for one frequency with:
%   1. Gradient step: Γ^tmp = Γ^(k) - α G
%   2. Complex amplitude soft thresholding for off-diagonal elements
%   3. Diagonal processing (real values)
%   4. Active set projection
%   5. Hermitian symmetrization
%   6. PD check with backtracking if needed
%
% Input Arguments:
%   Gamma_current - (double array, pxp) Current precision matrix Γ_ω^(k)
%   gradient      - (double array, pxp) Smooth gradient G_ω
%   alpha         - (double) Step size α
%   active_mask   - (logical array, pxp) Active set mask A_ω
%   aux_data      - (struct) Auxiliary data:
%     .lambda2    - (double) L1 penalty parameter
%   params        - (struct) Parameters:
%     .mode              - (string) 'simplified'|'joint'
%     .beta_backtrack    - (double) Backtracking shrinkage factor
%     .max_backtrack     - (integer) Maximum backtracking iterations
%     .psd_tolerance     - (double) PD check tolerance (default: 1e-12)
%     .verbose           - (logical) control warnings (default: false)
%
% Output Arguments:
%   Gamma_new - (double array, pxp) Updated precision matrix Γ_ω^(k+1)
%   step_info - (struct) Contains:
%     .backtrack_count       - (integer) Number of backtracking steps
%     .final_step_size       - (double) Final step size used
%     .objective_decrease    - (double) Objective function decrease (approx)
%     .psd_violations        - (integer) Number of PD violations encountered
%     .success               - (logical) Step success indicator

% ==================== Input Validation ====================
if nargin < 6
    error('module5_single_proximal_step:insufficient_input', ...
          'All 6 arguments are required');
end

% Validate matrix dimensions
p = size(Gamma_current, 1);
if ~isequal(size(Gamma_current), [p, p]) || ~isequal(size(gradient), [p, p])
    error('module5_single_proximal_step:dimension_mismatch', ...
          'Gamma_current and gradient must be same size');
end

if ~isequal(size(active_mask), [p, p])
    error('module5_single_proximal_step:active_mask_size', ...
          'active_mask must be %dx%d', p, p);
end

% Validate inputs are finite
if ~all(isfinite(Gamma_current(:))) || ~all(isfinite(gradient(:)))
    error('module5_single_proximal_step:invalid_values', ...
          'Gamma_current and gradient must contain finite values');
end

if alpha <= 0 || ~isfinite(alpha)
    error('module5_single_proximal_step:invalid_step_size', ...
          'Step size alpha must be positive and finite');
end

% Set defaults
if ~isfield(params, 'psd_tolerance')
    params.psd_tolerance = 1e-12;
end
if ~isfield(params, 'verbose')
    params.verbose = false;
end

% ==================== Initialize Step Info ====================
step_info = struct();
step_info.backtrack_count = 0;
step_info.psd_violations = 0;
step_info.success = false;

current_step_size = alpha;
lambda2 = aux_data.lambda2;

% ==================== Backtracking Loop ====================
for bt_iter = 0:params.max_backtrack
    
    % Step 1: Gradient step
    Gamma_temp = Gamma_current - current_step_size * gradient;
    
    % Step 2: Complex amplitude soft thresholding (off-diagonal only)
    tau_threshold = lambda2 * current_step_size;
    Gamma_prox = module5_soft_threshold_complex(Gamma_temp, tau_threshold, params.mode);
    
    % Step 3: Active set projection
    Gamma_projected = module5_active_set_projection(Gamma_prox, active_mask);
    
    % Step 4: Hermitian symmetrization
    Gamma_symmetric = module5_hermitian_symmetrize(Gamma_projected, params.mode);
    
    % Step 5: PD check
    [isPSD, ~] = module5_psd_check(Gamma_symmetric);
    
    if isPSD
        % Success: accept the step
        Gamma_new = Gamma_symmetric;
        step_info.final_step_size = current_step_size;
        step_info.success = true;
        break;
    else
        % PD violation: record and backtrack
        step_info.psd_violations = step_info.psd_violations + 1;
        step_info.backtrack_count = step_info.backtrack_count + 1;
        
        if bt_iter == params.max_backtrack
            % Maximum backtracking reached: use regularized fallback
            if params.verbose
                warning('module5_single_proximal_step:max_backtrack_reached', ...
                        'Maximum backtracking iterations reached, using regularized solution');
            end
            
            % Minimal regularization to ensure PD
            min_eigenval = 1e-8;
            [V, D] = eig(Gamma_symmetric);
            eigenvals = real(diag(D));
            eigenvals = max(eigenvals, min_eigenval);
            Gamma_new = V * diag(eigenvals) * V';
            Gamma_new = (Gamma_new + Gamma_new') / 2;  % Ensure numerical Hermitian
            
            step_info.final_step_size = current_step_size;
            step_info.success = true;  % Mark as success with regularization
            break;
        else
            % Reduce step size and retry
            current_step_size = current_step_size * params.beta_backtrack;
        end
    end
end

% ==================== Compute Step Quality Metrics ====================
if step_info.success
    % Compute objective decrease (approximate)
    norm_change = norm(Gamma_new - Gamma_current, 'fro');
    step_info.objective_decrease = current_step_size * real(trace(gradient' * (Gamma_new - Gamma_current)));
    
    % Additional quality metrics
    step_info.matrix_change_norm = norm_change;
    step_info.relative_change = norm_change / max(norm(Gamma_current, 'fro'), 1e-12);
    
    % Check final matrix properties
    step_info.final_is_hermitian = norm(Gamma_new - Gamma_new', 'fro') < 1e-12;
    step_info.final_diagonal_real = all(abs(imag(diag(Gamma_new))) < 1e-12);
    step_info.final_condition_number = cond(Gamma_new);
else
    % Failed step
    step_info.objective_decrease = 0;
    step_info.matrix_change_norm = 0;
    step_info.relative_change = 0;
    step_info.final_is_hermitian = false;
    step_info.final_diagonal_real = false;
    step_info.final_condition_number = Inf;
    
    % Fallback to previous matrix
    Gamma_new = Gamma_current;
end

end
