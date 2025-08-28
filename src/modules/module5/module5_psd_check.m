function [isPSD, chol_info] = module5_psd_check(Gamma, tolerance)
% MODULE5_PSD_CHECK - Check positive semidefiniteness via Cholesky decomposition
%
% MODULE5_PSD_CHECK - Check positive definiteness (PD) via Cholesky decomposition
%
% Syntax:
%   isPSD = module5_psd_check(Gamma)
%   [isPSD, chol_info] = module5_psd_check(Gamma, tolerance)
%
% Description:
%   Checks if a Hermitian matrix is positive definite using Cholesky
%   decomposition. This is stricter than “semidefinite” and matches the
%   requirement for log-determinant.
% ...
%   Uses MATLAB's chol() function which is optimized and numerically stable.
%   Also provides additional diagnostics for debugging numerical issues.
%
% Input Arguments:
%   Gamma     - (double array, pxp) Matrix to check (should be Hermitian)
%   tolerance - (double) Tolerance for numerical checks (default: 1e-12)
%
% Output Arguments:
%   isPSD     - (logical) True if matrix is positive semidefinite
%   chol_info - (struct) Detailed information:
%     .chol_success       - (logical) Whether Cholesky succeeded
%     .chol_flag          - (integer) Cholesky return flag (0 = success)
%     .pivot_column       - (integer) Failed pivot column (if any)
%     .min_diagonal       - (double) Minimum diagonal element
%     .condition_number   - (double) Condition number estimate
%     .hermitian_error    - (double) Hermitian symmetry violation
%     .determinant_sign   - (integer) Sign of determinant (1, 0, -1)
%
% Examples:
%   % Basic PSD check
%   isPSD = module5_psd_check(Gamma_matrix);
%   
%   % Detailed diagnostics
%   [isPSD, info] = module5_psd_check(Gamma_matrix, 1e-10);
%   if ~isPSD
%       fprintf('PSD check failed at column %d\n', info.pivot_column);
%       fprintf('Minimum diagonal: %.2e\n', info.min_diagonal);
%       fprintf('Condition number: %.2e\n', info.condition_number);
%   end
%   
%   % Check if matrix needs regularization
%   if info.condition_number > 1e12
%       fprintf('Matrix is ill-conditioned, consider regularization\n');
%   end
%
% Mathematical Background:
%   A Hermitian matrix Γ is positive semidefinite if and only if:
%   1. All eigenvalues are non-negative
%   2. Cholesky decomposition Γ = L*L^H exists
%   3. All leading principal minors are non-negative
%   
%   Cholesky is preferred over eig() for large matrices as it's O(n³/3)
%   versus O(n³) and more numerically stable.
%
% See also: MODULE5_SINGLE_PROXIMAL_STEP, MODULE5_BACKTRACKING
%
% Author: [Your Name]
% Date: [Current Date]
% Version: 1.0

% ==================== Input Validation ====================
if nargin < 1
    error('module5_psd_check:insufficient_input', ...
          'Matrix Gamma is required');
end

if nargin < 2
    tolerance = 1e-12;
end

% Validate input matrix
if ~isnumeric(Gamma) || ~ismatrix(Gamma)
    error('module5_psd_check:invalid_input_type', ...
          'Gamma must be a numeric matrix');
end

[p, q] = size(Gamma);
if p ~= q
    error('module5_psd_check:not_square', ...
          'Gamma must be square, got %dx%d', p, q);
end

if ~all(isfinite(Gamma(:)))
    error('module5_psd_check:invalid_values', ...
          'Gamma must contain finite values');
end

% Validate tolerance
if ~isscalar(tolerance) || ~isreal(tolerance) || tolerance <= 0
    error('module5_psd_check:invalid_tolerance', ...
          'tolerance must be a positive real scalar');
end

% ==================== Initialize Output Structure ====================
chol_info = struct();
chol_info.chol_success = false;
chol_info.chol_flag = -1;
chol_info.pivot_column = -1;
chol_info.min_diagonal = min(real(diag(Gamma)));
chol_info.condition_number = NaN;
chol_info.hermitian_error = norm(Gamma - Gamma', 'fro');
chol_info.determinant_sign = 0;

% ==================== Hermitian Check ====================
if chol_info.hermitian_error > tolerance
    warning('module5_psd_check:not_hermitian', ...
            'Matrix is not Hermitian (error: %.2e), PSD check may be unreliable', ...
            chol_info.hermitian_error);
end

% ==================== Cholesky Decomposition Test ====================
try
    [L, chol_flag] = chol(Gamma, 'lower');
    
    chol_info.chol_flag = chol_flag;
    
    if chol_flag == 0
        % Cholesky succeeded - matrix is PSD
        chol_info.chol_success = true;
        isPSD = true;
        
        % Additional diagnostics
        try
            chol_info.condition_number = cond(Gamma);
        catch
            % Condition number computation failed, use alternative
            chol_info.condition_number = Inf;
        end
        
        % Determinant sign (positive for PD)
        log_det = 2 * sum(log(real(diag(L))));
        if isfinite(log_det) && log_det > -50  % Avoid underflow
            chol_info.determinant_sign = 1;
        else
            chol_info.determinant_sign = 0;  % Numerically singular
        end
        
    else
        % Cholesky failed - matrix is not PSD
        chol_info.chol_success = false;
        chol_info.pivot_column = chol_flag;
        isPSD = false;
        
        % Try to estimate condition number via SVD (more expensive but robust)
        try
            s = svd(Gamma);
            if min(s) > tolerance
                chol_info.condition_number = max(s) / min(s);
                chol_info.determinant_sign = 1;
            else
                chol_info.condition_number = Inf;
                chol_info.determinant_sign = 0;
            end
        catch
            chol_info.condition_number = Inf;
            chol_info.determinant_sign = -1;  % Unknown/problematic
        end
    end
    
catch ME
    % Cholesky computation failed entirely
    isPSD = false;
    chol_info.chol_success = false;
    chol_info.condition_number = Inf;
    chol_info.determinant_sign = -1;
    
    warning('module5_psd_check:chol_computation_failed', ...
            'Cholesky computation failed: %s', ME.message);
end

% ==================== Additional Diagnostics ====================
% Quick diagonal check
if any(real(diag(Gamma)) <= 0)
    % Negative or zero diagonal elements - definitely not PSD
    isPSD = false;
    if chol_info.chol_success  % This shouldn't happen, but just in case
        warning('module5_psd_check:inconsistent_result', ...
                'Cholesky succeeded but diagonal has non-positive elements');
    end
end

% ==================== Summary for Debugging ====================
% Store summary information useful for debugging
chol_info.matrix_size = p;
chol_info.matrix_norm = norm(Gamma, 'fro');
chol_info.diagonal_sum = sum(real(diag(Gamma)));
chol_info.off_diagonal_norm = norm(Gamma - diag(diag(Gamma)), 'fro');

end