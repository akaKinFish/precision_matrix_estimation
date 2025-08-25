function T_xi_v = module2_residual_transfer_function(T_jv, L, varargin)
% MODULE2_RESIDUAL_TRANSFER_FUNCTION - Compute Residual Transfer Function (RTF)
%
% Syntax:
%   T_xi_v = module2_residual_transfer_function(T_jv, L)
%   T_xi_v = module2_residual_transfer_function(T_jv, L, 'Name', Value)
%
% Description:
%   Residual TF for the E-step:
%       T_ξv = I_p - L * T_jv
%   NOTE: RTF is a contraction operator in the ideal model (||T_ξv||_2 ≤ 1),
%   not an idempotent projector in general.
%
% Name-Value Arguments:
%   'validate_properties'  (logical)  Validate complementary identity. Default: true
%   'numerical_tolerance'  (double)   Tolerance for validations. Default: 1e-12
%   'verbose'              (logical)  Print diagnostics. Default: false
%   'enforce_hermitian'    (logical)  Post symmetrization (breaks identity). Default: false
%   'alternative_computation' (logical) If true and both 'Sigma_xi_xi' and 'Sigma_jj'
%                           are provided, compute T_ξv = Σ_ξξ (L Σ_jj L^H + Σ_ξξ)^(-1).
%                           Default: false
%   'Sigma_xi_xi'          (matrix)   Sensor noise covariance (p×p), used only when
%                           'alternative_computation' is true.
%   'Sigma_jj'             (matrix)   Source prior covariance (n×n), used only when
%                           'alternative_computation' is true.
%
% Output:
%   T_xi_v : (complex, p×p)

% ---- Input checks -------------------------------------------------------
if ~isnumeric(T_jv) || ndims(T_jv) ~= 2
    error('module2_residual_transfer_function:invalid_transfer_function', ...
          'T_jv must be a 2D numeric matrix');
end
[n, p_tjv] = size(T_jv);

if ~isnumeric(L) || ndims(L) ~= 2
    error('module2_residual_transfer_function:invalid_leadfield', ...
          'Leadfield matrix L must be a 2D numeric array');
end
[p, n_sources] = size(L);

if n_sources ~= n
    error('module2_residual_transfer_function:dimension_mismatch_sources', ...
          'Leadfield sources (%d) must match T_jv sources (%d)', n_sources, n);
end
if p_tjv ~= p
    error('module2_residual_transfer_function:dimension_mismatch_sensors', ...
          'T_jv sensors (%d) must match leadfield sensors (%d)', p_tjv, p);
end
if any(~isfinite(T_jv(:)))
    error('module2_residual_transfer_function:invalid_transfer_values', ...
          'T_jv contains NaN or Inf');
end
if any(~isfinite(L(:)))
    error('module2_residual_transfer_function:invalid_leadfield_values', ...
          'Leadfield contains NaN or Inf');
end

% ---- Parse options ------------------------------------------------------
ip = inputParser;
addParameter(ip, 'validate_properties', true,  @(x)islogical(x)&&isscalar(x));
addParameter(ip, 'numerical_tolerance', 1e-12, @(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(ip, 'verbose', false,            @(x)islogical(x)&&isscalar(x));
addParameter(ip, 'enforce_hermitian', false,  @(x)islogical(x)&&isscalar(x));
addParameter(ip, 'alternative_computation', false, @(x)islogical(x)&&isscalar(x));
addParameter(ip, 'Sigma_xi_xi', [], @(x)isnumeric(x) || isempty(x));
addParameter(ip, 'Sigma_jj',    [], @(x)isnumeric(x) || isempty(x));
parse(ip, varargin{:});
opt = ip.Results;

% ---- Main computation ---------------------------------------------------
try
    if opt.alternative_computation
        % Optional equivalent form (requires Σ_ξξ and Σ_jj)
        if isempty(opt.Sigma_xi_xi) || isempty(opt.Sigma_jj)
            % Fall back to the standard identity-preserving formula
            if opt.verbose
                fprintf(['[RTF] alternative_computation requested but Sigma_xi_xi/Sigma_jj ', ...
                         'are missing; falling back to I - L*T_jv.\n']);
            end
            T_xi_v = eye(p) - L * T_jv;
        else
            % Compute A = L Σ_jj L^H + Σ_ξξ and T_ξv = Σ_ξξ A^{-1}
            A = L * opt.Sigma_jj * L' + opt.Sigma_xi_xi;
            % Use a linear solve instead of explicit inverse for stability
            T_xi_v = opt.Sigma_xi_xi / A; % == opt.Sigma_xi_xi * inv(A)
        end
    else
        % Exact algebra that guarantees T_ξv + L*T_jv ≈ I_p (up to round-off)
        T_xi_v = eye(p) - L * T_jv;
    end

    % Optional post-symmetrization (for display only; breaks complementary identity)
    if opt.enforce_hermitian
        T_xi_v = (T_xi_v + T_xi_v') / 2;
    end

catch ME
    error('module2_residual_transfer_function:computation_failed', ...
          'RTF computation failed: %s', ME.message);
end

% ---- Validations (optional) --------------------------------------------
if opt.validate_properties
    tol = opt.numerical_tolerance;
    comp_err = norm(T_xi_v + L*T_jv - eye(p), 'fro') / sqrt(p);
    if comp_err > 10*tol && ~opt.alternative_computation
        % For the alternative path we do not know T_jv’s construction, so only warn when using the standard path
        warning('module2_residual_transfer_function:complementary_identity_failed', ...
                'Complementary identity deviation = %.3e', comp_err);
    end
end

% ---- Output checks ------------------------------------------------------
if any(~isfinite(T_xi_v(:)))
    error('module2_residual_transfer_function:invalid_output', ...
          'T_xi_v contains NaN or Inf');
end
if ~isequal(size(T_xi_v), [p p])
    error('module2_residual_transfer_function:output_dimension_error', ...
          'Output must be %dx%d', p, p);
end
end
