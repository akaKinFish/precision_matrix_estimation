function T_xi_v = module2_residual_transfer_function(T_jv, L, varargin)
% T_xi_v = I_p - L*T_jv    (default)
% Optional: T_xi_v = Σ_xi (L Σ_jj L' + Σ_xi)^{-1}  (if 'alternative_computation' true)

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

[n, p_tjv] = size(T_jv);
[p, nL]    = size(L);
if nL ~= n || p_tjv ~= p
    error('module2_residual_transfer_function:dimension_mismatch', 'L or T_jv size mismatch');
end

try
    if opt.alternative_computation
        if isempty(opt.Sigma_xi_xi) || isempty(opt.Sigma_jj)
            % Fallback to exact complementary identity
            T_xi_v = eye(p) - L * T_jv;
        else
            A = L * opt.Sigma_jj * L' + opt.Sigma_xi_xi;  % p×p
            % Use stable linear solve instead of explicit inverse:
            T_xi_v = (A \ opt.Sigma_xi_xi')';             % == Sigma_xi_xi / A
        end
    else
        T_xi_v = eye(p) - L * T_jv;
    end

    if opt.enforce_hermitian
        T_xi_v = (T_xi_v + T_xi_v')/2;
    end
catch ME
    error('module2_residual_transfer_function:computation_failed', ...
          'RTF computation failed: %s', ME.message);
end

if opt.validate_properties && ~opt.alternative_computation
    tol = opt.numerical_tolerance;
    comp_err = norm(T_xi_v + L*T_jv - eye(p), 'fro') / sqrt(p);
    if comp_err > 10*tol
        warning('module2_residual_transfer_function:complementary_identity_failed', ...
                'Complementary identity deviation = %.3e', comp_err);
    end
end

% Drop tiny imaginary leakage ONLY when L,T_jv are real
if isreal(L) && isreal(T_jv)
    if max(abs(imag(T_xi_v(:)))) < 1e-13
        T_xi_v = real(T_xi_v);
    end
end
end
