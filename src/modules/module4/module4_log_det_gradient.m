function [logdet_gradients, computation_stats] = module4_log_det_gradient(precision_matrices, params)
% MODULE4_LOG_DET_GRADIENT - Compute gradients of log-determinant terms (revised)
%
% ∇_{Γ_ω}[-log det(Γ_ω)] = -Γ_ω^{-1}
% Numerical scheme: symmetric projection + Cholesky with progressive
% diagonal regularization (no default pseudoinverse fallback).
%
% Params (with new defaults):
%   .chol_tolerance         (double)  default 1e-12
%   .force_hermitian        (logical) default true
%   .verbose                (logical) default false
%   .allow_pinv_fallback    (logical) default false  % diagnostics only
%
% Outputs:
%   logdet_gradients  - cell{F} of -Γ^{-1}
%   computation_stats - struct with fields:
%       .computation_time
%       .cholesky_failures
%       .condition_numbers
%       .method_used
%
% Author: (revised)

% -------- Input checks --------
if nargin < 1
    error('module4_log_det_gradient:insufficient_input', ...
          'precision_matrices is required');
end
if nargin < 2, params = struct(); end
if ~iscell(precision_matrices) || isempty(precision_matrices)
    error('module4_log_det_gradient:invalid_input', ...
          'precision_matrices must be a non-empty cell array');
end

% -------- Defaults --------
defaults = struct();
defaults.chol_tolerance      = 1e-12;
defaults.force_hermitian     = true;
defaults.verbose             = false;
defaults.allow_pinv_fallback = false;   % NEW: safer default

fn = fieldnames(defaults);
for i = 1:numel(fn)
    f = fn{i};
    if ~isfield(params, f), params.(f) = defaults.(f); end
end

% -------- Init --------
F = numel(precision_matrices);
p = size(precision_matrices{1}, 1);

logdet_gradients = cell(F, 1);
computation_stats = struct();
computation_stats.computation_time  = 0;
computation_stats.cholesky_failures = 0;
computation_stats.condition_numbers = zeros(F, 1);
computation_stats.method_used       = 'cholesky_solve';

t_total = tic;

% -------- Per-frequency --------
for f = 1:F
    Gamma = precision_matrices{f};

    if ~isnumeric(Gamma) || ~ismatrix(Gamma) || ~isequal(size(Gamma), [p,p])
        error('module4_log_det_gradient:invalid_matrix', ...
              'precision_matrices{%d} must be a %dx%d numeric matrix', f, p, p);
    end

    % Enforce exact Hermitian before factorization
    Gamma = (Gamma + Gamma')/2;

    try
        % Progressive regularization if needed
        tol = params.chol_tolerance;
        reg0 = max(tol, tol * mean(diag(Gamma)));
        success = false;

        for attempt = 0:4
            if attempt == 0
                A = Gamma;
            else
                A = Gamma + reg0 * (10^(attempt-1)) * eye(p);
            end

            [R, flag] = chol(A, 'upper');
            if flag == 0
                % Solve A^{-1} via triangular solves (no explicit inv)
                I = eye(p);
                X = R' \ I;         % R^H X = I
                invA = R \ X;       % R * invA = X
                logdet_grad = -invA;

                computation_stats.condition_numbers(f) = cond(A);
                success = true;
                break;
            end
        end

        if ~success
            computation_stats.cholesky_failures = computation_stats.cholesky_failures + 1;
            if params.allow_pinv_fallback
                if params.verbose
                    fprintf('  Using pseudoinverse for freq %d (diagnostic only)\n', f);
                end
                logdet_grad = -pinv(Gamma);
                computation_stats.condition_numbers(f) = Inf;
                computation_stats.method_used = 'pseudoinverse_fallback';
            else
                error('module4_log_det_gradient:spd_violation', ...
                      'Matrix at freq %d is not SPD after regularization.', f);
            end
        end

        % Optional: force Hermitian projection of the gradient
        if params.force_hermitian
            Gherm = (logdet_grad + logdet_grad')/2;
            logdet_grad = Gherm;
        end

        logdet_gradients{f} = logdet_grad;

    catch ME
        error('module4_log_det_gradient:computation_failed', ...
             'Failed to compute gradient for frequency %d: %s', f, ME.message);
    end
end

% -------- Stats --------
computation_stats.computation_time = toc(t_total);
finite_mask = isfinite(computation_stats.condition_numbers);
if any(finite_mask)
    computation_stats.mean_condition_number = mean(computation_stats.condition_numbers(finite_mask));
    computation_stats.max_condition_number  = max(computation_stats.condition_numbers(finite_mask));
else
    computation_stats.mean_condition_number = NaN;
    computation_stats.max_condition_number  = NaN;
end

if params.verbose && computation_stats.cholesky_failures > 0
    fprintf('  Log-det gradient: %d/%d Cholesky failures\n', ...
        computation_stats.cholesky_failures, F);
end
end
