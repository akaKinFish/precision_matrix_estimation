function objective_value = module5_objective(Gamma_cells, Sigma_tilde, aux_data, params)
% MODULE5_OBJECTIVE - Evaluate complete objective function for monitoring
%
% Syntax:
%   objective_value = module5_objective(Gamma_cells, Sigma_tilde, aux_data, params)
%
% Description:
%   Evaluates the complete objective function for the sparse precision 
%   matrix estimation problem. Used for convergence monitoring and 
%   backtracking line search.
%   
%   Objective: F̃ = Σ_ω[-log det(Γ̃_ω) + tr(Σ̃_ω Γ̃_ω)] + λ₁||smoothing||² + λ₂||Γ̃_ω||_L1
%
% Input Arguments:
%   Gamma_cells - (cell array, Fx1) Precision matrices {Γ̃_ω}
%   Sigma_tilde - (cell array, Fx1) Whitened covariances {Σ̃_ω}
%   aux_data    - (struct) Auxiliary data:
%     .smoothing_kernel - (double array, FxF) Kernel matrix k_{ω,ω'}
%     .weight_matrix   - (double array, pxp) Weight matrix W^Γ
%     .lambda1         - (double) Smoothing parameter
%     .lambda2         - (double) L1 penalty parameter
%   params      - (struct) Parameters:
%     .mode            - (string) 'simplified'|'joint' (default: 'simplified')
%     .penalize_diagonal - (logical) Include diagonal in L1 (default: false)
%     .safe_logdet     - (logical) Use safe log-det computation (default: true)
%
% Output Arguments:
%   objective_value - (double) Total objective function value
%
% Examples:
%   % Basic objective evaluation
%   aux_data.smoothing_kernel = K_matrix;
%   aux_data.weight_matrix = W_matrix;
%   aux_data.lambda1 = 0.01;
%   aux_data.lambda2 = 0.001;
%   F_val = module5_objective(Gamma_cells, Sigma_cells, aux_data, struct());
%   
%   % Monitor objective during optimization
%   objectives = zeros(max_iter, 1);
%   for iter = 1:max_iter
%       objectives(iter) = module5_objective(Gamma_current, Sigma_tilde, aux_data, params);
%       % ... perform proximal update ...
%   end
%   figure; plot(objectives); title('Objective Convergence');
%
% Mathematical Background:
%   The objective combines three terms:
%   1. Negative log-likelihood: -log det(Γ) + tr(ΣΓ)
%   2. Cross-frequency smoothing: λ₁ Σ_{ω,ω'} k_{ω,ω'} ||Γ_ω - Γ_{ω'}||²_{W^Γ}
%   3. Sparsity penalty: λ₂ ||Γ||_L1 (off-diagonal elements)
%
% See also: MODULE5_PROXIMAL_MAIN, MODULE4_OBJECTIVE_EVALUATION
%
% Author: [Your Name]
% Date: [Current Date]
% Version: 1.0

% ==================== Input Validation ====================
if nargin < 3
    error('module5_objective:insufficient_input', ...
          'Gamma_cells, Sigma_tilde, and aux_data are required');
end

if nargin < 4
    params = struct();
end

% Set default parameters
if ~isfield(params, 'mode'), params.mode = 'simplified'; end
if ~isfield(params, 'penalize_diagonal'), params.penalize_diagonal = false; end
if ~isfield(params, 'safe_logdet'), params.safe_logdet = true; end

% Validate dimensions
F = length(Gamma_cells);
if length(Sigma_tilde) ~= F
    error('module5_objective:length_mismatch', ...
          'Gamma_cells and Sigma_tilde must have same length');
end

p = size(Gamma_cells{1}, 1);

% Extract auxiliary data
K_smooth = aux_data.smoothing_kernel;
W_matrix = aux_data.weight_matrix;
lambda1 = aux_data.lambda1;
lambda2 = aux_data.lambda2;

% ==================== Term 1: Negative Log-Likelihood ====================
loglik_term = 0;

for f = 1:F
    Gamma_f = Gamma_cells{f};
    Sigma_f = Sigma_tilde{f};
    
    % Validate positive definiteness
    [isPSD, ~] = module5_psd_check(Gamma_f);
    if ~isPSD
        warning('module5_objective:not_psd', ...
                'Gamma_cells{%d} is not positive definite', f);
        objective_value = Inf;
        return;
    end
    
    % Compute -log det(Γ_f)
    if params.safe_logdet
        try
            [L, chol_flag] = chol(Gamma_f, 'lower');
            if chol_flag == 0
                log_det = 2 * sum(log(real(diag(L))));
            else
                log_det = -Inf;  % Not PSD
            end
        catch
            log_det = -Inf;
        end
    else
        % Direct computation (less safe)
        log_det = log(det(Gamma_f));
    end
    
    if ~isfinite(log_det)
        warning('module5_objective:infinite_logdet', ...
                'Log-determinant is not finite at frequency %d', f);
        objective_value = Inf;
        return;
    end
    
    % Compute tr(Σ_f Γ_f)
    trace_term = real(trace(Sigma_f * Gamma_f));
    
    % Add to log-likelihood
    loglik_term = loglik_term + (-log_det + trace_term);
end

% ==================== Term 2: Cross-Frequency Smoothing ====================
smoothing_term = 0;

if lambda1 > 0
    for f1 = 1:F
        for f2 = 1:F
            if K_smooth(f1, f2) ~= 0
                Gamma1 = Gamma_cells{f1};
                Gamma2 = Gamma_cells{f2};
                diff_matrix = Gamma1 - Gamma2;

                % ||Γ₁ - Γ₂||²_{W^Γ} = Re tr((Γ₁-Γ₂)^H W (Γ₁-Γ₂))
                weighted_norm_sq = real(trace(diff_matrix' * W_matrix * diff_matrix));

                % 关键：乘 0.5，避免对称配对被重复计数
                smoothing_term = smoothing_term + 0.5 * K_smooth(f1, f2) * weighted_norm_sq;
            end
        end
    end
    smoothing_term = lambda1 * smoothing_term;
end

% ==================== Term 3: L1 Sparsity Penalty ====================
l1_term = 0;

if lambda2 > 0
    for f = 1:F
        Gamma_f = Gamma_cells{f};

        if params.penalize_diagonal
            % 对角也惩罚：上三角×2 + 对角
            l1_term = l1_term + 2*sum(sum(abs(triu(Gamma_f,1)))) + sum(abs(real(diag(Gamma_f))));
        else
            % 默认只惩罚非对角：等价于上三角×2
            l1_term = l1_term + 2*sum(sum(abs(triu(Gamma_f,1))));
        end
    end
    l1_term = lambda2 * l1_term;
end

% ==================== Combine All Terms ====================
objective_value = loglik_term + smoothing_term + l1_term;

% ==================== Validity Check ====================
if ~isfinite(objective_value)
    warning('module5_objective:non_finite_objective', ...
            'Objective function is not finite (loglik=%.2e, smooth=%.2e, l1=%.2e)', ...
            loglik_term, smoothing_term, l1_term);
end

end