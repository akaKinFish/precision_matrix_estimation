function debug_em_convergence(Sigma_collection, varargin)
    % DEBUG_EM_CONVERGENCE - Diagnose convergence issues in EM algorithm
    %
    % Syntax:
    %   debug_em_convergence(Sigma_collection)
    %   debug_em_convergence(Sigma_collection, Name, Value)
    %
    % Description:
    %   Comprehensive diagnostic tool to identify numerical and algorithmic
    %   issues causing convergence problems in the EM pipeline.
    %
    % Input Arguments:
    %   Sigma_collection - (p×p×F complex) Covariance matrices
    %
    % Name-Value Arguments:
    %   weight_mode - (string) Weight mode: 'hadamard', 'uniform'. Default: 'hadamard'
    %   use_laplacian - (logical) Use Laplacian weights. Default: true
    %   verbose - (logical) Verbose output. Default: true
    %
    % Output:
    %   Diagnostic results printed to console and optionally saved
    %
    % Author: EM Debug Tool
    % Date: 2025-01-09
    % Version: 1.0
    
    % Parse inputs
    p = inputParser;
    addRequired(p, 'Sigma_collection', @(x) ndims(x) == 3);
    addParameter(p, 'weight_mode', 'hadamard', @ischar);
    addParameter(p, 'use_laplacian', true, @islogical);
    addParameter(p, 'verbose', true, @islogical);
    parse(p, Sigma_collection, varargin{:});
    
    opts = p.Results;
    [n, ~, F] = size(Sigma_collection);
    
    fprintf('\n========== EM CONVERGENCE DIAGNOSTIC ==========\n');
    fprintf('Problem size: n=%d, F=%d\n', n, F);
    fprintf('Configuration: weight_mode=%s, use_laplacian=%d\n', ...
            opts.weight_mode, opts.use_laplacian);
    
    %% 1. Check Hermitian property of input matrices
    fprintf('\n[1] Input Matrix Properties:\n');
    hermitian_errors = zeros(F, 1);
    eigenvalue_ranges = zeros(F, 2);
    condition_numbers = zeros(F, 1);
    
    for f = 1:F
        Sigma_f = Sigma_collection(:, :, f);
        
        % Hermitian check
        hermitian_errors(f) = norm(Sigma_f - Sigma_f', 'fro') / norm(Sigma_f, 'fro');
        
        % Eigenvalue analysis
        eigenvals = eig(Sigma_f);
        eigenvalue_ranges(f, :) = [min(real(eigenvals)), max(real(eigenvals))];
        condition_numbers(f) = cond(Sigma_f);
    end
    
    fprintf('  Hermitian errors: max=%.2e, mean=%.2e\n', ...
            max(hermitian_errors), mean(hermitian_errors));
    fprintf('  Eigenvalue range: [%.2e, %.2e]\n', ...
            min(eigenvalue_ranges(:,1)), max(eigenvalue_ranges(:,2)));
    fprintf('  Condition numbers: max=%.2e, mean=%.2e\n', ...
            max(condition_numbers), mean(condition_numbers));
    
    %% 2. Construct and analyze smoothing kernel K
    fprintf('\n[2] Smoothing Kernel Analysis:\n');
    
    % Build weight matrix
    if opts.use_laplacian
        adjacency = ones(n) - eye(n);  % Fully connected for now
        degree = sum(adjacency, 2);
        laplacian = diag(degree) - adjacency;
        weight_matrix = laplacian;
    else
        weight_matrix = ones(n);
    end
    
    % Construct K using different methods
    K_real = zeros(n, n);
    K_abs = zeros(n, n);
    K_complex = zeros(n, n);
    
    for f = 1:F
        Sigma_f = Sigma_collection(:, :, f);
        
        if strcmp(opts.weight_mode, 'hadamard')
            K_real = K_real + weight_matrix .* real(Sigma_f);
            K_abs = K_abs + weight_matrix .* abs(Sigma_f);
            K_complex = K_complex + weight_matrix .* Sigma_f;
        else
            K_real = K_real + real(Sigma_f);
            K_abs = K_abs + abs(Sigma_f);
            K_complex = K_complex + Sigma_f;
        end
    end
    
    % Analyze symmetry
    fprintf('  K_real symmetry error: %.2e\n', ...
            norm(K_real - K_real', 'fro') / norm(K_real, 'fro'));
    fprintf('  K_abs symmetry error: %.2e\n', ...
            norm(K_abs - K_abs', 'fro') / norm(K_abs, 'fro'));
    fprintf('  K_complex Hermitian error: %.2e\n', ...
            norm(K_complex - K_complex', 'fro') / norm(K_complex, 'fro'));
    
    % Check positive definiteness
    eig_real = eig(K_real);
    eig_abs = eig(K_abs);
    fprintf('  K_real: min_eig=%.2e, positive definite=%d\n', ...
            min(eig_real), all(eig_real > 0));
    fprintf('  K_abs: min_eig=%.2e, positive definite=%d\n', ...
            min(eig_abs), all(eig_abs > 0));
    
    %% 3. Hyperparameter scale analysis
    fprintf('\n[3] Hyperparameter Scale Analysis:\n');
    
    % Compute various norms and bounds
    K = K_real;  % Use real part for analysis
    K_max = max(sum(abs(K), 2));  % Row sum norm
    R_max = max(sum(abs(K - diag(diag(K))), 2));  % Off-diagonal row sum
    
    % Estimate appropriate scales
    data_scale = mean(diag(K));
    sparsity_scale = mean(abs(K(~eye(n))));
    
    fprintf('  Data scale (mean diagonal): %.2e\n', data_scale);
    fprintf('  Sparsity scale (mean off-diag): %.2e\n', sparsity_scale);
    fprintf('  K_max (row sum): %.2e\n', K_max);
    fprintf('  R_max (off-diag row sum): %.2e\n', R_max);
    
    % Suggest hyperparameters based on different heuristics
    lambda1_theory = 1 / (2 * K_max);
    lambda1_empirical = 0.1 * data_scale;
    alpha_theory = 1 / (2 * R_max);
    alpha_empirical = 0.01 * sparsity_scale;
    
    fprintf('\n  Suggested hyperparameters:\n');
    fprintf('    Theory-based: λ₁=%.2e, α=%.2e\n', lambda1_theory, alpha_theory);
    fprintf('    Empirical: λ₁=%.2e, α=%.2e\n', lambda1_empirical, alpha_empirical);
    
    %% 4. Gradient and objective function analysis
    fprintf('\n[4] Objective Function Analysis:\n');
    
    % Test gradient at random points
    n_test = 5;
    Gamma_test = eye(n) + 0.1 * randn(n);
    Gamma_test = (Gamma_test + Gamma_test') / 2;  % Ensure symmetric
    
    % Make positive definite
    [V, D] = eig(Gamma_test);
    D = diag(max(diag(D), 0.1));
    Gamma_test = V * D * V';
    
    % Compute objective and gradient (simplified version)
    obj = compute_simple_objective(Gamma_test, K, lambda1_theory);
    grad_numerical = compute_numerical_gradient(Gamma_test, K, lambda1_theory);
    
    fprintf('  Test objective value: %.2e\n', obj);
    fprintf('  Gradient norm: %.2e\n', norm(grad_numerical, 'fro'));
    fprintf('  Gradient symmetry: %.2e\n', ...
            norm(grad_numerical - grad_numerical', 'fro') / norm(grad_numerical, 'fro'));
    
    %% 5. Recommendations
    fprintf('\n[5] RECOMMENDATIONS:\n');
    fprintf('=====================================\n');
    
    if max(hermitian_errors) > 1e-10
        fprintf('⚠ Input matrices are not perfectly Hermitian\n');
        fprintf('  → Enforce Hermitian: Sigma = (Sigma + Sigma'')/2\n');
    end
    
    if norm(K_real - K_real', 'fro') / norm(K_real, 'fro') > 1e-10
        fprintf('⚠ Kernel matrix K is not symmetric\n');
        fprintf('  → Use real part explicitly: K = real(sum(Sigma))\n');
        fprintf('  → Or symmetrize: K = (K + K'')/2\n');
    end
    
    if min(eig_real) < 0
        fprintf('⚠ Kernel matrix is not positive definite\n');
        fprintf('  → Add regularization: K = K + epsilon*I\n');
        fprintf('  → Consider different weight scheme\n');
    end
    
    ratio = lambda1_empirical / lambda1_theory;
    if ratio > 100 || ratio < 0.01
        fprintf('⚠ Hyperparameter scales differ by %.0fx\n', max(ratio, 1/ratio));
        fprintf('  → Use adaptive scaling based on data\n');
        fprintf('  → Consider logarithmic grid search\n');
    end
    
    if max(condition_numbers) > 1e3
        fprintf('⚠ High condition number detected (%.2e)\n', max(condition_numbers));
        fprintf('  → Apply diagonal loading to input matrices\n');
        fprintf('  → Use regularized inversions\n');
    end
    
    fprintf('\n========== END DIAGNOSTIC ==========\n\n');
end

function obj = compute_simple_objective(Gamma, K, lambda1)
    % Simplified objective for testing
    n = size(Gamma, 1);
    
    % Ensure Gamma is positive definite for log-det
    [V, D] = eig(Gamma);
    D = diag(max(diag(D), 1e-6));
    Gamma_reg = V * D * V';
    
    % Compute terms
    logdet_term = -log(det(Gamma_reg));
    trace_term = trace(K * Gamma);
    l1_term = lambda1 * sum(abs(Gamma(:)));
    
    obj = logdet_term + trace_term + l1_term;
end

function grad = compute_numerical_gradient(Gamma, K, lambda1)
    % Numerical gradient for verification
    n = size(Gamma, 1);
    grad = zeros(n, n);
    epsilon = 1e-8;
    
    for i = 1:n
        for j = i:n  % Only upper triangle
            % Perturb element
            Gamma_plus = Gamma;
            Gamma_plus(i, j) = Gamma_plus(i, j) + epsilon;
            if i ~= j
                Gamma_plus(j, i) = Gamma_plus(j, i) + epsilon;  % Maintain symmetry
            end
            
            % Compute finite difference
            f_plus = compute_simple_objective(Gamma_plus, K, lambda1);
            f = compute_simple_objective(Gamma, K, lambda1);
            
            grad_ij = (f_plus - f) / epsilon;
            grad(i, j) = grad_ij;
            if i ~= j
                grad(j, i) = grad_ij;  % Symmetric gradient
            end
        end
    end
end