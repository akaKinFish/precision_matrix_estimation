function output = module6_hyperparameters(input, varargin)
% MODULE6_HYPERPARAMETERS - Automatic hyperparameter configuration.
%
% Mathematical Logic (Corrected):
%   1. Step Size (alpha): Depends on the Lipschitz constant L of the gradient.
%      L_total <= L_logdet + lambda1 * L_smooth
%      
%      L_logdet = || nabla^2 (-log det G) ||_2 = || G^{-1} ||_2^2
%      At initialization G_0 approx Sigma^{-1}, so G^{-1} approx Sigma.
%      Thus, L_logdet approx || Sigma ||_2^2 = lambda_max(Sigma)^2.
%
%   2. Smoothing (lambda1):
%      L_smooth <= 2 * K_max * R_max (Gershgorin bound of Laplacian term).
%      We set lambda1 such that lambda1 * L_smooth <= delta (safety margin).
%
%   3. Sparsity (lambda2):
%      Heuristic suggestion based on high-dimensional stats theory.
%
% Inputs:
%   .whitened_covariances : {F×1 cell} Sigma_tilde
%   .kernel_matrix        : [F×F] K
%   .weight_matrix        : [p×p] W
%
% Name-Value Pairs:
%   'safety_margin'  : delta (default 0.9)
%   'use_gershgorin' : true/false (default true). If true, uses Row-Sum bounds.
%                      If false, uses exact eigenvalues (slower but tighter).

%% 1. Parse & Validate
t0 = tic;
P = inputParser;
addParameter(P, 'safety_margin', 0.9, @isnumeric);
addParameter(P, 'lambda2_scale', 0.01, @isnumeric);
addParameter(P, 'use_gershgorin', true, @islogical);
addParameter(P, 'verbose', true, @islogical);
parse(P, varargin{:});

delta = P.Results.safety_margin;
use_gersh = P.Results.use_gershgorin;
verbose = P.Results.verbose;

Sigma = input.whitened_covariances;
K = input.kernel_matrix;
W = input.weight_matrix;

F = numel(Sigma);
p = size(Sigma{1}, 1);

% Ensure Symmetry of K and W
K = (K + K')/2;
W = (W + W')/2;

if verbose
    fprintf('Module 6: Hyperparameter Configuration\n');
    fprintf('=====================================\n');
    fprintf('Method: %s\n', ternary(use_gersh, 'Gershgorin Bounds (Fast, Conservative)', 'Exact Eigenvalues (Slow, Tight)'));
end

%% 2. Compute L_logdet (Curvature of Data Term)
% Theory: L_logdet approx lambda_max(Sigma)^2
if verbose, fprintf('Step 1: Computing L_logdet (Data Curvature)...\n'); end

L_vals = zeros(F, 1);

for f = 1:F
    S = Sigma{f};
    % Ensure Hermitian
    S = (S + S')/2;
    
    if use_gersh
        % Gershgorin Upper Bound for lambda_max(S)
        % Bound = max_i ( S_ii + sum_{j!=i} |S_ij| )
        % This is simply max(sum(abs(S), 2))
        row_sums = sum(abs(S), 2);
        lambda_max_bound = max(row_sums);
        
        % L = bound^2
        L_vals(f) = lambda_max_bound^2;
    else
        % Exact Spectral Norm
        % Use eig or svd. S is Hermitian, so eig is fine.
        % We want max eigenvalue.
        % 'approx' via eigs is faster for large p, but let's use eig for robustness here unless p is huge
        if p > 2000
            opts.issym = 1;
            l_max = eigs(S, 1, 'LA', opts);
        else
            l_max = max(eig(S));
        end
        L_vals(f) = l_max^2;
    end
end

L_logdet = max(L_vals);
if verbose, fprintf('  L_logdet = %.4e\n', L_logdet); end

%% 3. Compute L_smooth Constants (K_max, R_max)
% Theory: L_smooth <= 2 * ||K||_inf * ||W||_inf
if verbose, fprintf('Step 2: Computing Smoothing Bounds...\n'); end

% K_max: Max row sum of Kernel
K_max = max(sum(abs(K), 2));

% R_max: Max row sum of Weight Matrix (Gershgorin radius)
% Usually W has 0 on diagonal, but formula works regardless.
R_max = max(sum(abs(W), 2));

if verbose
    fprintf('  K_max = %.4f\n', K_max);
    fprintf('  R_max = %.4f\n', R_max);
end

%% 4. Calculate Lambda1 and Alpha
if verbose, fprintf('Step 3: Calculating Parameters...\n'); end

% Lambda1 Strategy:
% Limit the smoothing curvature contribution to 'delta'.
% lambda1 * (2 * K_max * R_max) = delta
denominator = 2 * K_max * R_max;

if denominator < 1e-12
    % No smoothing impact (e.g. K=0 or W=0)
    lambda1 = 0; 
    if verbose, fprintf('  Smoothing disabled (K or W is zero).\n'); end
else
    lambda1 = delta / denominator;
end

% Alpha Strategy:
% alpha <= 1 / L_total
% L_total = L_logdet + lambda1 * L_smooth
% With our choice of lambda1, lambda1 * L_smooth <= delta.
% So alpha = 1 / (L_logdet + delta).

alpha = 1 / (L_logdet + delta);

%% 5. Lambda2 Suggestion (Heuristic / Stats)
% We provide a statistical baseline: sqrt(log(p)/n) type scaling.
% Here we assume n is implicitly related to p or provided roughly.
% A robust baseline for unit-variance data is sqrt(2*log(p)/p) if n=p.
% Let's stick to the provided heuristic scale.

lambda2_suggested = P.Results.lambda2_scale * sqrt(log(p)/p);

if verbose
    fprintf('  Result: lambda1 = %.4e\n', lambda1);
    fprintf('  Result: alpha   = %.4e\n', alpha);
    fprintf('  Result: lambda2 = %.4e (Heuristic)\n', lambda2_suggested);
end

%% 6. Output
output.lambda1 = lambda1;
output.alpha = alpha;
output.lambda2_suggested = lambda2_suggested;

output.diagnostics.L_logdet = L_logdet;
output.diagnostics.K_max = K_max;
output.diagnostics.R_max = R_max;
output.diagnostics.method = ternary(use_gersh, 'gershgorin', 'exact');

end

function val = ternary(cond, a, b)
    if cond, val = a; else, val = b; end
end