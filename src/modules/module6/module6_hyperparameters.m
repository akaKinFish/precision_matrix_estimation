function output = module6_hyperparameters(input, varargin)
% MODULE6_HYPERPARAMETERS - Automatic hyperparameter configuration using safe bounds.
%
% Description:
%   Computes λ₁ (smoothing) and α (step size) using:
%     - L_logdet: max spectral norm of Γ^{-1} (log-det Hessian block bound)
%     - K_max    : max kernel row sum
%     - R_max    : Gershgorin off-diagonal row-sum bound of W^Γ
%   Then sets:
%     λ₁ = δ / (2 K_max R_max),  α = 1 / (L_logdet + δ)
%   with safeguards for degenerate cases.
%
% Important:
%   - 'use_gershgorin' ONLY controls how L_logdet is estimated:
%       true  : Gershgorin bound on ||Γ^{-1}||₂ (conservative)
%       false : exact via 1 / λ_min(Σ) (no explicit matrix inverse)
%   - R_max always uses Gershgorin off-diagonal row sums (diag ignored).
%
% Inputs (fields in 'input'):
%   .whitened_covariances : {F×1 cell} Σ̃_ω (Hermitian, size p×p)
%   .kernel_matrix        : [F×F double] smoothing kernel K
%   .weight_matrix        : [p×p double] edge weights W^Γ (diag ignored)
%   .active_set_mask      : {F×1 cell} optional logical masks (p×p)
%
% Name-Value:
%   'safety_margin'        δ ∈ (0,1), default 0.9
%   'lambda2_scale'        >0, default 0.01
%   'initialization_method' {'inverse','pseudoinverse','diagonal'}, default 'inverse'
%   'verbose'              logical, default true
%   'use_gershgorin'       logical, default true
%
% Outputs:
%   output.lambda1, output.alpha, output.lambda2_suggested
%   output.diagnostics (L_logdet, K_max, R_max, etc.)
%   output.computation_info, output.quality_assessment
%
% Author: Precision Matrix Project Team
% Version: 1.1

%% Parse & validate
t0 = tic;
if ~isstruct(input)
    error('module6_hyperparameters:invalid_input','Input must be a structure');
end
req = {'whitened_covariances','kernel_matrix','weight_matrix'};
for i=1:numel(req)
    if ~isfield(input, req{i})
        error('module6_hyperparameters:missing_field','Missing field "%s"', req{i});
    end
end

P = inputParser;
addParameter(P,'safety_margin',0.9,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);
addParameter(P,'lambda2_scale',0.01,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(P,'initialization_method','inverse',@(x)ischar(x)||isstring(x));
addParameter(P,'verbose',true,@islogical);
addParameter(P,'use_gershgorin',true,@islogical);
parse(P,varargin{:});

delta          = P.Results.safety_margin;
lambda2_scale  = P.Results.lambda2_scale;
init_method    = char(P.Results.initialization_method);
verbose        = P.Results.verbose;
use_gershgorin = P.Results.use_gershgorin;

Sigma = input.whitened_covariances;    % cell F×1
K     = input.kernel_matrix;           % F×F
Wg    = input.weight_matrix;           % p×p

F = numel(Sigma);
if F==0, error('module6_hyperparameters:empty_cov','Empty whitened_covariances'); end
[p1,p2] = size(Sigma{1});
if p1~=p2, error('module6_hyperparameters:dim','First covariance must be square'); end
p = p1;

% Check all Σ: size, finite, Hermitian (up to tolerance)
tolH = 1e-12;
for f=1:F
    A = Sigma{f};
    if ~isnumeric(A) || ndims(A)~=2
        error('module6_hyperparameters:invalid_cov','Covariance %d not numeric 2D', f);
    end
    [r,c] = size(A);
    if r~=p || c~=p
        error('module6_hyperparameters:dim_mismatch','Cov %d size %dx%d ≠ %dx%d', f, r,c,p,p);
    end
    if ~all(isfinite(A(:)))
        error('module6_hyperparameters:nonfinite','Covariance %d contains NaN/Inf', f);
    end
    % Hermitian check with relative tolerance; symmetrize if needed
    relErr = norm(A - A','fro')/max(1, norm(A,'fro'));
    if relErr > tolH
        warning('module6_hyperparameters:not_hermitian', ...
            'Covariance %d not Hermitian (rel.err=%.2e). Symmetrizing.', f, relErr);
        Sigma{f} = (A + A')/2;
    end
end

if size(K,1)~=F || size(K,2)~=F
    error('module6_hyperparameters:K_size','K must be F×F'); 
end
if size(Wg,1)~=p || size(Wg,2)~=p
    error('module6_hyperparameters:W_size','W^Γ must be p×p'); 
end

if verbose
    fprintf('Module 6: Hyperparameter Configuration\n');
    fprintf('=====================================\n');
    fprintf('Problem: p=%d nodes, F=%d frequencies\n', p, F);
    fprintf('Safety margin δ=%.2f\n', delta);
    fprintf('L_logdet method: %s\n', tern(use_gershgorin, ...
        'Gershgorin bound on ||Γ^{-1}||₂', 'exact via 1/λ_min(Σ)'));
end

%% Step 1: L_logdet
if verbose, fprintf('\nStep 1: Computing L_logdet...\n'); end
conds = zeros(F,1);
gersh_bounds = zeros(F,1);
Lvals = zeros(F,1);

Gamma_init = [];
if use_gershgorin
    Gamma_init = cell(F,1);
    for f=1:F
        A = Sigma{f};
        reg = 1e-10 * trace(real(A))/p;
        Areg = A + reg*eye(p);
        switch lower(init_method)
            case 'inverse'
                Gamma_init{f} = inv(Areg);
            case 'pseudoinverse'
                Gamma_init{f} = pinv(Areg);
            case 'diagonal'
                d = max(real(diag(Areg)), eps);
                Gamma_init{f} = diag(1./d);
            otherwise
                error('module6_hyperparameters:invalid_init_method','Unknown init: %s', init_method);
        end
    end
end

for f=1:F
    A = Sigma{f};
    conds(f) = cond(A);
    if use_gershgorin
        G = Gamma_init{f};
        dg = diag(G);
        offsum = sum(abs(G),2) - abs(dg);
        bnd = max(abs(dg) + offsum);
        gersh_bounds(f) = bnd;
        Lvals(f) = bnd^2;
    else
        reg = 1e-10 * trace(real(A))/p;
        Areg = (A + A')/2 + reg*eye(p);
        lam_min = max(min(eig(Areg,'vector')), eps);
        Lvals(f) = (1/lam_min)^2;
    end
end
L_logdet = max(Lvals);

if verbose
    fprintf('  L_logdet = %.4e\n', L_logdet);
    fprintf('  mean cond(Σ) = %.2e\n', mean(conds));
    if use_gershgorin
        fprintf('  Gershgorin ||Γ^{-1}||₂ bound range: [%.2e, %.2e]\n', min(gersh_bounds), max(gersh_bounds));
    end
end

%% Step 2: K_max
if verbose, fprintf('\nStep 2: Computing K_max...\n'); end
if ~issymmetric(K)
    if verbose
        warning('module6_hyperparameters:K_not_sym','K not symmetric; symmetrizing');
    end
    K = (K + K')/2;
end
K_rows = sum(abs(K),2);
K_max = max(K_rows);
if verbose
    fprintf('  K_max = %.4f (row-sum)\n', K_max);
end

%% Step 3: R_max
if verbose, fprintf('\nStep 3: Computing R_max...\n'); end
if ~issymmetric(Wg)
    if verbose
        warning('module6_hyperparameters:W_not_sym','W^Γ not symmetric; symmetrizing');
    end
    Wg = (Wg + Wg')/2;
end
off_row_sums = sum(abs(Wg),2) - abs(diag(Wg)); % diag ignored
R_max = max(off_row_sums);
if verbose
    fprintf('  R_max = %.4f (off-diagonal row-sum)\n', R_max);
end

%% Step 4: λ₁ and α
if verbose, fprintf('\nStep 4: Computing hyperparameters...\n'); end
if K_max * R_max < eps
    warning('module6_hyperparameters:zero_product','K_max*R_max≈0; using fallback λ₁');
    lambda1 = 1e-6;
else
    lambda1 = delta / (2 * K_max * R_max);
end
alpha = 1 / (max(L_logdet, eps) + delta);

% λ₂ suggestion
lambda2_policy = '';
active_density = NaN;
if isfield(input,'active_set_mask') && ~isempty(input.active_set_mask) ...
        && iscell(input.active_set_mask) && numel(input.active_set_mask)==F
    tot_act = 0; tot_poss = 0; ok = true;
    for f=1:F
        M = input.active_set_mask{f};
        if ~isnumeric(M) && ~islogical(M), ok=false; break; end
        [r,c] = size(M);
        if r~=p || c~=p, ok=false; break; end
        M = logical(M);
        tot_act = tot_act + sum(M(triu(true(p),1)));
        tot_poss = tot_poss + p*(p-1)/2;
    end
    if ok && tot_poss>0
        active_density = tot_act / tot_poss;
        lambda2_suggested = lambda2_scale * active_density;
        lambda2_policy = 'density_scaled';
    else
        lambda2_suggested = lambda2_scale * sqrt(log(p)/p);
        lambda2_policy = 'size_scaled';
    end
else
    lambda2_suggested = lambda2_scale * sqrt(log(p)/p);
    lambda2_policy = 'size_scaled';
end

if verbose
    fprintf('  λ₁ = %.6e, α = %.6e, λ₂(sugg) = %.6e\n', lambda1, alpha, lambda2_suggested);
end

%% Step 5: Diagnostics & output
theo_rate = alpha * L_logdet;     % should be < 1
ker_eigs = eig((K+K')/2);
ker_rank = sum(ker_eigs > 1e-10 * max(ker_eigs));
W_connect = (nnz(Wg) - p) / max(1, (p^2 - p));

stability_ok = alpha>0 && lambda1>0 && isfinite(alpha) && isfinite(lambda1);
if ~stability_ok
    error('module6_hyperparameters:stability','Unstable α=%.2e or λ₁=%.2e', alpha, lambda1);
end

output = struct();
output.lambda1 = lambda1;
output.alpha = alpha;
output.lambda2_suggested = lambda2_suggested;

output.diagnostics = struct();
output.diagnostics.L_logdet = L_logdet;
output.diagnostics.K_max = K_max;
output.diagnostics.R_max = R_max;
output.diagnostics.L_logdet_per_frequency = Lvals;
output.diagnostics.condition_numbers = conds;
output.diagnostics.gershgorin_bounds = gersh_bounds;
output.diagnostics.kernel_row_sums = K_rows;
output.diagnostics.weight_row_sums = off_row_sums;
output.diagnostics.theoretical_convergence_rate = theo_rate;
output.diagnostics.kernel_effective_rank = ker_rank;
output.diagnostics.weight_connectivity = W_connect;
output.diagnostics.lambda2_policy = lambda2_policy;
output.diagnostics.active_density = active_density;
output.diagnostics.stability_check = stability_ok;

output.computation_info = struct();
output.computation_info.total_time_seconds = toc(t0);
output.computation_info.initialization_method = init_method;
output.computation_info.eigenvalue_method = tern(use_gershgorin,'gershgorin','exact-min-eig');
output.computation_info.safety_margin_used = delta;
output.computation_info.n_nodes = p;
output.computation_info.n_frequencies = F;

qa = struct();
qa.lambda1_reasonable = (lambda1 > 1e-8) && (lambda1 < 1);
qa.alpha_reasonable   = (alpha   > 1e-6) && (alpha   < 1);
qa.condition_acceptable = all(conds < 1e12);
qa.overall_quality = qa.lambda1_reasonable && qa.alpha_reasonable && qa.condition_acceptable;
output.quality_assessment = qa;

if verbose
    fprintf('\nModule 6 finished in %.3f s\n', output.computation_info.total_time_seconds);
    fprintf('Quality: %s\n', tern(qa.overall_quality,'GOOD','NEEDS ATTENTION'));
end
end

%% Helper
function s = tern(c,a,b), if c, s=a; else, s=b; end, end
