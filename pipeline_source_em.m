%% PIPELINE_SOURCE_EM  (Source-domain EM pipeline, scalar estep_in version)
% Build estep_in as a *scalar struct* with per-frequency matrices stored as {F×1} cells.

%% 0) Simulation (Module 7)
n  = 128;  p  = 18;  F  = 16;  T  = 4096;

[Omega_true, Sigma_true, emp_covariance, sim] = module7_simulation_improved_complex( ...
    'n_nodes', n, 'n_sensors', p, 'n_freq', F, 'n_samples', T, ...
    'generate_leadfield', true, ...
    'leadfield_type', 'simple', ...
    'random_seed', 42);

L          = sim.leadfield_matrix;   % p×n
Sigma_xixi = sim.Sigma_xixi;         % p×p

%% 1) Build E-step input as a *scalar struct* (NOT a struct array)
emp_cov_cell = coerce_cov_cell(emp_covariance);     % supports cell / 3D / single
F            = numel(emp_cov_cell);
prior_cell   = repmat({eye(n)}, F, 1);

estep_in = struct();                                  % scalar struct
estep_in.leadfield_matrix         = L;                % p×n
estep_in.empirical_covariances    = emp_cov_cell;     % {F×1}, p×p
estep_in.source_prior_covariances = prior_cell;       % {F×1}, n×n
estep_in.noise_covariance         = Sigma_xixi;       % p×p
estep_in.frequencies              = (1:F);            % 1×F

assert(isstruct(estep_in) && isscalar(estep_in), 'estep_in must be a scalar struct');
assert(isa(estep_in.empirical_covariances,'cell') && numel(estep_in.empirical_covariances)==F);
assert(all(cellfun(@(A) isequal(size(A),[p p]), estep_in.empirical_covariances)));
assert(isa(estep_in.source_prior_covariances,'cell') && numel(estep_in.source_prior_covariances)==F);
assert(all(cellfun(@(S) isequal(size(S),[n n]), estep_in.source_prior_covariances)));

%% 2) E-step (Module 2 + wrapper)
estep_out = module2_estep(estep_in, struct( ...
    'ensure_hermitian', true, ...
    'ensure_real_diag', true, ...
    'ensure_psd',       false, ...
    'psd_tol',          1e-10, ...
    'diag_loading',     0 ...
));
Sjj_hat = estep_out.source_second_moments;   % {F×1}, n×n

%% 3) Preprocessing / Whitening in source domain (Module 1 wrapper)
pre = module1_preproc_from_covset(Sjj_hat, struct( ...
    'smoothing_method','moving_average', ...
    'loading_factor',  1e-6, ...
    'min_power',       1e-10, ...   % <- fixed name here
    'verbose',         false ...
));
D_src     = pre.D;            % {F×1}, n×n
Sjj_tilde = pre.Sigma_tilde;  % {F×1}, n×n

%% 4) Active set (Module 3)

input_data_m3 = struct();
input_data_m3.whitened_covariances = Sjj_tilde;  % 16x1 cell
input_data_m3.frequencies          = 1:F;        % 1xF double

act = module3_active_set(input_data_m3, struct( ...
        'proxy_method',   'correlation', ...
        'quantile_level', 0.25, ...
        'force_diagonal_active', true, ...
        'verbose', false ...
    ));
A_masks = arrayfun(@(f) logical(act.combined_active_mask(:,:,f)), 1:F, 'uni', 0);

%% 5) Hyperparameters (Module 6)
K = make_frequency_kernel(F, 3.0);
W = make_uniform_weight(n);
input_data_m6 = struct();
input_data_m6.whitened_covariances = Sjj_tilde;  % 16x1 cell
input_data_m6.kernel_matrix          = K;        % 1xF double
input_data_m6.weight_matrix          = W;        
input_data_m6.active_set_mask          = {A_masks{:}};        


hp = module6_hyperparameter_config(input_data_m6, struct('use_gershgorin', true));

lambda1 = hp.lambda1;  lambda2 = hp.lambda2_suggested;  alpha = hp.alpha;

%% 6) Proximal solver (Module 5) in whitened domain
Gamma0 = cellfun(@(S) diag(1 ./ max(real(diag(S)), 1e-8)), Sjj_tilde, 'uni', 0);

input_data_m5 = struct();
input_data_m5.whitened_covariances = Sjj_tilde;  % 16x1 cell
input_data_m5.initial_precision          = Gamma0;        % 1xF double
input_data_m5.smoothing_kernel          = K;        
input_data_m5.weight_matrix          = W;   
input_data_m5.active_set_mask          = {A_masks{:}};        


[Gamma_tilde_star, prox_res] = module5_proximal(input_data_m5, struct( ...
        'lambda1',         lambda1, ...
        'lambda2',         lambda2, ...
        'alpha',           alpha, ...
        'beta_backtrack',  0.5, ...
        'max_backtrack',   30, ...
        'max_iter',        300, ...
        'verbose',         true, ...
        'weight_mode',     'hadamard', ...
        'use_graph_laplacian', true ...
    ));

%% 7) Recolor to source scale (Module 8)
input_data_m8 = struct();
input_data_m8.whitened_precision_matrices = Gamma_tilde_star;  % 16x1 cell
input_data_m8.whitening_matrices          = D_src;        % 1xF double
recol = module8_recoloring(input_data_m8, struct());
Omega_src = recol.recolored_precision_matrices;   % {F×1}, n×n

%% 8) Simple readout (optional)
f_view = 1;
Om     = Omega_src{f_view};
pcorr  = abs(-Om) ./ sqrt((abs(diag(Om))+eps) * (abs(diag(Om))+eps)');
pcorr(1:n+1:end) = 0;
fprintf('Done. Example partial coherence at f=%d: max=%g, median=%g\n', ...
    f_view, max(pcorr(:)), median(pcorr(pcorr>0)));

viz_pipeline_summary(prox_res, Gamma_tilde_star, Sjj_tilde, K, A_masks, Omega_src, Omega_true);
metrics = gt_compare_and_plot(Omega_true, Omega_src, struct('mode','match_sparsity', 'f_view', 1, 'show_pr', true));
%% ===== Local helpers =====
function C = coerce_cov_cell(X, F_hint)
% Turn X into {F×1} cell of square matrices.
    if isa(X,'cell'), C = X(:); return; end
    if isnumeric(X) && ndims(X)==3 && size(X,1)==size(X,2)
        F = size(X,3); C = cell(F,1);
        for f = 1:F, C{f} = X(:,:,f); end, return;
    end
    if isnumeric(X) && ismatrix(X) && size(X,1)==size(X,2)
        if nargin>=2 && ~isempty(F_hint), C = repmat({X}, F_hint, 1);
        else, C = {X}; end, return;
    end
    error('coerce_cov_cell:unsupported','Expect cell{F,1} | p×p×F | single p×p.');
end

function K = make_frequency_kernel(F, sigma)
% Gaussian kernel along frequency index (1..F), row-normalized
    if nargin < 2, sigma = 3.0; end
    [I,J] = ndgrid(1:F,1:F);
    K = exp(-((I-J).^2)/(2*sigma^2));
    K = K ./ max(sum(K,2), eps);
end

function W = make_uniform_weight(n)
% Ones off-diagonal, zeros on diagonal
    W = ones(n); W(1:n+1:end) = 0;
end
