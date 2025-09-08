function estep_results = module2_estep_main(input_data, estep_params)
% MODULE2_ESTEP_MAIN - Complete E-step computation for EM algorithm (Full Source Model)
%
% CHANGELOG (2025-09):
%   * initial_precision_matrices are NOW in the SOURCE domain (n×n), computed as a
%     robust inverse of posterior source covariances Σ_{jj,post}^{(f)}.
%   * source_second_moments are returned:  Ŝ_{jj}^{(f)} = T_{jv}^{(f)} S_{vv}^{(f)} T_{jv}^{(f)H} + Σ_{jj,post}^{(f)}.
%   * For compatibility, residual-domain initial precisions ( (S_{ξξ}^{(f)})^{-1}, p×p ) are
%     provided under residual_initial_precision_matrices (NEW NAME).
%
% Inputs (fields of input_data):
%   .leadfield_matrix         (double, p×n)  : L
%   .empirical_covariances    (cell{F} or p×p): {S_vv^(f)}
%   .source_prior_covariances (cell{F} or n×n): {Σ_jj^(f)}
%   .noise_covariance         (complex, p×p)  : Σ_ξξ
%   .frequencies              (double, F×1 or 1×F)
%
% Name-Value in estep_params:
%   .regularization_factor (double) default 1e-8
%   .condition_threshold   (double) default 1e12
%   .min_eigenvalue_ratio  (double) default 1e-12
%   .verbose               (logical) default false
%
% Output (struct):
%   .transfer_functions                   {F}  (n×p)   : T_{jv}^{(f)}
%   .posterior_source_covariances         {F}  (n×n)   : Σ_{jj,post}^{(f)}
%   .source_second_moments                {F}  (n×n)   : Ŝ_{jj}^{(f)}
%   .residual_transfer_functions          {F}  (p×p)   : T_{ξv}^{(f)}
%   .residual_covariances                 {F}  (p×p)   : S_{ξξ}^{(f)}
%   .initial_precision_matrices           {F}  (n×n)   : Ω_{src,init}^{(f)} = inv_robust(Σ_{jj,post}^{(f)})
%   .residual_initial_precision_matrices  {F}  (p×p)   : Ω_{resid,init}^{(f)} = inv_robust(S_{ξξ}^{(f)})  [compat only]
%   .computation_stats / .success

% -------------------- Minimal input presence ----------------------------
if nargin < 1
    error('module2_estep_main:insufficient_input','input_data is required');
end
if nargin < 2
    estep_params = struct();
end

req = {'leadfield_matrix','empirical_covariances','source_prior_covariances','noise_covariance','frequencies'};
for k = 1:numel(req)
    if ~isfield(input_data, req{k})
        error('module2_estep_main:missing_field','Missing input_data.%s', req{k});
    end
end

% -------------------- Extract & robustify --------------------------------
L                = input_data.leadfield_matrix;        % p×n
Sigma_emp        = input_data.empirical_covariances;   % {F} or p×p
Sigma_prior_in   = input_data.source_prior_covariances;% {F} or n×n
Sigma_xi_xi      = input_data.noise_covariance;        % p×p
frequencies      = input_data.frequencies;

% Auto-wrap to cell for single-frequency usage
if ~iscell(Sigma_emp),      Sigma_emp      = {Sigma_emp};      end
if ~iscell(Sigma_prior_in), Sigma_prior_in = {Sigma_prior_in}; end

% Basic numeric checks
if ~isnumeric(L) || ndims(L) ~= 2
    error('module2_estep_main:invalid_leadfield','leadfield_matrix must be 2D numeric');
end
[p, n] = size(L);
if ~isnumeric(Sigma_xi_xi) || ~isequal(size(Sigma_xi_xi), [p p])
    error('module2_estep_main:invalid_noise_covariance','noise_covariance must be %dx%d', p, p);
end

F = numel(Sigma_emp);

% Frequencies length fix
if ~isvector(frequencies) || numel(frequencies) ~= F
    frequencies = 1:F;
end

% Broadcast a single prior to all F if needed
if numel(Sigma_prior_in) == 1 && F > 1
    Sigma_prior_in = repmat(Sigma_prior_in, F, 1);
end
if numel(Sigma_prior_in) ~= F
    error('module2_estep_main:source_covariance_mismatch', ...
          'Number of source priors (%d) must match number of empirical covariances (%d).', ...
          numel(Sigma_prior_in), F);
end

% -------------------- Params with defaults --------------------------------
def = struct('regularization_factor',1e-8, ...
             'condition_threshold', 1e12, ...
             'min_eigenvalue_ratio',1e-12, ...
             'verbose', false);
fn = fieldnames(def);
for i = 1:numel(fn)
    if ~isfield(estep_params, fn{i}), estep_params.(fn{i}) = def.(fn{i}); end
end

% -------------------- Prepare outputs/stats --------------------------------
estep_results = struct();
estep_results.transfer_functions                    = cell(F,1); % n×p
estep_results.posterior_source_covariances          = cell(F,1); % n×n
estep_results.source_second_moments                 = cell(F,1); % n×n  (NEW)
estep_results.residual_transfer_functions           = cell(F,1); % p×p
estep_results.residual_covariances                  = cell(F,1); % p×p
estep_results.initial_precision_matrices            = cell(F,1); % n×n  (NOW SOURCE DOMAIN)
estep_results.residual_initial_precision_matrices   = cell(F,1); % p×p  (compat/diagnostics)
estep_results.computation_stats                     = struct();
estep_results.success                               = false;

stats = struct();
stats.processing_times         = zeros(F,6); % +1 slot for source-moments
stats.condition_numbers        = zeros(F,3);
stats.regularization_applied   = false(F,1);
stats.eigenvalue_corrections   = zeros(F,4); % [neg/small for Σpost,  neg/small for S_xi_xi]
stats.total_computation_time   = 0;
stats.successful_frequencies   = 0;
stats.failed_frequencies       = 0;
stats.error_messages           = {};

overall_tic = tic;
if estep_params.verbose
    fprintf('========================================\n');
    fprintf('Module 2: E-Step (Full Source Model)\n');
    fprintf('Processing %d frequencies | sensors=%d, sources=%d\n', F, p, n);
    fprintf('========================================\n\n');
end

% -------------------- Main loop -------------------------------------------
for f = 1:F
    try
        freq_tic = tic;

        % --- Validate per-frequency inputs
        S_vv_f     = Sigma_emp{f};          % p×p empirical (sensor)
        Sigma_jj_f = Sigma_prior_in{f};     % n×n prior (source)
        if ~isnumeric(S_vv_f) || ~isequal(size(S_vv_f), [p p])
            error('Empirical covariance at f=%d must be %dx%d numeric', f, p, p);
        end
        if ~isnumeric(Sigma_jj_f) || ~isequal(size(Sigma_jj_f), [n n])
            error('Source prior at f=%d must be %dx%d numeric', f, n, n);
        end

        % 1) DSTF T_{jv}
        t1 = tic;
        opt_dstf = struct('regularization_factor', estep_params.regularization_factor, ...
                          'condition_threshold',   estep_params.condition_threshold, ...
                          'verbose', false);
        T_jv_f = module2_dstf_computation(L, Sigma_jj_f, Sigma_xi_xi, opt_dstf); % n×p
        estep_results.transfer_functions{f} = T_jv_f;
        stats.processing_times(f,1) = toc(t1);

        % 2) Posterior Σ_{jj,post}
        t2 = tic;
        opt_psc = struct('regularization_factor', estep_params.regularization_factor, ...
                         'condition_threshold',   estep_params.condition_threshold, ...
                         'min_eigenvalue_ratio',  estep_params.min_eigenvalue_ratio, ...
                         'verbose', false);
        Sigma_post_f = module2_posterior_source_covariance(Sigma_jj_f, L, Sigma_xi_xi, opt_psc); % n×n
        estep_results.posterior_source_covariances{f} = Sigma_post_f;
        stats.processing_times(f,2) = toc(t2);

        % 2.5) Source second moments  Ŝ_{jj} = T_jv S_vv T_jv' + Σ_{jj,post}
        t25 = tic;
        Sjj_hat_f = T_jv_f * S_vv_f * T_jv_f' + Sigma_post_f;      % n×n
        Sjj_hat_f = (Sjj_hat_f + Sjj_hat_f')/2;                    % ensure Hermitian
        estep_results.source_second_moments{f} = Sjj_hat_f;
        stats.processing_times(f,6) = toc(t25);

        % 3) Residual TF T_{ξv}
        t3 = tic;
        T_xi_v_f = module2_residual_transfer_function(T_jv_f, L, ...
                        'validate_properties', true, ...
                        'enforce_hermitian', false, ...
                        'verbose', false);
        estep_results.residual_transfer_functions{f} = T_xi_v_f;   % p×p
        stats.processing_times(f,3) = toc(t3);

        % 4) Residual empirical covariance S_{ξξ}
        t4 = tic;
        min_thr = estep_params.min_eigenvalue_ratio * max(real(eig(S_vv_f)));
        if ~isfinite(min_thr) || min_thr < 0, min_thr = 0; end
        opt_rec = struct('regularization_factor',  estep_params.regularization_factor, ...
                         'min_eigenvalue_threshold', max(min_thr, 0), ...
                         'enforce_hermitian', true, ...
                         'verbose', false);
        S_xi_xi_f = module2_residual_empirical_covariance(T_xi_v_f, S_vv_f, opt_rec); % p×p
        estep_results.residual_covariances{f} = S_xi_xi_f;
        stats.processing_times(f,4) = toc(t4);

        % 5) INITIAL PRECISIONS (NEW POLICY)
        % 5a) SOURCE-domain initial precision  Ω_src,init = inv_robust(Σ_{jj,post})
        eps_reg = estep_params.regularization_factor;
        min_rat = estep_params.min_eigenvalue_ratio;

        [Omega_src_init_f, neg1, small1] = inv_psd_robust(Sigma_post_f, eps_reg, min_rat); % n×n
        estep_results.initial_precision_matrices{f} = Omega_src_init_f;                    % << n×n

        % 5b) (Optional for diagnostics/compat) RESIDUAL-domain initial precision
        [Omega_resid_init_f, neg2, small2] = inv_psd_robust(S_xi_xi_f, eps_reg, min_rat);  % p×p
        estep_results.residual_initial_precision_matrices{f} = Omega_resid_init_f;

        stats.eigenvalue_corrections(f,:) = [neg1, small1, neg2, small2];

        % Diagnostics (unchanged)
        try, stats.condition_numbers(f,1) = cond(L * Sigma_jj_f * L' + Sigma_xi_xi); catch, stats.condition_numbers(f,1) = NaN; end
        try, stats.condition_numbers(f,2) = cond(Sigma_post_f);                       catch, stats.condition_numbers(f,2) = NaN; end
        try, stats.condition_numbers(f,3) = cond(S_xi_xi_f);                          catch, stats.condition_numbers(f,3) = NaN; end

        stats.regularization_applied(f) = true;
        stats.successful_frequencies = stats.successful_frequencies + 1;

        if estep_params.verbose
            fprintf('  ✓ f=%d/%d (%.2f Hz): done in %.3fs\n', ...
                f, F, frequencies(min(f,numel(frequencies))), toc(freq_tic));
        end

    catch ME
        stats.failed_frequencies = stats.failed_frequencies + 1;
        stats.error_messages{end+1} = sprintf('f=%d: %s', f, ME.message);

        estep_results.transfer_functions{f}                    = [];
        estep_results.posterior_source_covariances{f}          = [];
        estep_results.source_second_moments{f}                 = []; % new
        estep_results.residual_transfer_functions{f}           = [];
        estep_results.residual_covariances{f}                  = [];
        estep_results.initial_precision_matrices{f}            = []; % now n×n
        estep_results.residual_initial_precision_matrices{f}   = [];

        if estep_params.verbose
            fprintf('  ✗ f=%d failed: %s\n', f, ME.message);
        end
    end
end

stats.total_computation_time = toc(overall_tic);
estep_results.computation_stats = stats;
estep_results.success = (stats.successful_frequencies / F) >= 0.8;

if estep_params.verbose
    fprintf('\n========================================\n');
    fprintf('E-step summary: %d/%d succeeded (%.1f%%), total %.3fs\n', ...
        stats.successful_frequencies, F, 100*stats.successful_frequencies/max(F,1), stats.total_computation_time);
    if ~isempty(stats.error_messages)
        fprintf('Errors:\n'); fprintf('  - %s\n', stats.error_messages{:});
    end
end

% ---- Soft validation for the SOURCE-domain initial precisions ------------
if estep_results.success
    good = 0;
    for f = 1:F
        Om = estep_results.initial_precision_matrices{f}; % n×n
        if ~isempty(Om)
            try
                lam_min = min(real(eig((Om+Om')/2)));
                if lam_min > -1e-10, good = good + 1; end
            catch
            end
        end
    end
    if good < stats.successful_frequencies
        warning('module2_estep_main:invalid_source_precision', ...
            'Some SOURCE-domain initial precision matrices may not be positive definite.');
    end
end
end

% ========================= helpers =========================
function [Ainv, neg_cnt, small_cnt] = inv_psd_robust(A, eps_reg, min_ratio)
% Robust inverse for (near-)Hermitian PSD matrices.
% - Floors eigenvalues at max(min_ratio * max_eig, 0)
% - Optional ridge-like blend controlled by eps_reg
    Ah = (A + A')/2;
    [V,D] = eig(Ah);
    d = real(diag(D));
    dmax = max(d);
    floor_val = max(min_ratio * max(dmax, eps), 0);
    neg_cnt   = sum(d < 0);
    small_cnt = sum(d >= 0 & d < floor_val);
    d = max(d, floor_val);
    if eps_reg > 0
        d = (d + eps_reg * dmax) / (1 + eps_reg);
    end
    Ainv = V * diag(1./d) * V';
    Ainv = (Ainv + Ainv')/2;
end
