function estep_results = module2_estep_main(input_data, estep_params)
% MODULE2_ESTEP_MAIN - Complete E-step computation for EM algorithm (Full Source Model)
%
% Syntax:
%   estep_results = module2_estep_main(input_data, estep_params)
%
% Description:
%   Performs the complete E-step across F frequencies:
%     1) DSTF  T_jv
%     2) Posterior source covariance Σ_jj,post
%     3) Residual TF T_ξv
%     4) Residual empirical covariance S_ξξ
%     5) Initial precision estimate Ω with eigenvalue regularization
%
% Robustness additions in this revision:
%   - Accepts single matrices for empirical_covariances/source_prior_covariances by
%     wrapping them into 1×1 cell arrays automatically.
%   - If frequencies length does not match F, it is rebuilt as 1:F.
%   - If only one source prior covariance is provided, it is broadcast to all F.
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
%   .transfer_functions
%   .posterior_source_covariances
%   .residual_transfer_functions
%   .residual_covariances
%   .initial_precision_matrices
%   .computation_stats
%   .success

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
L                = input_data.leadfield_matrix;
Sigma_emp        = input_data.empirical_covariances;       % cell{F} or p×p
Sigma_prior_in   = input_data.source_prior_covariances;    % cell{F} or n×n
Sigma_xi_xi      = input_data.noise_covariance;
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
estep_results.transfer_functions             = cell(F,1);
estep_results.posterior_source_covariances   = cell(F,1);
estep_results.residual_transfer_functions    = cell(F,1);
estep_results.residual_covariances           = cell(F,1);
estep_results.initial_precision_matrices     = cell(F,1);
estep_results.computation_stats              = struct();
estep_results.success                        = false;

stats = struct();
stats.processing_times         = zeros(F,5);
stats.condition_numbers        = zeros(F,3);
stats.regularization_applied   = false(F,1);
stats.eigenvalue_corrections   = zeros(F,2);
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
        S_vv_f    = Sigma_emp{f};
        Sigma_jj_f= Sigma_prior_in{f};
        if ~isnumeric(S_vv_f) || ~isequal(size(S_vv_f), [p p])
            error('Empirical covariance at f=%d must be %dx%d numeric', f, p, p);
        end
        if ~isnumeric(Sigma_jj_f) || ~isequal(size(Sigma_jj_f), [n n])
            error('Source prior at f=%d must be %dx%d numeric', f, n, n);
        end

        % 1) DSTF
        t1 = tic;
        opt_dstf = struct('regularization_factor', estep_params.regularization_factor, ...
                          'condition_threshold',   estep_params.condition_threshold, ...
                          'verbose', false);
        T_jv_f = module2_dstf_computation(L, Sigma_jj_f, Sigma_xi_xi, opt_dstf);
        estep_results.transfer_functions{f} = T_jv_f;
        stats.processing_times(f,1) = toc(t1);

        % 2) Posterior Σ_jj,post
        t2 = tic;
        opt_psc = struct('regularization_factor', estep_params.regularization_factor, ...
                         'condition_threshold',   estep_params.condition_threshold, ...
                         'min_eigenvalue_ratio',  estep_params.min_eigenvalue_ratio, ...
                         'verbose', false);
        Sigma_post_f = module2_posterior_source_covariance(Sigma_jj_f, L, Sigma_xi_xi, opt_psc);
        estep_results.posterior_source_covariances{f} = Sigma_post_f;
        stats.processing_times(f,2) = toc(t2);

        % 3) Residual TF
        t3 = tic;
        T_xi_v_f = module2_residual_transfer_function(T_jv_f, L, ...
                        'validate_properties', true, ...
                        'enforce_hermitian', false, ...   % keep the complementary identity
                        'verbose', false);
        estep_results.residual_transfer_functions{f} = T_xi_v_f;
        stats.processing_times(f,3) = toc(t3);

        % 4) Residual empirical covariance S_ξξ
        t4 = tic;
        min_thr = estep_params.min_eigenvalue_ratio * max(real(eig(S_vv_f)));
        if ~isfinite(min_thr) || min_thr < 0, min_thr = 0; end
        opt_rec = struct('regularization_factor',  estep_params.regularization_factor, ...
                         'min_eigenvalue_threshold', max(min_thr, 0), ...
                         'enforce_hermitian', true, ...
                         'verbose', false);
        S_xi_xi_f = module2_residual_empirical_covariance(T_xi_v_f, S_vv_f, opt_rec);
        estep_results.residual_covariances{f} = S_xi_xi_f;
        stats.processing_times(f,4) = toc(t4);

        % 5) Initial precision Ω (eigenvalue floor + ridge-like smoothing)
        t5 = tic;
        [V, D]   = eig(S_xi_xi_f);
        evals    = real(diag(D));
        eps_reg  = estep_params.regularization_factor;
        dmax     = max(evals);
        evals_reg= (evals + eps_reg * dmax) / (1 + eps_reg);

        neg_cnt  = sum(evals < 0);
        small_cnt= sum(evals < estep_params.min_eigenvalue_ratio * max(dmax,eps) & evals >= 0);
        stats.eigenvalue_corrections(f,:) = [neg_cnt, small_cnt];

        Omega_f  = V * diag(1 ./ evals_reg) * V';
        Omega_f  = (Omega_f + Omega_f') / 2;   % keep Hermitian
        estep_results.initial_precision_matrices{f} = Omega_f;
        stats.processing_times(f,5) = toc(t5);

        % Diagnostics
        try
            stats.condition_numbers(f,1) = cond(L * Sigma_jj_f * L' + Sigma_xi_xi);
        catch, stats.condition_numbers(f,1) = NaN; end
        try
            stats.condition_numbers(f,2) = cond(Sigma_post_f);
        catch, stats.condition_numbers(f,2) = NaN; end
        try
            stats.condition_numbers(f,3) = cond(S_xi_xi_f);
        catch, stats.condition_numbers(f,3) = NaN; end

        stats.regularization_applied(f) = (eps_reg > 0) || (neg_cnt + small_cnt > 0);
        stats.successful_frequencies = stats.successful_frequencies + 1;

        if estep_params.verbose
            fprintf('  ✓ f=%d/%d (%.2f Hz): done in %.3fs\n', ...
                f, F, frequencies(min(f,numel(frequencies))), toc(freq_tic));
        end

    catch ME
        stats.failed_frequencies = stats.failed_frequencies + 1;
        stats.error_messages{end+1} = sprintf('f=%d: %s', f, ME.message);

        estep_results.transfer_functions{f}           = [];
        estep_results.posterior_source_covariances{f} = [];
        estep_results.residual_transfer_functions{f}  = [];
        estep_results.residual_covariances{f}         = [];
        estep_results.initial_precision_matrices{f}   = [];

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

% ---- Final validation (soft) --------------------------------------------
if estep_results.success
    good = 0;
    for f = 1:F
        Om = estep_results.initial_precision_matrices{f};
        if ~isempty(Om)
            try
                lam_min = min(real(eig(Om)));
                if ishermitian(Om) && lam_min > -1e-10
                    good = good + 1;
                end
            catch
                % ignore
            end
        end
    end
    if good < stats.successful_frequencies
        warning('module2_estep_main:invalid_precision_matrices', ...
            'Some precision matrices may not be positive definite.');
    end
end
end
