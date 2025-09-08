function estep_results = module2_estep_main(input_data, estep_params)
% MODULE2_ESTEP_MAIN - E-step (Full Source Model)
% Outputs:
%   .transfer_functions                    {F} (n×p)   T_{jv}
%   .posterior_source_covariances          {F} (n×n)   Σ_{jj,post}
%   .source_second_moments                 {F} (n×n)   Ŝ_{jj} = T S_{vv} T^H + Σ_{jj,post}
%   .residual_transfer_functions           {F} (p×p)   T_{ξv} = I - L T_{jv}
%   .residual_covariances                  {F} (p×p)   S_{ξξ} = T_{ξv} S_{vv} T_{ξv}^H
%   .initial_precision_matrices            {F} (n×n)   Ω_src,init = inv_robust(Σ_{jj,post})
%   .residual_initial_precision_matrices   {F} (p×p)   Ω_resid,init = inv_robust(S_{ξξ})  [diagnostic]

% -------- Input presence --------
if nargin < 1, error('module2_estep_main:insufficient_input','input_data is required'); end
if nargin < 2, estep_params = struct(); end
req = {'leadfield_matrix','empirical_covariances','source_prior_covariances','noise_covariance','frequencies'};
for k = 1:numel(req)
    if ~isfield(input_data, req{k})
        error('module2_estep_main:missing_field','Missing input_data.%s', req{k});
    end
end

% -------- Extract & robustify --------
L              = input_data.leadfield_matrix;        % p×n (typically real)
Sigma_emp      = input_data.empirical_covariances;   % {F} or p×p (can be complex)
Sigma_prior_in = input_data.source_prior_covariances;% {F} or n×n (typically real)
Sigma_xi_xi    = input_data.noise_covariance;        % p×p (typically real)
frequencies    = input_data.frequencies;

if ~iscell(Sigma_emp),      Sigma_emp      = {Sigma_emp};      end
if ~iscell(Sigma_prior_in), Sigma_prior_in = {Sigma_prior_in}; end

if ~isnumeric(L) || ndims(L) ~= 2, error('module2_estep_main:invalid_leadfield','L must be 2D numeric'); end
[p, n] = size(L);
if ~isnumeric(Sigma_xi_xi) || ~isequal(size(Sigma_xi_xi), [p p])
    error('module2_estep_main:invalid_noise_covariance','noise_covariance must be %dx%d', p, p);
end
F = numel(Sigma_emp);
if ~isvector(frequencies) || numel(frequencies) ~= F, frequencies = 1:F; end

if numel(Sigma_prior_in) == 1 && F > 1
    Sigma_prior_in = repmat(Sigma_prior_in, F, 1);
end
if numel(Sigma_prior_in) ~= F
    error('module2_estep_main:source_covariance_mismatch', ...
          'Number of source priors (%d) must match number of empirical covariances (%d).', ...
          numel(Sigma_prior_in), F);
end

% -------- Params with defaults --------
def = struct('regularization_factor',1e-8, ...
             'condition_threshold', 1e12, ...
             'min_eigenvalue_ratio',1e-12, ...
             'verbose', false);
fn = fieldnames(def);
for i = 1:numel(fn)
    if ~isfield(estep_params, fn{i}), estep_params.(fn{i}) = def.(fn{i}); end
end

% -------- Prepare outputs --------
estep_results = struct();
estep_results.transfer_functions                    = cell(F,1); % n×p
estep_results.posterior_source_covariances          = cell(F,1); % n×n
estep_results.source_second_moments                 = cell(F,1); % n×n
estep_results.residual_transfer_functions           = cell(F,1); % p×p
estep_results.residual_covariances                  = cell(F,1); % p×p
estep_results.initial_precision_matrices            = cell(F,1); % n×n
estep_results.residual_initial_precision_matrices   = cell(F,1); % p×p
estep_results.computation_stats                     = struct();
estep_results.success                               = false;

stats = struct();
stats.processing_times         = zeros(F,6);
stats.condition_numbers        = zeros(F,3);
stats.regularization_applied   = false(F,1);
stats.eigenvalue_corrections   = zeros(F,4); % [neg/small Σpost, neg/small S_xi_xi]
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

% -------- Main loop --------
for f = 1:F
    try
        tic_f = tic;

        % Validate per-frequency sizes
        S_vv_f     = Sigma_emp{f};      % p×p (can be complex, Hermitian)
        Sigma_jj_f = Sigma_prior_in{f}; % n×n (typically real)
        if ~isnumeric(S_vv_f)     || ~isequal(size(S_vv_f), [p p]), error('S_vv(f) must be %dx%d', p,p); end
        if ~isnumeric(Sigma_jj_f) || ~isequal(size(Sigma_jj_f), [n n]), error('Sigma_jj(f) must be %dx%d', n,n); end

        % (2) Posterior source covariance (information-form)
        t2 = tic;
        opt_psc = struct('regularization_factor', estep_params.regularization_factor, ...
                         'jitter_max_tries', 6, ...
                         'ensure_positive_definite', true, ...
                         'min_eigenvalue_ratio',  estep_params.min_eigenvalue_ratio, ...
                         'verbose', false);
        Sigma_post_f = module2_posterior_source_covariance(Sigma_jj_f, L, Sigma_xi_xi, opt_psc); % n×n
        estep_results.posterior_source_covariances{f} = Sigma_post_f;
        stats.processing_times(f,2) = toc(t2);

        % (1) T_{jv} from Σ_post:  T = Σ_post * L' * Σ_xi^{-1}
        t1 = tic;
        T_jv_f = Sigma_post_f * (L' / Sigma_xi_xi);             % n×p
        % Drop tiny imaginary leakage ONLY when all inputs for T are real
        if isreal(L) && isreal(Sigma_post_f) && isreal(Sigma_xi_xi)
            if max(abs(imag(T_jv_f(:)))) < 1e-13, T_jv_f = real(T_jv_f); end
        end
        estep_results.transfer_functions{f} = T_jv_f;
        stats.processing_times(f,1) = toc(t1);

        % (2.5) Source second moments: Ŝ_{jj} = T S_{vv} T^H + Σ_post
        t25 = tic;
        Sjj_hat_f = T_jv_f * S_vv_f * T_jv_f';                   % n×n (complex if S_vv is complex)
        Sjj_hat_f = (Sjj_hat_f + Sjj_hat_f')/2;                  % Hermitianize (does NOT drop complex off-diagonals)
        Sjj_hat_f = Sjj_hat_f + Sigma_post_f;
        estep_results.source_second_moments{f} = Sjj_hat_f;
        stats.processing_times(f,6) = toc(t25);

        % (3) Residual TF: T_{ξv} = I - L T_{jv}
        t3 = tic;
        T_xi_v_f = module2_residual_transfer_function(T_jv_f, L, ...
                        'validate_properties', true, ...
                        'numerical_tolerance', 1e-12, ...
                        'verbose', false, ...
                        'enforce_hermitian', false);
        estep_results.residual_transfer_functions{f} = T_xi_v_f; % p×p (typically real if L,T_jv real)
        stats.processing_times(f,3) = toc(t3);

        % (4) Residual empirical covariance: S_{ξξ} = T_{ξv} S_{vv} T_{ξv}^H
        t4 = tic;
        opt_rec = struct('enforce_hermitian', true, ...
                         'regularization_factor', 0, ...
                         'min_eigenvalue_threshold', 1e-12, ...
                         'psd_only', true, ...
                         'numerical_tolerance', 1e-12, ...
                         'verbose', false);
        S_xi_xi_f = module2_residual_empirical_covariance(T_xi_v_f, S_vv_f, opt_rec); % p×p (complex if S_vv complex)
        estep_results.residual_covariances{f} = S_xi_xi_f;
        stats.processing_times(f,4) = toc(t4);

        % (5) Initial precisions
        eps_reg = estep_params.regularization_factor;
        min_rat = estep_params.min_eigenvalue_ratio;
        [Omega_src_init_f, neg1, small1]   = inv_psd_robust(Sigma_post_f, eps_reg, min_rat); % n×n
        [Omega_resid_init_f, neg2, small2] = inv_psd_robust(S_xi_xi_f,  eps_reg, min_rat);   % p×p
        estep_results.initial_precision_matrices{f}          = Omega_src_init_f;   % source-domain
        estep_results.residual_initial_precision_matrices{f} = Omega_resid_init_f; % sensor-domain (diagnostic)

        stats.eigenvalue_corrections(f,:) = [neg1, small1, neg2, small2];

        % Diagnostics
        try, stats.condition_numbers(f,1) = cond(L * Sigma_jj_f * L' + Sigma_xi_xi); catch, stats.condition_numbers(f,1) = NaN; end
        try, stats.condition_numbers(f,2) = cond(Sigma_post_f);                       catch, stats.condition_numbers(f,2) = NaN; end
        try, stats.condition_numbers(f,3) = cond(S_xi_xi_f);                          catch, stats.condition_numbers(f,3) = NaN; end

        stats.regularization_applied(f) = true;
        stats.successful_frequencies = stats.successful_frequencies + 1;

        if estep_params.verbose
            fprintf('  ✓ f=%d/%d (%.2f Hz): done in %.3fs\n', ...
                f, F, frequencies(min(f,numel(frequencies))), toc(tic_f));
        end

    catch ME
        stats.failed_frequencies = stats.failed_frequencies + 1;
        stats.error_messages{end+1} = sprintf('f=%d: %s', f, ME.message);

        estep_results.transfer_functions{f}                    = [];
        estep_results.posterior_source_covariances{f}          = [];
        estep_results.source_second_moments{f}                 = [];
        estep_results.residual_transfer_functions{f}           = [];
        estep_results.residual_covariances{f}                  = [];
        estep_results.initial_precision_matrices{f}            = [];
        estep_results.residual_initial_precision_matrices{f}   = [];

        if estep_params.verbose
            fprintf('  ✗ f=%d failed: %s\n', f, ME.message);
        end
    end
end

stats.total_computation_time = toc(overall_tic);
estep_results.computation_stats = stats;
estep_results.success = (stats.successful_frequencies / max(F,1)) >= 0.8;

if estep_params.verbose
    fprintf('\n========================================\n');
    fprintf('E-step summary: %d/%d succeeded (%.1f%%), total %.3fs\n', ...
        stats.successful_frequencies, F, 100*stats.successful_frequencies/max(F,1), stats.total_computation_time);
    if ~isempty(stats.error_messages)
        fprintf('Errors:\n'); fprintf('  - %s\n', stats.error_messages{:});
    end
end

% ---- Soft validation for source-domain initial precisions ----
if estep_results.success
    good = 0;
    for f = 1:F
        Om = estep_results.initial_precision_matrices{f}; % n×n
        if ~isempty(Om)
            try
                lam_min = min(real(eig((Om+Om')/2)));
                if lam_min > -1e-10, good = good + 1; end
            catch, end
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
