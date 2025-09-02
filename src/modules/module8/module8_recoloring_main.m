function recoloring_results = module8_recoloring_main(input_data, recoloring_params)
% MODULE8_RECOLORING_MAIN - Inverse whitening transformation (recoloring) to original scale
%
% Mathematical transformation (correct):
%   Omega_ω = D(ω) * Gamma_tilde_ω * D(ω)
%
% Notes on correctness:
% - Whitening is defined as:  Sigma_tilde_ω = D(ω) * Sigma_ω * D(ω).
% - Therefore Sigma_ω = D(ω)^(-1) * Sigma_tilde_ω * D(ω)^(-1),
%   and Omega_ω = Sigma_ω^(-1) = D(ω) * Sigma_tilde_ω^(-1) * D(ω)
%                = D(ω) * Gamma_tilde_ω * D(ω).
%
% This file includes:
% - Vectorized recoloring (no inv(), element-wise scaling by d_i d_j)
% - Optional Hermitian enforcement
% - SPD validation (chol first, fallback eig) once per frequency
% - Quality metrics without redundant SPD checks
% - Robust sparsity mapping: thresholds scale by (d_i d_j), not division

% ==================== Input Validation ====================
if nargin < 1
    error('module8_recoloring_main:insufficient_input', ...
          'At least input_data is required');
end

if nargin < 2
    recoloring_params = struct();
end

if ~isstruct(input_data)
    error('module8_recoloring_main:invalid_input', 'input_data must be a structure');
end

required_fields = {'whitened_precision_matrices', 'whitening_matrices'};
for i = 1:numel(required_fields)
    if ~isfield(input_data, required_fields{i}) || isempty(input_data.(required_fields{i}))
        error('module8_recoloring_main:missing_field', ...
              'Required field %s is missing or empty', required_fields{i});
    end
end

Gamma_tilde = input_data.whitened_precision_matrices;
D = input_data.whitening_matrices;

if ~iscell(Gamma_tilde) || ~iscell(D)
    error('module8_recoloring_main:invalid_format', ...
          'Precision matrices and whitening matrices must be cell arrays');
end

F = numel(Gamma_tilde);
if numel(D) ~= F
    error('module8_recoloring_main:dimension_mismatch', ...
          'Number of precision matrices (%d) must match whitening matrices (%d)', ...
          numel(Gamma_tilde), numel(D));
end

if F == 0
    error('module8_recoloring_main:empty_input', 'Input data is empty');
end

p = size(Gamma_tilde{1}, 1);
for omega = 1:F
    % Precision matrix checks
    if ~isnumeric(Gamma_tilde{omega}) || ~ismatrix(Gamma_tilde{omega})
        error('module8_recoloring_main:invalid_precision', ...
              'Precision matrix at frequency %d must be numeric matrix', omega);
    end
    if ~isequal(size(Gamma_tilde{omega}), [p p])
        error('module8_recoloring_main:precision_size_mismatch', ...
              'All precision matrices must be %dx%d, got %dx%d at ω=%d', ...
              p, p, size(Gamma_tilde{omega}, 1), size(Gamma_tilde{omega}, 2), omega);
    end
    % Whitening matrix checks
    if ~isnumeric(D{omega}) || ~ismatrix(D{omega})
        error('module8_recoloring_main:invalid_whitening', ...
              'Whitening matrix at frequency %d must be numeric matrix', omega);
    end
    if ~isequal(size(D{omega}), [p p])
        error('module8_recoloring_main:whitening_size_mismatch', ...
              'All whitening matrices must be %dx%d, got %dx%d at ω=%d', ...
              p, p, size(D{omega}, 1), size(D{omega}, 2), omega);
    end
    % Diagonal warning (should be diagonal)
    off_diag_norm = norm(D{omega} - diag(diag(D{omega})), 'fro');
    if off_diag_norm > 1e-10
        warning('module8_recoloring_main:non_diagonal_whitening', ...
                'Whitening matrix at ω=%d appears non-diagonal (off-diag norm: %.2e)', ...
                omega, off_diag_norm);
    end
end

% ==================== Parameter Setup ====================
defaults = struct();
defaults.force_hermitian = true;
defaults.validate_spd = true;
defaults.g_min_threshold = 1e-12;  % floor for stability (on g_i -> on D diag)
defaults.inv_error_tolerance = []; % adaptive if empty
defaults.hermitian_tolerance = 1e-12;
defaults.compute_quality_metrics = true;
defaults.verbose = false;

fn = fieldnames(defaults);
for i = 1:numel(fn)
    f = fn{i};
    if ~isfield(recoloring_params, f)
        recoloring_params.(f) = defaults.(f);
    end
end

% ==================== Initialize Results Structure ====================
recoloring_results = struct();
recoloring_results.recolored_precision_matrices = cell(F, 1);
recoloring_results.recoloring_quality = cell(F, 1);
recoloring_results.transformation_stats = struct();
recoloring_results.computation_stats = struct();
recoloring_results.success = false;

comp_stats = struct();
comp_stats.processing_times = zeros(F, 1);
comp_stats.total_computation_time = 0;
comp_stats.successful_frequencies = 0;
comp_stats.failed_frequencies = 0;
comp_stats.error_messages = {};
comp_stats.hermitian_enforcements = 0;
comp_stats.spd_violations = 0;

overall_tic = tic;

if recoloring_params.verbose
    fprintf('========================================\n');
    fprintf('Module 8: Recoloring (Inverse Whitening)\n');
    fprintf('Processing %d frequencies | nodes=%d\n', F, p);
    fprintf('========================================\n\n');
end

% ==================== Main Recoloring Loop ====================
Omega = cell(F, 1);
quality_metrics = cell(F, 1);

for omega = 1:F
    freq_tic = tic;
    try
        if recoloring_params.verbose
            fprintf('Processing frequency %d/%d... ', omega, F);
        end

        Gamma_tilde_omega = Gamma_tilde{omega};
        D_omega = D{omega};

        validate_frequency_matrices(Gamma_tilde_omega, D_omega, omega, recoloring_params);

        % Core recoloring (vectorized) with clamped diagonal for stability
        Omega_omega = apply_recoloring_transformation(Gamma_tilde_omega, D_omega, recoloring_params);

        % Optional Hermitian enforcement (uses conjugate transpose ')
        if recoloring_params.force_hermitian
            Omega_omega = force_hermitian_symmetry(Omega_omega);
            comp_stats.hermitian_enforcements = comp_stats.hermitian_enforcements + 1;
        end

        % SPD validation once (reused by quality metrics)
        spd_status = struct('is_spd', true, 'details', []);
        if recoloring_params.validate_spd
            [is_spd, spd_info] = validate_spd_property(Omega_omega, omega);
            spd_status.is_spd = is_spd;
            spd_status.details = spd_info;
            if ~is_spd
                comp_stats.spd_violations = comp_stats.spd_violations + 1;
                if recoloring_params.verbose
                    fprintf('SPD violation at ω=%d: %s\n', omega, spd_info.message);
                end
            end
        end

        Omega{omega} = Omega_omega;

        if recoloring_params.compute_quality_metrics
            quality_metrics{omega} = compute_quality_metrics( ...
                Omega_omega, Gamma_tilde_omega, D_omega, omega, input_data, recoloring_params, spd_status);
        end

        comp_stats.processing_times(omega) = toc(freq_tic);
        comp_stats.successful_frequencies = comp_stats.successful_frequencies + 1;

        if recoloring_params.verbose
            fprintf('completed in %.3fs\n', comp_stats.processing_times(omega));
        end

    catch ME
        comp_stats.failed_frequencies = comp_stats.failed_frequencies + 1;
        comp_stats.error_messages{end+1} = sprintf('Frequency %d: %s', omega, ME.message);
        if recoloring_params.verbose
            fprintf('failed: %s\n', ME.message);
        end
        Omega{omega} = [];
        quality_metrics{omega} = struct('success', false, 'error', ME.message);
    end
end

% ==================== Finalize Results ====================
comp_stats.total_computation_time = toc(overall_tic);
comp_stats.average_time_per_frequency = comp_stats.total_computation_time / F;

recoloring_results.recolored_precision_matrices = Omega;
recoloring_results.recoloring_quality = quality_metrics;
recoloring_results.computation_stats = comp_stats;

if recoloring_params.compute_quality_metrics
    recoloring_results.transformation_stats = compute_transformation_stats(quality_metrics, F);
end

recoloring_results.success = (comp_stats.failed_frequencies == 0);

if recoloring_params.verbose
    fprintf('\n========================================\n');
    fprintf('Module 8 Recoloring Complete\n');
    fprintf('========================================\n');
    fprintf('Total time: %.3f seconds\n', comp_stats.total_computation_time);
    fprintf('Successful frequencies: %d/%d\n', comp_stats.successful_frequencies, F);
    if comp_stats.failed_frequencies > 0
        fprintf('Failed frequencies: %d\n', comp_stats.failed_frequencies);
    end
    if comp_stats.hermitian_enforcements > 0
        fprintf('Hermitian enforcements: %d\n', comp_stats.hermitian_enforcements);
    end
    if comp_stats.spd_violations > 0
        fprintf('SPD violations detected: %d\n', comp_stats.spd_violations);
    end
    fprintf('\n');
end

end

%% ==================== Helper Functions ====================

function validate_frequency_matrices(Gamma_tilde_omega, D_omega, omega, params)
% Validate per-frequency inputs and basic numerical health

if any(~isfinite(Gamma_tilde_omega(:)))
    error('module8_recoloring_main:invalid_precision_values', ...
          'Precision matrix at ω=%d contains NaN or Inf values', omega);
end

if any(~isfinite(D_omega(:)))
    error('module8_recoloring_main:invalid_whitening_values', ...
          'Whitening matrix at ω=%d contains NaN or Inf values', omega);
end

% Hermitian check (use conjugate transpose ')
herm_err = norm(Gamma_tilde_omega - Gamma_tilde_omega', 'fro');
rel_herm_err = herm_err / max(norm(Gamma_tilde_omega, 'fro'), eps);
if rel_herm_err > 1e-10
    warning('module8_recoloring_main:non_hermitian_input', ...
            'Input precision matrix at ω=%d not Hermitian (error: %.2e)', omega, rel_herm_err);
end

% D diagonal positivity
d_diag = diag(D_omega);
if any(real(d_diag) <= 0)
    error('module8_recoloring_main:invalid_whitening_diagonal', ...
          'Whitening matrix at ω=%d has non-positive diagonal elements', omega);
end

% Warn on extremely small diagonals (potential instability)
min_diag = min(real(d_diag));
if min_diag < max(sqrt(params.g_min_threshold), realmin)
    warning('module8_recoloring_main:small_whitening_diagonal', ...
            'Very small whitening diagonal at ω=%d (min: %.2e < %.2e)', ...
            omega, min_diag, sqrt(params.g_min_threshold));
end
end

function Omega_omega = apply_recoloring_transformation(Gamma_tilde_omega, D_omega, params)
% Vectorized recoloring: Omega = D * Gamma_tilde * D  <=>  Omega_ij = d_i * Gamma_ij * d_j
% Clamp D diagonals for stability (floor on g_i implies floor on D entries)

d_vec = diag(D_omega);
d_vec = max(d_vec, sqrt(params.g_min_threshold));   % clamp tiny scales
Douter = d_vec * d_vec.';                           % element-wise scaling matrix
Omega_omega = Gamma_tilde_omega .* Douter;          % no inv(), no explicit multiplications
end

function Omega_symmetric = force_hermitian_symmetry(Omega_omega)
% Enforce Hermitian symmetry (use conjugate transpose)
Omega_symmetric = (Omega_omega + Omega_omega') / 2;
end

function [is_spd, spd_info] = validate_spd_property(Omega_omega, omega)
% SPD validation: try chol first, fallback to eig (with small negative tolerance)

spd_info = struct();
spd_info.method_used = '';
spd_info.message = '';
spd_info.condition_number = NaN;
spd_info.min_eigenvalue = NaN;

try
    chol(Omega_omega);
    is_spd = true;
    spd_info.method_used = 'chol';
    spd_info.condition_number = cond(Omega_omega);
    spd_info.message = 'SPD verified by Cholesky';
catch
    try
        e = eig((Omega_omega + Omega_omega')/2);
        re = real(e);
        spd_info.min_eigenvalue = min(re);
        is_spd = all(re > -1e-10); % small negative tolerance
        spd_info.method_used = 'eig';
        if is_spd
            spd_info.condition_number = max(re) / max(min(re), eps);
            spd_info.message = sprintf('SPD verified by eigenvalues (lambda_min: %.2e)', spd_info.min_eigenvalue);
        else
            spd_info.message = sprintf('Not SPD: lambda_min = %.2e', spd_info.min_eigenvalue);
        end
    catch ME2
        is_spd = false;
        spd_info.method_used = 'failed';
        spd_info.message = sprintf('SPD validation failed: %s', ME2.message);
    end
end
end

function quality = compute_quality_metrics(Omega_omega, Gamma_tilde_omega, D_omega, omega, input_data, params, spd_status_precomputed)
% Compute per-frequency quality metrics. Reuse SPD status if already computed.

quality = struct();
quality.frequency_index = omega;
quality.matrix_size = size(Omega_omega, 1);

try
    % Hermitian error
    hermitian_diff = Omega_omega - Omega_omega';
    quality.hermitian_error = norm(hermitian_diff, 'fro') / max(norm(Omega_omega, 'fro'), eps);

    % SPD status (reuse precomputed)
    if nargin >= 7 && ~isempty(spd_status_precomputed)
        quality.spd_status = spd_status_precomputed;
        if isfield(spd_status_precomputed, 'details') && isfield(spd_status_precomputed.details, 'condition_number')
            quality.condition_number = spd_status_precomputed.details.condition_number;
        else
            quality.condition_number = NaN;
        end
    else
        [is_spd, spd_info] = validate_spd_property(Omega_omega, omega);
        quality.spd_status = struct('is_spd', is_spd, 'details', spd_info);
        quality.condition_number = spd_info.condition_number;
    end

    % Inverse relationship error if Sigma provided
    quality.inv_error = NaN;
    if isfield(input_data, 'original_covariances') && ~isempty(input_data.original_covariances) ...
       && omega <= numel(input_data.original_covariances)
        Sigma_original = input_data.original_covariances{omega};
        if ~isempty(Sigma_original)
            prodOS = Omega_omega * Sigma_original;
            I = eye(size(prodOS));
            quality.inv_error = norm(prodOS - I, 'fro') / sqrt(size(I,1));
        end
    end

    % Sparsity consistency (threshold scales by d_i d_j)
    if isfield(input_data, 'active_set_masks') && ~isempty(input_data.active_set_masks) ...
       && omega <= numel(input_data.active_set_masks)
        quality.sparsity_metrics = compute_sparsity_consistency( ...
            Omega_omega, Gamma_tilde_omega, D_omega, input_data.active_set_masks{omega});
    end

    % Scaling diagnostics
    quality.scaling_metrics = compute_scaling_metrics(Omega_omega, Gamma_tilde_omega, D_omega);

    quality.success = true;

catch ME
    quality.success = false;
    quality.error_message = ME.message;
    quality.hermitian_error = NaN;
    quality.condition_number = NaN;
    quality.inv_error = NaN;
end
end

function sparsity_metrics = compute_sparsity_consistency(Omega_omega, Gamma_tilde_omega, D_omega, active_mask)
% Compare sparsity before/after recoloring.
% IMPORTANT: Threshold in original domain must be scaled by (d_i d_j).

sparsity_metrics = struct();
try
    d_vec = diag(D_omega);
    p = length(d_vec);

    base_threshold = 1e-6;

    % Whitened-domain nnz by base threshold
    whitened_nonzeros = (abs(Gamma_tilde_omega) > base_threshold);

    % Original-domain threshold must MULTIPLY by (d_i d_j)
    scaled_threshold_matrix = base_threshold * (d_vec * d_vec.');  % <-- fix: use multiplication
    recolored_nonzeros = (abs(Omega_omega) > scaled_threshold_matrix);

    total_elements = p * p;
    consistent_elements = sum((whitened_nonzeros(:) & recolored_nonzeros(:)) | ...
                              (~whitened_nonzeros(:) & ~recolored_nonzeros(:)));

    sparsity_metrics.consistency_rate = consistent_elements / total_elements;
    sparsity_metrics.whitened_nnz = nnz(whitened_nonzeros);
    sparsity_metrics.recolored_nnz = nnz(recolored_nonzeros);
    sparsity_metrics.nnz_ratio = sparsity_metrics.recolored_nnz / max(sparsity_metrics.whitened_nnz, 1);

    if ~isempty(active_mask)
        active_elements = nnz(active_mask);
        recolored_active = recolored_nonzeros & active_mask;
        sparsity_metrics.active_set_consistency = nnz(recolored_active) / max(active_elements, 1);
    end

    sparsity_metrics.success = true;

catch ME
    sparsity_metrics.success = false;
    sparsity_metrics.error = ME.message;
end
end

function scaling_metrics = compute_scaling_metrics(Omega_omega, Gamma_tilde_omega, D_omega)
% Analyze scaling effect of recoloring (should match d_i^2 on diagonal, etc.)

scaling_metrics = struct();
try
    d_vec = diag(D_omega);
    scaling_matrix = d_vec * d_vec.';
    scaling_metrics.min_scaling = min(scaling_matrix(:));
    scaling_metrics.max_scaling = max(scaling_matrix(:));
    scaling_metrics.mean_scaling = mean(scaling_matrix(:));
    scaling_metrics.scaling_ratio = scaling_metrics.max_scaling / max(scaling_metrics.min_scaling, eps);

    scaling_metrics.frobenius_ratio = norm(Omega_omega, 'fro') / max(norm(Gamma_tilde_omega, 'fro'), eps);
    scaling_metrics.spectral_ratio = norm(Omega_omega, 2) / max(norm(Gamma_tilde_omega, 2), eps);

    diag_omega = diag(Omega_omega);
    diag_gamma = diag(Gamma_tilde_omega);
    diag_scaling_actual = diag_omega ./ max(diag_gamma, eps);
    diag_scaling_expected = d_vec .^ 2;

    scaling_metrics.diagonal_scaling_error = norm(diag_scaling_actual - diag_scaling_expected) / ...
                                             max(norm(diag_scaling_expected), eps);
    scaling_metrics.success = true;
catch ME
    scaling_metrics.success = false;
    scaling_metrics.error = ME.message;
end
end

function transformation_stats = compute_transformation_stats(quality_metrics, F)
% Aggregate per-frequency metrics

transformation_stats = struct();

successful = cellfun(@(q) isstruct(q) && isfield(q,'success') && q.success, quality_metrics);
n_successful = sum(successful);
if n_successful == 0
    transformation_stats.success = false;
    transformation_stats.message = 'No successful frequency transformations';
    return;
end

inv_errors = NaN(F, 1);
hermitian_errors = NaN(F, 1);
condition_numbers = NaN(F, 1);
spd_successes = false(F, 1);

for omega = 1:F
    if successful(omega)
        q = quality_metrics{omega};
        if isfield(q,'inv_error') && isfinite(q.inv_error), inv_errors(omega) = q.inv_error; end
        if isfield(q,'hermitian_error') && isfinite(q.hermitian_error), hermitian_errors(omega) = q.hermitian_error; end
        if isfield(q,'condition_number') && isfinite(q.condition_number), condition_numbers(omega) = q.condition_number; end
        if isfield(q,'spd_status') && isfield(q.spd_status,'is_spd'), spd_successes(omega) = q.spd_status.is_spd; end
    end
end

transformation_stats.n_successful_frequencies = n_successful;
transformation_stats.success_rate = n_successful / F;

valid_inv = inv_errors(~isnan(inv_errors));
if ~isempty(valid_inv)
    s = struct();
    s.mean = mean(valid_inv); s.median = median(valid_inv);
    s.max = max(valid_inv);  s.std = std(valid_inv);
    s.percentile_95 = prctile(valid_inv, 95);
    transformation_stats.inv_error_stats = s;
end

valid_herm = hermitian_errors(~isnan(hermitian_errors));
if ~isempty(valid_herm)
    s = struct();
    s.mean = mean(valid_herm);
    s.max = max(valid_herm);
    s.all_below_tolerance = all(valid_herm < 1e-12);
    transformation_stats.hermitian_error_stats = s;
end

valid_cond = condition_numbers(~isnan(condition_numbers));
if ~isempty(valid_cond)
    s = struct();
    s.mean = mean(valid_cond);
    s.median = median(valid_cond);
    s.max = max(valid_cond);
    s.geometric_mean = exp(mean(log(valid_cond)));
    transformation_stats.condition_stats = s;
end

transformation_stats.spd_success_rate = sum(spd_successes) / F;
transformation_stats.overall_quality_score = compute_overall_quality_score(transformation_stats);
transformation_stats.success = transformation_stats.overall_quality_score > 0.8;

end

function quality_score = compute_overall_quality_score(stats)
% Combine multiple stats into a [0,1] quality score (simple weighted rule)

quality_score = 0; total_weight = 0;

w_success = 0.30; quality_score = quality_score + w_success * stats.success_rate; total_weight = total_weight + w_success;
w_spd = 0.25;     quality_score = quality_score + w_spd * stats.spd_success_rate; total_weight = total_weight + w_spd;

w_inv = 0.25;
if isfield(stats,'inv_error_stats')
    max_acc = 0.1;
    inv_q = max(0, 1 - stats.inv_error_stats.median / max_acc);
    quality_score = quality_score + w_inv * inv_q; total_weight = total_weight + w_inv;
end

w_herm = 0.20;
if isfield(stats,'hermitian_error_stats')
    herm_q = double(stats.hermitian_error_stats.all_below_tolerance);
    quality_score = quality_score + w_herm * herm_q; total_weight = total_weight + w_herm;
end

if total_weight > 0
    quality_score = quality_score / total_weight;
else
    quality_score = 0;
end
end
