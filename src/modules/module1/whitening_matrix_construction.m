function [D, whitening_stats] = whitening_matrix_construction(g_smooth, varargin)
% WHITENING_MATRIX_CONSTRUCTION - Construct whitening matrices from smoothed diagonal powers
%
% This function constructs the whitening matrix D(f) from smoothed diagonal powers
% according to: D(f) = diag(1/sqrt(g_1(f)), ..., 1/sqrt(g_n(f)))
%
% Usage:
%   [D, stats] = whitening_matrix_construction(g_smooth)
%   [D, stats] = whitening_matrix_construction(g_smooth, 'param', value, ...)
%
% Inputs:
%   g_smooth - Cell array of smoothed diagonal powers {F x 1}
%             Each element is real vector of size [n x 1]
%   Optional parameters (name-value pairs):
%     'min_power'      - Minimum power threshold to prevent amplification (default: 1e-10)
%     'max_condition'  - Maximum condition number allowed (default: 1e8)
%     'regularization' - Regularization method: 'floor', 'shrink', 'adaptive' (default: 'adaptive')
%     'shrink_factor'  - Shrinkage factor for 'shrink' method (default: 0.1)
%     'adaptive_percentile' - Percentile for adaptive regularization (default: 0.05)
%     'verbose'        - Display detailed information (default: true)
%
% Outputs:
%   D              - Cell array of whitening matrices {F x 1}
%                   Each element is real diagonal matrix of size [n x n]
%   whitening_stats - Structure containing statistics:
%     .condition_numbers - Condition numbers of whitening matrices [F x 1]
%     .power_ranges     - Min/max powers before regularization [F x 2]
%     .regularized_count - Number of regularized entries per frequency [F x 1]
%     .final_power_ranges - Min/max powers after regularization [F x 2]
%
% File location: src/modules/module1_preprocessing/whitening_matrix_construction.m

    % Parse input arguments
    p = inputParser;
    addRequired(p, 'g_smooth', @(x) iscell(x) && ~isempty(x));
    addParameter(p, 'min_power', 1e-10, @(x) isscalar(x) && x > 0);
    addParameter(p, 'max_condition', 1e8, @(x) isscalar(x) && x > 1);
    addParameter(p, 'regularization', 'adaptive', ...
                @(x) ismember(x, {'floor', 'shrink', 'adaptive'}));
    addParameter(p, 'shrink_factor', 0.1, @(x) isscalar(x) && x > 0 && x < 1);
    addParameter(p, 'adaptive_percentile', 0.05, @(x) isscalar(x) && x > 0 && x < 0.5);
    addParameter(p, 'verbose', true, @islogical);
    
    parse(p, g_smooth, varargin{:});
    params = p.Results;
    
    % Validate input
    F = length(g_smooth);
    if F == 0
        error('whitening_matrix_construction:empty_input', 'Input power array is empty');
    end
    
    n = length(g_smooth{1});
    if n == 0
        error('whitening_matrix_construction:empty_powers', 'Power vectors are empty');
    end
    
    if params.verbose
        fprintf('Constructing whitening matrices for %d frequencies, %d nodes\n', F, n);
    end
    
    % Initialize outputs
    D = cell(F, 1);
    whitening_stats = initialize_stats_structure(F);
    
    % Compute global statistics for adaptive regularization
    global_stats = compute_global_power_statistics(g_smooth, params);
    
    % Process each frequency
    for omega = 1:F
        g_vals = g_smooth{omega};
        
        % Validate input powers
        validate_power_vector(g_vals, omega);
        
        % Store original power statistics
        whitening_stats.power_ranges(omega, :) = [min(g_vals), max(g_vals)];
        
        % Apply regularization if needed
        [g_regularized, n_regularized] = apply_power_regularization_enhanced(g_vals, params, global_stats);
        whitening_stats.regularized_count(omega) = n_regularized;
        whitening_stats.final_power_ranges(omega, :) = [min(g_regularized), max(g_regularized)];
        
        % Construct whitening matrix
        D_omega = construct_diagonal_matrix_safe(g_regularized, params);
        
        % Validate whitening matrix
        validate_whitening_matrix(D_omega, omega, params);
        
        % Store results
        D{omega} = D_omega;
        whitening_stats.condition_numbers(omega) = cond(D_omega);
        
        if params.verbose && mod(omega, max(1, floor(F/10))) == 0
            fprintf('Processed frequency %d/%d\n', omega, F);
        end
    end
    
    % Report summary statistics
    if params.verbose
        report_whitening_statistics(whitening_stats, params);
    end
    
    % Final validation
    validate_whitening_output(D, whitening_stats);
end

function stats = initialize_stats_structure(F)
% Initialize statistics structure
    stats = struct();
    stats.condition_numbers = zeros(F, 1);
    stats.power_ranges = zeros(F, 2);
    stats.regularized_count = zeros(F, 1);
    stats.final_power_ranges = zeros(F, 2);
end

function global_stats = compute_global_power_statistics(g_smooth, params)
% Compute global statistics across all frequencies for adaptive regularization
    
    F = length(g_smooth);
    all_powers = [];
    
    % Collect all power values
    for omega = 1:F
        all_powers = [all_powers; g_smooth{omega}(:)];
    end
    
    % Remove any problematic values before computing statistics
    valid_powers = all_powers(all_powers > 0 & isfinite(all_powers));
    
    global_stats = struct();
    if ~isempty(valid_powers)
        global_stats.median = median(valid_powers);
        global_stats.mean = mean(valid_powers);
        global_stats.std = std(valid_powers);
        global_stats.min_percentile = prctile(valid_powers, params.adaptive_percentile * 100);
        global_stats.max_percentile = prctile(valid_powers, (1 - params.adaptive_percentile) * 100);
        global_stats.robust_min = max(params.min_power, global_stats.min_percentile / 10);
        global_stats.robust_max = min(global_stats.max_percentile * 10, global_stats.median * 1000);
    else
        % Fallback values
        global_stats.median = 1.0;
        global_stats.mean = 1.0;
        global_stats.std = 0.1;
        global_stats.min_percentile = params.min_power;
        global_stats.max_percentile = 10.0;
        global_stats.robust_min = params.min_power;
        global_stats.robust_max = 10.0;
    end
end

function validate_power_vector(g_vals, omega)
% Validate input power vector
    
    if ~isvector(g_vals) || ~isreal(g_vals)
        error('whitening_matrix_construction:invalid_powers', ...
              'Powers at frequency %d must be real vector', omega);
    end
    
    if any(g_vals <= 0)
        error('whitening_matrix_construction:negative_powers', ...
              'Non-positive powers found at frequency %d', omega);
    end
    
    if any(~isfinite(g_vals))
        error('whitening_matrix_construction:infinite_powers', ...
              'Non-finite powers found at frequency %d', omega);
    end
end

function [g_regularized, n_regularized] = apply_power_regularization_enhanced(g_vals, params, global_stats)
% Apply enhanced regularization to power values
    
    g_regularized = g_vals(:); % Ensure column vector
    n_regularized = 0;
    
    switch params.regularization
        case 'adaptive'
            [g_regularized, n_regularized] = apply_adaptive_regularization(g_regularized, params, global_stats);
            
        case 'floor'
            % Simple floor regularization
            small_indices = g_regularized < params.min_power;
            g_regularized(small_indices) = params.min_power;
            n_regularized = sum(small_indices);
            
        case 'shrink'
            % Shrink extreme values towards robust center
            [g_regularized, n_regularized] = apply_shrinkage_regularization(g_regularized, params, global_stats);
            
        otherwise
            error('whitening_matrix_construction:invalid_regularization', ...
                  'Unknown regularization method: %s', params.regularization);
    end
    
    % Final safety check and condition number control
    [g_regularized, additional_reg] = apply_condition_number_control(g_regularized, params);
    n_regularized = n_regularized + additional_reg;
    
    % Ensure all values are within safe bounds
    g_regularized = max(g_regularized, params.min_power);
    g_regularized = min(g_regularized, global_stats.robust_max);
end

function [g_reg, n_reg] = apply_adaptive_regularization(g_vals, params, global_stats)
% Apply adaptive regularization based on global statistics
    
    g_reg = g_vals;
    n_reg = 0;
    
    % Identify outliers using robust statistics
    robust_center = global_stats.median;
    robust_scale = global_stats.std;
    
    % Lower bound: use global minimum percentile but not too aggressive
    lower_bound = max(params.min_power, global_stats.robust_min);
    
    % Upper bound: prevent extreme amplification
    upper_bound = global_stats.robust_max;
    
    % Apply bounds
    too_small = g_reg < lower_bound;
    too_large = g_reg > upper_bound;
    
    if any(too_small)
        % Use soft transition instead of hard threshold
        small_vals = g_reg(too_small);
        % Exponential transition to lower bound
        transition_factor = exp(-5 * (lower_bound - small_vals) / lower_bound);
        g_reg(too_small) = lower_bound + (small_vals - lower_bound) .* transition_factor;
        g_reg(too_small) = max(g_reg(too_small), lower_bound);
        n_reg = n_reg + sum(too_small);
    end
    
    if any(too_large)
        % Soft capping for large values
        large_vals = g_reg(too_large);
        % Logarithmic compression
        excess = (large_vals - upper_bound) / upper_bound;
        compressed_excess = upper_bound * log(1 + excess) / max(1, max(excess));
        g_reg(too_large) = upper_bound + compressed_excess;
        n_reg = n_reg + sum(too_large);
    end
end

function [g_reg, n_reg] = apply_shrinkage_regularization(g_vals, params, global_stats)
% Apply shrinkage-based regularization
    
    g_reg = g_vals;
    robust_center = global_stats.median;
    
    % Identify values needing shrinkage
    deviation_threshold = 3.0; % Standard deviations from median
    robust_mad = median(abs(g_vals - robust_center));
    
    if robust_mad > 0
        z_scores = abs(g_vals - robust_center) / robust_mad;
        outlier_mask = z_scores > deviation_threshold;
    else
        outlier_mask = false(size(g_vals));
    end
    
    % Also shrink values that are too extreme in absolute terms
    too_small = g_vals < params.min_power * 10;
    too_large = g_vals > robust_center * 100;
    
    shrink_mask = outlier_mask | too_small | too_large;
    n_reg = sum(shrink_mask);
    
    if n_reg > 0
        % Apply shrinkage towards robust center
        shrink_strength = params.shrink_factor;
        g_reg(shrink_mask) = (1 - shrink_strength) * g_vals(shrink_mask) + ...
                            shrink_strength * robust_center;
        
        % Ensure minimum threshold
        g_reg = max(g_reg, params.min_power);
    end
end

function [g_reg, n_additional] = apply_condition_number_control(g_vals, params)
% Control condition number by limiting power ratios
    
    g_reg = g_vals;
    n_additional = 0;
    
    % Compute current condition number (squared, since we'll take sqrt later)
    current_ratio = max(g_vals) / min(g_vals);
    max_allowed_ratio = sqrt(params.max_condition);
    
    if current_ratio > max_allowed_ratio
        % Need to compress the range
        
        % Use geometric mean as reference point
        geo_mean = exp(mean(log(g_vals)));
        
        % Compress towards geometric mean
        log_vals = log(g_vals);
        log_center = log(geo_mean);
        log_range = max(log_vals) - min(log_vals);
        max_allowed_log_range = log(max_allowed_ratio);
        
        if log_range > max_allowed_log_range
            compression_factor = max_allowed_log_range / log_range;
            compressed_log_vals = log_center + (log_vals - log_center) * compression_factor;
            g_reg = exp(compressed_log_vals);
            n_additional = length(g_vals); % All values were adjusted
        end
    end
end

function D_omega = construct_diagonal_matrix_safe(g_vals, params)
% Construct diagonal whitening matrix with safety checks
    
    n = length(g_vals);
    
    % Compute 1/sqrt(g_i) for diagonal entries with protection
    diagonal_entries = zeros(n, 1);
    
    for i = 1:n
        if g_vals(i) > 0
            diagonal_entries(i) = 1 / sqrt(g_vals(i));
        else
            % This should not happen after regularization, but safety first
            diagonal_entries(i) = 1 / sqrt(params.min_power);
        end
        
        % Additional safety: prevent extreme values
        if diagonal_entries(i) > 1e6
            diagonal_entries(i) = 1e6;
        end
        if diagonal_entries(i) < 1e-6
            diagonal_entries(i) = 1e-6;
        end
    end
    
    % Construct diagonal matrix
    D_omega = diag(diagonal_entries);
    
    % Ensure real-valued result
    if ~isreal(D_omega)
        warning('whitening_matrix_construction:complex_matrix', ...
                'Whitening matrix has complex entries, taking real part');
        D_omega = real(D_omega);
    end
end

function validate_whitening_matrix(D_omega, omega, params)
% Validate constructed whitening matrix
    
    % Check if matrix is diagonal
    off_diagonal_norm = norm(D_omega - diag(diag(D_omega)), 'fro');
    if off_diagonal_norm > 1e-14
        error('whitening_matrix_construction:not_diagonal', ...
              'Whitening matrix at frequency %d is not diagonal', omega);
    end
    
    % Check for positive diagonal entries
    diagonal_entries = diag(D_omega);
    if any(diagonal_entries <= 0)
        error('whitening_matrix_construction:negative_diagonal', ...
              'Negative diagonal entries in whitening matrix at frequency %d', omega);
    end
    
    % Check condition number
    condition_num = max(diagonal_entries) / min(diagonal_entries);
    if condition_num > params.max_condition
        warning('whitening_matrix_construction:high_condition', ...
                'High condition number (%.2e) at frequency %d', condition_num, omega);
    end
    
    % Check for finite entries
    if any(~isfinite(diagonal_entries))
        error('whitening_matrix_construction:infinite_entries', ...
              'Non-finite entries in whitening matrix at frequency %d', omega);
    end
end

function report_whitening_statistics(stats, params)
% Report summary statistics
    
    F = length(stats.condition_numbers);
    
    fprintf('\nWhitening matrix construction statistics:\n');
    fprintf('========================================\n');
    
    % Condition numbers
    fprintf('Condition numbers:\n');
    fprintf('  Mean: %.2e, Median: %.2e\n', ...
            mean(stats.condition_numbers), median(stats.condition_numbers));
    fprintf('  Min:  %.2e, Max:    %.2e\n', ...
            min(stats.condition_numbers), max(stats.condition_numbers));
    
    % Power ranges before regularization
    fprintf('\nOriginal power ranges:\n');
    fprintf('  Min powers - Mean: %.2e, Min: %.2e\n', ...
            mean(stats.power_ranges(:, 1)), min(stats.power_ranges(:, 1)));
    fprintf('  Max powers - Mean: %.2e, Max: %.2e\n', ...
            mean(stats.power_ranges(:, 2)), max(stats.power_ranges(:, 2)));
    
    % Power ranges after regularization
    fprintf('\nRegularized power ranges:\n');
    fprintf('  Min powers - Mean: %.2e, Min: %.2e\n', ...
            mean(stats.final_power_ranges(:, 1)), min(stats.final_power_ranges(:, 1)));
    fprintf('  Max powers - Mean: %.2e, Max: %.2e\n', ...
            mean(stats.final_power_ranges(:, 2)), max(stats.final_power_ranges(:, 2)));
    
    % Regularization statistics
    total_regularized = sum(stats.regularized_count);
    total_entries = F * size(stats.power_ranges, 1); % This should be F * n
    if total_entries > 0
        regularization_rate = total_regularized / total_entries * 100;
    else
        regularization_rate = 0;
    end
    
    fprintf('\nRegularization statistics:\n');
    fprintf('  Total regularized entries: %d (%.1f%%)\n', ...
            total_regularized, regularization_rate);
    fprintf('  Frequencies with regularization: %d/%d\n', ...
            sum(stats.regularized_count > 0), F);
    
    if total_regularized > 0
        fprintf('  Average regularized per frequency: %.1f\n', ...
                mean(stats.regularized_count(stats.regularized_count > 0)));
    end
    
    % Warnings
    high_condition_count = sum(stats.condition_numbers > params.max_condition * 0.1);
    if high_condition_count > 0
        fprintf('\nWarnings:\n');
        fprintf('  %d frequencies have condition numbers > %.1e\n', ...
                high_condition_count, params.max_condition * 0.1);
    end
    
    fprintf('\n');
end

function validate_whitening_output(D, stats)
% Final validation of whitening matrices
    
    F = length(D);
    
    if F ~= length(stats.condition_numbers)
        error('whitening_matrix_construction:inconsistent_output', ...
              'Number of matrices (%d) does not match statistics (%d)', ...
              F, length(stats.condition_numbers));
    end
    
    % Check matrix dimensions consistency
    if F > 0
        n = size(D{1}, 1);
        for omega = 1:F
            if size(D{omega}, 1) ~= n || size(D{omega}, 2) ~= n
                error('whitening_matrix_construction:inconsistent_dimensions', ...
                      'Matrix %d has inconsistent dimensions', omega);
            end
            
            % Check if truly diagonal
            if norm(D{omega} - diag(diag(D{omega})), 'fro') > 1e-12
                error('whitening_matrix_construction:not_diagonal_final', ...
                      'Final matrix %d is not diagonal', omega);
            end
            
            % Check positive definiteness
            if any(diag(D{omega}) <= 0)
                error('whitening_matrix_construction:not_positive_final', ...
                      'Final matrix %d has non-positive diagonal entries', omega);
            end
        end
    end
    
    fprintf('Whitening matrix construction validation passed\n');
end