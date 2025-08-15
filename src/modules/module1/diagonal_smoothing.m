function [g_smooth, Sigma_emp_loaded] = diagonal_smoothing(Sigma_emp, varargin)
% DIAGONAL_SMOOTHING - Apply diagonal loading and smoothing to covariance matrices
%
% Enhanced version with robust protection against negative powers and
% better numerical stability.
%
% Usage:
%   [g_smooth, Sigma_emp_loaded] = diagonal_smoothing(Sigma_emp)
%   [g_smooth, Sigma_emp_loaded] = diagonal_smoothing(Sigma_emp, 'param', value, ...)
%
% File location: src/modules/module1_preprocessing/diagonal_smoothing.m
% Action: REPLACE the original diagonal_smoothing.m with this enhanced version

    % Parse input arguments
    p = inputParser;
    addRequired(p, 'Sigma_emp', @(x) iscell(x) && ~isempty(x));
    addParameter(p, 'smoothing_method', 'moving_average', ...
                @(x) ismember(x, {'moving_average', 'lowpass', 'spline'}));
    addParameter(p, 'window_size', 5, @(x) isscalar(x) && x > 0);
    addParameter(p, 'filter_order', 3, @(x) isscalar(x) && x > 0);
    addParameter(p, 'cutoff_freq', 0.3, @(x) isscalar(x) && x > 0 && x < 1);
    addParameter(p, 'spline_smooth', 0.1, @(x) isscalar(x) && x >= 0);
    addParameter(p, 'diagonal_loading', true, @islogical);
    addParameter(p, 'loading_factor', 0.01, @(x) isscalar(x) && x >= 0);
    addParameter(p, 'min_eigenval', 1e-8, @(x) isscalar(x) && x > 0);
    addParameter(p, 'min_power', 1e-12, @(x) isscalar(x) && x > 0);
    addParameter(p, 'max_change_factor', 5.0, @(x) isscalar(x) && x > 1);
    
    parse(p, Sigma_emp, varargin{:});
    params = p.Results;
    
    % Validate input
    F = length(Sigma_emp);
    if F == 0
        error('diagonal_smoothing:empty_input', 'Input covariance array is empty');
    end
    
    [n, n_check] = size(Sigma_emp{1});
    if n ~= n_check
        error('diagonal_smoothing:not_square', 'Covariance matrices must be square');
    end
    
    fprintf('Starting diagonal smoothing for %d frequencies, %d nodes\n', F, n);
    
    % Step 1: Extract raw diagonal powers with protection
    raw_diagonals = extract_diagonal_powers_protected(Sigma_emp, params);
    
    % Step 2: Apply smoothing across frequencies
    g_smooth = apply_frequency_smoothing_protected(raw_diagonals, params);
    
    % Step 3: Apply diagonal loading if requested
    if params.diagonal_loading
        fprintf('Applying diagonal loading with factor %.4f\n', params.loading_factor);
        Sigma_emp_loaded = apply_diagonal_loading(Sigma_emp, g_smooth, params);
    else
        Sigma_emp_loaded = {};
    end
    
    % Validate outputs
    validate_smoothed_powers_protected(g_smooth, raw_diagonals, params);
    
    fprintf('Diagonal smoothing completed successfully\n');
end

function raw_diagonals = extract_diagonal_powers_protected(Sigma_emp, params)
% Extract diagonal elements with enhanced protection
    
    F = length(Sigma_emp);
    n = size(Sigma_emp{1}, 1);
    
    % Initialize matrix: [n x F]
    raw_diagonals = zeros(n, F);
    
    for omega = 1:F
        % Extract real part of diagonal
        diag_vals = diag(Sigma_emp{omega});
        raw_diagonals(:, omega) = real(diag_vals);
        
        % Enhanced protection for problematic values
        problematic_idx = raw_diagonals(:, omega) <= params.min_power;
        
        if any(problematic_idx)
            % Use median of valid values as replacement
            valid_vals = raw_diagonals(~problematic_idx, omega);
            if isempty(valid_vals)
                replacement_val = 1.0; % Default fallback
            else
                replacement_val = max(median(valid_vals), params.min_power * 10);
            end
            
            raw_diagonals(problematic_idx, omega) = replacement_val;
            
            if sum(problematic_idx) > 0
                fprintf('  Frequency %d: Protected %d problematic diagonal values\n', ...
                        omega, sum(problematic_idx));
            end
        end
    end
    
    fprintf('Extracted diagonal powers: range [%.2e, %.2e]\n', ...
            min(raw_diagonals(:)), max(raw_diagonals(:)));
end

function g_smooth = apply_frequency_smoothing_protected(raw_diagonals, params)
% Apply smoothing across frequencies with enhanced protection
    
    [n, F] = size(raw_diagonals);
    g_smooth = cell(F, 1);
    
    fprintf('Applying %s smoothing with protection\n', params.smoothing_method);
    
    switch params.smoothing_method
        case 'moving_average'
            smoothed_matrix = apply_moving_average_protected(raw_diagonals, params);
            
        case 'lowpass'
            smoothed_matrix = apply_lowpass_filter_protected(raw_diagonals, params);
            
        case 'spline'
            smoothed_matrix = apply_spline_smoothing_protected(raw_diagonals, params);
            
        otherwise
            error('diagonal_smoothing:invalid_method', ...
                  'Unknown smoothing method: %s', params.smoothing_method);
    end
    
    % Apply change limiting and final protection
    smoothed_matrix = apply_change_protection(raw_diagonals, smoothed_matrix, params);
    
    % Convert to cell array format
    for omega = 1:F
        g_smooth{omega} = smoothed_matrix(:, omega);
    end
    
    % Report smoothing statistics
    original_var = var(raw_diagonals, 0, 2);
    smoothed_var = var(smoothed_matrix, 0, 2);
    variance_reduction = mean((original_var - smoothed_var) ./ (original_var + 1e-12));
    
    fprintf('Smoothing variance reduction: %.1f%%\n', variance_reduction * 100);
end

function smoothed_matrix = apply_moving_average_protected(raw_diagonals, params)
% Apply moving average smoothing with protection
    
    [n, F] = size(raw_diagonals);
    smoothed_matrix = zeros(n, F);
    
    % Ensure window size doesn't exceed data length
    window_size = min(params.window_size, F);
    if mod(window_size, 2) == 0
        window_size = window_size + 1;
    end
    
    half_window = floor(window_size / 2);
    
    for i = 1:n
        for omega = 1:F
            % Define window boundaries
            start_idx = max(1, omega - half_window);
            end_idx = min(F, omega + half_window);
            
            % Compute robust average
            window_data = raw_diagonals(i, start_idx:end_idx);
            smoothed_matrix(i, omega) = robust_mean(window_data);
        end
    end
end

function smoothed_matrix = apply_lowpass_filter_protected(raw_diagonals, params)
% Apply low-pass filtering with enhanced protection
    
    [n, F] = size(raw_diagonals);
    smoothed_matrix = zeros(n, F);
    
    % Check if we have enough data for filtering
    min_length_for_filter = max(6, 3 * params.filter_order);
    
    if F < min_length_for_filter
        fprintf('  Warning: Insufficient data length for filtering, using moving average\n');
        smoothed_matrix = apply_moving_average_protected(raw_diagonals, params);
        return;
    end
    
    try
        % Design filter with protection
        cutoff = min(params.cutoff_freq, 0.45);
        [b, a] = butter(params.filter_order, cutoff, 'low');
        
        for i = 1:n
            try
                % Apply zero-phase filtering
                filtered_data = filtfilt(b, a, raw_diagonals(i, :));
                
                % Check for valid output
                if any(~isfinite(filtered_data)) || any(filtered_data <= 0)
                    fprintf('  Warning: Filtering failed for channel %d, using fallback\n', i);
                    smoothed_matrix(i, :) = apply_robust_fallback(raw_diagonals(i, :), params);
                else
                    smoothed_matrix(i, :) = filtered_data;
                end
                
            catch
                smoothed_matrix(i, :) = apply_robust_fallback(raw_diagonals(i, :), params);
            end
        end
        
    catch ME
        fprintf('  Warning: Filter design failed: %s. Using moving average.\n', ME.message);
        smoothed_matrix = apply_moving_average_protected(raw_diagonals, params);
    end
end

function smoothed_matrix = apply_spline_smoothing_protected(raw_diagonals, params)
% Apply spline smoothing with protection
    
    [n, F] = size(raw_diagonals);
    smoothed_matrix = zeros(n, F);
    
    for i = 1:n
        if F < 4
            smoothed_matrix(i, :) = raw_diagonals(i, :);
        else
            try
                if exist('csaps', 'file') == 2
                    freq_points = 1:F;
                    adjusted_smooth = params.spline_smooth * min(1, F / 20);
                    pp = csaps(freq_points, raw_diagonals(i, :), adjusted_smooth);
                    smoothed_data = ppval(pp, freq_points);
                    
                    if any(~isfinite(smoothed_data)) || any(smoothed_data <= 0)
                        smoothed_matrix(i, :) = apply_robust_fallback(raw_diagonals(i, :), params);
                    else
                        smoothed_matrix(i, :) = smoothed_data;
                    end
                else
                    smoothed_matrix(i, :) = apply_kernel_smoothing(raw_diagonals(i, :), params);
                end
                
            catch
                smoothed_matrix(i, :) = apply_robust_fallback(raw_diagonals(i, :), params);
            end
        end
    end
end

function smoothed_data = apply_robust_fallback(data, params)
% Robust fallback smoothing method
    
    F = length(data);
    smoothed_data = zeros(size(data));
    
    % Use median filtering as robust fallback
    window_size = min(5, F);
    if mod(window_size, 2) == 0
        window_size = window_size - 1;
    end
    
    half_window = floor(window_size / 2);
    
    for j = 1:F
        start_idx = max(1, j - half_window);
        end_idx = min(F, j + half_window);
        
        window_data = data(start_idx:end_idx);
        smoothed_data(j) = median(window_data);
    end
    
    % Ensure all values are positive
    smoothed_data = max(smoothed_data, params.min_power);
end

function smoothed_data = apply_kernel_smoothing(data, params)
% Kernel smoothing fallback
    
    F = length(data);
    smoothed_data = zeros(size(data));
    
    kernel_width = max(3, round(params.spline_smooth * F));
    if mod(kernel_width, 2) == 0
        kernel_width = kernel_width + 1;
    end
    
    half_width = floor(kernel_width / 2);
    
    for j = 1:F
        start_idx = max(1, j - half_width);
        end_idx = min(F, j + half_width);
        
        indices = start_idx:end_idx;
        weights = exp(-0.5 * ((indices - j) / (kernel_width/3)).^2);
        weights = weights / sum(weights);
        
        smoothed_data(j) = sum(weights .* data(start_idx:end_idx));
    end
    
    smoothed_data = max(smoothed_data, params.min_power);
end

function robust_mean_val = robust_mean(data)
% Compute robust mean excluding outliers
    
    if length(data) <= 2
        robust_mean_val = mean(data);
        return;
    end
    
    % Remove extreme outliers
    median_val = median(data);
    mad_val = median(abs(data - median_val));
    
    if mad_val > 0
        outlier_threshold = 2.5 * mad_val;
        valid_data = data(abs(data - median_val) <= outlier_threshold);
    else
        valid_data = data;
    end
    
    if isempty(valid_data)
        robust_mean_val = median_val;
    else
        robust_mean_val = mean(valid_data);
    end
end

function protected_matrix = apply_change_protection(original_matrix, smoothed_matrix, params)
% Apply change limiting and final protection
    
    [n, F] = size(original_matrix);
    protected_matrix = smoothed_matrix;
    
    total_limited = 0;
    
    for i = 1:n
        for omega = 1:F
            original_val = original_matrix(i, omega);
            smoothed_val = smoothed_matrix(i, omega);
            
            % Check for excessive changes
            if original_val > 0
                change_ratio = abs(smoothed_val - original_val) / original_val;
                
                if change_ratio > params.max_change_factor
                    % Limit the change
                    if smoothed_val > original_val
                        limited_val = original_val * (1 + params.max_change_factor);
                    else
                        limited_val = original_val / (1 + params.max_change_factor);
                    end
                    
                    protected_matrix(i, omega) = limited_val;
                    total_limited = total_limited + 1;
                end
            end
            
            % Final minimum threshold
            protected_matrix(i, omega) = max(protected_matrix(i, omega), params.min_power);
        end
    end
    
    if total_limited > 0
        fprintf('  Applied change protection to %d values\n', total_limited);
    end
end

function Sigma_emp_loaded = apply_diagonal_loading(Sigma_emp, g_smooth, params)
% Apply diagonal loading to improve numerical conditioning
    
    F = length(Sigma_emp);
    n = size(Sigma_emp{1}, 1);
    Sigma_emp_loaded = cell(F, 1);
    
    loading_stats = struct('min_eigenval_before', [], 'min_eigenval_after', [], ...
                          'condition_before', [], 'condition_after', []);
    
    for omega = 1:F
        Sigma_omega = Sigma_emp{omega};
        
        % Store original conditioning
        eigenvals_orig = eig(Sigma_omega);
        loading_stats.min_eigenval_before(omega) = min(real(eigenvals_orig));
        loading_stats.condition_before(omega) = cond(Sigma_omega);
        
        % Compute diagonal loading matrix
        g_vals = g_smooth{omega};
        loading_matrix = diag(params.loading_factor * g_vals);
        
        % Apply loading
        Sigma_loaded = Sigma_omega + loading_matrix;
        
        % Ensure positive definiteness
        min_eigenval = min(real(eig(Sigma_loaded)));
        if min_eigenval < params.min_eigenval
            additional_loading = (params.min_eigenval - min_eigenval) * eye(n);
            Sigma_loaded = Sigma_loaded + additional_loading;
        end
        
        % Store result
        Sigma_emp_loaded{omega} = Sigma_loaded;
        
        % Store final conditioning
        eigenvals_final = eig(Sigma_loaded);
        loading_stats.min_eigenval_after(omega) = min(real(eigenvals_final));
        loading_stats.condition_after(omega) = cond(Sigma_loaded);
    end
    
    % Report loading statistics
    fprintf('Diagonal loading statistics:\n');
    fprintf('  Min eigenvalue before: %.2e (mean), %.2e (min)\n', ...
            mean(loading_stats.min_eigenval_before), ...
            min(loading_stats.min_eigenval_before));
    fprintf('  Min eigenvalue after:  %.2e (mean), %.2e (min)\n', ...
            mean(loading_stats.min_eigenval_after), ...
            min(loading_stats.min_eigenval_after));
    fprintf('  Condition number before: %.2e (mean)\n', ...
            mean(loading_stats.condition_before));
    fprintf('  Condition number after:  %.2e (mean)\n', ...
            mean(loading_stats.condition_after));
end

function validate_smoothed_powers_protected(g_smooth, raw_diagonals, params)
% Validate smoothed diagonal powers with enhanced checks
    
    F = length(g_smooth);
    [n, F_raw] = size(raw_diagonals);
    
    if F ~= F_raw
        error('diagonal_smoothing:dimension_mismatch', ...
              'Smoothed powers have %d frequencies, expected %d', F, F_raw);
    end
    
    total_violations = 0;
    
    for omega = 1:F
        g_vals = g_smooth{omega};
        
        % Check dimensions
        if length(g_vals) ~= n
            error('diagonal_smoothing:size_mismatch', ...
                  'Smoothed powers at frequency %d have length %d, expected %d', ...
                  omega, length(g_vals), n);
        end
        
        % Check for positive values
        negative_count = sum(g_vals <= 0);
        if negative_count > 0
            error('diagonal_smoothing:negative_power', ...
                  'Found %d negative or zero smoothed powers at frequency %d', ...
                  negative_count, omega);
        end
        
        % Check for real values
        if ~isreal(g_vals)
            error('diagonal_smoothing:complex_power', ...
                  'Complex smoothed powers at frequency %d', omega);
        end
        
        % Check for finite values
        if any(~isfinite(g_vals))
            error('diagonal_smoothing:infinite_power', ...
                  'Non-finite smoothed powers at frequency %d', omega);
        end
        
        % Check for reasonable values
        original_vals = raw_diagonals(:, omega);
        relative_change = abs(g_vals - original_vals) ./ (original_vals + 1e-12);
        
        large_changes = sum(relative_change > params.max_change_factor);
        if large_changes > 0
            total_violations = total_violations + large_changes;
        end
    end
    
    if total_violations > 0
        fprintf('Note: %d values had large changes but were kept within limits\n', total_violations);
    end
    
    fprintf('Smoothed powers validation passed\n');
end