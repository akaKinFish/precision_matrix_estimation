function [threshold_value, threshold_stats] = module3_threshold_determination(edge_proxies, quantile_level, options)
% MODULE3_THRESHOLD_DETERMINATION - Determine threshold for active edge selection
%
% Syntax:
%   threshold_value = module3_threshold_determination(edge_proxies, quantile_level)
%   [threshold_value, threshold_stats] = module3_threshold_determination(edge_proxies, quantile_level, options)
%
% Description:
%   Determines the threshold τ for active edge selection by computing the
%   q-quantile of all edge proxy values across frequencies. Only considers
%   upper triangular elements (i < j) to avoid counting each edge twice.
%
%   Threshold computation:
%   τ = Quantile({c_ij(f) : i<j, f=1..F}, q)
%
% Input Arguments:
%   edge_proxies - (cell array, Fx1) Edge proxy matrices from edge proxy computation
%   quantile_level - (double) Quantile level q ∈ (0,1) for threshold determination
%
% Name-Value Arguments:
%   exclude_zeros     - (logical) Exclude zero values from threshold computation
%                       Default: true
%   min_threshold     - (double) Minimum allowable threshold value. Default: 0
%   max_threshold     - (double) Maximum allowable threshold value. Default: Inf
%   robust_quantile   - (logical) Use robust quantile estimation. Default: false
%   verbose           - (logical) Display detailed statistics. Default: false
%
% Output Arguments:
%   threshold_value - (double) Computed threshold τ
%   threshold_stats - (struct) Detailed statistics about threshold computation:
%     .total_proxy_values    - Total number of proxy values considered
%     .valid_proxy_values    - Number of finite, non-excluded values
%     .excluded_values       - Number of excluded values (zeros, NaN, Inf)
%     .proxy_range           - [min, max] of valid proxy values
%     .quantile_level_used   - Actual quantile level used
%     .threshold_percentile  - What percentile the threshold represents
%
% Examples:
%   % Basic threshold determination
%   tau = module3_threshold_determination(edge_proxies, 0.1);
%   
%   % Advanced usage with statistics and options
%   options.exclude_zeros = false;
%   options.verbose = true;
%   [tau, stats] = module3_threshold_determination(edge_proxies, 0.05, options);
%
% Mathematical Background:
%   The threshold τ is chosen to retain approximately (1-q)×100% of all edges
%   for further optimization. Lower quantile levels (e.g., q=0.05) result in
%   more aggressive pruning, while higher levels (e.g., q=0.2) are more
%   conservative.
%
% See also: MODULE3_EDGE_PROXY_COMPUTATION, MODULE3_ACTIVE_SET_MAIN
%
% Author: [Your Name]  
% Date: [Current Date]
% Version: 1.0

% ==================== Input Validation ====================
if nargin < 2
    error('module3_threshold_determination:insufficient_input', ...
          'edge_proxies and quantile_level are required');
end

if nargin < 3
    options = struct();
end

% Validate edge_proxies
if ~iscell(edge_proxies)
    error('module3_threshold_determination:invalid_input', ...
          'edge_proxies must be a cell array');
end

F = numel(edge_proxies);
if F == 0
    error('module3_threshold_determination:empty_input', ...
          'edge_proxies cannot be empty');
end

% Check that all matrices are valid
p = size(edge_proxies{1}, 1);
for f = 1:F
    if ~isnumeric(edge_proxies{f}) || ~isequal(size(edge_proxies{f}), [p p])
        error('module3_threshold_determination:invalid_matrix', ...
              'All proxy matrices must be %dx%d numeric, matrix %d is %dx%d', ...
              p, p, f, size(edge_proxies{f}, 1), size(edge_proxies{f}, 2));
    end
end

% Validate quantile_level
if ~isnumeric(quantile_level) || ~isscalar(quantile_level) || ...
   quantile_level <= 0 || quantile_level >= 1
    error('module3_threshold_determination:invalid_quantile', ...
          'quantile_level must be a scalar in (0,1), got %.3f', quantile_level);
end

% ==================== Parameter Setup ====================
defaults = struct();
defaults.exclude_zeros = true;
defaults.min_threshold = 0;
defaults.max_threshold = inf;
defaults.robust_quantile = false;
defaults.verbose = false;

field_names = fieldnames(defaults);
for i = 1:numel(field_names)
    fname = field_names{i};
    if ~isfield(options, fname)
        options.(fname) = defaults.(fname);
    end
end

% Validate parameters
if ~isnumeric(options.min_threshold) || ~isscalar(options.min_threshold)
    error('module3_threshold_determination:invalid_min_threshold', ...
          'min_threshold must be a numeric scalar');
end

if ~isnumeric(options.max_threshold) || ~isscalar(options.max_threshold)
    error('module3_threshold_determination:invalid_max_threshold', ...
          'max_threshold must be a numeric scalar');
end

if options.min_threshold >= options.max_threshold
    error('module3_threshold_determination:invalid_threshold_range', ...
          'min_threshold (%.3f) must be less than max_threshold (%.3f)', ...
          options.min_threshold, options.max_threshold);
end

% ==================== Extract Proxy Values ====================
if options.verbose
    fprintf('Extracting edge proxy values from %d frequencies (%dx%d matrices)\n', F, p, p);
end

all_proxies = [];
total_values = 0;
excluded_values = 0;

for f = 1:F
    proxy_matrix = edge_proxies{f};
    
    % Extract upper triangular part (i < j) to avoid double-counting edges
    [row_idx, col_idx] = find(triu(ones(p, p), 1));
    proxy_values = proxy_matrix(sub2ind([p, p], row_idx, col_idx));
    
    total_values = total_values + numel(proxy_values);
    
    % Apply exclusion criteria
    valid_mask = true(size(proxy_values));
    
    % Exclude non-finite values
    finite_mask = isfinite(proxy_values);
    if sum(~finite_mask) > 0
        if options.verbose
            fprintf('  f=%d: Excluding %d non-finite values\n', f, sum(~finite_mask));
        end
        excluded_values = excluded_values + sum(~finite_mask);
        valid_mask = valid_mask & finite_mask;
    end
    
    % Exclude zeros if requested
    if options.exclude_zeros
        zero_mask = (proxy_values ~= 0);
        if sum(~zero_mask) > 0
            if options.verbose
                fprintf('  f=%d: Excluding %d zero values\n', f, sum(~zero_mask));
            end
            excluded_values = excluded_values + sum(~zero_mask);
            valid_mask = valid_mask & zero_mask;
        end
    end
    
    % Apply range constraints
    range_mask = (proxy_values >= options.min_threshold) & (proxy_values <= options.max_threshold);
    if sum(~range_mask) > 0
        if options.verbose
            fprintf('  f=%d: Excluding %d values outside range [%.3f, %.3f]\n', ...
                    f, sum(~range_mask), options.min_threshold, options.max_threshold);
        end
        excluded_values = excluded_values + sum(~range_mask);
        valid_mask = valid_mask & range_mask;
    end
    
    % Add valid values to collection
    valid_values = proxy_values(valid_mask);
    all_proxies = [all_proxies; valid_values];
end

valid_count = numel(all_proxies);

if options.verbose
    fprintf('Proxy value extraction summary:\n');
    fprintf('  Total values: %d\n', total_values);
    fprintf('  Valid values: %d\n', valid_count);
    fprintf('  Excluded values: %d (%.1f%%)\n', excluded_values, ...
            100 * excluded_values / max(total_values, 1));
end

% ==================== Threshold Computation ====================
if valid_count == 0
    error('module3_threshold_determination:no_valid_values', ...
          'No valid proxy values found for threshold computation');
end

if valid_count == 1
    warning('module3_threshold_determination:single_value', ...
            'Only one valid proxy value found, using it as threshold');
    threshold_value = all_proxies(1);
else
    if options.robust_quantile
        % Use robust quantile estimation (median of bootstrap quantiles)
        n_bootstrap = 100;
        bootstrap_quantiles = zeros(n_bootstrap, 1);
        
        for b = 1:n_bootstrap
            % Bootstrap resample
            bootstrap_indices = randi(valid_count, valid_count, 1);
            bootstrap_sample = all_proxies(bootstrap_indices);
            bootstrap_quantiles(b) = quantile(bootstrap_sample, quantile_level);
        end
        
        threshold_value = median(bootstrap_quantiles);
        
        if options.verbose
            fprintf('Robust quantile estimation used (bootstrap median)\n');
            fprintf('  Bootstrap std: %.6f\n', std(bootstrap_quantiles));
        end
    else
        % Standard quantile computation
        threshold_value = quantile(all_proxies, quantile_level);
    end
end

% Apply final constraints
if threshold_value < options.min_threshold
    if options.verbose
        fprintf('Threshold %.6f below minimum, clamping to %.6f\n', ...
                threshold_value, options.min_threshold);
    end
    threshold_value = options.min_threshold;
elseif threshold_value > options.max_threshold
    if options.verbose
        fprintf('Threshold %.6f above maximum, clamping to %.6f\n', ...
                threshold_value, options.max_threshold);
    end
    threshold_value = options.max_threshold;
end

% ==================== Statistics Computation ====================
if nargout > 1 || options.verbose
    threshold_stats = struct();
    threshold_stats.total_proxy_values = total_values;
    threshold_stats.valid_proxy_values = valid_count;
    threshold_stats.excluded_values = excluded_values;
    
    if valid_count > 0
        threshold_stats.proxy_range = [min(all_proxies), max(all_proxies)];
        
        % Compute what percentile the threshold actually represents
        threshold_stats.threshold_percentile = mean(all_proxies <= threshold_value) * 100;
    else
        threshold_stats.proxy_range = [NaN, NaN];
        threshold_stats.threshold_percentile = NaN;
    end
    
    threshold_stats.quantile_level_used = quantile_level;
    
    % Additional statistics
    if valid_count > 1
        threshold_stats.proxy_mean = mean(all_proxies);
        threshold_stats.proxy_std = std(all_proxies);
        threshold_stats.proxy_median = median(all_proxies);
        
        % Compute some other useful quantiles
        other_quantiles = [0.01, 0.05, 0.25, 0.75, 0.95, 0.99];
        threshold_stats.proxy_quantiles = quantile(all_proxies, other_quantiles);
        threshold_stats.quantile_levels = other_quantiles;
    end
    
    if options.verbose
        fprintf('\nThreshold determination results:\n');
        fprintf('  Threshold: τ = %.6f\n', threshold_value);
        fprintf('  Target quantile level: %.3f\n', quantile_level);
        fprintf('  Actual percentile: %.1f%%\n', threshold_stats.threshold_percentile);
        fprintf('  Proxy range: [%.6f, %.6f]\n', threshold_stats.proxy_range(1), threshold_stats.proxy_range(2));
        
        if valid_count > 1
            fprintf('  Proxy statistics:\n');
            fprintf('    Mean: %.6f\n', threshold_stats.proxy_mean);
            fprintf('    Std:  %.6f\n', threshold_stats.proxy_std);
            fprintf('    Median: %.6f\n', threshold_stats.proxy_median);
        end
    end
end

% ==================== Final Validation ====================
if ~isfinite(threshold_value)
    error('module3_threshold_determination:invalid_threshold', ...
          'Computed threshold is not finite: %.6f', threshold_value);
end

if threshold_value < 0
    warning('module3_threshold_determination:negative_threshold', ...
            'Computed threshold is negative: %.6f', threshold_value);
end

end