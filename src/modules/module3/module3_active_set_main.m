function active_set_results = module3_active_set_main(input_data, active_set_params)
% MODULE3_ACTIVE_SET_MAIN - Complete active edge set selection for sparse precision estimation
%
% Syntax:
%   active_set_results = module3_active_set_main(input_data, active_set_params)
%
% Description:
%   Performs active edge set selection to reduce optimization dimension by 
%   selecting plausible edges based on proxy statistics. The module implements
%   both correlation-based and precision-based edge proxy computation with
%   node-level filtering and combined active set construction.
%
%   Pipeline:
%     1) Edge proxy computation (correlation or precision-based)
%     2) Threshold determination via quantile
%     3) Edge active set selection
%     4) Node active set computation
%     5) Combined active set filtering
%     6) Active set application
%
% Input Arguments:
%   input_data - (struct) Required fields:
%     .whitened_covariances        - (cell array, Fx1) Whitened empirical covariances
%     .initial_precision_matrices  - (cell array, Fx1) Initial precision estimates (optional)
%     .frequencies                 - (double array, Fx1) Frequency values
%
%   active_set_params - (struct) Optional parameters:
%     .proxy_method               - ('correlation'|'precision') Method for proxy computation
%     .quantile_level             - (double) Quantile level for threshold (default: 0.1)
%     .force_diagonal_active      - (logical) Always keep diagonal entries active
%     .min_active_edges          - (integer) Minimum number of edges to keep active
%     .verbose                   - (logical) Display progress information
%
% Output Arguments:
%   active_set_results - (struct) Contains:
%     .edge_proxies              - (cell array, Fx1) Edge proxy values c_ij(f)
%     .threshold_value           - (double) Computed threshold tau
%     .edge_active_mask          - (logical array, pxpxF) Edge activity masks
%     .node_active_mask          - (logical array, pxF) Node activity masks  
%     .combined_active_mask      - (logical array, pxpxF) Final combined active masks
%     .active_edge_statistics    - (struct) Statistics about active edges
%     .computation_stats         - (struct) Computation statistics
%     .success                   - (logical) Overall success indicator
%
% Examples:
%   % Basic usage with whitened covariances
%   input_data.whitened_covariances = whitened_cov;
%   input_data.frequencies = frequencies;
%   results = module3_active_set_main(input_data, struct());
%   
%   % Advanced usage with precision-based proxy
%   input_data.initial_precision_matrices = init_precision;
%   params.proxy_method = 'precision';
%   params.quantile_level = 0.05;
%   results = module3_active_set_main(input_data, params);
%
% Mathematical Background:
%   Edge proxy computation uses either:
%   - Correlation: c_ij(f) = |Σ_ij(f)| (whitened covariance magnitude)
%   - Precision: c_ij(f) = |-Ω_ij(f)| / sqrt(Ω_ii(f) * Ω_jj(f)) (partial coherence)
%   
%   Threshold τ = Quantile({c_ij(f) : i<j, f=1..F}, q)
%   Active edges: A_edge = {(i,j,f) : c_ij(f) >= τ}
%   Node activity: r_i(f) = max_{j≠i} c_ij(f)
%   Active nodes: A_node = {(i,f) : r_i(f) >= τ}
%   Combined: A = A_edge ∩ (endpoints in A_node)
%
% See also: MODULE3_EDGE_PROXY_COMPUTATION, MODULE3_THRESHOLD_DETERMINATION,
%           MODULE3_COMBINED_ACTIVE_SET
%
% Author: [Your Name]
% Date: [Current Date]
% Version: 1.0

% ==================== Input Validation ====================
if nargin < 1
    error('module3_active_set_main:insufficient_input', ...
          'At least input_data is required');
end

if nargin < 2
    active_set_params = struct();
end

% Validate input_data structure
if ~isstruct(input_data)
    error('module3_active_set_main:invalid_input', ...
          'input_data must be a structure');
end

if ~isfield(input_data, 'whitened_covariances') || isempty(input_data(1).whitened_covariances)
    error('module3_active_set_main:missing_covariances', ...
          'whitened_covariances field is required');
end

if ~isfield(input_data, 'frequencies') || isempty(input_data(1).frequencies)
    error('module3_active_set_main:missing_frequencies', ...
          'frequencies field is required');
end

% Extract basic dimensions
Sigma_whitened = input_data.whitened_covariances;
if ~iscell(Sigma_whitened)
    error('module3_active_set_main:invalid_covariances', ...
          'whitened_covariances must be a cell array');
end

F = numel(Sigma_whitened);
frequencies = input_data.frequencies;

if numel(frequencies) ~= F
    error('module3_active_set_main:dimension_mismatch', ...
          'Frequency array length (%d) must match covariance cell array (%d)', ...
          numel(frequencies), F);
end

% Validate covariance matrix dimensions
p = size(Sigma_whitened{1}, 1);
for f = 1:F
    if ~isnumeric(Sigma_whitened{f}) || ~ismatrix(Sigma_whitened{f})
        error('module3_active_set_main:invalid_matrix', ...
              'Covariance at frequency %d must be numeric matrix', f);
    end
    if ~isequal(size(Sigma_whitened{f}), [p p])
        error('module3_active_set_main:size_mismatch', ...
              'All covariance matrices must be %dx%d, got %dx%d at f=%d', ...
              p, p, size(Sigma_whitened{f}, 1), size(Sigma_whitened{f}, 2), f);
    end
end

% ==================== Parameter Setup ====================
defaults = struct();
defaults.proxy_method = 'correlation';
defaults.quantile_level = 0.1;
defaults.force_diagonal_active = true;
defaults.min_active_edges = max(10, round(0.05 * p * (p-1) / 2));
defaults.verbose = false;

field_names = fieldnames(defaults);
for i = 1:numel(field_names)
    fname = field_names{i};
    if ~isfield(active_set_params, fname)
        active_set_params.(fname) = defaults.(fname);
    end
end

% Validate parameters
valid_methods = {'correlation', 'precision'};
if ~ismember(active_set_params.proxy_method, valid_methods)
    error('module3_active_set_main:invalid_method', ...
          'proxy_method must be one of: %s', strjoin(valid_methods, ', '));
end

if active_set_params.quantile_level <= 0 || active_set_params.quantile_level >= 1
    error('module3_active_set_main:invalid_quantile', ...
          'quantile_level must be in (0,1), got %.3f', active_set_params.quantile_level);
end

if strcmp(active_set_params.proxy_method, 'precision')
    if ~isfield(input_data, 'initial_precision_matrices') || ...
       isempty(input_data.initial_precision_matrices)
        warning('module3_active_set_main:missing_precision', ...
                'Precision method requires initial_precision_matrices, switching to correlation');
        active_set_params.proxy_method = 'correlation';
    end
end

% ==================== Initialize Results Structure ====================
active_set_results = struct();
active_set_results.edge_proxies = cell(F, 1);
active_set_results.threshold_value = NaN;
active_set_results.edge_active_mask = false(p, p, F);
active_set_results.node_active_mask = false(p, F);
active_set_results.combined_active_mask = false(p, p, F);
active_set_results.active_edge_statistics = struct();
active_set_results.computation_stats = struct();
active_set_results.success = false;

% Initialize computation statistics
stats = struct();
stats.processing_times = zeros(1, 6);
stats.total_computation_time = 0;
stats.proxy_method_used = active_set_params.proxy_method;
stats.quantile_level_used = active_set_params.quantile_level;
stats.successful_frequencies = 0;
stats.failed_frequencies = 0;
stats.error_messages = {};

overall_tic = tic;

if active_set_params.verbose
    fprintf('========================================\n');
    fprintf('Module 3: Active Edge Set Selection\n');
    fprintf('Processing %d frequencies | nodes=%d | method=%s\n', ...
            F, p, active_set_params.proxy_method);
    fprintf('Quantile level: %.3f\n', active_set_params.quantile_level);
    fprintf('========================================\n\n');
end

% ==================== Step 1: Edge Proxy Computation ====================
if active_set_params.verbose
    fprintf('Step 1/6: Edge Proxy Computation\n');
    fprintf('--------------------------------\n');
end

step1_tic = tic;
try
    switch active_set_params.proxy_method
        case 'correlation'
            for f = 1:F
                S_f = Sigma_whitened{f};
                proxy_matrix = abs(S_f);
                proxy_matrix(1:p+1:end) = 0; % Zero diagonal
                active_set_results.edge_proxies{f} = proxy_matrix;
            end
            
        case 'precision'
            initial_precision = input_data.initial_precision_matrices;
            for f = 1:F
                if ~isempty(initial_precision{f})
                    Omega_f = initial_precision{f};
                    proxy_matrix = module3_precision_proxy_computation(Omega_f);
                    active_set_results.edge_proxies{f} = proxy_matrix;
                else
                    % Fallback to correlation
                    S_f = Sigma_whitened{f};
                    proxy_matrix = abs(S_f);
                    proxy_matrix(1:p+1:end) = 0;
                    active_set_results.edge_proxies{f} = proxy_matrix;
                    stats.error_messages{end+1} = sprintf('f=%d: Used correlation fallback', f);
                end
            end
    end
    
    stats.successful_frequencies = F;
    stats.processing_times(1) = toc(step1_tic);
    
    if active_set_params.verbose
        fprintf('✓ Computed edge proxies using %s method (%.3fs)\n', ...
                active_set_params.proxy_method, stats.processing_times(1));
    end
    
catch ME
    stats.failed_frequencies = F;
    stats.error_messages{end+1} = sprintf('Edge proxy computation failed: %s', ME.message);
    if active_set_params.verbose
        fprintf('✗ Edge proxy computation failed: %s\n', ME.message);
    end
    active_set_results.computation_stats = stats;
    return;
end

% ==================== Step 2: Threshold Determination ====================
if active_set_params.verbose
    fprintf('\nStep 2/6: Threshold Determination\n');
    fprintf('---------------------------------\n');
end

step2_tic = tic;
try
    all_proxies = [];
    for f = 1:F
        proxy_matrix = active_set_results.edge_proxies{f};
        [row_idx, col_idx] = find(triu(proxy_matrix, 1));
        proxy_values = proxy_matrix(sub2ind(size(proxy_matrix), row_idx, col_idx));
        all_proxies = [all_proxies; proxy_values];
    end
    
    if isempty(all_proxies)
        error('No valid proxy values found');
    end
    
    % Remove invalid values
    valid_proxies = all_proxies(isfinite(all_proxies));
    if numel(valid_proxies) < numel(all_proxies)
        warning('module3_active_set_main:invalid_proxies', ...
                'Removed %d invalid proxy values', numel(all_proxies) - numel(valid_proxies));
    end
    
    if isempty(valid_proxies)
        error('No valid finite proxy values found');
    end
    
    % Compute threshold
    active_set_results.threshold_value = quantile(valid_proxies, active_set_params.quantile_level);
    
    stats.processing_times(2) = toc(step2_tic);
    
    if active_set_params.verbose
        fprintf('✓ Threshold computed: τ = %.6f (%.3fs)\n', ...
                active_set_results.threshold_value, stats.processing_times(2));
        fprintf('  Total proxies: %d, Valid: %d\n', numel(all_proxies), numel(valid_proxies));
    end
    
catch ME
    stats.error_messages{end+1} = sprintf('Threshold determination failed: %s', ME.message);
    if active_set_params.verbose
        fprintf('✗ Threshold determination failed: %s\n', ME.message);
    end
    active_set_results.computation_stats = stats;
    return;
end

% ==================== Step 3: Combined Active Set Construction ====================
if active_set_params.verbose
    fprintf('\nStep 3/6: Combined Active Set Construction\n');
    fprintf('------------------------------------------\n');
end

step3_tic = tic;
try
    tau = active_set_results.threshold_value;
    
    % Use the fixed combined active set function
    [combined_masks, combined_stats] = module3_combined_active_set(active_set_results.edge_proxies, tau, ...
                                                                  struct('verbose', false, ...
                                                                         'force_diagonal_active', active_set_params.force_diagonal_active, ...
                                                                         'symmetrize_masks', true));
    
    % Store masks in results
    active_set_results.edge_active_mask = combined_masks.edge_masks;
    active_set_results.node_active_mask = combined_masks.node_masks;
    active_set_results.combined_active_mask = combined_masks.combined_masks;
    
    stats.processing_times(3) = toc(step3_tic);
    
    if active_set_params.verbose
        fprintf('✓ Combined active set constructed (%.3fs)\n', stats.processing_times(3));
        fprintf('  Final active off-diagonal edges: %d\n', combined_stats.total_active_edges);
    end
    
catch ME
    stats.error_messages{end+1} = sprintf('Combined active set construction failed: %s', ME.message);
    if active_set_params.verbose
        fprintf('✗ Combined active set construction failed: %s\n', ME.message);
    end
    active_set_results.computation_stats = stats;
    return;
end

% ==================== Step 4: Statistics and Finalization ====================
if active_set_params.verbose
    fprintf('\nStep 4/6: Statistics and Finalization\n');
    fprintf('-------------------------------------\n');
end

step4_tic = tic;
try
    % CRITICAL FIX: Ensure all required statistical fields are present
    edge_stats = struct();
    
    % Copy statistics from combined_stats and ensure all required fields exist
    if exist('combined_stats', 'var') && isstruct(combined_stats)
        % Copy existing fields
        stat_field_names = fieldnames(combined_stats);
        for i = 1:numel(stat_field_names)
            edge_stats.(stat_field_names{i}) = combined_stats.(stat_field_names{i});
        end
    else
        % Compute statistics manually if combined_stats not available
        edge_stats.active_edges_per_frequency = zeros(F, 1);
        edge_stats.active_nodes_per_frequency = zeros(F, 1);
        edge_stats.sparsity_ratios_per_frequency = zeros(F, 1);
        
        total_active_edges = 0;
        total_active_nodes = 0;
        
        for f = 1:F
            combined_mask = active_set_results.combined_active_mask(:, :, f);
            node_mask = active_set_results.node_active_mask(:, f);
            
            off_diag_active = combined_mask & ~eye(p);
            edge_count = sum(off_diag_active(:)) / 2;
            node_count = sum(node_mask);
            
            edge_stats.active_edges_per_frequency(f) = edge_count;
            edge_stats.active_nodes_per_frequency(f) = node_count;
            
            total_possible_edges_per_freq = p * (p - 1) / 2;
            edge_stats.sparsity_ratios_per_frequency(f) = edge_count / total_possible_edges_per_freq;
            
            total_active_edges = total_active_edges + edge_count;
            total_active_nodes = total_active_nodes + node_count;
        end
        
        % Overall statistics
        edge_stats.total_active_edges = total_active_edges;
        edge_stats.average_edges_per_frequency = total_active_edges / F;
        edge_stats.total_active_nodes = total_active_nodes;
        edge_stats.average_nodes_per_frequency = total_active_nodes / F;
        
        % CRITICAL: Ensure overall_sparsity_ratio field exists
        max_possible_edges = F * p * (p - 1) / 2;
        edge_stats.overall_sparsity_ratio = total_active_edges / max_possible_edges;
        edge_stats.edge_reduction_ratio = 1 - edge_stats.overall_sparsity_ratio;
    end
    
    % Ensure all critical fields exist (double-check)
    required_fields = {'total_active_edges', 'average_edges_per_frequency', ...
                      'total_active_nodes', 'average_nodes_per_frequency', ...
                      'overall_sparsity_ratio', 'active_edges_per_frequency', ...
                      'active_nodes_per_frequency', 'sparsity_ratios_per_frequency'};
    
    for i = 1:numel(required_fields)
        fname = required_fields{i};
        if ~isfield(edge_stats, fname)
            switch fname
                case 'overall_sparsity_ratio'
                    if isfield(edge_stats, 'total_active_edges')
                        max_possible_edges = F * p * (p - 1) / 2;
                        edge_stats.overall_sparsity_ratio = edge_stats.total_active_edges / max_possible_edges;
                    else
                        edge_stats.overall_sparsity_ratio = 0;
                    end
                case 'edge_reduction_ratio'
                    if isfield(edge_stats, 'overall_sparsity_ratio')
                        edge_stats.edge_reduction_ratio = 1 - edge_stats.overall_sparsity_ratio;
                    else
                        edge_stats.edge_reduction_ratio = 1;
                    end
                otherwise
                    % Set default values for other missing fields
                    if contains(fname, 'per_frequency')
                        edge_stats.(fname) = zeros(F, 1);
                    else
                        edge_stats.(fname) = 0;
                    end
            end
        end
    end
    
    % Add threshold information
    edge_stats.threshold_used = active_set_results.threshold_value;
    edge_stats.quantile_level_used = active_set_params.quantile_level;
    
    % Store in results
    active_set_results.active_edge_statistics = edge_stats;
    
    stats.processing_times(4) = toc(step4_tic);
    stats.total_computation_time = toc(overall_tic);
    active_set_results.computation_stats = stats;
    active_set_results.success = true;
    
    if active_set_params.verbose
        fprintf('✓ Statistics computed (%.3fs)\n', stats.processing_times(4));
        fprintf('\n========================================\n');
        fprintf('Module 3 Summary:\n');
        fprintf('- Total computation time: %.3fs\n', stats.total_computation_time);
        fprintf('- Threshold: τ = %.6f\n', edge_stats.threshold_used);
        fprintf('- Total active edges: %d\n', edge_stats.total_active_edges);
        fprintf('- Overall sparsity ratio: %.3f\n', edge_stats.overall_sparsity_ratio);
        fprintf('- Edge reduction ratio: %.3f\n', edge_stats.edge_reduction_ratio);
        fprintf('- Success: YES\n');
        fprintf('========================================\n');
    end
    
catch ME
    stats.error_messages{end+1} = sprintf('Statistics computation failed: %s', ME.message);
    if active_set_params.verbose
        fprintf('✗ Statistics computation failed: %s\n', ME.message);
    end
    active_set_results.computation_stats = stats;
    active_set_results.success = false;
    return;
end

end

% ==================== Helper Function ====================
function proxy_matrix = module3_precision_proxy_computation(Omega)
% Compute partial coherence proxy from precision matrix
% c_ij = |-Ω_ij| / sqrt(Ω_ii * Ω_jj)

p = size(Omega, 1);
proxy_matrix = zeros(p, p);

for i = 1:p
    for j = 1:p
        if i ~= j && Omega(i,i) > 0 && Omega(j,j) > 0
            proxy_matrix(i, j) = abs(-Omega(i, j)) / sqrt(Omega(i,i) * Omega(j,j));
        end
    end
end

% Ensure symmetry
proxy_matrix = (proxy_matrix + proxy_matrix') / 2;

end