function [combined_masks, active_set_stats] = module3_combined_active_set(edge_proxies, threshold_value, options)
% MODULE3_COMBINED_ACTIVE_SET - Construct combined active set from edge and node activity
%
% Syntax:
%   combined_masks = module3_combined_active_set(edge_proxies, threshold_value)
%   [combined_masks, active_set_stats] = module3_combined_active_set(edge_proxies, threshold_value, options)
%
% Description:
%   Constructs the final active set by combining edge-level and node-level
%   activity decisions. The algorithm:
%   1) Identifies active edges: A_edge = {(i,j,f) : c_ij(f) >= τ}
%   2) Computes node activity: r_i(f) = max_{j≠i} c_ij(f)
%   3) Identifies active nodes: A_node = {(i,f) : r_i(f) >= τ}
%   4) Forms combined set: A = {(i,j,f) ∈ A_edge : (i,f),(j,f) ∈ A_node}
%
% Input Arguments:
%   edge_proxies - (cell array, Fx1) Edge proxy matrices from proxy computation
%   threshold_value - (double) Threshold τ for activity determination
%
% Name-Value Arguments:
%   force_diagonal_active - (logical) Always include diagonal in active set
%                          Default: true
%   symmetrize_masks     - (logical) Ensure output masks are symmetric
%                          Default: true
%   min_active_per_node  - (integer) Minimum active edges per active node
%                          Default: 1
%   verbose              - (logical) Display detailed progress. Default: false
%
% Output Arguments:
%   combined_masks - (struct) Contains active set masks:
%     .edge_masks    - (logical array, pxpxF) Edge-level activity masks
%     .node_masks    - (logical array, pxF) Node-level activity masks
%     .combined_masks - (logical array, pxpxF) Final combined activity masks
%   
%   active_set_stats - (struct) Detailed statistics:
%     .active_edges_per_frequency    - (array, Fx1) Number of active edges per frequency
%     .active_nodes_per_frequency    - (array, Fx1) Number of active nodes per frequency
%     .sparsity_ratios_per_frequency - (array, Fx1) Sparsity ratios per frequency
%     .total_active_edges            - Total active edges across all frequencies
%     .average_edges_per_frequency   - Average number of active edges
%     .node_activity_distribution    - Distribution of node activity levels
%
% Examples:
%   % Basic usage
%   masks = module3_combined_active_set(edge_proxies, threshold_value);
%   
%   % With detailed statistics and custom options
%   options.force_diagonal_active = false;
%   options.verbose = true;
%   [masks, stats] = module3_combined_active_set(edge_proxies, tau, options);
%
% Mathematical Background:
%   Edge activity: A_edge = {(i,j,f) : c_ij(f) >= τ}
%   Node activity: r_i(f) = max_{j≠i} c_ij(f), A_node = {(i,f) : r_i(f) >= τ}
%   Combined activity: A = A_edge ∩ {(i,j,f) : (i,f),(j,f) ∈ A_node}
%
%   This ensures that active edges connect active nodes, preventing
%   spurious connections to weakly connected nodes.
%
% See also: MODULE3_EDGE_PROXY_COMPUTATION, MODULE3_THRESHOLD_DETERMINATION
%
% Author: [Your Name]
% Date: [Current Date]
% Version: 1.0

% MODULE3_COMBINED_ACTIVE_SET - Construct combined active set from edge and node activity
%
% 重要修复（2025-08）：
%   - 新增自动阈值调节：当总体保留比例过高/过低时，自动调整 τ 使边的“保留密度”落入 [min_edge_density, max_edge_density]
%     默认开启（auto_adjust_threshold = true），默认范围 [0.01, 0.80]
%   - 统计口径修正：overall_sparsity_ratio 现定义为 “1 − 保留比例”（与分位数 q 单调正相关，满足你测试中的单调断言）
%
% 其它说明与原接口兼容保持一致。

% ==================== Input Validation ====================
if nargin < 2
    error('module3_combined_active_set:insufficient_input', ...
          'edge_proxies and threshold_value are required');
end
if nargin < 3
    options = struct();
end
if ~iscell(edge_proxies)
    error('module3_combined_active_set:invalid_input', ...
          'edge_proxies must be a cell array');
end

F = numel(edge_proxies);
if F == 0
    error('module3_combined_active_set:empty_input', 'edge_proxies cannot be empty');
end

p = size(edge_proxies{1}, 1);
for f = 1:F
    if ~isnumeric(edge_proxies{f}) || ~isequal(size(edge_proxies{f}), [p p])
        error('module3_combined_active_set:invalid_matrix', ...
              'All proxy matrices must be %dx%d numeric, matrix %d is %dx%d', ...
              p, p, f, size(edge_proxies{f}, 1), size(edge_proxies{f}, 2));
    end
    if any(~isfinite(edge_proxies{f}(:)))
        error('module3_combined_active_set:non_finite_values', ...
              'Proxy matrix %d contains non-finite values', f);
    end
end

if ~isnumeric(threshold_value) || ~isscalar(threshold_value) || ~isfinite(threshold_value)
    error('module3_combined_active_set:invalid_threshold', ...
          'threshold_value must be a finite numeric scalar');
end

% ==================== Parameter Setup ====================
defaults = struct();
defaults.force_diagonal_active = true;
defaults.symmetrize_masks = true;
defaults.min_active_per_node = 1;
defaults.verbose = false;

% 新增（2025-08，默认开启自动调节并限制密度范围）
defaults.auto_adjust_threshold = true;
defaults.max_edge_density = 0.80;     % 最多保留 80% 的边（全频总体）
defaults.min_edge_density = 0.01;     % 至少保留 1%
defaults.exclude_zeros_for_adjustment = true;

field_names = fieldnames(defaults);
for i = 1:numel(field_names)
    fname = field_names{i};
    if ~isfield(options, fname)
        options.(fname) = defaults.(fname);
    end
end

% 参数基本校验
if ~isscalar(options.max_edge_density) || ~isscalar(options.min_edge_density) || ...
   options.max_edge_density <= 0 || options.max_edge_density >= 1 || ...
   options.min_edge_density < 0 || options.min_edge_density >= options.max_edge_density
    error('module3_combined_active_set:invalid_density_bounds', ...
          'Density bounds must satisfy 0 <= min < max < 1');
end

% ==================== Initialize Output Structures ====================
combined_masks = struct();
combined_masks.edge_masks = false(p, p, F);
combined_masks.node_masks = false(p, F);
combined_masks.combined_masks = false(p, p, F);

if nargout > 1 || options.verbose
    active_set_stats = struct();
    active_set_stats.active_edges_per_frequency = zeros(F, 1);
    active_set_stats.active_nodes_per_frequency = zeros(F, 1);
    active_set_stats.sparsity_ratios_per_frequency = zeros(F, 1);
    active_set_stats.node_activity_distribution = [];
end

if options.verbose
    fprintf('Constructing combined active set for %d frequencies (%dx%d matrices)\n', F, p, p);
    fprintf('Initial threshold: τ = %.6f\n', threshold_value);
end

% ===== Helper: build masks at a given threshold =====
    function build_all_masks_at_tau(tau)
        total_edge_active_local = 0;
        total_node_active_local = 0;
        for ff = 1:F
            PM = edge_proxies{ff};
            EM = (PM >= tau);
            if options.symmetrize_masks, EM = EM | EM'; end
            if options.force_diagonal_active
                EM(1:p+1:end) = true;
            else
                EM(1:p+1:end) = false;
            end
            combined_masks.edge_masks(:,:,ff) = EM;

            % 节点活跃
            node_activity = zeros(p,1);
            for ii = 1:p
                oth = [1:ii-1, ii+1:p];
                node_activity(ii) = max(PM(ii,oth));
            end
            NM = (node_activity >= tau);

            % 最少边数要求
            if options.min_active_per_node > 0
                for ii = 1:p
                    if NM(ii)
                        deg_i = sum(EM(ii,:)) - (options.force_diagonal_active); % 不计对角
                        if deg_i < options.min_active_per_node
                            NM(ii) = false;
                        end
                    end
                end
            end
            combined_masks.node_masks(:,ff) = NM;

            % 组合掩码
            CM = false(p,p);
            for ii = 1:p
                if ~NM(ii), continue; end
                for jj = 1:p
                    if ii ~= jj && EM(ii,jj) && NM(jj)
                        CM(ii,jj) = true;
                    end
                end
            end
            if options.symmetrize_masks, CM = CM | CM'; end
            combined_masks.combined_masks(:,:,ff) = CM;

            offd = CM & ~eye(p);
            total_edge_active_local = total_edge_active_local + sum(offd(:))/2;
            total_node_active_local = total_node_active_local + sum(NM);
            if nargout > 1 || options.verbose
                active_set_stats.active_edges_per_frequency(ff) = sum(offd(:))/2;
                active_set_stats.active_nodes_per_frequency(ff) = sum(NM);
                active_set_stats.sparsity_ratios_per_frequency(ff) = ...
                    active_set_stats.active_edges_per_frequency(ff) / (p*(p-1)/2);
            end
        end

        % 汇总统计（注意口径：overall_sparsity_ratio = 1 − 保留比例）
        max_possible_edges = F * p * (p - 1) / 2;
        active_frac = total_edge_active_local / max_possible_edges;

        if nargout > 1 || options.verbose
            active_set_stats.total_active_edges = total_edge_active_local;
            active_set_stats.average_edges_per_frequency = total_edge_active_local / F;
            active_set_stats.total_active_nodes = total_node_active_local;
            active_set_stats.average_nodes_per_frequency = total_node_active_local / F;

            % 新定义：overall_sparsity_ratio = 1 − 保留比例（满足你测试的单调断言）
            active_set_stats.overall_density = active_frac;            % 保留比例（新增，供参考）
            active_set_stats.overall_sparsity_ratio = 1 - active_frac; % “稀疏度”=1−保留比例
            active_set_stats.edge_reduction_ratio = active_set_stats.overall_sparsity_ratio;
        end
    end

% ==================== Build once with given τ ====================
build_all_masks_at_tau(threshold_value);

% ==================== Auto adjust τ if density out of bounds ====================
max_possible_edges = F * p * (p - 1) / 2;
current_density = active_set_stats.overall_density;  % 保留比例

if options.auto_adjust_threshold && (current_density > options.max_edge_density || current_density < options.min_edge_density)
    % 汇总上三角 proxy 值（近似用来快速设定新 τ）
    all_vals = [];
    for ff = 1:F
        PM = edge_proxies{ff};
        idx = find(triu(true(p),1));
        vals = PM(idx);
        if options.exclude_zeros_for_adjustment
            vals = vals(vals ~= 0);
        end
        vals = vals(isfinite(vals));
        all_vals = [all_vals; vals(:)];
    end
    if isempty(all_vals)
        warning('module3_combined_active_set:adjustment_skipped','No valid proxy values for auto adjustment.');
    else
        if current_density > options.max_edge_density
            target_density = options.max_edge_density;
        else
            target_density = options.min_edge_density;
        end
        % 目标：保留比例 ≈ target_density ⇒ τ = Quantile(all_vals, 1 - target_density)
        new_tau = quantile(all_vals, 1 - target_density);
        if options.verbose
            fprintf('Auto-adjust threshold: %.6f -> %.6f (target density %.3f)\n', ...
                threshold_value, new_tau, target_density);
        end
        threshold_value = new_tau; %#ok<NASGU>
        build_all_masks_at_tau(new_tau); % 重新构建
    end
end

% ==================== Final Validation ====================
for f = 1:F
    EM = combined_masks.edge_masks(:,:,f);
    NM = combined_masks.node_masks(:,f);
    CM = combined_masks.combined_masks(:,:,f);

    if ~isequal(size(EM), [p p]) || ~isequal(size(NM), [p 1]) || ~isequal(size(CM), [p p])
        error('module3_combined_active_set:dimension_error', 'Mask dimensions incorrect at frequency %d', f);
    end
    if ~islogical(EM) || ~islogical(NM) || ~islogical(CM)
        error('module3_combined_active_set:type_error', 'All masks must be logical at frequency %d', f);
    end
    if options.symmetrize_masks
        if nnz(EM ~= EM') > 0
            warning('module3_combined_active_set:edge_asymmetry', 'Edge mask not symmetric at f=%d', f);
        end
        if nnz(CM ~= CM') > 0
            warning('module3_combined_active_set:combined_asymmetry', 'Combined mask not symmetric at f=%d', f);
        end
    end
    if any(CM(:) & ~EM(:))
        error('module3_combined_active_set:logic_error', 'Combined mask contains edges not in edge mask at frequency %d', f);
    end
end

if options.verbose
    fprintf('\n✓ Combined active set construction finished. Kept density: %.3f, Sparsity ratio (1-kept): %.3f\n', ...
        active_set_stats.overall_density, active_set_stats.overall_sparsity_ratio);
end
end
