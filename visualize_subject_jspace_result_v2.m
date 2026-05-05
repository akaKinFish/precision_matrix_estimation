function out = visualize_subject_jspace_result_v2(result_mat_path, scalp_mat_path, opts)
%VISUALIZE_SUBJECT_JSPACE_RESULT Visualize one subject on scalp/source layers.
%   This version keeps the user's original workflow but adds:
%   1) Linear coherence (sqrt(MSC)) as the default display metric.
%   2) Robust color scaling for sparse PCOR matrices.
%   3) Optional peak-centered band averaging for peak matrices.
%   4) Printed summary statistics (p95 / p99 / max) for MSC, linear coherence, and PCOR.
%   5) NaN-Separation trick for plotting ALL 60k+ edges instantly.
%
% Required inputs:
%   result_mat_path : MAT file containing at least Sjj_est, and optionally Omega_est.
%   scalp_mat_path  : MAT file containing data_struct.CrossM (optional if recoverable from result file).
%   opts            : optional struct with fields:
%       alpha_band              : [f1 f2], default [8 12]
%       eps_floor               : numeric, default 1e-12
%       roi_order               : custom ROI order for matrices
%       save_dir                : directory to save figures
%       SC                      : Brainstorm cortex struct for region plots / 3D plots
%       atlas_name              : atlas name, default 'HCP_MMP1'
%       edge_keep_quantile      : default 0.998
%       edge_keep_max           : default 150
%       head_surface_view       : default [110 20]
%       roi_to_scout_map        : optional ROI->scout map
%       peak_avg_bins           : half-width in bins around detected peak (default 0)
%       coh_curve_metric        : 'linear' (default) or 'msc'
%       matrix_metric           : 'linear' (default) or 'msc'
%       print_stats             : true (default)
%
% Output:
%   out struct with figure handles, peak info, and computed matrices.
    if nargin < 1 || isempty(result_mat_path)
        error('Please provide result_mat_path.');
    end
    if nargin < 2
        scalp_mat_path = [];
    end
    if nargin < 3 || isempty(opts)
        opts = struct();
    end
    alpha_band = get_opt_(opts, 'alpha_band', [8 12]);
    eps_floor = get_opt_(opts, 'eps_floor', 1e-12);
    roi_order_in = get_opt_(opts, 'roi_order', []);
    save_dir = get_opt_(opts, 'save_dir', '');
    SC = get_opt_(opts, 'SC', []);
    atlas_name = get_opt_(opts, 'atlas_name', 'HCP_MMP1');
    edge_keep_quantile = get_opt_(opts, 'edge_keep_quantile', 0.998);
    edge_keep_max = get_opt_(opts, 'edge_keep_max', 150);
    head_surface_view = get_opt_(opts, 'head_surface_view', [110 20]);
    roi_to_scout_map = get_opt_(opts, 'roi_to_scout_map', []);
    peak_avg_bins = get_opt_(opts, 'peak_avg_bins', 0);
    coh_curve_metric = lower(strtrim(get_opt_(opts, 'coh_curve_metric', 'linear')));
    matrix_metric = lower(strtrim(get_opt_(opts, 'matrix_metric', 'linear')));
    print_stats = get_opt_(opts, 'print_stats', true);
    
    % 注意：使用了 NaN Trick 后，不再需要限制抽样条数，此参数保留仅为兼容性
    global_superimpose_max = get_opt_(opts, 'global_superimpose_max', 10000);
    
    % ------------------------------------------------------------
    % 1) Load result and scalp data
    % ------------------------------------------------------------
    dres = load(result_mat_path);
    if ~isfield(dres, 'Sjj_est')
        error('Result file missing Sjj_est: %s', result_mat_path);
    end
    Sjj = to_tensor_(dres.Sjj_est);
    Sjj = sanitize_cross_tensor_(Sjj, eps_floor);
    [Nr, ~, Fsrc] = size(Sjj);
    has_omega = isfield(dres, 'Omega_est') && ~isempty(dres.Omega_est);
    if has_omega
        Omega = to_tensor_(dres.Omega_est);
        Omega = sanitize_cross_tensor_(Omega, eps_floor);
    else
        warning('Omega_est not found in result MAT. PCOR plots will be skipped.');
    end
    scalp_file = resolve_scalp_path_(dres, scalp_mat_path);
    dscalp = load(scalp_file, 'data_struct');
    if ~isfield(dscalp, 'data_struct') || ~isfield(dscalp.data_struct, 'CrossM')
        error('Scalp file missing data_struct.CrossM: %s', scalp_file);
    end
    Svv = sanitize_cross_tensor_(dscalp.data_struct.CrossM, eps_floor);
    [Nc, ~, Fscalp] = size(Svv);
    freq = get_freq_(dscalp.data_struct, dres, Fscalp);
    Fuse = min([Fsrc, Fscalp, numel(freq)]);
    Sjj = Sjj(:, :, 1:Fuse);
    Svv = Svv(:, :, 1:Fuse);
    if has_omega
        Omega = Omega(:, :, 1:Fuse);
    end
    freq = freq(1:Fuse);
    alpha_idx = find(freq >= alpha_band(1) & freq <= alpha_band(2));
    if isempty(alpha_idx)
        alpha_idx = 1:numel(freq);
    end
    % ------------------------------------------------------------
    % 2) Power spectra
    % ------------------------------------------------------------
    scalp_power = real(reshape(diag3d_(Svv), Nc, Fuse));
    source_power = real(reshape(diag3d_(Sjj), Nr, Fuse));
    scalp_power_db = 10 * log10(max(scalp_power, eps_floor));
    source_power_db = 10 * log10(max(source_power, eps_floor));
    mean_src_pow = mean(source_power_db, 1);
    [~, peak_local_idx] = max(mean_src_pow(alpha_idx));
    peak_f_idx = alpha_idx(peak_local_idx);
    peak_freq = freq(peak_f_idx);
    fprintf('\n[Peak Detection] Found Alpha peak at %.2f Hz (Index: %d)\n', peak_freq, peak_f_idx);
    peak_bins = max(1, peak_f_idx - peak_avg_bins) : min(Fuse, peak_f_idx + peak_avg_bins);
    if numel(peak_bins) == 1
        peak_desc = sprintf('%.1f Hz', peak_freq);
    else
        peak_desc = sprintf('%.1f--%.1f Hz avg', freq(peak_bins(1)), freq(peak_bins(end)));
    end
    % ------------------------------------------------------------
    % 3) Coherence / PCOR tensors
    % ------------------------------------------------------------
    scalp_coh_msc = coherence_from_cross_(Svv, eps_floor);
    source_coh_msc = coherence_from_cross_(Sjj, eps_floor);
    scalp_coh_lin = sqrt(max(scalp_coh_msc, 0));
    source_coh_lin = sqrt(max(source_coh_msc, 0));
    scalp_alpha_msc = mean(scalp_coh_msc(:, :, alpha_idx), 3);
    scalp_alpha_lin = mean(scalp_coh_lin(:, :, alpha_idx), 3);
    source_alpha_msc = mean(source_coh_msc(:, :, alpha_idx), 3);
    source_alpha_lin = mean(source_coh_lin(:, :, alpha_idx), 3);
    source_msc_peak = mean(source_coh_msc(:, :, peak_bins), 3);
    source_lin_peak = mean(source_coh_lin(:, :, peak_bins), 3);
    source_msc_peak(1:Nr+1:end) = 0;
    source_lin_peak(1:Nr+1:end) = 0;
    source_pcor_peak = zeros(Nr, Nr);
    pcor_tensor = [];
    if has_omega
        pcor_tensor = pcor_from_precision_tensor_(Omega, eps_floor);
        source_pcor_peak = mean(pcor_tensor(:, :, peak_bins), 3);
        source_pcor_peak(1:Nr+1:end) = 0;
    end
    % Choose display metric for coherence-based curves / matrices
    switch coh_curve_metric
        case 'msc'
            source_curve_tensor = source_coh_msc;
            curve_tag = 'MSC coherence';
        otherwise
            source_curve_tensor = source_coh_lin;
            curve_tag = 'Linear coherence';
    end
    switch matrix_metric
        case 'msc'
            source_display_peak = source_msc_peak;
            source_alpha_display = source_alpha_msc;
            matrix_tag = 'MSC coherence';
        otherwise
            source_display_peak = source_lin_peak;
            source_alpha_display = source_alpha_lin;
            matrix_tag = 'Linear coherence';
    end
    if isempty(roi_order_in)
        off = source_display_peak;
        [~, roi_order] = sort(sum(off, 2), 'descend');
    else
        roi_order = roi_order_in(:);
    end
    % ------------------------------------------------------------
    % 4) Print useful statistics
    % ------------------------------------------------------------
    stats = struct();
    stats.peak.msc = summarize_upper_triangle_(source_msc_peak);
    stats.peak.linear = summarize_upper_triangle_(source_lin_peak);
    if has_omega
        stats.peak.pcor = summarize_upper_triangle_(source_pcor_peak);
    else
        stats.peak.pcor = struct('p95', NaN, 'p99', NaN, 'max', NaN, 'mean', NaN, 'nnz', NaN);
    end
    stats.alpha.msc = summarize_upper_triangle_(source_alpha_msc);
    stats.alpha.linear = summarize_upper_triangle_(source_alpha_lin);
    if has_omega
        source_alpha_pcor = mean(pcor_tensor(:, :, alpha_idx), 3);
        source_alpha_pcor(1:Nr+1:end) = 0;
        stats.alpha.pcor = summarize_upper_triangle_(source_alpha_pcor);
    else
        stats.alpha.pcor = struct('p95', NaN, 'p99', NaN, 'max', NaN, 'mean', NaN, 'nnz', NaN);
    end
    if print_stats
        fprintf('[Source edge stats @ %s]\n', peak_desc);
        fprintf('  MSC coherence    : mean=%.4g | p95=%.4g | p99=%.4g | max=%.4g\n', ...
            stats.peak.msc.mean, stats.peak.msc.p95, stats.peak.msc.p99, stats.peak.msc.max);
        fprintf('  Linear coherence : mean=%.4g | p95=%.4g | p99=%.4g | max=%.4g\n', ...
            stats.peak.linear.mean, stats.peak.linear.p95, stats.peak.linear.p99, stats.peak.linear.max);
        if has_omega
            fprintf('  PCOR             : mean=%.4g | p95=%.4g | p99=%.4g | max=%.4g\n', ...
                stats.peak.pcor.mean, stats.peak.pcor.p95, stats.peak.pcor.p99, stats.peak.pcor.max);
        end
        fprintf('[Alpha-band source edge stats]\n');
        fprintf('  MSC coherence    : mean=%.4g | p95=%.4g | p99=%.4g | max=%.4g\n', ...
            stats.alpha.msc.mean, stats.alpha.msc.p95, stats.alpha.msc.p99, stats.alpha.msc.max);
        fprintf('  Linear coherence : mean=%.4g | p95=%.4g | p99=%.4g | max=%.4g\n', ...
            stats.alpha.linear.mean, stats.alpha.linear.p95, stats.alpha.linear.p99, stats.alpha.linear.max);
        if has_omega
            fprintf('  PCOR             : mean=%.4g | p95=%.4g | p99=%.4g | max=%.4g\n', ...
                stats.alpha.pcor.mean, stats.alpha.pcor.p95, stats.alpha.pcor.p99, stats.alpha.pcor.max);
        end
    end
    % ------------------------------------------------------------
    % 5) Figures
    % ------------------------------------------------------------
    figs = struct();
    % 5.1 Power spectrum
    figs.spectrum = figure('Name', 'Power Spectrum & Peak Detection', 'Position', [100 100 1000 400]);
    subplot(1,2,1);
    plot(freq, scalp_power_db.', 'Color', [0.8 0.8 0.8]); hold on;
    plot(freq, mean(scalp_power_db, 1), 'k', 'LineWidth', 2);
    ylm = ylim;
    patch([alpha_band(1) alpha_band(2) alpha_band(2) alpha_band(1)], ...
          [ylm(1) ylm(1) ylm(2) ylm(2)], [0.9 0.95 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot([peak_freq peak_freq], ylm, 'r--', 'LineWidth', 1.5);
    xlabel('Frequency (Hz)'); ylabel('Power (dB)'); title('Scalp power spectrum'); grid on;
    subplot(1,2,2);
    med_src = median(source_power_db, 1);
    fill([freq fliplr(freq)], [percentile_cols_(source_power_db, 25) fliplr(percentile_cols_(source_power_db, 75))], ...
         [0.85 0.9 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.5); hold on;
    plot(freq, med_src, 'b', 'LineWidth', 2);
    ylm = ylim;
    patch([alpha_band(1) alpha_band(2) alpha_band(2) alpha_band(1)], ...
          [ylm(1) ylm(1) ylm(2) ylm(2)], [0.9 0.95 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot([peak_freq peak_freq], ylm, 'r--', 'LineWidth', 1.5);
    plot(peak_freq, med_src(peak_f_idx), 'rp', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    xlabel('Frequency (Hz)'); ylabel('Power (dB)'); title(sprintf('Source spectrum (Peak = %.1f Hz)', peak_freq)); grid on;
    % 5.2 Peak matrices
    figs.peak_matrices = figure('Name', 'Peak Source Matrices', 'Position', [150 150 900 400]);
    subplot(1,2,1);
    Csrc = source_display_peak(roi_order, roi_order);
    Csrc(1:Nr+1:end) = NaN;
    climC = robust_upper_clim_(source_display_peak, false);
    imagesc(Csrc, [0 max(max(source_display_peak))]); axis image; colorbar;
    title(sprintf('%s at %s', matrix_tag, peak_desc));
    if has_omega
        subplot(1,2,2);
        Psrc = source_pcor_peak(roi_order, roi_order);
        Psrc(1:Nr+1:end) = NaN;
        climP = robust_upper_clim_(source_pcor_peak, true);
        imagesc(Psrc, [0 max(max(source_pcor_peak))]); axis image; colorbar;
        title(sprintf('Partial Coherence (PCOR) at %s', peak_desc));
    end
    % 5.3 ROI coherence spectrum
    seed_info = select_occipital_seed_(SC, atlas_name, Nr);
    roi_vertices = seed_info.roi_vertices;
    seed_label = seed_info.seed_label;
    if isempty(roi_vertices)
        warning('No valid ROI vertices found for coherence spectrum figure.');
    else
        roi_tensor = source_curve_tensor(roi_vertices, :, :);
        vertex_mean_curve = squeeze(mean(roi_tensor, 2));
        if isvector(vertex_mean_curve)
            vertex_mean_curve = vertex_mean_curve(:).';
        end
        if size(vertex_mean_curve, 2) ~= Fuse
            vertex_mean_curve = vertex_mean_curve.';
        end
        figs.roi_coherence_spec = figure('Name', 'ROI Vertices Coherence Spectrum', 'Position', [200 200 700 450]);
        hold on;
        ylm_max = max(vertex_mean_curve(:)) * 1.1;
        if ~isfinite(ylm_max) || ylm_max <= 0
            ylm_max = 1;
        end
        patch([alpha_band(1) alpha_band(2) alpha_band(2) alpha_band(1)], ...
              [0 0 ylm_max ylm_max], [1 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        plot([peak_freq peak_freq], [0 ylm_max], 'k--', 'LineWidth', 1.5);
        colors = lines(size(vertex_mean_curve,1));
        for v = 1:size(vertex_mean_curve,1)
            plot(freq, vertex_mean_curve(v,:), 'Color', colors(v,:), 'LineWidth', 1.5);
        end
        ylim([0 ylm_max]); xlim([freq(1), freq(end)]);
        xlabel('Frequency (Hz)'); ylabel(sprintf('Mean Global %s', curve_tag));
        title(sprintf('Coherence Spectrum (ROI: %s, %d Vertices)', seed_label, size(vertex_mean_curve,1)), 'Interpreter', 'none');
        grid on;
        drawnow; 
    end
    
    % ------------------------------------------------------------
    % 5.4 Global coherence spectrum over all upper-triangular edges
    % ------------------------------------------------------------
    maskU = triu(true(Nr), 1);
    num_edges = nnz(maskU);
    all_curve_lines = zeros(num_edges, Fuse);
    for f = 1:Fuse
        Af = source_curve_tensor(:,:,f);
        all_curve_lines(:,f) = Af(maskU);
    end
    figs.global_coherence_spec = figure('Name', 'Global Coherence (Upper Triangular)', 'Position', [250 250 700 450]);
    
    % 强制关闭该坐标轴的默认交互和悬浮工具栏
    ax_global = gca;
    disableDefaultInteractivity(ax_global);
    ax_global.Toolbar.Visible = 'off';
    hold on;
    
    ylm_max = max(all_curve_lines(:));
    ylm_max = max(ylm_max * 1.05, 1e-6);
    
    % 1. 获取 Patch 句柄
    h_patch = patch([alpha_band(1) alpha_band(2) alpha_band(2) alpha_band(1)], ...
          [0 0 ylm_max ylm_max], [1 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
          
    % 2. 获取 Peak 虚线句柄
    h_peak = plot([peak_freq peak_freq], [0 ylm_max], 'k--', 'LineWidth', 1.5);
    
    % --- NaN 截断法绘制所有 Coherence 连边 ---
    Y_mat_coh = all_curve_lines.'; 
    X_mat_coh = repmat(freq(:), 1, num_edges); 
    Y_mat_coh(end+1, :) = NaN; 
    X_mat_coh(end+1, :) = NaN; 
    Y_vec_coh = Y_mat_coh(:);
    X_vec_coh = X_mat_coh(:);
    
    % 一次性画出所有线 (只生成 1 个 Line 对象)，设置极低透明度
    h_edges = plot(X_vec_coh, Y_vec_coh, 'Color', [0.75 0.75 0.75 0.02], 'LineWidth', 0.5);
    
    % 统计特性
    global_avg_curve = mean(all_curve_lines, 1, 'omitnan');
    global_p95_curve = prctile(all_curve_lines, 95, 1);
    
    % 4. 获取 Mean 和 95th 曲线的句柄
    h_mean = plot(freq, global_avg_curve, 'r', 'LineWidth', 2.5);
    h_p95 = plot(freq, global_p95_curve, 'b--', 'LineWidth', 1.5);
    
    ylim([0 ylm_max]); xlim([freq(1), freq(end)]);
    xlabel('Frequency (Hz)'); ylabel(curve_tag);
    title(sprintf('Global %s Spectrum (All %d Edges)', curve_tag, num_edges));
    
    % 5. 将指定的句柄数组传给 legend，彻底杜绝错乱
    legend([h_patch, h_peak, h_edges, h_mean, h_p95], ...
           {'Alpha band', 'Peak', 'All edges (Density)', 'Mean', '95th pct'}, ...
           'Location', 'best');
           
    grid on;
    drawnow;
    
    % 5.5 Region/chord/3D plots
    out_region_coh = struct();
    out_region_pcor = struct();
    if ~isempty(SC)
        if exist('visualize_source_coherence_brain_regions', 'file') == 2
            region_opts = struct('atlas_name', atlas_name, ...
                                 'edge_keep_quantile', edge_keep_quantile, ...
                                 'edge_keep_max', edge_keep_max, ...
                                 'chord_mode', 'both');
            region_opts.fig_prefix = sprintf('PEAK %s', upper(matrix_tag));
            region_out_coh = visualize_source_coherence_brain_regions(source_display_peak, SC, region_opts);
            if isfield(region_out_coh, 'figs')
                if isfield(region_out_coh.figs, 'region_chord'), figs.source_region_chord_coh = region_out_coh.figs.region_chord; end
                if isfield(region_out_coh.figs, 'roi_chord'), figs.source_roi_chord_coh = region_out_coh.figs.roi_chord; end
                if isfield(region_out_coh.figs, 'region_heatmap'), figs.source_region_heatmap_coh = region_out_coh.figs.region_heatmap; end
            end
            out_region_coh = region_out_coh;
            if has_omega
                region_opts.fig_prefix = 'PEAK PCOR';
                region_out_pcor = visualize_source_coherence_brain_regions(source_pcor_peak, SC, region_opts);
                if isfield(region_out_pcor, 'figs')
                    if isfield(region_out_pcor.figs, 'region_chord'), figs.source_region_chord_pcor = region_out_pcor.figs.region_chord; end
                    if isfield(region_out_pcor.figs, 'roi_chord'), figs.source_roi_chord_pcor = region_out_pcor.figs.roi_chord; end
                    if isfield(region_out_pcor.figs, 'region_heatmap'), figs.source_region_heatmap_pcor = region_out_pcor.figs.region_heatmap; end
                end
                out_region_pcor = region_out_pcor;
            end
        end
        figs.source_3d_network_coh = figure('Name', 'Source 3D Brain Network (Coherence)');
        plot_3d_brain_network_(SC, source_display_peak, atlas_name, roi_to_scout_map, ...
                               head_surface_view, edge_keep_quantile, edge_keep_max, 1e-6, ...
                               sprintf('Source %s Network at %s', matrix_tag, peak_desc));
        if has_omega
            figs.source_3d_network_pcor = figure('Name', 'Source 3D Brain Network (PCOR)');
            plot_3d_brain_network_(SC, source_pcor_peak, atlas_name, roi_to_scout_map, ...
                                   head_surface_view, edge_keep_quantile, edge_keep_max, 1e-6, ...
                                   sprintf('Source Partial Coherence Network at %s', peak_desc));
        end
    end
    
    % ------------------------------------------------------------
    % 5.4b Global Partial Coherence (PCOR) spectrum 
    % ------------------------------------------------------------
    if has_omega && ~isempty(pcor_tensor)
        all_pcor_lines = zeros(num_edges, Fuse);
        for f = 1:Fuse
            Pf_mat = pcor_tensor(:,:,f);
            all_pcor_lines(:,f) = Pf_mat(maskU);
        end
        
        figs.global_pcor_spec = figure('Name', 'Global PCOR (Upper Triangular)', 'Position', [300 300 700 450]);
        
        % 强制关闭该坐标轴的默认交互和悬浮工具栏
        ax_pcor = gca;
        disableDefaultInteractivity(ax_pcor);
        ax_pcor.Toolbar.Visible = 'off';
        hold on;
        
        % 获取 PCOR 的最大值范围
        ylm_max_pcor = max(all_pcor_lines(:));
        ylm_max_pcor = max(ylm_max_pcor * 1.05, 1e-6);
        
        % 1. 获取 Patch 句柄
        h_patch_pcor = patch([alpha_band(1) alpha_band(2) alpha_band(2) alpha_band(1)], ...
              [0 0 ylm_max_pcor ylm_max_pcor], [1 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
              
        % 2. 获取 Peak 虚线句柄
        h_peak_pcor = plot([peak_freq peak_freq], [0 ylm_max_pcor], 'k--', 'LineWidth', 1.5);
        
        % --- NaN 截断法绘制所有 PCOR 连边 ---
        Y_mat_pcor = all_pcor_lines.'; 
        X_mat_pcor = repmat(freq(:), 1, num_edges); 
        Y_mat_pcor(end+1, :) = NaN; 
        X_mat_pcor(end+1, :) = NaN; 
        Y_vec_pcor = Y_mat_pcor(:);
        X_vec_pcor = X_mat_pcor(:);
        
        % 一次性画出所有线 (只生成 1 个 Line 对象)，设置极低透明度
        h_edges_pcor = plot(X_vec_pcor, Y_vec_pcor, 'Color', [0.75 0.75 0.75 0.02], 'LineWidth', 0.5);
        
        % 统计特性
        global_avg_pcor = mean(all_pcor_lines, 1, 'omitnan');
        global_p95_pcor = prctile(all_pcor_lines, 95, 1);
        
        % 4. 获取 Mean 和 95th 曲线的句柄
        h_mean_pcor = plot(freq, global_avg_pcor, 'r', 'LineWidth', 2.5);
        h_p95_pcor = plot(freq, global_p95_pcor, 'b--', 'LineWidth', 1.5);
        
        ylim([0 ylm_max_pcor]); xlim([freq(1), freq(end)]);
        xlabel('Frequency (Hz)'); ylabel('Partial Coherence (PCOR)');
        title(sprintf('Global PCOR Spectrum (All %d Edges)', num_edges));
        
        % 5. 精确绑定图例
        legend([h_patch_pcor, h_peak_pcor, h_edges_pcor, h_mean_pcor, h_p95_pcor], ...
               {'Alpha band', 'Peak', 'All edges (Density)', 'Mean', '95th pct'}, ...
               'Location', 'best');
               
        grid on;
        drawnow;

        % ------------------------------------------------------------
        % 5.4c Global Partial Coherence (PCOR) Histogram
        % ------------------------------------------------------------
        figs.global_pcor_hist = figure('Name', 'Global PCOR Histogram', 'Position', [350 350 700 450]);
        
        % 提取用于画直方图的数据
        % 这里把所有边、所有频点的数据展平为一维向量。
        % (提示：如果你只想看 Alpha 峰值频率处的直方图，请改为 pcor_hist_data = all_pcor_lines(:, peak_f_idx);)
        pcor_hist_data = all_pcor_lines(:);
        
        % 计算整体的均值和 95% 分位数
        hist_mean = mean(pcor_hist_data, 'omitnan');
        hist_p95  = prctile(pcor_hist_data, 95);
        
        % 绘制非常密集的直方图 (Bins = 300，无边框)
        num_bins = 300; 
        h_hist = histogram(pcor_hist_data, num_bins, 'EdgeColor', 'none', 'FaceColor', [0.3 0.6 0.8], 'FaceAlpha', 0.8);
        hold on;
        
        % 获取当前 Y 轴的高度，用于画垂直辅助线
        ylm_hist = ylim;
        
        % 绘制均值线 (红色粗实线) 和 95% 分位数线 (紫色粗虚线)
        % h_mean_line = plot([hist_mean hist_mean], ylm_hist, 'r-', 'LineWidth', 2.5);
        % h_p95_line  = plot([hist_p95 hist_p95], ylm_hist, 'm--', 'LineWidth', 2.5);
        
        % 格式化图表
        xlabel('Partial Coherence (PCOR)');
        ylabel('Frequency (Counts)');
        title(sprintf('Global PCOR Distribution (All %d Edges \\times %d Freqs)', num_edges, Fuse));
        
        % 绑定图例，并在图例中直接显示具体的数值
        legend([h_hist, h_mean_line, h_p95_line], ...
               {'PCOR Distribution', ...
                sprintf('Mean = %.4f', hist_mean), ...
                sprintf('95th Pct = %.4f', hist_p95)}, ...
               'Location', 'best');
               
        grid on;
        drawnow;
    end
    
    out = struct();
    out.figs = figs;
    out.peak_freq = peak_freq;
    out.peak_f_idx = peak_f_idx;
    out.peak_bins = peak_bins;
    out.source_msc_peak = source_msc_peak;
    out.source_mag_peak = source_lin_peak;
    out.source_pcor_peak = source_pcor_peak;
    out.source_alpha_msc = source_alpha_msc;
    out.source_alpha_linear = source_alpha_lin;
    out.roi_order = roi_order;
    out.region_coh = out_region_coh;
    out.region_pcor = out_region_pcor;
    out.stats = stats;
end
% -------------------------------------------------------------------------
% Helpers
% -------------------------------------------------------------------------
function info = select_occipital_seed_(SC, atlas_name, Nr)
    info = struct('roi_vertices', [], 'seed_label', 'N/A');
    if isempty(SC) || ~isfield(SC, 'Atlas') || isempty(SC.Atlas)
        return;
    end
    ai = 1;
    for k = 1:numel(SC.Atlas)
        if isfield(SC.Atlas(k), 'Name') && strcmpi(char(SC.Atlas(k).Name), atlas_name)
            ai = k;
            break;
        end
    end
    if ~isfield(SC.Atlas(ai), 'Scouts') || isempty(SC.Atlas(ai).Scouts)
        return;
    end
    scouts = SC.Atlas(ai).Scouts;
    occ_keywords = {'Occipital', 'V1', 'O1', 'V2', 'cuneus', 'calcarine'};
    seed_idx = [];
    for i = 1:numel(scouts)
        lbl = '';
        if isfield(scouts(i), 'Label') && ~isempty(scouts(i).Label)
            lbl = char(scouts(i).Label);
        end
        is_left = contains(lower(lbl), ' left') || contains(lower(lbl), 'l_') || startsWith(lbl, 'L_') || endsWith(lbl, ' L') || startsWith(lbl, 'L');
        if any(contains(lbl, occ_keywords, 'IgnoreCase', true)) && is_left
            seed_idx = i;
            break;
        end
    end
    if isempty(seed_idx)
        seed_idx = 1;
    end
    info.seed_label = char(scouts(seed_idx).Label);
    if isfield(scouts(seed_idx), 'Vertices') && ~isempty(scouts(seed_idx).Vertices)
        roi_vertices = scouts(seed_idx).Vertices(:);
        roi_vertices = roi_vertices(roi_vertices > 0 & roi_vertices <= Nr);
        info.roi_vertices = roi_vertices;
    end
end
function plot_3d_brain_network_(SC, C, atlas_name, roi_to_scout_map, view_ang, q_keep, max_edges, min_edge_abs, plt_title)
    if ~isfield(SC, 'Vertices') || ~isfield(SC, 'Faces')
        text(0.5, 0.5, 'SC has no surface mesh.', 'HorizontalAlignment', 'center'); axis off; return;
    end
    V = SC.Vertices; N = size(C, 1);
    ai = 1;
    if isfield(SC, 'Atlas') && ~isempty(SC.Atlas)
        for k = 1:numel(SC.Atlas)
            if strcmpi(char(SC.Atlas(k).Name), atlas_name), ai = k; break; end
        end
    end
    scouts = SC.Atlas(ai).Scouts;
    if isempty(roi_to_scout_map), roi_to_scout_map = (1:N)'; end
    centers = zeros(N, 3); valid_node = false(N, 1);
    for i = 1:min(N, numel(roi_to_scout_map))
        si = roi_to_scout_map(i);
        if si > 0 && si <= numel(scouts) && isfield(scouts(si), 'Vertices') && ~isempty(scouts(si).Vertices)
            centers(i, :) = mean(V(scouts(si).Vertices, :), 1);
            valid_node(i) = true;
        end
    end
    allw = C(triu(true(N), 1));
    allw = allw(isfinite(allw) & allw > min_edge_abs);
    if isempty(allw)
        thr = inf;
    else
        thr = percentile_(allw, q_keep * 100);
    end
    edges = [];
    for i = 1:N
        for j = i+1:N
            if valid_node(i) && valid_node(j)
                w = C(i,j);
                if isfinite(w) && w >= thr && w > min_edge_abs
                    edges = [edges; i, j, w]; %#ok<AGROW>
                end
            end
        end
    end
    if ~isempty(max_edges) && size(edges,1) > max_edges
        [~, ord] = sort(edges(:,3), 'descend');
        edges = edges(ord(1:max_edges), :);
    end
    patch('Vertices', V, 'Faces', SC.Faces, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.15);
    hold on; axis equal off vis3d;
    if isempty(edges)
        title([plt_title ' (No edges)']);
        return;
    end
    cmap = colormap(gca, 'autumn');
    wmin = min(edges(:,3)); wmax = max(edges(:,3));
    for e = 1:size(edges,1)
        i = edges(e,1); j = edges(e,2); w = edges(e,3);
        if wmax > wmin, t = (w - wmin)/(wmax - wmin); else, t = 0.5; end
        col = cmap(max(1, min(size(cmap,1), round(t * size(cmap,1)))), :);
        p1 = centers(i,:); p2 = centers(j,:);
        plot3([p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)], '-', 'Color', col, 'LineWidth', 1 + 5 * t);
    end
    active_nodes = unique(edges(:, 1:2));
    plot3(centers(active_nodes, 1), centers(active_nodes, 2), centers(active_nodes, 3), ...
          'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 5, 'LineWidth', 1);
    view(view_ang(1), view_ang(2)); camlight('headlight'); lighting gouraud; title(plt_title); hold off;
end
function ai = get_atlas_index_(SC, atlas_name)
    ai = 1;
    if isfield(SC, 'Atlas') && ~isempty(SC.Atlas)
        for k = 1:numel(SC.Atlas)
            if isfield(SC.Atlas(k), 'Name') && strcmpi(char(SC.Atlas(k).Name), atlas_name)
                ai = k;
                return;
            end
        end
    end
end
function scalp_file = resolve_scalp_path_(dres, scalp_mat_path)
    if ~isempty(scalp_mat_path), scalp_file = scalp_mat_path; return; end
    if isfield(dres, 'site_dir') && isfield(dres, 'subject_id')
        scalp_file = fullfile(char(dres.site_dir), '2pre', 'crossSpec', char(dres.subject_id), [char(dres.subject_id) '.mat']);
        if isfile(scalp_file), return; end
    end
    error('Cannot locate scalp source file.');
end
function freq = get_freq_(data_struct, dres, F)
    if isfield(data_struct, 'freqrange') && ~isempty(data_struct.freqrange)
        freq = data_struct.freqrange(:)';
    elseif isfield(dres, 'freq') && ~isempty(dres.freq)
        freq = dres.freq(:)';
    else
        freq = 1:F;
    end
    if numel(freq) < F
        freq = [freq, (numel(freq)+1):F];
    end
    freq = freq(1:F);
end
function T = to_tensor_(X)
    if iscell(X), T = cat(3, X{:}); else, T = X; end
end
function T = sanitize_cross_tensor_(T, eps_floor)
    [N, ~, F] = size(T);
    for f = 1:F
        Sf = T(:, :, f);
        Sf = (Sf + Sf') / 2;
        Sf(1:N+1:end) = max(real(diag(Sf)), eps_floor);
        T(:, :, f) = Sf;
    end
end
function d = diag3d_(T)
    [N, ~, F] = size(T);
    d = zeros(N, F);
    for f = 1:F, d(:, f) = diag(T(:, :, f)); end
end
function Coh = coherence_from_cross_(S, eps_floor)
    [N, ~, F] = size(S);
    Coh = zeros(N, N, F);
    for f = 1:F
        Sf = S(:, :, f);
        p = max(real(diag(Sf)), eps_floor);
        C = (abs(Sf).^2) ./ (p * p.');
        C(1:N+1:end) = 1;
        Coh(:, :, f) = min(max(C, 0), 1);
    end
end
function P = pcor_from_precision_tensor_(Omega, eps_floor)
    [N, ~, F] = size(Omega);
    P = zeros(N, N, F);
    for f = 1:F
        Om = (Omega(:,:,f) + Omega(:,:,f)')/2;
        d = max(real(diag(Om)), eps_floor);
        Pf = abs(-Om ./ sqrt(d * d.'));
        Pf(1:N+1:end) = 0;
        P(:,:,f) = Pf;
    end
end
function stats = summarize_upper_triangle_(A)
    n = size(A,1);
    v = A(triu(true(n),1));
    v = v(isfinite(v));
    if isempty(v)
        stats = struct('mean', NaN, 'p95', NaN, 'p99', NaN, 'max', NaN, 'nnz', 0);
    else
        stats = struct('mean', mean(v), 'p95', prctile(v,95), 'p99', prctile(v,99), 'max', max(v), 'nnz', nnz(v > 0));
    end
end
function clim_hi = robust_upper_clim_(A, sparse_mode)
    n = size(A,1);
    v = A(triu(true(n),1));
    v = v(isfinite(v));
    if sparse_mode
        v = v(v > 1e-10);
    end
    if isempty(v)
        clim_hi = 1e-6;
    else
        clim_hi = max(prctile(v, 99), 1e-6);
    end
end
function q = percentile_cols_(X, p)
    [~, M] = size(X); q = zeros(1, M);
    for m = 1:M, q(m) = percentile_(X(:, m), p); end
end
function q = percentile_(X, p)
    x = sort(X(:)); x = x(isfinite(x));
    if isempty(x), q = NaN; return; end
    idx = 1 + (numel(x)-1)*(p/100);
    lo = floor(idx); hi = ceil(idx);
    if lo == hi
        q = x(lo);
    else
        q = x(lo) + (idx-lo)*(x(hi)-x(lo));
    end
end
function v = get_opt_(opts, key, default_v)
    if isfield(opts, key), v = opts.(key); else, v = default_v; end
end