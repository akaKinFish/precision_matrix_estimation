function out = visualize_subject_jspace_result(result_mat_path, scalp_mat_path, opts)
%VISUALIZE_SUBJECT_JSPACE_RESULT Visualize one subject on scalp/source layers.
%   Now includes Alpha Peak Detection and Partial Coherence (PCOR) from Omega and Sjj.

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

dres = load(result_mat_path);
if ~isfield(dres, 'Sjj_est')
    error('Result file missing Sjj_est: %s', result_mat_path);
end

Sjj = to_tensor_(dres.Sjj_est);
Sjj = sanitize_cross_tensor_(Sjj, eps_floor);
[Nr, ~, Fsrc] = size(Sjj);

% Load Omega_est for Partial Coherence
has_omega = isfield(dres, 'Omega_est');
if has_omega
    Omega = to_tensor_(dres.Omega_est);
    Omega = sanitize_cross_tensor_(Omega, eps_floor);
else
    warning('Omega_est not found in result mat. PCOR will not be computed.');
end

scalp_file = resolve_scalp_path_(dres, scalp_mat_path);
dscalp = load(scalp_file, 'data_struct');
Svv = sanitize_cross_tensor_(dscalp.data_struct.CrossM, eps_floor);
[Nc, ~, Fscalp] = size(Svv);

freq = get_freq_(dscalp.data_struct, dres, Fscalp);
alpha_idx = find(freq >= alpha_band(1) & freq <= alpha_band(2));
if isempty(alpha_idx), alpha_idx = 1:numel(freq); end

Fuse = min([Fsrc, Fscalp, numel(freq)]);
Sjj = Sjj(:, :, 1:Fuse); Svv = Svv(:, :, 1:Fuse);
if has_omega, Omega = Omega(:, :, 1:Fuse); end
freq = freq(1:Fuse); alpha_idx = intersect(alpha_idx, 1:Fuse);

scalp_power = real(reshape(diag3d_(Svv), Nc, Fuse));
source_power = real(reshape(diag3d_(Sjj), Nr, Fuse));
scalp_power_db = 10 * log10(max(scalp_power, eps_floor));
source_power_db = 10 * log10(max(source_power, eps_floor));
scalp_alpha_cross = mean(Svv(:, :, alpha_idx), 3);

% --- PEAK DETECTION IN ALPHA BAND ---
mean_src_pow = mean(source_power_db, 1);
[~, peak_local_idx] = max(mean_src_pow(alpha_idx));
peak_f_idx = alpha_idx(peak_local_idx);
peak_freq = freq(peak_f_idx);
fprintf('\n[Peak Detection] Found Alpha peak at %.2f Hz (Index: %d)\n', peak_freq, peak_f_idx);

% --- COHERENCE CALCULATIONS ---
% 1. Average MSC
scalp_coh_msc = coherence_from_cross_(Svv, eps_floor);
scalp_alpha_msc = mean(scalp_coh_msc(:, :, alpha_idx), 3);
source_coh_msc = coherence_from_cross_(Sjj, eps_floor);

% 2. Peak Coherence (Magnitude)
source_mag_peak = sqrt(source_coh_msc(:, :, peak_f_idx));
source_mag_peak(1:Nr+1:end) = 0; % Remove diag

% 3. Peak Partial Coherence (PCOR from Omega)
source_pcor_peak = zeros(Nr, Nr);
if has_omega
    Om_peak = Omega(:, :, peak_f_idx);
    d_om = max(real(diag(Om_peak)), eps_floor);
    denom_om = sqrt(d_om * d_om.');
    P_peak = -Om_peak ./ denom_om;
    P_peak(1:Nr+1:end) = 0;
    source_pcor_peak = abs(P_peak); % Magnitude of Partial Coherence
end

% 4. Peak Partial Coherence directly from Sjj
source_pcor_from_s_peak = zeros(Nr, Nr);
try
    Sf = Sjj(:, :, peak_f_idx);
    Sf = (Sf + Sf') / 2;
    Sf_spd = sanitize_cross_tensor_(Sf, eps_floor); 
    Of_from_s = inv(Sf_spd);
    Of_from_s = (Of_from_s + Of_from_s')/2;
    
    d_of = max(real(diag(Of_from_s)), eps_floor);
    denom_of = sqrt(d_of * d_of.');
    P_diag = -Of_from_s ./ denom_of;
    P_diag(1:Nr+1:end) = 0; % Remove diag
    source_pcor_from_s_peak = abs(P_diag);
catch
    warning('Inversion of Sjj failed at peak freq. Cannot compute PCOR-from-S.');
end
% ------------------------------

if isempty(roi_order_in)
    off = source_mag_peak;
    [~, roi_order] = sort(sum(off, 2), 'descend');
else
    roi_order = roi_order_in(:);
end

figs = struct();

% 1) Scalp & Source spectrum (Combined for cleaner view)
figs.spectrum = figure('Name', 'Power Spectrum & Peak Detection', 'Position', [100 100 1000 400]);
subplot(1,2,1);
plot(freq, scalp_power_db.', 'Color', [0.8 0.8 0.8]); hold on;
plot(freq, mean(scalp_power_db, 1), 'k', 'LineWidth', 2);
ylm = ylim; patch([alpha_band(1) alpha_band(2) alpha_band(2) alpha_band(1)], ...
    [ylm(1) ylm(1) ylm(2) ylm(2)], [0.9 0.95 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
plot([peak_freq peak_freq], ylm, 'r--', 'LineWidth', 1.5);
xlabel('Frequency (Hz)'); ylabel('Power (dB)'); title('Scalp power spectrum'); grid on;

subplot(1,2,2);
med_src = median(source_power_db, 1);
fill([freq fliplr(freq)], [percentile_cols_(source_power_db, 25) fliplr(percentile_cols_(source_power_db, 75))], ...
    [0.85 0.9 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.5); hold on;
plot(freq, med_src, 'b', 'LineWidth', 2);
ylm = ylim; patch([alpha_band(1) alpha_band(2) alpha_band(2) alpha_band(1)], ...
    [ylm(1) ylm(1) ylm(2) ylm(2)], [0.9 0.95 1.0], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
plot([peak_freq peak_freq], ylm, 'r--', 'LineWidth', 1.5);
plot(peak_freq, med_src(peak_f_idx), 'rp', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
xlabel('Frequency (Hz)'); ylabel('Power (dB)'); title(sprintf('Source spectrum (Peak = %.1f Hz)', peak_freq)); grid on;

% ============================================================
% 2) Matrices at Peak (Updated for dynamic 1x3 subplots)
% ============================================================
num_mats = 1;
if has_omega, num_mats = num_mats + 1; end
if any(source_pcor_from_s_peak(:) > 0), num_mats = num_mats + 1; end

figs.peak_matrices = figure('Name', 'Peak Source Matrices', 'Position', [150 150 450*num_mats 400]);
curr_sub = 1;

% [1] Linear Coherence
subplot(1, num_mats, curr_sub);
Csrc = source_mag_peak(roi_order, roi_order); Csrc(1:Nr+1:end) = NaN;
imagesc(Csrc, [0 max(percentile_(Csrc(isfinite(Csrc)), 99), 1e-6)]); axis image; colorbar;
title(sprintf('Linear Coherence at %.1f Hz', peak_freq));
curr_sub = curr_sub + 1;

% [2] JSPACE PCOR (from Omega)
if has_omega
    subplot(1, num_mats, curr_sub);
    Psrc = source_pcor_peak(roi_order, roi_order);
    Psrc(1:Nr+1:end) = NaN;
    maskU = triu(true(Nr), 1);
    v = Psrc(maskU);
    v = v(isfinite(v) & v > 1e-10);
    if isempty(v)
        clim_hi = 1e-6;
    else
        clim_hi = max(prctile(v, 99), 1e-6);
    end
    imagesc(Psrc, [0 clim_hi]); axis image; colorbar;
    title(sprintf('JSPACE PCOR at %.1f Hz', peak_freq));
    curr_sub = curr_sub + 1;
end

% [3] Native PCOR (from Sjj)
if any(source_pcor_from_s_peak(:) > 0)
    subplot(1, num_mats, curr_sub);
    Psrc_s = source_pcor_from_s_peak(roi_order, roi_order);
    Psrc_s(1:Nr+1:end) = NaN;
    maskU = triu(true(Nr), 1);
    v_s = Psrc_s(maskU);
    v_s = v_s(isfinite(v_s) & v_s > 1e-10);
    % if isempty(v_s)
    %     clim_s = 1e-6;
    % else
    %     clim_s = max(prctile(v_s, 99), 1e-6);
    % end
    imagesc(Psrc_s, [0 max(v_s)]); axis image; colorbar;
    title(sprintf('Native PCOR (from S^{-1}) at %.1f Hz', peak_freq));
end


% --- 新增：绘制 ROI 内所有顶点的相干性频谱图 ---
ai = 8;
if ai > numel(SC.Atlas), ai = 1; end
current_scouts = SC.Atlas(ai).Scouts;

seed_idx = [];
occ_keywords = {'Occipital', 'V1', 'O1', 'V2', 'cuneus', 'calcarine'};
for i = 1:numel(current_scouts)
    lbl = current_scouts(i).Label;
    if any(contains(lbl, occ_keywords, 'IgnoreCase', true)) && ...
            (contains(lbl, 'L') || contains(lbl, 'left') || startsWith(lbl, 'L_'))
        seed_idx = i;
        seed_label = lbl;
        break;
    end
end
if isempty(seed_idx)
    seed_idx = 1;
    seed_label = current_scouts(seed_idx).Label;
end

roi_vertices = current_scouts(seed_idx).Vertices;
roi_vertices = roi_vertices(roi_vertices > 0 & roi_vertices <= Nr);
num_vertices = numel(roi_vertices);

if num_vertices == 0
    warning('未找到有效的顶点，跳过绘制 ROI 相干性频谱图。');
else
    roi_coh_tensor = source_coh_msc(roi_vertices, :, :);
    vertex_mean_coh = squeeze(mean(roi_coh_tensor, 2));
    if size(vertex_mean_coh, 2) ~= Fuse
        vertex_mean_coh = vertex_mean_coh.';
    end
    
    figs.coherence_spec = figure('Name', 'ROI Vertices Coherence Spectrum', 'Position', [200 200 700 450]);
    hold on;
    ylm_max = max(vertex_mean_coh(:)) * 1.1;
    patch([alpha_band(1) alpha_band(2) alpha_band(2) alpha_band(1)], ...
        [0 0 ylm_max ylm_max], [1 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot([peak_freq peak_freq], [0 ylm_max], 'k--', 'LineWidth', 1.5);
    
    colors = lines(num_vertices);
    for v = 1:num_vertices
        plot(freq, vertex_mean_coh(v, :), 'Color', colors(v, :), 'LineWidth', 1.5);
    end
    ylim([0, ylm_max]); xlim([freq(1), freq(end)]);
    xlabel('Frequency (Hz)'); ylabel('Mean Global Coherence');
    title(sprintf('Coherence Spectrum (ROI: %s, %d Vertices)', seed_label, num_vertices), 'Interpreter', 'none');
    grid on;
end

% --- 新增：绘制全脑所有非冗余连边 (Upper Triangular) 的相干性频谱图 ---
maskU = triu(true(Nr), 1);
num_edges = nnz(maskU); 
all_coh_lines = zeros(num_edges, Fuse);
for f = 1:Fuse
    temp_mat = source_coh_msc(:, :, f);
    all_coh_lines(:, f) = temp_mat(maskU);
end

figs.global_coherence_spec = figure('Name', 'Global Coherence (Upper Triangular)', 'Position', [250 250 700 450]);
hold on;
ylm_max = 1.0; 
patch([alpha_band(1) alpha_band(2) alpha_band(2) alpha_band(1)], ...
    [0 0 ylm_max ylm_max], [1 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
plot([peak_freq peak_freq], [0 ylm_max], 'k--', 'LineWidth', 1.5);

if num_edges > 20000
    fprintf('全脑连边数量较多 (%d 条)，正在渲染叠加图，请稍候...\n', num_edges);
    plot(freq, all_coh_lines.', 'Color', [0.6 0.6 0.6 0.01]); 
else
    plot(freq, all_coh_lines.', 'Color', [0.7 0.7 0.7 0.05]); 
end

global_avg_coh = mean(all_coh_lines, 1);
plot(freq, global_avg_coh, 'r', 'LineWidth', 2.5);
ylim([0, 1]); xlim([freq(1), freq(end)]);
xlabel('Frequency (Hz)'); ylabel('MSC Coherence');
title(sprintf('Global Coherence Spectrum (Superimposed %d Unique Edges)', num_edges));
grid on;

% 3) Region Aggregation, Chords & 3D Network
out_region_coh = struct();
out_region_pcor = struct();
if ~isempty(SC)
    region_opts = struct('atlas_name', atlas_name, 'edge_keep_quantile', edge_keep_quantile, ...
        'edge_keep_max', edge_keep_max, 'chord_mode', 'both');
    
    % --- A. COHERENCE PLOTS ---
    region_opts.fig_prefix = 'PEAK COHERENCE';
    region_out_coh = visualize_source_coherence_brain_regions(source_mag_peak, SC, region_opts);
    figs.source_region_chord_coh = region_out_coh.figs.region_chord;
    figs.source_roi_chord_coh = region_out_coh.figs.roi_chord;
    out_region_coh = region_out_coh;
    
    figs.source_3d_network_coh = figure('Name', 'Source 3D Brain Network (Coherence)');
    plot_3d_brain_network_(SC, source_mag_peak, atlas_name, roi_to_scout_map, ...
        head_surface_view, edge_keep_quantile, edge_keep_max, 1e-6, ...
        sprintf('Source Coherence Network at %.1f Hz', peak_freq));
        
    % --- B. PARTIAL COHERENCE PLOTS (from JSPACE Omega) ---
    if has_omega
        region_opts.fig_prefix = 'PEAK PCOR';
        region_out_pcor = visualize_source_coherence_brain_regions(source_pcor_peak, SC, region_opts);
        figs.source_region_chord_pcor = region_out_pcor.figs.region_chord;
        figs.source_roi_chord_pcor = region_out_pcor.figs.roi_chord;
        out_region_pcor = region_out_pcor;
        
        figs.source_3d_network_pcor = figure('Name', 'Source 3D Brain Network (PCOR)');
        plot_3d_brain_network_(SC, source_pcor_peak, atlas_name, roi_to_scout_map, ...
            head_surface_view, edge_keep_quantile, edge_keep_max, 1e-6, ...
            sprintf('Source Partial Coherence (Direct Paths) at %.1f Hz', peak_freq));
    end
    
    % --- C. PCOR FROM Sjj PLOTS ---
    if any(source_pcor_from_s_peak(:) > 0)
        figs.source_3d_network_pcor_from_s = figure('Name', 'Source 3D Brain Network (PCOR from Sjj)');
        plot_3d_brain_network_(SC, source_pcor_from_s_peak, atlas_name, roi_to_scout_map, ...
            head_surface_view, edge_keep_quantile, edge_keep_max, 1e-6, ...
            sprintf('Source PCOR (from S_{jj}^{-1}) at %.1f Hz', peak_freq));
    end
end

% ============================================================
% --- [诊断模块]：基于 Omega 逆矩阵的分析 (不保存) ---
% ============================================================
if has_omega
    fprintf('\n[Diagnostics] Running Omega-inversion analysis...\n');
    Sjj_omega = zeros(Nr, Nr, Fuse);
    for f = 1:Fuse
        Of = Omega(:, :, f);
        Of_spd = sanitize_cross_tensor_(Of, eps_floor); 
        try
            Sf_inv = inv(Of_spd);
            Sf_inv = (Sf_inv + Sf_inv')/2; 
            Sjj_omega(:, :, f) = Sf_inv;
        catch
            Sjj_omega(:, :, f) = Of; 
            warning('Inversion failed at freq idx %d', f);
        end
    end
    
    source_power_om = real(reshape(diag3d_(Sjj_omega), Nr, Fuse));
    source_power_om_db = 10 * log10(max(source_power_om, eps_floor));
    
    figs.diag_spectrum_om = figure('Name', 'Implied Spectrum from Omega^-1', 'Position', [150 150 600 400]);
    hold on;
    fill([freq fliplr(freq)], [percentile_cols_(source_power_om_db, 25) fliplr(percentile_cols_(source_power_om_db, 75))], ...
        [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5); 
    
    plot(freq, mean(source_power_om_db, 1), 'k', 'LineWidth', 2);
    plot(freq, mean(source_power_db, 1), 'b--', 'LineWidth', 1.5); 
    
    ylm = ylim; patch([alpha_band(1) alpha_band(2) alpha_band(2) alpha_band(1)], ...
        [ylm(1) ylm(1) ylm(2) ylm(2)], [1 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    plot([peak_freq peak_freq], ylm, 'r--', 'LineWidth', 1.5);
    
    xlabel('Frequency (Hz)'); ylabel('Power (dB)'); 
    title('Implied Source Spectrum (\Omega^{-1} vs Original S)');
    legend({'IQR(\Omega^{-1})', 'Mean(\Omega^{-1})', 'Mean(Orig S)'}, 'Location', 'best');
    grid on;
    
    coh_msc_om = coherence_from_cross_(Sjj_omega, eps_floor);
    coh_lines_om = zeros(num_edges, Fuse);
    for f = 1:Fuse
        temp_mat = coh_msc_om(:, :, f);
        coh_lines_om(:, f) = temp_mat(maskU);
    end
    
    figs.diag_global_coh_om = figure('Name', 'Global Coherence from Omega^-1', 'Position', [200 200 700 450]);
    hold on;
    ylm_max = 1.0; 
    patch([alpha_band(1) alpha_band(2) alpha_band(2) alpha_band(1)], ...
        [0 0 ylm_max ylm_max], [1 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
    plot([peak_freq peak_freq], [0 ylm_max], 'k--', 'LineWidth', 1.5);
    
    plot(freq, coh_lines_om.', 'Color', [0.2 0.8 0.2 0.01]); 
    plot(freq, mean(coh_lines_om, 1), 'k', 'LineWidth', 2.5);
    plot(freq, global_avg_coh, 'r--', 'LineWidth', 2); 
    
    ylim([0, 1]); xlim([freq(1), freq(end)]);
    xlabel('Frequency (Hz)'); ylabel('MSC Coherence');
    title(sprintf('Implied Global Coherence from \\Omega^{-1} (%d Edges)', num_edges));
    legend({'Alpha Band', 'Peak', '\Omega^{-1} Edges', 'Mean(\Omega^{-1})', 'Mean(Orig S)'}, 'Location', 'bestoutside');
    grid on;
end

% ============================================================
% Save Figures
if ~isempty(save_dir)
    if ~isfolder(save_dir), mkdir(save_dir); end
    saveas(figs.spectrum, fullfile(save_dir, 'fig1_power_spectrum_with_peak.png'));
    saveas(figs.peak_matrices, fullfile(save_dir, 'fig2_peak_matrices.png'));
    if isfield(figs, 'source_3d_network_coh'), saveas(figs.source_3d_network_coh, fullfile(save_dir, 'fig3a_3d_network_coh.png')); end
    if isfield(figs, 'source_region_chord_coh'), saveas(figs.source_region_chord_coh, fullfile(save_dir, 'fig3b_region_chord_coh.png')); end
    if isfield(figs, 'source_roi_chord_coh'), saveas(figs.source_roi_chord_coh, fullfile(save_dir, 'fig3c_roi_chord_coh.png')); end
    if has_omega
        if isfield(figs, 'source_3d_network_pcor'), saveas(figs.source_3d_network_pcor, fullfile(save_dir, 'fig4a_3d_network_pcor.png')); end
        if isfield(figs, 'source_region_chord_pcor'), saveas(figs.source_region_chord_pcor, fullfile(save_dir, 'fig4b_region_chord_pcor.png')); end
        if isfield(figs, 'source_roi_chord_pcor'), saveas(figs.source_roi_chord_pcor, fullfile(save_dir, 'fig4c_roi_chord_pcor.png')); end
    end
    if isfield(figs, 'source_3d_network_pcor_from_s'), saveas(figs.source_3d_network_pcor_from_s, fullfile(save_dir, 'fig5_3d_network_pcor_s.png')); end
end

out = struct('figs', figs, 'peak_freq', peak_freq, 'peak_f_idx', peak_f_idx, ...
    'source_mag_peak', source_mag_peak, 'source_pcor_peak', source_pcor_peak, ...
    'source_pcor_from_s_peak', source_pcor_from_s_peak, ...
    'roi_order', roi_order, 'region_coh', out_region_coh, 'region_pcor', out_region_pcor);
end

% --- HELPER FUNCTIONS ---
function plot_3d_brain_network_(SC, C, atlas_name, roi_to_scout_map, view_ang, q_keep, max_edges, min_edge_abs, plt_title)
if ~isfield(SC, 'Vertices') || ~isfield(SC, 'Faces')
    text(0.5, 0.5, 'SC has no surface mesh.', 'HorizontalAlignment', 'center'); axis off; return;
end
V = SC.Vertices; N = size(C, 1);
ai = 1;
if isfield(SC, 'Atlas') && ~isempty(SC.Atlas)
    for k = 1:numel(SC.Atlas), if strcmpi(char(SC.Atlas(k).Name), atlas_name), ai = k; break; end; end
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
allw = C(triu(true(N), 1)); allw = allw(isfinite(allw) & allw > min_edge_abs);
if isempty(allw), thr = inf; else, thr = percentile_(allw, q_keep * 100); end
edges = [];
for i = 1:N
    for j = i+1:N
        if valid_node(i) && valid_node(j)
            w = C(i,j); if isfinite(w) && w >= thr && w > min_edge_abs, edges = [edges; i, j, w]; end
        end
    end
end
if ~isempty(max_edges) && size(edges,1) > max_edges
    [~, ord] = sort(edges(:,3), 'descend'); edges = edges(ord(1:max_edges), :);
end
patch('Vertices', V, 'Faces', SC.Faces, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.15);
hold on; axis equal off vis3d;
if isempty(edges), title([plt_title ' (No edges)']); return; end
cmap = colormap(gca, 'autumn'); 
wmin = min(edges(:,3)); wmax = max(edges(:,3));
for e = 1:size(edges,1)
    i = edges(e,1); j = edges(e,2); w = edges(e,3);
    if wmax > wmin, t = (w - wmin)/(wmax - wmin); else, t = 0.5; end
    col = cmap(max(1, min(size(cmap,1), round(t * size(cmap,1)))), :);
    p1 = centers(i,:); p2 = centers(j,:);
    plot3([p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)], '-', 'Color', [col 0.8], 'LineWidth', 1 + 5 * t);
end
active_nodes = unique(edges(:, 1:2));
plot3(centers(active_nodes, 1), centers(active_nodes, 2), centers(active_nodes, 3), ...
    'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 5, 'LineWidth', 1);
view(view_ang(1), view_ang(2)); camlight('headlight'); lighting gouraud; title(plt_title); hold off;
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
if isfield(data_struct, 'freqrange') && ~isempty(data_struct.freqrange), freq = data_struct.freqrange(:)';
elseif isfield(dres, 'freq') && ~isempty(dres.freq), freq = dres.freq(:)';
else, freq = 1:F; end
if numel(freq) < F, freq = [freq, (numel(freq)+1):F]; end
end

function T = to_tensor_(X), if iscell(X), T = cat(3, X{:}); else, T = X; end; end

function T = sanitize_cross_tensor_(T, eps_floor)
[N, ~, F] = size(T);
for f = 1:F, Sf = T(:, :, f); Sf = (Sf + Sf') / 2; Sf(1:N+1:end) = max(real(diag(Sf)), eps_floor); T(:, :, f) = Sf; end
end

function d = diag3d_(T)
[N, ~, F] = size(T); d = zeros(N, F);
for f = 1:F, d(:, f) = diag(T(:, :, f)); end
end

function Coh = coherence_from_cross_(S, eps_floor)
[N, ~, F] = size(S); Coh = zeros(N, N, F);
for f = 1:F, Sf = S(:, :, f); p = max(real(diag(Sf)), eps_floor);
    C = (abs(Sf).^2) ./ (p * p.'); C(1:N+1:end) = 1; Coh(:, :, f) = min(max(C, 0), 1); end
end

function q = percentile_cols_(X, p)
[~, M] = size(X); q = zeros(1, M);
for m = 1:M, q(m) = percentile_(X(:, m), p); end
end

function q = percentile_(X, p)
x = sort(X(:)); x = x(isfinite(x)); if isempty(x), q = NaN; return; end
idx = 1 + (numel(x)-1)*(p/100); lo = floor(idx); hi = ceil(idx);
if lo == hi, q = x(lo); else, q = x(lo) + (idx-lo)*(x(hi)-x(lo)); end
end

function v = get_opt_(opts, key, default_v), if isfield(opts, key), v = opts.(key); else, v = default_v; end; end