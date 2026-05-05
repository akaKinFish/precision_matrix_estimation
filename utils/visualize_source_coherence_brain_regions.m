function out = visualize_source_coherence_brain_regions(source_alpha_mag, SC, opts)
%VISUALIZE_SOURCE_COHERENCE_BRAIN_REGIONS Plot atlas-aware source coherence views.
%   Now supports fig_prefix to distinguish between Coherence and PCOR plots.

    if nargin < 3 || isempty(opts), opts = struct(); end
    
    atlas_name = get_opt_(opts, 'atlas_name', 'HCP_MMP1');
    hide_diagonal = get_opt_(opts, 'hide_diagonal', true);
    exclude_unassigned = get_opt_(opts, 'exclude_unassigned', true);
    
    roi_q   = get_opt_(opts, 'roi_edge_keep_quantile', 0.998);
    roi_max = get_opt_(opts, 'roi_edge_keep_max', 150);
    reg_q   = get_opt_(opts, 'region_edge_keep_quantile', 0.92);
    reg_max = get_opt_(opts, 'region_edge_keep_max', 15);
    
    region_stat = lower(strtrim(get_opt_(opts, 'region_stat', 'max'))); 
    chord_show_labels = get_opt_(opts, 'chord_show_labels', false);
    chord_mode = lower(strtrim(get_opt_(opts, 'chord_mode', 'both')));
    min_edge_abs = get_opt_(opts, 'min_edge_abs', 1e-6); 
    debug_print = get_opt_(opts, 'debug_print', true);
    debug_topK  = get_opt_(opts, 'debug_topK', 10);
    
    fig_prefix = get_opt_(opts, 'fig_prefix', '');

    if ~isempty(fig_prefix), fig_prefix = [fig_prefix ' - ']; end
    
    [roi_order, roi_labels, roi_region_id, region_labels, group_counts] = ...
        atlas_roi_order_(SC, size(source_alpha_mag, 1), atlas_name);
    region_labels_full = cellfun(@region_full_name_, region_labels, 'UniformOutput', false);
    
    Croi = source_alpha_mag(roi_order, roi_order);
    rid  = roi_region_id(roi_order);
    roi_labels_ord = roi_labels(roi_order);
    
    if exclude_unassigned
        keep_region = ~strcmp(region_labels, 'UNASSIGNED');
        if any(~keep_region)
            keep_roi = keep_region(rid);
            Croi = Croi(keep_roi, keep_roi); rid = rid(keep_roi);
            roi_labels_ord = roi_labels_ord(keep_roi);
            old2new = zeros(numel(region_labels), 1); old2new(keep_region) = 1:nnz(keep_region);
            rid = old2new(rid);
            region_labels = region_labels(keep_region); region_labels_full = region_labels_full(keep_region);
            group_counts = accumarray(rid, 1, [numel(region_labels), 1]);
        end
    end
    
    if hide_diagonal, Croi(1:size(Croi,1)+1:end) = NaN; end
    K = numel(region_labels);
    
    [Creg_mean_nz, Creg_max, region_counts] = region_aggregation_matrix_(Croi, rid, K, min_edge_abs);
    
    switch region_stat
        case 'mean', Creg_show = Creg_mean_nz; reg_tag = 'mean(non-zero)';
        case 'max',  Creg_show = Creg_max;     reg_tag = 'max';
        case 'both', Creg_show = Creg_max;     reg_tag = 'max'; 
        otherwise,   Creg_show = Creg_max;     reg_tag = 'max';
    end
    
    if debug_print
        fprintf('\n[DBG] === %sRegion/ROI chord settings ===\n', fig_prefix);
        print_top_pairs_(Creg_mean_nz, region_labels_full, debug_topK, [fig_prefix 'Region matrix (MEAN of NON-ZERO)']);
        print_top_pairs_(Creg_max, region_labels_full, debug_topK, [fig_prefix 'Region matrix (MAX)']);
    end
    
    figs = struct();
    
    figs.region_heatmap = figure('Name', [fig_prefix 'Source matrix reordered by atlas region']);
    lim_hi = percentile_(Croi(isfinite(Croi)), 99);
    imagesc(Croi, [0 max(lim_hi, 1e-6)]); axis image; colorbar;
    title([fig_prefix 'Source Matrix reordered by region']); hold on;
    boundaries = cumsum(group_counts);
    for k = 1:numel(boundaries)-1
        x = boundaries(k) + 0.5;
        plot([x x], [0.5 size(Croi,1)+0.5], 'w-', 'LineWidth', 0.5);
        plot([0.5 size(Croi,1)+0.5], [x x], 'w-', 'LineWidth', 0.5);
    end
    centers = boundaries - group_counts/2;
    set(gca, 'XTick', centers, 'YTick', centers, 'XTickLabel', region_labels_full, ...
        'YTickLabel', region_labels_full, 'TickLabelInterpreter', 'none', 'FontSize', 8);
    xtickangle(45); hold off;
    
    if strcmp(chord_mode, 'region') || strcmp(chord_mode, 'both')
        figs.region_chord = figure('Name', [fig_prefix 'Region-level chord']);
        plot_ribbon_chord_(Creg_show, region_labels_full, (1:K)', region_labels_full, reg_q, reg_max, min_edge_abs, true);
        title(sprintf('%sRegion chord (%s, q=%.3f, max=%d)', fig_prefix, reg_tag, reg_q, reg_max));
    end
    
    if strcmp(chord_mode, 'roi') || strcmp(chord_mode, 'both')
        figs.roi_chord = figure('Name', [fig_prefix 'ROI-level chord']);
        plot_ribbon_chord_(Croi, roi_labels_ord, rid, region_labels_full, roi_q, roi_max, min_edge_abs, chord_show_labels);
        title(sprintf('%sROI chord (q=%.3f, max=%d)', fig_prefix, roi_q, roi_max));
    end
    
    out = struct('figs', figs, 'roi_order', roi_order, 'region_matrix_mean_nz', Creg_mean_nz, 'region_matrix_max', Creg_max);
end

% --- HELPER FUNCTIONS ---
function [Cmean_nz, Cmax, counts] = region_aggregation_matrix_(C, rid, K, min_edge_abs)
    Cmean_nz = nan(K,K); Cmax = nan(K,K); counts = zeros(K,1);
    for a = 1:K
        ia = find(rid == a); counts(a) = numel(ia);
        for b = 1:K
            ib = find(rid == b); block = C(ia, ib);
            if a == b && numel(ia) > 0, block(1:numel(ia)+1:end) = NaN; end
            v = block(isfinite(block));
            if isempty(v)
                Cmean_nz(a,b) = NaN; Cmax(a,b) = NaN;
            else
                v_nz = v(abs(v) > min_edge_abs); 
                if isempty(v_nz), Cmean_nz(a,b) = 0; else, Cmean_nz(a,b) = mean(v_nz); end
                Cmax(a,b) = max(abs(v));
            end
        end
    end
    Cmean_nz = (Cmean_nz + Cmean_nz')/2; Cmax = (Cmax + Cmax')/2;
end

function plot_ribbon_chord_(C, labels, region_id, region_labels, q_keep, edge_keep_max, min_edge_abs, show_labels)
    N = size(C, 1); theta = linspace(0, 2*pi, N + 1); theta(end) = []; dth = 2*pi/N;
    cmap_region = lines(max(region_id)); clf; hold on; axis equal off;
    allw = C(triu(true(N), 1)); allw = allw(isfinite(allw) & allw > min_edge_abs);
    if isempty(allw), thr = inf; else, thr = percentile_(allw, q_keep * 100); end
    edges = [];
    for i = 1:N
        for j = i+1:N
            w = C(i,j); if isfinite(w) && w >= thr && w > min_edge_abs, edges = [edges; i, j, w]; end
        end
    end
    if ~isempty(edge_keep_max) && size(edges,1) > edge_keep_max
        [~, ord] = sort(edges(:,3), 'descend'); edges = edges(ord(1:edge_keep_max), :);
    end
    if ~isempty(edges)
        wmin = min(edges(:,3)); wmax = max(edges(:,3));
        for e = 1:size(edges,1)
            i = edges(e,1); j = edges(e,2); w = edges(e,3);
            if wmax > wmin, t = (w - wmin)/(wmax - wmin); else, t = 0.5; end
            halfw = (0.06 + 0.32*t) * dth;
            p1 = [cos(theta(i) + halfw), sin(theta(i) + halfw)]; p2 = [cos(theta(i) - halfw), sin(theta(i) - halfw)];
            q1 = [cos(theta(j) + halfw), sin(theta(j) + halfw)]; q2 = [cos(theta(j) - halfw), sin(theta(j) - halfw)];
            ctrl_radius = 0.35; 
            cp1 = ctrl_radius * p1; cp2 = ctrl_radius * q2; cp3 = ctrl_radius * q1; cp4 = ctrl_radius * p2;
            B1 = bezier_curve_(p1, cp1, cp2, q2, 24); B2 = bezier_curve_(q1, cp3, cp4, p2, 24);
            poly = [B1; B2]; c = 0.5*(cmap_region(region_id(i),:) + cmap_region(region_id(j),:));
            patch(poly(:,1), poly(:,2), c, 'EdgeColor','none', 'FaceAlpha', 0.2 + 0.6*t);
        end
    end
    r1 = 0.95; r2 = 1.05;
    for i = 1:N
        th = linspace(theta(i)-0.45*dth, theta(i)+0.45*dth, 14); c = cmap_region(region_id(i),:);
        patch([r2*cos(th), fliplr(r1*cos(th))], [r2*sin(th), fliplr(r1*sin(th))], c, 'EdgeColor','none');
        if show_labels
            text(1.12*cos(theta(i)), 1.12*sin(theta(i)), labels{i}, 'FontSize',6, 'HorizontalAlignment','center', 'Interpreter','none');
        end
    end
    lh = gobjects(numel(region_labels),1);
    for k = 1:numel(region_labels), lh(k) = plot(nan,nan,'o','MarkerFaceColor',cmap_region(k,:),'MarkerEdgeColor',cmap_region(k,:)); end
    legend(lh, region_labels, 'Interpreter','none', 'Location','eastoutside'); hold off;
end

function B = bezier_curve_(p0, p1, p2, p3, n)
    t = linspace(0,1,n)'; B = ((1-t).^3)*p0 + 3*((1-t).^2).*t.*p1 + 3*(1-t).*(t.^2).*p2 + (t.^3).*p3;
end
function print_top_pairs_(C, labels, topK, tag)
    K = size(C,1); M = C; M(1:K+1:end) = NaN; lin = find(triu(true(K),1)); vals = M(lin);
    good = isfinite(vals); lin = lin(good); vals = vals(good);
    if isempty(vals), fprintf('[DBG] %s: no finite edges.\n', tag); return; end
    [vals_sorted, ord] = sort(vals, 'descend'); nshow = min(topK, numel(vals_sorted));
    fprintf('\n[DBG] Top %d pairs: %s\n', nshow, tag);
    for k = 1:nshow, [i,j] = ind2sub([K K], lin(ord(k))); fprintf('  %2d) %s -- %s : %.6g\n', k, labels{i}, labels{j}, vals_sorted(k)); end
end
function [roi_order, roi_labels, roi_region_id, region_labels, group_counts] = atlas_roi_order_(SC, Nr, atlas_name)
    scouts = [];
    if isfield(SC, 'Atlas')
        idx = 1; for i = 1:numel(SC.Atlas), if strcmpi(char(SC.Atlas(i).Name), atlas_name), idx = i; break; end; end
        if isfield(SC.Atlas(idx), 'Scouts'), scouts = SC.Atlas(idx).Scouts; end
    end
    if isempty(scouts)
        roi_order = (1:Nr)'; roi_labels = arrayfun(@(k) sprintf('ROI_%d',k), 1:Nr, 'UniformOutput', false)';
        region_labels = {'ALL'}; roi_region_id = ones(Nr, 1); group_counts = Nr; return;
    end
    nScout = min(numel(scouts), Nr); roi_labels = cell(Nr, 1); region_name = cell(Nr, 1);
    for i = 1:nScout
        if isfield(scouts(i), 'Label') && ~isempty(scouts(i).Label), roi_labels{i} = char(scouts(i).Label); else, roi_labels{i} = sprintf('ROI_%d', i); end
        if isfield(scouts(i), 'Region') && ~isempty(scouts(i).Region), region_name{i} = strtrim(char(scouts(i).Region)); else, region_name{i} = 'UNASSIGNED'; end
    end
    for i = nScout+1:Nr, roi_labels{i} = sprintf('ROI_%d', i); region_name{i} = 'UNASSIGNED'; end
    [region_labels, ~, roi_region_id] = unique(region_name, 'stable');
    [~, roi_order] = sortrows([roi_region_id(:), (1:Nr)']);
    K = numel(region_labels); group_counts = zeros(K, 1);
    for k = 1:K, group_counts(k) = sum(roi_region_id == k); end
end
function full = region_full_name_(code)
    c = upper(strtrim(char(code)));
    switch c
        case 'LF', full = 'Left Frontal (LF)'; case 'LPF', full = 'Left Prefrontal (LPF)';
        case 'LC', full = 'Left Central (LC)'; case 'LL', full = 'Left Limbic (LL)';
        case 'LP', full = 'Left Parietal (LP)'; case 'LT', full = 'Left Temporal (LT)';
        case 'LO', full = 'Left Occipital (LO)'; case 'RF', full = 'Right Frontal (RF)';
        case 'RPF', full = 'Right Prefrontal (RPF)'; case 'RC', full = 'Right Central (RC)';
        case 'RL', full = 'Right Limbic (RL)'; case 'RP', full = 'Right Parietal (RP)';
        case 'RT', full = 'Right Temporal (RT)'; case 'RO', full = 'Right Occipital (RO)';
        case 'LU', full = 'Left Hippocampus (L_Hipp)'; case 'RU', full = 'Right Hippocampus (R_Hipp)';
        case 'UNASSIGNED', full = 'Unassigned (UNASSIGNED)'; otherwise, full = sprintf('%s (%s)', code, code);
    end
end
function v = get_opt_(opts, key, default_v), if isfield(opts, key), v = opts.(key); else, v = default_v; end; end
function q = percentile_(X, p)
    x = sort(X(:)); x = x(isfinite(x)); if isempty(x), q = NaN; return; end
    idx = 1 + (numel(x)-1)*(p/100); lo = floor(idx); hi = ceil(idx);
    if lo == hi, q = x(lo); else, q = x(lo) + (idx-lo)*(x(hi)-x(lo)); end
end