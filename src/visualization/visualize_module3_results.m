function visualize_module3_results(edge_proxies, combined_masks, active_set_stats, tau, varargin)
% VISUALIZE_MODULE3_RESULTS
% Compact dashboard for Module 3 deactivate-only active-set results.
%
% Usage:
%   visualize_module3_results(edge_proxies, combined_masks, active_set_stats, tau)
%   visualize_module3_results(..., 'frequencies', freq_vec, 'freq_index', k, 'title', txt)
%
% Inputs:
%   edge_proxies     : cell {F x 1} of p x p proxy matrices (real or complex)
%   combined_masks   : struct with fields:
%                       - edge_masks (p x p x F)
%                       - node_masks (p x F)
%                       - combined_masks (p x p x F)
%                       - deactivated_masks (p x p x F) [optional]
%   active_set_stats : struct (as emitted by module3_combined_active_set)
%   tau              : numeric threshold used
%
% Name-Value:
%   'frequencies'    : vector of F display frequencies (default: 1:F)
%   'freq_index'     : which frequency to zoom in for heatmaps (default: 1)
%   'title'          : figure suptitle string (default: 'Module 3 Active-Set Visualization')

    p = inputParser;
    addParameter(p, 'frequencies', [], @(x) isempty(x) || isvector(x));
    addParameter(p, 'freq_index', 1, @(x) isscalar(x) && x >= 1);
    addParameter(p, 'title', 'Module 3 Active-Set Visualization', @ischar);
    parse(p, varargin{:});
    opts = p.Results;

    % Basic dims and safety
    F = numel(edge_proxies);
    if isempty(opts.frequencies), opts.frequencies = 1:F; end
    freq_index = min(max(1, opts.freq_index), F);

    % Prepare counts
    E_keep = combined_masks.combined_masks;         % kept after node gating
    if isfield(combined_masks, 'deactivated_masks') && ~isempty(combined_masks.deactivated_masks)
        E_off = combined_masks.deactivated_masks;
    else
        % Fallback: anything not kept is considered deactivated (off-diagonal only)
        E_off = false(size(E_keep));
        for f = 1:F
            off = ~E_keep(:,:,f) & ~eye(size(E_keep,1));
            E_off(:,:,f) = off | off';
        end
    end

    kept_counts = zeros(F,1);
    off_counts  = zeros(F,1);
    for f = 1:F
        kept_off = E_keep(:,:,f) & ~eye(size(E_keep,1));
        off_off  = E_off(:,:,f)  & ~eye(size(E_off,1));
        kept_counts(f) = sum(kept_off(:))/2;
        off_counts(f)  = sum(off_off(:))/2;
    end

    % Pull stats if available
    if isfield(active_set_stats,'sparsity_ratios_per_frequency')
        sr = active_set_stats.sparsity_ratios_per_frequency(:);
    else
        % compute sparsity from kept_counts
        pN = size(E_keep,1);
        total_possible = pN*(pN-1)/2;
        sr = kept_counts / max(total_possible,1);
    end

    % Choose a frequency to show
    C = edge_proxies{freq_index};
    if ~isreal(C), Cshow = abs(C); else, Cshow = C; end
    K = E_keep(:,:,freq_index);
    OFF = E_off(:,:,freq_index);
    Nmask = combined_masks.node_masks(:,freq_index);

    % Histogram of proxy values (upper triangle only), with tau line
    ut_mask = triu(true(size(Cshow)),1);
    proxy_vals = Cshow(ut_mask);
    proxy_vals = proxy_vals(isfinite(proxy_vals));

    % ---- Figure layout ----
    figure('Name','Module 3 Active-Set Visualization','Color','w','Position',[100,60,1300,820]);

    % (1) Kept vs Deactivated per frequency
    subplot(2,3,1);
    b = bar([kept_counts, off_counts], 'stacked');
    title('Kept vs Deactivated Edges per Frequency');
    xlabel('Frequency Index'); ylabel('# Edges (off-diagonal)');
    legend({'Kept','Deactivated'}, 'Location','best');
    xticks(1:F); xticklabels(string(opts.frequencies));
    grid on;

    % (2) Sparsity (density) per frequency
    subplot(2,3,2);
    plot(1:F, sr, 'o-','LineWidth',1.8);
    title('Kept Edge Density per Frequency');
    xlabel('Frequency Index'); ylabel('Density (kept / possible)');
    xticks(1:F); xticklabels(string(opts.frequencies));
    ylim([0, min(1, max(sr)*1.1 + 0.02)]);
    grid on;

    % (3) Node activity mask across frequencies
    subplot(2,3,3);
    imagesc(double(combined_masks.node_masks));
    colormap(gca, gray);
    colorbar; caxis([0 1]);
    title('Node Activity Across Frequencies (1=active)');
    xlabel('Frequency Index'); ylabel('Node Index');
    xticks(1:F); xticklabels(string(opts.frequencies));

    % (4) Proxy heatmap at selected frequency (abs if complex), with τ line in colorbar label
    subplot(2,3,4);
    imagesc(Cshow);
    axis image; colorbar;
    title(sprintf('Proxy |C| at freq %d (\\tau=%.3g)', opts.frequencies(freq_index), tau));
    xlabel('j'); ylabel('i');

    % (5) Kept mask at selected frequency
    subplot(2,3,5);
    imagesc(double(K));
    axis image; colorbar; caxis([0 1]);
    title(sprintf('Kept Edges at freq %d', opts.frequencies(freq_index)));
    xlabel('j'); ylabel('i');

    % (6) Deactivated mask at selected frequency
    subplot(2,3,6);
    imagesc(double(OFF));
    axis image; colorbar; caxis([0 1]);
    title(sprintf('Deactivated Edges at freq %d', opts.frequencies(freq_index)));
    xlabel('j'); ylabel('i');

    sgtitle(opts.title, 'FontWeight', 'bold');

    % ---- Secondary quick views in separate figure (optional, helpful) ----
    figure('Name','Module 3 – Node degrees & Proxy histogram','Color','w','Position',[180,120,1100,420]);

    % Node degrees (kept) at selected frequency
    subplot(1,2,1);
    deg_kept = sum(K & ~eye(size(K)), 2);
    bar(deg_kept);
    title(sprintf('Kept Degree per Node (freq %d)', opts.frequencies(freq_index)));
    xlabel('Node'); ylabel('Degree'); grid on;

    % Proxy histogram (+ τ)
    subplot(1,2,2);
    if ~isempty(proxy_vals)
        histogram(proxy_vals, max(15, round(sqrt(numel(proxy_vals)))));
        hold on;
        yl = ylim;
        plot([tau tau], yl, 'r--', 'LineWidth', 2);
        hold off;
        title(sprintf('Proxy Distribution (upper-tri), freq %d', opts.frequencies(freq_index)));
        xlabel('|c_{ij}|'); ylabel('Count');
        legend({'Proxies',['\tau = ' num2str(tau)]}, 'Location','best');
        grid on;
    else
        text(0.5,0.5,'No finite proxy values','HorizontalAlignment','center');
        axis off;
    end
end
