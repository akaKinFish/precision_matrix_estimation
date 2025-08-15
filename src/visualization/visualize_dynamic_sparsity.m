function visualize_dynamic_sparsity(true_precision, params)
% VISUALIZE_DYNAMIC_SPARSITY - Visualize how sparsity pattern changes across frequencies
%
% This function creates visualizations specifically for understanding
% dynamic sparsity patterns in the improved simulation.
%
% Inputs:
%   true_precision - Cell array of precision matrices
%   params - Parameters including variable edge information

    n_freq = length(true_precision);
    n_nodes = size(true_precision{1}, 1);
    
    % Figure 1: Sparsity evolution
    figure('Name', 'Dynamic Sparsity Evolution', 'Position', [100, 100, 1200, 800]);
    
    % Select frequency points to display
    freq_points = round(linspace(1, n_freq, 6));
    
    % Plot sparsity patterns at different frequencies
    for i = 1:6
        subplot(2, 3, i);
        f = freq_points(i);
        
        % Create sparsity pattern (with magnitude information)
        P = true_precision{f};
        max_val = max(abs(P(:)));
        threshold = 1e-3 * max_val;
        
        % Use different colors for different magnitude ranges
        pattern = zeros(n_nodes, n_nodes, 3); % RGB
        
        for row = 1:n_nodes
            for col = 1:n_nodes
                val = abs(P(row, col));
                if val > threshold
                    % Color based on magnitude
                    if val > 0.5 * max_val
                        pattern(row, col, :) = [1, 0, 0]; % Red for strong
                    elseif val > 0.1 * max_val
                        pattern(row, col, :) = [1, 1, 0]; % Yellow for medium
                    else
                        pattern(row, col, :) = [0, 0, 1]; % Blue for weak
                    end
                end
            end
        end
        
        image(pattern);
        title(sprintf('Freq %d/%d', f, n_freq));
        axis square;
        set(gca, 'XTick', [], 'YTick', []);
    end
    
    % Add colorbar legend
    subplot(2, 3, 6);
    text(0.1, 0.8, 'Red: Strong edges', 'Color', 'r', 'FontWeight', 'bold');
    text(0.1, 0.6, 'Yellow: Medium edges', 'Color', [0.8, 0.8, 0], 'FontWeight', 'bold');
    text(0.1, 0.4, 'Blue: Weak edges', 'Color', 'b', 'FontWeight', 'bold');
    text(0.1, 0.2, 'Black: Zero/Near-zero', 'Color', 'k', 'FontWeight', 'bold');
    axis off;
    
    % Figure 2: Edge tracking
    figure('Name', 'Edge Evolution Tracking', 'Position', [100, 100, 1000, 600]);
    
    % Track specific edges
    if isfield(params, 'variable_edge_list') && ~isempty(params.variable_edge_list)
        % Track both base and variable edges
        n_base_track = min(3, size(params.base_edge_list, 1));
        n_var_track = min(3, size(params.variable_edge_list, 1));
        
        subplot(2, 1, 1);
        hold on;
        
        % Plot base edges (should be consistently active)
        for e = 1:n_base_track
            i = params.base_edge_list(e, 1);
            j = params.base_edge_list(e, 2);
            values = zeros(n_freq, 1);
            for f = 1:n_freq
                values(f) = true_precision{f}(i, j);
            end
            plot(1:n_freq, values, 'LineWidth', 2, ...
                 'DisplayName', sprintf('Base edge (%d,%d)', i, j));
        end
        
        title('Base Edges (Always Active)');
        xlabel('Frequency');
        ylabel('Edge Value');
        legend('Location', 'best');
        grid on;
        hold off;
        
        subplot(2, 1, 2);
        hold on;
        
        % Plot variable edges (should show activation/deactivation)
        for e = 1:n_var_track
            i = params.variable_edge_list(e, 1);
            j = params.variable_edge_list(e, 2);
            values = zeros(n_freq, 1);
            for f = 1:n_freq
                values(f) = true_precision{f}(i, j);
            end
            plot(1:n_freq, values, 'LineWidth', 2, ...
                 'DisplayName', sprintf('Variable edge (%d,%d)', i, j));
        end
        
        title('Variable Edges (Dynamic Activation)');
        xlabel('Frequency');
        ylabel('Edge Value');
        legend('Location', 'best');
        grid on;
        hold off;
    else
        % Fallback: just track some edges
        edges_to_track = min(6, params.n_edges);
        edge_list = params.edge_list(1:edges_to_track, :);
        
        hold on;
        colors = lines(edges_to_track);
        
        for e = 1:edges_to_track
            i = edge_list(e, 1);
            j = edge_list(e, 2);
            values = zeros(n_freq, 1);
            active = zeros(n_freq, 1);
            
            for f = 1:n_freq
                values(f) = true_precision{f}(i, j);
                active(f) = abs(values(f)) > 1e-10;
            end
            
            % Plot with different style for active/inactive
            plot(1:n_freq, values, 'Color', colors(e, :), 'LineWidth', 2, ...
                 'DisplayName', sprintf('Edge (%d,%d)', i, j));
            
            % Mark transitions
            transitions = find(diff(active) ~= 0) + 0.5;
            for t = transitions'
                plot([t, t], ylim, 'k:', 'LineWidth', 1, 'HandleVisibility', 'off');
            end
        end
        
        title('Edge Evolution with Activation Changes');
        xlabel('Frequency');
        ylabel('Edge Value');
        legend('Location', 'best');
        grid on;
        hold off;
    end
    
    % Figure 3: Sparsity statistics
    figure('Name', 'Sparsity Statistics', 'Position', [100, 100, 800, 600]);
    
    % Compute sparsity metrics
    n_edges_per_freq = zeros(n_freq, 1);
    total_variation = zeros(n_freq-1, 1);
    threshold = 1e-10;
    
    for f = 1:n_freq
        P = true_precision{f};
        mask = triu(true(n_nodes), 1);
        n_edges_per_freq(f) = nnz(abs(P(mask)) > threshold);
        
        if f > 1
            P_prev = true_precision{f-1};
            total_variation(f-1) = sum(sum(abs(P - P_prev))) / 2; % Symmetric
        end
    end
    
    subplot(2, 2, 1);
    plot(1:n_freq, n_edges_per_freq, 'b-', 'LineWidth', 2);
    xlabel('Frequency');
    ylabel('Number of Active Edges');
    title('Edge Count Evolution');
    grid on;
    
    % Add expected range
    if isfield(params, 'n_base_edges') && isfield(params, 'n_variable_edges')
        hold on;
        yline(params.n_base_edges, 'r--', 'Base edges');
        yline(params.n_base_edges + params.n_variable_edges, 'g--', 'Max possible');
        hold off;
    end
    
    subplot(2, 2, 2);
    if length(total_variation) > 0
        plot(2:n_freq, total_variation, 'r-', 'LineWidth', 2);
        xlabel('Frequency');
        ylabel('Total Variation');
        title('Matrix Change Between Adjacent Frequencies');
        grid on;
    end
    
    subplot(2, 2, 3);
    % Histogram of edge activation frequencies
    edge_activation_freq = zeros(n_nodes * (n_nodes - 1) / 2, 1);
    idx = 1;
    
    for i = 1:n_nodes
        for j = i+1:n_nodes
            activation_count = 0;
            for f = 1:n_freq
                if abs(true_precision{f}(i, j)) > threshold
                    activation_count = activation_count + 1;
                end
            end
            edge_activation_freq(idx) = activation_count / n_freq;
            idx = idx + 1;
        end
    end
    
    histogram(edge_activation_freq, 20);
    xlabel('Activation Frequency');
    ylabel('Number of Edges');
    title('Edge Activation Distribution');
    grid on;
    
    subplot(2, 2, 4);
    % Sparsity pattern stability
    pattern_changes = zeros(n_freq-1, 1);
    for f = 2:n_freq
        prev_pattern = abs(true_precision{f-1}) > threshold;
        curr_pattern = abs(true_precision{f}) > threshold;
        pattern_changes(f-1) = nnz(prev_pattern ~= curr_pattern) / 2; % Symmetric
    end
    
    bar(2:n_freq, pattern_changes);
    xlabel('Frequency');
    ylabel('Number of Pattern Changes');
    title('Sparsity Pattern Changes');
    grid on;
    
    % Add summary statistics
    fprintf('\n=== Dynamic Sparsity Summary ===\n');
    fprintf('Average active edges: %.1f (±%.1f)\n', mean(n_edges_per_freq), std(n_edges_per_freq));
    fprintf('Total pattern changes: %d\n', sum(pattern_changes));
    fprintf('Edges always active: %d\n', sum(edge_activation_freq == 1));
    fprintf('Edges never active: %d\n', sum(edge_activation_freq == 0));
    fprintf('Edges sometimes active: %d\n', sum(edge_activation_freq > 0 & edge_activation_freq < 1));
    fprintf('================================\n');
end