function visualize_simulation_data(true_precision, true_covariance, emp_covariance, params)
% VISUALIZE_SIMULATION_DATA - Comprehensive visualization of generated simulation data
%
% This function creates multiple figures to visualize different aspects of the
% simulated precision matrices, covariances, and their properties.
%
% Inputs:
%   true_precision  - Cell array of true precision matrices
%   true_covariance - Cell array of true covariance matrices  
%   emp_covariance  - Cell array of empirical covariance matrices
%   params          - Parameters structure from module7_simulation
%
% Outputs:
%   Creates 2-3 figures with multiple subplots
%
% Example:
%   [Omega, Sigma, Sigma_emp, params] = module7_simulation('n_nodes', 10);
%   visualize_simulation_data(Omega, Sigma, Sigma_emp, params);

    % Input validation
    if nargin < 4
        error('All four inputs required: true_precision, true_covariance, emp_covariance, params');
    end
    
    % Figure 1: Overview of data structure and matrices
    figure('Name', 'Simulation Data Overview', 'Position', [100, 100, 1200, 800]);
    
    % 1.1: Graph structure
    subplot(2, 3, 1);
    plot_graph_structure(params);
    title('Graph Structure', 'FontSize', 12);
    
    % 1.2-1.4: Precision matrices at different frequencies
    freq_indices = round(linspace(1, params.n_freq, 3));
    for i = 1:3
        subplot(2, 3, i+1);
        plot_precision_matrix(true_precision{freq_indices(i)});
        title(sprintf('|Precision| at freq %d/%d', freq_indices(i), params.n_freq), 'FontSize', 12);
    end
    
    % 1.5: Frequency smoothness
    subplot(2, 3, 5);
    plot_frequency_smoothness(true_precision, params);
    title('Edge Values Across Frequencies', 'FontSize', 12);
    
    % 1.6: Covariance comparison
    subplot(2, 3, 6);
    plot_covariance_comparison(true_covariance, emp_covariance);
    title('True vs Empirical Covariance Error', 'FontSize', 12);
    
        % Figure 2: Detailed analysis
    figure('Name', 'Detailed Analysis', 'Position', [100, 100, 1200, 600]);
    
    % 2.1: Basis functions
    subplot(1, 3, 1);
    plot_basis_functions(params);
    title('Smooth Basis Functions', 'FontSize', 12);
    
    % 2.2: Sparsity pattern
    subplot(1, 3, 2);
    plot_sparsity_pattern(true_precision);
    title('Average Sparsity Pattern', 'FontSize', 12);
    
    % 2.3: Eigenvalue distribution
    subplot(1, 3, 3);
    plot_eigenvalue_distribution(true_precision, true_covariance);
    title('Eigenvalue Distribution', 'FontSize', 12);
    
    % Figure 3: Statistical properties (optional, for detailed analysis)
    if params.n_freq >= 10
        figure('Name', 'Statistical Properties', 'Position', [100, 100, 1000, 600]);
        
        % 3.1: Condition numbers
        subplot(2, 2, 1);
        plot_condition_numbers(true_precision, true_covariance);
        title('Condition Numbers', 'FontSize', 12);
        
        % 3.2: Edge weight distribution
        subplot(2, 2, 2);
        plot_edge_distribution(true_precision);
        title('Edge Weight Distribution', 'FontSize', 12);
        
        % 3.3: Empirical vs true diagonal
        subplot(2, 2, 3);
        plot_diagonal_comparison(true_covariance, emp_covariance);
        title('Diagonal Elements: True vs Empirical', 'FontSize', 12);
        
        % 3.4: Frobenius norm errors
        subplot(2, 2, 4);
        plot_frobenius_errors(true_covariance, emp_covariance);
        title('Frobenius Norm Errors', 'FontSize', 12);
    end
end

%% Subfunctions for individual plots

function plot_graph_structure(params)
    % Create and visualize the graph structure
    n = params.n_nodes;
    A = zeros(n, n);
    
    % Build adjacency matrix from edge list
    for e = 1:size(params.edge_list, 1)
        i = params.edge_list(e, 1);
        j = params.edge_list(e, 2);
        A(i, j) = 1;
        A(j, i) = 1;
    end
    
    % Try to use graph plotting if available
    if exist('graph', 'class') == 8
        G = graph(A);
        h = plot(G, 'Layout', 'force', 'NodeColor', [0.2, 0.4, 0.8], ...
                 'MarkerSize', 8, 'LineWidth', 2);
        
        % Highlight hub nodes if applicable
        degrees = sum(A, 2);
        [max_degree, hub_idx] = max(degrees);
        if max_degree > mean(degrees) + 2*std(degrees)
            highlight(h, hub_idx, 'NodeColor', 'r', 'MarkerSize', 12);
        end
    else
        % Fallback: show adjacency matrix
        imagesc(A);
        colormap(flipud(gray));
        axis square;
        xlabel('Node index');
        ylabel('Node index');
    end
    
    % Add text information
    text(0.02, 0.98, sprintf('Type: %s\nNodes: %d\nEdges: %d', ...
         params.graph_type, params.n_nodes, params.n_edges), ...
         'Units', 'normalized', 'VerticalAlignment', 'top', ...
         'BackgroundColor', 'w', 'FontSize', 10);
end

function plot_precision_matrix(Omega)
    % Visualize a single precision matrix
    imagesc(log10(abs(Omega) + 1e-10));
    colormap(hot);
    colorbar;
    axis square;
    xlabel('Node index');
    ylabel('Node index');
    
    % Add grid for small matrices
    if size(Omega, 1) <= 20
        grid on;
        set(gca, 'GridColor', 'w', 'GridAlpha', 0.3);
    end
end

function plot_frequency_smoothness(true_precision, params)
    % Track specific edges across frequencies
    n_edges_to_plot = min(5, params.n_edges);
    edge_values = zeros(params.n_freq, n_edges_to_plot);
    
    % Select edges to track (first few from edge list)
    if n_edges_to_plot > 0
        edge_indices = params.edge_list(1:n_edges_to_plot, :);
        
        for f = 1:params.n_freq
            for e = 1:n_edges_to_plot
                i = edge_indices(e, 1);
                j = edge_indices(e, 2);
                edge_values(f, e) = true_precision{f}(i, j);
            end
        end
        
        % Plot with different line styles
        freq_axis = 1:params.n_freq;
        colors = lines(n_edges_to_plot);
        hold on;
        for e = 1:n_edges_to_plot
            plot(freq_axis, edge_values(:, e), 'LineWidth', 2, ...
                 'Color', colors(e, :), ...
                 'DisplayName', sprintf('Edge (%d,%d)', edge_indices(e,1), edge_indices(e,2)));
        end
        hold off;
        
        xlabel('Frequency index');
        ylabel('Edge value');
        legend('Location', 'best');
        grid on;
    else
        text(0.5, 0.5, 'No edges to display', 'HorizontalAlignment', 'center');
    end
end

function plot_covariance_comparison(true_covariance, emp_covariance)
    % Compare true and empirical covariances
    n_freq = length(true_covariance);
    relative_errors = zeros(n_freq, 1);
    absolute_errors = zeros(n_freq, 1);
    
    for f = 1:n_freq
        diff_mat = emp_covariance{f} - true_covariance{f};
        relative_errors(f) = norm(diff_mat, 'fro') / norm(true_covariance{f}, 'fro');
        absolute_errors(f) = norm(diff_mat, 'fro');
    end
    
    % Plot both error types
    yyaxis left;
    plot(1:n_freq, relative_errors, 'b-', 'LineWidth', 2);
    ylabel('Relative Error');
    ylim([0, max(relative_errors)*1.1]);
    
    yyaxis right;
    plot(1:n_freq, absolute_errors, 'r--', 'LineWidth', 2);
    ylabel('Absolute Error (Frobenius)');
    ylim([0, max(absolute_errors)*1.1]);
    
    xlabel('Frequency index');
    grid on;
    
    % Add mean error info
    text(0.02, 0.98, sprintf('Mean rel. error: %.3f\nMean abs. error: %.3f', ...
         mean(relative_errors), mean(absolute_errors)), ...
         'Units', 'normalized', 'VerticalAlignment', 'top', ...
         'BackgroundColor', 'w');
end

function plot_basis_functions(params)
    % Visualize the smooth basis functions
    if isfield(params, 'basis_functions') && ~isempty(params.basis_functions)
        basis = params.basis_functions;
        freq_axis = 1:size(basis, 1);
        
        % Plot each basis function
        colors = lines(params.n_basis);
        hold on;
        for b = 1:params.n_basis
            plot(freq_axis, basis(:, b), 'LineWidth', 2, ...
                 'Color', colors(b, :), ...
                 'DisplayName', sprintf('Basis %d', b));
        end
        hold off;
        
        xlabel('Frequency index');
        ylabel('Basis value');
        legend('Location', 'best');
        grid on;
        
        % Show sum of basis functions
        sum_basis = sum(basis, 2);
        if max(abs(sum_basis - mean(sum_basis))) > 0.01
            hold on;
            plot(freq_axis, sum_basis, 'k:', 'LineWidth', 2, ...
                 'DisplayName', 'Sum');
            hold off;
        end
    else
        text(0.5, 0.5, 'Basis functions not available', ...
             'HorizontalAlignment', 'center');
    end
end

function plot_sparsity_pattern(true_precision)
    % Show average sparsity pattern across frequencies
    n_freq = length(true_precision);
    n_nodes = size(true_precision{1}, 1);
    
    % Compute average absolute values
    avg_abs_precision = zeros(n_nodes, n_nodes);
    for f = 1:n_freq
        avg_abs_precision = avg_abs_precision + abs(true_precision{f});
    end
    avg_abs_precision = avg_abs_precision / n_freq;
    
    % Threshold to show pattern
    threshold = 0.01 * max(avg_abs_precision(:));
    sparsity_pattern = avg_abs_precision > threshold;
    
    % Plot
    imagesc(sparsity_pattern);
    colormap(flipud(gray));
    axis square;
    xlabel('Node index');
    ylabel('Node index');
    
    % Calculate and display sparsity
    off_diag_mask = ~eye(n_nodes);
    sparsity = nnz(sparsity_pattern(off_diag_mask)) / nnz(off_diag_mask);
    title(sprintf('Sparsity: %.1f%%', 100*sparsity), 'FontSize', 10);
    
    % Add grid for small matrices
    if n_nodes <= 20
        grid on;
        set(gca, 'GridColor', 'k', 'GridAlpha', 0.3);
    end
end

function plot_eigenvalue_distribution(true_precision, true_covariance)
    % Plot eigenvalues across frequencies
    n_freq = length(true_precision);
    n_nodes = size(true_precision{1}, 1);
    
    % Collect eigenvalues
    eigvals_prec = zeros(n_nodes, n_freq);
    eigvals_cov = zeros(n_nodes, n_freq);
    
    for f = 1:n_freq
        eigvals_prec(:, f) = sort(eig(true_precision{f}), 'descend');
        eigvals_cov(:, f) = sort(eig(true_covariance{f}), 'descend');
    end
    
    % Plot on log scale
    freq_axis = 1:n_freq;
    
    % Use color gradient for eigenvalues
    colors_prec = copper(n_nodes);
    colors_cov = winter(n_nodes);
    
    hold on;
    % Plot precision eigenvalues
    for i = 1:n_nodes
        semilogy(freq_axis, eigvals_prec(i, :), '-', ...
                'Color', colors_prec(i, :), 'LineWidth', 1.5);
    end
    
    % Plot covariance eigenvalues
    for i = 1:n_nodes
        semilogy(freq_axis, eigvals_cov(i, :), '--', ...
                'Color', colors_cov(i, :), 'LineWidth', 1.5);
    end
    hold off;
    
    xlabel('Frequency index');
    ylabel('Eigenvalue (log scale)');
    grid on;
    
    % Add legend for first and last
    legend({'Prec (max)', '', 'Prec (min)', 'Cov (max)', '', 'Cov (min)'}, ...
           'Location', 'best');
    
    % Add condition number info
    cond_prec = eigvals_prec(1, :) ./ eigvals_prec(end, :);
    cond_cov = eigvals_cov(1, :) ./ eigvals_cov(end, :);
    text(0.02, 0.02, sprintf('Mean κ: Prec=%.1f, Cov=%.1f', ...
         mean(cond_prec), mean(cond_cov)), ...
         'Units', 'normalized', 'VerticalAlignment', 'bottom', ...
         'BackgroundColor', 'w');
end

%% Additional plotting functions for Figure 3

function plot_condition_numbers(true_precision, true_covariance)
    % Plot condition numbers across frequencies
    n_freq = length(true_precision);
    cond_prec = zeros(n_freq, 1);
    cond_cov = zeros(n_freq, 1);
    
    for f = 1:n_freq
        cond_prec(f) = cond(true_precision{f});
        cond_cov(f) = cond(true_covariance{f});
    end
    
    freq_axis = 1:n_freq;
    semilogy(freq_axis, cond_prec, 'b-', 'LineWidth', 2, 'DisplayName', 'Precision');
    hold on;
    semilogy(freq_axis, cond_cov, 'r--', 'LineWidth', 2, 'DisplayName', 'Covariance');
    hold off;
    
    xlabel('Frequency index');
    ylabel('Condition number (log scale)');
    legend('Location', 'best');
    grid on;
    
    % Add statistics
    text(0.02, 0.98, sprintf('Precision: %.1f ± %.1f\nCovariance: %.1f ± %.1f', ...
         mean(cond_prec), std(cond_prec), mean(cond_cov), std(cond_cov)), ...
         'Units', 'normalized', 'VerticalAlignment', 'top', ...
         'BackgroundColor', 'w');
end

function plot_edge_distribution(true_precision)
    % Plot distribution of edge weights
    n_freq = length(true_precision);
    n_nodes = size(true_precision{1}, 1);
    
    % Collect all off-diagonal values
    all_edges = [];
    for f = 1:n_freq
        P = true_precision{f};
        mask = triu(true(n_nodes), 1); % Upper triangular mask
        all_edges = [all_edges; P(mask)];
    end
    
    % Remove near-zero values
    threshold = 1e-10;
    significant_edges = all_edges(abs(all_edges) > threshold);
    
    % Plot histogram
    if ~isempty(significant_edges)
        histogram(log10(abs(significant_edges)), 30, 'Normalization', 'probability');
        xlabel('log_{10}(|edge weight|)');
        ylabel('Probability');
        grid on;
        
        % Add statistics
        text(0.02, 0.98, sprintf('Non-zeros: %d\nMean: %.3f\nStd: %.3f', ...
             length(significant_edges), mean(abs(significant_edges)), ...
             std(abs(significant_edges))), ...
             'Units', 'normalized', 'VerticalAlignment', 'top', ...
             'BackgroundColor', 'w');
    else
        text(0.5, 0.5, 'No significant edges found', ...
             'HorizontalAlignment', 'center');
    end
end

function plot_diagonal_comparison(true_covariance, emp_covariance)
    % Compare diagonal elements of true vs empirical covariance
    n_freq = length(true_covariance);
    n_nodes = size(true_covariance{1}, 1);
    
    % Select a few nodes to track
    nodes_to_plot = min(5, n_nodes);
    node_indices = round(linspace(1, n_nodes, nodes_to_plot));
    
    colors = lines(nodes_to_plot);
    freq_axis = 1:n_freq;
    
    for i = 1:nodes_to_plot
        node = node_indices(i);
        true_diag = zeros(n_freq, 1);
        emp_diag = zeros(n_freq, 1);
        
        for f = 1:n_freq
            true_diag(f) = true_covariance{f}(node, node);
            emp_diag(f) = emp_covariance{f}(node, node);
        end
        
        plot(freq_axis, true_diag, '-', 'Color', colors(i,:), ...
             'LineWidth', 2, 'DisplayName', sprintf('True (node %d)', node));
        hold on;
        plot(freq_axis, emp_diag, ':', 'Color', colors(i,:), ...
             'LineWidth', 2, 'DisplayName', sprintf('Emp (node %d)', node));
    end
    hold off;
    
    xlabel('Frequency index');
    ylabel('Diagonal value');
    legend('Location', 'best', 'FontSize', 8);
    grid on;
end

function plot_frobenius_errors(true_covariance, emp_covariance)
    % Plot detailed error analysis
    n_freq = length(true_covariance);
    
    rel_errors = zeros(n_freq, 1);
    diag_errors = zeros(n_freq, 1);
    off_diag_errors = zeros(n_freq, 1);
    
    for f = 1:n_freq
        n = size(true_covariance{f}, 1);
        diff = emp_covariance{f} - true_covariance{f};
        
        % Overall relative error
        rel_errors(f) = norm(diff, 'fro') / norm(true_covariance{f}, 'fro');
        
        % Diagonal error
        diag_diff = diag(diag(diff));
        diag_errors(f) = norm(diag_diff, 'fro') / norm(diag(diag(true_covariance{f})), 'fro');
        
        % Off-diagonal error
        off_diag_mask = ~eye(n);
        off_diag_diff = diff .* off_diag_mask;
        off_diag_true = true_covariance{f} .* off_diag_mask;
        off_diag_errors(f) = norm(off_diag_diff, 'fro') / (norm(off_diag_true, 'fro') + eps);
    end
    
    freq_axis = 1:n_freq;
    plot(freq_axis, rel_errors, 'k-', 'LineWidth', 2, 'DisplayName', 'Overall');
    hold on;
    plot(freq_axis, diag_errors, 'b--', 'LineWidth', 2, 'DisplayName', 'Diagonal');
    plot(freq_axis, off_diag_errors, 'r:', 'LineWidth', 2, 'DisplayName', 'Off-diagonal');
    hold off;
    
    xlabel('Frequency index');
    ylabel('Relative Error');
    legend('Location', 'best');
    grid on;
    ylim([0, max([rel_errors; diag_errors; off_diag_errors])*1.1]);
end