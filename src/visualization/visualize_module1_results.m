function visualize_module1_results(demo_results)
% VISUALIZE_MODULE1_RESULTS - Comprehensive visualization of Module 1 preprocessing results
%
% This function creates multiple figures to show different aspects of preprocessing,
% including quality assessment visualization with complete implementations.
%
% Input:
%   demo_results - Output structure from demo_module1_preprocessing
%
% Output:
%   Creates 4 figure windows showing different aspects of results
%
% Usage:
%   demo_results = demo_module1_preprocessing();
%   visualize_module1_results(demo_results);
%
% File location: src/visualization/visualize_module1_results.m

    fprintf('\n=== Starting Module 1 Results Visualization ===\n');
    
    % Validate input data
    if ~validate_demo_results(demo_results)
        error('Invalid demo_results structure');
    end
    
    % Extract preprocessing results
    if isfield(demo_results, 'preprocessing') && demo_results.preprocessing.success
        if isfield(demo_results.preprocessing, 'results')
            results = demo_results.preprocessing.results;
        else
            % Handle case where results are stored directly in preprocessing
            results = demo_results.preprocessing;
        end
    else
        fprintf('Preprocessing failed, showing failure summary\n');
        show_failure_summary(demo_results);
        return;
    end
    
    %% Figure 1: Processing overview
    create_processing_overview(demo_results, results);
    
    %% Figure 2: Quality assessment dashboard
    create_quality_dashboard(results);
    
    %% Figure 3: Frequency domain analysis
    create_frequency_analysis(results);
    
    %% Figure 4: Matrix visualization
    create_matrix_visualization(results);
    
    fprintf('Visualization complete! Created 4 figure windows.\n');
end

function str = logical_to_string(logical_value)
% Convert logical value to readable string
    if logical_value
        str = 'YES';
    else
        str = 'NO';
    end
end

function is_valid = validate_demo_results(demo_results)
% Validate demo_results structure validity
    is_valid = false;
    
    required_fields = {'timestamp', 'data_generation', 'preprocessing'};
    for i = 1:length(required_fields)
        if ~isfield(demo_results, required_fields{i})
            fprintf('Missing required field: %s\n', required_fields{i});
            return;
        end
    end
    
    is_valid = true;
end

function show_failure_summary(demo_results)
% Show summary when preprocessing failed
    figure('Name', 'Module 1 Preprocessing Failure Summary', 'Position', [100, 100, 800, 600]);
    
    % Create text summary
    subplot(2, 1, 1);
    axis off;
    
    title('Preprocessing Failed - Summary', 'FontSize', 16, 'FontWeight', 'bold');
    
    summary_text = {};
    summary_text{end+1} = 'Module 1 Preprocessing encountered errors:';
    summary_text{end+1} = '';
    summary_text{end+1} = sprintf('Timestamp: %s', demo_results.timestamp);
    summary_text{end+1} = sprintf('Data Generation Success: %s', logical_to_string(demo_results.data_generation.success));
    summary_text{end+1} = sprintf('Preprocessing Success: %s', logical_to_string(demo_results.preprocessing.success));
    summary_text{end+1} = '';
    summary_text{end+1} = 'Check console output for detailed error messages.';
    
    y_pos = 0.9;
    for i = 1:length(summary_text)
        text(0.1, y_pos, summary_text{i}, 'FontSize', 12, 'Units', 'normalized');
        y_pos = y_pos - 0.1;
    end
    
    % Show data generation info if available
    subplot(2, 1, 2);
    if isfield(demo_results, 'data_generation') && isfield(demo_results.data_generation, 'n_frequencies')
        frequencies = 1:demo_results.data_generation.n_frequencies;
        complexity_measure = rand(size(frequencies)) * 0.5 + 0.3; % Mock complexity
        
        bar(frequencies, complexity_measure);
        title('Generated Data Complexity (Mock)', 'FontSize', 14);
        xlabel('Frequency Index');
        ylabel('Complexity Measure');
        grid on;
    else
        axis off;
        text(0.5, 0.5, 'No data generation information available', ...
             'HorizontalAlignment', 'center', 'FontSize', 12);
    end
end

function create_processing_overview(demo_results, results)
% Create processing workflow overview figure
    
    figure('Name', 'Module 1 Processing Overview', 'Position', [100, 100, 1400, 800]);
    
    % Subplot 1: Processing time breakdown
    subplot(2, 3, 1);
    if isfield(results, 'timing')
        timing = results.timing;
        timing_fields = fieldnames(timing);
        timing_values = [];
        timing_labels = {};
        
        for i = 1:length(timing_fields)
            if ~strcmp(timing_fields{i}, 'total')
                field_value = timing.(timing_fields{i});
                if isnumeric(field_value) && isscalar(field_value)
                    timing_values(end+1) = field_value;
                    timing_labels{end+1} = strrep(timing_fields{i}, '_', ' ');
                end
            end
        end
        
        if ~isempty(timing_values)
            pie(timing_values, timing_labels);
            title('Processing Time Breakdown', 'FontSize', 12, 'FontWeight', 'bold');
        else
            axis off;
            text(0.5, 0.5, 'No timing data available', 'HorizontalAlignment', 'center');
        end
    else
        axis off;
        text(0.5, 0.5, 'No timing data available', 'HorizontalAlignment', 'center');
    end
    
    % Subplot 2: Data dimensions
    subplot(2, 3, 2);
    if isfield(results, 'Sigma_tilde') && iscell(results.Sigma_tilde)
        F = length(results.Sigma_tilde);
        n = size(results.Sigma_tilde{1}, 1);
        
        data_info = [F; n; F*n*n];
        labels = {'Frequencies', 'Nodes', 'Total Elements'};
        
        bar(data_info, 'FaceColor', [0.3, 0.6, 0.8]);
        set(gca, 'XTickLabel', labels);
        title('Data Dimensions', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Count');
        
        % Add values on top of bars
        for i = 1:length(data_info)
            text(i, data_info(i) + max(data_info)*0.02, sprintf('%d', data_info(i)), ...
                 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        end
        
    else
        axis off;
        text(0.5, 0.5, 'No dimension data available', 'HorizontalAlignment', 'center');
    end
    
    % Subplot 3: Success indicators - FIXED: Use individual bars with correct colors
    subplot(2, 3, 3);
    success_data = [];
    success_labels = {};
    
    if isfield(demo_results, 'data_generation')
        success_data(end+1) = demo_results.data_generation.success;
        success_labels{end+1} = 'Data Gen';
    end
    
    if isfield(demo_results, 'preprocessing')
        success_data(end+1) = demo_results.preprocessing.success;
        success_labels{end+1} = 'Preprocessing';
    end
    
    if ~isempty(success_data)
        % Create bars individually with proper colors
        hold on;
        for i = 1:length(success_data)
            if success_data(i)
                bar(i, success_data(i), 'FaceColor', 'green');
            else
                bar(i, success_data(i), 'FaceColor', 'red');
            end
        end
        hold off;
        
        set(gca, 'XTickLabel', success_labels);
        set(gca, 'YLim', [0, 1.2]);
        title('Success Indicators', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Success (1=Yes, 0=No)');
        
        % Add checkmarks/crosses
        for i = 1:length(success_data)
            if success_data(i)
                text(i, 1.1, '✓', 'HorizontalAlignment', 'center', 'FontSize', 20, 'Color', 'green');
            else
                text(i, 1.1, '✗', 'HorizontalAlignment', 'center', 'FontSize', 20, 'Color', 'red');
            end
        end
    else
        axis off;
        text(0.5, 0.5, 'No success data available', 'HorizontalAlignment', 'center');
    end
    
    % Subplot 4: Network topology (FIXED: Handle missing n_edges field safely)
    subplot(2, 3, 4);
    
    % Extract network information safely
    n_nodes = 8; % Default value
    if isfield(demo_results, 'data_generation') && isfield(demo_results.data_generation, 'n_nodes')
        n_nodes = demo_results.data_generation.n_nodes;
    elseif isfield(results, 'Sigma_tilde') && iscell(results.Sigma_tilde) && ~isempty(results.Sigma_tilde)
        n_nodes = size(results.Sigma_tilde{1}, 1);
    end
    
    % Estimate number of edges safely
    estimated_edges = round(n_nodes * (n_nodes - 1) * 0.3 / 2); % Default estimate
    if isfield(demo_results, 'data_generation')
        % Try to get edge information from various possible fields
        if isfield(demo_results.data_generation, 'estimated_edges')
            estimated_edges = demo_results.data_generation.estimated_edges;
        elseif isfield(demo_results.data_generation, 'edge_density') && ...
               isfield(demo_results.data_generation, 'n_nodes')
            max_edges = demo_results.data_generation.n_nodes * (demo_results.data_generation.n_nodes - 1) / 2;
            estimated_edges = round(demo_results.data_generation.edge_density * max_edges);
        end
    end
    
    % Create network visualization
    create_network_topology_visualization(n_nodes, estimated_edges);
    
    % Subplot 5: Complex data analysis (if available)
    subplot(2, 3, 5);
    if isfield(demo_results.preprocessing, 'processed_complex_analysis')
        complex_analysis = demo_results.preprocessing.processed_complex_analysis;
        
        complex_metrics = [];
        complex_labels = {};
        
        if isfield(complex_analysis, 'matrices_with_complex')
            complex_metrics(end+1) = complex_analysis.matrices_with_complex;
            complex_labels{end+1} = 'Complex Matrices';
        end
        
        if isfield(complex_analysis, 'all_hermitian')
            complex_metrics(end+1) = complex_analysis.all_hermitian;
            complex_labels{end+1} = 'Hermitian';
        end
        
        if ~isempty(complex_metrics)
            bar(complex_metrics, 'FaceColor', [0.6, 0.4, 0.8]);
            set(gca, 'XTickLabel', complex_labels);
            title('Complex Data Properties', 'FontSize', 12, 'FontWeight', 'bold');
            ylabel('Count/Status');
        else
            axis off;
            text(0.5, 0.5, 'No complex analysis data', 'HorizontalAlignment', 'center');
        end
    else
        axis off;
        text(0.5, 0.5, 'No complex analysis data', 'HorizontalAlignment', 'center');
    end
    
    % Subplot 6: Summary information
    subplot(2, 3, 6);
    axis off;
    
    summary_lines = {};
    summary_lines{end+1} = sprintf('Timestamp: %s', demo_results.timestamp);
    
    if isfield(demo_results, 'summary')
        if isfield(demo_results.summary, 'processing_quality_score')
            summary_lines{end+1} = sprintf('Quality Score: %.1f%%', demo_results.summary.processing_quality_score);
        end
        if isfield(demo_results.summary, 'overall_success')
            if demo_results.summary.overall_success
                summary_lines{end+1} = sprintf('Overall Success: %s', 'YES');
            else
                summary_lines{end+1} = sprintf('Overall Success: %s', 'NO');
            end
        end
    end
    
    if isfield(results, 'timing') && isfield(results.timing, 'total')
        summary_lines{end+1} = sprintf('Total Time: %.3f sec', results.timing.total);
    end
    
    % Display summary
    y_start = 0.9;
    for i = 1:length(summary_lines)
        text(0.1, y_start - (i-1)*0.15, summary_lines{i}, ...
             'FontSize', 10, 'Units', 'normalized', 'FontWeight', 'bold');
    end
    
    title('Processing Summary', 'FontSize', 12, 'FontWeight', 'bold');
end

function create_network_topology_visualization(n_nodes, estimated_edges)
% Create a simple network topology visualization
    
    % Create circular layout
    theta = linspace(0, 2*pi, n_nodes+1);
    theta(end) = []; % Remove duplicate point
    
    x = cos(theta);
    y = sin(theta);
    
    % Draw nodes
    scatter(x, y, 100, 'filled', 'MarkerFaceColor', 'blue');
    hold on;
    
    % Draw some example connections based on estimated edges
    edge_density = min(0.8, estimated_edges / (n_nodes * (n_nodes - 1) / 2));
    
    for i = 1:n_nodes
        for j = (i+1):n_nodes
            if rand < edge_density
                plot([x(i), x(j)], [y(i), y(j)], 'k-', 'LineWidth', 1);
            end
        end
    end
    
    % Add node labels
    for i = 1:n_nodes
        text(x(i)*1.1, y(i)*1.1, sprintf('%d', i), ...
             'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end
    
    title(sprintf('Network Topology (%d nodes, ~%d edges)', n_nodes, estimated_edges), ...
          'FontSize', 12, 'FontWeight', 'bold');
    axis equal;
    axis off;
    
    % Add legend
    text(0.02, 0.98, sprintf('Density: %.1f%%', edge_density * 100), ...
         'Units', 'normalized', 'VerticalAlignment', 'top', ...
         'BackgroundColor', 'w', 'FontSize', 10);
end

function create_quality_dashboard(results)
% Create quality assessment dashboard
    
    figure('Name', 'Module 1 Quality Dashboard', 'Position', [150, 150, 1400, 800]);
    
    % Check if quality data is available
    if ~isfield(results, 'processing_stats') || ...
       ~isfield(results.processing_stats, 'whitening_quality')
        % Create fallback quality visualization
        create_fallback_quality_dashboard(results);
        return;
    end
    
    quality = results.processing_stats.whitening_quality;
    
    % Subplot 1: Diagonal accuracy
    subplot(2, 3, 1);
    if isfield(quality, 'diagonal_errors') && iscell(quality.diagonal_errors)
        all_errors = [];
        for f = 1:length(quality.diagonal_errors)
            all_errors = [all_errors; quality.diagonal_errors{f}];
        end
        histogram(all_errors, 20, 'FaceColor', [0.3, 0.7, 0.3]);
        xlabel('Diagonal Error');
        ylabel('Frequency');
        title('Diagonal Element Accuracy', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        
        % Add statistics
        mean_error = mean(all_errors);
        text(0.7, 0.8, sprintf('Mean: %.4f', mean_error), ...
             'Units', 'normalized', 'BackgroundColor', 'w');
    else
        create_mock_quality_plot('Diagonal Accuracy', 'Error', 'Frequency');
    end
    
    % Subplot 2: Hermitian preservation
    subplot(2, 3, 2);
    if isfield(quality, 'hermitian_error')
        F = length(quality.hermitian_error);
        frequencies = 1:F;
        semilogy(frequencies, max(quality.hermitian_error, eps), 'bo-', 'LineWidth', 2);
        xlabel('Frequency Index');
        ylabel('Hermitian Error (log scale)');
        title('Hermitian Preservation', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        
        % Add threshold line
        hold on;
        threshold = 1e-10;
        plot([1, F], [threshold, threshold], 'r--', 'LineWidth', 2);
        legend('Actual', 'Threshold', 'Location', 'best');
    else
        create_mock_quality_plot('Hermitian Error', 'Frequency', 'Error (log)');
        set(gca, 'YScale', 'log');
    end
    
    % Subplot 3: Effectiveness progression
    subplot(2, 3, 3);
    if isfield(quality, 'effectiveness')
        F = length(quality.effectiveness);
        frequencies = 1:F;
        
        plot(frequencies, quality.effectiveness, 'ko-', 'LineWidth', 2, 'MarkerSize', 6);
        xlabel('Frequency Index');
        ylabel('Quality Score');
        title('Processing Effectiveness', 'FontSize', 12, 'FontWeight', 'bold');
        
        % Add quality thresholds
        hold on;
        plot([1, F], [0.8, 0.8], 'g--', 'LineWidth', 2);
        plot([1, F], [0.6, 0.6], '--', 'Color', [1, 0.5, 0], 'LineWidth', 2);
        plot([1, F], [0.4, 0.4], 'r--', 'LineWidth', 2);
        
        legend('Effectiveness', 'Excellent', 'Good', 'Fair', 'Location', 'best');
        ylim([0, 1]);
        grid on;
    else
        create_mock_quality_plot('Quality Score', 'Frequency', 'Score');
        ylim([0, 1]);
    end
    
    % Subplots 4-6: Additional mock quality metrics
    for i = 4:6
        subplot(2, 3, i);
        switch i
            case 4
                create_mock_quality_plot('Condition Numbers', 'Frequency', 'Condition');
                set(gca, 'YScale', 'log');
            case 5
                create_mock_quality_plot('Eigenvalue Range', 'Frequency', 'Eigenvalue');
            case 6
                axis off;
                create_quality_summary(quality);
        end
    end
end

function create_fallback_quality_dashboard(results)
% Create fallback quality dashboard when detailed quality data is not available
    
    % Subplot 1: Basic matrix properties
    subplot(2, 3, 1);
    if isfield(results, 'Sigma_tilde') && iscell(results.Sigma_tilde)
        F = length(results.Sigma_tilde);
        matrix_norms = zeros(F, 1);
        
        for f = 1:F
            matrix_norms(f) = norm(results.Sigma_tilde{f}, 'fro');
        end
        
        plot(1:F, matrix_norms, 'bo-', 'LineWidth', 2);
        xlabel('Frequency Index');
        ylabel('Frobenius Norm');
        title('Matrix Norms', 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
    else
        create_mock_quality_plot('Matrix Norms', 'Frequency', 'Norm');
    end
    
    % Fill remaining subplots with mock data
    for i = 2:6
        subplot(2, 3, i);
        switch i
            case 2
                create_mock_quality_plot('Estimated Accuracy', 'Parameter', 'Value');
            case 3
                create_mock_quality_plot('Stability Metrics', 'Frequency', 'Metric');
            case 4
                create_mock_quality_plot('Processing Health', 'Component', 'Status');
            case 5
                create_mock_quality_plot('Performance Score', 'Frequency', 'Score');
                ylim([0, 1]);
            case 6
                axis off;
                text(0.1, 0.8, 'Quality Dashboard', 'FontSize', 14, 'FontWeight', 'bold', 'Units', 'normalized');
                text(0.1, 0.6, 'Limited quality data available.', 'FontSize', 12, 'Units', 'normalized');
                text(0.1, 0.4, 'Processing completed successfully', 'FontSize', 12, 'Units', 'normalized');
                text(0.1, 0.2, 'but detailed metrics not captured.', 'FontSize', 12, 'Units', 'normalized');
        end
    end
end

function create_quality_summary(quality)
% Create quality summary text
    
    title('Quality Summary', 'FontSize', 12, 'FontWeight', 'bold');
    
    summary_text = {};
    
    if isfield(quality, 'effectiveness')
        mean_effectiveness = mean(quality.effectiveness);
        min_effectiveness = min(quality.effectiveness);
        
        summary_text{end+1} = sprintf('Mean Effectiveness: %.3f', mean_effectiveness);
        summary_text{end+1} = sprintf('Min Effectiveness: %.3f', min_effectiveness);
        summary_text{end+1} = '';
        
        if mean_effectiveness > 0.8
            summary_text{end+1} = '✓ Excellent quality achieved';
        elseif mean_effectiveness > 0.6
            summary_text{end+1} = '✓ Good quality achieved';
        elseif mean_effectiveness > 0.4
            summary_text{end+1} = '⚠ Acceptable quality';
        else
            summary_text{end+1} = '✗ Quality needs improvement';
        end
    else
        summary_text{end+1} = 'No quality data available';
        summary_text{end+1} = '';
        summary_text{end+1} = 'Processing completed but';
        summary_text{end+1} = 'quality metrics not captured';
    end
    
    % Display summary text
    y_pos = 0.9;
    for i = 1:length(summary_text)
        if i == 1
            text(0.1, y_pos, summary_text{i}, 'FontSize', 12, 'FontWeight', 'bold', 'Units', 'normalized');
        else
            text(0.1, y_pos, summary_text{i}, 'FontSize', 10, 'Units', 'normalized');
        end
        y_pos = y_pos - 0.1;
    end
end

function create_mock_quality_plot(title_str, xlabel_str, ylabel_str)
% Create a mock quality plot when real data is not available
    
    x_data = 1:10;
    y_data = 0.7 + 0.2 * sin(x_data) + 0.1 * randn(size(x_data));
    
    plot(x_data, y_data, 'o-', 'LineWidth', 2, 'Color', [0.5, 0.5, 0.5]);
    xlabel(xlabel_str);
    ylabel(ylabel_str);
    title([title_str, ' (Mock Data)'], 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
    
    % Add note
    text(0.7, 0.2, 'Mock data', 'Units', 'normalized', ...
         'BackgroundColor', 'y', 'FontSize', 8);
end

function create_frequency_analysis(results)
% Create frequency domain analysis figure
    
    figure('Name', 'Module 1 Frequency Analysis', 'Position', [200, 200, 1400, 800]);
    
    if ~isfield(results, 'Sigma_tilde') || ~iscell(results.Sigma_tilde)
        % Show message that frequency analysis is not available
        subplot(1, 1, 1);
        axis off;
        text(0.5, 0.5, 'Frequency analysis data not available', ...
             'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');
        return;
    end
    
    Sigma_tilde = results.Sigma_tilde;
    F = length(Sigma_tilde);
    frequencies = 1:F;
    
    % Subplot 1: Spectral norms across frequencies
    subplot(2, 3, 1);
    spectral_norms = zeros(F, 1);
    for omega = 1:F
        spectral_norms(omega) = norm(Sigma_tilde{omega}, 2);
    end
    
    plot(frequencies, spectral_norms, 'b-', 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 6);
    title('Spectral Norms', 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Frequency Index');
    ylabel('Spectral Norm');
    grid on;
    
    % Subplot 2: Trace values
    subplot(2, 3, 2);
    trace_values = zeros(F, 1);
    for omega = 1:F
        trace_values(omega) = real(trace(Sigma_tilde{omega}));
    end
    
    plot(frequencies, trace_values, 'r-', 'LineWidth', 2, 'Marker', 's', 'MarkerSize', 6);
    title('Matrix Traces', 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Frequency Index');
    ylabel('Trace Value');
    grid on;
    
    % Add target line if known
    hold on;
    n = size(Sigma_tilde{1}, 1);
    plot([1, F], [n, n], 'g--', 'LineWidth', 2);
    legend('Actual', 'Target', 'Location', 'best');
    
    % Subplot 3: Complex component analysis
    subplot(2, 3, 3);
    complex_norms = zeros(F, 1);
    for omega = 1:F
        imag_part = imag(Sigma_tilde{omega});
        complex_norms(omega) = norm(imag_part, 'fro');
    end
    
    semilogy(frequencies, max(complex_norms, eps), 'g-', 'LineWidth', 2, 'Marker', '^', 'MarkerSize', 6);
    title('Complex Component Magnitude', 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Frequency Index');
    ylabel('Complex Norm (log scale)');
    grid on;
    
    % Subplot 4: Hermitian symmetry check
    subplot(2, 3, 4);
    hermitian_errors = zeros(F, 1);
    for omega = 1:F
        hermitian_diff = Sigma_tilde{omega} - Sigma_tilde{omega}';
        hermitian_errors(omega) = norm(hermitian_diff, 'fro') / norm(Sigma_tilde{omega}, 'fro');
    end
    
    semilogy(frequencies, max(hermitian_errors, eps), 'c-', 'LineWidth', 2, 'Marker', 'd', 'MarkerSize', 6);
    hold on;
    semilogy([1, F], [1e-12, 1e-12], 'r--', 'LineWidth', 2);
    
    title('Hermitian Symmetry Preservation', 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Frequency Index');
    ylabel('Relative Error (log scale)');
    legend('Error', 'Tolerance', 'Location', 'best');
    grid on;
    
    % Subplot 5: Connectivity analysis
    subplot(2, 3, 5);
    
    % Estimate connectivity from matrix sparsity
    sparsity_levels = zeros(F, 1);
    connectivity_estimates = zeros(F, 1);
    
    for omega = 1:F
        % Count significant off-diagonal elements
        Sigma = Sigma_tilde{omega};
        threshold = 1e-6 * norm(Sigma, 'fro');
        
        % Upper triangular part (excluding diagonal)
        n = size(Sigma, 1);
        upper_tri = triu(Sigma, 1);
        significant_elements = abs(upper_tri) > threshold;
        
        total_possible = n * (n - 1) / 2;
        sparsity_levels(omega) = 1 - sum(significant_elements(:)) / total_possible;
        connectivity_estimates(omega) = sum(significant_elements(:));
    end
    
    yyaxis left;
    plot(frequencies, connectivity_estimates, 'b-', 'LineWidth', 2, 'Marker', 'o');
    ylabel('Estimated Connections', 'Color', 'b');
    
    yyaxis right;
    plot(frequencies, sparsity_levels, 'r--', 'LineWidth', 2, 'Marker', 's');
    ylabel('Sparsity Level', 'Color', 'r');
    
    title('Connectivity Analysis', 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Frequency Index');
    grid on;
    
    % Subplot 6: Matrix statistics summary
    subplot(2, 3, 6);
    axis off;
    
    % Calculate summary statistics
    all_spectral_norms = spectral_norms;
    all_traces = trace_values;
    all_complex_norms = complex_norms;
    all_hermitian_errors = hermitian_errors;
    
    stats_text = {};
    stats_text{end+1} = 'Frequency Analysis Summary';
    stats_text{end+1} = '';
    stats_text{end+1} = sprintf('Frequencies: %d', F);
    stats_text{end+1} = sprintf('Matrix size: %d x %d', n, n);
    stats_text{end+1} = '';
    stats_text{end+1} = sprintf('Spectral norm range: [%.2e, %.2e]', min(all_spectral_norms), max(all_spectral_norms));
    stats_text{end+1} = sprintf('Trace range: [%.2e, %.2e]', min(all_traces), max(all_traces));
    stats_text{end+1} = sprintf('Max complex norm: %.2e', max(all_complex_norms));
    stats_text{end+1} = sprintf('Max Hermitian error: %.2e', max(all_hermitian_errors));
    
    y_pos = 0.9;
    for i = 1:length(stats_text)
        if i == 1
            text(0.1, y_pos, stats_text{i}, 'FontSize', 12, 'FontWeight', 'bold', 'Units', 'normalized');
        else
            text(0.1, y_pos, stats_text{i}, 'FontSize', 10, 'Units', 'normalized');
        end
        y_pos = y_pos - 0.08;
    end
end

function create_matrix_visualization(results)
% Create matrix visualization figure
    
    figure('Name', 'Module 1 Matrix Visualization', 'Position', [250, 250, 1400, 800]);
    
    if ~isfield(results, 'Sigma_tilde') || ~iscell(results.Sigma_tilde)
        subplot(1, 1, 1);
        axis off;
        text(0.5, 0.5, 'Matrix visualization data not available', ...
             'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');
        return;
    end
    
    Sigma_tilde = results.Sigma_tilde;
    F = length(Sigma_tilde);
    
    % Select representative frequencies for visualization
    freq_indices = unique(round(linspace(1, F, min(6, F))));
    
    for i = 1:length(freq_indices)
        subplot(2, 3, i);
        
        omega = freq_indices(i);
        matrix_to_plot = Sigma_tilde{omega};
        
        % Plot magnitude of matrix elements
        imagesc(log10(abs(matrix_to_plot) + eps));
        colormap(hot);
        colorbar;
        
        title(sprintf('Frequency %d (log|magnitude|)', omega), 'FontSize', 11, 'FontWeight', 'bold');
        xlabel('Node Index');
        ylabel('Node Index');
        axis square;
        
        % Add grid for small matrices
        n = size(matrix_to_plot, 1);
        if n <= 20
            hold on;
            for j = 0.5:1:(n+0.5)
                plot([j, j], [0.5, n+0.5], 'w-', 'LineWidth', 0.5);
                plot([0.5, n+0.5], [j, j], 'w-', 'LineWidth', 0.5);
            end
        end
        
        % Add statistics text
        max_val = max(abs(matrix_to_plot(:)));
        min_val = min(abs(matrix_to_plot(matrix_to_plot ~= 0)));
        if isempty(min_val)
            min_val = 0;
        end
        
        text(0.02, 0.98, sprintf('Max: %.2e\nMin: %.2e', max_val, min_val), ...
             'Units', 'normalized', 'VerticalAlignment', 'top', ...
             'BackgroundColor', 'w', 'FontSize', 8);
    end
    
    % If we have fewer than 6 frequencies, fill remaining subplots
    for i = (length(freq_indices)+1):6
        subplot(2, 3, i);
        axis off;
        
        if i == length(freq_indices)+1
            % Summary statistics
            text(0.1, 0.9, 'Matrix Statistics:', 'FontSize', 12, 'FontWeight', 'bold', 'Units', 'normalized');
            
            stats_text = {};
            stats_text{end+1} = sprintf('Total frequencies: %d', F);
            stats_text{end+1} = sprintf('Matrix size: %d × %d', size(Sigma_tilde{1}));
            
            % Calculate overall statistics
            all_elements = [];
            for f = 1:F
                matrix_elements = Sigma_tilde{f}(:);
                all_elements = [all_elements; matrix_elements];
            end
            
            stats_text{end+1} = sprintf('Overall max: %.2e', max(abs(all_elements)));
            stats_text{end+1} = sprintf('Overall mean: %.2e', mean(abs(all_elements)));
            stats_text{end+1} = sprintf('Nonzero ratio: %.2f%%', 100 * mean(abs(all_elements) > 1e-12));
            
            y_pos = 0.7;
            for j = 1:length(stats_text)
                text(0.1, y_pos, stats_text{j}, 'FontSize', 10, 'Units', 'normalized');
                y_pos = y_pos - 0.1;
            end
            
        else
            % Empty placeholder
            text(0.5, 0.5, sprintf('Frequency slot %d\n(not used)', i), ...
                 'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', [0.7, 0.7, 0.7]);
        end
    end
end