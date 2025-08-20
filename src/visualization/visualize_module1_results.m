function visualize_module1_results(demo_results)
% VISUALIZE_MODULE1_RESULTS - Comprehensive visualization of Module 1 preprocessing results
%
% This function creates multiple figures to show different aspects of preprocessing,
% including quality assessment visualization
%
% Input:
%   demo_results - Output structure from demo_module1_preprocessing
%
% Output:
%   Creates 3-4 figure windows showing different aspects of results
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
        fprintf('Preprocessing failed, can only show basic information\n');
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

function is_valid = validate_demo_results(demo_results)
% Validate demo_results structure validity
    is_valid = false;
    
    required_fields = {'timestamp', 'data_generation', 'preprocessing', 'summary'};
    for i = 1:length(required_fields)
        if ~isfield(demo_results, required_fields{i})
            fprintf('Missing required field: %s\n', required_fields{i});
            return;
        end
    end
    
    is_valid = true;
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
                if isnumeric(field_value) && field_value > 0
                    timing_values(end+1) = field_value;
                    timing_labels{end+1} = strrep(timing_fields{i}, '_', ' ');
                end
            end
        end
        
        if ~isempty(timing_values)
            pie(timing_values, timing_labels);
            title('Processing Time Breakdown', 'FontSize', 12, 'FontWeight', 'bold');
        else
            text(0.5, 0.5, 'No timing data available', 'HorizontalAlignment', 'center');
            title('Processing Time Breakdown');
        end
    else
        text(0.5, 0.5, 'Timing information not available', 'HorizontalAlignment', 'center');
        title('Processing Time Breakdown');
    end
    
    % Subplot 2: Data dimension information (FIXED - safe field access)
    subplot(2, 3, 2);
    if demo_results.data_generation.success
        params = demo_results.data_generation.params;
        
        % Safe field access with proper defaults - FIXED
        n_nodes = get_safe_field(params, 'n_nodes', 'Unknown');
        n_freq = get_safe_field(params, 'n_freq', 'Unknown');
        n_samples = get_safe_field(params, 'n_samples', 'Unknown');
        graph_type = get_safe_field(params, 'graph_type', 'Unknown');
        
        % Safe field access for edges - FIXED
        if isfield(params, 'n_edges')
            n_edges = params.n_edges;
            edges_text = sprintf('Edges: %d', n_edges);
            
            % Calculate sparsity if possible
            if isnumeric(n_nodes) && n_nodes > 1
                max_edges = n_nodes * (n_nodes-1) / 2;
                sparsity = (1 - n_edges/max_edges) * 100;
                sparsity_text = sprintf('Sparsity: %.1f%%', sparsity);
            else
                sparsity_text = 'Sparsity: Unknown';
            end
        else
            edges_text = 'Edges: Not available';
            sparsity_text = 'Sparsity: Not available';
        end
        
        % Create information text
        info_text = {
            sprintf('Nodes: %s', convert_to_string(n_nodes)),
            sprintf('Frequencies: %s', convert_to_string(n_freq)),
            sprintf('Samples: %s', convert_to_string(n_samples)),
            sprintf('Graph type: %s', convert_to_string(graph_type)),
            edges_text,
            sparsity_text
        };
        
        text(0.1, 0.9, info_text, 'VerticalAlignment', 'top', 'FontSize', 10);
        axis off;
        title('Data Characteristics', 'FontSize', 12, 'FontWeight', 'bold');
    else
        text(0.5, 0.5, 'Data generation failed', 'HorizontalAlignment', 'center');
        title('Data Characteristics');
    end
    
    % Subplot 3: Processing status
    subplot(2, 3, 3);
    create_status_indicator(demo_results);
    
    % Subplot 4: Diagonal smoothing effects
    subplot(2, 3, 4);
    plot_diagonal_smoothing(results);
    
    % Subplot 5: Condition number improvement
    subplot(2, 3, 5);
    plot_condition_improvement(results);
    
    % Subplot 6: Memory usage estimation
    subplot(2, 3, 6);
    plot_memory_usage(results);
end

function create_status_indicator(demo_results)
% Create visual status indicator
    
    % Status indicators
    data_gen_ok = demo_results.data_generation.success;
    preprocessing_ok = demo_results.preprocessing.success;
    
    status_text = {
        sprintf('Data Generation: %s', get_status_text(data_gen_ok)),
        sprintf('Preprocessing: %s', get_status_text(preprocessing_ok)),
        '',
        sprintf('Overall: %s', get_status_text(data_gen_ok && preprocessing_ok))
    };
    
    text(0.1, 0.9, status_text, 'VerticalAlignment', 'top', 'FontSize', 12);
    axis off;
    title('Processing Status', 'FontSize', 12, 'FontWeight', 'bold');
end

function status_str = get_status_text(is_success)
% Get status text with emoji
    if is_success
        status_str = '✓ SUCCESS';
    else
        status_str = '✗ FAILED';
    end
end

function plot_diagonal_smoothing(results)
% Plot diagonal smoothing effects with REAL implementation
    
    % Try to get actual data for smoothing visualization
    if isfield(results, 'Sigma_tilde') && iscell(results.Sigma_tilde)
        matrices = results.Sigma_tilde;
        n_freq = length(matrices);
        n_nodes = size(matrices{1}, 1);
        
        % Extract diagonal elements
        diag_elements = zeros(n_freq, n_nodes);
        for f = 1:n_freq
            diag_elements(f, :) = real(diag(matrices{f}));
        end
        
        % Plot diagonal power evolution
        freq_indices = 1:n_freq;
        plot(freq_indices, diag_elements(:, 1:min(3, n_nodes)), 'LineWidth', 2);
        xlabel('Frequency Index');
        ylabel('Diagonal Elements');
        title('Diagonal Element Evolution');
        legend(arrayfun(@(x) sprintf('Node %d', x), 1:min(3, n_nodes), 'UniformOutput', false));
        grid on;
    else
        text(0.5, 0.5, 'Matrix data not available for smoothing analysis', 'HorizontalAlignment', 'center');
        title('Diagonal Smoothing Effects');
        axis off;
    end
end

function plot_condition_improvement(results)
% Plot condition number improvement with REAL implementation
    
    if isfield(results, 'Sigma_tilde') && iscell(results.Sigma_tilde)
        matrices = results.Sigma_tilde;
        n_freq = length(matrices);
        
        % Calculate condition numbers
        condition_numbers = zeros(n_freq, 1);
        for f = 1:n_freq
            try
                eigvals = eig(matrices{f});
                eigvals_real = real(eigvals(abs(eigvals) > 1e-12));
                if length(eigvals_real) > 1
                    condition_numbers(f) = max(eigvals_real) / min(eigvals_real);
                else
                    condition_numbers(f) = 1;
                end
            catch
                condition_numbers(f) = NaN;
            end
        end
        
        % Plot condition numbers
        semilogy(1:n_freq, condition_numbers, 'b-o', 'LineWidth', 2, 'MarkerSize', 6);
        xlabel('Frequency Index');
        ylabel('Condition Number (log scale)');
        title('Matrix Condition Numbers');
        grid on;
        
        % Add reference lines
        hold on;
        ylims = ylim;
        plot([1, n_freq], [100, 100], 'g--', 'LineWidth', 1);
        plot([1, n_freq], [1000, 1000], 'y--', 'LineWidth', 1);
        plot([1, n_freq], [10000, 10000], 'r--', 'LineWidth', 1);
        legend('Condition Numbers', 'Good (100)', 'Fair (1K)', 'Poor (10K)', 'Location', 'best');
        hold off;
    else
        text(0.5, 0.5, 'Matrix data not available for condition analysis', 'HorizontalAlignment', 'center');
        title('Condition Number Improvement');
        axis off;
    end
end

function plot_memory_usage(results)
% Plot memory usage estimation with REAL implementation
    
    if isfield(results, 'Sigma_tilde') && iscell(results.Sigma_tilde)
        matrices = results.Sigma_tilde;
        n_matrices = length(matrices);
        matrix_size = size(matrices{1});
        
        % Calculate memory usage
        elements_per_matrix = prod(matrix_size);
        total_elements = elements_per_matrix * n_matrices;
        
        % Check if complex
        is_complex = ~isreal(matrices{1});
        bytes_per_element = 16; % Complex double
        if ~is_complex
            bytes_per_element = 8; % Real double
        end
        
        memory_bytes = total_elements * bytes_per_element;
        memory_mb = memory_bytes / (1024^2);
        
        % Create bar chart of memory usage breakdown
        categories = {'Matrices', 'Processing', 'Temporary'};
        memory_breakdown = [memory_mb, memory_mb * 0.3, memory_mb * 0.2];
        
        bar(memory_breakdown, 'FaceColor', [0.2, 0.6, 0.8]);
        set(gca, 'XTickLabel', categories);
        ylabel('Memory Usage (MB)');
        title(sprintf('Memory Usage Estimation\nTotal: %.1f MB', sum(memory_breakdown)));
        grid on;
        
        % Add text annotations
        for i = 1:length(memory_breakdown)
            text(i, memory_breakdown(i) + max(memory_breakdown)*0.02, ...
                 sprintf('%.1f MB', memory_breakdown(i)), ...
                 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        end
    else
        text(0.5, 0.5, 'Memory estimation unavailable', 'HorizontalAlignment', 'center');
        title('Memory Usage');
        axis off;
    end
end

function create_quality_dashboard(results)
% Create quality assessment dashboard
    
    figure('Name', 'Quality Dashboard', 'Position', [200, 200, 1200, 800]);
    
    % Subplot 1: Overall quality gauge
    subplot(2, 3, 1);
    create_quality_gauge(results);
    
    % Subplot 2: Diagonal error distribution
    subplot(2, 3, 2);
    create_diagonal_analysis(results);
    
    % Subplot 3: Whitening effectiveness
    subplot(2, 3, 3);
    create_effectiveness_plot(results);
    
    % Subplot 4: Condition number analysis
    subplot(2, 3, 4);
    create_condition_analysis(results);
    
    % Subplot 5: Error trends across frequencies
    subplot(2, 3, 5);
    create_error_trends(results);
    
    % Subplot 6: Quality summary
    subplot(2, 3, 6);
    create_quality_summary_text(results);
end

function create_quality_gauge(results)
% Create overall quality gauge
    
    overall_score = 0.5; % Default value
    
    if isfield(results, 'processing_stats') && ...
       isfield(results.processing_stats, 'whitening_quality')
        wq = results.processing_stats.whitening_quality;
        
        if isfield(wq, 'overall_score')
            overall_score = wq.overall_score / 100;  % Convert to 0-1 range
        elseif isfield(wq, 'whitening_effectiveness')
            effectiveness = wq.whitening_effectiveness;
            if isnumeric(effectiveness) && ~isempty(effectiveness)
                overall_score = mean(effectiveness(~isnan(effectiveness)));
            end
        end
    end
    
    % Ensure score is in valid range
    overall_score = max(0, min(1, overall_score));
    
    % Create simple gauge using bar
    theta = linspace(0, pi, 100);
    r_outer = 1;
    r_inner = 0.7;
    
    % Background arc
    fill([r_inner*cos(theta), fliplr(r_outer*cos(theta))], ...
         [r_inner*sin(theta), fliplr(r_outer*sin(theta))], ...
         [0.9 0.9 0.9], 'EdgeColor', 'none');
    hold on;
    
    % Quality arc
    theta_score = pi * (1 - overall_score);
    theta_fill = linspace(theta_score, pi, 50);
    if ~isempty(theta_fill)
        color = [1-overall_score, overall_score, 0];
        fill([r_inner*cos(theta_fill), fliplr(r_outer*cos(theta_fill))], ...
             [r_inner*sin(theta_fill), fliplr(r_outer*sin(theta_fill))], ...
             color, 'EdgeColor', 'none');
    end
    
    % Needle
    needle_angle = pi * (1 - overall_score);
    plot([0, 0.9*cos(needle_angle)], [0, 0.9*sin(needle_angle)], ...
         'k-', 'LineWidth', 3);
    
    % Labels
    text(0, -0.3, sprintf('%.2f', overall_score), ...
         'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
    text(-1, 0, 'Poor', 'HorizontalAlignment', 'center');
    text(1, 0, 'Excellent', 'HorizontalAlignment', 'center');
    
    axis equal;
    axis off;
    title('Overall Quality Score');
end

function create_diagonal_analysis(results)
% Create diagonal error analysis with real implementation
    
    if isfield(results, 'processing_stats') && ...
       isfield(results.processing_stats, 'whitening_quality') && ...
       isfield(results.processing_stats.whitening_quality, 'diagonal_errors')
        
        errors = results.processing_stats.whitening_quality.diagonal_errors;
        errors = errors(~isnan(errors) & isfinite(errors));
        
        if ~isempty(errors)
            % Create histogram of diagonal errors
            histogram(errors, min(20, length(errors)));
            xlabel('Diagonal Error');
            ylabel('Frequency');
            title(sprintf('Diagonal Errors Distribution\n(Mean: %.3f, Max: %.3f)', mean(errors), max(errors)));
            grid on;
            
            % Add vertical lines for thresholds
            hold on;
            ylims = ylim;
            plot([0.05, 0.05], ylims, 'g--', 'LineWidth', 2); % Good threshold
            plot([0.10, 0.10], ylims, 'y--', 'LineWidth', 2); % Acceptable threshold
            plot([0.20, 0.20], ylims, 'r--', 'LineWidth', 2); % Poor threshold
            legend('Errors', 'Good (5%)', 'Acceptable (10%)', 'Poor (20%)', 'Location', 'best');
            hold off;
        else
            text(0.5, 0.5, 'No valid diagonal error data', 'HorizontalAlignment', 'center');
            title('Diagonal Error Analysis');
        end
    else
        text(0.5, 0.5, 'Diagonal error data not available', 'HorizontalAlignment', 'center');
        title('Diagonal Error Analysis');
    end
    axis off;
end

function create_effectiveness_plot(results)
% Create whitening effectiveness plot with real implementation
    
    if isfield(results, 'processing_stats') && ...
       isfield(results.processing_stats, 'whitening_quality')
        
        wq = results.processing_stats.whitening_quality;
        
        % Extract effectiveness data
        if isfield(wq, 'effectiveness') && ~isempty(wq.effectiveness)
            effectiveness = wq.effectiveness(~isnan(wq.effectiveness) & isfinite(wq.effectiveness));
            
            if ~isempty(effectiveness)
                freq_indices = 1:length(effectiveness);
                
                % Plot effectiveness vs frequency
                plot(freq_indices, effectiveness, 'bo-', 'LineWidth', 2, 'MarkerSize', 6);
                xlabel('Frequency Index');
                ylabel('Whitening Effectiveness');
                title('Whitening Effectiveness by Frequency');
                grid on;
                ylim([0, 1.1]);
                
                % Add horizontal reference lines
                hold on;
                plot(xlim, [0.9, 0.9], 'g--', 'LineWidth', 1); % Excellent
                plot(xlim, [0.8, 0.8], 'y--', 'LineWidth', 1); % Good
                plot(xlim, [0.7, 0.7], 'r--', 'LineWidth', 1); % Acceptable
                legend('Effectiveness', 'Excellent (90%)', 'Good (80%)', 'Acceptable (70%)', 'Location', 'best');
                hold off;
                
                return;
            end
        end
    end
    
    % Fallback display
    text(0.5, 0.5, {'Whitening Effectiveness', '', 'Data not available or invalid'}, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    title('Whitening Effectiveness');
    axis off;
end

function create_condition_analysis(results)
% Create condition number analysis with real implementation
    
    if isfield(results, 'processing_stats') && ...
       isfield(results.processing_stats, 'whitening_quality') && ...
       isfield(results.processing_stats.whitening_quality, 'condition_numbers')
        
        cond_nums = results.processing_stats.whitening_quality.condition_numbers;
        cond_nums = cond_nums(~isnan(cond_nums) & isfinite(cond_nums));
        
        if ~isempty(cond_nums)
            freq_indices = 1:length(cond_nums);
            
            % Plot condition numbers (log scale)
            semilogy(freq_indices, cond_nums, 'ro-', 'LineWidth', 2, 'MarkerSize', 6);
            xlabel('Frequency Index');
            ylabel('Condition Number (log scale)');
            title('Condition Numbers by Frequency');
            grid on;
            
            % Add horizontal reference lines
            hold on;
            plot(xlim, [1e2, 1e2], 'g--', 'LineWidth', 1); % Good
            plot(xlim, [1e4, 1e4], 'y--', 'LineWidth', 1); % Acceptable
            plot(xlim, [1e6, 1e6], 'r--', 'LineWidth', 1); % Poor
            legend('Condition Numbers', 'Good (<100)', 'Acceptable (<10K)', 'Poor (<1M)', 'Location', 'best');
            hold off;
            
            return;
        end
    end
    
    % Fallback display
    text(0.5, 0.5, {'Condition Number Analysis', '', 'Data not available or invalid'}, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    title('Condition Number Analysis');
    axis off;
end

function create_error_trends(results)
% Create error trends across frequencies with real implementation
    
    if isfield(results, 'processing_stats') && ...
       isfield(results.processing_stats, 'whitening_quality')
        
        wq = results.processing_stats.whitening_quality;
        
        % Extract error data
        has_diagonal = isfield(wq, 'diagonal_errors') && ~isempty(wq.diagonal_errors);
        has_hermitian = isfield(wq, 'hermitian_errors') && ~isempty(wq.hermitian_errors);
        
        if has_diagonal || has_hermitian
            hold on;
            legends = {};
            
            if has_diagonal
                diag_errors = wq.diagonal_errors(~isnan(wq.diagonal_errors) & isfinite(wq.diagonal_errors));
                if ~isempty(diag_errors)
                    freq_indices = 1:length(diag_errors);
                    semilogy(freq_indices, diag_errors, 'b-o', 'LineWidth', 2, 'MarkerSize', 4);
                    legends{end+1} = 'Diagonal Errors';
                end
            end
            
            if has_hermitian
                herm_errors = wq.hermitian_errors(~isnan(wq.hermitian_errors) & isfinite(wq.hermitian_errors));
                if ~isempty(herm_errors) && any(herm_errors > 0)
                    freq_indices = 1:length(herm_errors);
                    semilogy(freq_indices, herm_errors, 'r-s', 'LineWidth', 2, 'MarkerSize', 4);
                    legends{end+1} = 'Hermitian Errors';
                end
            end
            
            xlabel('Frequency Index');
            ylabel('Error Magnitude (log scale)');
            title('Error Evolution Across Frequencies');
            grid on;
            
            if ~isempty(legends)
                legend(legends, 'Location', 'best');
            end
            hold off;
            
            return;
        end
    end
    
    % Fallback display
    text(0.5, 0.5, {'Error Trends Analysis', '', 'Error data not available'}, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    title('Error Trends');
    axis off;
end

function create_quality_summary_text(results)
% Create quality summary text
    
    summary_text = {'Quality Summary:'};
    
    if isfield(results, 'processing_stats') && ...
       isfield(results.processing_stats, 'whitening_quality')
        
        wq = results.processing_stats.whitening_quality;
        
        if isfield(wq, 'diagonal_errors') && ~isempty(wq.diagonal_errors)
            errors = wq.diagonal_errors(~isnan(wq.diagonal_errors));
            if ~isempty(errors)
                summary_text{end+1} = sprintf('Mean diagonal error: %.3f', mean(errors));
                summary_text{end+1} = sprintf('Max diagonal error: %.3f', max(errors));
            end
        end
        
        if isfield(wq, 'success_rates')
            if isfield(wq.success_rates, 'good_diagonal')
                summary_text{end+1} = sprintf('Good diagonal rate: %.1f%%', ...
                                             wq.success_rates.good_diagonal * 100);
            end
        end
        
        if isfield(wq, 'overall_score')
            summary_text{end+1} = sprintf('Overall score: %.1f', wq.overall_score);
        end
    else
        summary_text{end+1} = 'Quality data not available';
    end
    
    text(0.1, 0.9, summary_text, 'VerticalAlignment', 'top', 'FontSize', 10);
    axis off;
    title('Quality Summary');
end

function create_frequency_analysis(results)
% Create frequency domain analysis
    
    figure('Name', 'Frequency Analysis', 'Position', [300, 300, 800, 600]);
    
    text(0.5, 0.5, {'Frequency Domain Analysis', '', 'Advanced frequency analysis', 'features coming soon'}, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 14);
    title('Frequency Domain Analysis');
    axis off;
end

function create_matrix_visualization(results)
% Create matrix visualization
    
    figure('Name', 'Matrix Visualization', 'Position', [400, 400, 800, 600]);
    
    text(0.5, 0.5, {'Matrix Visualization', '', 'Advanced matrix visualization', 'features coming soon'}, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 14);
    title('Matrix Visualization');
    axis off;
end

function show_failure_summary(demo_results)
% Show failure information
    
    fprintf('Processing failure summary:\n');
    if ~demo_results.data_generation.success
        fprintf('  Data generation failed\n');
    end
    if ~demo_results.preprocessing.success
        fprintf('  Preprocessing failed: %s\n', demo_results.preprocessing.error);
    end
end

function value = get_safe_field(struct_data, field_name, default_value)
% Safely get field value with default
    if isfield(struct_data, field_name)
        value = struct_data.(field_name);
    else
        value = default_value;
    end
end

function str_result = convert_to_string(input_value)
% Convert input to string safely
    if ischar(input_value) || isstring(input_value)
        str_result = char(input_value);
    elseif isnumeric(input_value)
        str_result = num2str(input_value);
    else
        str_result = 'Unknown';
    end
end