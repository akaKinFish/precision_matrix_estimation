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
    if demo_results.preprocessing.success
        results = demo_results.preprocessing.results;
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
                timing_values(end+1) = timing.(timing_fields{i});
                timing_labels{end+1} = strrep(timing_fields{i}, '_', ' ');
            end
        end
        
        if ~isempty(timing_values)
            pie(timing_values, timing_labels);
            title('Processing Time Breakdown', 'FontSize', 12, 'FontWeight', 'bold');
        end
    else
        text(0.5, 0.5, 'Timing information not available', 'HorizontalAlignment', 'center');
        title('Processing Time Breakdown');
    end
    
    % Subplot 2: Data dimension information
    subplot(2, 3, 2);
    if demo_results.data_generation.success
        params = demo_results.data_generation.params;
        
        % Create information text
        info_text = {
            sprintf('Nodes: %d', params.n_nodes),
            sprintf('Frequencies: %d', params.n_freq),
            sprintf('Samples: %d', params.n_samples),
            sprintf('Graph type: %s', params.graph_type),
            sprintf('Edges: %d', params.n_edges),
            sprintf('Sparsity: %.1f%%', (1 - params.n_edges/(params.n_nodes*(params.n_nodes-1)/2))*100)
        };
        
        text(0.1, 0.9, info_text, 'VerticalAlignment', 'top', 'FontSize', 10);
        axis off;
        title('Data Characteristics', 'FontSize', 12, 'FontWeight', 'bold');
    end
    
    % Subplot 3: Processing status
    subplot(2, 3, 3);
    create_status_indicator(demo_results);
    
    % Subplot 4: Diagonal smoothing effects
    subplot(2, 3, 4);
    if isfield(results, 'g_smooth') && iscell(results.g_smooth)
        plot_diagonal_smoothing(results);
    else
        text(0.5, 0.5, 'Diagonal data not available', 'HorizontalAlignment', 'center');
        title('Diagonal Smoothing Effects');
    end
    
    % Subplot 5: Condition number improvement
    subplot(2, 3, 5);
    if isfield(results, 'processing_stats')
        plot_condition_improvement(results);
    else
        text(0.5, 0.5, 'Condition number data not available', 'HorizontalAlignment', 'center');
        title('Condition Number Improvement');
    end
    
    % Subplot 6: Memory usage estimation
    subplot(2, 3, 6);
    if isfield(results, 'Sigma_tilde') && iscell(results.Sigma_tilde)
        plot_memory_usage(results);
    else
        text(0.5, 0.5, 'Memory info not available', 'HorizontalAlignment', 'center');
        title('Memory Usage');
    end
end

function create_quality_dashboard(results)
% Create quality assessment dashboard
    
    figure('Name', 'Module 1 Quality Assessment Dashboard', 'Position', [150, 150, 1200, 900]);
    
    % Extract quality metrics
    if isfield(results, 'processing_stats') && ...
       isfield(results.processing_stats, 'whitening_quality')
        wq = results.processing_stats.whitening_quality;
    else
        fprintf('Warning: Whitening quality data not available\n');
        wq = struct();
    end
    
    % Subplot 1: Diagonal normalization quality
    subplot(3, 3, 1);
    plot_diagonal_quality(wq);
    
    % Subplot 2: Whitening effectiveness distribution
    subplot(3, 3, 2);
    plot_whitening_effectiveness(wq);
    
    % Subplot 3: Condition number distribution
    subplot(3, 3, 3);
    plot_condition_distribution(wq);
    
    % Subplot 4: Quality score gauge
    subplot(3, 3, 4);
    plot_quality_gauge(wq);
    
    % Subplot 5: Success rate bar chart
    subplot(3, 3, 5);
    plot_success_rates(wq);
    
    % Subplot 6: Error distribution heatmap
    subplot(3, 3, 6);
    plot_error_heatmap(wq);
    
    % Subplot 7: Numerical stability indicators
    subplot(3, 3, 7);
    plot_numerical_stability(wq);
    
    % Subplot 8: Quality trend (across frequencies)
    subplot(3, 3, 8);
    plot_quality_trend(wq);
    
    % Subplot 9: Problem diagnosis
    subplot(3, 3, 9);
    plot_problem_diagnosis(wq);
end

function create_frequency_analysis(results)
% Create frequency domain analysis figure
    
    figure('Name', 'Module 1 Frequency Domain Analysis', 'Position', [200, 200, 1300, 700]);
    
    if ~isfield(results, 'Sigma_tilde') || ~iscell(results.Sigma_tilde)
        text(0.5, 0.5, 'Whitened data not available', 'HorizontalAlignment', 'center');
        return;
    end
    
    Sigma_tilde = results.Sigma_tilde;
    F = length(Sigma_tilde);
    
    % Compute frequency domain metrics
    diagonal_means = zeros(F, 1);
    frobenius_norms = zeros(F, 1);
    condition_numbers = zeros(F, 1);
    eigenvalue_spreads = zeros(F, 1);
    
    for f = 1:F
        S = Sigma_tilde{f};
        diagonal_means(f) = mean(real(diag(S)));
        frobenius_norms(f) = norm(S, 'fro');
        
        % Compute condition number and eigenvalue distribution
        try
            eigs_S = eig(S);
            eigs_S = real(eigs_S(eigs_S > 1e-12));
            if length(eigs_S) > 1
                condition_numbers(f) = max(eigs_S) / min(eigs_S);
                eigenvalue_spreads(f) = std(eigs_S) / mean(eigs_S);
            else
                condition_numbers(f) = 1;
                eigenvalue_spreads(f) = 0;
            end
        catch
            condition_numbers(f) = NaN;
            eigenvalue_spreads(f) = NaN;
        end
    end
    
    % Subplot 1: Diagonal mean trend
    subplot(2, 3, 1);
    plot(1:F, diagonal_means, 'b-o', 'LineWidth', 2, 'MarkerSize', 6);
    hold on;
    yline(1.0, 'r--', 'LineWidth', 2, 'DisplayName', 'Target');
    xlabel('Frequency Index');
    ylabel('Diagonal Mean');
    title('Post-Whitening Diagonal Mean');
    grid on;
    legend;
    
    % Subplot 2: Frobenius norm
    subplot(2, 3, 2);
    plot(1:F, frobenius_norms, 'g-s', 'LineWidth', 2, 'MarkerSize', 6);
    xlabel('Frequency Index');
    ylabel('Frobenius Norm');
    title('Covariance Matrix Norms');
    grid on;
    
    % Subplot 3: Condition numbers
    subplot(2, 3, 3);
    semilogy(1:F, condition_numbers, 'r-^', 'LineWidth', 2, 'MarkerSize', 6);
    xlabel('Frequency Index');
    ylabel('Condition Number (log)');
    title('Matrix Condition Numbers');
    grid on;
    
    % Subplot 4: Eigenvalue distribution heatmap
    subplot(2, 3, 4);
    eigenvalue_matrix = zeros(F, size(Sigma_tilde{1}, 1));
    for f = 1:F
        try
            eigs_f = sort(real(eig(Sigma_tilde{f})), 'descend');
            eigenvalue_matrix(f, :) = eigs_f;
        catch
            eigenvalue_matrix(f, :) = NaN;
        end
    end
    imagesc(eigenvalue_matrix);
    colorbar;
    xlabel('Eigenvalue Index');
    ylabel('Frequency Index');
    title('Eigenvalue Heatmap');
    
    % Subplot 5: Eigenvalue coefficient of variation
    subplot(2, 3, 5);
    plot(1:F, eigenvalue_spreads, 'm-d', 'LineWidth', 2, 'MarkerSize', 6);
    xlabel('Frequency Index');
    ylabel('Coefficient of Variation');
    title('Eigenvalue Distribution Variation');
    grid on;
    
    % Subplot 6: Quality assessment radar chart
    subplot(2, 3, 6);
    create_quality_radar(diagonal_means, condition_numbers, frobenius_norms);
end

function create_matrix_visualization(results)
% Create matrix visualization figure
    
    figure('Name', 'Module 1 Matrix Visualization', 'Position', [250, 250, 1400, 600]);
    
    if ~isfield(results, 'Sigma_tilde') || ~iscell(results.Sigma_tilde)
        text(0.5, 0.5, 'Matrix data not available', 'HorizontalAlignment', 'center');
        return;
    end
    
    Sigma_tilde = results.Sigma_tilde;
    F = length(Sigma_tilde);
    
    % Select frequencies to display
    freq_indices = round(linspace(1, F, min(6, F)));
    
    for i = 1:length(freq_indices)
        f = freq_indices(i);
        
        % Real part
        subplot(2, length(freq_indices), i);
        imagesc(real(Sigma_tilde{f}));
        colorbar;
        title(sprintf('Real Part - Freq %d', f));
        axis square;
        
        % Imaginary part (if exists)
        subplot(2, length(freq_indices), i + length(freq_indices));
        if any(any(imag(Sigma_tilde{f}) ~= 0))
            imagesc(imag(Sigma_tilde{f}));
            title(sprintf('Imaginary Part - Freq %d', f));
        else
            imagesc(abs(Sigma_tilde{f}));
            title(sprintf('Magnitude - Freq %d', f));
        end
        colorbar;
        axis square;
    end
end

%% Helper plotting functions

function create_status_indicator(demo_results)
% Create status indicator
    
    % Status information
    data_status = demo_results.data_generation.success;
    proc_status = demo_results.preprocessing.success;
    overall_status = strcmp(demo_results.summary.overall_status, 'SUCCESS');
    
    statuses = [data_status, proc_status, overall_status];
    labels = {'Data Generation', 'Preprocessing', 'Overall'};
    colors = [0.2 0.8 0.2; 0.8 0.2 0.2]; % Green for success, red for failure
    
    for i = 1:3
        if statuses(i)
            rectangle('Position', [i-0.4, 0.4, 0.8, 0.2], ...
                     'FaceColor', colors(1,:), 'EdgeColor', 'k');
            text(i, 0.5, 'SUCCESS', 'HorizontalAlignment', 'center', ...
                 'FontWeight', 'bold', 'Color', 'w');
        else
            rectangle('Position', [i-0.4, 0.4, 0.8, 0.2], ...
                     'FaceColor', colors(2,:), 'EdgeColor', 'k');
            text(i, 0.5, 'FAILED', 'HorizontalAlignment', 'center', ...
                 'FontWeight', 'bold', 'Color', 'w');
        end
        text(i, 0.2, labels{i}, 'HorizontalAlignment', 'center', ...
             'FontSize', 10);
    end
    
    xlim([0.5, 3.5]);
    ylim([0, 1]);
    axis off;
    title('Processing Status', 'FontSize', 12, 'FontWeight', 'bold');
end

function plot_diagonal_smoothing(results)
% Plot diagonal smoothing effects
    
    if isfield(results, 'Sigma_emp') && isfield(results, 'g_smooth')
        F = length(results.g_smooth);
        
        % Extract original and smoothed diagonals
        orig_diag = zeros(F, size(results.Sigma_emp{1}, 1));
        smooth_diag = zeros(F, size(results.g_smooth{1}, 1));
        
        for f = 1:F
            orig_diag(f, :) = real(diag(results.Sigma_emp{f}));
            smooth_diag(f, :) = results.g_smooth{f};
        end
        
        % Plot diagonal evolution for a few nodes
        n_nodes_show = min(3, size(orig_diag, 2));
        for i = 1:n_nodes_show
            plot(1:F, orig_diag(:, i), '--', 'LineWidth', 1.5, ...
                 'DisplayName', sprintf('Original-Node%d', i));
            hold on;
            plot(1:F, smooth_diag(:, i), '-', 'LineWidth', 2, ...
                 'DisplayName', sprintf('Smoothed-Node%d', i));
        end
        
        xlabel('Frequency Index');
        ylabel('Diagonal Power');
        title('Diagonal Smoothing Effects');
        legend('Location', 'best');
        grid on;
    else
        text(0.5, 0.5, 'Diagonal data missing', 'HorizontalAlignment', 'center');
        title('Diagonal Smoothing Effects');
    end
end

function plot_diagonal_quality(wq)
% Plot diagonal quality metrics
    
    if isfield(wq, 'diagonal_normalization')
        dn = wq.diagonal_normalization;
        
        if isfield(dn, 'max_errors') && ~isempty(dn.max_errors)
            histogram(dn.max_errors, 10, 'FaceColor', [0.7 0.7 1]);
            xlabel('Maximum Diagonal Error');
            ylabel('Frequency Count');
            title('Diagonal Normalization Error Distribution');
            
            % Add threshold lines
            hold on;
            xline(0.05, 'r--', '5% Threshold', 'LineWidth', 2);
            xline(0.1, 'y--', '10% Threshold', 'LineWidth', 2);
        else
            text(0.5, 0.5, 'Diagonal error data not available', 'HorizontalAlignment', 'center');
        end
    else
        text(0.5, 0.5, 'Diagonal normalization data not available', 'HorizontalAlignment', 'center');
    end
    grid on;
    title('Diagonal Quality');
end

function plot_quality_gauge(wq)
% Create quality score gauge
    
    % Try to get overall quality score from different fields
    overall_score = 0;
    if isfield(wq, 'overall_score')
        overall_score = wq.overall_score;
    elseif isfield(wq, 'whitening_effectiveness')
        effectiveness = wq.whitening_effectiveness;
        if isnumeric(effectiveness) && ~isempty(effectiveness)
            overall_score = mean(effectiveness(~isnan(effectiveness)));
        end
    end
    
    % Create simple gauge
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
        fill([r_inner*cos(theta_fill), fliplr(r_outer*cos(theta_fill))], ...
             [r_inner*sin(theta_fill), fliplr(r_outer*sin(theta_fill))], ...
             [1-overall_score, overall_score, 0], 'EdgeColor', 'none');
    end
    
    % Needle
    needle_angle = pi * (1 - overall_score);
    plot([0, 0.9*cos(needle_angle)], [0, 0.9*sin(needle_angle)], ...
         'k-', 'LineWidth', 3);
    
    % Labels
    text(0, -0.3, sprintf('%.3f', overall_score), ...
         'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
    text(-1, 0, 'Poor', 'HorizontalAlignment', 'center');
    text(1, 0, 'Excellent', 'HorizontalAlignment', 'center');
    
    axis equal;
    axis off;
    title('Overall Quality Score');
end

function create_quality_radar(diagonal_means, condition_numbers, frobenius_norms)
% Create quality assessment radar chart
    
    % Normalize metrics (0-1, 1 is best)
    diag_score = 1 - abs(mean(diagonal_means) - 1.0);  % Closer to 1 is better
    valid_conds = condition_numbers(~isnan(condition_numbers));
    if ~isempty(valid_conds)
        cond_score = 1 / (1 + log10(mean(valid_conds)));  % Lower condition number is better
    else
        cond_score = 0;
    end
    norm_score = 1 - (std(frobenius_norms) / mean(frobenius_norms));  % Lower norm variation is better
    
    scores = [diag_score, cond_score, norm_score];
    labels = {'Diagonal Quality', 'Condition Number', 'Norm Stability'};
    
    % Create simplified radar chart using plot instead of polar
    angles = [0, 2*pi/3, 4*pi/3, 0];  % Closed triangle
    r = [scores, scores(1)];  % Closed
    
    % Convert to cartesian coordinates
    x = r .* cos(angles);
    y = r .* sin(angles);
    
    plot(x, y, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
    hold on;
    
    % Reference circle
    theta = linspace(0, 2*pi, 100);
    plot(cos(theta), sin(theta), 'r--', 'LineWidth', 1);
    
    % Add labels
    label_r = 1.2;
    text(label_r * cos(0), label_r * sin(0), labels{1}, 'HorizontalAlignment', 'center');
    text(label_r * cos(2*pi/3), label_r * sin(2*pi/3), labels{2}, 'HorizontalAlignment', 'center');
    text(label_r * cos(4*pi/3), label_r * sin(4*pi/3), labels{3}, 'HorizontalAlignment', 'center');
    
    axis equal;
    xlim([-1.5, 1.5]);
    ylim([-1.5, 1.5]);
    title('Quality Radar Chart');
    legend({'Actual', 'Reference'}, 'Location', 'best');
    grid on;
end

function plot_whitening_effectiveness(wq)
% Plot whitening effectiveness
    
    % Create placeholder with informative message
    text(0.5, 0.5, {'Whitening Effectiveness Analysis', '', ...
                   'This plot would show:', ...
                   '- Effectiveness scores per frequency', ...
                   '- Distribution of quality metrics', ...
                   '- Threshold comparisons'}, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    xlim([0, 1]);
    ylim([0, 1]);
    title('Whitening Effectiveness');
    axis off;
end

function plot_condition_distribution(wq)
% Plot condition number distribution
    
    text(0.5, 0.5, {'Condition Number Distribution', '', ...
                   'This plot would show:', ...
                   '- Histogram of condition numbers', ...
                   '- Statistical summary', ...
                   '- Stability indicators'}, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    xlim([0, 1]);
    ylim([0, 1]);
    title('Condition Number Distribution');
    axis off;
end

function plot_success_rates(wq)
% Plot success rates
    
    text(0.5, 0.5, {'Success Rate Statistics', '', ...
                   'This plot would show:', ...
                   '- Pass/fail rates by threshold', ...
                   '- Quality level distribution', ...
                   '- Performance metrics'}, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    xlim([0, 1]);
    ylim([0, 1]);
    title('Success Rates');
    axis off;
end

function plot_error_heatmap(wq)
% Plot error distribution heatmap
    
    text(0.5, 0.5, {'Error Distribution Heatmap', '', ...
                   'This plot would show:', ...
                   '- Error patterns across frequencies', ...
                   '- Spatial error distribution', ...
                   '- Problem area identification'}, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    xlim([0, 1]);
    ylim([0, 1]);
    title('Error Distribution');
    axis off;
end

function plot_numerical_stability(wq)
% Plot numerical stability indicators
    
    text(0.5, 0.5, {'Numerical Stability Indicators', '', ...
                   'This plot would show:', ...
                   '- Eigenvalue analysis', ...
                   '- Conditioning metrics', ...
                   '- Stability warnings'}, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    xlim([0, 1]);
    ylim([0, 1]);
    title('Numerical Stability');
    axis off;
end

function plot_quality_trend(wq)
% Plot quality trend analysis
    
    text(0.5, 0.5, {'Quality Trend Analysis', '', ...
                   'This plot would show:', ...
                   '- Quality evolution across frequencies', ...
                   '- Trend identification', ...
                   '- Pattern analysis'}, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    xlim([0, 1]);
    ylim([0, 1]);
    title('Quality Trends');
    axis off;
end

function plot_problem_diagnosis(wq)
% Plot problem diagnosis
    
    text(0.5, 0.5, {'Problem Diagnosis', '', ...
                   'This plot would show:', ...
                   '- Issue categorization', ...
                   '- Root cause analysis', ...
                   '- Fix recommendations'}, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    xlim([0, 1]);
    ylim([0, 1]);
    title('Problem Diagnosis');
    axis off;
end

function plot_condition_improvement(results)
% Plot condition number improvement
    
    text(0.5, 0.5, {'Condition Number Improvement', '', ...
                   'This plot would show:', ...
                   '- Before vs after comparison', ...
                   '- Improvement statistics', ...
                   '- Processing effectiveness'}, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    xlim([0, 1]);
    ylim([0, 1]);
    title('Condition Number Improvement');
    axis off;
end

function plot_memory_usage(results)
% Plot memory usage estimation
    if iscell(results.Sigma_tilde)
        F = length(results.Sigma_tilde);
        n = size(results.Sigma_tilde{1}, 1);
        
        % Estimate memory usage
        bytes_per_complex = 16;  % Complex double precision
        memory_per_matrix = n * n * bytes_per_complex / 1024^2;  % MB
        total_memory = F * memory_per_matrix;
        
        % Display memory information
        text(0.1, 0.8, sprintf('Total Memory: %.2f MB', total_memory), 'FontSize', 10);
        text(0.1, 0.6, sprintf('Per Matrix: %.2f MB', memory_per_matrix), 'FontSize', 10);
        text(0.1, 0.4, sprintf('Matrix Count: %d', F), 'FontSize', 10);
        text(0.1, 0.2, sprintf('Matrix Size: %dx%d', n, n), 'FontSize', 10);
        
        axis off;
        title('Memory Usage Estimation');
    end
end

function show_failure_summary(demo_results)
% Show failure summary
    figure('Name', 'Processing Failure Summary', 'Position', [300, 300, 800, 600]);
    
    if ~demo_results.data_generation.success
        text(0.5, 0.8, 'Data Generation Failed:', 'HorizontalAlignment', 'center', ...
             'FontSize', 14, 'FontWeight', 'bold');
        text(0.5, 0.7, demo_results.data_generation.error, ...
             'HorizontalAlignment', 'center', 'FontSize', 12);
    end
    
    if ~demo_results.preprocessing.success
        text(0.5, 0.5, 'Preprocessing Failed:', 'HorizontalAlignment', 'center', ...
             'FontSize', 14, 'FontWeight', 'bold');
        text(0.5, 0.4, demo_results.preprocessing.error, ...
             'HorizontalAlignment', 'center', 'FontSize', 12);
    end
    
    axis off;
    title('Error Summary');
end