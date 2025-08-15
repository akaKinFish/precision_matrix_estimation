function figure_handles = create_synthetic_visualizations(preprocessing_results, input_data)
% CREATE_SYNTHETIC_VISUALIZATIONS - Create visualizations for synthetic data demo
%
% File location: demos/visualization/create_synthetic_visualizations.m

    figure_handles = struct();
    
    % Extract data
    Sigma_emp = preprocessing_results.Sigma_emp;
    Sigma_tilde = preprocessing_results.Sigma_tilde;
    D = preprocessing_results.D;
    g_smooth = preprocessing_results.g_smooth;
    
    F = length(Sigma_emp);
    n = size(Sigma_emp{1}, 1);
    
    %% Figure 1: Diagonal Power Evolution
    figure_handles.fig1 = figure('Name', 'Diagonal Power Evolution', 'Position', [100, 600, 1200, 400]);
    
    % Extract diagonal powers across frequencies
    raw_powers = zeros(n, F);
    smoothed_powers = zeros(n, F);
    whitened_powers = zeros(n, F);
    
    for omega = 1:F
        raw_powers(:, omega) = real(diag(Sigma_emp{omega}));
        smoothed_powers(:, omega) = g_smooth{omega};
        whitened_powers(:, omega) = real(diag(Sigma_tilde{omega}));
    end
    
    subplot(1, 3, 1);
    imagesc(raw_powers);
    colorbar;
    title('Raw Diagonal Powers');
    xlabel('Frequency Index');
    ylabel('Channel Index');
    
    subplot(1, 3, 2);
    imagesc(smoothed_powers);
    colorbar;
    title('Smoothed Diagonal Powers');
    xlabel('Frequency Index');
    ylabel('Channel Index');
    
    subplot(1, 3, 3);
    imagesc(whitened_powers);
    colorbar;
    title('Whitened Diagonal Powers');
    xlabel('Frequency Index');
    ylabel('Channel Index');
    caxis([0, 2]); % Focus on range around 1
    
    %% Figure 2: Whitening Quality Assessment
    figure_handles.fig2 = figure('Name', 'Whitening Quality Assessment', 'Position', [150, 500, 1200, 500]);
    
    quality = preprocessing_results.processing_stats.whitening_quality;
    
    subplot(2, 3, 1);
    plot(quality.diagonal_errors');
    title('Diagonal Errors per Frequency');
    xlabel('Channel Index');
    ylabel('Error from Target (1.0)');
    legend('Freq 1', 'Freq 2', '...', 'Location', 'best');
    
    subplot(2, 3, 2);
    plot(max(quality.diagonal_errors, [], 2), 'b-o');
    title('Max Diagonal Error per Frequency');
    xlabel('Frequency Index');
    ylabel('Max Error');
    
    subplot(2, 3, 3);
    semilogy(quality.condition_numbers, 'r-s');
    title('Condition Numbers');
    xlabel('Frequency Index');
    ylabel('Condition Number');
    
    subplot(2, 3, 4);
    plot(quality.min_eigenvalues, 'g-^');
    title('Minimum Eigenvalues');
    xlabel('Frequency Index');
    ylabel('Min Eigenvalue');
    
    subplot(2, 3, 5);
    plot(quality.whitening_effectiveness, 'k-d');
    title('Whitening Effectiveness');
    xlabel('Frequency Index');
    ylabel('Effectiveness Score');
    ylim([0, 1.1]);
    
    subplot(2, 3, 6);
    histogram(quality.whitening_effectiveness, 20);
    title('Effectiveness Distribution');
    xlabel('Effectiveness Score');
    ylabel('Frequency Count');
    
    %% Figure 3: Covariance Matrix Comparison
    figure_handles.fig3 = figure('Name', 'Covariance Matrix Comparison', 'Position', [200, 400, 1500, 400]);
    
    % Show middle frequency
    mid_freq = round(F/2);
    
    subplot(1, 4, 1);
    imagesc(real(Sigma_emp{mid_freq}));
    colorbar;
    title(sprintf('Original Σ (Freq %d)', mid_freq));
    axis square;
    
    subplot(1, 4, 2);
    imagesc(real(Sigma_tilde{mid_freq}));
    colorbar;
    title(sprintf('Whitened Σ̃ (Freq %d)', mid_freq));
    axis square;
    
    subplot(1, 4, 3);
    imagesc(abs(Sigma_emp{mid_freq} - Sigma_tilde{mid_freq}));
    colorbar;
    title('Absolute Difference');
    axis square;
    
    subplot(1, 4, 4);
    % Correlation comparison
    orig_corr = cov2corr(real(Sigma_emp{mid_freq}));
    whitened_corr = cov2corr(real(Sigma_tilde{mid_freq}));
    imagesc(abs(orig_corr - whitened_corr));
    colorbar;
    title('Correlation Change');
    axis square;
    
    %% Figure 4: Processing Pipeline Visualization
    figure_handles.fig4 = figure('Name', 'Processing Pipeline', 'Position', [250, 300, 1400, 600]);
    
    % Show evolution for a few selected channels
    selected_channels = [1, round(n/4), round(n/2), round(3*n/4), n];
    colors = lines(length(selected_channels));
    
    subplot(2, 2, 1);
    hold on;
    for i = 1:length(selected_channels)
        ch = selected_channels(i);
        plot(raw_powers(ch, :), 'Color', colors(i, :), 'LineStyle', '-', 'LineWidth', 2);
    end
    title('Raw Diagonal Powers');
    xlabel('Frequency Index');
    ylabel('Power');
    legend(arrayfun(@(x) sprintf('Ch %d', x), selected_channels, 'UniformOutput', false), ...
           'Location', 'best');
    
    subplot(2, 2, 2);
    hold on;
    for i = 1:length(selected_channels)
        ch = selected_channels(i);
        plot(smoothed_powers(ch, :), 'Color', colors(i, :), 'LineStyle', '--', 'LineWidth', 2);
    end
    title('Smoothed Diagonal Powers');
    xlabel('Frequency Index');
    ylabel('Power');
    
    subplot(2, 2, 3);
    hold on;
    for i = 1:length(selected_channels)
        ch = selected_channels(i);
        whitening_vals = zeros(1, F);
        for omega = 1:F
            whitening_vals(omega) = D{omega}(ch, ch);
        end
        plot(whitening_vals, 'Color', colors(i, :), 'LineStyle', ':', 'LineWidth', 2);
    end
    title('Whitening Matrix Diagonal');
    xlabel('Frequency Index');
    ylabel('1/√(power)');
    
    subplot(2, 2, 4);
    hold on;
    for i = 1:length(selected_channels)
        ch = selected_channels(i);
        plot(whitened_powers(ch, :), 'Color', colors(i, :), 'LineStyle', '-.', 'LineWidth', 2);
    end
    title('Final Whitened Powers');
    xlabel('Frequency Index');
    ylabel('Power');
    yline(1, 'k--', 'Target', 'LineWidth', 1);
    
    fprintf('Created %d figures for synthetic visualization\n', length(fieldnames(figure_handles)));
end

function figure_handles = create_eeg_visualizations(preprocessing_results, input_data)
% CREATE_EEG_VISUALIZATIONS - Create EEG-specific visualizations
%
% File location: demos/visualization/create_eeg_visualizations.m

    figure_handles = struct();
    
    % Extract data
    Sigma_emp = preprocessing_results.Sigma_emp;
    Sigma_tilde = preprocessing_results.Sigma_tilde;
    g_smooth = preprocessing_results.g_smooth;
    
    F = length(Sigma_emp);
    n = size(Sigma_emp{1}, 1);
    frequencies = linspace(1, 40, F); % 1-40 Hz
    
    %% Figure 1: EEG Power Spectrum Analysis
    figure_handles.fig1 = figure('Name', 'EEG Power Spectrum Analysis', 'Position', [100, 600, 1400, 500]);
    
    % Extract average power across channels for each frequency
    avg_power_orig = zeros(1, F);
    avg_power_whitened = zeros(1, F);
    
    for omega = 1:F
        avg_power_orig(omega) = mean(real(diag(Sigma_emp{omega})));
        avg_power_whitened(omega) = mean(real(diag(Sigma_tilde{omega})));
    end
    
    subplot(2, 3, 1);
    semilogy(frequencies, avg_power_orig, 'b-', 'LineWidth', 2);
    hold on;
    semilogy(frequencies, avg_power_whitened, 'r--', 'LineWidth', 2);
    title('Average Power Spectrum');
    xlabel('Frequency (Hz)');
    ylabel('Average Power');
    legend('Original', 'Whitened', 'Location', 'best');
    grid on;
    
    % Highlight frequency bands
    alpha_band = [8, 12];
    beta_band = [13, 30];
    
    subplot(2, 3, 2);
    plot(frequencies, avg_power_orig, 'b-', 'LineWidth', 2);
    hold on;
    xline(alpha_band(1), 'g--', 'Alpha start');
    xline(alpha_band(2), 'g--', 'Alpha end');
    xline(beta_band(1), 'r--', 'Beta start');
    xline(beta_band(2), 'r--', 'Beta end');
    title('Original Power with Frequency Bands');
    xlabel('Frequency (Hz)');
    ylabel('Average Power');
    grid on;
    
    % Channel-wise power distribution
    subplot(2, 3, 3);
    alpha_idx = find(frequencies >= alpha_band(1) & frequencies <= alpha_band(2));
    beta_idx = find(frequencies >= beta_band(1) & frequencies <= beta_band(2));
    
    alpha_power = zeros(n, 1);
    beta_power = zeros(n, 1);
    
    for ch = 1:n
        alpha_power(ch) = mean(arrayfun(@(f) real(Sigma_emp{f}(ch, ch)), alpha_idx));
        beta_power(ch) = mean(arrayfun(@(f) real(Sigma_emp{f}(ch, ch)), beta_idx));
    end
    
    scatter(alpha_power, beta_power, 60, 1:n, 'filled');
    colorbar;
    xlabel('Alpha Power');
    ylabel('Beta Power');
    title('Channel Alpha vs Beta Power');
    
    % Power smoothing effectiveness
    subplot(2, 3, 4);
    for ch = 1:min(5, n) % Show first 5 channels
        raw_ch_power = arrayfun(@(f) real(Sigma_emp{f}(ch, ch)), 1:F);
        smooth_ch_power = arrayfun(@(f) g_smooth{f}(ch), 1:F);
        
        plot(frequencies, raw_ch_power, '-', 'Color', [0.7, 0.7, 0.7]);
        hold on;
        plot(frequencies, smooth_ch_power, '-', 'LineWidth', 2);
    end
    title('Power Smoothing (First 5 Channels)');
    xlabel('Frequency (Hz)');
    ylabel('Power');
    legend('Raw', 'Smoothed', 'Location', 'best');
    
    % Whitening effectiveness across frequency bands
    subplot(2, 3, 5);
    quality = preprocessing_results.processing_stats.whitening_quality;
    plot(frequencies, quality.whitening_effectiveness, 'k-o', 'LineWidth', 2);
    hold on;
    fill([alpha_band(1), alpha_band(2), alpha_band(2), alpha_band(1)], ...
         [0, 0, 1, 1], 'g', 'Alpha', 0.2);
    fill([beta_band(1), beta_band(2), beta_band(2), beta_band(1)], ...
         [0, 0, 1, 1], 'r', 'Alpha', 0.2);
    title('Whitening Effectiveness by Frequency');
    xlabel('Frequency (Hz)');
    ylabel('Effectiveness Score');
    ylim([0, 1.1]);
    
    % Condition number analysis
    subplot(2, 3, 6);
    semilogy(frequencies, quality.condition_numbers, 'r-s', 'LineWidth', 2);
    title('Condition Numbers vs Frequency');
    xlabel('Frequency (Hz)');
    ylabel('Condition Number');
    grid on;
    
    %% Figure 2: EEG Spatial Correlation Analysis
    figure_handles.fig2 = figure('Name', 'EEG Spatial Correlation Analysis', 'Position', [150, 500, 1200, 800]);
    
    % Select representative frequencies
    low_freq_idx = find(frequencies <= 8, 1, 'last');   % <= 8 Hz
    alpha_freq_idx = alpha_idx(round(length(alpha_idx)/2)); % Middle of alpha
    beta_freq_idx = beta_idx(round(length(beta_idx)/2));    % Middle of beta
    high_freq_idx = find(frequencies >= 30, 1, 'first');   % >= 30 Hz
    
    freq_indices = [low_freq_idx, alpha_freq_idx, beta_freq_idx, high_freq_idx];
    freq_labels = {'Low (≤8Hz)', 'Alpha (~10Hz)', 'Beta (~20Hz)', 'High (≥30Hz)'};
    
    for i = 1:4
        freq_idx = freq_indices(i);
        
        % Original correlation
        subplot(2, 4, i);
        orig_corr = cov2corr(real(Sigma_emp{freq_idx}));
        imagesc(orig_corr);
        colorbar;
        title(sprintf('Original Corr - %s', freq_labels{i}));
        axis square;
        caxis([-1, 1]);
        
        % Whitened correlation
        subplot(2, 4, i+4);
        whitened_corr = cov2corr(real(Sigma_tilde{freq_idx}));
        imagesc(whitened_corr);
        colorbar;
        title(sprintf('Whitened Corr - %s', freq_labels{i}));
        axis square;
        caxis([-1, 1]);
    end
    
    %% Figure 3: EEG Topographical Power Maps (Simulated)
    figure_handles.fig3 = figure('Name', 'EEG Topographical Analysis', 'Position', [200, 400, 1400, 600]);
    
    % Create simulated electrode positions (circular layout)
    angles = linspace(0, 2*pi, n+1);
    angles = angles(1:end-1);
    x_pos = cos(angles);
    y_pos = sin(angles);
    
    % Power in different frequency bands
    delta_power = zeros(n, 1);  % 1-4 Hz
    theta_power = zeros(n, 1);  % 4-8 Hz
    alpha_power = zeros(n, 1);  % 8-12 Hz
    beta_power = zeros(n, 1);   % 13-30 Hz
    
    delta_idx = find(frequencies >= 1 & frequencies <= 4);
    theta_idx = find(frequencies >= 4 & frequencies <= 8);
    
    for ch = 1:n
        if ~isempty(delta_idx)
            delta_power(ch) = mean(arrayfun(@(f) real(Sigma_emp{f}(ch, ch)), delta_idx));
        end
        if ~isempty(theta_idx)
            theta_power(ch) = mean(arrayfun(@(f) real(Sigma_emp{f}(ch, ch)), theta_idx));
        end
        alpha_power(ch) = mean(arrayfun(@(f) real(Sigma_emp{f}(ch, ch)), alpha_idx));
        beta_power(ch) = mean(arrayfun(@(f) real(Sigma_emp{f}(ch, ch)), beta_idx));
    end
    
    % Plot topographical maps
    subplot(2, 4, 1);
    scatter(x_pos, y_pos, 100, delta_power, 'filled');
    colorbar;
    title('Delta Power (1-4 Hz)');
    axis equal; axis off;
    
    subplot(2, 4, 2);
    scatter(x_pos, y_pos, 100, theta_power, 'filled');
    colorbar;
    title('Theta Power (4-8 Hz)');
    axis equal; axis off;
    
    subplot(2, 4, 3);
    scatter(x_pos, y_pos, 100, alpha_power, 'filled');
    colorbar;
    title('Alpha Power (8-12 Hz)');
    axis equal; axis off;
    
    subplot(2, 4, 4);
    scatter(x_pos, y_pos, 100, beta_power, 'filled');
    colorbar;
    title('Beta Power (13-30 Hz)');
    axis equal; axis off;
    
    % After whitening
    delta_power_w = zeros(n, 1);
    theta_power_w = zeros(n, 1);
    alpha_power_w = zeros(n, 1);
    beta_power_w = zeros(n, 1);
    
    for ch = 1:n
        if ~isempty(delta_idx)
            delta_power_w(ch) = mean(arrayfun(@(f) real(Sigma_tilde{f}(ch, ch)), delta_idx));
        end
        if ~isempty(theta_idx)
            theta_power_w(ch) = mean(arrayfun(@(f) real(Sigma_tilde{f}(ch, ch)), theta_idx));
        end
        alpha_power_w(ch) = mean(arrayfun(@(f) real(Sigma_tilde{f}(ch, ch)), alpha_idx));
        beta_power_w(ch) = mean(arrayfun(@(f) real(Sigma_tilde{f}(ch, ch)), beta_idx));
    end
    
    subplot(2, 4, 5);
    scatter(x_pos, y_pos, 100, delta_power_w, 'filled');
    colorbar;
    title('Delta Power - Whitened');
    axis equal; axis off;
    
    subplot(2, 4, 6);
    scatter(x_pos, y_pos, 100, theta_power_w, 'filled');
    colorbar;
    title('Theta Power - Whitened');
    axis equal; axis off;
    
    subplot(2, 4, 7);
    scatter(x_pos, y_pos, 100, alpha_power_w, 'filled');
    colorbar;
    title('Alpha Power - Whitened');
    axis equal; axis off;
    
    subplot(2, 4, 8);
    scatter(x_pos, y_pos, 100, beta_power_w, 'filled');
    colorbar;
    title('Beta Power - Whitened');
    axis equal; axis off;
    
    fprintf('Created %d figures for EEG visualization\n', length(fieldnames(figure_handles)));
end

function figure_handles = create_parameter_sensitivity_plots(results)
% CREATE_PARAMETER_SENSITIVITY_PLOTS - Visualize parameter sensitivity analysis
%
% File location: demos/visualization/create_parameter_sensitivity_plots.m

    figure_handles = struct();
    
    %% Figure 1: Smoothing Method Comparison
    figure_handles.fig1 = figure('Name', 'Smoothing Method Comparison', 'Position', [100, 600, 1200, 400]);
    
    methods = fieldnames(results.parameter_tests.smoothing);
    colors = lines(length(methods));
    
    subplot(1, 3, 1);
    hold on;
    for i = 1:length(methods)
        method = methods{i};
        quality = results.parameter_tests.smoothing.(method).processing_stats.whitening_quality;
        plot(quality.whitening_effectiveness, 'Color', colors(i, :), 'LineWidth', 2, ...
             'DisplayName', strrep(method, '_', ' '));
    end
    title('Whitening Effectiveness by Method');
    xlabel('Frequency Index');
    ylabel('Effectiveness Score');
    legend('Location', 'best');
    grid on;
    
    subplot(1, 3, 2);
    method_scores = zeros(length(methods), 1);
    for i = 1:length(methods)
        method = methods{i};
        quality = results.parameter_tests.smoothing.(method).processing_stats.whitening_quality;
        method_scores(i) = mean(quality.whitening_effectiveness);
    end
    bar(method_scores);
    set(gca, 'XTickLabel', strrep(methods, '_', ' '));
    title('Average Effectiveness by Method');
    ylabel('Mean Effectiveness Score');
    
    subplot(1, 3, 3);
    timing_scores = zeros(length(methods), 1);
    for i = 1:length(methods)
        method = methods{i};
        timing_scores(i) = results.parameter_tests.smoothing.(method).timing.total;
    end
    bar(timing_scores);
    set(gca, 'XTickLabel', strrep(methods, '_', ' '));
    title('Processing Time by Method');
    ylabel('Time (seconds)');
    
    %% Figure 2: Window Size Analysis
    figure_handles.fig2 = figure('Name', 'Window Size Analysis', 'Position', [150, 500, 1200, 400]);
    
    window_fields = fieldnames(results.parameter_tests.window_size);
    window_sizes = zeros(length(window_fields), 1);
    window_effectiveness = zeros(length(window_fields), 1);
    window_smoothness = zeros(length(window_fields), 1);
    
    for i = 1:length(window_fields)
        field = window_fields{i};
        % Extract window size from field name (e.g., 'size_5' -> 5)
        window_sizes(i) = str2double(regexp(field, '\d+', 'match', 'once'));
        
        quality = results.parameter_tests.window_size.(field).processing_stats.whitening_quality;
        window_effectiveness(i) = mean(quality.whitening_effectiveness);
        
        % Measure smoothness as inverse of variance in effectiveness
        window_smoothness(i) = 1 / (var(quality.whitening_effectiveness) + 1e-6);
    end
    
    [~, sort_idx] = sort(window_sizes);
    window_sizes = window_sizes(sort_idx);
    window_effectiveness = window_effectiveness(sort_idx);
    window_smoothness = window_smoothness(sort_idx);
    
    subplot(1, 3, 1);
    plot(window_sizes, window_effectiveness, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
    title('Effectiveness vs Window Size');
    xlabel('Window Size');
    ylabel('Mean Effectiveness');
    grid on;
    
    subplot(1, 3, 2);
    plot(window_sizes, window_smoothness, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
    title('Smoothness vs Window Size');
    xlabel('Window Size');
    ylabel('Smoothness (1/variance)');
    grid on;
    
    subplot(1, 3, 3);
    % Combined score
    normalized_eff = (window_effectiveness - min(window_effectiveness)) / ...
                    (max(window_effectiveness) - min(window_effectiveness) + 1e-6);
    normalized_smooth = (window_smoothness - min(window_smoothness)) / ...
                       (max(window_smoothness) - min(window_smoothness) + 1e-6);
    combined_score = 0.7 * normalized_eff + 0.3 * normalized_smooth;
    
    plot(window_sizes, combined_score, 'go-', 'LineWidth', 2, 'MarkerSize', 8);
    title('Combined Score vs Window Size');
    xlabel('Window Size');
    ylabel('Combined Score');
    grid on;
    
    %% Figure 3: Loading Factor Analysis
    figure_handles.fig3 = figure('Name', 'Loading Factor Analysis', 'Position', [200, 400, 1400, 500]);
    
    loading_fields = fieldnames(results.parameter_tests.loading_factor);
    loading_factors = zeros(length(loading_fields), 1);
    loading_effectiveness = zeros(length(loading_fields), 1);
    loading_condition = zeros(length(loading_fields), 1);
    loading_min_eigenvals = zeros(length(loading_fields), 1);
    
    for i = 1:length(loading_fields)
        field = loading_fields{i};
        % Extract loading factor from field name
        factor_str = regexp(field, '[\d\.]+', 'match', 'once');
        loading_factors(i) = str2double(factor_str);
        
        quality = results.parameter_tests.loading_factor.(field).processing_stats.whitening_quality;
        loading_effectiveness(i) = mean(quality.whitening_effectiveness);
        loading_condition(i) = mean(quality.condition_numbers);
        loading_min_eigenvals(i) = min(quality.min_eigenvalues);
    end
    
    [~, sort_idx] = sort(loading_factors);
    loading_factors = loading_factors(sort_idx);
    loading_effectiveness = loading_effectiveness(sort_idx);
    loading_condition = loading_condition(sort_idx);
    loading_min_eigenvals = loading_min_eigenvals(sort_idx);
    
    subplot(2, 2, 1);
    semilogx(loading_factors, loading_effectiveness, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
    title('Effectiveness vs Loading Factor');
    xlabel('Loading Factor');
    ylabel('Mean Effectiveness');
    grid on;
    
    subplot(2, 2, 2);
    loglog(loading_factors, loading_condition, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
    title('Condition Number vs Loading Factor');
    xlabel('Loading Factor');
    ylabel('Mean Condition Number');
    grid on;
    
    subplot(2, 2, 3);
    semilogx(loading_factors, loading_min_eigenvals, 'go-', 'LineWidth', 2, 'MarkerSize', 8);
    title('Min Eigenvalue vs Loading Factor');
    xlabel('Loading Factor');
    ylabel('Min Eigenvalue');
    grid on;
    
    subplot(2, 2, 4);
    % Trade-off analysis
    % Normalize metrics
    norm_eff = loading_effectiveness / max(loading_effectiveness);
    norm_cond = min(loading_condition) ./ loading_condition; % Smaller is better
    norm_eigen = (loading_min_eigenvals - min(loading_min_eigenvals)) / ...
                 (max(loading_min_eigenvals) - min(loading_min_eigenvals) + 1e-6);
    
    combined = 0.5 * norm_eff + 0.3 * norm_cond + 0.2 * norm_eigen;
    
    semilogx(loading_factors, combined, 'ko-', 'LineWidth', 2, 'MarkerSize', 8);
    title('Combined Performance Score');
    xlabel('Loading Factor');
    ylabel('Combined Score');
    grid on;
    
    fprintf('Created %d figures for parameter sensitivity analysis\n', length(fieldnames(figure_handles)));
end

function figure_handles = create_noise_robustness_plots(results)
% CREATE_NOISE_ROBUSTNESS_PLOTS - Visualize noise robustness analysis
%
% File location: demos/visualization/create_noise_robustness_plots.m

    figure_handles = struct();
    
    noise_levels = results.noise_levels;
    n_noise = length(noise_levels);
    
    % Extract metrics for each noise level
    effectiveness_scores = zeros(n_noise, 1);
    condition_numbers = zeros(n_noise, 1);
    min_eigenvalues = zeros(n_noise, 1);
    diagonal_errors = zeros(n_noise, 1);
    processing_times = zeros(n_noise, 1);
    
    noise_fields = fieldnames(results.noise_tests);
    
    for i = 1:n_noise
        field = noise_fields{i};
        test_result = results.noise_tests.(field);
        
        quality = test_result.processing_stats.whitening_quality;
        effectiveness_scores(i) = mean(quality.whitening_effectiveness);
        condition_numbers(i) = mean(quality.condition_numbers);
        min_eigenvalues(i) = min(quality.min_eigenvalues);
        diagonal_errors(i) = mean(max(quality.diagonal_errors, [], 2));
        processing_times(i) = test_result.timing.total;
    end
    
    %% Figure 1: Robustness Overview
    figure_handles.fig1 = figure('Name', 'Noise Robustness Overview', 'Position', [100, 600, 1400, 500]);
    
    subplot(2, 3, 1);
    semilogx(noise_levels, effectiveness_scores, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
    title('Whitening Effectiveness vs Noise');
    xlabel('Noise Level');
    ylabel('Mean Effectiveness Score');
    grid on;
    ylim([0, 1.1]);
    
    subplot(2, 3, 2);
    loglog(noise_levels, condition_numbers, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
    title('Condition Numbers vs Noise');
    xlabel('Noise Level');
    ylabel('Mean Condition Number');
    grid on;
    
    subplot(2, 3, 3);
    semilogx(noise_levels, min_eigenvalues, 'go-', 'LineWidth', 2, 'MarkerSize', 8);
    title('Min Eigenvalues vs Noise');
    xlabel('Noise Level');
    ylabel('Min Eigenvalue');
    grid on;
    
    subplot(2, 3, 4);
    semilogx(noise_levels, diagonal_errors, 'mo-', 'LineWidth', 2, 'MarkerSize', 8);
    title('Diagonal Errors vs Noise');
    xlabel('Noise Level');
    ylabel('Mean Max Diagonal Error');
    grid on;
    
    subplot(2, 3, 5);
    semilogx(noise_levels, processing_times, 'co-', 'LineWidth', 2, 'MarkerSize', 8);
    title('Processing Time vs Noise');
    xlabel('Noise Level');
    ylabel('Processing Time (s)');
    grid on;
    
    subplot(2, 3, 6);
    % Robustness score (higher is better)
    robustness_score = effectiveness_scores .* (1 ./ (1 + diagonal_errors)) .* ...
                      (1 ./ (1 + log10(condition_numbers)));
    
    semilogx(noise_levels, robustness_score, 'ko-', 'LineWidth', 2, 'MarkerSize', 8);
    title('Overall Robustness Score');
    xlabel('Noise Level');
    ylabel('Robustness Score');
    grid on;
    
    %% Figure 2: Detailed Noise Analysis
    figure_handles.fig2 = figure('Name', 'Detailed Noise Analysis', 'Position', [150, 500, 1200, 800]);
    
    % Show effectiveness evolution for different noise levels
    selected_noise_idx = [1, round(n_noise/3), round(2*n_noise/3), n_noise];
    colors = lines(length(selected_noise_idx));
    
    subplot(2, 2, 1);
    hold on;
    for i = 1:length(selected_noise_idx)
        idx = selected_noise_idx(i);
        field = noise_fields{idx};
        quality = results.noise_tests.(field).processing_stats.whitening_quality;
        
        plot(quality.whitening_effectiveness, 'Color', colors(i, :), 'LineWidth', 2, ...
             'DisplayName', sprintf('Noise %.3f', noise_levels(idx)));
    end
    title('Effectiveness vs Frequency');
    xlabel('Frequency Index');
    ylabel('Effectiveness Score');
    legend('Location', 'best');
    grid on;
    
    subplot(2, 2, 2);
    hold on;
    for i = 1:length(selected_noise_idx)
        idx = selected_noise_idx(i);
        field = noise_fields{idx};
        quality = results.noise_tests.(field).processing_stats.whitening_quality;
        
        semilogy(quality.condition_numbers, 'Color', colors(i, :), 'LineWidth', 2, ...
                'DisplayName', sprintf('Noise %.3f', noise_levels(idx)));
    end
    title('Condition Numbers vs Frequency');
    xlabel('Frequency Index');
    ylabel('Condition Number');
    legend('Location', 'best');
    grid on;
    
    % Degradation analysis
    subplot(2, 2, 3);
    degradation = (effectiveness_scores(1) - effectiveness_scores) / effectiveness_scores(1) * 100;
    semilogx(noise_levels, degradation, 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
    title('Performance Degradation');
    xlabel('Noise Level');
    ylabel('Degradation (%)');
    grid on;
    
    % Critical noise level analysis
    subplot(2, 2, 4);
    % Find noise level where effectiveness drops below 80% of clean performance
    threshold = 0.8 * effectiveness_scores(1);
    critical_idx = find(effectiveness_scores < threshold, 1, 'first');
    
    semilogx(noise_levels, effectiveness_scores, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
    hold on;
    yline(threshold, 'r--', 'Critical Threshold', 'LineWidth', 2);
    if ~isempty(critical_idx)
        xline(noise_levels(critical_idx), 'g--', 'Critical Noise Level', 'LineWidth', 2);
        plot(noise_levels(critical_idx), effectiveness_scores(critical_idx), 'ro', ...
             'MarkerSize', 12, 'LineWidth', 3);
    end
    title('Critical Noise Level Analysis');
    xlabel('Noise Level');
    ylabel('Effectiveness Score');
    grid on;
    
    fprintf('Created %d figures for noise robustness analysis\n', length(fieldnames(figure_handles)));
end

function corr_matrix = cov2corr(cov_matrix)
% Convert covariance matrix to correlation matrix
    d = sqrt(diag(cov_matrix));
    d(d == 0) = 1; % Avoid division by zero
    corr_matrix = cov_matrix ./ (d * d');
end