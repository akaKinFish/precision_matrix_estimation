function visualize_nmm_data(Svv_cell, Sjj_true_cell, freqs, target_freq_hz)
% VISUALIZE_NMM_DATA - Plot PSD and Connectivity from NMM Simulation
%
% Inputs:
%   Svv_cell      : {F x 1} Sensor Cross-Spectra
%   Sjj_true_cell : {F x 1} Source Cross-Spectra (GT)
%   freqs         : Vector of frequencies
%   target_freq_hz: (Optional) Frequency to plot connectivity matrix (default 10)

    if nargin < 4, target_freq_hz = 10; end

    F = length(freqs);
    Nr = size(Sjj_true_cell{1}, 1);
    Ne = size(Svv_cell{1}, 1);

    % 1. Extract Power Spectra (Diagonals)
    % Convert Cell to Matrix [Nodes x Freqs]
    power_source = zeros(Nr, F);
    power_sensor = zeros(Ne, F);

    for f = 1:F
        power_source(:, f) = real(diag(Sjj_true_cell{f}));
        power_sensor(:, f) = real(diag(Svv_cell{f}));
    end

    % 2. Find Target Frequency Index (e.g., 10Hz)
    [~, f_idx] = min(abs(freqs - target_freq_hz));
    actual_freq = freqs(f_idx);

    % 3. Plotting
    figure('Name', 'NMM Simulation Quality Check', 'Color', 'w', 'Position', [100, 100, 1200, 800]);

    % --- Subplot 1: Source Power Spectrum ---
    subplot(2, 2, 1);
    plot(freqs, 10*log10(power_source + eps), 'LineWidth', 0.5, 'Color', [0.8 0.2 0.2 0.3]); % Semi-transparent red
    hold on;
    plot(freqs, mean(10*log10(power_source + eps), 1), 'k', 'LineWidth', 2); % Mean in black
    title('Source Power Spectrum (Cortex)');
    xlabel('Frequency (Hz)'); ylabel('Power (dB)');
    xlim([min(freqs), max(freqs)]); grid on;

    % --- Subplot 2: Sensor Power Spectrum ---
    subplot(2, 2, 2);
    plot(freqs, 10*log10(power_sensor + eps), 'LineWidth', 0.5, 'Color', [0.2 0.4 0.8 0.3]); % Semi-transparent blue
    hold on;
    plot(freqs, mean(10*log10(power_sensor + eps), 1), 'k', 'LineWidth', 2); % Mean in black
    title('Sensor Power Spectrum (Scalp)');
    xlabel('Frequency (Hz)'); ylabel('Power (dB)');
    xlim([min(freqs), max(freqs)]); grid on;

    % --- Subplot 3: Source Connectivity Matrix at ~10Hz ---
    subplot(2, 2, 3);
    S_src = Sjj_true_cell{f_idx};
    % Compute Coherence: |S_ij| / sqrt(S_ii * S_jj)
    d_src = sqrt(diag(S_src));
    Coh_src = abs(S_src) ./ (d_src * d_src');
    imagesc(Coh_src);
    colormap(gca, 'jet'); colorbar; axis square;
    title(sprintf('Source Coherence (GT) at %.1f Hz', actual_freq));
    xlabel('Regions'); ylabel('Regions');

    % --- Subplot 4: Sensor Connectivity Matrix at ~10Hz ---
    subplot(2, 2, 4);
    S_sens = Svv_cell{f_idx};
    d_sens = sqrt(diag(S_sens));
    Coh_sens = abs(S_sens) ./ (d_sens * d_sens');
    imagesc(Coh_sens);
    colormap(gca, 'jet'); colorbar; axis square;
    title(sprintf('Sensor Coherence at %.1f Hz', actual_freq));
    xlabel('Sensors'); ylabel('Sensors');

end