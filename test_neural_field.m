%% Main_NMM_JSPACE.m - Validate J-SPACE on Neural Mass Models
clear; clc;

% ---------------------------------------------------------
% 1. Load Parameters
% ---------------------------------------------------------
% Replace with your actual path
param_file = 'G:\OneDrive_4_24-11-2025\parameters.mat'; 

if ~exist(param_file, 'file')
    error('Parameters file not found: %s', param_file);
end
fprintf('Loading parameters from: %s\n', param_file);
a = load(param_file);
parameters = a; % Keep consistent naming

% ---------------------------------------------------------
% 2. Run Simulation (Using Wrapper)
% ---------------------------------------------------------
[Svv_cell, Sjj_true_cell, L, freqs] = simulation_nmm_wrapper(parameters);

fprintf('Visualizing NMM Data...\n');
visualize_nmm_data(Svv_cell, Sigma_src_est, freqs, 10); % 检查 10Hz 附近的 alpha 波
drawnow;

% Parameters for J-SPACE
F = length(Svv_cell);
[Ns, Nr] = size(L);

% ---------------------------------------------------------
% 3. Run J-SPACE Solver
% ---------------------------------------------------------
fprintf('\n>>> Running J-SPACE Solver <<<\n');

% Configuration
cfg = struct();
cfg.max_em_iter   = 10;
cfg.grid_size     = 10;
cfg.m_samples     = 1000; % Effective samples (simulated)
cfg.verbose       = true;
cfg.use_gpu       = true; % Turn on if you have GPU
cfg.lambda1_ratio = 1.0;  % Standard smoothing
cfg.lambda3_ratio = 0.1;

% Spatial Laplacian (Optional: construct from Cortex geometry if available)
% For now, we leave it empty or use identity
GraphLaplacian = []; 

tic;
[Omega_est, Sigma_src_est, outs] = solver_jspace(Svv_cell, L, GraphLaplacian, cfg);
time_jspace = toc;

fprintf('J-SPACE completed in %.2f seconds.\n', time_jspace);

% ---------------------------------------------------------
% 4. Performance Evaluation (Metrics)
% ---------------------------------------------------------
fprintf('\n[Step 4] Calculating Metrics...\n');

err_fro_coh = zeros(F, 1);
err_l1_coh  = zeros(F, 1);

for f = 1:F
    % --- Ground Truth Coherence ---
    S_true = Sjj_true_cell{f};
    % Force diagonal real just in case
    P_true_diag = real(diag(S_true));
    % Outer product for denominator
    Denom_true = sqrt(P_true_diag * P_true_diag');
    
    % Calculate Coherence: |S| / sqrt(P_i * P_j)
    Coh_true = abs(S_true) ./ (Denom_true + eps);
    Coh_true(1:Nr+1:end) = 0; % Remove diagonal (self-loop)
    
    % --- Estimated Coherence ---
    S_est = Sigma_src_est{f};
    
    % [FIX] Ensure P_est is strictly real and positive
    P_est_diag = real(diag(S_est));
    P_est_diag(P_est_diag < 0) = 0; % Safety clip
    
    Denom_est = sqrt(P_est_diag * P_est_diag');
    
    % Calculate Coherence
    Coh_est = abs(S_est) ./ (Denom_est + eps);
    Coh_est(1:Nr+1:end) = 0; 
    
    % --- Metrics ---
    diff = Coh_est - Coh_true;
    
    % Normalized Error
    norm_true_fro = norm(Coh_true, 'fro');
    norm_true_l1  = sum(abs(Coh_true(:)));
    
    if norm_true_fro > 0
        err_fro_coh(f) = norm(diff, 'fro') / norm_true_fro;
    else
        err_fro_coh(f) = 0;
    end
    
    if norm_true_l1 > 0
        err_l1_coh(f)  = sum(abs(diff(:))) / norm_true_l1;
    else
        err_l1_coh(f) = 0;
    end
end

fprintf('  Mean Relative Frobenius Error (Coherence): %.4f\n', mean(err_fro_coh));
fprintf('  Mean Relative L1 Error (Coherence):        %.4f\n', mean(err_l1_coh));

% ---------------------------------------------------------
% 5. Visualization (Plotting)
% ---------------------------------------------------------
fprintf('\n[Step 5] Visualizing Results...\n');

figure('Name', 'J-SPACE Performance Analysis', 'Color', 'w', 'Position', [50, 50, 1400, 800]);

% Pick Alpha band peak (~10Hz) or first freq if low
[~, f_idx] = min(abs(freqs - 10)); 
if isempty(f_idx), f_idx = 1; end
target_freq = freqs(f_idx);

% Re-calculate matrices for the specific frequency for plotting
% GT
S_gt = Sjj_true_cell{f_idx};
P_gt = real(diag(S_gt));
C_gt = abs(S_gt) ./ (sqrt(P_gt * P_gt') + eps);
C_gt(1:Nr+1:end) = 0;

% Est
S_est = Sigma_src_est{f_idx};
P_est = real(diag(S_est));
C_est = abs(S_est) ./ (sqrt(P_est * P_est') + eps);
C_est(1:Nr+1:end) = 0;

% [FIX] Ensure data is real for imagesc
C_gt = double(real(C_gt));
C_est = double(real(C_est));

% --- Plotting ---

% 1. GT
subplot(2, 3, 1);
imagesc(C_gt); 
axis square; colorbar; title(sprintf('Ground Truth (%.1f Hz)', target_freq));
xlabel('ROI'); ylabel('ROI'); caxis([0 1]);

% 2. Est
subplot(2, 3, 2);
imagesc(C_est); 
axis square; colorbar; title('J-SPACE Estimated');
xlabel('ROI'); ylabel('ROI'); caxis([0 1]);

% 3. Error
subplot(2, 3, 3);
imagesc(abs(C_est - C_gt));
axis square; colorbar; title('Absolute Error');
xlabel('ROI'); ylabel('ROI'); colormap(gca, 'jet');

% 4. Convergence
subplot(2, 3, 4);
if ~isempty(outs.loglik)
    plot(outs.loglik, '-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
    title('Log-Likelihood'); xlabel('Iteration'); grid on;
    xlim([1, length(outs.loglik)]);
end

% 5. Best Lambda
subplot(2, 3, 5);
if ~isempty(outs.best_lambdas)
    plot(outs.best_lambdas, '-s', 'LineWidth', 2, 'MarkerFaceColor', 'r', 'Color', 'r');
    title('Sparsity (\lambda_2)'); xlabel('Iteration'); grid on;
    xlim([1, length(outs.best_lambdas)]);
end

% 6. Spectrum Error
subplot(2, 3, 6);
plot(freqs, err_fro_coh, 'b-', 'LineWidth', 1.5);
title('Error Spectrum'); xlabel('Frequency (Hz)'); ylabel('Rel. Fro. Error');
grid on; xlim([min(freqs), max(freqs)]);