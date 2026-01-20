%% J-SPACE Spline Comparison Script
% Import required helper functions
clear all;
close all;

% --- Add J-SPACE Spline Paths ---
% Assuming the folder structure is maintained relative to this script
addpath(genpath('core'));
addpath(genpath('utils'));
addpath(genpath('helpers'));

import functions.auxx.ModelVectorization.*;
import guide.Visualization.*;
import functions.auxx.ZeroInflatedModels.*;
import functions.auxx.Refine_Solution.*;
import functions.auxx.OptimizedOperations.*;

Cortex = load("templates/Cortex.mat");
rng(2026);

% Path to the JSON file with model result metadata. Modify this directions manually
json_path = 'G:\OneDrive - CCLAB\New_Data_XiAlphaNET\xialphanet_newresults22\XIALPHANET.json';
dir_data = 'G:\BC-V_result_12_6_part1\CHBMP\2pre\crossSpec';

[dataset_dir, ~, ~] = fileparts(json_path);
dataset = jsondecode(fileread(json_path));
dataset.Location = dataset_dir;

parameters = load(fullfile(dataset_dir, 'structural', 'parameters.mat'));
parameters.Model = parameters.Compact_Model;
parameters.Compact_Model = parameters.Compact_Model;
parameters.Dimensions = parameters.Dimensions;
parameters.Dimensions.Nv = parameters.Dimensions.Nr;

% Clear temporary variables to save memory
clear Compact_Model Model Dimensions;

% Import additional functions
import functions.*;
import functions.auxx.*;
import functions.auxx.BayesOptimization.*;
import functions.auxx.CrossSpectrum.*;
import functions.auxx.ModelVectorization.*;
import functions.auxx.Simulations.*;
import functions.auxx.TOperator.*;
import tools.*;
import functions.auxx.DataPreprocessing.*;
import functions.auxx.Simulations.inverse.*;
import functions.auxx.Simulations.private.*;

% Set up simulation parameters
Ne = parameters.Dimensions.Ne;  % Number of electrodes
Nr = parameters.Dimensions.Nr;  % Number of ROIs
Nv = parameters.Dimensions.Nr;  % Voxel dimension = ROI dimension
Nw = parameters.Dimensions.Nw;  % Number of frequency bins
Nsub = 5;                       % Number of simulations
N_wishart = 1;

conn_spec = norm(parameters.Model.C, 'fro');
disp("--> Estimating source cross-spectrum");

% Xi-AlphaNET properties
properties.model_params.nFreqs = Nw;
properties.model_params.BayesIter_Delay = 50;
properties.model_params.BayesIter_Reg1 = 20;
properties.model_params.BayesIter_Reg2 = 100;
properties.model_params.Nrand1 = 10;
properties.model_params.Nrand2 = 50;
properties.model_params.delay.lambda_space_cd = [[0.4, 1.6]; [10^(-10), 1/conn_spec]];
properties.general_params.parallel.conn_delay = 1;
properties.model_params.stoch1 = 1;
properties.model_params.stoch2 = 1;
properties.model_params.tensor_field.default = 0;

% Load model and transformation matrices
L = parameters.Model.K;  % Transformation matrix for cross-spectrum
K = parameters.Model.K;
R = parameters.Compact_Model.R;

% Select subjects for simulation
subject_folders = dir(fullfile(dir_data, '*'));
subject_folders = subject_folders([subject_folders.isdir] & ~startsWith({subject_folders.name}, '.'));
if length(subject_folders) < Nsub
    Nsub = length(subject_folders);
end
selected_folders = subject_folders(randperm(length(subject_folders), Nsub));

% Define frequency bands
bands = {
    'Delta', 1, 4;
    'Theta', 4, 8;
    'Alpha', 8, 13;
    'Beta', 13, 30;
    'Gamma', 30, 50
    };
Nbands = size(bands, 1);

% Pre-allocate variables
Nw = 47; % Assuming freq is defined elsewhere
Nroi = 360;        % Adjust if needed
xl = zeros(Nroi, Nw, Nsub);
l  = zeros(Nroi, Nw, Nsub);
mf = zeros(Nroi, Nw, Nsub);
js = zeros(Nroi, Nw, Nsub); % J-SPACE Power

% Frobenius and L1 norm benchmark errors for each subject
XA_Spectrum_Fro = zeros(Nsub, 1); XA_Spectrum_FroR = zeros(Nsub, 1);
XA_Spectrum_L1 = zeros(Nsub, 1);  XA_Spectrum_L1R = zeros(Nsub, 1);
mn_Spectrum_Fro = zeros(Nsub, 1); mn_Spectrum_FroR = zeros(Nsub, 1);
mn_Spectrum_L1 = zeros(Nsub, 1);  mn_Spectrum_L1R = zeros(Nsub, 1);

% J-SPACE Spectrum Metrics
JS_Spectrum_Fro = zeros(Nsub, 1); JS_Spectrum_FroR = zeros(Nsub, 1);
JS_Spectrum_L1 = zeros(Nsub, 1);  JS_Spectrum_L1R = zeros(Nsub, 1);

XA_Coherence_Fro = zeros(Nsub, 1); XA_Coherence_FroR = zeros(Nsub, 1);
XA_Coherence_L1 = zeros(Nsub, 1);  XA_Coherence_L1R = zeros(Nsub, 1);
mn_Coherence_Fro = zeros(Nsub, 1); mn_Coherence_FroR = zeros(Nsub, 1);
mn_Coherence_L1 = zeros(Nsub, 1);  mn_Coherence_L1R = zeros(Nsub, 1);

% J-SPACE Coherence Metrics
JS_Coherence_Fro = zeros(Nsub, 1); JS_Coherence_FroR = zeros(Nsub, 1);
JS_Coherence_L1 = zeros(Nsub, 1);  JS_Coherence_L1R = zeros(Nsub, 1);

XA_Phase_Fro = zeros(Nsub, 1); XA_Phase_FroR = zeros(Nsub, 1);
XA_Phase_L1 = zeros(Nsub, 1);  XA_Phase_L1R = zeros(Nsub, 1);
mn_Phase_Fro = zeros(Nsub, 1); mn_Phase_FroR = zeros(Nsub, 1);
mn_Phase_L1 = zeros(Nsub, 1);  mn_Phase_L1R = zeros(Nsub, 1);

% J-SPACE Phase Metrics
JS_Phase_Fro = zeros(Nsub, 1); JS_Phase_FroR = zeros(Nsub, 1);
JS_Phase_L1 = zeros(Nsub, 1);  JS_Phase_L1R = zeros(Nsub, 1);


% Simulation loop
for j = 1:Nsub
    tic;
    fprintf('\n=== Processing Subject %d / %d ===\n', j, Nsub);
    
    % Load subject data
    subject_folder = selected_folders(j).name;
    mat_file_path = fullfile(dir_data, subject_folder, [subject_folder, '.mat']);
    data_struct = load(mat_file_path);
    freq = data_struct.data_struct.freqrange(1:Nw);
    parameters.Data.freq = freq;
    
    % Source cross-spectrum simulation
    disp('->> Estimating Source Cross Using a Neural Mass Simulation');
    [Sjj, ~] = functions.auxx.Simulations.neural_mass_simulation(0, json_path);
    
    % Generate scalp-level cross-spectrum with noise
    N_wishart_sim = 2 * size(L, 1); 
    
    parfor i = 1:Nw
        Sjj_cross(:, :, i) = generate_complex_wishart(Sjj(:, :, i), N_wishart_sim);
        
        % Projection
        Svv_clean = L * Sjj_cross(:, :, i) * L';
        
        % Add sensor noise to ensure full rank
        noise_level = 0.05 * trace(Svv_clean) / size(L, 1);
        Svv_cross(:, :, i) = Svv_clean + noise_level * eye(size(L, 1)); 
    end
    
    [Log_Spec,freq]= log_spectrum(Svv_cross,freq);
    [n, ~, F] = size(Svv_cross);
    diag_mask = repmat(logical(eye(n)), [1, 1, F]);
    Svv_cross(diag_mask) = real(Svv_cross(diag_mask));
    toc;
    
    % ---------------------------------------------------------------------
    % Xi-AlphaNET estimation
    % ---------------------------------------------------------------------
    data.Cross = Svv_cross;
    data.age = 25;  % Subject age
    data.freq = freq;
    disp('->> Xi-AlphaNeT Inverse Solution');
    [x, ~, G, x0] = Xi_ALphaNET(properties, data, parameters);
    source_act_cross = functions.auxx.CrossSpectrum.eval_source_conn(x.Solution, data.freq, parameters.Model.R, properties, parameters);
    XA_Sjj_cross = source_act_cross.Cross.Full;
    
    % ---------------------------------------------------------------------
    % Mean cross-spectrum
    % ---------------------------------------------------------------------
    mn_Sjj_cross = mn_cross(Svv_cross, K, 0);
    
    % ---------------------------------------------------------------------
    % eLORETA
    % ---------------------------------------------------------------------
    disp('->> eLORETA Processing...');
    eL_Sjj_cross = zeros(Nroi, Nroi, Nw);
    parfor i = 1:Nw
        source = inverse(Svv_cross(:,:,i), L);
        eL_Sjj_cross(:,:,i) = source.eloreata.Sjj;
    end
    
    % ---------------------------------------------------------------------
    % [CHANGED] J-SPACE SPLINE Processing (Route A)
    % ---------------------------------------------------------------------
    disp('->> J-SPACE SPLINE Processing (Group-Lasso + P-Spline)...');
    
    % --- Configure J-SPACE Spline ---
    cfg_spline = [];
    
    % 1. Input Processing
    cfg_spline.svv.jitter = 1e-10; % Robustness
    cfg_spline.noise.sigma2 = 1e-3; % Sensor noise assumption
    
    % 2. Spline Settings (Model Complexity)
    cfg_spline.spline.K = 8;          % Number of B-spline bases (e.g. 6-10)
    cfg_spline.spline.degree = 3;     % Cubic spline
    cfg_spline.spline.diff_order = 2; % P-spline penalty order
    
    % 3. M-step Hyperparameters (Tuning Knobs)
    % lambda1: Group Lasso (Structure Sparsity) - Controls edge density
    % lambda_ps: P-spline (Frequency Smoothness) - Controls smoothness
    cfg_spline.mstep.lambda1   = 0.05;  
    cfg_spline.mstep.lambda_ps = 0.1;   
    cfg_spline.mstep.lambda2   = 0.0;   % Spatial L2 (optional, off)
    
    % 4. Optimization Settings (FISTA)
    cfg_spline.mstep.fista.max_iter = 50;
    cfg_spline.mstep.fista.verbose = false;
    
    % 5. EM Settings
    cfg_spline.em.max_iter = 15;
    cfg_spline.em.update_rate = 0.3; % Inertia update
    cfg_spline.em.verbose = true;
    
    % 6. Structural Prior (DWI)
    % Use the loaded compact model connectivity
    dwi_C_prior = parameters.Compact_Model.C;
    cfg_spline.dwi.alpha = 1.0;
    
    % --- Run Solver ---
    % Inputs: Svv (Ne x Ne x F), L (Ne x N), freq, dwi_C, cfg
    [~, JS_Sjj_cell, outs_js] = run_jspace_pspline_real(Svv_cross, L, freq, dwi_C_prior, cfg_spline);
    
    % Map output to variable name expected by metrics (JS_Sjj_cross can be a cell)
    JS_Sjj_cross = JS_Sjj_cell; 
    
    % ---------------------------------------------------------------------
    % Metric Collection
    % ---------------------------------------------------------------------
    
    % Initialize mean power array for each ROI and frequency
    mn_power = zeros(Nroi, Nw);
    
    % Extract and log-transform diagonal power for each ROI
    for i = 1:Nw
        xl(:, i, j)= log10(real(diag(XA_Sjj_cross(:, :, i))));
        
        % Original cross-spectrum power (log-transformed)
        l(:, i, j) = log10(real(diag(Sjj_cross(:, :, i))));
        
        % eLORETA power (not log-transformed)
        lo(:, i) = real(diag(eL_Sjj_cross(:, :, i)));
        
        % J-SPACE power [UPDATED for Cell Array]
        js(:, i, j) = log10(real(diag(JS_Sjj_cross{i})) + eps);
        
        % Mean cross-spectrum power
        mn_power(:, i) = real(diag(mn_Sjj_cross(:, :, i)));
    end
    
    % Fit FOOOF and xi-alpha models per ROI
    parfor roi = 1:Nroi
        % Fit the FOOOF model to the mean power spectrum
        mn_fooof = fooof_matlab(freq, mn_power(roi, :), [min(freq), max(freq)], 1);
        mf(roi, :, j) = mn_fooof.model_fit(:).';
        % Fit the xi-alpha model to eLORETA diagonal power
        [params_xa_lo, fit_xa_lo] = functions.auxx.GenerateSourceSample.fit_xi_alpha_multi(lo(roi,:), freq,10, 0);
        el_fit_xa(roi,:,j) = fit_xa_lo;
    end
    
    % Store ground-truth cross-spectrum (used for comparisons)
    ground_truth_Sjj = Sjj_cross;
    
    % Consolidate J-SPACE result cell array into 3-D tensor for metrics
    if iscell(JS_Sjj_cross)
        JS_tensor = cat(3, JS_Sjj_cross{:});
    else
        JS_tensor = JS_Sjj_cross;
    end
    
    % ----- Spectral Comparison -----
    % Frobenius Norm Comparison
    XA_Spectrum_Fro(j)   = tensor_norm(XA_Sjj_cross - ground_truth_Sjj, 2);
    mn_Spectrum_Fro(j)   = tensor_norm(mn_Sjj_cross - ground_truth_Sjj, 2);
    eL_Spectrum_Fro(j)   = tensor_norm(eL_Sjj_cross - ground_truth_Sjj, 2);
    JS_Spectrum_Fro(j)   = tensor_norm(JS_tensor - ground_truth_Sjj, 2); 
    
    % Relative Frobenius Norm (normalized by ground-truth)
    XA_Spectrum_FroR(j)  = XA_Spectrum_Fro(j) / tensor_norm(ground_truth_Sjj, 2);
    mn_Spectrum_FroR(j)  = mn_Spectrum_Fro(j) / tensor_norm(ground_truth_Sjj, 2);
    eL_Spectrum_FroR(j)  = eL_Spectrum_Fro(j) / tensor_norm(ground_truth_Sjj, 2);
    JS_Spectrum_FroR(j)  = JS_Spectrum_Fro(j) / tensor_norm(ground_truth_Sjj, 2); 
    
    % L1 Norm Comparison
    XA_Spectrum_L1(j)    = tensor_norm(XA_Sjj_cross - ground_truth_Sjj, 1);
    mn_Spectrum_L1(j)    = tensor_norm(mn_Sjj_cross - ground_truth_Sjj, 1);
    eL_Spectrum_L1(j)    = tensor_norm(eL_Sjj_cross - ground_truth_Sjj, 1);
    JS_Spectrum_L1(j)    = tensor_norm(JS_tensor - ground_truth_Sjj, 1); 
    
    % Relative L1 Norm (normalized by ground-truth)
    XA_Spectrum_L1R(j)   = XA_Spectrum_L1(j) / tensor_norm(ground_truth_Sjj, 1);
    mn_Spectrum_L1R(j)   = mn_Spectrum_L1(j) / tensor_norm(ground_truth_Sjj, 1);
    eL_Spectrum_L1R(j)   = eL_Spectrum_L1(j) / tensor_norm(ground_truth_Sjj, 1);
    JS_Spectrum_L1R(j)   = JS_Spectrum_L1(j) / tensor_norm(ground_truth_Sjj, 1); 
    
    % ----- Coherence Comparison -----
    % Compare coherence matrices for each model
    coherence_Sjj = coherence(ground_truth_Sjj);
    coherence_XA  = coherence(XA_Sjj_cross);
    coherence_mn  = coherence(mn_Sjj_cross);
    coherence_eL  = coherence(eL_Sjj_cross);
    coherence_JS  = coherence(JS_tensor); 
    
    % Frobenius Norm Comparison for Coherence
    XA_Coherence_Fro(j)  = tensor_norm(coherence_XA - coherence_Sjj, 2);
    mn_Coherence_Fro(j)  = tensor_norm(coherence_mn - coherence_Sjj, 2);
    eL_Coherence_Fro(j)  = tensor_norm(coherence_eL - coherence_Sjj, 2);
    JS_Coherence_Fro(j)  = tensor_norm(coherence_JS - coherence_Sjj, 2); 
    
    % Relative Frobenius Norm
    XA_Coherence_FroR(j) = XA_Coherence_Fro(j) / tensor_norm(coherence_Sjj, 2);
    mn_Coherence_FroR(j) = mn_Coherence_Fro(j) / tensor_norm(coherence_Sjj, 2);
    eL_Coherence_FroR(j) = eL_Coherence_Fro(j) / tensor_norm(coherence_Sjj, 2);
    JS_Coherence_FroR(j) = JS_Coherence_Fro(j) / tensor_norm(coherence_Sjj, 2); 
    
    % L1 Norm Comparison for Coherence
    XA_Coherence_L1(j)   = tensor_norm(coherence_XA - coherence_Sjj, 1);
    mn_Coherence_L1(j)   = tensor_norm(coherence_mn - coherence_Sjj, 1);
    eL_Coherence_L1(j)   = tensor_norm(coherence_eL - coherence_Sjj, 1);
    JS_Coherence_L1(j)   = tensor_norm(coherence_JS - coherence_Sjj, 1); 
    
    % Relative L1 Norm
    XA_Coherence_L1R(j)  = XA_Coherence_L1(j) / tensor_norm(coherence_Sjj, 1);
    mn_Coherence_L1R(j)  = mn_Coherence_L1(j) / tensor_norm(coherence_Sjj, 1);
    eL_Coherence_L1R(j)  = eL_Coherence_L1(j) / tensor_norm(coherence_Sjj, 1);
    JS_Coherence_L1R(j)  = JS_Coherence_L1(j) / tensor_norm(coherence_Sjj, 1); 
    
    % ----- Phase Comparison -----
    % Compare phase angles across models
    angle_Sjj = angle(ground_truth_Sjj);
    angle_XA  = angle(XA_Sjj_cross);
    angle_mn  = angle(mn_Sjj_cross);
    angle_eL  = angle(eL_Sjj_cross);
    angle_JS  = angle(JS_tensor); 
    
    % Frobenius Norm Comparison for Phase
    XA_Phase_Fro(j)      = tensor_norm(angle_XA - angle_Sjj, 2);
    mn_Phase_Fro(j)      = tensor_norm(angle_mn - angle_Sjj, 2);
    eL_Phase_Fro(j)      = tensor_norm(angle_eL - angle_Sjj, 2);
    JS_Phase_Fro(j)      = tensor_norm(angle_JS - angle_Sjj, 2); 
    
    % Relative Frobenius Norm
    XA_Phase_FroR(j)     = XA_Phase_Fro(j) / tensor_norm(angle_Sjj, 2);
    mn_Phase_FroR(j)     = mn_Phase_Fro(j) / tensor_norm(angle_Sjj, 2);
    eL_Phase_FroR(j)     = eL_Phase_Fro(j) / tensor_norm(angle_Sjj, 2);
    JS_Phase_FroR(j)     = JS_Phase_Fro(j) / tensor_norm(angle_Sjj, 2); 
    
    % L1 Norm Comparison for Phase
    XA_Phase_L1(j)       = tensor_norm(angle_XA - angle_Sjj, 1);
    mn_Phase_L1(j)       = tensor_norm(angle_mn - angle_Sjj, 1);
    eL_Phase_L1(j)       = tensor_norm(angle_eL - angle_Sjj, 1);
    JS_Phase_L1(j)       = tensor_norm(angle_JS - angle_Sjj, 1); 
    
    % Relative L1 Norm
    XA_Phase_L1R(j)      = XA_Phase_L1(j) / tensor_norm(angle_Sjj, 1);
    mn_Phase_L1R(j)      = mn_Phase_L1(j) / tensor_norm(angle_Sjj, 1);
    eL_Phase_L1R(j)      = eL_Phase_L1(j) / tensor_norm(angle_Sjj, 1);
    JS_Phase_L1R(j)      = JS_Phase_L1(j) / tensor_norm(angle_Sjj, 1); 
end


%% Visualization (Violin Plots)
import guide.Visualization.*
import guide.Visualization.DataVizm.*
import guide.Visualization.DataVizm.daviolinplot.*

colors = [0.2 0.6 0.8;   % eLORETA
          0.8 0.4 0.2;   % LCMV / MNE
          0.6 0.8 0.2;   % Xi-AlphaNET
          0.6 0.2 0.8;   % J-SPACE SPLINE [ADDED]
          0.5 0.5 0.5];  % gray

color_xa = colors(3, :);
color_foof = colors(1, :);
color_elxa = colors(2,:);
color_jspace = colors(4, :); 

% Overall MSE
MSE_xl = zeros(Nsub, Nroi);
MSE_mf = zeros(Nsub, Nroi);
MSE_lo = zeros(Nsub, Nroi);
MSE_js = zeros(Nsub, Nroi); 

for j = 1:Nsub
    for roi = 1:Nroi
        l_roi  = squeeze(l(roi,:,j));
        xl_roi = squeeze(xl(roi,:,j));
        mf_roi = squeeze(mf(roi,:,j));
        lo_roi = squeeze(el_fit_xa(roi,:,j));
        js_roi = squeeze(js(roi,:,j)); 
        
        MSE_xl(j, roi) = mean((xl_roi - l_roi).^2, 'omitnan');
        MSE_mf(j, roi) = mean((mf_roi - l_roi).^2, 'omitnan');
        MSE_lo(j, roi) = mean((lo_roi - l_roi).^2, 'omitnan');
        MSE_js(j, roi) = mean((js_roi - l_roi).^2, 'omitnan'); 
    end
end

mean_MSE_xl = mean(MSE_xl, 2, 'omitnan');
mean_MSE_mf = mean(MSE_mf, 2, 'omitnan');
mean_MSE_lo = mean(MSE_lo, 2, 'omitnan');
mean_MSE_js = mean(MSE_js, 2, 'omitnan'); 

% Paired t-test
[d_h, p_ttest, ci, stats] = ttest(mean_MSE_xl, mean_MSE_js);  
d = mean_MSE_xl - mean_MSE_js;
fprintf('\n[Paired t-test] H=%d p=%.6g | t=%.3f | df=%d\n', d_h, p_ttest, stats.tstat, stats.df);
fprintf('[Mean diff] mean(xl - js)=%.3e | std=%.3e\n', mean(d,'omitnan'), std(d,'omitnan'));

if mean(d,'omitnan') < 0
    fprintf('[Direction] Xi-AlphaNET has LOWER MSE than J-SPACE on average.\n');
else
    fprintf('[Direction] J-SPACE has LOWER MSE than Xi-AlphaNET on average.\n');
end

% Band-specific MSEs
freq_band_indices = cell(Nbands,1);
for b = 1:Nbands
    freq_band_indices{b} = find(freq >= bands{b,2} & freq <= bands{b,3});
end

MSE_xl_band = zeros(Nsub, Nroi, Nbands);
MSE_mf_band = zeros(Nsub, Nroi, Nbands);
MSE_eL_band = zeros(Nsub, Nroi, Nbands);
MSE_js_band = zeros(Nsub, Nroi, Nbands); 

for j = 1:Nsub
    for roi = 1:Nroi
        l_roi  = squeeze(l(roi,:,j));
        xl_roi = squeeze(xl(roi,:,j));
        mf_roi = squeeze(mf(roi,:,j));
        lo_roi = squeeze(el_fit_xa(roi,:,j));
        js_roi = squeeze(js(roi,:,j)); 
        
        for b = 1:Nbands
            idx = freq_band_indices{b};
            MSE_xl_band(j, roi, b) = mean((xl_roi(idx)-l_roi(idx)).^2, 'omitnan');
            MSE_mf_band(j, roi, b) = mean((mf_roi(idx)-l_roi(idx)).^2, 'omitnan');
            MSE_eL_band(j, roi, b) = mean((lo_roi(idx)-l_roi(idx)).^2, 'omitnan');
            MSE_js_band(j, roi, b) = mean((js_roi(idx)-l_roi(idx)).^2, 'omitnan'); 
        end
    end
end

mean_MSE_xl_band = squeeze(mean(MSE_xl_band,2,'omitnan'));
mean_MSE_mf_band = squeeze(mean(MSE_mf_band,2,'omitnan'));
mean_MSE_eL_band = squeeze(mean(MSE_eL_band,2,'omitnan'));
mean_MSE_js_band = squeeze(mean(MSE_js_band,2,'omitnan'));

% Overall violin plot
errors_combined = {mean_MSE_xl, mean_MSE_mf, mean_MSE_lo, mean_MSE_js};
methods = {'\xi-\alphaNET', 'MNE+FOOOF', 'eLORETA+\xi-\alpha', 'J-SPACE Spline'};
colors_overall = [color_xa; color_foof; color_elxa; color_jspace];

figure('Color','w');
daviolinplot(errors_combined, 'violin', 'half', 'box', 3, ...
    'xtlabels', methods, 'scatter', 0, 'violinalpha', 0.7, ...
    'colors', colors_overall, 'boxwidth', 1.5, 'violinwidth', 1.2, ...
    'boxcolors', 'w', 'boxspacing', 1.2);
ylabel('Mean Squared Error','FontSize',14,'FontWeight','bold');
title('Overall Reconstruction Errors','FontSize',18,'FontWeight','bold');
grid on;
set(gca,'FontSize',14);

% Plot band-specific violin plots
selected_bands = 1:(Nbands-1); % Exclude Gamma if desired
errors_band = {};
xtlabels_band = {};
colors_band = [];

for idx = 1:numel(selected_bands)
    b = selected_bands(idx);
    errors_band{end+1} = log10(mean_MSE_xl_band(:,b));
    errors_band{end+1} = log10(mean_MSE_mf_band(:,b));
    errors_band{end+1} = log10(mean_MSE_eL_band(:,b));
    errors_band{end+1} = log10(mean_MSE_js_band(:,b)); 
    
    band_name = bands{b,1};
    if iscell(band_name), band_name = band_name{1}; end
    if isstring(band_name), band_name = char(band_name); end
    
    xtlabels_band{end+1} = [band_name, ' (\xi-\alpha)'];
    xtlabels_band{end+1} = [band_name, ' (MNE)'];
    xtlabels_band{end+1} = [band_name, ' (eL+\xi)'];
    xtlabels_band{end+1} = [band_name, ' (J-SPACE)']; 
    
    colors_band = [colors_band; color_xa; color_foof; color_elxa; color_jspace];
end

figure('Color','w');
daviolinplot(errors_band, 'violin', 'half', 'box', 3, ...
    'xtlabels', xtlabels_band, 'scatter', 0, 'violinalpha', 0.7, ...
    'colors', colors_band, 'boxwidth', 1.5, 'violinwidth', 1.2, ...
    'boxcolors', 'w', 'boxspacing', 1.2);
ylabel('Mean Squared Error (log10)','FontSize',14,'FontWeight','bold');
title('Band-Specific Reconstruction Errors','FontSize',18,'FontWeight','bold');
grid on;
set(gca,'FontSize',14);
xtickangle(45);

%% Coherence Performance
if size(XA_Coherence_FroR,2) > 1, XA_Coherence_FroR_mean = mean(XA_Coherence_FroR, 2, 'omitnan'); else, XA_Coherence_FroR_mean = XA_Coherence_FroR; end
if size(mn_Coherence_FroR,2) > 1, mn_Coherence_FroR_mean = mean(mn_Coherence_FroR, 2, 'omitnan'); else, mn_Coherence_FroR_mean = mn_Coherence_FroR; end
if size(eL_Coherence_FroR,2) > 1, eL_Coherence_FroR_mean = mean(eL_Coherence_FroR, 2, 'omitnan'); else, eL_Coherence_FroR_mean = eL_Coherence_FroR; end
if size(JS_Coherence_FroR,2) > 1, JS_Coherence_FroR_mean = mean(JS_Coherence_FroR, 2, 'omitnan'); else, JS_Coherence_FroR_mean = JS_Coherence_FroR; end

errors_coherence = {log10(XA_Coherence_FroR_mean), log10(mn_Coherence_FroR_mean), log10(eL_Coherence_FroR_mean), log10(JS_Coherence_FroR_mean)};
methods_coherence = {'\xi-\alphaNET', 'MNE+FOOOF', 'eLORETA+\xi-\alpha', 'J-SPACE Spline'};
colors_coherence = [color_xa; color_foof; color_elxa; color_jspace];

figure('Color','w');
daviolinplot(errors_coherence, 'violin', 'half', 'box', 3, ...
    'xtlabels', methods_coherence, 'scatter', 0, 'violinalpha', 0.7, ...
    'colors', colors_coherence, 'boxwidth', 1.5, 'violinwidth', 1.2, ...
    'boxcolors', 'w', 'boxspacing', 1.2);
ylabel('log_{10}(Relative Frobenius Norm)','FontSize',14,'FontWeight','bold');
title('Performance Comparison (Coherence)','FontSize',18,'FontWeight','bold');
grid on;
set(gca,'FontSize',14);

%% Phase Performance
if size(XA_Phase_FroR,2) > 1, XA_Phase_FroR_mean = mean(XA_Phase_FroR, 2, 'omitnan'); else, XA_Phase_FroR_mean = XA_Phase_FroR; end
if size(mn_Phase_FroR,2) > 1, mn_Phase_FroR_mean = mean(mn_Phase_FroR, 2, 'omitnan'); else, mn_Phase_FroR_mean = mn_Phase_FroR; end
if size(eL_Phase_FroR,2) > 1, eL_Phase_FroR_mean = mean(eL_Phase_FroR, 2, 'omitnan'); else, eL_Phase_FroR_mean = eL_Phase_FroR; end
if size(JS_Phase_FroR,2) > 1, JS_Phase_FroR_mean = mean(JS_Phase_FroR, 2, 'omitnan'); else, JS_Phase_FroR_mean = JS_Phase_FroR; end

errors_phase = {log10(XA_Phase_FroR_mean), log10(mn_Phase_FroR_mean), log10(eL_Phase_FroR_mean), log10(JS_Phase_FroR_mean)};

figure('Color','w');
daviolinplot(errors_phase, 'violin', 'half', 'box', 3, ...
    'xtlabels', methods_coherence, 'scatter', 0, 'violinalpha', 0.7, ...
    'colors', colors_coherence, 'boxwidth', 1.5, 'violinwidth', 1.2, ...
    'boxcolors', 'w', 'boxspacing', 1.2);
ylabel('log_{10}(Relative Frobenius Norm)','FontSize',14,'FontWeight','bold');
title('Performance Comparison (Phase)','FontSize',18,'FontWeight','bold');
grid on;
set(gca,'FontSize',14);