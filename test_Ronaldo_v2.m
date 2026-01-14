


%%
% Import required helper functions
clear all;
close all;
import functions.auxx.ModelVectorization.*;
import guide.Visualization.*;
import functions.auxx.ZeroInflatedModels.*;
import functions.auxx.Refine_Solution.*;
import functions.auxx.OptimizedOperations.*;
Cortex = load("templates/Cortex.mat");
rng(2026);
% Path to the JSON file with model result metadata. Modify this directions
% manually acording to the location of the downloaded data 
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
js = zeros(Nroi, Nw, Nsub); % [ADDED] J-SPACE Power

% Frobenius and L1 norm benchmark errors for each subject
XA_Spectrum_Fro = zeros(Nsub, 1);
XA_Spectrum_FroR = zeros(Nsub, 1);
XA_Spectrum_L1 = zeros(Nsub, 1);
XA_Spectrum_L1R = zeros(Nsub, 1);

mn_Spectrum_Fro = zeros(Nsub, 1);
mn_Spectrum_FroR = zeros(Nsub, 1);
mn_Spectrum_L1 = zeros(Nsub, 1);
mn_Spectrum_L1R = zeros(Nsub, 1);

% [ADDED] J-SPACE Spectrum Metrics
JS_Spectrum_Fro = zeros(Nsub, 1);
JS_Spectrum_FroR = zeros(Nsub, 1);
JS_Spectrum_L1 = zeros(Nsub, 1);
JS_Spectrum_L1R = zeros(Nsub, 1);

XA_Coherence_Fro = zeros(Nsub, 1);
XA_Coherence_FroR = zeros(Nsub, 1);
XA_Coherence_L1 = zeros(Nsub, 1);
XA_Coherence_L1R = zeros(Nsub, 1);

mn_Coherence_Fro = zeros(Nsub, 1);
mn_Coherence_FroR = zeros(Nsub, 1);
mn_Coherence_L1 = zeros(Nsub, 1);
mn_Coherence_L1R = zeros(Nsub, 1);

% [ADDED] J-SPACE Coherence Metrics
JS_Coherence_Fro = zeros(Nsub, 1);
JS_Coherence_FroR = zeros(Nsub, 1);
JS_Coherence_L1 = zeros(Nsub, 1);
JS_Coherence_L1R = zeros(Nsub, 1);

XA_Phase_Fro = zeros(Nsub, 1);
XA_Phase_FroR = zeros(Nsub, 1);
XA_Phase_L1 = zeros(Nsub, 1);
XA_Phase_L1R = zeros(Nsub, 1);

mn_Phase_Fro = zeros(Nsub, 1);
mn_Phase_FroR = zeros(Nsub, 1);
mn_Phase_L1 = zeros(Nsub, 1);
mn_Phase_L1R = zeros(Nsub, 1);

% [ADDED] J-SPACE Phase Metrics
JS_Phase_Fro = zeros(Nsub, 1);
JS_Phase_FroR = zeros(Nsub, 1);
JS_Phase_L1 = zeros(Nsub, 1);
JS_Phase_L1R = zeros(Nsub, 1);

% ============================================================
% [J-SPACE] Warm-start state across subjects (DO NOT TOUCH OTHERS)
% ============================================================
eta_prev = [];          % 1x3 best eta from previous subject
eta_bank = [];          % Kx3 bank of good etas
eta_span_decades = 1.0; % +/- decades around eta_prev (1.0 => ±1 decade)
eta_bank_keep = 8;      % keep last K etas (simple & robust)


% Simulation loop
for j = 1:Nsub
    tic;
    % Load subject data
    subject_folder = selected_folders(j).name;
    mat_file_path = fullfile(dir_data, subject_folder, [subject_folder, '.mat']);
    data_struct = load(mat_file_path);

    freq = data_struct.data_struct.freqrange(1:Nw);
    parameters.Data.freq = freq;

    % Source cross-spectrum simulation
    %%
    disp('->> Estimating Source Cross Using a Neural Mass Simulation');
    [Sjj, ~] = functions.auxx.Simulations.neural_mass_simulation(0, json_path);

    % Generate scalp-level cross-spectrum with noise

    % N_wishart = 4;
    % disp('->> Simulating Scalp Cross Wishart Noise + FModel');
    % parfor i = 1:Nw
    %     Sjj_cross(:, :, i) = generate_complex_wishart(Sjj(:, :, i), N_wishart);
    %     Svv_cross(:, :, i) = L * Sjj_cross(:, :, i) * L';  % Apply transformation
    % end
    N_wishart = 2 * size(L, 1); % 自由度至少要大于感应器数量，保证满秩
% 或者保持 N_wishart 小，但必须加噪声

parfor i = 1:Nw
    Sjj_cross(:, :, i) = generate_complex_wishart(Sjj(:, :, i), N_wishart);
    
    % 投影
    Svv_clean = L * Sjj_cross(:, :, i) * L';
    
    % 【关键】加上感应器白噪声 (Sensor Noise)
    % 这一步能保证矩阵满秩，白化后不会全为 1
    noise_level = 0.05 * trace(Svv_clean) / size(L, 1);
    Svv_cross(:, :, i) = Svv_clean + noise_level * eye(size(L, 1)); 
end
    [Log_Spec,freq]= log_spectrum(Svv_cross,freq);
    [n, ~, F] = size(Svv_cross);

diag_mask = repmat(logical(eye(n)), [1, 1, F]);

Svv_cross(diag_mask) = real(Svv_cross(diag_mask));
    % plot(freq,Log_Spec')

    %
    toc;

    % Xi-AlphaNET estimation
    data.Cross = Svv_cross;
    data.age = 25;  % Subject age
    data.freq = freq;

    disp('->> Xi-AlphaNeT Inverse Solution');
    [x, ~, G, x0] = Xi_ALphaNET(properties, data, parameters);
    source_act_cross = functions.auxx.CrossSpectrum.eval_source_conn(x.Solution, data.freq, parameters.Model.R, properties, parameters);
    XA_Sjj_cross = source_act_cross.Cross.Full;

    % Additional processing (mean cross-spectrum)
    % Calculate the mean cross-spectrum across frequencies
    mn_Sjj_cross = mn_cross(Svv_cross, K, 0);

    disp('->> eLORETA Processing...');
    % Compute eLORETA source cross-spectra for each frequency bin
    parfor i = 1:Nw
        source = inverse(Svv_cross(:,:,i), L);
        eL_Sjj_cross(:,:,i) = source.eloreata.Sjj;
    end

% [ADDED] J-SPACE Processing (Plan-1: DWI soft prior + eta-search + dynamic lambdas)
disp('->> J-SPACE Processing (Global 3D Surrogate Opt)...');

% --- Build cfg for J-SPACE ---
cfg_js = struct();
cfg_js.max_em_iter   = 10;
cfg_js.verbose       = true;
cfg_js.use_gpu       = false;
cfg_js.debug_print   = false;

% --- Optimizer settings ---
cfg_js.opt_max_evals = 50;
cfg_js.opt_min_points = 10;
cfg_js.opt_use_parallel = true;

% --- Objective / EBIC sample size ---
cfg_js.m_samples     = 1000;

% --- M-step solver choice ---
cfg_js.optimizer     = 'stoch_fista';   % [CHANGED] new solver uses stochastic FISTA
cfg_js.use_fista     = false;           % [ADDED] avoid legacy override if someone sets it elsewhere

% --- (PLAN-1) Pass DWI connectivity as soft prior ---
cfg_js.dwi_C = parameters.Compact_Model.C;
cfg_js.freq  = parameters.Data.freq;

cfg_js.dwi_weight_mode   = 'power';
cfg_js.dwi_weight_alpha  = 1.0;
cfg_js.dwi_weight_clip   = [0.25, 4];
cfg_js.dwi_weight_normalize_median = true;

cfg_js.postprocess_enable = false;

% ============================================================
% [STOCH CONFIG] (Added, no deletions)
%   This matches helper functions:
%   - stoch_prepare_freq_sampler_
%   - stoch_sample_freq_indices_
%   - stoch_sub_kernel_rescale_
% ============================================================
if ~isfield(cfg_js,'stoch') || isempty(cfg_js.stoch)
    cfg_js.stoch = struct();
end

cfg_js.stoch.mode = 'B';                 % keep your setting ('B' = banded)
cfg_js.stoch.seed = 0;                   % [ADDED] RNG seed base (used with iter to be reproducible)

% How many freqs sampled per band per stochastic M-step
% (your helper supports either 'n_per_band' or 'Nsfreq'; we set both for compatibility)
cfg_js.stoch.n_per_band = 2;             % [ADDED] suggested default (tune later)
cfg_js.stoch.Nsfreq     = cfg_js.stoch.n_per_band;  % [ADDED] alias

% Optional cap on total sampled freqs each iteration
cfg_js.stoch.max_total  = 16;            % [ADDED] safety cap (<=47)

% EEG-like band edges (Hz). You can customize.
% For your freq range ~[1.17, 19.14], these edges work well.
cfg_js.stoch.band_edges = [0 4 8 13 20]; % [ADDED] last edge <= max(freq)

% Stochastic FISTA inner loop (fixed-iter as you requested)
cfg_js.stoch.max_iter   = 40;            % [ADDED] fixed iter count for stochastic FISTA
cfg_js.stoch.tol        = 0;             % [ADDED] disable tol-stop (fixed iter)
cfg_js.stoch.backtracking_beta = 0.5;    % [ADDED] match your fista settings style
cfg_js.stoch.max_backtracking  = 20;     % [ADDED]
cfg_js.stoch.monotone          = true;   % [ADDED]
cfg_js.stoch.use_restart       = true;   % [ADDED]

% Optional: kernel mass rescale when using subset of freqs
cfg_js.stoch.rescale_kernel_mass = true; % [ADDED]

% ============================================================
% [J-SPACE MOD-1] Tighten eta bounds (hard safety + warm-start bounds)
% ============================================================
eta_global_lb = [3e-4, 3e-4, 3e-4];
eta_global_ub = [3e-2, 3e-1, 5e-1];     % NOTE: eta2_ub tightened to 0.3

if ~isempty(eta_prev)
    cen_log = log10(eta_prev(:)');
    lb = 10.^(cen_log - eta_span_decades);
    ub = 10.^(cen_log + eta_span_decades);

    lb = max(lb, eta_global_lb);
    ub = min(ub, eta_global_ub);
else
    lb = eta_global_lb;
    ub = eta_global_ub;
end

cfg_js.eta1_lb = lb(1);  cfg_js.eta1_ub = ub(1);
cfg_js.eta2_lb = lb(2);  cfg_js.eta2_ub = ub(2);
cfg_js.eta3_lb = lb(3);  cfg_js.eta3_ub = ub(3);

% ============================================================
% [J-SPACE MOD-2] Provide MULTI initial points to surrogateopt
% ============================================================
init_eta = [];

eta0_default = [0.1, 0.2, 0.1];
eta0_default = min(max(eta0_default, lb), ub);
init_eta = [init_eta; eta0_default];

if ~isempty(eta_prev)
    init_eta = [init_eta; min(max(eta_prev, lb), ub)];
end

if ~isempty(eta_bank)
    k = min(size(eta_bank,1), 6);
    cand = eta_bank(end-k+1:end, :);
    cand = min(max(cand, lb), ub);
    init_eta = [init_eta; cand];
end

need_extra = max(0, cfg_js.opt_min_points - size(init_eta,1));
if need_extra > 0
    try
        X = lhsdesign(need_extra, 3);
    catch
        X = rand(need_extra, 3);
    end
    lb_log = log10(lb); ub_log = log10(ub);
    extra_eta = 10.^(lb_log + X .* (ub_log - lb_log));
    init_eta = [init_eta; extra_eta];
end

init_eta = unique(round(init_eta, 12), 'rows', 'stable');

cfg_js.opt_initial_points_log = log10(init_eta);

cfg_js.eta1_init = init_eta(1,1);
cfg_js.eta2_init = init_eta(1,2);
cfg_js.eta3_init = init_eta(1,3);

% --- Call the solver ---
% [Omega_est, JS_Sjj_cross, outs_js] = solver_jspace_3d_opt(Svv_cross, L, [], cfg_js);
% [Omega_est, JS_Sjj_cross, outs_js] = solver_jspace_3d_opt_bayesopt(Svv_cross, L, [], cfg_js);
[Omega_est, JS_Sjj_cross, outs_js] = solver_jspace_3d_opt_stoch(Svv_cross, L, [], cfg_js);

% ============================================================
% [J-SPACE] Update warm-start state for next subject
% ============================================================
if isfield(outs_js, 'best_eta')
    eta_prev = [outs_js.best_eta.eta1, outs_js.best_eta.eta2, outs_js.best_eta.eta3];
else
    % fallback if field name differs
    try
        eta_prev = [outs_js.global_hyperparams.eta1, outs_js.global_hyperparams.eta2, outs_js.global_hyperparams.eta3];
    catch
        eta_prev = [];
    end
end
if ~isempty(eta_prev)
    eta_bank = [eta_bank; eta_prev];
    if size(eta_bank,1) > eta_bank_keep
        eta_bank = eta_bank(end-eta_bank_keep+1:end, :);
    end
end

% --- (Optional) Print results ---
if cfg_js.verbose
    fprintf('\n[J-SPACE] 3D Optimization Summary:\n');
    if isfield(outs_js,'best_eta')
        fprintf('  > BEST eta: eta1=%.3e | eta2=%.3e | eta3=%.3e\n', ...
            outs_js.best_eta.eta1, outs_js.best_eta.eta2, outs_js.best_eta.eta3);
    end
    if isfield(outs_js,'em_lambdas') && ~isempty(outs_js.em_lambdas)
        lam_last = outs_js.em_lambdas(end,:);
        fprintf('  > Final EM lambdas: [%.3e %.3e %.3e]\n', lam_last(1), lam_last(2), lam_last(3));
    end
    if isfield(outs_js,'opt_results') && isfield(outs_js.opt_results,'fval')
        fprintf('  > Best AIC+Penalty: %.3e\n', outs_js.opt_results.fval);
    end
end


    % Initialize mean power array for each ROI and frequency
    mn_power = zeros(Nroi, Nw);

    % Extract and log-transform diagonal power for each ROI
    for i = 1:Nw
        xl(:, i, j)= log10(real(diag(XA_Sjj_cross(:, :, i))));
        % Original cross-spectrum power (log-transformed)
        l(:, i, j) = log10(real(diag(Sjj_cross(:, :, i))));
        % eLORETA power (not log-transformed)
        lo(:, i) = real(diag(eL_Sjj_cross(:, :, i)));
        % J-SPACE power [ADDED]
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
    % Compare estimated cross-spectra to ground-truth using Frobenius and L1 norms

    % Notation:
    % - "XA" corresponds to your Xi-AlphaNET or a different model; rename as needed
    % - "mn" corresponds to the mean cross-spectrum
    % - "eL" corresponds to the eLORETA estimates

    % Frobenius Norm Comparison
    XA_Spectrum_Fro(j)   = tensor_norm(XA_Sjj_cross - ground_truth_Sjj, 2);
    mn_Spectrum_Fro(j)   = tensor_norm(mn_Sjj_cross - ground_truth_Sjj, 2);
    eL_Spectrum_Fro(j)   = tensor_norm(eL_Sjj_cross - ground_truth_Sjj, 2);
    JS_Spectrum_Fro(j)   = tensor_norm(JS_tensor - ground_truth_Sjj, 2); % [ADDED]

    % Relative Frobenius Norm (normalized by ground-truth)
    XA_Spectrum_FroR(j)  = XA_Spectrum_Fro(j) / tensor_norm(ground_truth_Sjj, 2);
    mn_Spectrum_FroR(j)  = mn_Spectrum_Fro(j) / tensor_norm(ground_truth_Sjj, 2);
    eL_Spectrum_FroR(j)  = eL_Spectrum_Fro(j) / tensor_norm(ground_truth_Sjj, 2);
    JS_Spectrum_FroR(j)  = JS_Spectrum_Fro(j) / tensor_norm(ground_truth_Sjj, 2); % [ADDED]

    % L1 Norm Comparison
    XA_Spectrum_L1(j)    = tensor_norm(XA_Sjj_cross - ground_truth_Sjj, 1);
    mn_Spectrum_L1(j)    = tensor_norm(mn_Sjj_cross - ground_truth_Sjj, 1);
    eL_Spectrum_L1(j)    = tensor_norm(eL_Sjj_cross - ground_truth_Sjj, 1);
    JS_Spectrum_L1(j)    = tensor_norm(JS_tensor - ground_truth_Sjj, 1); % [ADDED]

    % Relative L1 Norm (normalized by ground-truth)
    XA_Spectrum_L1R(j)   = XA_Spectrum_L1(j) / tensor_norm(ground_truth_Sjj, 1);
    mn_Spectrum_L1R(j)   = mn_Spectrum_L1(j) / tensor_norm(ground_truth_Sjj, 1);
    eL_Spectrum_L1R(j)   = eL_Spectrum_L1(j) / tensor_norm(ground_truth_Sjj, 1);
    JS_Spectrum_L1R(j)   = JS_Spectrum_L1(j) / tensor_norm(ground_truth_Sjj, 1); % [ADDED]

    % ----- Coherence Comparison -----
    % Compare coherence matrices for each model
    coherence_Sjj = coherence(ground_truth_Sjj);
    coherence_XA  = coherence(XA_Sjj_cross);
    coherence_mn  = coherence(mn_Sjj_cross);
    coherence_eL  = coherence(eL_Sjj_cross);
    coherence_JS  = coherence(JS_tensor); % [ADDED]

    % Frobenius Norm Comparison for Coherence
    XA_Coherence_Fro(j)  = tensor_norm(coherence_XA - coherence_Sjj, 2);
    mn_Coherence_Fro(j)  = tensor_norm(coherence_mn - coherence_Sjj, 2);
    eL_Coherence_Fro(j)  = tensor_norm(coherence_eL - coherence_Sjj, 2);
    JS_Coherence_Fro(j)  = tensor_norm(coherence_JS - coherence_Sjj, 2); % [ADDED]

    % Relative Frobenius Norm
    XA_Coherence_FroR(j) = XA_Coherence_Fro(j) / tensor_norm(coherence_Sjj, 2);
    mn_Coherence_FroR(j) = mn_Coherence_Fro(j) / tensor_norm(coherence_Sjj, 2);
    eL_Coherence_FroR(j) = eL_Coherence_Fro(j) / tensor_norm(coherence_Sjj, 2);
    JS_Coherence_FroR(j) = JS_Coherence_Fro(j) / tensor_norm(coherence_Sjj, 2); % [ADDED]

    % L1 Norm Comparison for Coherence
    XA_Coherence_L1(j)   = tensor_norm(coherence_XA - coherence_Sjj, 1);
    mn_Coherence_L1(j)   = tensor_norm(coherence_mn - coherence_Sjj, 1);
    eL_Coherence_L1(j)   = tensor_norm(coherence_eL - coherence_Sjj, 1);
    JS_Coherence_L1(j)   = tensor_norm(coherence_JS - coherence_Sjj, 1); % [ADDED]

    % Relative L1 Norm
    XA_Coherence_L1R(j)  = XA_Coherence_L1(j) / tensor_norm(coherence_Sjj, 1);
    mn_Coherence_L1R(j)  = mn_Coherence_L1(j) / tensor_norm(coherence_Sjj, 1);
    eL_Coherence_L1R(j)  = eL_Coherence_L1(j) / tensor_norm(coherence_Sjj, 1);
    JS_Coherence_L1R(j)  = JS_Coherence_L1(j) / tensor_norm(coherence_Sjj, 1); % [ADDED]

    % ----- Phase Comparison -----
    % Compare phase angles across models
    angle_Sjj = angle(ground_truth_Sjj);
    angle_XA  = angle(XA_Sjj_cross);
    angle_mn  = angle(mn_Sjj_cross);
    angle_eL  = angle(eL_Sjj_cross);
    angle_JS  = angle(JS_tensor); % [ADDED]

    % Frobenius Norm Comparison for Phase
    XA_Phase_Fro(j)      = tensor_norm(angle_XA - angle_Sjj, 2);
    mn_Phase_Fro(j)      = tensor_norm(angle_mn - angle_Sjj, 2);
    eL_Phase_Fro(j)      = tensor_norm(angle_eL - angle_Sjj, 2);
    JS_Phase_Fro(j)      = tensor_norm(angle_JS - angle_Sjj, 2); % [ADDED]

    % Relative Frobenius Norm
    XA_Phase_FroR(j)     = XA_Phase_Fro(j) / tensor_norm(angle_Sjj, 2);
    mn_Phase_FroR(j)     = mn_Phase_Fro(j) / tensor_norm(angle_Sjj, 2);
    eL_Phase_FroR(j)     = eL_Phase_Fro(j) / tensor_norm(angle_Sjj, 2);
    JS_Phase_FroR(j)     = JS_Phase_Fro(j) / tensor_norm(angle_Sjj, 2); % [ADDED]

    % L1 Norm Comparison for Phase
    XA_Phase_L1(j)       = tensor_norm(angle_XA - angle_Sjj, 1);
    mn_Phase_L1(j)       = tensor_norm(angle_mn - angle_Sjj, 1);
    eL_Phase_L1(j)       = tensor_norm(angle_eL - angle_Sjj, 1);
    JS_Phase_L1(j)       = tensor_norm(angle_JS - angle_Sjj, 1); % [ADDED]

    % Relative L1 Norm
    XA_Phase_L1R(j)      = XA_Phase_L1(j) / tensor_norm(angle_Sjj, 1);
    mn_Phase_L1R(j)      = mn_Phase_L1(j) / tensor_norm(angle_Sjj, 1);
    eL_Phase_L1R(j)      = eL_Phase_L1(j) / tensor_norm(angle_Sjj, 1);
    JS_Phase_L1R(j)      = JS_Phase_L1(j) / tensor_norm(angle_Sjj, 1); % [ADDED]
end

% -----------------------------------
%% Colors
import guide.Visualization.*
import guide.Visualization.DataVizm.*
import guide.Visualization.DataVizm.daviolinplot.*
colors = [0.2 0.6 0.8;   % eLORETA
          0.8 0.4 0.2;   % LCMV / MNE
          0.6 0.8 0.2;   % Xi-AlphaNET
          0.6 0.2 0.8;   % J-SPACE [ADDED]
          0.5 0.5 0.5];  % gray

color_xa = colors(3, :);
color_foof = colors(1, :);
color_elxa = colors(2,:);
color_jspace = colors(4, :); % [ADDED]

% Overall MSE
MSE_xl = zeros(Nsub, Nroi);
MSE_mf = zeros(Nsub, Nroi);
MSE_lo = zeros(Nsub, Nroi);
MSE_js = zeros(Nsub, Nroi); % [ADDED]

for j = 1:Nsub
    for roi = 1:Nroi
        l_roi  = squeeze(l(roi,:,j));
        xl_roi = squeeze(xl(roi,:,j));
        mf_roi = squeeze(mf(roi,:,j));
        lo_roi = squeeze(el_fit_xa(roi,:,j));
        js_roi = squeeze(js(roi,:,j)); % [ADDED]

        MSE_xl(j, roi) = mean((xl_roi - l_roi).^2, 'omitnan');
        MSE_mf(j, roi) = mean((mf_roi - l_roi).^2, 'omitnan');
        MSE_lo(j, roi) = mean((lo_roi - l_roi).^2, 'omitnan');
        MSE_js(j, roi) = mean((js_roi - l_roi).^2, 'omitnan'); % [ADDED]
    end
end

mean_MSE_xl = mean(MSE_xl, 2, 'omitnan');
mean_MSE_mf = mean(MSE_mf, 2, 'omitnan');
mean_MSE_lo = mean(MSE_lo, 2, 'omitnan');
mean_MSE_js = mean(MSE_js, 2, 'omitnan'); % [ADDED]

% Paired t-test
[d_h, p_ttest, ci, stats] = ttest(mean_MSE_xl, mean_MSE_js);  % paired t-test
d = mean_MSE_xl - mean_MSE_js;

fprintf('\n[Paired t-test] H=%d p=%.6g | t=%.3f | df=%d\n', d_h, p_ttest, stats.tstat, stats.df);
fprintf('[Mean diff] mean(xl - mf)=%.3e | std=%.3e\n', mean(d,'omitnan'), std(d,'omitnan'));
fprintf('[95%% CI] [%.3e, %.3e] for mean(xl - mf)\n', ci(1), ci(2));

if mean(d,'omitnan') < 0
    fprintf('[Direction] Xi-AlphaNET (xl) has LOWER MSE than JSAPCE on average.\n');
else
    fprintf('[Direction] jsapce has LOWER MSE than Xi-AlphaNET (xl) on average.\n');
end


% Band-specific MSEs
freq_band_indices = cell(Nbands,1);
for b = 1:Nbands
    freq_band_indices{b} = find(freq >= bands{b,2} & freq <= bands{b,3});
end

MSE_xl_band = zeros(Nsub, Nroi, Nbands);
MSE_mf_band = zeros(Nsub, Nroi, Nbands);
MSE_eL_band = zeros(Nsub, Nroi, Nbands);
MSE_js_band = zeros(Nsub, Nroi, Nbands); % [ADDED]

for j = 1:Nsub
    for roi = 1:Nroi
        l_roi  = squeeze(l(roi,:,j));
        xl_roi = squeeze(xl(roi,:,j));
        mf_roi = squeeze(mf(roi,:,j));
        lo_roi = squeeze(el_fit_xa(roi,:,j));
        js_roi = squeeze(js(roi,:,j)); % [ADDED]
        
        for b = 1:Nbands
            idx = freq_band_indices{b};
            MSE_xl_band(j, roi, b) = mean((xl_roi(idx)-l_roi(idx)).^2, 'omitnan');
            MSE_mf_band(j, roi, b) = mean((mf_roi(idx)-l_roi(idx)).^2, 'omitnan');
            MSE_eL_band(j, roi, b) = mean((lo_roi(idx)-l_roi(idx)).^2, 'omitnan');
            MSE_js_band(j, roi, b) = mean((js_roi(idx)-l_roi(idx)).^2, 'omitnan'); % [ADDED]
        end
    end
end

mean_MSE_xl_band = squeeze(mean(MSE_xl_band,2,'omitnan'));
mean_MSE_mf_band = squeeze(mean(MSE_mf_band,2,'omitnan'));
mean_MSE_eL_band = squeeze(mean(MSE_eL_band,2,'omitnan'));
mean_MSE_js_band = squeeze(mean(MSE_js_band,2,'omitnan')); % [ADDED]

% Overall violin plot [MODIFIED]
errors_combined = {mean_MSE_xl, mean_MSE_mf, mean_MSE_lo, mean_MSE_js};
methods = {'\xi-\alphaNET', 'MNE+FOOOF', 'eLORETA+\xi-\alpha', 'J-SPACE'};
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

% Plot band-specific violin plots [MODIFIED]
selected_bands = 1:(Nbands-1); % Exclude Gamma
errors_band = {};
xtlabels_band = {};
colors_band = [];

for idx = 1:numel(selected_bands)
    b = selected_bands(idx);

    % Prepare cell array for daviolinplot
    errors_band{end+1} = log10(mean_MSE_xl_band(:,b));
    errors_band{end+1} = log10(mean_MSE_mf_band(:,b));
    errors_band{end+1} = log10(mean_MSE_eL_band(:,b));
    errors_band{end+1} = log10(mean_MSE_js_band(:,b)); % [ADDED]
    
    band_name = bands{b,1};
    if iscell(band_name), band_name = band_name{1}; end
    if isstring(band_name), band_name = char(band_name); end

    xtlabels_band{end+1} = [band_name, ' (\xi-\alpha)'];
    xtlabels_band{end+1} = [band_name, ' (MNE)'];
    xtlabels_band{end+1} = [band_name, ' (eL+\xi)'];
    xtlabels_band{end+1} = [band_name, ' (J-SPACE)']; % [ADDED]

    colors_band = [colors_band; color_xa; color_foof; color_elxa; color_jspace];
end

figure('Color','w');
daviolinplot(errors_band, 'violin', 'half', 'box', 3, ...
    'xtlabels', xtlabels_band, 'scatter', 0, 'violinalpha', 0.7, ...
    'colors', colors_band, 'boxwidth', 1.5, 'violinwidth', 1.2, ...
    'boxcolors', 'w', 'boxspacing', 1.2);
ylabel('Mean Squared Error (log10)','FontSize',14,'FontWeight','bold');
title('Band-Specific Reconstruction Errors (Gamma Excluded)','FontSize',18,'FontWeight','bold');
grid on;
set(gca,'FontSize',14);
xtickangle(45);
%%

% ----- Coherence Performance [MODIFIED] -----
if size(XA_Coherence_FroR,2) > 1
    XA_Coherence_FroR_mean = mean(XA_Coherence_FroR, 2, 'omitnan');
else
    XA_Coherence_FroR_mean = XA_Coherence_FroR;
end
XA_Coherence_FroR_log = log10(XA_Coherence_FroR_mean);

if size(mn_Coherence_FroR,2) > 1
    mn_Coherence_FroR_mean = mean(mn_Coherence_FroR, 2, 'omitnan');
else
    mn_Coherence_FroR_mean = mn_Coherence_FroR;
end
mn_Coherence_FroR_log = log10(mn_Coherence_FroR_mean);

if size(eL_Coherence_FroR,2) > 1
    eL_Coherence_FroR_mean = mean(eL_Coherence_FroR, 2, 'omitnan');
else
    eL_Coherence_FroR_mean = eL_Coherence_FroR;
end
eL_Coherence_FroR_log = log10(eL_Coherence_FroR_mean);

% [ADDED] J-SPACE Coherence
if size(JS_Coherence_FroR,2) > 1
    JS_Coherence_FroR_mean = mean(JS_Coherence_FroR, 2, 'omitnan');
else
    JS_Coherence_FroR_mean = JS_Coherence_FroR;
end
JS_Coherence_FroR_log = log10(JS_Coherence_FroR_mean);

errors_coherence = {XA_Coherence_FroR_log, mn_Coherence_FroR_log, eL_Coherence_FroR_log, JS_Coherence_FroR_log};
methods_coherence = {'\xi-\alphaNET', 'MNE+FOOOF', 'eLORETA+\xi-\alpha', 'J-SPACE'};
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

% ----- Phase Performance [MODIFIED] -----
if size(XA_Phase_FroR,2) > 1
    XA_Phase_FroR_mean = mean(XA_Phase_FroR, 2, 'omitnan');
else
    XA_Phase_FroR_mean = XA_Phase_FroR;
end
XA_Phase_FroR_log = log10(XA_Phase_FroR_mean);

if size(mn_Phase_FroR,2) > 1
    mn_Phase_FroR_mean = mean(mn_Phase_FroR, 2, 'omitnan');
else
    mn_Phase_FroR_mean = mn_Phase_FroR;
end
mn_Phase_FroR_log = log10(mn_Phase_FroR_mean);

if size(eL_Phase_FroR,2) > 1
    eL_Phase_FroR_mean = mean(eL_Phase_FroR, 2, 'omitnan');
else
    eL_Phase_FroR_mean = eL_Phase_FroR;
end
eL_Phase_FroR_log = log10(eL_Phase_FroR_mean);

% [ADDED] J-SPACE Phase
if size(JS_Phase_FroR,2) > 1
    JS_Phase_FroR_mean = mean(JS_Phase_FroR, 2, 'omitnan');
else
    JS_Phase_FroR_mean = JS_Phase_FroR;
end
JS_Phase_FroR_log = log10(JS_Phase_FroR_mean);

errors_phase = {XA_Phase_FroR_log, mn_Phase_FroR_log, eL_Phase_FroR_log, JS_Phase_FroR_log};
methods_phase = {'\xi-\alphaNET', 'MNE+FOOOF', 'eLORETA+\xi-\alpha', 'J-SPACE'};
colors_phase = [color_xa; color_foof; color_elxa; color_jspace];

figure('Color','w');
daviolinplot(errors_phase, 'violin', 'half', 'box', 3, ...
    'xtlabels', methods_phase, 'scatter', 0, 'violinalpha', 0.7, ...
    'colors', colors_phase, 'boxwidth', 1.5, 'violinwidth', 1.2, ...
    'boxcolors', 'w', 'boxspacing', 1.2);
ylabel('log_{10}(Relative Frobenius Norm)','FontSize',14,'FontWeight','bold');
title('Performance Comparison (Phase)','FontSize',18,'FontWeight','bold');
grid on;
set(gca,'FontSize',14);
