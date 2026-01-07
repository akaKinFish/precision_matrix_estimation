function Sjj_jspace = run_jspace_for_comparison(Svv_cross, L, freq)
% RUN_JSPACE_FOR_COMPARISON - Adapter to run J-SPACE in the comparison framework
%
% Inputs:
%   Svv_cross : [Ne x Ne x Nw] Sensor Cross-Spectrum
%   L         : [Ne x Nr] Leadfield Matrix
%   freq      : [1 x Nw] Frequency Vector
%
% Output:
%   Sjj_jspace : [Nr x Nr x Nw] Estimated Source Cross-Spectrum

[Ne, ~, Nw] = size(Svv_cross);
[~, Nr] = size(L);

% 1. Prepare Inputs for J-SPACE (Cell Array Format)
Svv_cell = cell(Nw, 1);
for k = 1:Nw
    Svv_cell{k} = Svv_cross(:, :, k);
end

% 2. Configure J-SPACE
cfg = struct();
cfg.max_em_iter   = 5;      % Fast convergence
cfg.grid_size     = 100;     % Optimal sparsity search
cfg.m_samples     = 1000;   % Effective samples heuristic
cfg.verbose       = true;  % Silent mode for batch run
cfg.use_gpu       = true;   % Enable GPU
cfg.lambda1_ratio = 1.0;
cfg.lambda3_ratio = 0.1;
cfg.enable_debias = true;
cfg.enable_rayleigh_search = true;
cfg.selection_metric = 'ebic';
cfg.ebic_gamma = 0;
cfg.hyper_grid_mode = 'scale_from_sbar';
cfg.alpha1_grid = [0.1 0.3 1];
cfg.alpha3_grid = [0 0.1];
cfg.inner_verbose = false;
cfg.plot = false;
% 建议包含 3.16 (sqrt(10)) 在内
cfg.rayleigh_range = 2.5:0.25:4.5;
% 3. Run Solver
% Note: J-SPACE returns Sigma_src_est (Covariance/Cross-Spectrum) and Omega_est (Precision)
try
    [~, Sigma_src_est, ~] = solver_jspace_hypersearch_global(Svv_cell, L, [], cfg);

    % 4. Convert Output back to 3D Matrix
    Sjj_jspace = zeros(Nr, Nr, Nw);
    for k = 1:Nw
        Sjj_jspace(:, :, k) = Sigma_src_est{k};
    end

catch ME
    warning('J-SPACE failed: %s. Returning zeros.', ME.message);
    Sjj_jspace = zeros(Nr, Nr, Nw);
end
end
