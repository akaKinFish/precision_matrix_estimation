function [] = InverseSolvers(output_source, paths, sens_system)
%% Run Inverse Solvers: HG-LASSO (Baseline) vs J-SPACE (Proposed)
%
% This script iterates through simulations, applies the inverse solvers,
% and saves the results.
%
% Solvers:
%   1. Higgs (h-hggm): The benchmark method.
%   2. J-SPACE: Joint Spatial-Spectral Precision And Connectivity Estimation.

%% 1. Loading Simulation Substrate or Real Data
% Assuming this file contains 'Svv_sim', 'Seeders_sim', 'LeadFields', 'Nsamp'
load('./data/Pseudorand_Net.mat')

%% 2. Prepare Head Model (Real Head Case)
if strcmp(sens_system,'real_head') == 1
    headmodel_data = load(paths.headmodel_mat);
    % Extract Leadfield (Gain)
    if isfield(headmodel_data, 'HeadModel') && isfield(headmodel_data.HeadModel, 'Gain')
        K = headmodel_data.HeadModel.Gain;
    elseif isfield(headmodel_data, 'Gain')
        K = headmodel_data.Gain;
    else
        error('Gain matrix not found in headmodel file');
    end
    % In simulation mode, LeadFields are usually loaded from .mat
end

%% 3. Select Subject / Leadfield
subject    = [1]; % User defined
LeadFields = LeadFields(1, subject);
Lvj0       = LeadFields{1};

%% ---------------------------------------------------------
%% Solver 1: Run h-hggm (Higgs) - Benchmark
%% ---------------------------------------------------------
fprintf('Running Higgs Solver...\n');
sol_higgs = InverseSolver_higgs(Svv_sim, LeadFields, Seeders_sim, Nsamp, sens_system);

% Create result directory if not exists
if ~exist('result', 'dir'), mkdir('result'); end
save(fullfile('result', 'Solutions_higgs.mat'), 'sol_higgs', '-v7.3');

%% ---------------------------------------------------------
%% Solver 2: Run J-SPACE (Proposed)
%% ---------------------------------------------------------
fprintf('Running J-SPACE Solver...\n');

try
    % Check if solver exists
    if exist('solver_jspace', 'file') == 2
        
        Nsim = size(Svv_sim, 2);
        q    = size(Seeders_sim, 1);
        p    = size(Lvj0, 1);
        
        % Initialize storage
        Theta_jspace_all = nan(q, q, Nsim);
        
        % J-SPACE Configuration
        cfg_jspace = struct();
        cfg_jspace.max_em_iter = 10;      % Sufficient for convergence
        cfg_jspace.grid_size   = 10;      % Lambda2 search resolution
        cfg_jspace.m_samples   = Nsamp;   % Crucial for EBIC
        cfg_jspace.verbose     = false;   % Keep loop clean
        cfg_jspace.use_gpu     = 1;       % Enable GPU if available
        
        % Optional: Graph Laplacian for Spatial Smoothing (lambda3)
        % If you have a spatial graph, load it here. Otherwise empty.
        GraphLaplacian = []; 

        % Loop over simulations
        for sim = 1:Nsim
            fprintf('  Sim %d/%d...', sim, Nsim);
            tic_sim = tic;
            
            try
                % 1. Prepare Data
                Svv_this = Svv_sim{1, sim}{1}; % Extract sensor covariance
                seeders  = Seeders_sim(:, sim);
                L_sel    = Lvj0(:, seeders);   % Selected Leadfield for active sources
                
                % J-SPACE expects cell array for Frequency (even if F=1)
                Svv_cell = {Svv_this};
                
                % 2. Run J-SPACE Solver
                % This handles E-step, Whitening, Auto-Tuning, PGD, Recoloring internally
                [Omega_est, ~, outs] = solver_jspace(Svv_cell, L_sel, GraphLaplacian, cfg_jspace);
                
                % 3. Extract Result
                % Omega_est is {F x 1}. We take the first freq.
                Theta_hat = Omega_est{1};
                
                % Force symmetry (numerical safety)
                Theta_hat = (Theta_hat + Theta_hat') / 2;
                
                % Store
                Theta_jspace_all(:, :, sim) = Theta_hat;
                
                % Optional: Log selected lambda
                best_lam = outs.best_lambdas(end);
                fprintf(' Done (%.2fs, lambda=%.2e)\n', toc(tic_sim), best_lam);
                
            catch MEi
                warning('\n  J-SPACE failed on sim %d: %s\n', sim, MEi.message);
            end
        end

        % Pack Results
        sol_jspace = struct();
        sol_jspace.name        = 'J-SPACE';
        sol_jspace.Theta       = Theta_jspace_all;
        sol_jspace.Seeders_sim = Seeders_sim;
        sol_jspace.Nsamp       = Nsamp;
        sol_jspace.config      = cfg_jspace;
        
        % Save
        save(fullfile('result', 'Solutions_jspace.mat'), 'sol_jspace', '-v7.3');
        fprintf('J-SPACE completed. Results saved.\n');
        
    else
        warning('solver_jspace.m not found on path. Skipping.');
    end
catch ME
    warning('J-SPACE run encountered an error: %s', ME.message);
    % Print stack trace for debugging
    for k = 1:length(ME.stack)
        fprintf('File: %s, Line: %d\n', ME.stack(k).file, ME.stack(k).line);
    end
end

end