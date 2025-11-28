function [Svv_cell, Sjj_true_cell, L, freqs] = simulation_nmm_wrapper(parameters)
% SIMULATION_NMM_WRAPPER - Force Sparse Generation
%
% Key Fix: Drastically reduce structural connectivity density and strength
% to prevent whole-brain synchronization (yellow matrix issue).

    % 1. Extract Parameters
    Ne = parameters.Dimensions.Ne; 
    Nr = parameters.Dimensions.Nr; 
    Nw = parameters.Dimensions.Nw; 
    
    % === FIX 1: Enforce Sparse Connectivity ===
    C_raw = parameters.Compact_Model.C; 
    D = parameters.Compact_Model.D; 
    L = parameters.Compact_Model.K; 
    
    % Keep only top 2% edges
    if max(C_raw(:)) > 0
        thresh = quantile(C_raw(C_raw>0), 0.98);
        C = zeros(Nr, Nr);
        % Set coupling strength to 0.5 (weak enough to avoid lock-step)
        C(C_raw > thresh) = 0.5; 
    else
        C = zeros(Nr, Nr); % Should not happen
    end
    
    % Simulation Constants
    T = 1; dt = 0.001; 
    n_steps = round(T / dt);
    D_steps = round(D / dt);
    
    % === FIX 2: Increase Independent Noise ===
    noise_std = 2.0; % High independent noise
    measurement_noise_std = 0.1;
    M = 10; 

    % Oscillator Params
    tau = 0.02; 
    omega1 = 0; 
    rng('shuffle');
    freq0 = 10; % Fixed 10Hz for clear Alpha test
    omega2 = 2 * pi * freq0; 
    zeta1 = 0.00001;
    zeta2 = 0.00000001;

    % === FIX 3: Sparse Activation ===
    mask1 = zeros(Nr, 1); 
    mask2 = zeros(Nr, 1);
    
    % Activate 5 random nodes strongly
    seeds = randperm(Nr, 5);
    mask1(seeds) = 10.0; % Strong kick
    
    fprintf('[NMM Sim] Running Sparse Simulation (M=%d)...\n', M);
    
    % 2. Run Simulation
    u_total_ensemble = zeros(M, Nr, n_steps);
    if isempty(gcp('nocreate')), pool_size = 0; else, pool_size = Inf; end
    
    parfor (m = 1:M, pool_size)
        u1 = zeros(Nr, n_steps); v1 = zeros(Nr, n_steps);
        u2 = zeros(Nr, n_steps); v2 = zeros(Nr, n_steps);
        
        u1(:,1) = mask1; u2(:,1) = mask2;
        I1 = noise_std * randn(Nr, n_steps);
        I2 = noise_std * randn(Nr, n_steps);
        
        for t = 2:n_steps
            input_sum = zeros(Nr, 1);
            % Naive loop is slow but correct. 
            % Optimization: Pre-calculate indices of C > 0
            % For now, trust parfor to handle it.
            for i = 1:Nr
                for j = 1:Nr
                    if C(i,j) > 0
                        d_idx = t - D_steps(i,j);
                        if d_idx > 0
                            delayed = u1(j, d_idx) + u2(j, d_idx);
                            % Sigmoid coupling
                            input_sum(i) = input_sum(i) + C(i,j) * (1 / (1 + exp(-delayed)));
                        end
                    end
                end
            end
            
            % Oscillator 1
            du1 = v1(:, t-1);
            dv1 = -2*zeta1*omega1*v1(:,t-1) - omega1^2*u1(:,t-1) + input_sum + I1(:,t);
            u1(:,t) = u1(:,t-1) + dt*du1; 
            v1(:,t) = v1(:,t-1) + dt*dv1;
            
            % Oscillator 2
            du2 = v2(:, t-1);
            dv2 = -2*zeta2*omega2*v2(:,t-1) - omega2^2*u2(:,t-1) + input_sum + I2(:,t);
            u2(:,t) = u2(:,t-1) + dt*du2;
            v2(:,t) = v2(:,t-1) + dt*dv2;
        end
        u_total_ensemble(m, :, :) = u1 + u2;
    end
    
    % 3. FFT & Cross-Spectra
    freqs_full = (0:(n_steps/2)) / (n_steps*dt);
    target_freqs = linspace(1, 40, Nw); 
    fft_indices = zeros(1, Nw);
    for k=1:Nw, [~, fft_indices(k)] = min(abs(freqs_full - target_freqs(k))); end
    freqs = freqs_full(fft_indices);
    
    Sjj_avg = zeros(Nr, Nr, Nw);
    Svv_avg = zeros(Ne, Ne, Nw);
    
    for m = 1:M
        u_total = squeeze(u_total_ensemble(m, :, :));
        
        U_fft = fft(u_total, [], 2);
        U_sel = U_fft(:, fft_indices);
        
        % Normalize source power to avoid numerical explosion
        U_sel = U_sel ./ (mean(abs(U_sel(:))) + eps);
        
        v_total = L * u_total + measurement_noise_std * randn(Ne, n_steps);
        V_fft = fft(v_total, [], 2);
        V_sel = V_fft(:, fft_indices);
        V_sel = V_sel ./ (mean(abs(V_sel(:))) + eps);
        
        for k = 1:Nw
            uk = U_sel(:, k);
            vk = V_sel(:, k);
            Sjj_avg(:,:,k) = Sjj_avg(:,:,k) + (uk * uk') / n_steps;
            Svv_avg(:,:,k) = Svv_avg(:,:,k) + (vk * vk') / n_steps;
        end
    end
    
    Sjj_avg = Sjj_avg / M;
    Svv_avg = Svv_avg / M;
    
    Svv_cell = cell(Nw, 1);
    Sjj_true_cell = cell(Nw, 1);
    for k = 1:Nw
        Svv_cell{k} = Svv_avg(:,:,k);
        Sjj_true_cell{k} = Sjj_avg(:,:,k);
    end
    
    fprintf('[NMM Sim] Sparse Data generated.\n');
end