function [Psijj_cell, stats] = module2_estep(Svv_cell, L, Sigmajj_cell, Sigma_noise, varargin)
% MODULE2_ESTEP - Expectation Step (Source Posterior Moments)
%
% Purpose:
%   Calculates the expected second moment of the sources (Psijj), which serves
%   as the input for the M-step (Module 1 & PGD).
%
% Mathematical Equivalence:
%   This function implements the exact same logic as 'higgs_expectation',
%   but uses the numerically stable Kalman Gain formulation:
%
%   1. Model Covariance: Q = L * Sigmajj * L' + Sigma_noise
%   2. Inverse Operator: T = Sigmajj * L' * inv(Q)
%   3. Posterior Cov :   Sigma_post = Sigmajj - T * L * Sigmajj
%   4. Empirical Cov :   Sjj = T * Svv * T'
%   5. Second Moment :   Psijj = Sigma_post + Sjj
%
% Inputs:
%   Svv_cell      : {F x 1} Sensor covariance matrices (v * v')
%   L             : (Ns x Nr) Leadfield matrix
%   Sigmajj_cell  : {F x 1} Current source covariance estimates (D^(k))
%   Sigma_noise   : (Ns x Ns) Noise covariance matrix (R)
%
% Outputs:
%   Psijj_cell    : {F x 1} The "Psi" matrix for M-step (E[ss'|v])
%   stats         : Struct with T (Transfer function) and Likelihoods

    % 1. Input Parsing
    if ~iscell(Svv_cell), Svv_cell = {Svv_cell}; end
    if ~iscell(Sigmajj_cell), Sigmajj_cell = {Sigmajj_cell}; end
    
    F = numel(Svv_cell);
    [Ns, Nr] = size(L);
    
    Psijj_cell = cell(F, 1);
    
    % Storage for transfer functions (optional, for activation analysis)
    T_cell = cell(F, 1);
    log_lik_sum = 0;
    
    % 2. Loop over frequencies
    for f = 1:F
        Svv = Svv_cell{f};          % Sensor Data Covariance
        Sigmajj = Sigmajj_cell{f};  % Current Source Covariance
        
        % --- Step A: Model Covariance (Q) ---
        % Q = L * D * L' + R
        % This represents the covariance we EXPECT to see at sensors
        Q = (L * Sigmajj * L') + Sigma_noise;
        Q = utils_math.make_hermitian(Q); % Numerical stability
        
        % --- Step B: Inverse Operator / DSTF (T) ---
        % higgs code: Tjv = Sigmajj_post * L' * inv(R)
        % Our code:   T   = Sigmajj * L' * inv(Q)
        % These are mathematically IDENTICAL via Matrix Inversion Lemma,
        % but inv(Q) is safer when R is small.
        
        % Use Cholesky solve for speed and stability: T = (Q \ (L * Sigmajj'))'
        % Note: Sigmajj is Hermitian, so Sigmajj' = Sigmajj
        K_gain = Q \ (L * Sigmajj); 
        T = K_gain'; 
        
        % --- Step C: Posterior Covariance (Sigma_post) ---
        % Sigma_post = Sigmajj - T * L * Sigmajj
        % "How much uncertainty is reduced by the data?"
        Sigma_post = Sigmajj - T * (L * Sigmajj);
        Sigma_post = utils_math.make_hermitian(Sigma_post);
        
        % --- Step D: Empirical Source Covariance (Sjj) ---
        % Sjj = T * Svv * T'
        % "What the data actually says the sources look like"
        Sjj = T * Svv * T';
        Sjj = utils_math.make_hermitian(Sjj);
        
        % --- Step E: Effective Second Moment (Psijj) ---
        % Psijj = Sigma_post + Sjj
        % This is the target for the M-step optimization
        Psijj = Sigma_post + Sjj;
        
        % Store results
        Psijj_cell{f} = Psijj;
        T_cell{f} = T;
        
        % (Optional) Log-Likelihood calculation for monitoring
        % LL = -log|Q| - tr(Q^-1 * Svv)
        [L_chol, is_spd] = utils_math.safe_log_det(Q);
        if is_spd
            term1 = L_chol; % log|Q|
            % Use only the stable solve against Q to avoid rank warnings from dividing by L'
            term2 = real(trace(Q \ Svv));
            log_lik_sum = log_lik_sum - 0.5 * (term1 + term2);
        end
    end
    
    % 3. Pack stats
    stats.T_transfer = T_cell;
    stats.log_likelihood = log_lik_sum;
    stats.description = 'Psijj is the E-step Second Moment (Input to M-step)';

end
