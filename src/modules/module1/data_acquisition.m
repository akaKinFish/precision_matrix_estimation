function Sigma_emp = data_acquisition(varargin)
% DATA_ACQUISITION - Load or generate empirical covariance matrices
%
% This function handles both real EEG data loading and simulation data
% integration. It serves as the entry point for Module 1.
%
% Usage:
%   Sigma_emp = data_acquisition('simulation', sim_results)
%   Sigma_emp = data_acquisition('eeg_data', data_path, params)
%
% Inputs:
%   varargin - Variable arguments:
%     For simulation mode:
%       'simulation', sim_results - Results from Module 7
%     For real data mode:
%       'eeg_data', data_path, params - EEG data path and parameters
%
% Outputs:
%   Sigma_emp - Cell array of empirical covariance matrices {F x 1}
%              Each element is complex matrix of size [n x n]
%
% File location: src/modules/module1_preprocessing/data_acquisition.m

    % Parse input arguments
    if nargin < 2
        error('data_acquisition:insufficient_args', ...
              'At least 2 arguments required');
    end
    
    mode = varargin{1};
    
    switch lower(mode)
        case 'simulation'
            % Load from simulation results (Module 7 output)
            if nargin < 2
                error('data_acquisition:missing_sim_results', ...
                      'Simulation results required for simulation mode');
            end
            sim_results = varargin{2};
            Sigma_emp = load_simulation_data(sim_results);
            
        case 'eeg_data'
            % Load from real EEG data files
            if nargin < 3
                error('data_acquisition:missing_eeg_params', ...
                      'Data path and parameters required for EEG mode');
            end
            data_path = varargin{2};
            params = varargin{3};
            Sigma_emp = load_eeg_data(data_path, params);
            
        otherwise
            error('data_acquisition:invalid_mode', ...
                  'Mode must be ''simulation'' or ''eeg_data''');
    end
    
    % Validate output
    validate_covariance_matrices(Sigma_emp);
    
    fprintf('Data acquisition completed. Loaded %d frequency points.\n', ...
            length(Sigma_emp));
end

function Sigma_emp = load_simulation_data(sim_results)
% Load empirical covariances from simulation results
    
    % Validate simulation results structure
    required_fields = {'Sigma_emp', 'F', 'n', 'T'};
    for i = 1:length(required_fields)
        if ~isfield(sim_results, required_fields{i})
            error('data_acquisition:invalid_sim_results', ...
                  'Missing field: %s', required_fields{i});
        end
    end
    
    % Extract empirical covariances
    Sigma_emp = sim_results.Sigma_emp;
    
    % Validate dimensions
    F = sim_results.F;
    n = sim_results.n;
    
    if length(Sigma_emp) ~= F
        error('data_acquisition:dimension_mismatch', ...
              'Number of covariance matrices (%d) does not match F (%d)', ...
              length(Sigma_emp), F);
    end
    
    for omega = 1:F
        [rows, cols] = size(Sigma_emp{omega});
        if rows ~= n || cols ~= n
            error('data_acquisition:size_mismatch', ...
                  'Covariance matrix %d has size [%d x %d], expected [%d x %d]', ...
                  omega, rows, cols, n, n);
        end
    end
    
    fprintf('Loaded simulation data: %d nodes, %d frequencies, %d samples per frequency\n', ...
            n, F, sim_results.T);
end

function Sigma_emp = load_eeg_data(data_path, params)
% Load empirical covariances from real EEG data
    
    % Validate parameters
    default_params = struct('n', [], 'F', [], 'T', [], ...
                           'fs', 500, 'freq_range', [1, 50], ...
                           'overlap', 0.5, 'window_type', 'hann');
    params = merge_structs(default_params, params);
    
    if isempty(params.n) || isempty(params.F) || isempty(params.T)
        error('data_acquisition:missing_params', ...
              'Parameters n, F, and T must be specified for EEG data');
    end
    
    % Check if data file exists
    if ~exist(data_path, 'file')
        error('data_acquisition:file_not_found', ...
              'EEG data file not found: %s', data_path);
    end
    
    try
        % Load EEG data (assuming .mat file with variable 'eeg_data')
        fprintf('Loading EEG data from: %s\n', data_path);
        data_struct = load(data_path);
        
        % Try to find EEG data in common variable names
        possible_names = {'eeg_data', 'data', 'X', 'signal'};
        eeg_data = [];
        for i = 1:length(possible_names)
            if isfield(data_struct, possible_names{i})
                eeg_data = data_struct.(possible_names{i});
                fprintf('Found EEG data in variable: %s\n', possible_names{i});
                break;
            end
        end
        
        if isempty(eeg_data)
            error('data_acquisition:no_eeg_data', ...
                  'Could not find EEG data in file. Expected variables: %s', ...
                  strjoin(possible_names, ', '));
        end
        
        % Convert time-domain data to frequency-domain covariances
        Sigma_emp = compute_frequency_covariances(eeg_data, params);
        
    catch ME
        error('data_acquisition:load_failed', ...
              'Failed to load EEG data: %s', ME.message);
    end
end

function Sigma_emp = compute_frequency_covariances(eeg_data, params)
% Compute frequency-domain covariance matrices from time-domain EEG data
    
    [n_channels, n_samples] = size(eeg_data);
    
    if n_channels ~= params.n
        warning('data_acquisition:channel_mismatch', ...
                'EEG data has %d channels, expected %d. Adjusting...', ...
                n_channels, params.n);
        if n_channels > params.n
            eeg_data = eeg_data(1:params.n, :);
        else
            error('data_acquisition:insufficient_channels', ...
                  'EEG data has only %d channels, need %d', ...
                  n_channels, params.n);
        end
    end
    
    % Compute frequency axis
    freq_min = params.freq_range(1);
    freq_max = params.freq_range(2);
    frequencies = linspace(freq_min, freq_max, params.F);
    
    % Initialize output
    Sigma_emp = cell(params.F, 1);
    
    % Window parameters for short-time FFT
    window_length = round(params.fs / freq_min * 4); % 4 cycles at min frequency
    window_length = min(window_length, floor(n_samples / 4));
    overlap_samples = round(window_length * params.overlap);
    
    fprintf('Computing frequency-domain covariances...\n');
    fprintf('Window length: %d samples, Overlap: %d samples\n', ...
            window_length, overlap_samples);
    
    for omega = 1:params.F
        freq = frequencies(omega);
        
        % Extract complex-valued data at this frequency using short-time FFT
        X_omega = compute_stft_at_frequency(eeg_data, freq, params.fs, ...
                                           window_length, overlap_samples, ...
                                           params.window_type);
        
        % Compute empirical covariance
        T_freq = size(X_omega, 2);
        if T_freq < params.T
            warning('data_acquisition:insufficient_samples', ...
                    'Frequency %d has only %d samples, expected %d', ...
                    omega, T_freq, params.T);
        end
        
        % Use first T samples or all available samples
        T_use = min(T_freq, params.T);
        X_omega = X_omega(:, 1:T_use);
        
        % Compute empirical covariance: Sigma = (1/T) * X * X^H
        Sigma_emp{omega} = (X_omega * X_omega') / T_use;
        
        if mod(omega, max(1, floor(params.F/10))) == 0
            fprintf('Processed frequency %d/%d (%.1f Hz)\n', ...
                    omega, params.F, freq);
        end
    end
    
    fprintf('Frequency-domain covariance computation completed.\n');
end

function X_freq = compute_stft_at_frequency(data, target_freq, fs, ...
                                           window_length, overlap, window_type)
% Compute short-time Fourier transform at specific frequency
    
    [n_channels, n_samples] = size(data);
    
    % Create window
    switch lower(window_type)
        case 'hann'
            window = hann(window_length)';
        case 'hamming'
            window = hamming(window_length)';
        case 'blackman'
            window = blackman(window_length)';
        otherwise
            window = ones(1, window_length); % rectangular
    end
    
    % Calculate number of segments
    step = window_length - overlap;
    n_segments = floor((n_samples - overlap) / step);
    
    % Initialize output
    X_freq = complex(zeros(n_channels, n_segments));
    
    % Frequency index for target frequency
    freq_axis = (0:window_length-1) * fs / window_length;
    [~, freq_idx] = min(abs(freq_axis - target_freq));
    
    % Compute STFT at target frequency for each segment
    for seg = 1:n_segments
        start_idx = (seg - 1) * step + 1;
        end_idx = start_idx + window_length - 1;
        
        if end_idx > n_samples
            break;
        end
        
        % Extract windowed data
        windowed_data = data(:, start_idx:end_idx) .* window;
        
        % Compute FFT and extract target frequency
        fft_data = fft(windowed_data, [], 2);
        X_freq(:, seg) = fft_data(:, freq_idx);
    end
    
    % Trim to actual number of computed segments
    X_freq = X_freq(:, 1:min(seg, n_segments));
end

function validate_covariance_matrices(Sigma_emp)
% Validate empirical covariance matrices
    
    if ~iscell(Sigma_emp)
        error('data_acquisition:invalid_format', ...
              'Sigma_emp must be a cell array');
    end
    
    F = length(Sigma_emp);
    if F == 0
        error('data_acquisition:empty_data', ...
              'No covariance matrices found');
    end
    
    % Check first matrix to determine expected size
    if isempty(Sigma_emp{1})
        error('data_acquisition:empty_matrix', ...
              'First covariance matrix is empty');
    end
    
    [n, n_check] = size(Sigma_emp{1});
    if n ~= n_check
        error('data_acquisition:not_square', ...
              'Covariance matrices must be square');
    end
    
    % Validate each matrix
    for omega = 1:F
        Sigma_omega = Sigma_emp{omega};
        
        % Check dimensions
        [rows, cols] = size(Sigma_omega);
        if rows ~= n || cols ~= n
            error('data_acquisition:inconsistent_size', ...
                  'Matrix %d has size [%d x %d], expected [%d x %d]', ...
                  omega, rows, cols, n, n);
        end
        
        % Check for complex values
        if ~isnumeric(Sigma_omega)
            error('data_acquisition:non_numeric', ...
                  'Matrix %d contains non-numeric values', omega);
        end
        
        % Check for Hermitian property (within tolerance)
        hermitian_error = max(max(abs(Sigma_omega - Sigma_omega')));
        if hermitian_error > 1e-12
            warning('data_acquisition:not_hermitian', ...
                    'Matrix %d is not Hermitian (error: %.2e)', ...
                    omega, hermitian_error);
        end
        
        % Check positive semi-definiteness
        eigenvals = eig(Sigma_omega);
        min_eigenval = min(real(eigenvals));
        if min_eigenval < -1e-12
            warning('data_acquisition:not_psd', ...
                    'Matrix %d may not be positive semi-definite (min eigenvalue: %.2e)', ...
                    omega, min_eigenval);
        end
    end
    
    fprintf('Validation completed: %d matrices of size [%d x %d]\n', F, n, n);
end

function result = merge_structs(default_struct, user_struct)
% Merge user parameters with defaults
    result = default_struct;
    if ~isempty(user_struct)
        user_fields = fieldnames(user_struct);
        for i = 1:length(user_fields)
            result.(user_fields{i}) = user_struct.(user_fields{i});
        end
    end
end