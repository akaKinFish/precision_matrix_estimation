function [Sigma_tilde, D, stats] = module1_data_whitening(Sigma_emp, varargin)
% MODULE1_DATA_WHITENING - Streamlined covariance whitening module.
%
% Purpose:
%   Transforms empirical covariance matrices to "whitened" space where
%   diagonal elements are normalized to approximately 1.
%   Optionally applies frequency-domain smoothing to the diagonal powers.
%
% Mathematical Operation:
%   1. Extract diagonal powers: g(w) = diag(Sigma(w))
%   2. Smooth powers (optional): g_smooth(w) = MovingAverage(g(w))
%   3. Construct whitening matrix: D(w) = diag( 1 ./ sqrt(g_smooth(w)) )
%   4. Whiten: Sigma_tilde(w) = D(w) * Sigma(w) * D(w)
%
% Usage:
%   [S_tilde, D] = module1_data_whitening(Sigma_cell);
%   [S_tilde, D] = module1_data_whitening(Sigma_cell, 'smoothing_window', 5);
%
% Inputs:
%   Sigma_emp : {F x 1} Cell array of covariance matrices (p x p)
%               OR (p x p x F) 3D array
%               OR (p x p) single matrix
%
% Parameters (Name-Value):
%   'smoothing_window' : (int) Window size for frequency smoothing. 
%                        Set to 1 or 0 to disable smoothing. Default: 1 (No smoothing)
%   'min_power'        : (double) Floor for diagonal values to avoid div-by-zero. Default: 1e-12.
%   'force_hermitian'  : (bool) Force output to be Hermitian. Default: true.
%
% Outputs:
%   Sigma_tilde : {F x 1} Whitened covariance matrices.
%   D           : {F x 1} Diagonal whitening matrices.
%   stats       : Struct with processing statistics.

    % ============================================================
    % 1. Input Parsing & Validation
    % ============================================================
    p_parser = inputParser;
    addRequired(p_parser, 'Sigma_emp');
    addParameter(p_parser, 'smoothing_window', 1, @(x) isscalar(x) && x >= 0);
    addParameter(p_parser, 'min_power', 1e-12, @(x) isscalar(x) && x > 0);
    addParameter(p_parser, 'force_hermitian', true, @islogical);
    parse(p_parser, Sigma_emp, varargin{:});
    
    opts = p_parser.Results;
    
    % Normalize Input to Cell Array {F x 1}
    if iscell(Sigma_emp)
        Sigma_in = Sigma_emp;
    elseif isnumeric(Sigma_emp) && ndims(Sigma_emp) == 3
        [p1, p2, F_in] = size(Sigma_emp);
        Sigma_in = cell(F_in, 1);
        for i = 1:F_in, Sigma_in{i} = Sigma_emp(:,:,i); end
    elseif isnumeric(Sigma_emp) && ismatrix(Sigma_emp)
        Sigma_in = {Sigma_emp};
    else
        error('Module1:InvalidInput', 'Input must be Cell, 3D array, or Matrix.');
    end
    
    F = numel(Sigma_in);
    p = size(Sigma_in{1}, 1);
    
    % Check dimensions
    if size(Sigma_in{1}, 2) ~= p
        error('Module1:NonSquare', 'Covariance matrices must be square.');
    end

    % ============================================================
    % 2. Extract & Smooth Diagonals
    % ============================================================
    % Extract raw diagonals: Matrix [p x F]
    raw_diagonals = zeros(p, F);
    for f = 1:F
        % Ensure real diagonal (power spectrum is real)
        raw_diagonals(:, f) = real(diag(Sigma_in{f}));
    end
    
    % Apply Smoothing (if F > 1 and window > 1)
    smooth_diagonals = raw_diagonals;
    if F > 1 && opts.smoothing_window > 1
        % Apply simple moving average across frequency dimension
        w_size = opts.smoothing_window;
        for i = 1:p
            % Using built-in smooth or custom filter
            % Custom implementation to avoid Toolbox dependency:
            smooth_diagonals(i, :) = simple_moving_average(raw_diagonals(i, :), w_size);
        end
    end
    
    % Apply Min Power Floor (Regularization)
    smooth_diagonals = max(smooth_diagonals, opts.min_power);

    % ============================================================
    % 3. Construct D and Whiten
    % ============================================================
    Sigma_tilde = cell(F, 1);
    D = cell(F, 1);
    
    stats.avg_diagonal_error = 0;
    stats.hermitian_enforced = 0;
    
    for f = 1:F
        % 3a. Construct Whitening Matrix D
        % D = diag(1 ./ sqrt(g))
        d_vec = 1 ./ sqrt(smooth_diagonals(:, f));
        D_mat = diag(d_vec);
        D{f} = D_mat;
        
        % 3b. Apply Whitening: D * Sigma * D'
        % Optimization: Element-wise multiplication is faster for diagonal D
        % Sigma_new(i,j) = Sigma(i,j) * d(i) * d(j)
        S_raw = Sigma_in{f};
        
        % Use broadcasting (outer product of d_vec)
        Scaling = d_vec * d_vec'; 
        S_white = S_raw .* Scaling;
        
        % 3c. Enforce Hermitian (if requested)
        if opts.force_hermitian
            if norm(S_white - S_white', 'fro') > 1e-12
                stats.hermitian_enforced = stats.hermitian_enforced + 1;
            end
            S_white = (S_white + S_white') / 2;
        end
        
        Sigma_tilde{f} = S_white;
        
        % Stats: How close is diagonal to 1?
        current_diag = real(diag(S_white));
        stats.avg_diagonal_error = stats.avg_diagonal_error + mean(abs(current_diag - 1.0));
    end
    
    stats.avg_diagonal_error = stats.avg_diagonal_error / F;
    stats.n_nodes = p;
    stats.n_freq = F;
    
end

% ============================================================
% Helper: Simple Moving Average (Dependency Free)
% ============================================================
function y = simple_moving_average(x, window)
    % Handles vector x with window size. 
    % Uses 'valid' convolution logic with padding for boundary handling.
    
    N = length(x);
    if window >= N
        y = ones(1, N) * mean(x);
        return;
    end
    
    % Kernel
    k = ones(1, window) / window;
    
    % Convolve (same size output)
    % Reflective padding is usually better for spectra than zero padding
    pad_len = floor(window/2);
    x_pad = [x(pad_len:-1:1), x, x(end:-1:end-pad_len+1)];
    
    y_conv = conv(x_pad, k, 'valid');
    
    % Trim if size mismatch due to odd/even window
    if length(y_conv) > N
        y_conv = y_conv(1:N);
    elseif length(y_conv) < N
        % Should not happen with proper padding logic, but safety:
        y_conv = [y_conv, x(end-(N-length(y_conv))+1:end)];
    end
    y = y_conv;
end