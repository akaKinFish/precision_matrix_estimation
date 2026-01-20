function [Psi_cell, stats] = estep_posterior_psi( ...
    Svv_in, L, Sigma_jj_cell, Sigma_ee_cell, m_samples)
%ESTEP_POSTERIOR_PSI  E-step: compute posterior second moments.
%
% Math per frequency:
%   Sigma_vv = L * Sigma_jj * L' + Sigma_ee
%   A        = Sigma_jj * L' * inv(Sigma_vv)
%   Sigma_j|v= Sigma_jj - A * L * Sigma_jj
%   Psi      = Sigma_j|v + A * Svv * A'
%
% Inputs:
%   Svv_in        : (Ne x Ne x F) tensor or {F x 1} cell
%   L             : (Ne x N) leadfield
%   Sigma_jj_cell : {F x 1} or (N x N) source covariance
%   Sigma_ee_cell : {F x 1} or (Ne x Ne) sensor noise covariance
%   m_samples     : scalar or (F x 1) vector, used for log-likelihood
%
% Outputs:
%   Psi_cell      : {F x 1} posterior second moments
%   stats         : struct with log-likelihood and transfer matrices

    if nargin < 5
        m_samples = 1;
    end

    [Ne, N] = size(L);

    % Normalize Svv to cell format
    if isnumeric(Svv_in)
        if ndims(Svv_in) == 2
            Svv_in = reshape(Svv_in, size(Svv_in, 1), size(Svv_in, 2), 1);
        end
        F = size(Svv_in, 3);
        Svv_cell = cell(F, 1);
        for t = 1:F
            Svv_cell{t} = Svv_in(:, :, t);
        end
    elseif iscell(Svv_in)
        Svv_cell = Svv_in(:);
        F = numel(Svv_cell);
    else
        error('estep_posterior_psi:InvalidSvv', 'Svv_in must be numeric or a cell.');
    end

    % Normalize Sigma_jj_cell and Sigma_ee_cell to cell format
    Sigma_jj_cell = normalize_cell_matrix(Sigma_jj_cell, F, 'Sigma_jj_cell');
    Sigma_ee_cell = normalize_cell_matrix(Sigma_ee_cell, F, 'Sigma_ee_cell');

    % Expand scalar noise to scaled identity if needed
    for t = 1:F
        if isscalar(Sigma_ee_cell{t})
            Sigma_ee_cell{t} = Sigma_ee_cell{t} * eye(Ne, 'like', L);
        end
    end

    % Sample counts
    if isempty(m_samples)
        m_vec = ones(F, 1);
    elseif isscalar(m_samples)
        m_vec = repmat(m_samples, F, 1);
    else
        m_vec = m_samples(:);
        if numel(m_vec) ~= F
            error('estep_posterior_psi:SizeMismatch', ...
                'm_samples must be scalar or have length F.');
        end
    end

    Psi_cell = cell(F, 1);
    A_cell = cell(F, 1);
    loglik_total = 0;

    for t = 1:F
        Svv_data = Svv_cell{t};
        Sjj_curr = Sigma_jj_cell{t};
        See_curr = Sigma_ee_cell{t};

        if size(Svv_data, 1) ~= Ne || size(Svv_data, 2) ~= Ne
            error('estep_posterior_psi:SizeMismatch', ...
                'Svv size does not match leadfield at t=%d.', t);
        end
        if size(Sjj_curr, 1) ~= N || size(Sjj_curr, 2) ~= N
            error('estep_posterior_psi:SizeMismatch', ...
                'Sigma_jj size does not match leadfield at t=%d.', t);
        end
        if size(See_curr, 1) ~= Ne || size(See_curr, 2) ~= Ne
            error('estep_posterior_psi:SizeMismatch', ...
                'Sigma_ee size does not match leadfield at t=%d.', t);
        end

        % Enforce Hermitian on input covariance
        Svv_data = 0.5 * (Svv_data + Svv_data');
        Svv_data(1:Ne+1:end) = real(diag(Svv_data));

        % Step A: Model covariance
        Sigma_vv = L * Sjj_curr * L' + See_curr;
        Sigma_vv = 0.5 * (Sigma_vv + Sigma_vv');
        Sigma_vv(1:Ne+1:end) = real(diag(Sigma_vv));

        % Step B: Robust inverse and logdet
        [inv_Svv, logdet_Svv, is_spd] = spd_inv_and_logdet(Sigma_vv);
        if ~is_spd
            warning('estep_posterior_psi:NonSPD', ...
                'Frequency %d: Sigma_vv is not SPD after jitter.', t);
        end

        % Step C: Inverse operator
        A = Sjj_curr * L' * inv_Svv;

        % Step D: Posterior covariance
        Sigma_post = Sjj_curr - A * (L * Sjj_curr);
        Sigma_post = 0.5 * (Sigma_post + Sigma_post');

        % Step E: Posterior second moment
        Sjj_empirical = A * Svv_data * A';
        Psi = Sigma_post + Sjj_empirical;
        Psi = 0.5 * (Psi + Psi');
        Psi(1:N+1:end) = real(diag(Psi));

        Psi_cell{t} = Psi;
        A_cell{t} = A;

        if is_spd
            tr_val = real(trace(inv_Svv * Svv_data));
            loglik_total = loglik_total - 0.5 * m_vec(t) * (logdet_Svv + tr_val);
        end
    end

    stats = struct();
    stats.loglik = loglik_total;
    stats.A_cell = A_cell;
    stats.description = 'Psi = E[jj'' | Svv]';
end

function C = normalize_cell_matrix(x, F, name)
    if nargin < 3
        name = 'input';
    end
    if isempty(x)
        error('estep_posterior_psi:MissingInput', '%s is empty.', name);
    end
    if iscell(x)
        if numel(x) == 1 && F > 1
            C = repmat(x, F, 1);
        elseif numel(x) == F
            C = x(:);
        else
            error('estep_posterior_psi:SizeMismatch', ...
                '%s must have length 1 or F.', name);
        end
    else
        C = repmat({x}, F, 1);
    end
end
