function Theta_cell = enforce_diag_policy(Theta_cell, policy_mode)
%ENFORCE_DIAG_POLICY  Constrain the diagonal structure of spline coefficients.
%
% Modes:
%   'theta0_only' (Recommended):
%       Theta_0 keeps its diagonal.
%       Theta_1...Theta_K diagonals are forced to 0.
%       This separates base precision from frequency fluctuations.
%
%   'free':
%       All Theta_k can have diagonals.

    if nargin < 2 || isempty(policy_mode)
        policy_mode = 'theta0_only';
    end

    K = numel(Theta_cell);
    N = size(Theta_cell{1}, 1);

    mode = lower(policy_mode);
    if strcmp(mode, 'theta0_only')
        Theta_cell{1}(1:N+1:end) = real(diag(Theta_cell{1}));
        for k = 2:K
            Theta_cell{k}(1:N+1:end) = 0;
        end
    elseif strcmp(mode, 'free')
        for k = 1:K
            Theta_cell{k}(1:N+1:end) = real(diag(Theta_cell{k}));
        end
    else
        % Fallback to safe default
        Theta_cell{1}(1:N+1:end) = real(diag(Theta_cell{1}));
        for k = 2:K
            Theta_cell{k}(1:N+1:end) = 0;
        end
    end
end
