function [Sigma_tilde, quality] = apply_whitening_with_fix(Sigma_emp, D, varargin)
% APPLY_WHITENING_WITH_FIX - Fixed version of whitening transformation
%
% This function applies whitening with proper array dimension handling
%
% Usage:
%   [Sigma_tilde, quality] = apply_whitening_with_fix(Sigma_emp, D)

    % Parse parameters
    p = inputParser;
    addParameter(p, 'verbose', true, @islogical);
    addParameter(p, 'target_diagonal', 1.0, @isnumeric);
    addParameter(p, 'diagonal_tolerance', 0.1, @isnumeric);
    parse(p, varargin{:});
    params = p.Results;
    
    F = length(Sigma_emp);
    n = size(Sigma_emp{1}, 1);
    
    Sigma_tilde = cell(F, 1);
    quality = struct();
    quality.diagonal_errors = zeros(F, 1);
    quality.whitening_effectiveness = zeros(F, 1);
    
    for omega = 1:F
        % Get whitening matrix diagonal
        d_omega = diag(D{omega});
        
        % Ensure positive values
        d_omega = max(d_omega, 1e-10);
        
        % Create whitening matrix
        D_sqrt = diag(sqrt(d_omega));
        
        % Apply whitening
        Sigma_tilde{omega} = D_sqrt * Sigma_emp{omega} * D_sqrt;
        
        % Force Hermitian
        Sigma_tilde{omega} = (Sigma_tilde{omega} + Sigma_tilde{omega}') / 2;
        
        % Calculate quality metrics
        diag_elements = real(diag(Sigma_tilde{omega}));
        quality.diagonal_errors(omega) = max(abs(diag_elements - params.target_diagonal));
        quality.whitening_effectiveness(omega) = 1 / (1 + quality.diagonal_errors(omega));
    end
    
    % Overall quality score
    quality.overall_score = mean(quality.whitening_effectiveness) * 100;
    
    if params.verbose
        fprintf('Whitening completed: Overall quality = %.1f%%\n', quality.overall_score);
    end
end
