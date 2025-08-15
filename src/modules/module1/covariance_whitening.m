function [Sigma_tilde, whitening_quality] = covariance_whitening(Sigma_emp, D, varargin)
% COVARIANCE_WHITENING - Legacy wrapper for CovarianceWhitening.whiten
%
% This function provides backward compatibility for existing code that calls
% the original covariance_whitening function. All functionality is preserved.
% 
% IMPORTANT: This is exactly the same as your original function, just moved to a class!
%
% Usage (SAME as before):
%   [Sigma_tilde, quality] = covariance_whitening(Sigma_emp, D)
%   [Sigma_tilde, quality] = covariance_whitening(Sigma_emp, D, 'param', value, ...)
%
% This is equivalent to:
%   [Sigma_tilde, quality] = CovarianceWhitening.whiten(Sigma_emp, D, ...)
%
% NEW capabilities available in the class (not available in this wrapper):
%   - CovarianceWhitening.compute_quality_metrics() - Independent quality assessment
%   - CovarianceWhitening.report_quality() - Standalone quality reporting
%
% File location: src/modules/module1/covariance_whitening.m
% Status: COMPATIBLE - Works exactly like your original function

    % Call the new class method with all the same parameters
    [Sigma_tilde, whitening_quality] = CovarianceWhitening.whiten(Sigma_emp, D, varargin{:});
    
end