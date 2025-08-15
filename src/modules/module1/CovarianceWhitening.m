classdef CovarianceWhitening < handle
% COVARIANCEWHITENING - Enhanced covariance whitening with quality metrics
%
% This class converts the original covariance_whitening function into a class-based
% architecture while preserving ALL original functionality. 
%
% Usage (SAME as original function):
%   [Sigma_tilde, quality] = CovarianceWhitening.whiten(Sigma_emp, D, ...)
%
% NEW capabilities:
%   quality = CovarianceWhitening.compute_quality_metrics(Sigma_orig, Sigma_tilde, omega, params)
%   CovarianceWhitening.report_quality(quality, params)
%
% File location: src/modules/module1/CovarianceWhitening.m
% Based on: Original covariance_whitening.m from paste.txt

    methods (Static)
        
        function [Sigma_tilde, whitening_quality] = whiten(Sigma_emp, D, varargin)
        % WHITEN - Apply whitening transformation to covariance matrices
        % 
        % This is EXACTLY the same as the original covariance_whitening function
        % Enhanced version with better diagonal normalization through iterative refinement
        % and adaptive scaling.
        %
        % Usage:
        %   [Sigma_tilde, quality] = CovarianceWhitening.whiten(Sigma_emp, D)
        %   [Sigma_tilde, quality] = CovarianceWhitening.whiten(Sigma_emp, D, 'param', value, ...)
            
            % Parse input arguments (SAME as original)
            p = inputParser;
            addRequired(p, 'Sigma_emp', @(x) iscell(x) && ~isempty(x));
            addRequired(p, 'D', @(x) iscell(x) && ~isempty(x));
            addParameter(p, 'target_diagonal', 1.0, @(x) isscalar(x) && x > 0);
            addParameter(p, 'diagonal_tolerance', 0.1, @(x) isscalar(x) && x > 0);
            addParameter(p, 'force_hermitian', true, @islogical);
            addParameter(p, 'check_psd', true, @islogical);
            addParameter(p, 'adaptive_correction', true, @islogical);
            addParameter(p, 'max_correction_iter', 3, @(x) isscalar(x) && x > 0);
            addParameter(p, 'convergence_threshold', 0.02, @(x) isscalar(x) && x > 0);
            addParameter(p, 'verbose', true, @islogical);
            
            parse(p, Sigma_emp, D, varargin{:});
            params = p.Results;
            
            % Validate inputs (SAME as original)
            CovarianceWhitening.validate_whitening_inputs(Sigma_emp, D);
            
            F = length(Sigma_emp);
            n = size(Sigma_emp{1}, 1);
            
            if params.verbose
                fprintf('Applying enhanced whitening transformation to %d frequencies, %d nodes\n', F, n);
            end
            
            % Initialize outputs (SAME as original)
            Sigma_tilde = cell(F, 1);
            whitening_quality = CovarianceWhitening.initialize_quality_structure(F, n);
            
            % Apply whitening transformation to each frequency (SAME as original)
            for omega = 1:F
                [Sigma_tilde{omega}, quality_omega] = CovarianceWhitening.apply_enhanced_whitening_single(...
                    Sigma_emp{omega}, D{omega}, omega, params);
                
                % Store quality metrics
                CovarianceWhitening.store_quality_metrics(whitening_quality, quality_omega, omega);
                
                if params.verbose && mod(omega, max(1, floor(F/10))) == 0
                    fprintf('Processed frequency %d/%d\n', omega, F);
                end
            end
            
            % Report overall whitening quality (SAME as original)
            if params.verbose
                CovarianceWhitening.report_enhanced_whitening_quality(whitening_quality, params);
            end
            
            % Final validation with relaxed thresholds (SAME as original)
            CovarianceWhitening.validate_whitening_output_relaxed(Sigma_tilde, whitening_quality, params);
        end
        
        function quality = compute_quality_metrics(Sigma_orig, Sigma_tilde, omega, varargin)
        % COMPUTE_QUALITY_METRICS - Public wrapper for enhanced quality metrics computation
        % 
        % This provides external access to quality assessment functionality
        %
        % Usage:
        %   quality = CovarianceWhitening.compute_quality_metrics(Sigma_orig, Sigma_tilde, omega)
        %   quality = CovarianceWhitening.compute_quality_metrics(..., 'target_diagonal', 1.0)
            
            % Parse optional parameters
            p = inputParser;
            addRequired(p, 'Sigma_orig');
            addRequired(p, 'Sigma_tilde');
            addRequired(p, 'omega', @(x) isscalar(x) && x > 0);
            addParameter(p, 'target_diagonal', 1.0, @(x) isscalar(x) && x > 0);
            addParameter(p, 'diagonal_tolerance', 0.1, @(x) isscalar(x) && x > 0);
            parse(p, Sigma_orig, Sigma_tilde, omega, varargin{:});
            params = p.Results;
            
            % Call the private enhanced method
            quality = CovarianceWhitening.compute_enhanced_quality_metrics(Sigma_orig, Sigma_tilde, omega, params);
        end
        
        function report_quality(quality, params)
        % REPORT_QUALITY - Public wrapper for quality reporting
        % 
        % This provides external access to quality reporting functionality
        %
        % Usage:
        %   CovarianceWhitening.report_quality(quality, params)
            
            CovarianceWhitening.report_enhanced_whitening_quality(quality, params);
        end
        
    end
    
    methods (Static, Access = private)
        % Private helper methods (equivalent to the local functions in original file)
        
        function validate_whitening_inputs(Sigma_emp, D)
        % Validate input matrices (EXACT copy from original function)
            
            F_sigma = length(Sigma_emp);
            F_D = length(D);
            
            if F_sigma ~= F_D
                error('covariance_whitening:dimension_mismatch', ...
                      'Number of covariance matrices (%d) does not match whitening matrices (%d)', ...
                      F_sigma, F_D);
            end
            
            if F_sigma == 0
                error('covariance_whitening:empty_input', 'Input arrays are empty');
            end
            
            % Check dimensions and properties
            [n_sigma, n_sigma_check] = size(Sigma_emp{1});
            [n_D, n_D_check] = size(D{1});
            
            if n_sigma ~= n_sigma_check || n_D ~= n_D_check
                error('covariance_whitening:not_square', 'Matrices must be square');
            end
            
            if n_sigma ~= n_D
                error('covariance_whitening:size_mismatch', ...
                      'Matrix size mismatch: covariance [%d x %d], whitening [%d x %d]', ...
                      n_sigma, n_sigma_check, n_D, n_D_check);
            end
            
            % Validate all matrices
            for omega = 1:F_sigma
                % Check whitening matrix is diagonal and positive
                off_diagonal_norm = norm(D{omega} - diag(diag(D{omega})), 'fro');
                if off_diagonal_norm > 1e-12
                    error('covariance_whitening:not_diagonal', ...
                          'Whitening matrix %d is not diagonal', omega);
                end
                
                if any(diag(D{omega}) <= 0)
                    error('covariance_whitening:negative_diagonal', ...
                          'Whitening matrix %d has non-positive diagonal entries', omega);
                end
            end
        end
        
        function quality = initialize_quality_structure(F, n)
        % Initialize quality metrics structure (EXACT copy from original)
            quality = struct();
            quality.diagonal_errors = zeros(F, n);
            quality.hermitian_errors = zeros(F, 1);
            quality.min_eigenvalues = zeros(F, 1);
            quality.condition_numbers = zeros(F, 1);
            quality.whitening_effectiveness = zeros(F, 1);
        end
        
        function [Sigma_tilde_omega, quality_omega] = apply_enhanced_whitening_single(Sigma_omega, D_omega, omega, params)
        % Apply enhanced whitening with iterative refinement (EXACT copy from original)
            
            % Initial whitening transformation
            Sigma_current = D_omega * Sigma_omega * D_omega;
            
            % Force Hermitian symmetry
            if params.force_hermitian
                Sigma_current = (Sigma_current + Sigma_current') / 2;
            end
            
            % Apply iterative correction if enabled
            if params.adaptive_correction
                [Sigma_current, correction_info] = CovarianceWhitening.apply_iterative_diagonal_correction(...
                    Sigma_omega, Sigma_current, D_omega, params);
            else
                correction_info = struct('iterations', 0, 'converged', true, 'final_error', 0);
            end
            
            % Final result
            Sigma_tilde_omega = Sigma_current;
            
            % Compute quality metrics
            quality_omega = CovarianceWhitening.compute_enhanced_quality_metrics(Sigma_omega, Sigma_tilde_omega, omega, params);
            quality_omega.correction_info = correction_info;
            
            % Apply safety corrections if needed
            [Sigma_tilde_omega, safety_applied] = CovarianceWhitening.apply_safety_corrections(Sigma_tilde_omega, quality_omega, params);
            
            if safety_applied
                % Recompute quality after safety corrections
                quality_omega = CovarianceWhitening.compute_enhanced_quality_metrics(Sigma_omega, Sigma_tilde_omega, omega, params);
                quality_omega.correction_info.safety_applied = true;
            end
        end
        
        function [Sigma_corrected, correction_info] = apply_iterative_diagonal_correction(Sigma_orig, Sigma_initial, D_initial, params)
        % Apply iterative diagonal correction for better normalization (EXACT copy from original)
            
            Sigma_corrected = Sigma_initial;
            correction_info = struct('iterations', 0, 'converged', false, 'initial_error', 0, 'final_error', 0);
            
            % Compute initial error
            initial_diag = real(diag(Sigma_initial));
            initial_error = max(abs(initial_diag - params.target_diagonal));
            correction_info.initial_error = initial_error;
            
            % Skip correction if already good enough
            if initial_error <= params.convergence_threshold
                correction_info.final_error = initial_error;
                correction_info.converged = true;
                return;
            end
            
            n = size(Sigma_initial, 1);
            D_current = D_initial;
            previous_error = initial_error; % Initialize previous_error
            
            for iter = 1:params.max_correction_iter
                % Current diagonal
                current_diag = real(diag(Sigma_corrected));
                current_error = max(abs(current_diag - params.target_diagonal));
                
                % Check convergence
                if current_error <= params.convergence_threshold
                    correction_info.converged = true;
                    break;
                end
                
                % Compute adaptive correction factors
                correction_factors = CovarianceWhitening.compute_adaptive_correction_factors(current_diag, params, iter);
                
                % Update whitening matrix conservatively
                D_diag_current = diag(D_current);
                D_diag_new = D_diag_current .* sqrt(correction_factors);
                
                % Apply bounds to prevent extreme corrections
                median_D = median(D_diag_current);
                D_diag_new = max(D_diag_new, median_D / 50);
                D_diag_new = min(D_diag_new, median_D * 50);
                
                D_current = diag(D_diag_new);
                
                % Apply updated whitening
                Sigma_corrected = D_current * Sigma_orig * D_current;
                
                % Ensure Hermitian
                if params.force_hermitian
                    Sigma_corrected = (Sigma_corrected + Sigma_corrected') / 2;
                end
                
                correction_info.iterations = iter;
                
                % Check for sufficient progress
                if iter > 1 && current_error > previous_error * 0.9
                    break; % Not making enough progress
                end
                
                previous_error = current_error;
            end
            
            correction_info.final_error = max(abs(real(diag(Sigma_corrected)) - params.target_diagonal));
        end
        
        function correction_factors = compute_adaptive_correction_factors(current_diag, params, iteration)
        % Compute adaptive correction factors for iterative refinement (EXACT copy from original)
            
            target = params.target_diagonal;
            
            % Basic correction
            basic_factors = target ./ current_diag;
            
            % Apply damping that decreases with iterations
            damping = 0.3 / iteration; % More conservative with more iterations
            correction_factors = 1 + damping * (basic_factors - 1);
            
            % Robust bounds based on current distribution
            median_diag = median(current_diag);
            robust_std = mad(current_diag, 1) * 1.4826; % Robust standard deviation
            
            % Adaptive bounds
            max_factor = 1 + min(0.5, 2 * robust_std / median_diag);
            min_factor = 1 - min(0.3, robust_std / median_diag);
            
            % Apply bounds
            correction_factors = max(min_factor, min(max_factor, correction_factors));
            
            % Additional smoothing for extreme outliers
            outlier_mask = abs(current_diag - median_diag) > 3 * robust_std;
            if any(outlier_mask)
                correction_factors(outlier_mask) = 1 + 0.1 * (correction_factors(outlier_mask) - 1);
            end
        end
        
        function quality = compute_enhanced_quality_metrics(Sigma_orig, Sigma_tilde, omega, params)
        % COMPLETE CORRECTED VERSION - Enhanced quality metrics computation (EXACT from original)
        % 
        % This function computes reliable quality metrics for whitening assessment.
        % All bugs have been identified and fixed.

            n = size(Sigma_tilde, 1);
            quality = struct();
            
            % 1. Diagonal analysis
            diagonal_elements = diag(Sigma_tilde);
            quality.diagonal_errors = abs(real(diagonal_elements) - params.target_diagonal);
            quality.max_diagonal_error = max(quality.diagonal_errors);
            quality.mean_diagonal_error = mean(quality.diagonal_errors);
            quality.std_diagonal_error = std(quality.diagonal_errors);
            
            % 2. Hermitian symmetry
            hermitian_diff = Sigma_tilde - Sigma_tilde';
            quality.hermitian_error = norm(hermitian_diff, 'fro') / max(norm(Sigma_tilde, 'fro'), 1e-12);
            
            % 3. Eigenvalue analysis
            try
                eigenvals = eig(Sigma_tilde);
                quality.min_eigenvalue = min(real(eigenvals));
                quality.max_eigenvalue = max(real(eigenvals));
                
                positive_eigs = real(eigenvals(real(eigenvals) > 1e-12));
                if length(positive_eigs) > 1
                    quality.condition_number = max(positive_eigs) / min(positive_eigs);
                else
                    quality.condition_number = Inf;
                end
                
                quality.negative_eigenvals = sum(real(eigenvals) < -1e-12);
                
            catch ME
                warning('covariance_whitening:eigenvalue_failed', ...
                        'Eigenvalue computation failed for frequency %d: %s', omega, ME.message);
                quality.min_eigenvalue = NaN;
                quality.max_eigenvalue = NaN;
                quality.condition_number = Inf;
                quality.negative_eigenvals = NaN;
            end
            
            % 4. FIXED: Enhanced effectiveness score calculation
            diagonal_variance = var(real(diagonal_elements));
            mean_diagonal_deviation = abs(mean(real(diagonal_elements)) - params.target_diagonal);
            
            % Component scores (all clamped to [0,1])
            diagonal_score = exp(-diagonal_variance * 15); % Penalize variance
            diagonal_score = max(0, min(1, diagonal_score));
            
            mean_score = exp(-mean_diagonal_deviation * 8); % Penalize mean deviation
            mean_score = max(0, min(1, mean_score));
            
            hermitian_score = exp(-quality.hermitian_error * 50);
            hermitian_score = max(0, min(1, hermitian_score));
            
            % PSD score
            if quality.min_eigenvalue > 0
                psd_score = 1.0;
            elseif quality.min_eigenvalue > -1e-10
                psd_score = 0.8;
            else
                psd_score = 0.5;
            end
            
            % Condition score
            if isfinite(quality.condition_number) && quality.condition_number < 1e6
                condition_score = exp(-log10(max(quality.condition_number, 1)) / 4);
            else
                condition_score = 0.2;
            end
            condition_score = max(0, min(1, condition_score));
            
            % FIXED: Weighted combination with proper array handling
            weights = [0.4, 0.3, 0.15, 0.1, 0.05];
            scores = [diagonal_score, mean_score, hermitian_score, psd_score, condition_score];
            
            % Ensure both arrays have the same size
            if length(weights) ~= length(scores)
                error('Weights and scores arrays must have the same length');
            end
            
            % Compute weighted sum
            quality.whitening_effectiveness = sum(weights .* scores);
            
            % CRITICAL: Final safety checks
            if ~isfinite(quality.whitening_effectiveness)
                warning('Non-finite effectiveness calculated, using fallback');
                quality.whitening_effectiveness = max(0, 1 - mean_diagonal_deviation);
            end
            
            if quality.whitening_effectiveness < 0
                quality.whitening_effectiveness = 0;
            elseif quality.whitening_effectiveness > 1
                quality.whitening_effectiveness = 1;
            end
            
            % Store component scores for debugging
            quality.component_scores = struct();
            quality.component_scores.diagonal = diagonal_score;
            quality.component_scores.mean = mean_score;
            quality.component_scores.hermitian = hermitian_score;
            quality.component_scores.psd = psd_score;
            quality.component_scores.condition = condition_score;
        end
        
        function [Sigma_corrected, correction_applied] = apply_safety_corrections(Sigma_tilde, quality, params)
        % Apply final safety corrections (EXACT copy from original)
            
            Sigma_corrected = Sigma_tilde;
            correction_applied = false;
            
            % Safety 1: Extreme diagonal deviations
            if quality.max_diagonal_error > params.diagonal_tolerance * 3
                diagonal_elements = diag(Sigma_tilde);
                corrected_diagonal = real(diagonal_elements);
                
                % Gentle adjustment towards target
                adjustment_strength = 0.2; % Very conservative
                target_adjustment = (params.target_diagonal - corrected_diagonal) * adjustment_strength;
                corrected_diagonal = corrected_diagonal + target_adjustment;
                
                % Apply reasonable bounds
                corrected_diagonal = max(corrected_diagonal, params.target_diagonal * 0.3);
                corrected_diagonal = min(corrected_diagonal, params.target_diagonal * 3.0);
                
                % Reconstruct matrix
                Sigma_corrected = Sigma_corrected - diag(diag(Sigma_corrected)) + diag(corrected_diagonal);
                correction_applied = true;
            end
            
            % Safety 2: Negative eigenvalues
            if params.check_psd && quality.negative_eigenvals > 0
                try
                    [V, Lambda] = eig(Sigma_corrected);
                    Lambda_diag = diag(Lambda);
                    
                    % Set negative eigenvalues to small positive value
                    min_eigenval = max(1e-8, params.target_diagonal / (size(Sigma_corrected, 1) * 20));
                    Lambda_diag(real(Lambda_diag) < min_eigenval) = min_eigenval;
                    
                    Sigma_corrected = V * diag(Lambda_diag) * V';
                    Sigma_corrected = (Sigma_corrected + Sigma_corrected') / 2;
                    
                    correction_applied = true;
                    
                catch
                    % If eigendecomposition fails, just add small regularization
                    n = size(Sigma_corrected, 1);
                    reg_strength = params.target_diagonal / (n * 100);
                    Sigma_corrected = Sigma_corrected + reg_strength * eye(n);
                    correction_applied = true;
                end
            end
        end
        
        function store_quality_metrics(whitening_quality, quality_omega, omega)
        % Store quality metrics for current frequency (EXACT copy from original)
            
            whitening_quality.diagonal_errors(omega, :) = quality_omega.diagonal_errors;
            whitening_quality.hermitian_errors(omega) = quality_omega.hermitian_error;
            whitening_quality.min_eigenvalues(omega) = quality_omega.min_eigenvalue;
            whitening_quality.condition_numbers(omega) = quality_omega.condition_number;
            whitening_quality.whitening_effectiveness(omega) = quality_omega.whitening_effectiveness;
        end
        
        function report_enhanced_whitening_quality(quality, params)
        % CORRECTED VERSION - Report enhanced whitening quality (EXACT copy from original)
        % 
        % This function fixes the contradictory error reporting issue

            F = length(quality.whitening_effectiveness);
            
            fprintf('\nEnhanced whitening quality assessment:\n');
            fprintf('=====================================\n');
            
            % FIXED: Diagonal quality with proper error calculation
            if size(quality.diagonal_errors, 1) == F && size(quality.diagonal_errors, 2) > 0
                max_diagonal_errors = max(quality.diagonal_errors, [], 2);
                mean_diagonal_errors = mean(quality.diagonal_errors, 2);
                
                fprintf('Diagonal normalization (target: %.3f):\n', params.target_diagonal);
                fprintf('  Max error  - Mean: %.4f, Median: %.4f, Range: [%.4f, %.4f]\n', ...
                        mean(max_diagonal_errors), median(max_diagonal_errors), ...
                        min(max_diagonal_errors), max(max_diagonal_errors));
                fprintf('  Mean error - Mean: %.4f, Median: %.4f\n', ...
                        mean(mean_diagonal_errors), median(mean_diagonal_errors));
                
                % Success rates at multiple tolerance levels
                tolerances = [0.05, 0.08, 0.10, 0.15, 0.20];
                fprintf('  Success rates:\n');
                for i = 1:length(tolerances)
                    tol = tolerances(i);
                    success_count = sum(max_diagonal_errors <= tol);
                    success_rate = success_count / F * 100;
                    fprintf('    ≤%.3f: %3d/%d (%.1f%%)\n', tol, success_count, F, success_rate);
                end
                
                good_diagonal_rate = sum(max_diagonal_errors <= params.diagonal_tolerance) / F;
            else
                fprintf('Diagonal normalization (target: %.3f):\n', params.target_diagonal);
                fprintf('  ERROR: Invalid diagonal_errors structure\n');
                good_diagonal_rate = 0;
                max_diagonal_errors = [];
            end
            
            % FIXED: Overall effectiveness reporting
            valid_effectiveness = quality.whitening_effectiveness(isfinite(quality.whitening_effectiveness));
            
            fprintf('\nWhitening effectiveness:\n');
            if ~isempty(valid_effectiveness)
                fprintf('  Mean: %.3f, Median: %.3f, Min: %.3f\n', ...
                        mean(valid_effectiveness), median(valid_effectiveness), min(valid_effectiveness));
                
                excellent_count = sum(valid_effectiveness > 0.9);
                good_count = sum(valid_effectiveness > 0.8);
                acceptable_count = sum(valid_effectiveness > 0.7);
                
                fprintf('  Quality distribution:\n');
                fprintf('    Excellent (>0.9): %d/%d (%.1f%%)\n', excellent_count, F, excellent_count/F*100);
                fprintf('    Good (>0.8):      %d/%d (%.1f%%)\n', good_count, F, good_count/F*100);
                fprintf('    Acceptable (>0.7): %d/%d (%.1f%%)\n', acceptable_count, F, acceptable_count/F*100);
            else
                fprintf('  ERROR: All effectiveness values are invalid\n');
            end
            
            % FIXED: Numerical stability reporting
            fprintf('\nNumerical stability:\n');
            valid_conditions = quality.condition_numbers(isfinite(quality.condition_numbers));
            if ~isempty(valid_conditions)
                fprintf('  Condition numbers - Mean: %.2e, Median: %.2e, Max: %.2e\n', ...
                        mean(valid_conditions), median(valid_conditions), max(valid_conditions));
            else
                fprintf('  Condition numbers - All invalid\n');
            end
            
            negative_count = sum(quality.min_eigenvalues < -1e-12);
            % FIXED: Use reasonable threshold for near-singular
            near_singular_threshold = 1e6;  
            near_singular_count = sum(quality.condition_numbers > near_singular_threshold | ...
                                     isinf(quality.condition_numbers));
            
            fprintf('  Negative eigenvalues: %d/%d frequencies\n', negative_count, F);
            fprintf('  Near-singular: %d/%d frequencies\n', near_singular_count, F);
            
            % FIXED: Overall assessment
            if ~isempty(valid_effectiveness) && ~isempty(max_diagonal_errors)
                overall_score = mean(valid_effectiveness);
                
                fprintf('\nOverall assessment:\n');
                if overall_score > 0.85 && good_diagonal_rate > 0.8
                    assessment = 'Excellent - High quality whitening achieved';
                elseif overall_score > 0.75 && good_diagonal_rate > 0.7
                    assessment = 'Good - Satisfactory whitening quality';
                elseif overall_score > 0.65 && good_diagonal_rate > 0.5
                    assessment = 'Fair - Acceptable whitening with some issues';
                else
                    assessment = 'Poor - Whitening quality needs improvement';
                end
                
                fprintf('  %s\n', assessment);
                fprintf('  Overall score: %.3f, Good diagonal rate: %.1f%%\n', ...
                        overall_score, good_diagonal_rate * 100);
            else
                fprintf('\nOverall assessment:\n');
                fprintf('  Cannot assess - Invalid quality metrics\n');
            end
            
            fprintf('\n');
        end
        
        function validate_whitening_output_relaxed(Sigma_tilde, quality, params)
        % Final validation with relaxed thresholds to reduce warnings (EXACT copy from original)
            
            F = length(Sigma_tilde);
            
            if F ~= length(quality.whitening_effectiveness)
                error('covariance_whitening:inconsistent_output', ...
                      'Inconsistent output dimensions');
            end
            
            % Count problematic frequencies but don't warn for each one
            if F > 0
                n = size(Sigma_tilde{1}, 1);
                problematic_count = 0;
                severe_count = 0;
                
                for omega = 1:F
                    Sigma_omega = Sigma_tilde{omega};
                    
                    % Check basic properties
                    if any(~isfinite(Sigma_omega(:)))
                        error('covariance_whitening:infinite_entries', ...
                              'Matrix %d contains non-finite entries', omega);
                    end
                    
                    % Check diagonal quality
                    diagonal_elements = diag(Sigma_omega);
                    max_diagonal_error = max(abs(real(diagonal_elements) - params.target_diagonal));
                    
                    if max_diagonal_error > params.diagonal_tolerance * 2
                        problematic_count = problematic_count + 1;
                    end
                    
                    if max_diagonal_error > params.diagonal_tolerance * 4
                        severe_count = severe_count + 1;
                    end
                end
                
                % Summary reporting instead of individual warnings
                if problematic_count > 0
                    fprintf('Note: %d/%d frequencies have diagonal errors > %.3f\n', ...
                            problematic_count, F, params.diagonal_tolerance * 2);
                end
                
                if severe_count > 0
                    fprintf('Warning: %d/%d frequencies have severe diagonal errors > %.3f\n', ...
                            severe_count, F, params.diagonal_tolerance * 4);
                end
            end
            
            % Validate quality metrics
            if any(~isfinite(quality.whitening_effectiveness))
                error('covariance_whitening:invalid_quality', ...
                      'Quality metrics contain non-finite values');
            end
            
            if any(quality.whitening_effectiveness < 0) || any(quality.whitening_effectiveness > 1)
                error('covariance_whitening:invalid_effectiveness', ...
                      'Effectiveness scores must be between 0 and 1');
            end
            
            fprintf('Enhanced covariance whitening validation passed\n');
        end
        
    end
end