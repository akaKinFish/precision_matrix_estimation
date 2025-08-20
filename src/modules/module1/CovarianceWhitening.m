classdef CovarianceWhitening < handle
    % COVARIANCE_WHITENING - Enhanced covariance whitening with quality control
    %
    % This class performs covariance matrix whitening with comprehensive
    % quality assessment and diagonal loading capabilities.
    %
    % Key features:
    % - Enhanced diagonal normalization
    % - Comprehensive quality metrics
    % - Safety corrections and validation
    % - Support for complex-valued matrices
    %
    % Usage:
    %   whitener = CovarianceWhitening();
    %   [Sigma_tilde, quality] = whitener.apply_whitening(Sigma_emp, D);
    %
    % See also: module1_preprocessing_main, construct_whitening_matrices
    %
    % Author: Enhanced Module 1 Team
    % Date: August 2025
    % Version: 3.0 (Fully Fixed)
    
    properties (Access = private)
        verbose = true;
    end
    
    methods (Access = public)
        function obj = CovarianceWhitening(varargin)
            % Constructor with parameter parsing
            p = inputParser;
            addParameter(p, 'verbose', true, @islogical);
            parse(p, varargin{:});
            
            obj.verbose = p.Results.verbose;
        end
        
        function [Sigma_tilde, quality_stats] = apply_whitening(obj, Sigma_emp, D, varargin)
            % APPLY_WHITENING - Apply covariance whitening with quality assessment
            %
            % Syntax:
            %   [Sigma_tilde, quality] = apply_whitening(Sigma_emp, D)
            %   [Sigma_tilde, quality] = apply_whitening(..., Name, Value)
            %
            % Input Arguments:
            %   Sigma_emp - Cell array of empirical covariance matrices
            %   D - Cell array of whitening matrices (diagonal)
            %
            % Name-Value Arguments:
            %   target_diagonal - Target diagonal value (default: 1.0)
            %   diagonal_tolerance - Tolerance for diagonal quality (default: 0.1)
            %
            % Output Arguments:
            %   Sigma_tilde - Cell array of whitened covariance matrices
            %   quality_stats - Structure with quality metrics
            
            % Parse parameters
            p = inputParser;
            addParameter(p, 'target_diagonal', 1.0, @(x) isnumeric(x) && isscalar(x) && x > 0);
            addParameter(p, 'diagonal_tolerance', 0.1, @(x) isnumeric(x) && isscalar(x) && x > 0);
            parse(p, varargin{:});
            
            params = p.Results;
            
            % Validate inputs
            obj.validate_whitening_inputs(Sigma_emp, D);
            
            F = length(Sigma_emp);
            n = size(Sigma_emp{1}, 1);
            
            % Initialize outputs
            Sigma_tilde = cell(F, 1);
            quality_stats = struct();
            
            % Initialize quality tracking arrays - FIXED
            diagonal_errors = cell(F, 1);
            max_diagonal_errors = zeros(F, 1);
            mean_diagonal_errors = zeros(F, 1);
            hermitian_errors = zeros(F, 1);
            min_eigenvalues = zeros(F, 1);
            condition_numbers = zeros(F, 1);
            whitening_effectiveness = zeros(F, 1);
            
            if obj.verbose
                fprintf('Applying enhanced whitening transformation to %d frequencies, %d nodes\n', F, n);
            end
            
            % Process each frequency
            for omega = 1:F
                try
                    % Apply whitening transformation
                    Sigma_tilde{omega} = D{omega} * Sigma_emp{omega} * D{omega};
                    
                    % Compute quality metrics for this frequency
                    freq_quality = obj.compute_quality_metrics(Sigma_tilde{omega}, params, omega);
                    
                    % Store quality metrics - FIXED array assignments
                    diagonal_errors{omega} = freq_quality.diagonal_errors;
                    max_diagonal_errors(omega) = freq_quality.max_diagonal_error;
                    mean_diagonal_errors(omega) = freq_quality.mean_diagonal_error;
                    hermitian_errors(omega) = freq_quality.hermitian_error;
                    min_eigenvalues(omega) = freq_quality.min_eigenvalue;
                    condition_numbers(omega) = freq_quality.condition_number;
                    whitening_effectiveness(omega) = freq_quality.whitening_effectiveness;
                    
                    % Apply safety corrections if needed
                    [Sigma_tilde{omega}, correction_applied] = obj.apply_safety_corrections(...
                        Sigma_tilde{omega}, freq_quality, params);
                    
                    if obj.verbose && mod(omega, max(1, floor(F/5))) == 0
                        fprintf('Processed frequency %d/%d\n', omega, F);
                    end
                    
                catch ME
                    error('covariance_whitening:processing_failed', ...
                          'Failed to process frequency %d: %s', omega, ME.message);
                end
            end
            
            % Compile overall quality statistics - ALL FIXED
            quality_stats.diagonal_errors = diagonal_errors;
            quality_stats.max_diagonal_errors = max_diagonal_errors;
            quality_stats.mean_diagonal_errors = mean_diagonal_errors;
            quality_stats.hermitian_errors = hermitian_errors;
            quality_stats.min_eigenvalues = min_eigenvalues;
            quality_stats.condition_numbers = condition_numbers;
            quality_stats.whitening_effectiveness = whitening_effectiveness;
            
            % Compute summary statistics - FIXED to handle arrays properly
            quality_stats.overall_max_error = max(max_diagonal_errors);
            quality_stats.overall_mean_error = mean(mean_diagonal_errors);
            quality_stats.success_rates = struct();
            
            % Success rates for different tolerances - FIXED
            tolerances = [0.050, 0.080, 0.100, 0.150, 0.200];
            for i = 1:length(tolerances)
                tol = tolerances(i);
                success_count = sum(max_diagonal_errors <= tol);
                field_name = sprintf('tol_%03d', round(tol * 1000));
                quality_stats.success_rates.(field_name) = success_count / F;
            end
            
            % Overall assessment - FIXED
            mean_effectiveness = mean(whitening_effectiveness);
            good_diagonal_rate = quality_stats.success_rates.tol_200;
            
            if mean_effectiveness > 0.8 && good_diagonal_rate > 0.9
                assessment = 'Excellent';
            elseif mean_effectiveness > 0.6 && good_diagonal_rate > 0.8
                assessment = 'Good';
            elseif mean_effectiveness > 0.4 && good_diagonal_rate > 0.7
                assessment = 'Acceptable';
            else
                assessment = 'Poor';
            end
            
            quality_stats.overall_assessment = assessment;
            quality_stats.mean_effectiveness = mean_effectiveness;
            quality_stats.good_diagonal_rate = good_diagonal_rate;
            
            % Validation and reporting
            if obj.verbose
                obj.perform_enhanced_validation(Sigma_tilde, quality_stats, params);
                obj.report_enhanced_whitening_quality(quality_stats, params);
            end
        end
        
    end
    
    methods (Access = private)
        
        function validate_whitening_inputs(obj, Sigma_emp, D)
            % Validate input matrices
            
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
                          'Whitening matrix %d has non-positive diagonal elements', omega);
                end
            end
        end
        
        function quality = compute_quality_metrics(obj, Sigma_tilde, params, omega)
            % Compute comprehensive quality metrics for a single matrix
            
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
                if obj.verbose
                    warning('covariance_whitening:eigenvalue_failed', ...
                            'Eigenvalue computation failed for frequency %d: %s', omega, ME.message);
                end
                quality.min_eigenvalue = NaN;
                quality.max_eigenvalue = NaN;
                quality.condition_number = Inf;
                quality.negative_eigenvals = NaN;
            end
            
            % 4. Whitening effectiveness calculation - COMPLETELY FIXED
            diagonal_variance = var(real(diagonal_elements));
            mean_diagonal_deviation = abs(mean(real(diagonal_elements)) - params.target_diagonal);
            
            % Component scores (all clamped to [0,1])
            diagonal_score = exp(-diagonal_variance * 15);
            diagonal_score = max(0, min(1, diagonal_score));
            
            mean_score = exp(-mean_diagonal_deviation * 8);
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
            
            % COMPLETELY FIXED: Proper scalar computation
            weights = [0.4, 0.3, 0.15, 0.1, 0.05];
            scores = [diagonal_score, mean_score, hermitian_score, psd_score, condition_score];
            
            % Compute weighted sum as scalar
            quality.whitening_effectiveness = 0;
            for i = 1:length(weights)
                quality.whitening_effectiveness = quality.whitening_effectiveness + weights(i) * scores(i);
            end
            
            % Final safety checks
            if ~isfinite(quality.whitening_effectiveness)
                if obj.verbose
                    warning('Non-finite effectiveness calculated, using fallback');
                end
                quality.whitening_effectiveness = max(0, 1 - mean_diagonal_deviation);
            end
            
            quality.whitening_effectiveness = max(0, min(1, quality.whitening_effectiveness));
            
            % Store component scores for debugging
            quality.component_scores = struct();
            quality.component_scores.diagonal = diagonal_score;
            quality.component_scores.mean = mean_score;
            quality.component_scores.hermitian = hermitian_score;
            quality.component_scores.psd = psd_score;
            quality.component_scores.condition = condition_score;
        end
        
        function [Sigma_corrected, correction_applied] = apply_safety_corrections(obj, Sigma_tilde, quality, params)
            % Apply final safety corrections
            
            Sigma_corrected = Sigma_tilde;
            correction_applied = false;
            
            % Safety 1: Extreme diagonal deviations
            if quality.max_diagonal_error > params.diagonal_tolerance * 3
                diagonal_elements = diag(Sigma_tilde);
                corrected_diagonal = real(diagonal_elements);
                
                % Gentle adjustment towards target
                adjustment_strength = 0.2;
                target_adjustment = (params.target_diagonal - corrected_diagonal) * adjustment_strength;
                corrected_diagonal = corrected_diagonal + target_adjustment;
                
                % Apply reasonable bounds
                corrected_diagonal = max(corrected_diagonal, params.target_diagonal * 0.3);
                corrected_diagonal = min(corrected_diagonal, params.target_diagonal * 3.0);
                
                % Reconstruct matrix
                off_diagonal = Sigma_tilde - diag(diag(Sigma_tilde));
                Sigma_corrected = off_diagonal + diag(corrected_diagonal);
                correction_applied = true;
            end
            
            % Safety 2: Enforce Hermitian property
            if quality.hermitian_error > 1e-10
                Sigma_corrected = (Sigma_corrected + Sigma_corrected') / 2;
                correction_applied = true;
            end
        end
        
        function perform_enhanced_validation(obj, Sigma_tilde, quality, params)
            % Perform comprehensive validation
            
            F = length(Sigma_tilde);
            problematic_count = 0;
            severe_count = 0;
            
            % Batch validation to reduce output
            for omega = 1:F
                Sigma_omega = Sigma_tilde{omega};
                
                % Check basic properties
                if any(~isfinite(Sigma_omega(:)))
                    error('covariance_whitening:infinite_entries', ...
                          'Matrix %d contains non-finite entries', omega);
                end
                
                % Check diagonal quality
                max_diagonal_error = quality.max_diagonal_errors(omega);
                
                if max_diagonal_error > params.diagonal_tolerance * 2
                    problematic_count = problematic_count + 1;
                end
                
                if max_diagonal_error > params.diagonal_tolerance * 4
                    severe_count = severe_count + 1;
                end
            end
            
            % Summary reporting
            if problematic_count > 0
                fprintf('Note: %d/%d frequencies have diagonal errors > %.3f\n', ...
                        problematic_count, F, params.diagonal_tolerance * 2);
            end
            
            if severe_count > 0
                fprintf('Warning: %d/%d frequencies have severe diagonal errors > %.3f\n', ...
                        severe_count, F, params.diagonal_tolerance * 4);
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
        
        function report_enhanced_whitening_quality(obj, quality, params)
            % Generate comprehensive quality report
            
            F = length(quality.whitening_effectiveness);
            
            fprintf('\nEnhanced whitening quality assessment:\n');
            fprintf('=====================================\n');
            
            % Diagonal normalization analysis
            fprintf('Diagonal normalization (target: %.3f):\n', params.target_diagonal);
            fprintf('  Max error  - Mean: %.4f, Median: %.4f, Range: [%.4f, %.4f]\n', ...
                    mean(quality.max_diagonal_errors), median(quality.max_diagonal_errors), ...
                    min(quality.max_diagonal_errors), max(quality.max_diagonal_errors));
            fprintf('  Mean error - Mean: %.4f, Median: %.4f\n', ...
                    mean(quality.mean_diagonal_errors), median(quality.mean_diagonal_errors));
            
            % Success rates
            fprintf('  Success rates:\n');
            rate_fields = fieldnames(quality.success_rates);
            for i = 1:length(rate_fields)
                field = rate_fields{i};
                tolerance = str2double(field(5:end)) / 1000;
                rate = quality.success_rates.(field);
                count = round(rate * F);
                fprintf('    ≤%.3f:  %d/%d (%.1f%%)\n', tolerance, count, F, rate * 100);
            end
            
            % Whitening effectiveness reporting - FIXED
            valid_effectiveness = quality.whitening_effectiveness;
            finite_mask = isfinite(valid_effectiveness);
            valid_effectiveness = valid_effectiveness(finite_mask);
            
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
            
            % Numerical stability reporting - FIXED
            fprintf('\nNumerical stability:\n');
            valid_conditions = quality.condition_numbers(isfinite(quality.condition_numbers));
            if ~isempty(valid_conditions)
                fprintf('  Condition numbers - Mean: %.2e, Median: %.2e, Max: %.2e\n', ...
                        mean(valid_conditions), median(valid_conditions), max(valid_conditions));
            else
                fprintf('  Condition numbers - All invalid\n');
            end
            
            negative_count = sum(quality.min_eigenvalues < -1e-12);
            near_singular_threshold = 1e6;
            near_singular_count = sum(quality.condition_numbers > near_singular_threshold | ...
                                     isinf(quality.condition_numbers));
            
            fprintf('  Negative eigenvalues: %d/%d frequencies\n', negative_count, F);
            fprintf('  Near-singular: %d/%d frequencies\n', near_singular_count, F);
            
            % Overall assessment
            fprintf('\nOverall assessment:\n');
            fprintf('  %s - %s\n', quality.overall_assessment, obj.get_assessment_description(quality.overall_assessment));
            fprintf('  Overall score: %.3f, Good diagonal rate: %.1f%%\n', ...
                    quality.mean_effectiveness, quality.good_diagonal_rate * 100);
        end
        
        function description = get_assessment_description(obj, assessment)
            % Get description for assessment level
            switch assessment
                case 'Excellent'
                    description = 'Whitening quality is excellent';
                case 'Good'
                    description = 'Whitening quality is good';
                case 'Acceptable'
                    description = 'Whitening quality is acceptable';
                case 'Poor'
                    description = 'Whitening quality needs improvement';
                otherwise
                    description = 'Unknown assessment level';
            end
        end
        
    end
end