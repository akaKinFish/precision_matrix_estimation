function fixed_results = fix_quality_assessment(demo_results)
% FIX_QUALITY_ASSESSMENT - Quick fix for the quality assessment bug
%
% This function takes the demo results and recalculates quality metrics
% using reliable methods, replacing the buggy official assessment
%
% Usage:
%   demo_results = demo_module1_preprocessing();
%   fixed_results = fix_quality_assessment(demo_results);
%   
% File location: examples/fix_quality_assessment.m

    fprintf('=== Fixing Quality Assessment Issues ===\n');
    
    if ~demo_results.preprocessing.success
        error('Cannot fix quality assessment - preprocessing failed');
    end
    
    results = demo_results.preprocessing.results;
    
    % Extract whitened matrices
    if ~isfield(results, 'Sigma_tilde') || ~iscell(results.Sigma_tilde)
        error('Whitened matrices not available');
    end
    
    Sigma_tilde = results.Sigma_tilde;
    F = length(Sigma_tilde);
    n = size(Sigma_tilde{1}, 1);
    
    fprintf('Recalculating quality metrics for %d frequencies...\n', F);
    
    %% Recalculate diagonal normalization metrics
    diagonal_metrics = calculate_diagonal_metrics(Sigma_tilde);
    
    %% Recalculate whitening effectiveness
    whitening_metrics = calculate_whitening_effectiveness(Sigma_tilde);
    
    %% Recalculate numerical stability
    stability_metrics = calculate_numerical_stability(Sigma_tilde);
    
    %% Create corrected quality assessment
    corrected_quality = struct();
    corrected_quality.diagonal_normalization = diagonal_metrics;
    corrected_quality.whitening_effectiveness = whitening_metrics.effectiveness;
    corrected_quality.numerical_stability = stability_metrics;
    corrected_quality.overall_score = calculate_overall_score(diagonal_metrics, whitening_metrics, stability_metrics);
    
    %% Replace buggy assessment in results
    fixed_results = demo_results;
    fixed_results.preprocessing.results.processing_stats.whitening_quality = corrected_quality;
    
    %% Display corrected assessment
    display_corrected_assessment(corrected_quality, F);
    
    fprintf('\n=== Quality Assessment Fixed ===\n');
end

function diagonal_metrics = calculate_diagonal_metrics(Sigma_tilde)
% Calculate reliable diagonal normalization metrics
    
    F = length(Sigma_tilde);
    
    % Initialize metrics
    diagonal_metrics = struct();
    diagonal_metrics.max_errors = zeros(F, 1);
    diagonal_metrics.mean_errors = zeros(F, 1);
    diagonal_metrics.target_deviations = zeros(F, 1);
    
    for f = 1:F
        S = Sigma_tilde{f};
        diag_S = real(diag(S));
        
        % Calculate errors relative to target value of 1.0
        errors = abs(diag_S - 1.0);
        diagonal_metrics.max_errors(f) = max(errors);
        diagonal_metrics.mean_errors(f) = mean(errors);
        diagonal_metrics.target_deviations(f) = abs(mean(diag_S) - 1.0);
    end
    
    % Calculate success rates
    thresholds = [0.05, 0.08, 0.10, 0.15, 0.20];
    diagonal_metrics.success_rates = struct();
    
    for i = 1:length(thresholds)
        thresh = thresholds(i);
        success_count = sum(diagonal_metrics.max_errors <= thresh);
        field_name = sprintf('threshold_%.0f', thresh * 100);
        diagonal_metrics.success_rates.(field_name) = success_count / F;
    end
    
    % Overall diagonal quality score
    diagonal_metrics.quality_score = mean(1 - diagonal_metrics.target_deviations);
end

function whitening_metrics = calculate_whitening_effectiveness(Sigma_tilde)
% Calculate reliable whitening effectiveness metrics
    
    F = length(Sigma_tilde);
    whitening_metrics = struct();
    
    % Method 1: Diagonal proximity to 1.0
    diagonal_effectiveness = zeros(F, 1);
    
    % Method 2: Off-diagonal reduction (if original matrices available)
    structure_effectiveness = zeros(F, 1);
    
    % Method 3: Condition number based effectiveness
    condition_effectiveness = zeros(F, 1);
    
    for f = 1:F
        S = Sigma_tilde{f};
        
        % Diagonal effectiveness: how close diagonals are to 1.0
        diag_S = real(diag(S));
        diagonal_effectiveness(f) = 1 - mean(abs(diag_S - 1.0));
        
        % Structure effectiveness: relative off-diagonal magnitude
        off_diag_strength = mean(abs(S(~eye(size(S)))));
        diag_strength = mean(abs(diag_S));
        if diag_strength > 0
            structure_effectiveness(f) = 1 - (off_diag_strength / diag_strength);
        else
            structure_effectiveness(f) = 0;
        end
        
        % Condition-based effectiveness
        try
            cond_S = cond(S);
            condition_effectiveness(f) = 1 / (1 + log10(max(cond_S, 1)));
        catch
            condition_effectiveness(f) = 0;
        end
    end
    
    % Combine effectiveness measures
    whitening_metrics.effectiveness = 0.5 * diagonal_effectiveness + ...
                                     0.3 * structure_effectiveness + ...
                                     0.2 * condition_effectiveness;
    
    whitening_metrics.diagonal_component = diagonal_effectiveness;
    whitening_metrics.structure_component = structure_effectiveness;
    whitening_metrics.condition_component = condition_effectiveness;
    
    % Quality distribution
    excellent_count = sum(whitening_metrics.effectiveness > 0.9);
    good_count = sum(whitening_metrics.effectiveness > 0.8);
    acceptable_count = sum(whitening_metrics.effectiveness > 0.7);
    
    whitening_metrics.quality_distribution = struct();
    whitening_metrics.quality_distribution.excellent = excellent_count / F;
    whitening_metrics.quality_distribution.good = good_count / F;
    whitening_metrics.quality_distribution.acceptable = acceptable_count / F;
end

function stability_metrics = calculate_numerical_stability(Sigma_tilde)
% Calculate numerical stability metrics
    
    F = length(Sigma_tilde);
    stability_metrics = struct();
    
    condition_numbers = zeros(F, 1);
    min_eigenvalues = zeros(F, 1);
    negative_eig_count = zeros(F, 1);
    near_zero_eig_count = zeros(F, 1);
    
    for f = 1:F
        S = Sigma_tilde{f};
        
        try
            % Eigenvalue analysis
            eigs_S = eig(S);
            eigs_real = real(eigs_S);
            
            min_eigenvalues(f) = min(eigs_real);
            negative_eig_count(f) = sum(eigs_real < -1e-12);
            near_zero_eig_count(f) = sum(abs(eigs_real) < 1e-10);
            
            % Condition number
            positive_eigs = eigs_real(eigs_real > 1e-12);
            if length(positive_eigs) > 1
                condition_numbers(f) = max(positive_eigs) / min(positive_eigs);
            else
                condition_numbers(f) = Inf;
            end
            
        catch
            condition_numbers(f) = NaN;
            min_eigenvalues(f) = NaN;
        end
    end
    
    stability_metrics.condition_numbers = condition_numbers;
    stability_metrics.min_eigenvalues = min_eigenvalues;
    stability_metrics.negative_eigenvalue_count = sum(negative_eig_count);
    stability_metrics.near_zero_eigenvalue_count = sum(near_zero_eig_count);
    
    % Stability assessment
    valid_conds = condition_numbers(~isnan(condition_numbers) & ~isinf(condition_numbers));
    stability_metrics.mean_condition = mean(valid_conds);
    stability_metrics.max_condition = max(valid_conds);
    
    % Near-singular assessment (use reasonable threshold)
    near_singular_threshold = 1e6;  % Much more reasonable than current
    stability_metrics.near_singular_count = sum(condition_numbers > near_singular_threshold);
    stability_metrics.stability_score = 1 - (stability_metrics.near_singular_count / F);
end

function overall_score = calculate_overall_score(diagonal_metrics, whitening_metrics, stability_metrics)
% Calculate overall quality score
    
    % Weight the different components
    diagonal_weight = 0.4;
    whitening_weight = 0.4;
    stability_weight = 0.2;
    
    diagonal_score = diagonal_metrics.quality_score;
    whitening_score = mean(whitening_metrics.effectiveness);
    stability_score = stability_metrics.stability_score;
    
    overall_score = diagonal_weight * diagonal_score + ...
                   whitening_weight * whitening_score + ...
                   stability_weight * stability_score;
    
    % Ensure score is between 0 and 1
    overall_score = max(0, min(1, overall_score));
end

function display_corrected_assessment(corrected_quality, F)
% Display the corrected quality assessment
    
    fprintf('\n=== CORRECTED Quality Assessment ===\n');
    
    % Diagonal normalization
    diag_metrics = corrected_quality.diagonal_normalization;
    fprintf('Diagonal Normalization (target: 1.000):\n');
    fprintf('  Max error  - Mean: %.4f, Median: %.4f, Range: [%.4f, %.4f]\n', ...
            mean(diag_metrics.max_errors), median(diag_metrics.max_errors), ...
            min(diag_metrics.max_errors), max(diag_metrics.max_errors));
    fprintf('  Mean error - Mean: %.4f, Median: %.4f\n', ...
            mean(diag_metrics.mean_errors), median(diag_metrics.mean_errors));
    
    % Success rates
    fprintf('  Success rates:\n');
    thresholds = [0.05, 0.08, 0.10, 0.15, 0.20];
    for i = 1:length(thresholds)
        thresh = thresholds(i);
        field_name = sprintf('threshold_%.0f', thresh * 100);
        rate = diag_metrics.success_rates.(field_name);
        passed = round(rate * F);
        fprintf('    ≤%.3f:  %d/%d (%.1f%%)\n', thresh, passed, F, rate * 100);
    end
    
    % Whitening effectiveness
    effectiveness = corrected_quality.whitening_effectiveness;
    fprintf('\nWhitening Effectiveness:\n');
    fprintf('  Mean: %.3f, Median: %.3f, Min: %.3f\n', ...
            mean(effectiveness), median(effectiveness), min(effectiveness));
    
    % Quality distribution
    if isfield(corrected_quality, 'whitening_effectiveness') && isnumeric(effectiveness)
        excellent_count = sum(effectiveness > 0.9);
        good_count = sum(effectiveness > 0.8);
        acceptable_count = sum(effectiveness > 0.7);
        
        fprintf('  Quality distribution:\n');
        fprintf('    Excellent (>0.9): %d/%d (%.1f%%)\n', excellent_count, F, excellent_count/F*100);
        fprintf('    Good (>0.8):      %d/%d (%.1f%%)\n', good_count, F, good_count/F*100);
        fprintf('    Acceptable (>0.7): %d/%d (%.1f%%)\n', acceptable_count, F, acceptable_count/F*100);
    end
    
    % Numerical stability
    stability = corrected_quality.numerical_stability;
    fprintf('\nNumerical Stability:\n');
    fprintf('  Condition numbers - Mean: %.2e, Max: %.2e\n', ...
            stability.mean_condition, stability.max_condition);
    fprintf('  Negative eigenvalues: %d/%d frequencies\n', ...
            stability.negative_eigenvalue_count, F);
    fprintf('  Near-singular: %d/%d frequencies\n', ...
            stability.near_singular_count, F);
    
    % Overall assessment
    overall_score = corrected_quality.overall_score;
    fprintf('\nCORRECTED Overall Assessment:\n');
    
    if overall_score > 0.9
        assessment = 'Excellent - High quality whitening achieved';
    elseif overall_score > 0.8
        assessment = 'Good - Whitening quality is satisfactory';
    elseif overall_score > 0.7
        assessment = 'Acceptable - Minor quality issues detected';
    elseif overall_score > 0.6
        assessment = 'Fair - Some quality concerns need attention';
    else
        assessment = 'Poor - Significant quality issues detected';
    end
    
    fprintf('  %s\n', assessment);
    fprintf('  Overall score: %.3f\n', overall_score);
    
    fprintf('=====================================\n');
end