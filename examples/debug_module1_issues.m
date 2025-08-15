function debug_module1_issues()
% DEBUG_MODULE1_ISSUES - Diagnose Module 1 whitening quality issues
%
% This script helps diagnose why whitening quality assessment is 0 
% and near-singular issues
%
% Usage:
%   debug_module1_issues();
%
% File location: examples/debug_module1_issues.m

    fprintf('=== Module 1 Issue Debugging Script ===\n\n');
    
    %% 1. Re-run demo and collect detailed information
    fprintf('Step 1: Re-running demo and collecting debug info...\n');
    
    try
        demo_results = demo_module1_preprocessing();
        
        if ~demo_results.preprocessing.success
            fprintf('ERROR: Preprocessing failed, cannot continue debugging\n');
            fprintf('Error message: %s\n', demo_results.preprocessing.error);
            return;
        end
        
        results = demo_results.preprocessing.results;
        fprintf('SUCCESS: Demo ran successfully, starting detailed analysis...\n\n');
        
    catch ME
        fprintf('ERROR: Demo execution failed: %s\n', ME.message);
        return;
    end
    
    %% 2. Check data structure integrity
    fprintf('Step 2: Checking data structure integrity...\n');
    check_data_structure(results);
    
    %% 3. Analyze whitening matrix quality
    fprintf('\nStep 3: Analyzing whitening matrix quality...\n');
    analyze_whitening_matrices(results);
    
    %% 4. Diagnose conditioning issues
    fprintf('\nStep 4: Diagnosing conditioning and singularity issues...\n');
    diagnose_conditioning_issues(results);
    
    %% 5. Check quality assessment algorithm
    fprintf('\nStep 5: Checking quality assessment algorithm...\n');
    verify_quality_assessment(results);
    
    %% 6. Provide fix recommendations
    fprintf('\nStep 6: Generating fix recommendations...\n');
    generate_fix_recommendations(results);
    
    %% 7. Create diagnostic plots
    fprintf('\nStep 7: Creating diagnostic visualization...\n');
    create_diagnostic_plots(results);
    
    fprintf('\n=== Debugging Complete ===\n');
end

function check_data_structure(results)
% Check data structure integrity
    
    fprintf('Checking required fields...\n');
    
    required_fields = {'Sigma_emp', 'Sigma_tilde', 'g_smooth', 'D', 'processing_stats'};
    missing_fields = {};
    
    for i = 1:length(required_fields)
        if ~isfield(results, required_fields{i})
            missing_fields{end+1} = required_fields{i};
        else
            fprintf('  OK %s: exists\n', required_fields{i});
        end
    end
    
    if ~isempty(missing_fields)
        fprintf('  ERROR Missing fields: %s\n', strjoin(missing_fields, ', '));
    end
    
    % Check cell array dimensions
    if isfield(results, 'Sigma_tilde') && iscell(results.Sigma_tilde)
        F = length(results.Sigma_tilde);
        n = size(results.Sigma_tilde{1}, 1);
        fprintf('  INFO Data dimensions: %d frequencies x %dx%d matrices\n', F, n, n);
        
        % Check all matrix dimensions are consistent
        dim_consistent = true;
        for f = 1:F
            if any(size(results.Sigma_tilde{f}) ~= [n, n])
                dim_consistent = false;
                break;
            end
        end
        
        if dim_consistent
            fprintf('  OK Matrix dimensions consistent\n');
        else
            fprintf('  ERROR Matrix dimensions inconsistent\n');
        end
    end
end

function analyze_whitening_matrices(results)
% Analyze whitening matrices in detail
    
    if ~isfield(results, 'Sigma_tilde') || ~iscell(results.Sigma_tilde)
        fprintf('ERROR: Whitening matrix data not available\n');
        return;
    end
    
    Sigma_tilde = results.Sigma_tilde;
    F = length(Sigma_tilde);
    n = size(Sigma_tilde{1}, 1);
    
    fprintf('Analyzing %d whitening matrices...\n', F);
    
    % Detailed statistics
    diagonal_stats = struct();
    diagonal_stats.means = zeros(F, 1);
    diagonal_stats.stds = zeros(F, 1);
    diagonal_stats.mins = zeros(F, 1);
    diagonal_stats.maxs = zeros(F, 1);
    diagonal_stats.target_errors = zeros(F, 1);
    
    eigenvalue_stats = struct();
    eigenvalue_stats.min_eig = zeros(F, 1);
    eigenvalue_stats.max_eig = zeros(F, 1);
    eigenvalue_stats.condition_nums = zeros(F, 1);
    eigenvalue_stats.num_negative = zeros(F, 1);
    eigenvalue_stats.num_near_zero = zeros(F, 1);
    
    for f = 1:F
        S = Sigma_tilde{f};
        
        % Diagonal analysis
        diag_S = real(diag(S));
        diagonal_stats.means(f) = mean(diag_S);
        diagonal_stats.stds(f) = std(diag_S);
        diagonal_stats.mins(f) = min(diag_S);
        diagonal_stats.maxs(f) = max(diag_S);
        diagonal_stats.target_errors(f) = abs(diagonal_stats.means(f) - 1.0);
        
        % Eigenvalue analysis
        try
            eigs_S = eig(S);
            eigs_real = real(eigs_S);
            
            eigenvalue_stats.min_eig(f) = min(eigs_real);
            eigenvalue_stats.max_eig(f) = max(eigs_real);
            
            positive_eigs = eigs_real(eigs_real > 1e-12);
            if length(positive_eigs) > 1
                eigenvalue_stats.condition_nums(f) = max(positive_eigs) / min(positive_eigs);
            else
                eigenvalue_stats.condition_nums(f) = Inf;
            end
            
            eigenvalue_stats.num_negative(f) = sum(eigs_real < -1e-12);
            eigenvalue_stats.num_near_zero(f) = sum(abs(eigs_real) < 1e-10);
            
        catch ME
            fprintf('  WARNING: Frequency %d eigenvalue computation failed: %s\n', f, ME.message);
            eigenvalue_stats.condition_nums(f) = NaN;
        end
    end
    
    % Print statistical summary
    fprintf('\nDiagonal Statistics Summary:\n');
    fprintf('  Mean value range: [%.4f, %.4f]\n', min(diagonal_stats.means), max(diagonal_stats.means));
    fprintf('  Target error average: %.4f\n', mean(diagonal_stats.target_errors));
    fprintf('  Target error maximum: %.4f\n', max(diagonal_stats.target_errors));
    
    large_error_freqs = find(diagonal_stats.target_errors > 0.1);
    if ~isempty(large_error_freqs)
        fprintf('  ERROR Large error frequencies (>10%%): %s\n', mat2str(large_error_freqs));
    else
        fprintf('  OK All frequency diagonal errors <10%%\n');
    end
    
    fprintf('\nEigenvalue Statistics Summary:\n');
    valid_conds = eigenvalue_stats.condition_nums(~isnan(eigenvalue_stats.condition_nums) & ~isinf(eigenvalue_stats.condition_nums));
    if ~isempty(valid_conds)
        fprintf('  Condition number range: [%.2e, %.2e]\n', min(valid_conds), max(valid_conds));
        fprintf('  Condition number average: %.2e\n', mean(valid_conds));
    end
    
    fprintf('  Total negative eigenvalues: %d\n', sum(eigenvalue_stats.num_negative));
    fprintf('  Total near-zero eigenvalues: %d\n', sum(eigenvalue_stats.num_near_zero));
    
    % Check if truly near-singular
    problematic_freqs = find(eigenvalue_stats.condition_nums > 1e10 | isinf(eigenvalue_stats.condition_nums));
    if ~isempty(problematic_freqs)
        fprintf('  ERROR Truly near-singular frequencies: %s\n', mat2str(problematic_freqs));
    else
        fprintf('  OK No truly near-singular issues\n');
    end
end

function diagnose_conditioning_issues(results)
% Diagnose conditioning and numerical stability issues
    
    fprintf('Diagnosing conditioning and numerical stability issues...\n');
    
    % Compare original and whitened condition numbers
    if isfield(results, 'Sigma_emp') && isfield(results, 'Sigma_tilde')
        Sigma_emp = results.Sigma_emp;
        Sigma_tilde = results.Sigma_tilde;
        F = length(Sigma_emp);
        
        cond_before = zeros(F, 1);
        cond_after = zeros(F, 1);
        
        for f = 1:F
            try
                cond_before(f) = cond(Sigma_emp{f});
                cond_after(f) = cond(Sigma_tilde{f});
            catch
                cond_before(f) = NaN;
                cond_after(f) = NaN;
            end
        end
        
        fprintf('\nCondition Number Comparison:\n');
        fprintf('  Original matrix condition: %.2e (mean), %.2e (max)\n', ...
                nanmean(cond_before), nanmax(cond_before));
        fprintf('  Whitened matrix condition: %.2e (mean), %.2e (max)\n', ...
                nanmean(cond_after), nanmax(cond_after));
        
        improvement = cond_before ./ cond_after;
        valid_improvement = improvement(~isnan(improvement) & ~isinf(improvement));
        
        if ~isempty(valid_improvement)
            if mean(valid_improvement) > 1
                fprintf('  OK Condition number improved by: %.2fx (average)\n', mean(valid_improvement));
            else
                fprintf('  ERROR Condition number degraded by: %.2fx (average)\n', mean(valid_improvement));
            end
        end
    end
    
    % Check whitening matrix D properties
    if isfield(results, 'D') && iscell(results.D)
        fprintf('\nWhitening Matrix D Diagnosis:\n');
        D = results.D;
        
        for f = 1:min(3, length(D))  % Only check first 3
            D_f = D{f};
            fprintf('  Frequency %d - D matrix condition: %.2e\n', f, cond(D_f));
            fprintf('  Frequency %d - D diagonal range: [%.4f, %.4f]\n', f, min(real(diag(D_f))), max(real(diag(D_f))));
        end
    end
end

function verify_quality_assessment(results)
% Verify the quality assessment algorithm itself
    
    fprintf('Verifying quality assessment algorithm...\n');
    
    if ~isfield(results, 'processing_stats') || ...
       ~isfield(results.processing_stats, 'whitening_quality')
        fprintf('ERROR: Quality assessment data does not exist\n');
        return;
    end
    
    wq = results.processing_stats.whitening_quality;
    
    % Check critical fields
    critical_fields = {'whitening_effectiveness', 'diagonal_normalization', 'overall_score'};
    
    for i = 1:length(critical_fields)
        field = critical_fields{i};
        if isfield(wq, field)
            fprintf('  OK %s: exists\n', field);
            
            if strcmp(field, 'whitening_effectiveness')
                effectiveness = wq.(field);
                if isnumeric(effectiveness)
                    fprintf('    Values: [%.6f, %.6f] (range)\n', min(effectiveness), max(effectiveness));
                    fprintf('    Average: %.6f\n', mean(effectiveness));
                    
                    if all(effectiveness == 0)
                        fprintf('    ERROR: All whitening effectiveness values are 0 - this is the problem!\n');
                    end
                end
            end
        else
            fprintf('  ERROR %s: missing\n', field);
        end
    end
    
    % Re-compute a simple quality metric for comparison
    fprintf('\nRe-computing simple quality metric:\n');
    if isfield(results, 'Sigma_tilde')
        Sigma_tilde = results.Sigma_tilde;
        F = length(Sigma_tilde);
        
        simple_quality = zeros(F, 1);
        for f = 1:F
            S = Sigma_tilde{f};
            diag_S = real(diag(S));
            
            % Simple quality = 1 - deviation from target diagonal values
            simple_quality(f) = 1 - mean(abs(diag_S - 1.0));
        end
        
        fprintf('  Simple quality assessment average: %.6f\n', mean(simple_quality));
        fprintf('  Simple quality assessment range: [%.6f, %.6f]\n', min(simple_quality), max(simple_quality));
        
        if mean(simple_quality) > 0.5
            fprintf('  OK Simple assessment shows acceptable quality\n');
        else
            fprintf('  ERROR Simple assessment also shows poor quality\n');
        end
    end
end

function generate_fix_recommendations(results)
% Generate specific fix recommendations
    
    fprintf('Generating fix recommendations...\n\n');
    
    fprintf('IMMEDIATE FIX RECOMMENDATIONS:\n');
    fprintf('==============================\n');
    
    % Give recommendations based on analysis results
    if isfield(results, 'Sigma_tilde')
        Sigma_tilde = results.Sigma_tilde;
        F = length(Sigma_tilde);
        
        % Check diagonal issues
        diagonal_issues = false;
        for f = 1:F
            diag_mean = mean(real(diag(Sigma_tilde{f})));
            if abs(diag_mean - 1.0) > 0.1
                diagonal_issues = true;
                break;
            end
        end
        
        if diagonal_issues
            fprintf('1. Fix diagonal normalization:\n');
            fprintf('   - Check whitening matrix D construction algorithm\n');
            fprintf('   - Verify D = diag(g_smooth)^(-1/2) computation\n');
            fprintf('   - Ensure Sigma_tilde = D * Sigma_emp_loaded * D implementation is correct\n\n');
        end
        
        % Check quality assessment algorithm
        fprintf('2. Fix quality assessment algorithm:\n');
        fprintf('   - Re-implement whitening_effectiveness calculation\n');
        fprintf('   - Check numerical stability of condition number computation\n');
        fprintf('   - Adjust near-singular judgment threshold\n\n');
        
        fprintf('3. Improve numerical stability:\n');
        fprintf('   - Increase diagonal loading strength (current 0.02 -> 0.05)\n');
        fprintf('   - Use more stable matrix decomposition methods\n');
        fprintf('   - Add positive definiteness checks and corrections\n\n');
        
        fprintf('4. Enhance debugging output:\n');
        fprintf('   - Add numerical range checks at key steps\n');
        fprintf('   - Record intermediate computation results\n');
        fprintf('   - Add more detailed error diagnostics\n\n');
    end
    
    fprintf('SUGGESTED CODE MODIFICATIONS:\n');
    fprintf('=============================\n');
    fprintf('1. In module1_preprocessing_main.m add:\n');
    fprintf('   fprintf(''Pre-whitening diagonal range: [%%.3f, %%.3f]\\n'', min(diag(S)), max(diag(S)));\n');
    fprintf('   fprintf(''Post-whitening diagonal range: [%%.3f, %%.3f]\\n'', min(diag(S_tilde)), max(diag(S_tilde)));\n\n');
    
    fprintf('2. Modify quality assessment function to use simpler reliable metrics\n');
    fprintf('3. Adjust parameters: diagonal_loading_factor from 0.02 to 0.05\n\n');
end

function create_diagnostic_plots(results)
% Create diagnostic plots
    
    if ~isfield(results, 'Sigma_tilde') || ~iscell(results.Sigma_tilde)
        fprintf('ERROR: Cannot create diagnostic plots - data not available\n');
        return;
    end
    
    Sigma_tilde = results.Sigma_tilde;
    F = length(Sigma_tilde);
    
    figure('Name', 'Module 1 Debug Diagnostics', 'Position', [100, 100, 1200, 800]);
    
    % Subplot 1: Diagonal distribution
    subplot(2, 3, 1);
    all_diagonals = [];
    for f = 1:F
        all_diagonals = [all_diagonals; real(diag(Sigma_tilde{f}))];
    end
    histogram(all_diagonals, 20);
    xline(1.0, 'r--', 'Target', 'LineWidth', 2);
    xlabel('Diagonal Values');
    ylabel('Frequency');
    title('All Diagonal Elements Distribution');
    grid on;
    
    % Subplot 2: Diagonal mean trend
    subplot(2, 3, 2);
    diagonal_means = zeros(F, 1);
    for f = 1:F
        diagonal_means(f) = mean(real(diag(Sigma_tilde{f})));
    end
    plot(1:F, diagonal_means, 'b-o', 'LineWidth', 2);
    hold on;
    yline(1.0, 'r--', 'Target', 'LineWidth', 2);
    xlabel('Frequency Index');
    ylabel('Diagonal Mean');
    title('Diagonal Mean Across Frequencies');
    grid on;
    
    % Subplot 3: Condition numbers
    subplot(2, 3, 3);
    condition_nums = zeros(F, 1);
    for f = 1:F
        try
            condition_nums(f) = cond(Sigma_tilde{f});
        catch
            condition_nums(f) = NaN;
        end
    end
    semilogy(1:F, condition_nums, 'g-s', 'LineWidth', 2);
    xlabel('Frequency Index');
    ylabel('Condition Number');
    title('Condition Number Evolution');
    grid on;
    
    % Subplot 4: Matrix heatmap example (first frequency)
    subplot(2, 3, 4);
    imagesc(real(Sigma_tilde{1}));
    colorbar;
    title('Frequency 1 - Real Part Heatmap');
    axis square;
    
    % Subplot 5: Eigenvalue distribution
    subplot(2, 3, 5);
    eig_example = real(eig(Sigma_tilde{1}));
    stem(1:length(eig_example), sort(eig_example, 'descend'), 'filled');
    xlabel('Eigenvalue Index');
    ylabel('Eigenvalue');
    title('Frequency 1 - Eigenvalue Distribution');
    grid on;
    
    % Subplot 6: Summary statistics
    subplot(2, 3, 6);
    axis off;
    
    % Calculate key statistics
    diag_error = mean(abs(diagonal_means - 1.0));
    mean_cond = nanmean(condition_nums);
    num_problematic = sum(condition_nums > 1000 | isnan(condition_nums));
    
    % Create status indicators using if-else instead of ternary operator
    if diag_error < 0.1
        diag_status = 'PASS';
    else
        diag_status = 'FAIL';
    end
    
    if mean_cond < 1000
        cond_status = 'PASS';
    else
        cond_status = 'FAIL';
    end
    
    if num_problematic == 0
        stability_status = 'PASS';
    else
        stability_status = 'FAIL';
    end
    
    summary_text = {
        sprintf('Frequencies: %d', F),
        sprintf('Matrix size: %dx%d', size(Sigma_tilde{1})),
        sprintf('Avg diagonal error: %.4f', diag_error),
        sprintf('Avg condition number: %.2e', mean_cond),
        sprintf('Problematic freqs: %d', num_problematic),
        '',
        'Status Assessment:',
        sprintf('Diagonal: %s', diag_status),
        sprintf('Condition: %s', cond_status),
        sprintf('Stability: %s', stability_status)
    };
    
    text(0.1, 0.9, summary_text, 'VerticalAlignment', 'top', 'FontSize', 10);
    title('Diagnostic Summary');
    
    fprintf('OK Diagnostic plots created\n');
end