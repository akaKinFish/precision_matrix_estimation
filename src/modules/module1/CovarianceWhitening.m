classdef CovarianceWhitening < handle
    % COVARIANCE_WHITENING - Enhanced covariance whitening with quality control
    % Version: 3.1 (Compatibility & PSD-safe) — 2025-08

    properties (Access = private)
        verbose = true;
    end

    methods (Access = public)
        function obj = CovarianceWhitening(varargin)
            p = inputParser;
            addParameter(p, 'verbose', true, @islogical);
            parse(p, varargin{:});
            obj.verbose = p.Results.verbose;
        end

        function [Sigma_tilde, quality_stats] = apply_whitening(obj, Sigma_emp, D, varargin)
            % APPLY_WHITENING — accepts Name-Value:
            %   'target_diagonal' (default 1.0)
            %   'diagonal_tolerance' (default 0.1)
            %   'force_hermitian' (default true)     <-- 新增兼容
            %   'check_psd'       (default true)     <-- 新增兼容
            %   'verbose'         (default obj.verbose)

            p = inputParser;
            addParameter(p, 'target_diagonal', 1.0, @(x) isnumeric(x) && isscalar(x) && x > 0);
            addParameter(p, 'diagonal_tolerance', 0.1, @(x) isnumeric(x) && isscalar(x) && x > 0);
            addParameter(p, 'force_hermitian', true, @islogical);  % 新增
            addParameter(p, 'check_psd', true, @islogical);        % 新增
            addParameter(p, 'verbose', obj.verbose, @islogical);   % 新增
            parse(p, varargin{:});
            params = p.Results;
            obj.verbose = params.verbose;

            obj.validate_whitening_inputs(Sigma_emp, D);
            F = numel(Sigma_emp);
            n = size(Sigma_emp{1},1);

            Sigma_tilde = cell(F,1);
            quality_stats = struct();

            diagonal_errors = cell(F, 1);
            max_diagonal_errors = zeros(F, 1);
            mean_diagonal_errors = zeros(F, 1);
            hermitian_errors = zeros(F, 1);
            min_eigenvalues = zeros(F, 1);
            condition_numbers = zeros(F, 1);
            whitening_effectiveness = zeros(F, 1);

            if obj.verbose
                fprintf('Applying enhanced whitening to %d frequencies, %d nodes\n', F, n);
            end

            for omega = 1:F
                try
                    % 白化（使用 D * Σ * D'，更稳健）
                    Sigma_t = D{omega} * Sigma_emp{omega} * D{omega}';
                    % 对称化（若需要）
                    if params.force_hermitian
                        Sigma_t = (Sigma_t + Sigma_t') / 2;
                    end
                    % PSD 保护（若需要）
                    if params.check_psd
                        % 简单快速的 PSD 保护：若最小特征值 < 0，则做微小对角加载
                        mineig = min(real(eig(Sigma_t)));
                        if ~isfinite(mineig)
                            % 极端数值问题，先强制对称再退避
                            Sigma_t = (Sigma_t + Sigma_t')/2;
                            mineig = min(real(eig(Sigma_t)));
                        end
                        if mineig < -1e-10
                            eps_load = (-mineig) * 1.05 + 1e-12;
                            Sigma_t = Sigma_t + eps_load * eye(n);
                            Sigma_t = (Sigma_t + Sigma_t')/2;
                        end
                    end

                    Sigma_tilde{omega} = Sigma_t;

                    % 质量度量
                    fq = obj.compute_quality_metrics(Sigma_t, params, omega);
                    diagonal_errors{omega} = fq.diagonal_errors;
                    max_diagonal_errors(omega) = fq.max_diagonal_error;
                    mean_diagonal_errors(omega) = fq.mean_diagonal_error;
                    hermitian_errors(omega) = fq.hermitian_error;
                    min_eigenvalues(omega) = fq.min_eigenvalue;
                    condition_numbers(omega) = fq.condition_number;
                    whitening_effectiveness(omega) = fq.whitening_effectiveness;

                    if obj.verbose && mod(omega, max(1, floor(F/5))) == 0
                        fprintf('Processed frequency %d/%d\n', omega, F);
                    end

                catch ME
                    error('covariance_whitening:processing_failed', ...
                          'Failed to process frequency %d: %s', omega, ME.message);
                end
            end

            % 汇总
            quality_stats.diagonal_errors = diagonal_errors;
            quality_stats.max_diagonal_errors = max_diagonal_errors;
            quality_stats.mean_diagonal_errors = mean_diagonal_errors;
            quality_stats.hermitian_errors = hermitian_errors;
            quality_stats.min_eigenvalues = min_eigenvalues;
            quality_stats.condition_numbers = condition_numbers;
            quality_stats.whitening_effectiveness = whitening_effectiveness;

            quality_stats.overall_max_error = max(max_diagonal_errors);
            quality_stats.overall_mean_error = mean(mean_diagonal_errors);
            quality_stats.success_rates = struct();

            tolerances = [0.050, 0.080, 0.100, 0.150, 0.200];
            for i = 1:numel(tolerances)
                tol = tolerances(i);
                success_count = sum(max_diagonal_errors <= tol);
                field_name = sprintf('tol_%03d', round(tol * 1000));
                quality_stats.success_rates.(field_name) = success_count / F;
            end

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

            if obj.verbose
                obj.perform_enhanced_validation(Sigma_tilde, quality_stats, params);
                obj.report_enhanced_whitening_quality(quality_stats, params);
            end
        end
    end

    methods (Access = private)
        function validate_whitening_inputs(~, Sigma_emp, D)
            F_sigma = numel(Sigma_emp);
            F_D = numel(D);
            if F_sigma ~= F_D
                error('covariance_whitening:dimension_mismatch', ...
                    'Number of covariance matrices (%d) does not match whitening matrices (%d)', F_sigma, F_D);
            end
            if F_sigma == 0, error('covariance_whitening:empty_input','Input arrays are empty'); end

            [n_sigma, m_sigma] = size(Sigma_emp{1});
            [n_D, m_D] = size(D{1});
            if n_sigma ~= m_sigma || n_D ~= m_D
                error('covariance_whitening:not_square','Matrices must be square');
            end
            if n_sigma ~= n_D
                error('covariance_whitening:size_mismatch', ...
                    'Matrix size mismatch: covariance [%d x %d], whitening [%d x %d]', n_sigma, m_sigma, n_D, m_D);
            end
            for k = 1:F_sigma
                off_diag = norm(D{k} - diag(diag(D{k})), 'fro');
                if off_diag > 1e-12
                    error('covariance_whitening:not_diagonal','Whitening matrix %d is not diagonal', k);
                end
                if any(diag(D{k}) <= 0)
                    error('covariance_whitening:negative_diagonal','Whitening matrix %d has non-positive diagonal', k);
                end
            end
        end

        function quality = compute_quality_metrics(~, Sigma_tilde, params, ~)
            n = size(Sigma_tilde,1);
            quality = struct();

            diag_elems = diag(Sigma_tilde);
            quality.diagonal_errors = abs(real(diag_elems) - params.target_diagonal);
            quality.max_diagonal_error = max(quality.diagonal_errors);
            quality.mean_diagonal_error = mean(quality.diagonal_errors);
            quality.std_diagonal_error = std(quality.diagonal_errors);

            H = Sigma_tilde - Sigma_tilde';
            denom = max(norm(Sigma_tilde,'fro'), 1e-12);
            quality.hermitian_error = norm(H,'fro')/denom;

            try
                ev = eig(Sigma_tilde);
                re = real(ev);
                quality.min_eigenvalue = min(re);
                quality.max_eigenvalue = max(re);
                pos = re(re > 1e-12);
                if numel(pos) > 1
                    quality.condition_number = max(pos)/min(pos);
                else
                    quality.condition_number = Inf;
                end
                quality.negative_eigenvals = sum(re < -1e-12);
            catch
                quality.min_eigenvalue = NaN;
                quality.max_eigenvalue = NaN;
                quality.condition_number = Inf;
                quality.negative_eigenvals = NaN;
            end

            diag_var = var(real(diag_elems));
            mean_dev = abs(mean(real(diag_elems)) - params.target_diagonal);
            diag_score = exp(-15*diag_var); diag_score = max(0,min(1,diag_score));
            mean_score = exp(-8*mean_dev);  mean_score = max(0,min(1,mean_score));
            herm_score = exp(-50*quality.hermitian_error); herm_score = max(0,min(1,herm_score));
            if quality.min_eigenvalue > 0
                psd_score = 1.0;
            elseif quality.min_eigenvalue > -1e-10
                psd_score = 0.8;
            else
                psd_score = 0.5;
            end
            if isfinite(quality.condition_number) && quality.condition_number < 1e6
                cond_score = exp(-log10(max(quality.condition_number,1))/4);
            else
                cond_score = 0.2;
            end
            cond_score = max(0,min(1,cond_score));

            w = [0.4, 0.3, 0.15, 0.1, 0.05];
            s = [diag_score, mean_score, herm_score, psd_score, cond_score];
            q = 0; for i=1:numel(w), q = q + w(i)*s(i); end
            if ~isfinite(q), q = max(0,1-mean_dev); end
            quality.whitening_effectiveness = max(0,min(1,q));
            quality.component_scores = struct('diagonal',diag_score,'mean',mean_score, ...
                                              'hermitian',herm_score,'psd',psd_score,'condition',cond_score);
        end

        function perform_enhanced_validation(~, Sigma_tilde, quality, params)
            F = numel(Sigma_tilde);
            problematic = 0; severe = 0;
            for k = 1:F
                S = Sigma_tilde{k};
                if any(~isfinite(S(:)))
                    error('covariance_whitening:infinite_entries','Matrix %d contains non-finite entries', k);
                end
                max_err = quality.max_diagonal_errors(k);
                if max_err > params.diagonal_tolerance * 2, problematic = problematic + 1; end
                if max_err > params.diagonal_tolerance * 4, severe = severe + 1; end
            end
            if problematic > 0
                fprintf('Note: %d/%d frequencies have diagonal errors > %.3f\n', problematic, F, params.diagonal_tolerance*2);
            end
            if severe > 0
                fprintf('Warning: %d/%d frequencies have severe diagonal errors > %.3f\n', severe, F, params.diagonal_tolerance*4);
            end
            if any(~isfinite(quality.whitening_effectiveness))
                error('covariance_whitening:invalid_quality','Quality metrics contain non-finite values');
            end
            if any(quality.whitening_effectiveness < 0 | quality.whitening_effectiveness > 1)
                error('covariance_whitening:invalid_effectiveness','Effectiveness scores must be in [0,1]');
            end
            fprintf('Enhanced covariance whitening validation passed\n');
        end

        function report_enhanced_whitening_quality(~, quality, params)
            F = numel(quality.whitening_effectiveness);
            fprintf('\nEnhanced whitening quality assessment:\n');
            fprintf('=====================================\n');
            fprintf('Diagonal normalization (target: %.3f):\n', params.target_diagonal);
            fprintf('  Max error  - Mean: %.4f, Median: %.4f, Range: [%.4f, %.4f]\n', ...
                mean(quality.max_diagonal_errors), median(quality.max_diagonal_errors), ...
                min(quality.max_diagonal_errors), max(quality.max_diagonal_errors));
            fprintf('  Mean error - Mean: %.4f, Median: %.4f\n', ...
                mean(quality.mean_diagonal_errors), median(quality.mean_diagonal_errors));
            fprintf('  Success rates:\n');
            rf = fieldnames(quality.success_rates);
            for i = 1:numel(rf)
                field = rf{i};
                tol = str2double(field(5:end))/1000;
                rate = quality.success_rates.(field);
                cnt = round(rate*F);
                fprintf('    ≤%.3f:  %d/%d (%.1f%%)\n', tol, cnt, F, rate*100);
            end

            ve = quality.whitening_effectiveness;
            mask = isfinite(ve);
            ve = ve(mask);
            fprintf('\nWhitening effectiveness:\n');
            if ~isempty(ve)
                fprintf('  Mean: %.3f, Median: %.3f, Min: %.3f\n', mean(ve), median(ve), min(ve));
                ex = sum(ve > 0.9); gd = sum(ve > 0.8); ac = sum(ve > 0.7);
                fprintf('  Quality distribution:\n');
                fprintf('    Excellent (>0.9): %d/%d (%.1f%%)\n', ex, F, ex/F*100);
                fprintf('    Good (>0.8):      %d/%d (%.1f%%)\n', gd, F, gd/F*100);
                fprintf('    Acceptable (>0.7): %d/%d (%.1f%%)\n', ac, F, ac/F*100);
            else
                fprintf('  ERROR: All effectiveness values are invalid\n');
            end

            fprintf('\nNumerical stability:\n');
            vc = quality.condition_numbers(isfinite(quality.condition_numbers));
            if ~isempty(vc)
                fprintf('  Condition numbers - Mean: %.2e, Median: %.2e, Max: %.2e\n', ...
                    mean(vc), median(vc), max(vc));
            else
                fprintf('  Condition numbers - All invalid\n');
            end
            negc = sum(quality.min_eigenvalues < -1e-12);
            near_sing = sum(quality.condition_numbers > 1e6 | isinf(quality.condition_numbers));
            fprintf('  Negative eigenvalues: %d/%d frequencies\n', negc, F);
            fprintf('  Near-singular: %d/%d frequencies\n', near_sing, F);

            fprintf('\nOverall assessment:\n');
            fprintf('  %s\n', quality.overall_assessment);
            fprintf('  Overall score: %.3f, Good diagonal rate: %.1f%%\n', ...
                quality.mean_effectiveness, quality.good_diagonal_rate*100);
        end
    end

    methods (Static)
        function [Sigma_tilde, quality_stats] = whiten(Sigma_emp, D, varargin)
            % 兼容 Module 1：允许传入 'force_hermitian' / 'check_psd' / 'verbose'
            obj = CovarianceWhitening('verbose', true);
            [Sigma_tilde, quality_stats] = obj.apply_whitening(Sigma_emp, D, varargin{:});
        end
    end
end
