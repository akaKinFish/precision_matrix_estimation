classdef ComplexDataVisualizer < handle
    % COMPLEXDATAVISUALIZER - Enhanced visualization class for complex data
    %
    % This class provides comprehensive visualization capabilities for Module 1
    % preprocessing results that include complex matrix support, Hermitian 
    % property verification, and complex-specific quality assessments.
    %
    % Usage:
    %   visualizer = ComplexDataVisualizer();
    %   visualizer.visualize_results(demo_results);
    %
    % Public Methods:
    %   visualize_results(demo_results) - Create all visualization figures
    %   create_overview_figure(demo_results) - Complex data overview
    %   create_matrix_analysis_figure(demo_results) - Detailed matrix analysis
    %   create_hermitian_verification_figure(demo_results) - Hermitian verification
    %   create_whitening_dashboard_figure(demo_results) - Whitening quality dashboard
    %   create_comparison_figure(demo_results) - Before/after comparison
    %   create_pipeline_figure(demo_results) - Processing pipeline analysis
    %
    % File location: src/visualization/ComplexDataVisualizer.m
    
    properties (Access = private)
        figure_handles
        current_demo_results
        visualization_options
    end
    
    methods (Access = public)
        
        function obj = ComplexDataVisualizer(varargin)
            % Constructor
            %
            % Usage:
            %   visualizer = ComplexDataVisualizer()
            %   visualizer = ComplexDataVisualizer('option1', value1, ...)
            %
            % Optional parameters:
            %   'figure_position' - Default figure position [x, y, width, height]
            %   'color_scheme'    - Color scheme ('default', 'dark', 'colorful')
            %   'save_figures'    - Save figures to files (true/false)
            %   'output_dir'      - Directory for saving figures
            
            % Parse input arguments
            p = inputParser;
            addParameter(p, 'figure_position', [100, 100, 1400, 800], @isnumeric);
            addParameter(p, 'color_scheme', 'default', @ischar);
            addParameter(p, 'save_figures', false, @islogical);
            addParameter(p, 'output_dir', './figures/', @ischar);
            parse(p, varargin{:});
            
            obj.visualization_options = p.Results;
            obj.figure_handles = {};
            obj.current_demo_results = [];
            
            fprintf('ComplexDataVisualizer initialized\n');
        end
        
        function visualize_results(obj, demo_results)
            % Main visualization method - creates all figures
            %
            % Input:
            %   demo_results - Output structure from demo_module1_preprocessing_updated
            %
            % Output:
            %   Creates 6 figure windows with comprehensive analysis
            
            fprintf('\n=== Starting Enhanced Complex Data Visualization ===\n');
            
            % Validate input data
            if ~obj.validate_demo_results(demo_results)
                error('Invalid demo_results structure for complex visualization');
            end
            
            obj.current_demo_results = demo_results;
            
            % Check if preprocessing was successful
            if ~demo_results.preprocessing.success
                fprintf('Preprocessing failed, showing failure analysis\n');
                obj.show_failure_analysis();
                return;
            end
            
            % Create all visualization figures
            try
                obj.create_overview_figure();
                obj.create_matrix_analysis_figure();
                obj.create_hermitian_verification_figure();
                obj.create_whitening_dashboard_figure();
                obj.create_comparison_figure();
                obj.create_pipeline_figure();
                
                fprintf('Enhanced complex visualization complete! Created 6 figure windows.\n');
                
                % Save figures if requested
                if obj.visualization_options.save_figures
                    obj.save_all_figures();
                end
                
            catch ME
                fprintf('Visualization error: %s\n', ME.message);
                rethrow(ME);
            end
        end
        
        function fig_handle = create_overview_figure(obj)
            % Create complex data overview figure
            
            demo_results = obj.current_demo_results;
            
            fig_handle = figure('Name', 'Complex Data Overview', ...
                               'Position', obj.visualization_options.figure_position);
            obj.figure_handles{end+1} = fig_handle;
            
            % Get complex analysis data
            if isfield(demo_results.data_generation, 'complex_analysis')
                complex_analysis = demo_results.data_generation.complex_analysis;
            else
                complex_analysis = struct('matrices_with_complex', 0, 'avg_complex_fraction', 0);
            end
            
            % Subplot 1: Complex data statistics
            subplot(2, 4, 1);
            obj.plot_complex_statistics(complex_analysis);
            
            % Subplot 2: Complex fraction distribution
            subplot(2, 4, 2);
            obj.plot_complex_fraction_distribution();
            
            % Subplot 3: Imaginary component evolution
            subplot(2, 4, 3);
            obj.plot_imaginary_component_evolution();
            
            % Subplot 4: Processing success indicators
            subplot(2, 4, 4);
            obj.plot_success_indicators();
            
            % Subplot 5: Hermitian error evolution
            subplot(2, 4, 5);
            obj.plot_hermitian_error_evolution();
            
            % Subplot 6: Processing timeline
            subplot(2, 4, 6);
            obj.plot_processing_timeline();
            
            % Subplot 7: Quality metrics summary
            subplot(2, 4, 7);
            obj.plot_quality_metrics_summary();
            
            % Subplot 8: Summary statistics
            subplot(2, 4, 8);
            obj.create_summary_text_panel(complex_analysis);
            
            sgtitle('Module 1 Complex Data Processing Overview', 'FontSize', 16, 'FontWeight', 'bold');
        end
        
        function fig_handle = create_matrix_analysis_figure(obj)
            % Create detailed complex matrix analysis figure
            
            demo_results = obj.current_demo_results;
            results = demo_results.preprocessing.results;
            
            fig_handle = figure('Name', 'Complex Matrix Analysis', ...
                               'Position', obj.visualization_options.figure_position);
            obj.figure_handles{end+1} = fig_handle;
            
            % Get original and processed data
            original_data = obj.get_original_data();
            processed_data = results.Sigma_tilde;
            
            if ~isempty(original_data) && ~isempty(processed_data)
                n_freq = min(length(original_data), length(processed_data));
                
                % Select representative frequencies for detailed analysis
                selected_freq = obj.select_representative_frequencies(n_freq);
                
                for i = 1:length(selected_freq)
                    f = selected_freq(i);
                    
                    % Get matrices
                    if iscell(original_data)
                        matrix_orig = original_data{f};
                    else
                        matrix_orig = original_data;
                    end
                    matrix_proc = processed_data{f};
                    
                    % Original matrix magnitude
                    subplot(3, 6, (i-1)*6 + 1);
                    obj.plot_matrix_magnitude(matrix_orig, sprintf('Original |M_%d|', f));
                    
                    % Original matrix phase
                    subplot(3, 6, (i-1)*6 + 2);
                    obj.plot_matrix_phase(matrix_orig, sprintf('Original ∠M_%d', f));
                    
                    % Processed matrix magnitude
                    subplot(3, 6, (i-1)*6 + 3);
                    obj.plot_matrix_magnitude(matrix_proc, sprintf('Processed |M_%d|', f));
                    
                    % Processed matrix phase
                    subplot(3, 6, (i-1)*6 + 4);
                    obj.plot_matrix_phase(matrix_proc, sprintf('Processed ∠M_%d', f));
                    
                    % Magnitude difference
                    subplot(3, 6, (i-1)*6 + 5);
                    obj.plot_magnitude_difference(matrix_orig, matrix_proc);
                    
                    % Complex component analysis
                    subplot(3, 6, (i-1)*6 + 6);
                    obj.plot_complex_component_analysis(matrix_orig, matrix_proc);
                end
            else
                text(0.5, 0.5, 'Matrix data not available for analysis', ...
                     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                     'FontSize', 14, 'FontWeight', 'bold');
            end
            
            sgtitle('Detailed Complex Matrix Analysis', 'FontSize', 16, 'FontWeight', 'bold');
        end
        
        function fig_handle = create_hermitian_verification_figure(obj)
            % Create Hermitian property verification figure
            
            demo_results = obj.current_demo_results;
            results = demo_results.preprocessing.results;
            
            fig_handle = figure('Name', 'Hermitian Property Verification', ...
                               'Position', [150, 150, 1400, 800]);
            obj.figure_handles{end+1} = fig_handle;
            
            if isfield(results, 'Sigma_tilde') && iscell(results.Sigma_tilde)
                [hermitian_errors, diagonal_imaginary, off_diagonal_symmetry] = ...
                    obj.analyze_hermitian_properties(results.Sigma_tilde);
                
                % Plot 1: Hermitian error evolution
                subplot(2, 3, 1);
                obj.plot_hermitian_errors(hermitian_errors);
                
                % Plot 2: Diagonal imaginary components
                subplot(2, 3, 2);
                obj.plot_diagonal_imaginary(diagonal_imaginary);
                
                % Plot 3: Off-diagonal symmetry
                subplot(2, 3, 3);
                obj.plot_off_diagonal_symmetry(off_diagonal_symmetry);
                
                % Plot 4: Hermitian verification heatmap
                subplot(2, 3, 4);
                obj.plot_hermitian_heatmap(hermitian_errors, diagonal_imaginary, off_diagonal_symmetry);
                
                % Plot 5: Overall Hermitian quality score
                subplot(2, 3, 5);
                obj.plot_hermitian_quality(hermitian_errors, diagonal_imaginary);
                
                % Plot 6: Summary statistics
                subplot(2, 3, 6);
                obj.create_hermitian_summary_panel(hermitian_errors, diagonal_imaginary, off_diagonal_symmetry);
                
            else
                text(0.5, 0.5, 'Processed matrix data not available', ...
                     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                     'FontSize', 14, 'FontWeight', 'bold');
            end
            
            sgtitle('Hermitian Property Verification Analysis', 'FontSize', 16, 'FontWeight', 'bold');
        end
        
        function fig_handle = create_whitening_dashboard_figure(obj)
            % Create complex whitening quality dashboard
            
            demo_results = obj.current_demo_results;
            results = demo_results.preprocessing.results;
            
            fig_handle = figure('Name', 'Complex Whitening Quality Dashboard', ...
                               'Position', [200, 200, 1500, 800]);
            obj.figure_handles{end+1} = fig_handle;
            
            if isfield(results, 'processing_stats') && isfield(results.processing_stats, 'whitening_quality')
                quality = results.processing_stats.whitening_quality;
                
                % Plot 1: Diagonal target achievement
                subplot(2, 4, 1);
                obj.plot_diagonal_achievement(quality);
                
                % Plot 2: Condition number analysis
                subplot(2, 4, 2);
                obj.plot_condition_numbers(quality);
                
                % Plot 3: Whitening effectiveness
                subplot(2, 4, 3);
                obj.plot_whitening_effectiveness(quality);
                
                % Plot 4: Complex preservation analysis
                subplot(2, 4, 4);
                obj.plot_complex_preservation(results);
                
                % Plot 5: Processing time breakdown
                subplot(2, 4, 5);
                obj.plot_timing_breakdown(results);
                
                % Plot 6: Quality metrics evolution
                subplot(2, 4, 6);
                obj.plot_quality_evolution(quality);
                
                % Plot 7: Error correlation analysis
                subplot(2, 4, 7);
                obj.plot_error_correlation(quality);
                
                % Plot 8: Overall quality assessment
                subplot(2, 4, 8);
                obj.plot_overall_quality(quality);
                
            else
                text(0.5, 0.5, 'Whitening quality data not available', ...
                     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                     'FontSize', 14, 'FontWeight', 'bold');
            end
            
            sgtitle('Complex Whitening Quality Dashboard', 'FontSize', 16, 'FontWeight', 'bold');
        end
        
        function fig_handle = create_comparison_figure(obj)
            % Create before/after comparison figure
            
            demo_results = obj.current_demo_results;
            results = demo_results.preprocessing.results;
            
            fig_handle = figure('Name', 'Before/After Complex Comparison', ...
                               'Position', [250, 250, 1600, 900]);
            obj.figure_handles{end+1} = fig_handle;
            
            % Get before and after data
            before_data = obj.get_original_data();
            after_data = results.Sigma_tilde;
            
            if ~isempty(before_data) && ~isempty(after_data)
                % Statistical comparison
                before_stats = obj.compute_complex_statistics(before_data);
                after_stats = obj.compute_complex_statistics(after_data);
                
                % Plot 1: Complex fraction comparison
                subplot(2, 4, 1);
                obj.plot_complex_fraction_comparison(before_stats, after_stats);
                
                % Plot 2: Imaginary component magnitude
                subplot(2, 4, 2);
                obj.plot_imaginary_magnitude_comparison(before_stats, after_stats);
                
                % Plot 3: Matrix condition numbers
                subplot(2, 4, 3);
                obj.plot_condition_comparison(before_stats, after_stats);
                
                % Plot 4: Diagonal elements comparison
                subplot(2, 4, 4);
                obj.plot_diagonal_comparison(before_data, after_data);
                
                % Plot 5: Eigenvalue analysis
                subplot(2, 4, 5);
                obj.plot_eigenvalue_comparison(before_data, after_data);
                
                % Plot 6: Spectral properties
                subplot(2, 4, 6);
                obj.plot_spectral_comparison(before_data, after_data);
                
                % Plot 7: Phase distribution
                subplot(2, 4, 7);
                obj.plot_phase_distribution_comparison(before_data, after_data);
                
                % Plot 8: Summary metrics
                subplot(2, 4, 8);
                obj.plot_comparison_summary(before_stats, after_stats);
                
            else
                text(0.5, 0.5, 'Before/after data not available for comparison', ...
                     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                     'FontSize', 14, 'FontWeight', 'bold');
            end
            
            sgtitle('Before/After Complex Data Comparison', 'FontSize', 16, 'FontWeight', 'bold');
        end
        
        function fig_handle = create_pipeline_figure(obj)
            % Create processing pipeline visualization
            
            demo_results = obj.current_demo_results;
            results = demo_results.preprocessing.results;
            
            fig_handle = figure('Name', 'Complex Processing Pipeline', ...
                               'Position', [300, 300, 1600, 900]);
            obj.figure_handles{end+1} = fig_handle;
            
            % Pipeline flow diagram
            subplot(3, 4, [1, 2]);
            obj.create_pipeline_flowchart();
            
            % Step-by-step analysis
            if isfield(results, 'Sigma_emp') && isfield(results, 'Sigma_tilde')
                
                % Step 1: Input data analysis
                subplot(3, 4, 3);
                obj.plot_pipeline_step(results.Sigma_emp, 'Input Data');
                
                % Step 2: Smoothed data analysis
                subplot(3, 4, 4);
                if isfield(results, 'g_smooth')
                    obj.plot_smoothing_step(results.g_smooth);
                else
                    obj.plot_unavailable_step('Smoothing Data');
                end
                
                % Step 3: Whitening matrix analysis
                subplot(3, 4, 7);
                if isfield(results, 'D')
                    obj.plot_pipeline_step(results.D, 'Whitening Matrices');
                else
                    obj.plot_unavailable_step('Whitening Matrices');
                end
                
                % Step 4: Final output analysis
                subplot(3, 4, 8);
                obj.plot_pipeline_step(results.Sigma_tilde, 'Output Data');
                
                % Processing quality evolution
                subplot(3, 4, [9, 10]);
                obj.plot_processing_quality_evolution();
                
                % Error propagation analysis
                subplot(3, 4, [11, 12]);
                obj.plot_error_propagation();
                
            else
                text(0.5, 0.5, 'Processing pipeline data not available', ...
                     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                     'FontSize', 14, 'FontWeight', 'bold');
            end
            
            sgtitle('Complex Data Processing Pipeline Analysis', 'FontSize', 16, 'FontWeight', 'bold');
        end
        
    end
    
    methods (Access = private)
        
        function is_valid = validate_demo_results(obj, demo_results)
            % Validate demo_results structure
            is_valid = false;
            
            required_fields = {'timestamp', 'data_generation', 'preprocessing'};
            for i = 1:length(required_fields)
                if ~isfield(demo_results, required_fields{i})
                    fprintf('Missing required field: %s\n', required_fields{i});
                    return;
                end
            end
            
            % Check for complex-specific fields
            if isfield(demo_results, 'data_generation') && demo_results.data_generation.success
                if ~isfield(demo_results.data_generation, 'complex_analysis')
                    fprintf('Warning: Missing complex_analysis in data_generation\n');
                end
            end
            
            is_valid = true;
        end
        
        function show_failure_analysis(obj)
            % Show analysis when processing fails
            demo_results = obj.current_demo_results;
            
            fig_handle = figure('Name', 'Complex Processing Failure Analysis', ...
                               'Position', [400, 400, 1200, 600]);
            obj.figure_handles{end+1} = fig_handle;
            
            subplot(1, 2, 1);
            axis off;
            
            failure_text = {'PROCESSING FAILURE ANALYSIS', '', ''};
            
            if isfield(demo_results.preprocessing, 'error')
                error_msg = demo_results.preprocessing.error;
                failure_text{end+1} = 'Error Message:';
                failure_text{end+1} = error_msg;
                failure_text{end+1} = '';
                
                % Analyze error type and provide recommendations
                obj.add_error_recommendations(failure_text, error_msg);
            end
            
            for i = 1:length(failure_text)
                if i <= 3
                    text(0.05, 0.95 - (i-1)*0.08, failure_text{i}, 'FontSize', 14, ...
                         'FontWeight', 'bold', 'VerticalAlignment', 'top');
                else
                    text(0.05, 0.95 - (i-1)*0.08, failure_text{i}, 'FontSize', 10, ...
                         'VerticalAlignment', 'top');
                end
            end
            
            subplot(1, 2, 2);
            if isfield(demo_results, 'data_generation') && demo_results.data_generation.success
                obj.show_successful_data_generation(demo_results);
            else
                text(0.5, 0.5, 'No data generated', 'HorizontalAlignment', 'center');
            end
        end
        
        function add_error_recommendations(obj, failure_text, error_msg)
            % Add specific recommendations based on error message
            if contains(error_msg, 'complex')
                failure_text{end+1} = 'Issue Type: Complex data handling';
                failure_text{end+1} = 'Recommendations:';
                failure_text{end+1} = '• Check CovarianceWhitening class complex support';
                failure_text{end+1} = '• Verify Hermitian matrix operations';
                failure_text{end+1} = '• Update diagonal smoothing for complex data';
            elseif contains(error_msg, 'dimension')
                failure_text{end+1} = 'Issue Type: Dimension mismatch';
                failure_text{end+1} = 'Recommendations:';
                failure_text{end+1} = '• Check matrix size consistency';
                failure_text{end+1} = '• Verify frequency dimension alignment';
            else
                failure_text{end+1} = 'Issue Type: General processing error';
                failure_text{end+1} = 'Recommendations:';
                failure_text{end+1} = '• Review preprocessing parameters';
                failure_text{end+1} = '• Check input data validity';
            end
        end
        
        function show_successful_data_generation(obj, demo_results)
            % Show successfully generated data information
            complex_analysis = demo_results.data_generation.complex_analysis;
            
            bar([complex_analysis.matrices_with_complex, ...
                 complex_analysis.n_frequencies - complex_analysis.matrices_with_complex]);
            set(gca, 'XTickLabel', {'Complex', 'Real'});
            ylabel('Number of Matrices');
            title('Successfully Generated Data');
            grid on;
        end
        
        function original_data = get_original_data(obj)
            % Get original data from demo results
            demo_results = obj.current_demo_results;
            
            if isfield(demo_results.data_generation, 'params') && ...
               isfield(demo_results.data_generation.params, 'Sigma_emp')
                original_data = demo_results.data_generation.params.Sigma_emp;
            else
                original_data = {};
            end
        end
        
        function selected_freq = select_representative_frequencies(obj, n_freq)
            % Select representative frequencies for analysis
            if n_freq >= 3
                selected_freq = [1, round(n_freq/2), n_freq];
            elseif n_freq == 2
                selected_freq = [1, 2];
            else
                selected_freq = 1;
            end
        end
        
        function save_all_figures(obj)
            % Save all figures to files
            output_dir = obj.visualization_options.output_dir;
            
            if ~exist(output_dir, 'dir')
                mkdir(output_dir);
            end
            
            for i = 1:length(obj.figure_handles)
                fig = obj.figure_handles{i};
                filename = sprintf('complex_visualization_%d.png', i);
                filepath = fullfile(output_dir, filename);
                saveas(fig, filepath);
            end
            
            fprintf('Saved %d figures to %s\n', length(obj.figure_handles), output_dir);
        end
        
        % Placeholder methods for specific plotting functions
        % (These would be implemented with actual plotting code)
        
        function plot_complex_statistics(obj, complex_analysis)
            if complex_analysis.matrices_with_complex > 0
                stats_data = [complex_analysis.matrices_with_complex, ...
                             complex_analysis.n_frequencies - complex_analysis.matrices_with_complex];
                pie(stats_data, {'Complex Matrices', 'Real Matrices'});
                title('Matrix Type Distribution');
                colormap([0.2 0.6 1; 0.8 0.8 0.8]);
            else
                text(0.5, 0.5, 'No Complex Matrices', 'HorizontalAlignment', 'center', ...
                     'VerticalAlignment', 'middle', 'FontSize', 12, 'FontWeight', 'bold');
                title('Matrix Type Distribution');
            end
        end
        
        function plot_complex_fraction_distribution(obj)
            text(0.5, 0.5, 'Complex Fraction\nDistribution', 'HorizontalAlignment', 'center');
            title('Complex Fraction Distribution');
        end
        
        function plot_imaginary_component_evolution(obj)
            text(0.5, 0.5, 'Imaginary Component\nEvolution', 'HorizontalAlignment', 'center');
            title('Imaginary Component Evolution');
        end
        
        function plot_success_indicators(obj)
            text(0.5, 0.5, 'Success\nIndicators', 'HorizontalAlignment', 'center');
            title('Processing Success');
        end
        
        function plot_hermitian_error_evolution(obj)
            text(0.5, 0.5, 'Hermitian Error\nEvolution', 'HorizontalAlignment', 'center');
            title('Hermitian Error Evolution');
        end
        
        function plot_processing_timeline(obj)
            text(0.5, 0.5, 'Processing\nTimeline', 'HorizontalAlignment', 'center');
            title('Processing Timeline');
        end
        
        function plot_quality_metrics_summary(obj)
            text(0.5, 0.5, 'Quality Metrics\nSummary', 'HorizontalAlignment', 'center');
            title('Quality Metrics');
        end
        
        function create_summary_text_panel(obj, complex_analysis)
            axis off;
            summary_text = {
                'COMPLEX DATA SUMMARY',
                '',
                sprintf('Total Frequencies: %d', complex_analysis.n_frequencies),
                sprintf('Complex Matrices: %d', complex_analysis.matrices_with_complex),
                sprintf('Avg Complex Fraction: %.3f', complex_analysis.avg_complex_fraction),
                sprintf('Max Imaginary: %.4f', complex_analysis.max_imag_component),
                sprintf('All Hermitian: %s', obj.logical_to_string(complex_analysis.all_hermitian))
            };
            
            for i = 1:length(summary_text)
                if i == 1
                    text(0.05, 0.95 - (i-1)*0.12, summary_text{i}, 'FontSize', 12, ...
                         'FontWeight', 'bold', 'VerticalAlignment', 'top');
                else
                    text(0.05, 0.95 - (i-1)*0.12, summary_text{i}, 'FontSize', 10, ...
                         'VerticalAlignment', 'top');
                end
            end
        end
        
        % Additional plotting method placeholders
        function plot_matrix_magnitude(obj, matrix, title_str)
            imagesc(abs(matrix));
            colorbar;
            title(title_str);
        end
        
        function plot_matrix_phase(obj, matrix, title_str)
            imagesc(angle(matrix));
            colorbar;
            title(title_str);
        end
        
        function plot_magnitude_difference(obj, matrix_orig, matrix_proc)
            mag_diff = abs(matrix_proc) - abs(matrix_orig);
            imagesc(mag_diff);
            colorbar;
            title('Magnitude Difference');
        end
        
        function plot_complex_component_analysis(obj, matrix_orig, matrix_proc)
            imag_orig = mean(abs(imag(matrix_orig(:))));
            imag_proc = mean(abs(imag(matrix_proc(:))));
            
            bar([imag_orig, imag_proc]);
            set(gca, 'XTickLabel', {'Original', 'Processed'});
            ylabel('Mean |Imaginary|');
            title('Complex Component');
            grid on;
        end
        
        function [hermitian_errors, diagonal_imaginary, off_diagonal_symmetry] = analyze_hermitian_properties(obj, matrix_cell_array)
            n_freq = length(matrix_cell_array);
            hermitian_errors = zeros(n_freq, 1);
            diagonal_imaginary = zeros(n_freq, 1);
            off_diagonal_symmetry = zeros(n_freq, 1);
            
            for f = 1:n_freq
                matrix = matrix_cell_array{f};
                hermitian_errors(f) = max(abs(matrix - matrix'));
                diagonal_imaginary(f) = max(abs(imag(diag(matrix))));
                upper_tri = triu(matrix, 1);
                lower_tri = tril(matrix, -1)';
                off_diagonal_symmetry(f) = max(abs(upper_tri(:) - conj(lower_tri(:))));
            end
        end
        
        function plot_hermitian_errors(obj, hermitian_errors)
            semilogy(1:length(hermitian_errors), hermitian_errors, 'b.-', 'LineWidth', 2);
            xlabel('Frequency Index');
            ylabel('Hermitian Error (log)');
            title('Hermitian Error');
            grid on;
        end
        
        function plot_diagonal_imaginary(obj, diagonal_imaginary)
            semilogy(1:length(diagonal_imaginary), diagonal_imaginary, 'r.-', 'LineWidth', 2);
            xlabel('Frequency Index');
            ylabel('Max Diagonal Imaginary (log)');
            title('Diagonal Imaginary');
            grid on;
        end
        
        function plot_off_diagonal_symmetry(obj, off_diagonal_symmetry)
            semilogy(1:length(off_diagonal_symmetry), off_diagonal_symmetry, 'g.-', 'LineWidth', 2);
            xlabel('Frequency Index');
            ylabel('Off-diagonal Asymmetry (log)');
            title('Off-diagonal Symmetry');
            grid on;
        end
        
        function plot_hermitian_heatmap(obj, hermitian_errors, diagonal_imaginary, off_diagonal_symmetry)
            verification_matrix = [hermitian_errors, diagonal_imaginary, off_diagonal_symmetry];
            imagesc(log10(verification_matrix'));
            colorbar;
            xlabel('Frequency Index');
            ylabel('Error Type');
            set(gca, 'YTick', 1:3, 'YTickLabel', {'Hermitian', 'Diagonal Imag', 'Off-diag Sym'});
            title('Hermitian Verification Heatmap (log10)');
        end
        
        function plot_hermitian_quality(obj, hermitian_errors, diagonal_imaginary)
            hermitian_quality = 100 * (hermitian_errors < 1e-10 & diagonal_imaginary < 1e-12);
            bar(1:length(hermitian_quality), hermitian_quality);
            xlabel('Frequency Index');
            ylabel('Hermitian Quality (%)');
            title('Hermitian Quality Score');
            ylim([0, 110]);
            grid on;
        end
        
        function create_hermitian_summary_panel(obj, hermitian_errors, diagonal_imaginary, off_diagonal_symmetry)
            n_freq = length(hermitian_errors);
            
            stats_text = {
                sprintf('Total Frequencies: %d', n_freq),
                sprintf('Perfect Hermitian: %d (%.1f%%)', sum(hermitian_errors < 1e-10), ...
                        100*sum(hermitian_errors < 1e-10)/n_freq),
                sprintf('Max Hermitian Error: %.2e', max(hermitian_errors)),
                sprintf('Max Diagonal Imag: %.2e', max(diagonal_imaginary)),
                sprintf('Max Off-diag Error: %.2e', max(off_diagonal_symmetry)),
                '',
                'Quality Assessment:',
                sprintf('Hermitian: %s', obj.pass_fail_string(all(hermitian_errors < 1e-10))),
                sprintf('Real Diagonal: %s', obj.pass_fail_string(all(diagonal_imaginary < 1e-12))),
                sprintf('Conjugate Sym: %s', obj.pass_fail_string(all(off_diagonal_symmetry < 1e-10)))
            };
            
            axis off;
            for i = 1:length(stats_text)
                text(0.05, 0.95 - (i-1)*0.08, stats_text{i}, 'FontSize', 10, ...
                     'VerticalAlignment', 'top', 'FontWeight', 'bold');
            end
        end
        
        function plot_diagonal_achievement(obj, quality)
            if isfield(quality, 'diagonal_errors')
                histogram(quality.diagonal_errors, 20);
                xlabel('Diagonal Error');
                ylabel('Frequency Count');
                title('Diagonal Target Achievement');
                grid on;
                
                hold on;
                target_line = 0.1;
                plot([target_line, target_line], ylim, 'r--', 'LineWidth', 2);
                legend('Error Distribution', 'Target Tolerance', 'Location', 'best');
                hold off;
            else
                text(0.5, 0.5, 'Diagonal error\ndata not available', 'HorizontalAlignment', 'center');
            end
        end
        
        function plot_condition_numbers(obj, quality)
            if isfield(quality, 'condition_numbers')
                semilogy(1:length(quality.condition_numbers), quality.condition_numbers, 'b.-', ...
                         'LineWidth', 2, 'MarkerSize', 8);
                xlabel('Frequency Index');
                ylabel('Condition Number (log scale)');
                title('Matrix Conditioning');
                grid on;
                
                hold on;
                plot([1, length(quality.condition_numbers)], [1e12, 1e12], 'r--', 'LineWidth', 2);
                legend('Condition Numbers', 'Ill-conditioning Threshold', 'Location', 'best');
                hold off;
            else
                text(0.5, 0.5, 'Condition number\ndata not available', 'HorizontalAlignment', 'center');
            end
        end
        
        function plot_whitening_effectiveness(obj, quality)
            if isfield(quality, 'whitening_effectiveness')
                bar(1:length(quality.whitening_effectiveness), quality.whitening_effectiveness);
                xlabel('Frequency Index');
                ylabel('Effectiveness Score');
                title('Whitening Effectiveness');
                ylim([0, 1.1]);
                grid on;
            else
                text(0.5, 0.5, 'Effectiveness data\nnot available', 'HorizontalAlignment', 'center');
            end
        end
        
        function plot_complex_preservation(obj, results)
            text(0.5, 0.5, 'Complex\nPreservation\nAnalysis', 'HorizontalAlignment', 'center');
            title('Complex Preservation');
        end
        
        function plot_timing_breakdown(obj, results)
            if isfield(results, 'timing')
                timing = results.timing;
                timing_fields = fieldnames(timing);
                timing_values = [];
                timing_labels = {};
                
                for i = 1:length(timing_fields)
                    if ~strcmp(timing_fields{i}, 'total')
                        timing_values(end+1) = timing.(timing_fields{i});
                        timing_labels{end+1} = strrep(timing_fields{i}, '_', ' ');
                    end
                end
                
                if ~isempty(timing_values)
                    pie(timing_values, timing_labels);
                    title('Processing Time Breakdown');
                end
            else
                text(0.5, 0.5, 'Timing data\nnot available', 'HorizontalAlignment', 'center');
            end
        end
        
        function plot_quality_evolution(obj, quality)
            text(0.5, 0.5, 'Quality\nEvolution', 'HorizontalAlignment', 'center');
            title('Quality Evolution');
        end
        
        function plot_error_correlation(obj, quality)
            text(0.5, 0.5, 'Error\nCorrelation', 'HorizontalAlignment', 'center');
            title('Error Correlation');
        end
        
        function plot_overall_quality(obj, quality)
            text(0.5, 0.5, 'Overall\nQuality', 'HorizontalAlignment', 'center');
            title('Overall Quality');
        end
        
        function stats = compute_complex_statistics(obj, matrix_cell_array)
            n_freq = length(matrix_cell_array);
            complex_fractions = zeros(n_freq, 1);
            imag_magnitudes = zeros(n_freq, 1);
            condition_numbers = zeros(n_freq, 1);
            
            for f = 1:n_freq
                if iscell(matrix_cell_array)
                    matrix = matrix_cell_array{f};
                else
                    matrix = matrix_cell_array;
                    break;
                end
                
                complex_elements = sum(abs(imag(matrix(:))) > 1e-12);
                complex_fractions(f) = complex_elements / numel(matrix);
                imag_magnitudes(f) = mean(abs(imag(matrix(:))));
                condition_numbers(f) = cond(matrix);
            end
            
            stats = struct();
            stats.complex_fractions = complex_fractions;
            stats.avg_complex_fraction = mean(complex_fractions);
            stats.imag_magnitudes = imag_magnitudes;
            stats.avg_imag_magnitude = mean(imag_magnitudes);
            stats.condition_numbers = condition_numbers;
            stats.avg_condition_number = mean(condition_numbers);
        end
        
        function plot_complex_fraction_comparison(obj, before_stats, after_stats)
            bar([before_stats.avg_complex_fraction, after_stats.avg_complex_fraction]);
            set(gca, 'XTickLabel', {'Before', 'After'});
            ylabel('Average Complex Fraction');
            title('Complex Element Fraction');
            grid on;
        end
        
        function plot_imaginary_magnitude_comparison(obj, before_stats, after_stats)
            bar([before_stats.avg_imag_magnitude, after_stats.avg_imag_magnitude]);
            set(gca, 'XTickLabel', {'Before', 'After'});
            ylabel('Average |Imaginary|');
            title('Imaginary Component Magnitude');
            grid on;
        end
        
        function plot_condition_comparison(obj, before_stats, after_stats)
            semilogy(1:length(before_stats.condition_numbers), before_stats.condition_numbers, 'b.-', ...
                     1:length(after_stats.condition_numbers), after_stats.condition_numbers, 'r.-', ...
                     'LineWidth', 2, 'MarkerSize', 6);
            xlabel('Frequency Index');
            ylabel('Condition Number (log scale)');
            title('Condition Number Comparison');
            legend('Before', 'After', 'Location', 'best');
            grid on;
        end
        
        function plot_diagonal_comparison(obj, before_data, after_data)
            text(0.5, 0.5, 'Diagonal\nComparison', 'HorizontalAlignment', 'center');
            title('Diagonal Comparison');
        end
        
        function plot_eigenvalue_comparison(obj, before_data, after_data)
            text(0.5, 0.5, 'Eigenvalue\nComparison', 'HorizontalAlignment', 'center');
            title('Eigenvalue Comparison');
        end
        
        function plot_spectral_comparison(obj, before_data, after_data)
            text(0.5, 0.5, 'Spectral\nProperties', 'HorizontalAlignment', 'center');
            title('Spectral Properties');
        end
        
        function plot_phase_distribution_comparison(obj, before_data, after_data)
            text(0.5, 0.5, 'Phase\nDistribution', 'HorizontalAlignment', 'center');
            title('Phase Distribution');
        end
        
        function plot_comparison_summary(obj, before_stats, after_stats)
            text(0.5, 0.5, 'Summary\nMetrics', 'HorizontalAlignment', 'center');
            title('Summary Metrics');
        end
        
        function create_pipeline_flowchart(obj)
            text(0.5, 0.5, 'Processing\nPipeline\nFlowchart', 'HorizontalAlignment', 'center');
            title('Processing Pipeline');
        end
        
        function plot_pipeline_step(obj, data, step_name)
            text(0.5, 0.5, step_name, 'HorizontalAlignment', 'center');
            title(step_name);
        end
        
        function plot_smoothing_step(obj, data)
            text(0.5, 0.5, 'Smoothing\nAnalysis', 'HorizontalAlignment', 'center');
            title('Smoothed Data');
        end
        
        function plot_unavailable_step(obj, step_name)
            text(0.5, 0.5, sprintf('%s\nnot available', step_name), 'HorizontalAlignment', 'center');
            title(step_name);
        end
        
        function plot_processing_quality_evolution(obj)
            text(0.5, 0.5, 'Quality\nEvolution', 'HorizontalAlignment', 'center');
            title('Processing Quality Evolution');
        end
        
        function plot_error_propagation(obj)
            text(0.5, 0.5, 'Error\nPropagation', 'HorizontalAlignment', 'center');
            title('Error Propagation Analysis');
        end
        
        function str = logical_to_string(obj, logical_value)
            if logical_value
                str = 'YES';
            else
                str = 'NO';
            end
        end
        
        function str = pass_fail_string(obj, condition)
            if condition
                str = 'PASS';
            else
                str = 'FAIL';
            end
        end
        
    end
    
    methods (Static)
        
        function demo()
            % Static demo method to test the visualizer
            %
            % Usage:
            %   ComplexDataVisualizer.demo()
            
            fprintf('Running ComplexDataVisualizer demo...\n');
            
            try
                % Run the updated demo
                demo_results = demo_module1_preprocessing_updated();
                
                % Create visualizer and run visualization
                visualizer = ComplexDataVisualizer('save_figures', false);
                visualizer.visualize_results(demo_results);
                
                fprintf('Demo completed successfully!\n');
                
            catch ME
                fprintf('Demo failed: %s\n', ME.message);
                rethrow(ME);
            end
        end
        
    end
    
end