classdef Module1Visualizer < handle
% MODULE1VISUALIZER - Visualization utilities for Module 1 preprocessing results
%
% This class provides comprehensive visualization capabilities for Module 1
% preprocessing results, including quality assessment, diagnostic plots,
% and comparative analysis.
%
% Usage:
%   viz = Module1Visualizer();
%   viz.plot_preprocessing_overview(preprocessing_results);
%   viz.plot_quality_assessment(preprocessing_results);
%   viz.plot_eeg_analysis(preprocessing_results, frequencies);
%
% File location: src/utils/Module1Visualizer.m
% Action: CREATE new file (utility class)

    properties (Access = private)
        default_colors
        figure_counter
        save_figures
        output_dir
    end
    
    methods
        function obj = Module1Visualizer(varargin)
            % Constructor
            p = inputParser;
            addParameter(p, 'save_figures', false, @islogical);
            addParameter(p, 'output_dir', './figures/', @ischar);
            parse(p, varargin{:});
            
            obj.save_figures = p.Results.save_figures;
            obj.output_dir = p.Results.output_dir;
            obj.figure_counter = 0;
            
            % Initialize default colors
            obj.default_colors = struct();
            obj.default_colors.primary = [0.2, 0.4, 0.8];
            obj.default_colors.secondary = [0.8, 0.4, 0.2];
            obj.default_colors.success = [0.2, 0.8, 0.4];
            obj.default_colors.warning = [0.9, 0.7, 0.1];
            obj.default_colors.error = [0.8, 0.2, 0.2];
            
            if obj.save_figures && ~exist(obj.output_dir, 'dir')
                mkdir(obj.output_dir);
            end
        end
        
        function fig_handles = plot_preprocessing_overview(obj, preprocessing_results)
            % Create comprehensive overview of preprocessing results
            
            fprintf('Creating preprocessing overview visualization...\n');
            
            % Extract data
            data = obj.extract_preprocessing_data(preprocessing_results);
            
            % Create main overview figure
            fig_handles = struct();
            fig_handles.overview = obj.create_figure('Preprocessing Overview', [100, 600, 1400, 700]);
            
            % Subplot 1: Diagonal Power Evolution
            subplot(2, 4, 1);
            obj.plot_diagonal_power_evolution(data);
            
            % Subplot 2: Smoothing Effectiveness
            subplot(2, 4, 2);
            obj.plot_smoothing_effectiveness(data);
            
            % Subplot 3: Whitening Quality
            subplot(2, 4, 3);
            obj.plot_whitening_quality_summary(data);
            
            % Subplot 4: Diagonal Error Distribution
            subplot(2, 4, 4);
            obj.plot_diagonal_error_distribution(data);
            
            % Subplot 5: Effectiveness vs Frequency
            subplot(2, 4, 5);
            obj.plot_effectiveness_vs_frequency(data);
            
            % Subplot 6: Condition Numbers
            subplot(2, 4, 6);
            obj.plot_condition_numbers(data);
            
            % Subplot 7: Success Rate Analysis
            subplot(2, 4, 7);
            obj.plot_success_rate_analysis(data);
            
            % Subplot 8: Processing Time Breakdown
            subplot(2, 4, 8);
            obj.plot_timing_breakdown(data);
            
            if obj.save_figures
                obj.save_figure(fig_handles.overview, 'preprocessing_overview');
            end
            
            fprintf('Created preprocessing overview visualization\n');
        end
        
        function fig_handles = plot_quality_assessment(obj, preprocessing_results)
            % Create detailed quality assessment visualization
            
            fprintf('Creating quality assessment visualization...\n');
            
            data = obj.extract_preprocessing_data(preprocessing_results);
            fig_handles = struct();
            
            % Detailed quality analysis
            fig_handles.quality = obj.create_figure('Quality Assessment', [150, 500, 1200, 800]);
            
            % Diagonal error heatmap
            subplot(2, 3, 1);
            obj.plot_diagonal_error_heatmap(data);
            
            % Quality metrics correlation
            subplot(2, 3, 2);
            obj.plot_quality_metrics_correlation(data);
            
            % Frequency-wise quality trend
            subplot(2, 3, 3);
            obj.plot_frequency_quality_trend(data);
            
            % Channel-wise quality distribution
            subplot(2, 3, 4);
            obj.plot_channel_quality_distribution(data);
            
            % Quality improvement tracking
            subplot(2, 3, 5);
            obj.plot_quality_improvement_tracking(data);
            
            % Overall quality score
            subplot(2, 3, 6);
            obj.plot_overall_quality_score(data);
            
            if obj.save_figures
                obj.save_figure(fig_handles.quality, 'quality_assessment');
            end
            
            fprintf('Created quality assessment visualization\n');
        end
        
        function fig_handles = plot_eeg_analysis(obj, preprocessing_results, frequencies)
            % Create EEG-specific analysis visualization
            
            if nargin < 3
                F = length(preprocessing_results.Sigma_emp);
                frequencies = linspace(1, 40, F); % Default EEG range
            end
            
            fprintf('Creating EEG analysis visualization...\n');
            
            data = obj.extract_preprocessing_data(preprocessing_results);
            data.frequencies = frequencies;
            
            fig_handles = struct();
            
            % EEG power spectrum analysis
            fig_handles.eeg_spectrum = obj.create_figure('EEG Power Spectrum Analysis', [200, 600, 1400, 600]);
            
            subplot(2, 3, 1);
            obj.plot_eeg_power_spectrum(data);
            
            subplot(2, 3, 2);
            obj.plot_eeg_frequency_bands(data);
            
            subplot(2, 3, 3);
            obj.plot_eeg_channel_power_distribution(data);
            
            subplot(2, 3, 4);
            obj.plot_eeg_spatial_correlation(data);
            
            subplot(2, 3, 5);
            obj.plot_eeg_whitening_effect(data);
            
            subplot(2, 3, 6);
            obj.plot_eeg_topography_simulation(data);
            
            if obj.save_figures
                obj.save_figure(fig_handles.eeg_spectrum, 'eeg_analysis');
            end
            
            fprintf('Created EEG analysis visualization\n');
        end
        
        function fig_handles = plot_comparison(obj, results1, results2, labels)
            % Compare two preprocessing results
            
            if nargin < 4
                labels = {'Result 1', 'Result 2'};
            end
            
            fprintf('Creating comparison visualization...\n');
            
            data1 = obj.extract_preprocessing_data(results1);
            data2 = obj.extract_preprocessing_data(results2);
            
            fig_handles = struct();
            fig_handles.comparison = obj.create_figure('Preprocessing Comparison', [250, 500, 1200, 700]);
            
            % Diagonal error comparison
            subplot(2, 3, 1);
            obj.plot_diagonal_error_comparison(data1, data2, labels);
            
            % Effectiveness comparison
            subplot(2, 3, 2);
            obj.plot_effectiveness_comparison(data1, data2, labels);
            
            % Success rate comparison
            subplot(2, 3, 3);
            obj.plot_success_rate_comparison(data1, data2, labels);
            
            % Timing comparison
            subplot(2, 3, 4);
            obj.plot_timing_comparison(data1, data2, labels);
            
            % Quality distribution comparison
            subplot(2, 3, 5);
            obj.plot_quality_distribution_comparison(data1, data2, labels);
            
            % Overall assessment
            subplot(2, 3, 6);
            obj.plot_overall_assessment_comparison(data1, data2, labels);
            
            if obj.save_figures
                obj.save_figure(fig_handles.comparison, 'preprocessing_comparison');
            end
            
            fprintf('Created comparison visualization\n');
        end
        
        function fig_handles = plot_diagnostic(obj, preprocessing_results)
            % Create diagnostic plots for troubleshooting
            
            fprintf('Creating diagnostic visualization...\n');
            
            data = obj.extract_preprocessing_data(preprocessing_results);
            fig_handles = struct();
            
            fig_handles.diagnostic = obj.create_figure('Diagnostic Analysis', [300, 400, 1400, 800]);
            
            % Raw vs processed data
            subplot(2, 4, 1);
            obj.plot_raw_vs_processed_powers(data);
            
            % Smoothing artifacts
            subplot(2, 4, 2);
            obj.plot_smoothing_artifacts(data);
            
            % Whitening matrix properties
            subplot(2, 4, 3);
            obj.plot_whitening_matrix_properties(data);
            
            % Eigenvalue analysis
            subplot(2, 4, 4);
            obj.plot_eigenvalue_analysis(data);
            
            % Outlier detection
            subplot(2, 4, 5);
            obj.plot_outlier_detection(data);
            
            % Convergence analysis
            subplot(2, 4, 6);
            obj.plot_convergence_analysis(data);
            
            % Numerical stability
            subplot(2, 4, 7);
            obj.plot_numerical_stability(data);
            
            % Recommendations
            subplot(2, 4, 8);
            obj.plot_recommendations(data);
            
            if obj.save_figures
                obj.save_figure(fig_handles.diagnostic, 'diagnostic_analysis');
            end
            
            fprintf('Created diagnostic visualization\n');
        end
    end
    
    methods (Access = private)
        function data = extract_preprocessing_data(obj, preprocessing_results)
            % Extract and organize data for visualization
            
            data = struct();
            
            % Basic dimensions
            data.F = length(preprocessing_results.Sigma_emp);
            data.n = size(preprocessing_results.Sigma_emp{1}, 1);
            
            % Extract matrices
            data.Sigma_emp = preprocessing_results.Sigma_emp;
            data.Sigma_tilde = preprocessing_results.Sigma_tilde;
            if isfield(preprocessing_results, 'g_smooth')
                data.g_smooth = preprocessing_results.g_smooth;
            end
            if isfield(preprocessing_results, 'D')
                data.D = preprocessing_results.D;
            end
            
            % Extract quality metrics
            if isfield(preprocessing_results, 'processing_stats') && ...
               isfield(preprocessing_results.processing_stats, 'whitening_quality')
                data.quality = preprocessing_results.processing_stats.whitening_quality;
            end
            
            % Extract timing information
            if isfield(preprocessing_results, 'timing')
                data.timing = preprocessing_results.timing;
            end
            
            % Extract parameters
            if isfield(preprocessing_results, 'parameters')
                data.parameters = preprocessing_results.parameters;
            end
            
            % Compute derived quantities
            data = obj.compute_derived_quantities(data);
        end
        
        function data = compute_derived_quantities(obj, data)
            % Compute derived quantities for visualization
            
            % Diagonal powers across frequencies
            if isfield(data, 'Sigma_emp') && isfield(data, 'Sigma_tilde')
                data.raw_powers = zeros(data.n, data.F);
                data.whitened_powers = zeros(data.n, data.F);
                
                for omega = 1:data.F
                    data.raw_powers(:, omega) = real(diag(data.Sigma_emp{omega}));
                    data.whitened_powers(:, omega) = real(diag(data.Sigma_tilde{omega}));
                end
            end
            
            % Smoothed powers
            if isfield(data, 'g_smooth')
                data.smoothed_powers = zeros(data.n, data.F);
                for omega = 1:data.F
                    data.smoothed_powers(:, omega) = data.g_smooth{omega};
                end
            end
            
            % Quality-derived metrics
            if isfield(data, 'quality')
                if isfield(data.quality, 'diagonal_errors')
                    data.max_diagonal_errors = max(data.quality.diagonal_errors, [], 2);
                    data.mean_diagonal_errors = mean(data.quality.diagonal_errors, 2);
                end
                
                if isfield(data.quality, 'whitening_effectiveness')
                    data.overall_effectiveness = mean(data.quality.whitening_effectiveness);
                end
            end
        end
        
        function fig = create_figure(obj, title_str, position)
            % Create figure with consistent styling
            
            obj.figure_counter = obj.figure_counter + 1;
            
            fig = figure('Name', title_str, 'Position', position, ...
                        'Color', 'white', 'NumberTitle', 'off');
            
            % Set default properties
            set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
            set(groot, 'defaultLegendInterpreter', 'latex');
            set(groot, 'defaultTextInterpreter', 'latex');
        end
        
        function save_figure(obj, fig, filename)
            % Save figure with consistent naming
            
            if obj.save_figures
                full_filename = fullfile(obj.output_dir, ...
                    sprintf('%03d_%s.png', obj.figure_counter, filename));
                saveas(fig, full_filename);
                fprintf('Saved figure: %s\n', full_filename);
            end
        end
        
        % Individual plotting methods would continue here...
        % (Implementation of specific plotting functions)
        
    end
end