classdef Module1Utils < handle
% MODULE1UTILS - Utility functions for Module 1 preprocessing
%
% This class contains static methods for supporting module1 preprocessing
% operations including data validation, error handling, and helper functions.
%
% Usage:
%   data = Module1Utils.validate_data_structure(input_data);
%   metrics = Module1Utils.compute_quality_metrics(results);
%
% File location: src/modules/module1/Module1Utils.m

    methods (Static)
        
        function is_valid = validate_data_structure(input_data)
            % VALIDATE_DATA_STRUCTURE - Validate input data structure
            %
            % Input:
            %   input_data - Input data structure
            %
            % Output:
            %   is_valid - Boolean indicating if structure is valid
            
            is_valid = false;
            
            try
                % Check required fields
                if ~isfield(input_data, 'mode')
                    warning('Missing required field: mode');
                    return;
                end
                
                switch input_data.mode
                    case 'simulation'
                        if ~isfield(input_data, 'sim_results')
                            warning('Missing sim_results field for simulation mode');
                            return;
                        end
                        
                        sim_results = input_data.sim_results;
                        required_fields = {'Sigma_emp', 'F', 'n', 'T'};
                        
                        for i = 1:length(required_fields)
                            if ~isfield(sim_results, required_fields{i})
                                warning('Missing field in sim_results: %s', required_fields{i});
                                return;
                            end
                        end
                        
                        % Validate Sigma_emp is cell array
                        if ~iscell(sim_results.Sigma_emp)
                            warning('Sigma_emp must be a cell array');
                            return;
                        end
                        
                        % Check dimensions
                        if length(sim_results.Sigma_emp) ~= sim_results.F
                            warning('Sigma_emp length does not match F');
                            return;
                        end
                        
                    case 'eeg_data'
                        if ~isfield(input_data, 'data_path')
                            warning('Missing data_path field for eeg_data mode');
                            return;
                        end
                        
                        if ~isfield(input_data, 'params')
                            warning('Missing params field for eeg_data mode');
                            return;
                        end
                        
                    otherwise
                        warning('Unknown mode: %s', input_data.mode);
                        return;
                end
                
                is_valid = true;
                
            catch ME
                warning(ME.identifier, '%s', ME.message);
                is_valid = false;
            end
        end
        
        function metrics = compute_quality_metrics(preprocessing_results)
            % COMPUTE_QUALITY_METRICS - Compute preprocessing quality metrics
            %
            % Input:
            %   preprocessing_results - Results from preprocessing
            %
            % Output:
            %   metrics - Structure containing quality metrics
            
            metrics = struct();
            
            try
                % Check if required fields exist
                if ~isfield(preprocessing_results, 'Sigma_tilde') || ...
                   ~iscell(preprocessing_results.Sigma_tilde)
                    warning('Missing or invalid Sigma_tilde field');
                    metrics.success = false;
                    return;
                end
                
                Sigma_tilde = preprocessing_results.Sigma_tilde;
                F = length(Sigma_tilde);
                n = size(Sigma_tilde{1}, 1);
                
                % Initialize metric arrays
                diagonal_means = zeros(F, 1);
                condition_numbers = zeros(F, 1);
                frobenius_norms = zeros(F, 1);
                
                % Compute metrics for each frequency
                for f = 1:F
                    S = Sigma_tilde{f};
                    
                    % Diagonal statistics
                    diagonal_means(f) = mean(real(diag(S)));
                    
                    % Condition number (with safety check)
                    try
                        eigs_S = eig(S);
                        eigs_S = real(eigs_S(eigs_S > 1e-12));
                        if length(eigs_S) > 1
                            condition_numbers(f) = max(eigs_S) / min(eigs_S);
                        else
                            condition_numbers(f) = 1;
                        end
                    catch
                        condition_numbers(f) = NaN;
                    end
                    
                    % Frobenius norm
                    frobenius_norms(f) = norm(S, 'fro');
                end
                
                % Summary statistics
                metrics.diagonal_stats = struct();
                metrics.diagonal_stats.mean = mean(diagonal_means);
                metrics.diagonal_stats.std = std(diagonal_means);
                metrics.diagonal_stats.target_deviation = abs(mean(diagonal_means) - 1.0);
                
                metrics.conditioning = struct();
                valid_cond = condition_numbers(~isnan(condition_numbers));
                if ~isempty(valid_cond)
                    metrics.conditioning.mean = mean(valid_cond);
                    metrics.conditioning.max = max(valid_cond);
                    metrics.conditioning.min = min(valid_cond);
                else
                    metrics.conditioning.mean = NaN;
                    metrics.conditioning.max = NaN;
                    metrics.conditioning.min = NaN;
                end
                
                metrics.norms = struct();
                metrics.norms.mean = mean(frobenius_norms);
                metrics.norms.std = std(frobenius_norms);
                
                % Overall quality assessment
                quality_score = 1.0;
                
                % Penalize deviation from target diagonal
                if metrics.diagonal_stats.target_deviation > 0.2
                    quality_score = quality_score * 0.8;
                end
                
                % Penalize high condition numbers
                if ~isnan(metrics.conditioning.mean) && metrics.conditioning.mean > 100
                    quality_score = quality_score * 0.7;
                end
                
                metrics.overall_quality = quality_score;
                metrics.success = true;
                
            catch ME
                warning(ME.identifier, '%s', ME.message);
                metrics.success = false;
                metrics.error = ME.message;
            end
        end
        
        function result = safe_field_access(structure, field_path, default_value)
            % SAFE_FIELD_ACCESS - Safely access nested structure fields
            %
            % Inputs:
            %   structure - Input structure
            %   field_path - String or cell array of field names
            %   default_value - Value to return if field doesn't exist
            %
            % Output:
            %   result - Field value or default_value
            
            if nargin < 3
                default_value = [];
            end
            
            result = default_value;
            
            try
                if ischar(field_path)
                    field_path = {field_path};
                end
                
                current = structure;
                for i = 1:length(field_path)
                    if ~isstruct(current) || ~isfield(current, field_path{i})
                        return;
                    end
                    current = current.(field_path{i});
                end
                
                result = current;
                
            catch
                % Return default value on any error
                result = default_value;
            end
        end
        
        function summary_text = create_processing_summary(results)
            % CREATE_PROCESSING_SUMMARY - Create human-readable summary
            %
            % Input:
            %   results - Processing results structure
            %
            % Output:
            %   summary_text - Multi-line string summary
            
            summary_lines = {};
            
            try
                % Basic information
                if Module1Utils.safe_field_access(results, 'success', false)
                    summary_lines{end+1} = 'Processing Status: SUCCESS';
                else
                    summary_lines{end+1} = 'Processing Status: FAILED';
                    error_msg = Module1Utils.safe_field_access(results, 'error', 'Unknown error');
                    summary_lines{end+1} = sprintf('Error: %s', error_msg);
                    summary_text = strjoin(summary_lines, '\n');
                    return;
                end
                
                % Dimensions
                F = Module1Utils.safe_field_access(results, 'F', 'Unknown');
                n = Module1Utils.safe_field_access(results, 'n', 'Unknown');
                summary_lines{end+1} = sprintf('Dimensions: %s frequencies, %s nodes', ...
                                              num2str(F), num2str(n));
                
                % Quality metrics
                quality_metrics = Module1Utils.safe_field_access(results, 'quality_metrics', struct());
                if isstruct(quality_metrics) && Module1Utils.safe_field_access(quality_metrics, 'success', false)
                    overall_quality = Module1Utils.safe_field_access(quality_metrics, 'overall_quality', NaN);
                    if ~isnan(overall_quality)
                        summary_lines{end+1} = sprintf('Overall Quality: %.3f', overall_quality);
                    end
                    
                    diagonal_mean = Module1Utils.safe_field_access(quality_metrics, ...
                                                                  {'diagonal_stats', 'mean'}, NaN);
                    if ~isnan(diagonal_mean)
                        summary_lines{end+1} = sprintf('Mean Diagonal: %.3f', diagonal_mean);
                    end
                    
                    cond_mean = Module1Utils.safe_field_access(quality_metrics, ...
                                                              {'conditioning', 'mean'}, NaN);
                    if ~isnan(cond_mean)
                        summary_lines{end+1} = sprintf('Mean Condition Number: %.1f', cond_mean);
                    end
                end
                
                % Timing information
                total_time = Module1Utils.safe_field_access(results, 'total_time', NaN);
                if ~isnan(total_time)
                    summary_lines{end+1} = sprintf('Processing Time: %.2f seconds', total_time);
                end
                
            catch ME
                summary_lines{end+1} = sprintf('Summary generation failed: %s', ME.message);
            end
            
            summary_text = strjoin(summary_lines, '\n');
        end
        
        function data = convert_simulation_data(sim_results)
            % CONVERT_SIMULATION_DATA - Convert module7 output to module1 input
            %
            % Input:
            %   sim_results - Output from module7_simulation
            %
            % Output:
            %   data - Properly formatted data for module1
            
            data = struct();
            
            try
                % Direct field mapping
                data.Sigma_emp = sim_results.Sigma_emp;
                data.F = sim_results.F;
                data.n = sim_results.n;
                data.T = sim_results.T;
                
                % Validate conversion
                if ~iscell(data.Sigma_emp)
                    error('Sigma_emp must be a cell array');
                end
                
                if length(data.Sigma_emp) ~= data.F
                    error('Inconsistent dimensions: Sigma_emp length vs F');
                end
                
                if size(data.Sigma_emp{1}, 1) ~= data.n
                    error('Inconsistent dimensions: matrix size vs n');
                end
                
            catch ME
                error('Data conversion failed: %s', ME.message);
            end
        end
        
    end
end