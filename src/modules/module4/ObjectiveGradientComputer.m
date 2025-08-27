classdef ObjectiveGradientComputer < handle
    % OBJECTIVE_GRADIENT_COMPUTER - Encapsulates objective and gradient computations
    %
    % Description:
    %   Provides a unified interface for computing objectives and gradients in the
    %   sparse precision matrix estimation problem. Supports both direct and
    %   graph Laplacian methods with comprehensive validation and statistics.
    %
    %   Key Features:
    %   - Numerically stable Cholesky-based computations
    %   - Automatic Hermitian symmetry enforcement
    %   - Flexible L1 penalty configuration
    %   - Performance monitoring and error handling
    %   - Cross-method validation capabilities
    %
    % Properties:
    %   lambda1 - Smoothing regularization parameter
    %   lambda2 - L1 sparsity parameter
    %   computation_options - Configuration options structure
    %   last_computation_stats - Statistics from most recent computation
    %
    % Methods:
    %   compute_smooth_gradients - Compute smooth part gradients only
    %   compute_full_objective - Compute complete objective function
    %   validate_inputs - Check input data consistency
    %   compare_gradient_methods - Compare direct vs Laplacian implementations
    %
    % Static Methods:
    %   create_default_options - Generate default configuration
    %   verify_hermitian_properties - Validate matrix Hermitian properties
    %
    % Examples:
    %   % Create computer with default settings
    %   computer = ObjectiveGradientComputer();
    %
    %   % Configure parameters
    %   computer.lambda1 = 0.05;
    %   computer.lambda2 = 0.02;
    %   computer.computation_options.verbose = true;
    %
    %   % Compute gradients
    %   gradients = computer.compute_smooth_gradients(input_data);
    %
    %   % Evaluate objective
    %   objective = computer.compute_full_objective(input_data);
    %
    %   % Performance analysis
    %   stats = computer.last_computation_stats;
    %   fprintf('Computation time: %.3fs\n', stats.total_time);
    %
    % See also: MODULE4_OBJECTIVE_GRADIENT_MAIN, MODULE4_OBJECTIVE_EVALUATION

    properties (Access = public)
        lambda1                 % Smoothing parameter
        lambda2                 % L1 penalty parameter
        computation_options     % Configuration options
        last_computation_stats  % Statistics from last computation
    end

    properties (Access = private)
        input_cache            % Cached input data for efficiency
        gradient_cache         % Cached gradient results
        objective_cache        % Cached objective results
        cache_valid           % Cache validity flags
    end

    methods (Access = public)

        function obj = ObjectiveGradientComputer(varargin)
            % OBJECTIVE_GRADIENT_COMPUTER - Constructor
            %
            % Syntax:
            %   computer = ObjectiveGradientComputer()
            %   computer = ObjectiveGradientComputer(Name, Value, ...)
            %
            % Name-Value Arguments:
            %   lambda1 - (double) Smoothing parameter (default: 0.01)
            %   lambda2 - (double) L1 penalty parameter (default: 0.01)
            %   options - (struct) Computation options (default: created automatically)

            % Parse input arguments
            p = inputParser;
            addParameter(p, 'lambda1', 0.01, @(x) isnumeric(x) && isscalar(x) && x >= 0);
            addParameter(p, 'lambda2', 0.01, @(x) isnumeric(x) && isscalar(x) && x >= 0);
            addParameter(p, 'options', [], @(x) isstruct(x) || isempty(x));
            parse(p, varargin{:});

            % Set properties
            obj.lambda1 = p.Results.lambda1;
            obj.lambda2 = p.Results.lambda2;

            if isempty(p.Results.options)
                obj.computation_options = ObjectiveGradientComputer.create_default_options();
            else
                obj.computation_options = p.Results.options;
            end

            % Initialize private properties
            obj.input_cache = struct();
            obj.gradient_cache = struct();
            obj.objective_cache = struct();
            obj.cache_valid = struct('gradients', false, 'objective', false);
            obj.last_computation_stats = struct();

            % Validate configuration
            obj.validate_configuration();
        end

        function gradients = compute_smooth_gradients(obj, input_data, options)
            % COMPUTE_SMOOTH_GRADIENTS - Compute gradients of smooth objective parts
            %
            % Syntax:
            %   gradients = compute_smooth_gradients(obj, input_data)
            %   gradients = compute_smooth_gradients(obj, input_data, options)

            if nargin < 3
                options = obj.computation_options;
            else
                % Merge with default options
                default_opts = obj.computation_options;
                merged_opts = obj.merge_options(default_opts, options);
                options = merged_opts;
            end

            % Validate inputs
            obj.validate_inputs(input_data);

            % Check cache validity
            if obj.cache_valid.gradients && obj.input_data_unchanged(input_data)
                if options.verbose
                    fprintf('Using cached gradient results\n');
                end
                gradients = obj.gradient_cache.smooth_gradients;
                obj.last_computation_stats = obj.gradient_cache.computation_stats;
                return;
            end

            % Prepare parameters for gradient computation
            gradient_params = options;
            gradient_params.lambda1 = obj.lambda1;

            % Compute gradients using main module
            gradient_results = module4_objective_gradient_main(input_data, gradient_params);

            % Extract results
            gradients = gradient_results.smooth_gradients;

            % Update cache
            obj.gradient_cache = gradient_results;
            obj.input_cache = input_data;
            obj.cache_valid.gradients = true;
            obj.last_computation_stats = gradient_results.computation_stats;

            % Validate results
            obj.validate_gradient_results(gradients, input_data);
        end

        function objective_values = compute_full_objective(obj, input_data, options)
            % COMPUTE_FULL_OBJECTIVE - Compute complete objective function value
            %
            % Syntax:
            %   objective_values = compute_full_objective(obj, input_data)
            %   objective_values = compute_full_objective(obj, input_data, options)

            if nargin < 3
                options = obj.computation_options;
            else
                % Merge with default options
                default_opts = obj.computation_options;
                merged_opts = obj.merge_options(default_opts, options);
                options = merged_opts;
            end

            % Validate inputs
            obj.validate_inputs(input_data);

            % Check cache validity
            if obj.cache_valid.objective && obj.input_data_unchanged(input_data)
                if options.verbose
                    fprintf('Using cached objective results\n');
                end
                objective_values = obj.objective_cache.objective_values;
                obj.last_computation_stats = obj.objective_cache.computation_stats;
                return;
            end

            % Prepare parameters for objective computation
            evaluation_params = options;
            evaluation_params.lambda1 = obj.lambda1;
            evaluation_params.lambda2 = obj.lambda2;

            % Compute objective using evaluation module
            [objective_vals, computation_stats] = module4_objective_evaluation(input_data, evaluation_params);

            % Extract results
            objective_values = objective_vals;

            % Update cache
            obj.objective_cache.objective_values = objective_values;
            obj.objective_cache.computation_stats = computation_stats;
            obj.input_cache = input_data;
            obj.cache_valid.objective = true;
            obj.last_computation_stats = computation_stats;
        end

        function comparison_results = compare_gradient_methods(obj, input_data)
            % COMPARE_GRADIENT_METHODS - Compare direct vs Laplacian gradient methods
            %
            % Syntax:
            %   comparison_results = compare_gradient_methods(obj, input_data)
            %
            % Output:
            %   comparison_results - Structure with comparison statistics

            fprintf('Comparing gradient computation methods...\n');

            % Compute gradients with direct method
            opts_direct = obj.computation_options;
            opts_direct.use_graph_laplacian = false;
            opts_direct.verbose = false;

            tic_direct = tic;
            gradients_direct = obj.compute_smooth_gradients(input_data, opts_direct);
            time_direct = toc(tic_direct);

            % Clear cache to force recomputation
            obj.invalidate_cache();

            % Compute gradients with Laplacian method
            opts_laplacian = obj.computation_options;
            opts_laplacian.use_graph_laplacian = true;
            opts_laplacian.verbose = false;

            tic_laplacian = tic;
            gradients_laplacian = obj.compute_smooth_gradients(input_data, opts_laplacian);
            time_laplacian = toc(tic_laplacian);

            % Compare results
            F = length(gradients_direct);
            max_difference = 0;
            mean_difference = 0;
            relative_errors = zeros(F, 1);

            for f = 1:F
                diff_f = norm(gradients_direct{f} - gradients_laplacian{f}, 'fro');
                norm_f = norm(gradients_direct{f}, 'fro');

                max_difference = max(max_difference, diff_f);
                mean_difference = mean_difference + diff_f;

                if norm_f > 1e-12
                    relative_errors(f) = diff_f / norm_f;
                else
                    relative_errors(f) = diff_f;
                end
            end
            mean_difference = mean_difference / F;

            % Create comparison results
            comparison_results = struct();
            comparison_results.time_direct = time_direct;
            comparison_results.time_laplacian = time_laplacian;
            comparison_results.speedup_ratio = time_direct / time_laplacian;
            comparison_results.max_absolute_difference = max_difference;
            comparison_results.mean_absolute_difference = mean_difference;
            comparison_results.max_relative_error = max(relative_errors);
            comparison_results.mean_relative_error = mean(relative_errors);
            comparison_results.methods_consistent = max_difference < 1e-10;

            % Display results
            fprintf('====================================\n');
            fprintf('Gradient Method Comparison Results:\n');
            fprintf('  Direct method time:     %.3fs\n', time_direct);
            fprintf('  Laplacian method time:  %.3fs\n', time_laplacian);
            fprintf('  Speedup ratio:          %.2fx\n', comparison_results.speedup_ratio);
            fprintf('  Max absolute diff:      %.2e\n', max_difference);
            fprintf('  Max relative error:     %.2e\n', comparison_results.max_relative_error);
            if comparison_results.methods_consistent
                fprintf('  Methods consistent:     YES\n');
            else
                fprintf('  Methods consistent:     NO\n');
            end
            fprintf('====================================\n');
        end

        function validate_inputs(obj, input_data)
            % VALIDATE_INPUTS - Comprehensive input validation
            %
            % Validates all aspects of input_data structure for consistency
            % and mathematical properties.

            % Check required fields
            required_fields = {'precision_matrices', 'whitened_covariances', ...
                'smoothing_kernel', 'weight_matrix'};
            for i = 1:length(required_fields)
                if ~isfield(input_data, required_fields{i})
                    error('ObjectiveGradientComputer:missing_field', ...
                        'Required field "%s" not found', required_fields{i});
                end
            end

            % Extract data
            Gammas = input_data.precision_matrices;
            Sigmas = input_data.whitened_covariances;
            K = input_data.smoothing_kernel;
            W = input_data.weight_matrix;

            % Validate basic properties
            if ~iscell(Gammas) || ~iscell(Sigmas) || isempty(Gammas) || isempty(Sigmas)
                error('ObjectiveGradientComputer:invalid_cell_arrays', ...
                    'precision_matrices and whitened_covariances must be non-empty cell arrays');
            end

            F = length(Gammas);
            p = size(Gammas{1}, 1);

            if length(Sigmas) ~= F
                error('ObjectiveGradientComputer:length_mismatch', ...
                    'precision_matrices and whitened_covariances must have same length');
            end

            % Validate matrix dimensions and properties
            for f = 1:F
                % Check dimensions
                if ~isequal(size(Gammas{f}), [p, p]) || ~isequal(size(Sigmas{f}), [p, p])
                    error('ObjectiveGradientComputer:dimension_mismatch', ...
                        'All matrices must be %dx%d at frequency %d', p, p, f);
                end

                % Check positive definiteness of precision matrices
                try
                    chol(Gammas{f});
                catch
                    error('ObjectiveGradientComputer:not_positive_definite', ...
                        'precision_matrices{%d} is not positive definite', f);
                end

                % Check Hermitian properties
                if norm(Gammas{f} - Gammas{f}', 'fro') > 1e-10
                    warning('ObjectiveGradientComputer:not_hermitian', ...
                        'precision_matrices{%d} is not Hermitian', f);
                end

                if norm(Sigmas{f} - Sigmas{f}', 'fro') > 1e-10
                    warning('ObjectiveGradientComputer:covariance_not_hermitian', ...
                        'whitened_covariances{%d} is not Hermitian', f);
                end
            end

            % Validate kernel matrix
            if ~isnumeric(K) || ~isequal(size(K), [F, F])
                error('ObjectiveGradientComputer:invalid_kernel', ...
                    'smoothing_kernel must be %dx%d numeric matrix', F, F);
            end

            if norm(K - K', 'fro') > 1e-10
                warning('ObjectiveGradientComputer:kernel_not_symmetric', ...
                    'smoothing_kernel should be symmetric');
            end

            % Validate weight matrix
            if ~isnumeric(W) || ~isequal(size(W), [p, p])
                error('ObjectiveGradientComputer:invalid_weight', ...
                    'weight_matrix must be %dx%d numeric matrix', p, p);
            end

            ObjectiveGradientComputer.verify_hermitian_properties(W, 'weight_matrix');
        end

        function invalidate_cache(obj)
            % INVALIDATE_CACHE - Clear all cached results
            obj.cache_valid.gradients = false;
            obj.cache_valid.objective = false;
        end

    end

    methods (Access = private)

        function validate_configuration(obj)
            % Validate object configuration
            if obj.lambda1 < 0 || obj.lambda2 < 0
                error('ObjectiveGradientComputer:invalid_lambda', ...
                    'lambda1 and lambda2 must be non-negative');
            end

            if ~isstruct(obj.computation_options)
                error('ObjectiveGradientComputer:invalid_options', ...
                    'computation_options must be a struct');
            end
        end

        function unchanged = input_data_unchanged(obj, input_data)
            % Check if input data is the same as cached version
            unchanged = false;

            if isempty(obj.input_cache)
                return;
            end

            % Simple check - compare data structures
            try
                % Check precision matrices
                if length(input_data.precision_matrices) ~= length(obj.input_cache.precision_matrices)
                    return;
                end

                for f = 1:length(input_data.precision_matrices)
                    if ~isequal(input_data.precision_matrices{f}, obj.input_cache.precision_matrices{f})
                        return;
                    end
                end

                % Check other fields
                if ~isequal(input_data.smoothing_kernel, obj.input_cache.smoothing_kernel) || ...
                        ~isequal(input_data.weight_matrix, obj.input_cache.weight_matrix)
                    return;
                end

                unchanged = true;
            catch
                unchanged = false;
            end
        end

        function merged_options = merge_options(obj, default_opts, user_opts)
            % Merge user options with defaults
            merged_options = default_opts;

            if isstruct(user_opts)
                field_names = fieldnames(user_opts);
                for i = 1:numel(field_names)
                    fname = field_names{i};
                    merged_options.(fname) = user_opts.(fname);
                end
            end
        end

        function validate_gradient_results(obj, gradients, input_data)
            % Validate computed gradient results
            F = length(gradients);
            p = size(gradients{1}, 1);

            for f = 1:F
                % Check dimensions
                if ~isequal(size(gradients{f}), [p, p])
                    error('ObjectiveGradientComputer:invalid_gradient_size', ...
                        'Gradient %d has incorrect size', f);
                end

                % Check Hermitian property
                if obj.computation_options.force_hermitian
                    hermitian_error = norm(gradients{f} - gradients{f}', 'fro');
                    if hermitian_error > 1e-8
                        warning('ObjectiveGradientComputer:gradient_not_hermitian', ...
                            'Gradient %d is not Hermitian (error: %.2e)', f, hermitian_error);
                    end
                end

                % Check for NaN or Inf
                if any(~isfinite(gradients{f}(:)))
                    error('ObjectiveGradientComputer:gradient_not_finite', ...
                        'Gradient %d contains NaN or Inf values', f);
                end
            end
        end

    end

    methods (Static)

        function options = create_default_options()
            % CREATE_DEFAULT_OPTIONS - Unified defaults for all computations
            options = struct();
            options.penalize_diagonal         = false;
            options.use_graph_laplacian       = true;
            options.chol_tolerance            = 1e-12;
            options.symmetrization_tolerance  = 1e-10;
            options.force_hermitian           = true;
            options.verbose                   = false;

            % NEW unified defaults to keep interfaces consistent
            options.kernel_zero_tol           = 1e-12;   % drop tiny kernel entries
            options.strict_input_validation   = true;    % allow loosening if desired
            options.allow_pinv_fallback       = false;   % diagnostics only; default off
        end




        function verify_hermitian_properties(matrix, matrix_name)
            % VERIFY_HERMITIAN_PROPERTIES - Validate Hermitian matrix properties
            %
            % Checks if a matrix is Hermitian and positive semi-definite

            if nargin < 2
                matrix_name = 'matrix';
            end

            % Check Hermitian property
            hermitian_error = norm(matrix - matrix', 'fro');
            if hermitian_error > 1e-10
                warning('ObjectiveGradientComputer:not_hermitian', ...
                    '%s is not Hermitian (error: %.2e)', matrix_name, hermitian_error);
            end

            % Check positive semi-definiteness
            try
                eigenvals = eig(matrix);
                min_eigenval = min(real(eigenvals));

                if min_eigenval < -1e-10
                    warning('ObjectiveGradientComputer:not_psd', ...
                        '%s may not be positive semi-definite (min eigenvalue: %.2e)', ...
                        matrix_name, min_eigenval);
                end
            catch ME
                warning('ObjectiveGradientComputer:eigenvalue_computation_failed', ...
                    'Could not compute eigenvalues of %s: %s', matrix_name, ME.message);
            end
        end

    end

end