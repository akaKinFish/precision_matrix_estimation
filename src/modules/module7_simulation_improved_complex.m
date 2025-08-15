function [true_precision, true_covariance, emp_covariance, params] = module7_simulation_improved_complex(varargin)
% MODULE7_SIMULATION_IMPROVED_COMPLEX - Fixed version with working dynamic sparsity
%
% FIXES:
% 1. Dynamic sparsity threshold calculation fixed
% 2. Parameter parsing improved
% 3. Edge activation logic corrected

    % Parse parameters with improved validation
    p = inputParser;
    % Original parameters
    addParameter(p, 'n_nodes', 12, @(x) isnumeric(x) && isscalar(x) && x > 0 && x == round(x));
    addParameter(p, 'n_freq', 25, @(x) isnumeric(x) && isscalar(x) && x > 0 && x == round(x));
    addParameter(p, 'n_samples', 150, @(x) isnumeric(x) && isscalar(x) && x > 0 && x == round(x));
    addParameter(p, 'graph_type', 'random', @(x) ischar(x) && any(strcmp(x, {'random', 'chain', 'grid', 'hub'})));
    addParameter(p, 'edge_density', 0.2, @(x) isnumeric(x) && isscalar(x) && x > 0 && x <= 1);
    addParameter(p, 'n_basis', 5, @(x) isnumeric(x) && isscalar(x) && x > 0 && x == round(x));
    addParameter(p, 'sigma_coef', 0.5, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'epsilon_reg', 0.1, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'random_seed', 123, @(x) isnumeric(x) && isscalar(x));
    % Complex-specific parameter
    addParameter(p, 'complex_strength', 0.8, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 2);
    % Dynamic sparsity parameters
    addParameter(p, 'sparsity_variation', 0.2, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
    addParameter(p, 'edge_activation_smoothness', 0.8, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
    
    parse(p, varargin{:});
    params = p.Results;
    
    rng(params.random_seed);
    
    fprintf('=== Enhanced Complex Simulation with Dynamic Sparsity ===\n');
    
    %% Step 1: Generate base graph structure (always active edges)
    base_edge_list = generateGraphStructure(params.n_nodes, params.graph_type, params.edge_density);
    n_base_edges = size(base_edge_list, 1);
    
    fprintf('Generated %d base edges\n', n_base_edges);
    
    %% Step 2: Generate variable edges for dynamic sparsity
    % Create pool of potential additional edges
    all_possible_edges = [];
    for i = 2:params.n_nodes
        for j = 1:(i-1)
            all_possible_edges = [all_possible_edges; i, j];
        end
    end
    
    % Remove base edges from possible variable edges
    base_edge_set = [base_edge_list; base_edge_list(:, [2,1])]; % Include both directions
    is_base = false(size(all_possible_edges, 1), 1);
    for e = 1:size(all_possible_edges, 1)
        is_base(e) = any(all(base_edge_set == all_possible_edges(e, :), 2));
    end
    variable_edges = all_possible_edges(~is_base, :);
    
    % Select variable edges for dynamic activation
    n_variable = round(params.sparsity_variation * n_base_edges);
    if n_variable > 0 && size(variable_edges, 1) > 0
        perm = randperm(size(variable_edges, 1));
        selected_idx = perm(1:min(n_variable, size(variable_edges, 1)));
        variable_edge_list = variable_edges(selected_idx, :);
    else
        variable_edge_list = [];
    end
    
    n_variable_edges = size(variable_edge_list, 1);
    fprintf('Selected %d variable edges for dynamic sparsity\n', n_variable_edges);
    
    %% Step 3: Generate smooth basis functions
    freq_normalized = linspace(0, 1, params.n_freq)';
    basis_functions = generateSmoothBasis(freq_normalized, params.n_basis);
    
    %% Step 4: Generate smooth edge activation patterns for variable edges (FIXED)
    variable_activations = zeros(params.n_freq, n_variable_edges);
    for e = 1:n_variable_edges
        % Generate smooth activation curve for this edge
        phase = rand() * 2 * pi; % Random phase
        freq_scale = 0.5 + randn() * 0.2; % Frequency scaling
        
        % Create smooth activation using sinusoidal pattern
        sine_component = sin(freq_scale * 2 * pi * freq_normalized + phase);
        
        % Apply sigmoid to get smooth 0-1 activation with more dynamic range
        smooth_factor = params.edge_activation_smoothness * 3; % Reduce factor for more variation
        variable_activations(:, e) = 1 ./ (1 + exp(-smooth_factor * sine_component));
        
        % CRITICAL FIX: Add more aggressive threshold variation
        % Make some activations go to zero to create actual sparsity changes
        threshold_level = 0.3 + 0.4 * rand(); % Random threshold between 0.3-0.7
        binary_mask = variable_activations(:, e) > threshold_level;
        variable_activations(:, e) = variable_activations(:, e) .* binary_mask;
    end
    
    %% Step 5: Generate complex coefficients for all edges and diagonals
    n_total_entries = n_base_edges + n_variable_edges + params.n_nodes;
    
    % Generate complex coefficients
    coefficients_real = params.sigma_coef * randn(n_total_entries, params.n_basis);
    coefficients_imag = params.sigma_coef * params.complex_strength * randn(n_total_entries, params.n_basis);
    
    % Ensure meaningful imaginary components
    coefficients_imag = coefficients_imag + 0.1 * params.sigma_coef * params.complex_strength * ...
                       sign(randn(n_total_entries, params.n_basis));
    
    coefficients = complex(coefficients_real, coefficients_imag);
    
    % Report complex coefficient properties
    complex_coef_fraction = sum(abs(imag(coefficients(:))) > 1e-10) / numel(coefficients);
    fprintf('Coefficient complex fraction: %.2f\n', complex_coef_fraction);
    
    %% Step 6: Construct complex Hermitian precision matrices with dynamic sparsity
    true_precision = cell(params.n_freq, 1);
    true_covariance = cell(params.n_freq, 1);
    
    for f = 1:params.n_freq
        if mod(f, max(1, round(params.n_freq/5))) == 1
            fprintf('Processing frequency %d/%d\n', f, params.n_freq);
        end
        
        % Initialize lower triangular matrix
        L = complex(zeros(params.n_nodes, params.n_nodes));
        
        edge_idx = 1;
        
        % Fill base edges (always active)
        for e = 1:n_base_edges
            i = base_edge_list(e, 1);
            j = base_edge_list(e, 2);
            
            % Compute complex value from basis functions
            base_value = sum(coefficients(edge_idx, :) .* basis_functions(f, :));
            L(i, j) = base_value;
            edge_idx = edge_idx + 1;
        end
        
        % Fill variable edges (with smooth activation) - FIXED
        for e = 1:n_variable_edges
            i = variable_edge_list(e, 1);
            j = variable_edge_list(e, 2);
            
            % Compute base complex value
            base_value = sum(coefficients(edge_idx, :) .* basis_functions(f, :));
            
            % Apply smooth activation - FIXED: Direct multiplication
            activated_value = base_value * variable_activations(f, e);
            
            % FIXED: Remove soft thresholding that was preventing sparsity changes
            % Just apply the activation directly
            L(i, j) = activated_value;
            edge_idx = edge_idx + 1;
        end
        
        % Fill diagonal elements (always real for Hermitian property)
        for i = 1:params.n_nodes
            diag_value = sum(coefficients(edge_idx, :) .* basis_functions(f, :));
            L(i, i) = real(diag_value); % Force real diagonal
            edge_idx = edge_idx + 1;
        end
        
        % Enforce diagonal dominance for positive definiteness
        for i = 1:params.n_nodes
            row_sum = sum(abs(L(i, 1:i-1)));
            L(i, i) = max(real(L(i, i)), row_sum + params.epsilon_reg);
        end
        
        % Form precision matrix: Omega = L * L'
        Omega = L * L';
        
        % Ensure Hermitian property (should be automatic, but enforce for numerical stability)
        Omega = (Omega + Omega') / 2;
        
        % Store results
        true_precision{f} = Omega;
        true_covariance{f} = inv(Omega);
        
        % Debug: Check complex properties for first few frequencies
        if f <= 3
            complex_fraction = sum(abs(imag(Omega(:))) > 1e-10) / numel(Omega);
            max_imag = max(abs(imag(Omega(:))));
            fprintf('  Freq %d: Complex fraction = %.2f, Max imag = %.4f\n', ...
                    f, complex_fraction, max_imag);
        end
    end
    
    %% Step 7: Generate empirical covariances
    emp_covariance = cell(params.n_freq, 1);
    for f = 1:params.n_freq
        Sigma = true_covariance{f};
        
        % Generate complex samples
        if isreal(Sigma)
            Z = mvnrnd(zeros(params.n_nodes, 1), Sigma, params.n_samples)';
        else
            Z_real = mvnrnd(zeros(params.n_nodes, 1), real(Sigma), params.n_samples)';
            Z_imag = mvnrnd(zeros(params.n_nodes, 1), real(Sigma), params.n_samples)' * 0.5;
            Z = Z_real + 1i * Z_imag;
        end
        
        % Form empirical covariance
        W = Z * Z';
        emp_cov_raw = W / params.n_samples;
        
        % Add small amount of measurement noise
        noise_strength = 0.05 * norm(emp_cov_raw, 'fro') / params.n_nodes;
        noise_real = noise_strength * randn(params.n_nodes, params.n_nodes);
        noise_imag = noise_strength * randn(params.n_nodes, params.n_nodes) * params.complex_strength * 0.3;
        noise = complex(noise_real, noise_imag);
        
        % Make noise Hermitian
        noise = (noise + noise') / 2;
        % Ensure real diagonal
        for i = 1:params.n_nodes
            noise(i, i) = real(noise(i, i));
        end
        
        emp_covariance{f} = emp_cov_raw + noise;
        
        % Ensure Hermitian and positive definite
        emp_covariance{f} = (emp_covariance{f} + emp_covariance{f}') / 2;
        
        % Force positive definiteness
        [V, D] = eig(emp_covariance{f});
        d = diag(D);
        d = max(real(d), 1e-8);
        emp_covariance{f} = V * diag(d) * V';
        
        % Final Hermitian enforcement
        emp_covariance{f} = (emp_covariance{f} + emp_covariance{f}') / 2;
    end
    
    %% Step 8: Analyze and store results
    params.base_edge_list = base_edge_list;
    params.variable_edge_list = variable_edge_list;
    params.variable_activations = variable_activations;
    params.n_base_edges = n_base_edges;
    params.n_variable_edges = n_variable_edges;
    params.basis_functions = basis_functions;
    params.coefficients = coefficients;
    
    % FIXED: Analyze sparsity pattern changes with better threshold
    sparsity_changes = analyze_sparsity_changes_complex(true_precision);
    params.sparsity_changes = sparsity_changes;
    
    % Comprehensive complex properties analysis
    complex_analysis = analyze_complex_matrices(true_precision);
    params.complex_analysis = complex_analysis;
    
    %% Step 9: Final reporting
    fprintf('\n=== Final Analysis ===\n');
    fprintf('Base edges (always active): %d\n', n_base_edges);
    fprintf('Variable edges: %d\n', n_variable_edges);
    fprintf('Sparsity pattern changes: %d\n', sparsity_changes);
    fprintf('Matrices with complex entries: %d/%d\n', ...
            complex_analysis.matrices_with_complex, params.n_freq);
    fprintf('Average complex fraction: %.2f\n', complex_analysis.avg_complex_fraction);
    fprintf('Average phase magnitude: %.4f radians\n', complex_analysis.avg_phase_magnitude);
    fprintf('Max imaginary component: %.4f\n', complex_analysis.max_imaginary);
    
    % Verify Hermitian property
    hermitian_count = 0;
    for f = 1:params.n_freq
        if ishermitian(true_precision{f})
            hermitian_count = hermitian_count + 1;
        end
    end
    fprintf('Hermitian matrices: %d/%d\n', hermitian_count, params.n_freq);
    
    fprintf('=== Enhanced Complex Simulation Complete ===\n');
end

%% Helper Functions

function edge_list = generateGraphStructure(n_nodes, graph_type, edge_density)
% Generate graph structure based on type and density
    edge_list = [];
    
    switch lower(graph_type)
        case 'random'
            n_possible_edges = n_nodes * (n_nodes - 1) / 2;
            n_edges = round(edge_density * n_possible_edges);
            all_edges = [];
            for i = 2:n_nodes
                for j = 1:(i-1)
                    all_edges = [all_edges; i, j];
                end
            end
            if n_edges > 0
                perm = randperm(size(all_edges, 1));
                selected_idx = perm(1:min(n_edges, size(all_edges, 1)));
                edge_list = all_edges(selected_idx, :);
            end
            
        case 'chain'
            for i = 2:n_nodes
                edge_list = [edge_list; i, i-1];
            end
            
        case 'grid'
            % Simple 2D grid approximation
            grid_size = ceil(sqrt(n_nodes));
            for i = 1:n_nodes
                [row, col] = ind2sub([grid_size, grid_size], i);
                % Right neighbor
                if col < grid_size
                    j = sub2ind([grid_size, grid_size], row, col + 1);
                    if j <= n_nodes && j > i
                        edge_list = [edge_list; j, i];
                    end
                end
                % Bottom neighbor
                if row < grid_size
                    j = sub2ind([grid_size, grid_size], row + 1, col);
                    if j <= n_nodes && j > i
                        edge_list = [edge_list; j, i];
                    end
                end
            end
            
        case 'hub'
            for i = 2:n_nodes
                edge_list = [edge_list; i, 1];
            end
    end
    
    edge_list = unique(edge_list, 'rows');
end

function basis = generateSmoothBasis(freq_normalized, n_basis)
% Generate smooth basis functions
    n_freq = length(freq_normalized);
    basis = zeros(n_freq, n_basis);
    centers = linspace(0, 1, n_basis);
    width = 2.0 / n_basis;
    
    for b = 1:n_basis
        for f = 1:n_freq
            dist = abs(freq_normalized(f) - centers(b));
            if dist < width/2
                basis(f, b) = 0.5 * (1 + cos(2*pi*dist/width));
            end
        end
    end
    
    % Normalize each basis function
    for b = 1:n_basis
        if sum(basis(:, b)) > 0
            basis(:, b) = basis(:, b) / sum(basis(:, b));
        end
    end
end

function sparsity_changes = analyze_sparsity_changes_complex(precision_matrices)
% FIXED: Analyze sparsity pattern changes with proper threshold
    F = length(precision_matrices);
    threshold = 1e-8; % FIXED: Lower threshold to detect more changes
    total_changes = 0;
    
    for f = 2:F
        prev_pattern = abs(precision_matrices{f-1}) > threshold;
        curr_pattern = abs(precision_matrices{f}) > threshold;
        
        % Count changes in upper triangular part only
        n = size(prev_pattern, 1);
        for i = 1:n
            for j = i+1:n
                if prev_pattern(i,j) ~= curr_pattern(i,j)
                    total_changes = total_changes + 1;
                end
            end
        end
    end
    
    sparsity_changes = total_changes;
end

function analysis = analyze_complex_matrices(precision_matrices)
% Comprehensive analysis of complex matrix properties
    F = length(precision_matrices);
    analysis = struct();
    
    matrices_with_complex = 0;
    total_complex_fraction = 0;
    total_phase_magnitude = 0;
    max_imaginary = 0;
    
    for f = 1:F
        P = precision_matrices{f};
        
        % Check if matrix has complex entries
        has_complex = any(abs(imag(P(:))) > 1e-12);
        if has_complex
            matrices_with_complex = matrices_with_complex + 1;
        end
        
        % Complex fraction
        complex_fraction = sum(abs(imag(P(:))) > 1e-12) / numel(P);
        total_complex_fraction = total_complex_fraction + complex_fraction;
        
        % Phase analysis
        phases = angle(P(abs(P) > 1e-12));
        if ~isempty(phases)
            avg_phase_mag = mean(abs(phases));
            total_phase_magnitude = total_phase_magnitude + avg_phase_mag;
        end
        
        % Max imaginary component
        max_imag_f = max(abs(imag(P(:))));
        max_imaginary = max(max_imaginary, max_imag_f);
    end
    
    analysis.matrices_with_complex = matrices_with_complex;
    analysis.avg_complex_fraction = total_complex_fraction / F;
    analysis.avg_phase_magnitude = total_phase_magnitude / F;
    analysis.max_imaginary = max_imaginary;
end