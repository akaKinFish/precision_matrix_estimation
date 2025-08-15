function [true_precision, true_covariance, emp_covariance, params] = module7_simulation_improved(varargin)
% MODULE7_SIMULATION_IMPROVED - Enhanced version with dynamic sparsity patterns
%
% This version generates precision matrices where:
% 1. Values change smoothly across frequencies
% 2. Sparsity patterns can evolve (edges can appear/disappear smoothly)
% 3. The overall structure remains coherent
%
% Additional parameters:
%   'sparsity_variation' - Controls how much sparsity pattern changes (0-1, default: 0.2)
%   'edge_activation_smoothness' - Controls edge activation smoothness (default: 0.8)

    % Parse all original parameters plus new ones with simplified validation
    p = inputParser;
    % Original parameters - using validateattributes for consistency
    addParameter(p, 'n_nodes', 20, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
    addParameter(p, 'n_freq', 50, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
    addParameter(p, 'n_samples', 100, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
    addParameter(p, 'graph_type', 'random', @(x) any(validatestring(x, {'random', 'chain', 'grid', 'hub'})));
    addParameter(p, 'edge_density', 0.1, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', '<=', 1}));
    addParameter(p, 'n_basis', 5, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
    addParameter(p, 'sigma_coef', 0.5, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
    addParameter(p, 'epsilon_reg', 0.1, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
    addParameter(p, 'random_seed', 42, @(x) validateattributes(x, {'numeric'}, {'scalar'}));
    % New parameters
    addParameter(p, 'sparsity_variation', 0.2, @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0, '<=', 1}));
    addParameter(p, 'edge_activation_smoothness', 0.8, @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0, '<=', 1}));
    
    parse(p, varargin{:});
    params = p.Results;
    
    % Set random seed
    rng(params.random_seed);
    
    % Step 1: Generate base graph structure (core edges that are always active)
    base_edge_list = generateGraphStructure(params.n_nodes, params.graph_type, params.edge_density);
    n_base_edges = size(base_edge_list, 1);
    
    % Step 2: Generate potential additional edges for variation
    % These edges can be activated/deactivated across frequencies
    all_possible_edges = [];
    for i = 2:params.n_nodes
        for j = 1:(i-1)
            all_possible_edges = [all_possible_edges; i, j];
        end
    end
    
    % Remove base edges from possible edges
    base_edge_set = [base_edge_list; base_edge_list(:, [2,1])]; % Include both directions
    is_base = false(size(all_possible_edges, 1), 1);
    for e = 1:size(all_possible_edges, 1)
        is_base(e) = any(all(base_edge_set == all_possible_edges(e, :), 2));
    end
    variable_edges = all_possible_edges(~is_base, :);
    
    % Select some variable edges
    n_variable = round(params.sparsity_variation * n_base_edges);
    if n_variable > 0 && size(variable_edges, 1) > 0
        perm = randperm(size(variable_edges, 1));
        selected_idx = perm(1:min(n_variable, size(variable_edges, 1)));
        variable_edge_list = variable_edges(selected_idx, :);
    else
        variable_edge_list = [];
    end
    
    % Step 3: Generate smooth activation patterns for variable edges
    freq_normalized = linspace(0, 1, params.n_freq)';
    
    % Generate smooth activation functions for each variable edge
    variable_activations = zeros(params.n_freq, size(variable_edge_list, 1));
    for e = 1:size(variable_edge_list, 1)
        % Create smooth activation pattern using sigmoid-like functions
        center = rand(); % Random activation center
        width = params.edge_activation_smoothness * 0.5;
        
        % Smooth activation function
        activation = 1 ./ (1 + exp(-10*(freq_normalized - center)/width));
        
        % Add some randomness
        activation = activation + 0.1*randn(params.n_freq, 1);
        activation = smooth(activation, 5); % Smooth the randomness
        
        % Threshold to binary (but keep some values near threshold)
        threshold = 0.5 + 0.2*(rand() - 0.5);
        soft_threshold = 0.1;
        activation_prob = 1 ./ (1 + exp(-(activation - threshold)/soft_threshold));
        
        variable_activations(:, e) = activation_prob;
    end
    
    % Step 4: Generate basis functions for coefficient values
    basis_functions = generateSmoothBasis(freq_normalized, params.n_basis);
    
    % Step 5: Generate coefficients for all edges (base + variable)
    total_entries = n_base_edges + size(variable_edge_list, 1) + params.n_nodes;
    coefficients = params.sigma_coef * randn(total_entries, params.n_basis);
    
    % Step 6: Construct precision matrices
    true_precision = cell(params.n_freq, 1);
    true_covariance = cell(params.n_freq, 1);
    
    for f = 1:params.n_freq
        % Initialize L
        L = zeros(params.n_nodes, params.n_nodes);
        edge_idx = 1;
        
        % Fill base edges (always active)
        for e = 1:n_base_edges
            i = base_edge_list(e, 1);
            j = base_edge_list(e, 2);
            value = sum(coefficients(edge_idx, :) .* basis_functions(f, :));
            L(i, j) = value;
            edge_idx = edge_idx + 1;
        end
        
        % Fill variable edges (with activation)
        for e = 1:size(variable_edge_list, 1)
            i = variable_edge_list(e, 1);
            j = variable_edge_list(e, 2);
            
            % Base value from basis functions
            base_value = sum(coefficients(edge_idx, :) .* basis_functions(f, :));
            
            % Apply smooth activation
            activated_value = base_value * variable_activations(f, e);
            
            % Use soft thresholding instead of hard thresholding
            if abs(activated_value) > 0.01 * params.sigma_coef
                L(i, j) = activated_value;
            end
            edge_idx = edge_idx + 1;
        end
        
        % Fill diagonal
        for i = 1:params.n_nodes
            diag_value = sum(coefficients(edge_idx, :) .* basis_functions(f, :));
            L(i, i) = diag_value;
            edge_idx = edge_idx + 1;
        end
        
        % Enforce diagonal dominance
        for i = 1:params.n_nodes
            row_sum = sum(abs(L(i, 1:i-1)));
            L(i, i) = max(L(i, i), row_sum + params.epsilon_reg);
        end
        
        % Form precision matrix
        Omega = L * L';
        true_precision{f} = Omega;
        true_covariance{f} = inv(Omega);
    end
    
    % Step 7: Generate empirical covariances (same as before)
    emp_covariance = cell(params.n_freq, 1);
    for f = 1:params.n_freq
        Sigma = true_covariance{f};
        Z = mvnrnd(zeros(params.n_nodes, 1), Sigma, params.n_samples)';
        W = Z * Z';
        emp_covariance{f} = W / params.n_samples;
        emp_covariance{f} = (emp_covariance{f} + emp_covariance{f}') / 2;
    end
    
    % Store additional information
    params.base_edge_list = base_edge_list;
    params.variable_edge_list = variable_edge_list;
    params.variable_activations = variable_activations;
    params.n_base_edges = n_base_edges;
    params.n_variable_edges = size(variable_edge_list, 1);
    params.basis_functions = basis_functions;
    params.coefficients = coefficients;
    
    % Analyze actual sparsity pattern variation
    sparsity_changes = analyze_sparsity_changes_improved(true_precision);
    params.sparsity_changes = sparsity_changes;
    
    fprintf('Module 7 (Improved): Simulation data generated\n');
    fprintf('  Base edges (always active): %d\n', n_base_edges);
    fprintf('  Variable edges: %d\n', size(variable_edge_list, 1));
    fprintf('  Sparsity pattern changes: %d\n', sparsity_changes);
end

%% Helper: Analyze sparsity pattern changes
function n_changes = analyze_sparsity_changes_improved(true_precision)
    n_freq = length(true_precision);
    n_nodes = size(true_precision{1}, 1);
    threshold = 1e-10;
    
    n_changes = 0;
    for f = 2:n_freq
        prev_pattern = abs(true_precision{f-1}) > threshold;
        curr_pattern = abs(true_precision{f}) > threshold;
        
        % Count edge changes (upper triangular only)
        for i = 1:n_nodes
            for j = i+1:n_nodes
                if prev_pattern(i,j) ~= curr_pattern(i,j)
                    n_changes = n_changes + 1;
                end
            end
        end
    end
end

%% Keep original helper functions
function edge_list = generateGraphStructure(n_nodes, graph_type, edge_density)
    % Same as original implementation
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
            perm = randperm(size(all_edges, 1));
            selected_idx = perm(1:n_edges);
            edge_list = all_edges(selected_idx, :);
            
        case 'chain'
            for i = 2:n_nodes
                edge_list = [edge_list; i, i-1];
            end
            
        case 'grid'
            grid_size = round(sqrt(n_nodes));
            for row = 1:grid_size
                for col = 1:(grid_size-1)
                    node1 = (row-1)*grid_size + col;
                    node2 = node1 + 1;
                    if node2 <= n_nodes && node1 <= n_nodes
                        edge_list = [edge_list; max(node1,node2), min(node1,node2)];
                    end
                end
            end
            for row = 1:(grid_size-1)
                for col = 1:grid_size
                    node1 = (row-1)*grid_size + col;
                    node2 = node1 + grid_size;
                    if node2 <= n_nodes && node1 <= n_nodes
                        edge_list = [edge_list; max(node1,node2), min(node1,node2)];
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
    % Same as original
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
    
    for b = 1:n_basis
        if sum(basis(:, b)) > 0
            basis(:, b) = basis(:, b) / sum(basis(:, b));
        end
    end
end