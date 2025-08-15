function [true_precision, true_covariance, emp_covariance, params] = module7_simulation(varargin)
% MODULE7_SIMULATION - Generate synthetic data for testing the precision matrix estimation algorithm
%
% This module implements the Wishart simulation method described in the algorithm document.
% It generates sparse, smoothly-varying precision matrices across frequencies.
%
% Syntax:
%   [true_precision, true_covariance, emp_covariance, params] = module7_simulation(...)
%
% Inputs (name-value pairs):
%   'n_nodes'      - Number of nodes/channels (default: 20)
%   'n_freq'       - Number of frequency points (default: 50)
%   'n_samples'    - Number of samples per frequency (default: 100)
%   'graph_type'   - Type of graph: 'random', 'chain', 'grid', 'hub' (default: 'random')
%   'edge_density' - Edge density for random graph (default: 0.1)
%   'n_basis'      - Number of smooth basis functions (default: 5)
%   'sigma_coef'   - Std of basis coefficients (default: 0.5)
%   'epsilon_reg'  - Diagonal dominance parameter (default: 0.1)
%   'random_seed'  - Random seed for reproducibility (default: 42)
%
% Outputs:
%   true_precision  - Cell array {n_freq x 1} of true precision matrices
%   true_covariance - Cell array {n_freq x 1} of true covariance matrices
%   emp_covariance  - Cell array {n_freq x 1} of empirical covariance matrices
%   params          - Structure containing all parameters used
%
% Example:
%   % Generate data with chain graph structure
%   [Omega, Sigma, Sigma_emp, params] = module7_simulation(...
%       'n_nodes', 10, 'graph_type', 'chain', 'n_freq', 30);
%
% See also: generateGraphStructure, generateSmoothBasis

    % Parse input arguments with simplified validation
    p = inputParser;
    % Use simpler validation functions to avoid issues
    addParameter(p, 'n_nodes', 20, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
    addParameter(p, 'n_freq', 50, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
    addParameter(p, 'n_samples', 100, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
    addParameter(p, 'graph_type', 'random', @(x) any(validatestring(x, {'random', 'chain', 'grid', 'hub'})));
    addParameter(p, 'edge_density', 0.1, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', '<=', 1}));
    addParameter(p, 'n_basis', 5, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
    addParameter(p, 'sigma_coef', 0.5, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
    addParameter(p, 'epsilon_reg', 0.1, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
    addParameter(p, 'random_seed', 42, @(x) validateattributes(x, {'numeric'}, {'scalar'}));
    parse(p, varargin{:});
    params = p.Results;
    
    % Set random seed for reproducibility
    rng(params.random_seed);
    
    % Step 1: Generate graph structure
    edge_list = generateGraphStructure(params.n_nodes, params.graph_type, params.edge_density);
    params.edge_list = edge_list;
    params.n_edges = size(edge_list, 1);
    
    % Step 2: Generate smooth basis functions
    freq_normalized = linspace(0, 1, params.n_freq)';
    basis_functions = generateSmoothBasis(freq_normalized, params.n_basis);
    
    % Step 3: Generate basis coefficients for each edge
    % Total entries = edges + diagonal elements
    n_total_entries = params.n_edges + params.n_nodes;
    coefficients = params.sigma_coef * randn(n_total_entries, params.n_basis);
    
    % Step 4: Construct precision matrices for each frequency
    true_precision = cell(params.n_freq, 1);
    true_covariance = cell(params.n_freq, 1);
    
    for f = 1:params.n_freq
        % Initialize lower triangular factor L
        L = zeros(params.n_nodes, params.n_nodes);
        
        % Fill in edge values using basis functions
        edge_idx = 1;
        for e = 1:params.n_edges
            i = edge_list(e, 1);
            j = edge_list(e, 2);
            % Compute value as weighted sum of basis functions
            value = sum(coefficients(edge_idx, :) .* basis_functions(f, :));
            L(i, j) = value;
            edge_idx = edge_idx + 1;
        end
        
        % Fill in diagonal values
        for i = 1:params.n_nodes
            diag_value = sum(coefficients(edge_idx, :) .* basis_functions(f, :));
            L(i, i) = diag_value;
            edge_idx = edge_idx + 1;
        end
        
        % Enforce diagonal dominance for positive definiteness
        for i = 1:params.n_nodes
            row_sum = sum(abs(L(i, 1:i-1)));
            L(i, i) = max(L(i, i), row_sum + params.epsilon_reg);
        end
        
        % Form precision matrix: Omega = L * L'
        Omega = L * L';
        true_precision{f} = Omega;
        
        % Compute true covariance: Sigma = Omega^(-1)
        true_covariance{f} = inv(Omega);
    end
    
    % Step 5: Generate empirical covariances via Wishart distribution
    emp_covariance = cell(params.n_freq, 1);
    
    for f = 1:params.n_freq
        Sigma = true_covariance{f};
        
        % Generate Wishart sample: W ~ Wishart(T, Sigma)
        % Method: W = Z * Z' where columns of Z ~ N(0, Sigma)
        Z = mvnrnd(zeros(params.n_nodes, 1), Sigma, params.n_samples)';
        W = Z * Z';
        
        % Empirical covariance = W / T
        emp_covariance{f} = W / params.n_samples;
        
        % Ensure numerical symmetry
        emp_covariance{f} = (emp_covariance{f} + emp_covariance{f}') / 2;
    end
    
    % Add diagnostic information to params
    params.basis_functions = basis_functions;
    params.coefficients = coefficients;
    params.creation_date = datestr(now);
    
    % Calculate sparsity pattern changes
    params.sparsity_changes = calculate_sparsity_changes(true_precision);
    
    % Print summary
    fprintf('Module 7: Simulation data generated successfully\n');
    fprintf('  - Nodes: %d\n', params.n_nodes);
    fprintf('  - Frequencies: %d\n', params.n_freq);
    fprintf('  - Samples per frequency: %d\n', params.n_samples);
    fprintf('  - Graph type: %s\n', params.graph_type);
    fprintf('  - Number of edges: %d\n', params.n_edges);
    fprintf('  - Sparsity: %.2f%%\n', 100 * params.n_edges / (params.n_nodes * (params.n_nodes - 1) / 2));
    fprintf('  - Sparsity pattern changes: %d\n', params.sparsity_changes);
end

%% Helper function to calculate sparsity changes
function n_changes = calculate_sparsity_changes(true_precision)
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

%% Local function: Generate graph structure
function edge_list = generateGraphStructure(n_nodes, graph_type, edge_density)
% Generate graph structure as edge list
% Returns edge_list as [i, j] pairs where i > j (lower triangular)

    edge_list = [];
    
    switch lower(graph_type)
        case 'random'
            % Random graph with specified density
            n_possible_edges = n_nodes * (n_nodes - 1) / 2;
            n_edges = round(edge_density * n_possible_edges);
            
            % Generate all possible edges
            all_edges = [];
            for i = 2:n_nodes
                for j = 1:(i-1)
                    all_edges = [all_edges; i, j];
                end
            end
            
            % Randomly select edges
            perm = randperm(size(all_edges, 1));
            selected_idx = perm(1:n_edges);
            edge_list = all_edges(selected_idx, :);
            
        case 'chain'
            % Chain graph: 1-2-3-...-n
            for i = 2:n_nodes
                edge_list = [edge_list; i, i-1];
            end
            
        case 'grid'
            % 2D grid graph (assumes n_nodes is approximately square)
            grid_size = round(sqrt(n_nodes));
            if grid_size^2 ~= n_nodes
                warning('n_nodes is not perfect square, using %d x %d grid', grid_size, grid_size);
            end
            
            % Add horizontal edges
            for row = 1:grid_size
                for col = 1:(grid_size-1)
                    node1 = (row-1)*grid_size + col;
                    node2 = node1 + 1;
                    if node2 <= n_nodes && node1 <= n_nodes
                        edge_list = [edge_list; max(node1,node2), min(node1,node2)];
                    end
                end
            end
            
            % Add vertical edges
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
            % Hub graph: node 1 connected to all others
            for i = 2:n_nodes
                edge_list = [edge_list; i, 1];
            end
            
        otherwise
            error('Unknown graph type: %s', graph_type);
    end
    
    % Sort edge list for consistency
    edge_list = unique(edge_list, 'rows');
end

%% Local function: Generate smooth basis functions
function basis = generateSmoothBasis(freq_normalized, n_basis)
% Generate smooth basis functions using raised cosines
% freq_normalized: normalized frequency points in [0, 1]
% Returns: [n_freq x n_basis] matrix

    n_freq = length(freq_normalized);
    basis = zeros(n_freq, n_basis);
    
    % Create raised cosine basis functions with overlapping support
    centers = linspace(0, 1, n_basis);
    width = 2.0 / n_basis; % Ensures overlap between adjacent basis functions
    
    for b = 1:n_basis
        for f = 1:n_freq
            dist = abs(freq_normalized(f) - centers(b));
            if dist < width/2
                % Raised cosine formula
                basis(f, b) = 0.5 * (1 + cos(2*pi*dist/width));
            end
        end
    end
    
    % Normalize each basis function to sum to 1
    for b = 1:n_basis
        if sum(basis(:, b)) > 0
            basis(:, b) = basis(:, b) / sum(basis(:, b));
        end
    end
end