function [true_precision, true_covariance, emp_covariance, params] = generateSimulationData(varargin)
% generateSimulationData - Generate synthetic data for testing the precision matrix estimation algorithm
%
% Syntax:
%   [true_precision, true_covariance, emp_covariance, params] = generateSimulationData(...)
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
%   true_precision - Cell array {n_freq x 1} of true precision matrices
%   true_covariance - Cell array {n_freq x 1} of true covariance matrices
%   emp_covariance - Cell array {n_freq x 1} of empirical covariance matrices
%   params - Structure containing all parameters used

    % Parse input arguments
    p = inputParser;
    addParameter(p, 'n_nodes', 20, @(x) x > 0);
    addParameter(p, 'n_freq', 50, @(x) x > 0);
    addParameter(p, 'n_samples', 100, @(x) x > 0);
    addParameter(p, 'graph_type', 'random', @ischar);
    addParameter(p, 'edge_density', 0.1, @(x) x > 0 && x <= 1);
    addParameter(p, 'n_basis', 5, @(x) x > 0);
    addParameter(p, 'sigma_coef', 0.5, @(x) x > 0);
    addParameter(p, 'epsilon_reg', 0.1, @(x) x > 0);
    addParameter(p, 'random_seed', 42, @isnumeric);
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
    % Coefficients for edges (including diagonal)
    n_total_entries = params.n_edges + params.n_nodes; % edges + diagonal
    coefficients = params.sigma_coef * randn(n_total_entries, params.n_basis);
    
    % Step 4: Construct precision matrices for each frequency
    true_precision = cell(params.n_freq, 1);
    true_covariance = cell(params.n_freq, 1);
    
    for f = 1:params.n_freq
        % Initialize lower triangular factor
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
        
        % Form precision matrix
        Omega = L * L';
        true_precision{f} = Omega;
        
        % Compute true covariance
        true_covariance{f} = inv(Omega);
    end
    
    % Step 5: Generate empirical covariances via Wishart distribution
    emp_covariance = cell(params.n_freq, 1);
    
    for f = 1:params.n_freq
        Sigma = true_covariance{f};
        
        % Generate Wishart sample
        % W ~ Wishart(T, Sigma) => W = Z * Z' where Z columns ~ N(0, Sigma)
        Z = mvnrnd(zeros(params.n_nodes, 1), Sigma, params.n_samples)';
        W = Z * Z';
        
        % Empirical covariance
        emp_covariance{f} = W / params.n_samples;
    end
    
    % Add diagnostic information to params
    params.basis_functions = basis_functions;
    params.coefficients = coefficients;
    
    fprintf('Simulation data generated successfully:\n');
    fprintf('  - Nodes: %d\n', params.n_nodes);
    fprintf('  - Frequencies: %d\n', params.n_freq);
    fprintf('  - Samples per frequency: %d\n', params.n_samples);
    fprintf('  - Graph type: %s\n', params.graph_type);
    fprintf('  - Number of edges: %d\n', params.n_edges);
    fprintf('  - Sparsity: %.2f%%\n', 100 * params.n_edges / (params.n_nodes * (params.n_nodes - 1) / 2));
end

function edge_list = generateGraphStructure(n_nodes, graph_type, edge_density)
% Generate graph structure as edge list
% Returns edge_list as [i, j] pairs where i > j

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
            % 2D grid graph (assumes n_nodes is perfect square)
            grid_size = round(sqrt(n_nodes));
            if grid_size^2 ~= n_nodes
                warning('n_nodes is not perfect square, using approximate grid');
            end
            
            % Add horizontal edges
            for row = 1:grid_size
                for col = 1:(grid_size-1)
                    node1 = (row-1)*grid_size + col;
                    node2 = node1 + 1;
                    if node2 <= n_nodes
                        edge_list = [edge_list; max(node1,node2), min(node1,node2)];
                    end
                end
            end
            
            % Add vertical edges
            for row = 1:(grid_size-1)
                for col = 1:grid_size
                    node1 = (row-1)*grid_size + col;
                    node2 = node1 + grid_size;
                    if node2 <= n_nodes
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
    edge_list = sortrows(edge_list);
end

function basis = generateSmoothBasis(freq_normalized, n_basis)
% Generate smooth basis functions (raised cosines)
% freq_normalized: normalized frequency points in [0, 1]
% Returns: [n_freq x n_basis] matrix

    n_freq = length(freq_normalized);
    basis = zeros(n_freq, n_basis);
    
    % Create raised cosine basis functions
    centers = linspace(0, 1, n_basis);
    width = 2.0 / n_basis; % Overlapping support
    
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
        basis(:, b) = basis(:, b) / sum(basis(:, b));
    end
end