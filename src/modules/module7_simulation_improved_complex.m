function [true_precision, true_covariance, emp_covariance, params] = module7_simulation_improved_complex(varargin)
% MODULE7_SIMULATION_IMPROVED_COMPLEX - Enhanced complex Hermitian simulation
% 
% MAJOR UPDATE: Fixed to ensure frequent small changes between adjacent frequencies
% 
% Key improvements:
% 1. More variable edges selected for dynamic patterns
% 2. Multi-layer sinusoidal activation ensuring every frequency has changes
% 3. Neighborhood correlation for smooth but frequent transitions
% 4. Post-processing to guarantee minimum change levels

%% Parse input parameters (same as before)
p = inputParser;
addParameter(p, 'n_nodes', 8, @(x) isscalar(x) && x > 2);
addParameter(p, 'n_freq', 12, @(x) isscalar(x) && x > 1);
addParameter(p, 'n_samples', 100, @(x) isscalar(x) && x > 0);
addParameter(p, 'graph_type', 'random', @(x) ismember(x, {'random', 'chain', 'hub'}));
addParameter(p, 'edge_density', 0.4, @(x) isscalar(x) && x > 0 && x < 1);
addParameter(p, 'sparsity_variation', 0.3, @(x) isscalar(x) && x >= 0 && x <= 1);
addParameter(p, 'edge_activation_smoothness', 0.8, @(x) isscalar(x) && x > 0);
addParameter(p, 'n_basis', 4, @(x) isscalar(x) && x > 0);
addParameter(p, 'sigma_coef', 0.5, @(x) isscalar(x) && x > 0);
addParameter(p, 'complex_strength', 1.0, @(x) isscalar(x) && x >= 0 && x <= 2);
addParameter(p, 'epsilon_reg', 0.1, @(x) isscalar(x) && x > 0);
addParameter(p, 'random_seed', [], @(x) isempty(x) || (isscalar(x) && x >= 0));
addParameter(p, 'coefficient_complex_fraction', 1.0, @(x) isscalar(x) && x >= 0 && x <= 1);

parse(p, varargin{:});
params = p.Results;

if ~isempty(params.random_seed)
    rng(params.random_seed);
end

fprintf('=== Enhanced Complex Simulation with Dynamic Sparsity ===\n');

%% Step 1: Generate base graph structure (same as before)
base_edge_list = generateGraphStructure(params.n_nodes, params.graph_type, params.edge_density);
n_base_edges = size(base_edge_list, 1);
fprintf('Generated %d base edges\n', n_base_edges);

%% Step 2: Generate MORE variable edges for frequent changes - ENHANCED
all_possible_edges = [];
for i = 2:params.n_nodes
    for j = 1:(i-1)
        all_possible_edges = [all_possible_edges; i, j];
    end
end

% Remove base edges from possible variable edges
base_edge_set = [base_edge_list; base_edge_list(:, [2,1])];
is_base = false(size(all_possible_edges, 1), 1);
for e = 1:size(all_possible_edges, 1)
    is_base(e) = any(all(base_edge_set == all_possible_edges(e, :), 2));
end
variable_edges = all_possible_edges(~is_base, :);

% ENHANCED: Select MORE variable edges for frequent changes
% Increase the number significantly to allow more variation
min_variable_edges = max(3, round(n_base_edges * 0.4)); % At least 40% of base edges
n_variable = round(params.sparsity_variation * n_base_edges) + min_variable_edges;

if n_variable > 0 && size(variable_edges, 1) > 0
    n_variable = min(n_variable, size(variable_edges, 1)); % Don't exceed available edges
    perm = randperm(size(variable_edges, 1));
    selected_idx = perm(1:n_variable);
    variable_edge_list = variable_edges(selected_idx, :);
else
    variable_edge_list = [];
end

n_variable_edges = size(variable_edge_list, 1);
fprintf('Selected %d variable edges for dynamic sparsity\n', n_variable_edges);

%% Step 3: Generate smooth basis functions (same as before)
freq_normalized = linspace(0, 1, params.n_freq)';
basis_functions = generateSmoothBasis(freq_normalized, params.n_basis);

%% Step 4: Generate complex coefficients (same as before)
fprintf('Coefficient complex fraction: %.2f\n', params.coefficient_complex_fraction);

n_total_edges = n_base_edges + n_variable_edges;
n_total_coefficients = n_total_edges + params.n_nodes;

coefficients = zeros(n_total_coefficients, params.n_basis);

for c = 1:n_total_coefficients
    for b = 1:params.n_basis
        is_complex = rand() < params.coefficient_complex_fraction;
        
        if is_complex
            real_part = randn() * params.sigma_coef;
            imag_part = randn() * params.sigma_coef * params.complex_strength;
            coefficients(c, b) = complex(real_part, imag_part);
        else
            coefficients(c, b) = randn() * params.sigma_coef;
        end
    end
end

%% Step 5: Generate variable edge activation patterns - MAJOR ENHANCEMENT FOR FREQUENT CHANGES
variable_activations = ones(params.n_freq, n_variable_edges);

if n_variable_edges > 0
    fprintf('Designing frequent change patterns for %d variable edges\n', n_variable_edges);
    
    for e = 1:n_variable_edges
        % Strategy: Multi-layer activation for frequent but smooth changes
        freq_indices = (0:params.n_freq-1)' / max(1, params.n_freq-1);
        
        % Layer 1: Slow background variation (creates overall structure)
        slow_freq = 0.8 + 0.4 * rand(); % Random frequency [0.8, 1.2]
        slow_phase = rand() * 2 * pi;
        slow_component = sin(2 * pi * freq_indices * slow_freq + slow_phase);
        
        % Layer 2: Medium frequency changes (creates groups of similar frequencies)  
        medium_freq = 2.5 + 1.5 * rand(); % Random frequency [2.5, 4.0]
        medium_phase = rand() * 2 * pi;
        medium_component = 0.7 * sin(2 * pi * freq_indices * medium_freq + medium_phase);
        
        % Layer 3: High frequency changes (ensures adjacent frequencies differ)
        high_freq = 6 + 6 * rand(); % Random frequency [6, 12] 
        high_phase = rand() * 2 * pi;
        high_component = 0.5 * sin(2 * pi * freq_indices * high_freq + high_phase);
        
        % Layer 4: Very high frequency for fine-scale changes
        very_high_freq = 15 + 10 * rand(); % Random frequency [15, 25]
        very_high_phase = rand() * 2 * pi;
        very_high_component = 0.3 * sin(2 * pi * freq_indices * very_high_freq + very_high_phase);
        
        % Combine all layers
        combined_signal = slow_component + medium_component + high_component + very_high_component;
        
        % Apply smoothness parameter to control how sharp the transitions are
        smoothness_factor = params.edge_activation_smoothness * 1.5;
        activation_probability = 1 ./ (1 + exp(-smoothness_factor * combined_signal));
        
        % Add neighborhood smoothing to ensure gradual changes between adjacent frequencies
        smoothed_probability = activation_probability;
        alpha = 0.3; % Smoothing strength
        for f = 2:params.n_freq-1
            smoothed_probability(f) = (1-alpha) * activation_probability(f) + ...
                                    alpha * 0.5 * (activation_probability(f-1) + activation_probability(f+1));
        end
        
        % Convert to binary with dynamic thresholds
        base_threshold = 0.4 + 0.3 * rand(); % Random base threshold [0.4, 0.7]
        
        % Add threshold variation to create more diverse patterns
        threshold_variation = 0.2 * sin(2 * pi * freq_indices * 4 + rand() * 2 * pi);
        dynamic_threshold = base_threshold + threshold_variation;
        
        % Apply thresholding with small random perturbations
        activation = zeros(params.n_freq, 1);
        for f = 1:params.n_freq
            noise = 0.03 * randn(); % Small random noise
            if smoothed_probability(f) + noise > dynamic_threshold(f)
                activation(f) = 1;
            else
                activation(f) = 0;
            end
        end
        
        % Post-processing: Ensure sufficient changes
        n_changes = sum(abs(diff(activation)));
        min_required_changes = max(2, round(params.n_freq * 0.25)); % At least 25% of frequencies should have changes
        
        if n_changes < min_required_changes
            % Force additional changes by flipping some states
            additional_changes_needed = min_required_changes - n_changes;
            change_positions = randperm(params.n_freq-1, additional_changes_needed);
            
            for pos = change_positions
                activation(pos+1) = 1 - activation(pos);
            end
        end
        
        % Final check: Ensure not all frequencies are the same
        if all(activation == activation(1))
            % If somehow all are the same, create a more diverse pattern
            n_to_flip = max(3, round(params.n_freq * 0.4));
            flip_indices = randperm(params.n_freq, n_to_flip);
            activation(flip_indices) = 1 - activation(1);
        end
        
        variable_activations(:, e) = activation;
        
        % Report change statistics for this edge
        edge_changes = sum(abs(diff(activation)));
        active_frequencies = sum(activation);
        fprintf('  Edge %d: %d changes, active in %d/%d frequencies\n', ...
                e, edge_changes, active_frequencies, params.n_freq);
    end
    
    % Report overall statistics
    total_changes = 0;
    for e = 1:n_variable_edges
        total_changes = total_changes + sum(abs(diff(variable_activations(:, e))));
    end
    fprintf('Total activation changes across all variable edges: %d\n', total_changes);
    
    % Ensure we have changes in the sparsity pattern
    if total_changes == 0
        fprintf('WARNING: No activation changes detected, forcing changes...\n');
        % Force at least one edge to change
        for e = 1:min(2, n_variable_edges)
            mid_point = round(params.n_freq / 2);
            variable_activations(mid_point:end, e) = 1 - variable_activations(1, e);
        end
    end
end

%% Step 6: Generate precision matrices (same structure, but will now have more changes)
true_precision = cell(params.n_freq, 1);
true_covariance = cell(params.n_freq, 1);

% Track changes for analysis
sparsity_changes = 0;
prev_pattern = [];
matrices_with_complex = 0;
total_complex_fraction = 0;
max_imag_component = 0;
total_phase_magnitude = 0;
hermitian_count = 0;

for f = 1:params.n_freq
    % Print progress for some frequencies
    if mod(f-1, round(params.n_freq/6)) == 0
        fprintf('Processing frequency %d/%d\n', f, params.n_freq);
    end
    
    % Initialize Cholesky factor
    L = zeros(params.n_nodes);
    edge_idx = 1;
    
    % Fill base edges (always active)
    for e = 1:n_base_edges
        i = base_edge_list(e, 1);
        j = base_edge_list(e, 2);
        
        if i < j
            temp = i; i = j; j = temp;
        end
        
        base_value = sum(coefficients(edge_idx, :) .* basis_functions(f, :));
        L(i, j) = base_value;
        edge_idx = edge_idx + 1;
    end
    
    % Fill variable edges (dynamic activation)
    for e = 1:n_variable_edges
        i = variable_edge_list(e, 1);
        j = variable_edge_list(e, 2);
        
        if i < j
            temp = i; i = j; j = temp;
        end
        
        base_value = sum(coefficients(edge_idx, :) .* basis_functions(f, :));
        activated_value = base_value * variable_activations(f, e);
        
        threshold = 0.05 * params.sigma_coef;
        if abs(activated_value) > threshold
            L(i, j) = activated_value;
        end
        edge_idx = edge_idx + 1;
    end
    
    % Fill diagonal elements
    for i = 1:params.n_nodes
        diag_value = sum(coefficients(edge_idx, :) .* basis_functions(f, :));
        L(i, i) = diag_value;
        edge_idx = edge_idx + 1;
    end
    
    % Ensure diagonal dominance for positive definiteness
    for i = 1:params.n_nodes
        off_diag_sum = sum(abs(L(i, 1:i-1))) + sum(abs(L(i+1:end, i)));
        min_diag = off_diag_sum + params.epsilon_reg;
        
        if abs(L(i, i)) < min_diag
            if abs(L(i, i)) > 1e-12
                phase_factor = L(i, i) / abs(L(i, i));
                L(i, i) = min_diag * phase_factor;
            else
                L(i, i) = min_diag;
            end
        end
    end
    
    % Form precision matrix
    Omega = L * L';
    Omega = (Omega + Omega') / 2; % Force Hermitian
    
    true_precision{f} = Omega;
    true_covariance{f} = inv(Omega);
    
    % Track properties
    if any(abs(imag(Omega(:))) > 1e-12)
        matrices_with_complex = matrices_with_complex + 1;
        complex_fraction = sum(abs(imag(Omega(:))) > 1e-12) / numel(Omega);
        total_complex_fraction = total_complex_fraction + complex_fraction;
        
        max_imag = max(abs(imag(Omega(:))));
        max_imag_component = max(max_imag_component, max_imag);
        
        complex_elements = Omega(abs(imag(Omega)) > 1e-12);
        if ~isempty(complex_elements)
            phases = angle(complex_elements);
            avg_phase = mean(abs(phases));
            total_phase_magnitude = total_phase_magnitude + avg_phase;
        end
        
        fprintf('  Freq %d: Complex fraction = %.2f, Max imag = %.4f\n', ...
                f, complex_fraction, max_imag);
    end
    
    if ishermitian(Omega)
        hermitian_count = hermitian_count + 1;
    end
    
    % Track sparsity pattern changes
    current_pattern = abs(triu(Omega, 1)) > 1e-10;
    if ~isempty(prev_pattern) && ~isequal(current_pattern, prev_pattern)
        sparsity_changes = sparsity_changes + 1;
    end
    prev_pattern = current_pattern;
end

%% Step 7: Generate empirical covariances (same as before)
emp_covariance = cell(params.n_freq, 1);
for f = 1:params.n_freq
    Sigma = true_covariance{f};
    
    if any(abs(imag(Sigma(:))) > 1e-12)
        Z_real = mvnrnd(zeros(params.n_nodes, 1), real(Sigma), params.n_samples)';
        Z_imag = mvnrnd(zeros(params.n_nodes, 1), real(Sigma), params.n_samples)';
        Z = Z_real + 1i * Z_imag;
    else
        Z = mvnrnd(zeros(params.n_nodes, 1), real(Sigma), params.n_samples)';
    end
    
    W = Z * Z';
    emp_covariance{f} = W / params.n_samples;
    emp_covariance{f} = (emp_covariance{f} + emp_covariance{f}') / 2;
end

%% Step 8: Store results and analysis
fprintf('\n=== Final Analysis ===\n');
fprintf('Base edges (always active): %d\n', n_base_edges);
fprintf('Variable edges: %d\n', n_variable_edges);
fprintf('Sparsity pattern changes: %d\n', sparsity_changes);
fprintf('Matrices with complex entries: %d/%d\n', matrices_with_complex, params.n_freq);

if matrices_with_complex > 0
    avg_complex_fraction = total_complex_fraction / matrices_with_complex;
    avg_phase_magnitude = total_phase_magnitude / matrices_with_complex;
    fprintf('Average complex fraction: %.2f\n', avg_complex_fraction);
    fprintf('Average phase magnitude: %.4f radians\n', avg_phase_magnitude);
else
    avg_complex_fraction = 0;
    avg_phase_magnitude = 0;
end

fprintf('Max imaginary component: %.4f\n', max_imag_component);
fprintf('Hermitian matrices: %d/%d\n', hermitian_count, params.n_freq);

% Store parameters
params.base_edge_list = base_edge_list;
params.variable_edge_list = variable_edge_list;
params.variable_activations = variable_activations;
params.n_base_edges = n_base_edges;
params.n_variable_edges = n_variable_edges;
params.basis_functions = basis_functions;
params.coefficients = coefficients;
params.sparsity_changes = sparsity_changes;
params.matrices_with_complex = matrices_with_complex;
params.avg_complex_fraction = avg_complex_fraction;
params.max_imag_component = max_imag_component;
params.avg_phase_magnitude = avg_phase_magnitude;
params.hermitian_count = hermitian_count;

fprintf('=== Enhanced Complex Simulation Complete ===\n');

end

% 辅助函数保持不变
function edge_list = generateGraphStructure(n_nodes, graph_type, edge_density)
% Same as before...
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
        
        if n_edges > 0 && n_edges <= size(all_edges, 1)
            perm = randperm(size(all_edges, 1));
            edge_list = all_edges(perm(1:n_edges), :);
        else
            edge_list = [];
        end
        
    case 'chain'
        edge_list = [];
        for i = 2:n_nodes
            edge_list = [edge_list; i, i-1];
        end
        
        remaining_density = edge_density - (n_nodes-1)/(n_nodes*(n_nodes-1)/2);
        if remaining_density > 0
            all_other_edges = [];
            for i = 3:n_nodes
                for j = 1:(i-2)
                    all_other_edges = [all_other_edges; i, j];
                end
            end
            
            if ~isempty(all_other_edges)
                n_additional = round(remaining_density * n_nodes * (n_nodes-1) / 2);
                n_additional = min(n_additional, size(all_other_edges, 1));
                if n_additional > 0
                    perm = randperm(size(all_other_edges, 1));
                    edge_list = [edge_list; all_other_edges(perm(1:n_additional), :)];
                end
            end
        end
        
    case 'hub'
        hub_node = ceil(n_nodes / 2);
        edge_list = [];
        
        for i = 1:n_nodes
            if i ~= hub_node
                if i > hub_node
                    edge_list = [edge_list; i, hub_node];
                else
                    edge_list = [edge_list; hub_node, i];
                end
            end
        end
        
        remaining_density = edge_density - (n_nodes-1)/(n_nodes*(n_nodes-1)/2);
        if remaining_density > 0
            all_other_edges = [];
            for i = 2:n_nodes
                for j = 1:(i-1)
                    if ~((i == hub_node && j ~= hub_node) || (j == hub_node && i ~= hub_node))
                        all_other_edges = [all_other_edges; i, j];
                    end
                end
            end
            
            if ~isempty(all_other_edges)
                n_additional = round(remaining_density * n_nodes * (n_nodes-1) / 2);
                n_additional = min(n_additional, size(all_other_edges, 1));
                if n_additional > 0
                    perm = randperm(size(all_other_edges, 1));
                    edge_list = [edge_list; all_other_edges(perm(1:n_additional), :)];
                end
            end
        end
        
    otherwise
        error('Unknown graph_type: %s. Use ''random'', ''chain'', or ''hub''.', graph_type);
end

if isempty(edge_list)
    edge_list = [2, 1];
end
end

function basis_functions = generateSmoothBasis(freq_normalized, n_basis)
% GENERATESMOOTHBASIS - Generate smooth basis functions for frequency evolution

basis_functions = zeros(length(freq_normalized), n_basis);

for b = 1:n_basis
    if b == 1
        % Constant term
        basis_functions(:, b) = ones(size(freq_normalized));
    elseif b == 2
        % Linear term
        basis_functions(:, b) = freq_normalized;
    else
        % Sinusoidal terms with different frequencies
        freq_mult = floor((b-2)/2) + 1;
        if mod(b, 2) == 1
            % Sine component
            basis_functions(:, b) = sin(freq_mult * 2 * pi * freq_normalized);
        else
            % Cosine component
            basis_functions(:, b) = cos(freq_mult * 2 * pi * freq_normalized);
        end
    end
end

end