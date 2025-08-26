function edge_proxies = module3_edge_proxy_computation(input_matrices, proxy_method, options)
% MODULE3_EDGE_PROXY_COMPUTATION - Compute edge proxy values for active set selection
%
% Syntax:
%   edge_proxies = module3_edge_proxy_computation(input_matrices, proxy_method)
%   edge_proxies = module3_edge_proxy_computation(input_matrices, proxy_method, options)
%
% Description:
%   Computes edge proxy values using either correlation-based or precision-based
%   methods. These proxy values are used to identify potentially active edges
%   in sparse precision matrix estimation.
%
%   Methods:
%   - 'correlation': c_ij(f) = |Σ_ij(f)| (magnitude of whitened covariance)
%   - 'precision': c_ij(f) = |-Ω_ij(f)| / sqrt(Ω_ii(f) * Ω_jj(f)) (partial coherence)
%
% Input Arguments:
%   input_matrices - (cell array, Fx1) Input matrices for proxy computation
%                    For 'correlation': whitened covariance matrices
%                    For 'precision': initial precision matrices
%   proxy_method - ('correlation'|'precision') Method for proxy computation
%
% Name-Value Arguments:
%   regularization_factor - (double, >=0) Small value added to diagonal for stability
%                          Default: 1e-8
%   symmetrize_output     - (logical) Force output symmetry. Default: true
%   zero_diagonal         - (logical) Set diagonal elements to zero. Default: true
%   verbose               - (logical) Display computation details. Default: false
%
% Output Arguments:
%   edge_proxies - (cell array, Fx1) Edge proxy matrices, one per frequency
%                  Each matrix is pxp with proxy values c_ij(f)
%
% Examples:
%   % Correlation-based proxies
%   proxies = module3_edge_proxy_computation(whitened_cov, 'correlation');
%   
%   % Precision-based proxies with options
%   options.regularization_factor = 1e-6;
%   proxies = module3_edge_proxy_computation(init_precision, 'precision', options);
%
% Mathematical Background:
%   Correlation method: c_ij(f) = |Σ̃_ij(f)| where Σ̃ is whitened covariance
%   Precision method: c_ij(f) = |-Ω_ij(f)| / sqrt(Ω_ii(f) * Ω_jj(f))
%   
%   The precision method computes a normalized partial coherence measure
%   that is scale-invariant and bounded between 0 and 1.
%
% See also: MODULE3_ACTIVE_SET_MAIN, MODULE3_THRESHOLD_DETERMINATION
%
% Author: [Your Name]
% Date: [Current Date]
% Version: 1.0

% ==================== Input Validation ====================
if nargin < 2
    error('module3_edge_proxy_computation:insufficient_input', ...
          'input_matrices and proxy_method are required');
end

if nargin < 3
    options = struct();
end

% Validate input_matrices
if ~iscell(input_matrices)
    error('module3_edge_proxy_computation:invalid_input', ...
          'input_matrices must be a cell array');
end

F = numel(input_matrices);
if F == 0
    error('module3_edge_proxy_computation:empty_input', ...
          'input_matrices cannot be empty');
end

% Get dimensions from first matrix
first_matrix = input_matrices{1};
if ~isnumeric(first_matrix) || ~ismatrix(first_matrix)
    error('module3_edge_proxy_computation:invalid_matrix', ...
          'Each input matrix must be numeric');
end

p = size(first_matrix, 1);
if size(first_matrix, 2) ~= p
    error('module3_edge_proxy_computation:not_square', ...
          'Input matrices must be square');
end

% Validate all matrices have consistent dimensions
for f = 1:F
    if ~isnumeric(input_matrices{f}) || ~isequal(size(input_matrices{f}), [p p])
        error('module3_edge_proxy_computation:dimension_mismatch', ...
              'All matrices must be %dx%d, matrix %d has size %dx%d', ...
              p, p, f, size(input_matrices{f}, 1), size(input_matrices{f}, 2));
    end
end

% Validate proxy_method
valid_methods = {'correlation', 'precision'};
if ~ischar(proxy_method) && ~isstring(proxy_method)
    error('module3_edge_proxy_computation:invalid_method_type', ...
          'proxy_method must be a string or char');
end

proxy_method = char(proxy_method);
if ~ismember(proxy_method, valid_methods)
    error('module3_edge_proxy_computation:invalid_method', ...
          'proxy_method must be one of: %s', strjoin(valid_methods, ', '));
end

% ==================== Parameter Setup ====================
defaults = struct();
defaults.regularization_factor = 1e-8;
defaults.symmetrize_output = true;
defaults.zero_diagonal = true;
defaults.verbose = false;

field_names = fieldnames(defaults);
for i = 1:numel(field_names)
    fname = field_names{i};
    if ~isfield(options, fname)
        options.(fname) = defaults.(fname);
    end
end

% Validate parameters
if ~isnumeric(options.regularization_factor) || options.regularization_factor < 0
    error('module3_edge_proxy_computation:invalid_regularization', ...
          'regularization_factor must be non-negative scalar');
end

% ==================== Computation ====================
edge_proxies = cell(F, 1);

if options.verbose
    fprintf('Computing edge proxies using %s method for %d frequencies\n', ...
            proxy_method, F);
end

switch proxy_method
    case 'correlation'
        % Correlation-based: c_ij(f) = |Σ_ij(f)|
        for f = 1:F
            S_f = input_matrices{f};
            
            % Input validation for this frequency
            if any(isnan(S_f(:))) || any(isinf(S_f(:)))
                error('module3_edge_proxy_computation:invalid_values', ...
                      'Matrix %d contains NaN or Inf values', f);
            end
            
            % Compute absolute values
            proxy_matrix = abs(S_f);
            
            % Apply regularization if requested
            if options.regularization_factor > 0
                proxy_matrix = proxy_matrix + options.regularization_factor * eye(p);
            end
            
            % Zero diagonal if requested
            if options.zero_diagonal
                proxy_matrix(1:p+1:end) = 0;
            end
            
            % Symmetrize if requested
            if options.symmetrize_output
                proxy_matrix = (proxy_matrix + proxy_matrix') / 2;
            end
            
            edge_proxies{f} = proxy_matrix;
        end
        
    case 'precision'
        % Precision-based: c_ij(f) = |-Ω_ij(f)| / sqrt(Ω_ii(f) * Ω_jj(f))
        for f = 1:F
            Omega_f = input_matrices{f};
            
            % Input validation for this frequency
            if any(isnan(Omega_f(:))) || any(isinf(Omega_f(:)))
                error('module3_edge_proxy_computation:invalid_values', ...
                      'Matrix %d contains NaN or Inf values', f);
            end
            
            % Check positive definiteness (diagonal elements should be positive)
            diag_elements = diag(Omega_f);
            if any(real(diag_elements) <= 0)
                warning('module3_edge_proxy_computation:non_positive_diagonal', ...
                        'Matrix %d has non-positive diagonal elements, regularizing', f);
                Omega_f = Omega_f + options.regularization_factor * eye(p);
                diag_elements = diag(Omega_f);
            end
            
            % Compute partial coherence proxy
            proxy_matrix = zeros(p, p);
            
            for i = 1:p
                for j = 1:p
                    if i ~= j
                        denominator = sqrt(abs(diag_elements(i) * diag_elements(j)));
                        if denominator > options.regularization_factor
                            proxy_matrix(i, j) = abs(-Omega_f(i, j)) / denominator;
                        else
                            proxy_matrix(i, j) = 0;
                        end
                    end
                end
            end
            
            % Zero diagonal if requested
            if options.zero_diagonal
                proxy_matrix(1:p+1:end) = 0;
            end
            
            % Symmetrize if requested  
            if options.symmetrize_output
                proxy_matrix = (proxy_matrix + proxy_matrix') / 2;
            end
            
            edge_proxies{f} = proxy_matrix;
        end
        
    otherwise
        error('module3_edge_proxy_computation:unknown_method', ...
              'Unknown proxy method: %s', proxy_method);
end

% ==================== Output Validation ====================
for f = 1:F
    proxy_matrix = edge_proxies{f};
    
    % Check for finite values
    if any(~isfinite(proxy_matrix(:)))
        warning('module3_edge_proxy_computation:infinite_values', ...
                'Proxy matrix %d contains non-finite values', f);
    end
    
    % Check for negative values (should not occur for these methods)
    if any(proxy_matrix(:) < 0)
        warning('module3_edge_proxy_computation:negative_values', ...
                'Proxy matrix %d contains negative values', f);
    end
    
    % Check symmetry if requested
    if options.symmetrize_output
        asymmetry = norm(proxy_matrix - proxy_matrix', 'fro');
        if asymmetry > sqrt(eps) * norm(proxy_matrix, 'fro')
            warning('module3_edge_proxy_computation:asymmetry', ...
                    'Proxy matrix %d is not symmetric (error: %.2e)', f, asymmetry);
        end
    end
end

if options.verbose
    % Report statistics
    total_proxies = 0;
    max_proxy = 0;
    min_proxy = inf;
    
    for f = 1:F
        proxy_matrix = edge_proxies{f};
        if options.zero_diagonal
            relevant_values = proxy_matrix(~eye(p));
        else
            relevant_values = proxy_matrix(:);
        end
        
        valid_values = relevant_values(isfinite(relevant_values));
        if ~isempty(valid_values)
            total_proxies = total_proxies + numel(valid_values);
            max_proxy = max(max_proxy, max(valid_values));
            min_proxy = min(min_proxy, min(valid_values));
        end
    end
    
    fprintf('Edge proxy computation completed:\n');
    fprintf('  Method: %s\n', proxy_method);
    fprintf('  Total proxy values: %d\n', total_proxies);
    fprintf('  Range: [%.6f, %.6f]\n', min_proxy, max_proxy);
end

end