function test_module4_gradient()
% TEST_MODULE4_GRADIENT - Unit tests for module4_gradient
%
% Run all tests with: test_module4_gradient()

fprintf('\n=== Testing module4_gradient ===\n');

% Test 1: Basic functionality
test_basic_functionality();

% Test 2: Parameter validation
test_parameter_validation();

% Test 3: Wrapper compatibility
test_wrapper_compatibility();

fprintf('\nAll basic tests passed!\n');
fprintf('For comprehensive testing, run: test_module4_comprehensive()\n');

end

function test_basic_functionality()
fprintf('Test 1: Basic functionality... ');

try
    % Create simple test data
    p = 4;  % 4x4 matrices
    F = 2;  % 2 frequencies
    
    % Create test precision matrices (positive definite)
    Gammas = cell(F, 1);
    for f = 1:F
        A = randn(p, p);
        Gammas{f} = A' * A + eye(p);  % Ensure positive definite
    end
    
    % Create test covariance matrices
    Sigmas = cell(F, 1);
    for f = 1:F
        B = randn(p, p);
        Sigmas{f} = B * B' + 0.1 * eye(p);  % Positive definite
        Sigmas{f} = (Sigmas{f} + Sigmas{f}') / 2;  % Ensure Hermitian
    end
    
    % Create smoothing kernel
    K = 0.2 * ones(F, F) - 0.2 * eye(F);  % Off-diagonal coupling
    
    % Create weight matrix
    W = eye(p);
    
    % Prepare input data
    input_data = struct();
    input_data.precision_matrices = Gammas;
    input_data.whitened_covariances = Sigmas;
    input_data.smoothing_kernel = K;
    input_data.weight_matrix = W;
    
    % Test with default parameters
    results = module4_gradient(input_data);
    
    % Verify results structure
    assert(isstruct(results), 'Results should be a struct');
    assert(isfield(results, 'smooth_gradients'), 'Results should have smooth_gradients field');
    assert(isfield(results, 'success'), 'Results should have success field');
    assert(results.success, 'Computation should succeed');
    assert(length(results.smooth_gradients) == F, 'Should have gradients for all frequencies');
    
    % Check gradient dimensions
    for f = 1:F
        assert(isequal(size(results.smooth_gradients{f}), [p, p]), ...
               'Gradient should have correct dimensions');
        assert(all(isfinite(results.smooth_gradients{f}(:))), ...
               'Gradient should be finite');
    end
    
    fprintf('PASSED\n');
    
catch ME
    fprintf('FAILED: %s\n', ME.message);
    rethrow(ME);
end

end

function test_parameter_validation()
fprintf('Test 2: Parameter validation... ');

try
    % Create minimal valid input
    p = 3;
    F = 2;
    
    input_data = struct();
    input_data.precision_matrices = {eye(p), eye(p)};
    input_data.whitened_covariances = {eye(p), eye(p)};
    input_data.smoothing_kernel = zeros(F, F);
    input_data.weight_matrix = eye(p);
    
    % Test with custom parameters
    params = struct();
    params.lambda1 = 0.05;
    params.lambda2 = 0.02;
    params.verbose = false;
    params.use_graph_laplacian = true;
    
    results = module4_gradient(input_data, params);
    
    assert(results.success, 'Should succeed with valid parameters');
    
    % Test invalid lambda1 (should be handled gracefully)
    params.lambda1 = -0.01;  % Negative
    
    try
        module4_gradient(input_data, params);
        error('Should have failed with negative lambda1');
    catch ME
        % Expected to fail
        assert(contains(ME.message, 'lambda1') || contains(ME.message, 'negative'), ...
               'Should give appropriate error for negative lambda1');
    end
    
    fprintf('PASSED\n');
    
catch ME
    fprintf('FAILED: %s\n', ME.message);
    rethrow(ME);
end

end

function test_wrapper_compatibility()
fprintf('Test 3: Wrapper compatibility... ');

try
    % Create test data
    p = 3;
    F = 2;
    
    input_data = struct();
    input_data.precision_matrices = {2*eye(p), 1.5*eye(p)};
    input_data.whitened_covariances = {0.6*eye(p), 0.8*eye(p)};
    input_data.smoothing_kernel = 0.1 * ones(F, F);
    input_data.weight_matrix = eye(p);
    
    % Test wrapper function
    results_wrapper = module4_gradient(input_data);
    
    % Test main implementation directly
    results_main = module4_objective_gradient_main(input_data, struct());
    
    % Compare results
    assert(isequal(size(results_wrapper.smooth_gradients), size(results_main.smooth_gradients)), ...
           'Wrapper and main should return same structure');
    
    for f = 1:F
        diff = norm(results_wrapper.smooth_gradients{f} - results_main.smooth_gradients{f}, 'fro');
        assert(diff < 1e-12, 'Wrapper and main should give identical results');
    end
    
    fprintf('PASSED\n');
    
catch ME
    fprintf('FAILED: %s\n', ME.message);
    rethrow(ME);
end

end