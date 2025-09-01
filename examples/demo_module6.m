%% Module 6 Demonstration: Hyperparameter Configuration
% Demonstrates λ₁ and α computed via safe bounds (no changes to Module 7).

clear; clc;
fprintf('========================================\n');
fprintf('Module 6: Hyperparameter Configuration Demo\n');
fprintf('========================================\n\n');

%% Step 1: Basic hyperparameter computation
fprintf('Step 1: Basic hyperparameter computation\n');
fprintf('----------------------------------------\n');

n_nodes = 6; n_freq = 4;
test_input = create_demo_input(n_nodes, n_freq);

tic; cfg = module6_hyperparameters(test_input); t = toc;
fprintf('Computed hyperparameters:\n');
fprintf('  λ₁ = %.6e\n', cfg.lambda1);
fprintf('  α  = %.6e\n', cfg.alpha);
fprintf('  λ₂ = %.6e\n', cfg.lambda2_suggested);
fprintf('  Time = %.3f s\n', t);

%% Step 2: Safety margin sensitivity
fprintf('\nStep 2: Safety margin sensitivity analysis\n');
fprintf('-----------------------------------------\n');
margins = [0.5, 0.7, 0.9, 0.95];
l1 = zeros(size(margins)); a = zeros(size(margins));
fprintf('Safety Margin | λ₁        | α         \n');
fprintf('-------------|-----------|----------\n');
for i=1:numel(margins)
    c = module6_hyperparameters(test_input,'safety_margin',margins(i),'verbose',false);
    l1(i)=c.lambda1; a(i)=c.alpha;
    fprintf('    %.2f     | %.3e | %.3e\n', margins(i), l1(i), a(i));
end
fprintf('\nMonotonicity:\n');
fprintf('  λ₁ decreases with δ: %s\n', tern(all(diff(l1)<=0),'YES','NO'));
fprintf('  α  decreases with δ: %s\n', tern(all(diff(a )<=0),'YES','NO'));

%% Step 3: Scaling with problem size
fprintf('\nStep 3: Problem size scaling behavior\n');
fprintf('------------------------------------\n');
sizes = [4,6,8,10]; times = zeros(size(sizes));
fprintf('Problem Size | λ₁        | α         | Time (s)\n');
fprintf('------------|-----------|-----------|--------\n');
for i=1:numel(sizes)
    p = sizes(i);
    in = create_demo_input(p,3);
    t0=tic; c = module6_hyperparameters(in,'verbose',false); times(i)=toc(t0);
    fprintf('     %2d     | %.3e | %.3e | %.4f\n', p, c.lambda1, c.alpha, times(i));
end
fprintf('\nTime increase 4→10 nodes: %.2fx\n', times(end)/max(1e-12,times(1)));

%% Step 4: Gershgorin vs exact
fprintf('\nStep 4: Method comparison (Gershgorin vs exact)\n');
fprintf('----------------------------------------------\n');
cmp_in = create_demo_input(8,4);
t0=tic; ce = module6_hyperparameters(cmp_in,'use_gershgorin',false,'verbose',false); te=toc(t0);
t0=tic; cg = module6_hyperparameters(cmp_in,'use_gershgorin',true, 'verbose',false); tg=toc(t0);
Le = ce.diagnostics.L_logdet; Lg = cg.diagnostics.L_logdet; ratio = Lg/max(1e-12,Le);
fprintf('                  | Exact      | Gershgorin | Ratio\n');
fprintf('------------------|------------|------------|------\n');
fprintf('L_logdet          | %.3e | %.3e | %.2f\n', Le, Lg, ratio);
fprintf('λ₁                | %.3e | %.3e | %.2f\n', ce.lambda1, cg.lambda1, cg.lambda1/max(1e-12,ce.lambda1));
fprintf('α                 | %.3e | %.3e | %.2f\n', ce.alpha  , cg.alpha  , cg.alpha  /max(1e-12,ce.alpha));
fprintf('Time              | %.4f s  | %.4f s  | %.2fx\n', te, tg, max(1e-12,tg)/max(1e-12,te));

%% Step 5: Integration with Module 7 (no changes to module7)
fprintf('\nStep 5: Integration with Module 7 (no changes)\n');
fprintf('----------------------------------------------\n');
try
    % NOTE: do NOT pass 'verbose' to module7 to keep its original signature
    [~, ~, Sigma_emp, ~] = module7_simulation_improved_complex('n_nodes',6,'n_freq',4,'n_samples',80);
    integ_in = struct();
    integ_in.whitened_covariances = Sigma_emp;
    integ_in.kernel_matrix = create_rbf_kernel(4, 0.8);
    integ_in.weight_matrix = create_correlation_weights(6);

    c = module6_hyperparameters(integ_in,'verbose',false);
    fprintf('Module 7 integration successful:\n');
    fprintf('  λ₁ = %.6e\n', c.lambda1);
    fprintf('  α  = %.6e\n', c.alpha);
    fprintf('  λ₂ = %.6e\n', c.lambda2_suggested);
catch ME
    fprintf('Module 7 integration failed: %s\n', ME.message);
    fprintf('Using synthetic data instead...\n');
    c = module6_hyperparameters(create_demo_input(6,4),'verbose',false);
    fprintf('Synthetic results: λ₁=%.3e, α=%.3e\n', c.lambda1, c.alpha);
end

%% Final summary
final = module6_hyperparameters(create_demo_input(6,4),'verbose',false);
fprintf('\nTypical parameters (6 nodes, 4 freqs): λ₁=%.3e, α=%.3e, λ₂=%.3e\n', ...
    final.lambda1, final.alpha, final.lambda2_suggested);
fprintf('\n=== Module 6 Demo Complete ===\n');

%% Local helpers
function input_struct = create_demo_input(p, F)
% Create SPD covariances; zero-diagonal weights to match theory.
input_struct = struct();
input_struct.whitened_covariances = cell(F,1);
for f=1:F
    A = randn(p); A = A + A';
    m = min(eig(A));
    if m<=0, A = A + (abs(m)+0.1)*eye(p); end
    input_struct.whitened_covariances{f} = (A + A')/2;
end
input_struct.kernel_matrix = create_rbf_kernel(F, 1.0);
input_struct.weight_matrix = create_correlation_weights(p);
end

function K = create_rbf_kernel(F, bw)
x = linspace(0,1,F); K = zeros(F);
for i=1:F, for j=1:F
    K(i,j) = exp(-0.5*((x(i)-x(j))/bw)^2);
end, end
K = K ./ sum(K,2);
end

function W = create_correlation_weights(p)
B = randn(p); B = B + B';
W = abs(B); W = W./max(W(:));
W(1:p+1:end) = 0; % zero diagonal
end

function s = tern(c,a,b), if c, s=a; else, s=b; end, end
