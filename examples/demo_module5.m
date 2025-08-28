% DEMO_MODULE5_INTEGRATED - Integrated demonstration of Modules 1-5 pipeline
%
% Pipeline: Module 7 (data) -> Module 1 (whitening) -> Module 4 (grad) -> Module 5 (prox)
%
% Author: [Your Name]
% Date: [Current Date]
% Version: 1.1

clear; clc; close all;

fprintf('=========================================\n');
fprintf('Module 5 Integrated Pipeline Demo\n');
fprintf('=========================================\n');

rng(42); % reproducibility

%% Step 1: Generate Test Data using Module 7
fprintf('\nStep 1: Generating test data using Module 7...\n');

p = 6;   % nodes
F = 4;   % frequencies
n_samples = 150;

try
    [Gamma_true, Sigma_true, Sigma_emp, sim_results] = ...
        module7_simulation_improved_complex('n_nodes', p, 'n_freq', F, ...
        'n_samples', n_samples, 'graph_type', 'random', ...
        'edge_density', 0.3, 'complex_strength', 0.7, ...
        'random_seed', 42);
    fprintf('  ✓ Using module7_simulation_improved_complex\n');
catch
    try
        [Gamma_true, Sigma_true, Sigma_emp, sim_results] = ...
            module7_simulation_improved('n_nodes', p, 'n_freq', F, ...
            'n_samples', n_samples, 'graph_type', 'random', ...
            'edge_density', 0.3, 'random_seed', 42);
        fprintf('  ✓ Using module7_simulation_improved\n');
    catch
        [Gamma_true, Sigma_true, Sigma_emp, sim_results] = ...
            module7_simulation('n_nodes', p, 'n_freq', F, ...
            'n_samples', n_samples, 'graph_type', 'random', ...
            'edge_density', 0.3, 'random_seed', 42);
        fprintf('  ✓ Using module7_simulation\n');
    end
end

fprintf('  Problem size: %d nodes, %d frequencies\n', p, F);
fprintf('  Generated %d sparse precision matrices\n', F);

% True sparsity (upper-tri off-diagonal)
true_sparsity_total = 0; true_elements_total = 0;
for f = 1:F
    elements = Gamma_true{f}(triu(true(p),1));
    true_sparsity_total = true_sparsity_total + nnz(abs(elements) < 1e-10);
    true_elements_total = true_elements_total + numel(elements);
end
fprintf('  Average sparsity: %.1f%%\n', 100 * true_sparsity_total / true_elements_total);

%% Step 2: Module 1 - Covariance Whitening
fprintf('\nStep 2: Applying Module 1 covariance whitening...\n');
whitening_operators = cell(F,1);
for f = 1:F
    S = Sigma_emp{f};
    d = sqrt(1 ./ max(abs(diag(S)),1e-12));  % target diag 1
    whitening_operators{f} = diag(d);
end

% Try class/function; fallback to manual
try
    [Sigma_whitened, ~] = CovarianceWhitening.whiten( ...
        Sigma_emp, whitening_operators, ...
        'target_diagonal', 1.0, 'force_hermitian', true, ...
        'check_psd', true, 'verbose', false);
    fprintf('  ✓ Using CovarianceWhitening class\n');
catch
    try
        [Sigma_whitened, ~] = covariance_whitening( ...
            Sigma_emp, whitening_operators, ...
            'target_diagonal', 1.0, 'force_hermitian', true, ...
            'verbose', false);
        fprintf('  ✓ Using covariance_whitening function\n');
    catch
        Sigma_whitened = cell(F,1);
        for f = 1:F
            D = whitening_operators{f};
            S_white = D * Sigma_emp{f} * D';
            S_white = (S_white + S_white')/2;
            % normalize small diag drift
            dfix = real(diag(S_white));
            D2 = diag(1 ./ sqrt(max(dfix,1e-12)));
            Sigma_whitened{f} = (D2*S_white*D2');
        end
        fprintf('  ✓ Manual whitening completed\n');
    end
end

% Quality
whitening_errors = zeros(F,1);
for f = 1:F
    whitening_errors(f) = norm(real(diag(Sigma_whitened{f})) - 1);
end
fprintf('  Whitening quality: max diagonal error = %.3e\n', max(whitening_errors));

%% Step 3: Module 2 - E-step (optional init)
fprintf('\nStep 3: Running Module 2 E-step computation...\n');
use_estep_initialization = false;
try
    n_sources = p*2;
    L = randn(p,n_sources)/sqrt(n_sources);
    module2_input = struct();
    module2_input.leadfield_matrix = L;
    module2_input.empirical_covariances = Sigma_whitened;
    module2_input.source_prior_covariances = repmat({0.3*eye(n_sources)},F,1);
    module2_input.frequencies = linspace(8,12,F);
    module2_input.noise_covariance = 0.05*eye(p);
    estep_params = struct('verbose', false, 'regularization_factor', 1e-6);
    estep_results = module2_estep_main(module2_input, estep_params);

    if estep_results.success
        fprintf('  ✓ E-step computation successful\n');
        Gamma_init = estep_results.initial_precision_matrices;
        for f = 1:F
            if isempty(Gamma_init{f})
                S = (Sigma_whitened{f}+Sigma_whitened{f}')/2;
                rr = 1e-8*trace(S)/p; Sreg = S + rr*eye(p);
                G0 = (Sreg\eye(p)); G0 = (G0+G0')/2; G0(1:p+1:end)=real(diag(G0));
                Gamma_init{f} = G0;
            end
        end
        use_estep_initialization = true;
    else
        fprintf('  ⚠ E-step computation failed, using fallback init\n');
    end
catch ME
    fprintf('  ⚠ Module 2 E-step failed: %s\n', ME.message);
end

if ~use_estep_initialization
    Gamma_init = cell(F,1);
    for f = 1:F
        S = (Sigma_whitened{f}+Sigma_whitened{f}')/2;
        rr = 1e-8*trace(S)/p; Sreg = S + rr*eye(p);
        G0 = (Sreg\eye(p)); G0 = (G0+G0')/2; G0(1:p+1:end)=real(diag(G0));
        % safety shift
        mineig = min(real(eig(G0)));
        if mineig <= 0, G0 = G0 + (abs(mineig)+1e-6)*eye(p); end
        Gamma_init{f} = G0;
    end
end
for f = 1:F
    if isempty(Gamma_init{f}) || any(~isfinite(Gamma_init{f}(:)))
        Gamma_init{f} = eye(p);
    end
    Gamma_init{f} = (Gamma_init{f} + Gamma_init{f}')/2;          % Hermitian
    Gamma_init{f}(1:p+1:end) = real(diag(Gamma_init{f}));        % real diag
end
init_err = zeros(F,1);
for f = 1:F
    init_err(f) = norm(Gamma_init{f}-Gamma_true{f},'fro')/norm(Gamma_true{f},'fro');
end
fprintf('  Initialization error range: [%.3f, %.3f]\n', min(init_err), max(init_err));

%% Step 4: Smoothing kernel & weight matrix
fprintf('\nStep 4: Setting up smoothing kernel and weight matrix...\n');

K_smooth = zeros(F,F);
for f = 1:F-1, K_smooth(f,f+1)=0.3; K_smooth(f+1,f)=0.3; end
for f = 1:F-2, K_smooth(f,f+2)=0.1; K_smooth(f+2,f)=0.1; end
fprintf('  Smoothing kernel: %d non-zero entries\n', nnz(K_smooth));

% PSD weight matrix (correlation-like)
base_corr = 0.2;
W_matrix = base_corr*ones(p) + (1.5-base_corr)*eye(p);
if isfield(sim_results,'graph_structure')
    Gs = sim_results.graph_structure;
    Wp = zeros(p);
    for i=1:p
        for j=i+1:p
            Wp(i,j) = (Gs(i,j) ~= 0) * (-0.1) + (Gs(i,j)==0) * (0.1);
            Wp(j,i) = Wp(i,j);
        end
    end
    maxneg = abs(min(real(eig(Wp))));
    alpha = 0.3 / max(maxneg,1e-10);
    W_matrix = W_matrix + alpha*Wp;
end
if min(real(eig(W_matrix))) <= 1e-12
    error('Weight matrix construction failed - not PSD!');
end
fprintf('  Weight matrix condition number: %.2f\n', cond(W_matrix));

%% Step 5: Module 3 - Active set
fprintf('\nStep 5: Running Module 3 active set selection...\n');
use_active_set = false;
try
    m3in = struct();
    m3in.whitened_covariances = Sigma_whitened;
    m3in.initial_precision_matrices = Gamma_init;
    m3in.frequencies = linspace(8,12,F);
    m3p = struct('proxy_method','correlation','quantile_level',0.15, ...
                 'force_diagonal_active',true,'verbose',false);
    m3res = module3_active_set_main(m3in,m3p);
    if m3res.success
        active_masks = cell(F,1);
        for f=1:F, active_masks{f} = m3res.combined_active_mask(:,:,f); end
        use_active_set = true;
        total_active = sum(cellfun(@(x) sum(sum(triu(x,1))), active_masks));
        total_possible = F*p*(p-1)/2;
        fprintf('  ✓ Active edges: %d / %d (%.1f%%)\n', ...
            total_active, total_possible, 100*total_active/total_possible);
    else
        fprintf('  ⚠ Active set selection failed\n');
    end
catch ME
    fprintf('  ⚠ Module 3 failed: %s\n', ME.message);
end
if ~use_active_set
    active_masks = cell(F,1);
    for f=1:F
        M = rand(p)>0.15; M = M|M'; M(1:p+1:end)=true;
        active_masks{f}=M;
    end
    total_active = sum(cellfun(@(x) sum(sum(triu(x,1))), active_masks));
    total_possible = F*p*(p-1)/2;
    fprintf('  Fallback active edges: %d / %d (%.1f%%)\n', ...
        total_active, total_possible, 100*total_active/total_possible);
end

%% Step 6: Module 4 gradient test (robust fields)
fprintf('\nStep 6: Testing Module 4 gradient computation...\n');
try
    m4in = struct('precision_matrices',Gamma_init, ...
                  'whitened_covariances',Sigma_whitened, ...
                  'smoothing_kernel',K_smooth, ...
                  'weight_matrix',W_matrix);
    m4p  = struct('lambda1',0.01,'penalize_diagonal',false,'use_graph_laplacian',true,'verbose',false);
    gout = module4_objective_gradient_main(m4in,m4p);
    Gcells = pick_gradients(gout);
    gnorms = cellfun(@(G) norm(G,'fro'), Gcells);
    fprintf('  ✓ Module 4 gradient norms: %s\n', mat2str(gnorms,3));
catch ME
    fprintf('  ⚠ Module 4 gradient computation failed: %s\n', ME.message);
    fprintf('  Continuing with Module 5...\n');
end

%% Step 7: Package input for Module 5
fprintf('\nStep 7: Preparing Module 5 input data structure...\n');
input_data = struct();
input_data.whitened_covariances = Sigma_whitened;
input_data.initial_precision    = Gamma_init;
input_data.smoothing_kernel     = K_smooth;
input_data.weight_matrix        = W_matrix;
input_data.active_set_masks     = active_masks;
fprintf('  Input data structure ready\n');

%% Step 8: Run Module 5
fprintf('\nStep 8: Running Module 5 proximal gradient optimization...\n');
params = struct();
params.lambda1  = [];     % auto
params.lambda2  = 0.02;   % show visible sparsity in demo
params.max_iter = 100;
params.eps_x    = 1e-3;
params.eps_f    = 1e-4;
params.verbose  = true;
params.use_parfor = false;

tic;
[Gamma_opt, res] = module5_proximal_main(input_data, params);
t_opt = toc;

fprintf('\n  Module 5 optimization: %s\n', tern(res.convergence_info.converged,'SUCCESS','DONE'));
fprintf('  Iterations: %d, time: %.2fs\n', res.convergence_info.iterations, t_opt);
fprintf('  Final objective: %.6e\n', res.convergence_info.final_objective);

%% Step 9: Evaluate in whitened coordinates
fprintf('\nStep 9: Evaluating end-to-end pipeline performance...\n');

% Build the effective whitening transform T_f so that S_white = T S T'
T_whiten = whitening_operators;
for f = 1:F
    S_tmp = T_whiten{f} * Sigma_emp{f} * T_whiten{f}';
    S_tmp = (S_tmp + S_tmp')/2;
    d2 = real(diag(S_tmp));
    D2 = diag(1./sqrt(max(d2,1e-12)));
    T_whiten{f} = D2 * T_whiten{f};
end

% Transform true Gamma to whitened coordinates: T^{-T} * Gamma * T^{-1}
Gamma_true_white = cell(F,1);
for f=1:F
    Tf = T_whiten{f}; invT = inv(Tf);
    Gamma_true_white{f} = invT' * Gamma_true{f} * invT;
end

recovery_errors = zeros(F,1);
achieved_sparsity = zeros(F,1);
for f=1:F
    recovery_errors(f) = norm(Gamma_opt{f}-Gamma_true_white{f},'fro')/ ...
                         norm(Gamma_true_white{f},'fro');
    offd = Gamma_opt{f}(triu(true(p),1));
    achieved_sparsity(f) = nnz(abs(offd)<1e-6) / numel(offd);
end
fprintf('  Recovery mean error: %.3f  (min=%.3f, max=%.3f)\n', ...
    mean(recovery_errors), min(recovery_errors), max(recovery_errors));
fprintf('  Achieved sparsity: mean %.1f%%\n', 100*mean(achieved_sparsity));

% PD & Hermitian checks
all_pd = true; all_herm = true; maxH = 0;
for f=1:F
    [isPD,~] = module5_psd_check(Gamma_opt{f}); % 接口不变，仅注释统一为 PD
    all_pd = all_pd && isPD;
    Herr = norm(Gamma_opt{f}-Gamma_opt{f}','fro'); maxH = max(maxH,Herr);
    all_herm = all_herm && (Herr<1e-12);
end
fprintf('  All matrices PD: %s, Hermitian: %s (max err=%.2e)\n', ...
    tern(all_pd,'YES','NO'), tern(all_herm,'YES','NO'), maxH);

%% Step 10: Sparsity sweep over lambda2
fprintf('\nStep 10: Demonstrating sparsity control across λ₂ values...\n');
lambda2_values = [1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 2e-2, 5e-2];
sp_out = zeros(size(lambda2_values));
err_out = zeros(size(lambda2_values));
obj_out = zeros(size(lambda2_values));

for i=1:numel(lambda2_values)
    p2 = params; p2.lambda2 = lambda2_values(i); p2.verbose=false; p2.max_iter=50;
    [G2, r2] = module5_proximal_main(input_data, p2);
    tot_off = 0; tot_zero = 0; tot_err2 = 0;
    for f=1:F
        el = G2{f}(triu(true(p),1));
        tot_off  = tot_off + numel(el);
        tot_zero = tot_zero + nnz(abs(el)<1e-6);
        tot_err2 = tot_err2 + norm(G2{f}-Gamma_true_white{f},'fro')^2;
    end
    sp_out(i)  = tot_zero / tot_off;
    err_out(i) = sqrt(tot_err2/F);
    obj_out(i) = r2.convergence_info.final_objective;
    fprintf('  λ₂=%6.1e → sparsity=%5.1f%%, error=%.3f\n', lambda2_values(i), 100*sp_out(i), err_out(i));
end

%% Step 11: Cross-module validation (random/hub/chain)
fprintf('\nStep 11: Cross-module validation testing...\n');
validation_configs = {
    {'graph_type','random','edge_density',0.20}
    {'graph_type','hub',   'edge_density',0.30}
    {'graph_type','chain', 'edge_density',0.25}
};
val_success = false(numel(validation_configs),1);

for k=1:numel(validation_configs)
    cfg = validation_configs{k};
    try
        [~,~,Sigma_emp_val,~] = module7_simulation_improved_complex( ...
            'n_nodes',p,'n_freq',F,'n_samples',n_samples, cfg{:}, 'random_seed',100*k);

        % quick whitening: diagonal scaling to unit variance
        Sigma_white_val = cell(F,1);
        for f=1:F
            S = Sigma_emp_val{f};
            d = sqrt(1 ./ max(abs(diag(S)),1e-12));
            D = diag(d);
            S2 = (D*S*D'); Sigma_white_val{f} = (S2+S2')/2;
        end

        in2 = input_data; in2.whitened_covariances = Sigma_white_val;
        p2  = params; p2.max_iter = 30; p2.verbose=false;
        [~,~] = module5_proximal_main(in2, p2);
        val_success(k) = true;
        fprintf('  %s graph: SUCCESS\n', cfg{2});
    catch ME
        fprintf('  %s graph: FAILED (%s)\n', cfg{2}, ME.message);
    end
end

%% Step 12: Visualization
fprintf('\nStep 12: Creating comprehensive visualization...\n');
try
    figure('Position',[100,100,1400,1000]);

    % Sparsity vs lambda2
    subplot(2,3,1);
    semilogx(lambda2_values, 100*sp_out, 'o-','LineWidth',2);
    xlabel('\lambda_2 (log scale)'); ylabel('Sparsity Achieved (%)');
    title('Sparsity Control'); grid on;

    % Recovery error vs lambda2
    subplot(2,3,2);
    semilogx(lambda2_values, err_out, 's-','LineWidth',2);
    xlabel('\lambda_2 (log scale)'); ylabel('Recovery Error');
    title('Recovery vs Sparsity Trade-off'); grid on;

    % Convergence history
    subplot(2,3,3);
    if isfield(res,'objective_history') && ~isempty(res.objective_history)
        plot(res.objective_history,'LineWidth',2);
        xlabel('Iteration'); ylabel('Objective Value'); title('Convergence History'); grid on;
    else
        text(0.5,0.5,'Objective history not available','HorizontalAlignment','center');
        title('Convergence History (N/A)'); axis off;
    end

    % True vs estimated (freq 1)
    subplot(2,3,4);
    imagesc(abs(Gamma_true_white{1})); colorbar; title('True Structure (|Γ_1|)');
    xlabel('Node j'); ylabel('Node i');

    subplot(2,3,5);
    imagesc(abs(Gamma_opt{1})); colorbar; title('Estimated Structure (|Γ̂_1|)');
    xlabel('Node j'); ylabel('Node i');

    % Cross-validation success bars
    subplot(2,3,6);
    names = {'Random','Hub','Chain'};
    vals = double(val_success);
    bar(vals,'FaceColor',[0.2 0.6 0.9]); ylim([0 1.2]); grid on;
    set(gca,'XTick',1:numel(names),'XTickLabel',names); ylabel('Success (0/1)');
    title('Cross-Validation Results');

    sgtitle('Module 5 Integrated Pipeline Results');
    fprintf('  Comprehensive visualization created\n');
catch ME
    fprintf('  Visualization failed: %s\n', ME.message);
end

%% Step 13: Baseline comparison (ridge)
fprintf('\nStep 13: Module performance comparison...\n');
ridge_err = zeros(F,1);
for f=1:F
    rr = 0.01*trace(Sigma_whitened{f})/p;
    G_r = (Sigma_whitened{f} + rr*eye(p)) \ eye(p);
    ridge_err(f) = norm(G_r - Gamma_true_white{f},'fro') / norm(Gamma_true_white{f},'fro');
end
fprintf('  Ridge baseline mean error: %.3f\n', mean(ridge_err));
fprintf('  Our pipeline mean error:   %.3f\n', mean(recovery_errors));
fprintf('  Improvement factor (Ridge/Our): %.2fx\n', mean(ridge_err)/mean(recovery_errors));

%% Step 14: Summary
fprintf('\nStep 14: Summary and recommendations\n');
fprintf('\n=== Integrated Pipeline Demo Summary ===\n');
fprintf('✓ Data: %dx%d (Module 7)\n', p, F);
fprintf('✓ Whitening max diag error: %.2e\n', max(whitening_errors));
fprintf('✓ Module 5: %d iters, %.2fs, PD=%s, Hermitian=%s\n', ...
    res.convergence_info.iterations, t_opt, tern(all_pd,'YES','NO'), tern(all_herm,'YES','NO'));
fprintf('✓ Sparsity (λ₂=%.3g): %.1f%%\n', params.lambda2, 100*mean(achieved_sparsity));
fprintf('✓ Recovery mean error (whitened coords): %.3f\n', mean(recovery_errors));
fprintf('\nBest practices:\n');
fprintf('• Always evaluate in the whitened coordinate system.\n');
fprintf('• Use auto λ₁/α; tune λ₂∈[1e-4,5e-2] to control sparsity.\n');
fprintf('• Ensure weight matrix W is PD (e.g., correlation-like + safe perturbation).\n');
fprintf('• Active-set can speed up & stabilize on sparse problems.\n');
fprintf('\n=== Integrated Demo Complete ===\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local helpers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = tern(cond,a,b), if cond, s=a; else, s=b; end, end

function Gcells = pick_gradients(gout)
% Robustly extract gradient cell array from Module 4 output
if isfield(gout,'gradients'), Gcells = gout.gradients; return; end
if isfield(gout,'smooth_gradients'), Gcells = gout.smooth_gradients; return; end
if isfield(gout,'gradient_components') && isfield(gout.gradient_components,'smoothing_gradients')
    Gcells = gout.gradient_components.smoothing_gradients; return;
end
error('Module 4: gradients not found in output struct.');
end
