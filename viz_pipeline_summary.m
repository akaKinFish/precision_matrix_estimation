function viz_pipeline_summary(prox_res, Gamma_tilde_star, Sjj_tilde, K, A_masks, Omega_src, Omega_true)
% VIZ_PIPELINE_SUMMARY - Post-run visualization & quick diagnostics
% Inputs:
%   prox_res         : struct from module5_proximal (objective_history, gradient_norm_history, etc.)
%   Gamma_tilde_star : {F×1} final precision in whitened domain
%   Sjj_tilde        : {F×1} whitened covariances
%   K                : F×F smoothing kernel
%   A_masks          : {F×1} logical masks for active set (p×p)
%   Omega_src        : {F×1} recolored precision (source scale)
%   Omega_true       : {F×1} or p×p×F (optional, ground-truth from sim)

F = numel(Gamma_tilde_star);
p = size(Gamma_tilde_star{1},1);
f_view = 1;  % default frequency to visualize

%% 0) Figure layout
figure('Name','Module5 Summary','NumberTitle','off'); clf;

%% 1) Convergence curves
subplot(2,3,1);
yyaxis left;
semilogy(max(prox_res.objective_history - min(prox_res.objective_history), eps), '-o'); grid on;
ylabel('Objective (shifted, log)');
yyaxis right;
semilogy(max(prox_res.gradient_norm_history, eps), '-s');
ylabel('||grad||_F (log)');
xlabel('Iteration');
title('Convergence');

%% 2) Step size & backtracks
subplot(2,3,2);
stairs(prox_res.step_size_history, '-'); grid on; hold on;
yyaxis right;
if ~isempty(prox_res.backtrack_counts)
    stem(1:numel(prox_res.backtrack_counts), prox_res.backtrack_counts, '.');
end
xlabel('Iteration'); ylabel('Backtracks');
yyaxis left; ylabel('Step size');
title('Step size & backtracks');

%% 3) Active-set dynamics
subplot(2,3,3);
if isfield(prox_res,'active_set_changes') && ~isempty(prox_res.active_set_changes)
    bar(prox_res.active_set_changes);
    title('Active-set changes / iter'); xlabel('Iteration'); ylabel('# flips');
else
    text(0.1,0.5,'No active-set change history','FontSize',12);
    axis off
end
grid on;

%% 4) Fit residual in whitened domain: ||S_f * Gamma_f - I||_F / ||I||
res = zeros(F,1);
for f=1:F
    Af = Sjj_tilde{f} * Gamma_tilde_star{f};
    res(f) = norm(Af - eye(p),'fro') / sqrt(p);
end
subplot(2,3,4);
plot(res,'-o'); grid on;
xlabel('Frequency'); ylabel('Fro residual');
title('Whitened fit: ||SΓ - I||_F');

%% 5) Sparsity per frequency (off-diagonal nnz)
nnz_off = zeros(F,1);
for f=1:F
    G = Gamma_tilde_star{f};
    nnz_off(f) = nnz(abs(G) > 1e-10) - p;
end
subplot(2,3,5);
plot(nnz_off,'-o'); grid on;
xlabel('Frequency'); ylabel('# offdiag nz');
title('Estimated sparsity (whitened Γ)');

%% 6) Heatmap: partial coherence at a chosen frequency (source scale if available)
subplot(2,3,6);
if nargin>=6 && ~isempty(Omega_src)
    Om = Omega_src{f_view};
    pcorr = abs(-Om) ./ sqrt((abs(diag(Om))+eps) * (abs(diag(Om))+eps)');
    pcorr(1:p+1:end) = 0;
    imagesc(pcorr); colorbar; title(sprintf('Partial coherence (f=%d)', f_view));
else
    G = Gamma_tilde_star{f_view};
    pc = abs(-G) ./ sqrt((abs(diag(G))+eps) * (abs(diag(G))+eps)');
    pc(1:p+1:end) = 0;
    imagesc(pc); colorbar; title(sprintf('Partial coherence (whitened, f=%d)', f_view));
end
axis square; xlabel('j'); ylabel('i');

drawnow;

%% 7) If ground-truth exists: precision/recall & F1 (off-diagonal support)
if nargin>=7 && ~isempty(Omega_true)
    Omega_true_cell = coerce_cov_cell_local(Omega_true, F);
    thr = 1e-8;
    tp=0; fp=0; fn=0;
    for f=1:F
        T = Omega_true_cell{f}; T(1:p+1:end)=0;
        Etrue = abs(T) > thr;

        if nargin>=6 && ~isempty(Omega_src)
            Est = Omega_src{f};
        else
            Est = Gamma_tilde_star{f};
        end
        Est(1:p+1:end)=0;
        Eest = abs(Est) > thr;

        % undirected counting (upper triangular)
        T_ut = triu(Etrue,1); E_ut = triu(Eest,1);
        tp = tp + nnz(T_ut & E_ut);
        fp = fp + nnz(~T_ut & E_ut);
        fn = fn + nnz(T_ut & ~E_ut);
    end
    prec = tp / max(tp+fp,1);
    rec  = tp / max(tp+fn,1);
    F1   = 2*prec*rec / max(prec+rec,eps);
    fprintf('[GT] Precision=%.3f  Recall=%.3f  F1=%.3f  (tp=%d, fp=%d, fn=%d)\n',prec,rec,F1,tp,fp,fn);
end

%% 8) Optional: kernel sanity
fprintf('Kernel symmetry error ||K-K''||_F = %.2e\n', norm(K-K','fro'));
end

function C = coerce_cov_cell_local(X, F_hint)
% Accepts cell{F,1} | p×p×F | single p×p (replicate)
if isa(X,'cell'), C = X(:); return; end
if isnumeric(X) && ndims(X)==3 && size(X,1)==size(X,2)
    F = size(X,3); C = cell(F,1); for f=1:F, C{f}=X(:,:,f); end; return;
end
if isnumeric(X) && ismatrix(X) && size(X,1)==size(X,2)
    F = F_hint; C = repmat({X}, F, 1); return;
end
error('coerce_cov_cell_local:unsupported','Expect cell{F,1} | p×p×F | single p×p.');
end
