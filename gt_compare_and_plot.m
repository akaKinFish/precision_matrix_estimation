function metrics = gt_compare_and_plot(Omega_true, Omega_est, opts, D_src, Gamma_tilde_star)
% GT_COMPARE_AND_PLOT - Compare estimated precision with ground truth.
% Inputs:
%   Omega_true : {F×1} or p×p×F   (ground truth, source domain)
%   Omega_est  : {F×1}            (estimate, source domain; e.g., Omega_src)
%   opts       : struct with fields:
%                  .mode        : 'absolute'|'match_sparsity' (default 'match_sparsity')
%                  .thr         : absolute threshold when mode='absolute' (default 1e-8)
%                  .f_view      : frequency index for heatmaps (default 1)
%                  .show_pr     : show aggregated PR curve (default true)
%                  .show_roc    : show aggregated ROC curve (default true)
%                  .show_curves : per-fig plotting switch for curves (default false)
%                  .Gamma       : (opt) {F×1} or p×p weight/prior Γ for correlation with |Ω_true|
%   D_src (opt): {F×1} whitening matrices; if given且提供 Gamma_tilde_star，
%                会计算“白化域”真值用于和 Gamma_tilde_star 对比。
%   Gamma_tilde_star (opt): {F×1} estimated precision in whitened domain.
%
% Output:
%   metrics : struct with fields:
%       .per_freq.precision / .recall / .F1
%       .per_freq.relFrob_mag / .fro_rel_ut / .fro_rel_off
%       .per_freq.kl / .stein
%       .overall.precision / .recall / .F1  (micro)
%       .agg_pr.(precision,recall,aupr) / .agg_roc.(fpr,tpr,auroc)
%       .corr.(pearson,spearman)  (per-freq & mean)
%       .whitened_relFrob (if D_src & Gamma_tilde_star provided)

if nargin < 3, opts = struct(); end
if ~isfield(opts,'mode'),        opts.mode = 'match_sparsity'; end
if ~isfield(opts,'thr'),         opts.thr  = 1e-8; end
if ~isfield(opts,'f_view'),      opts.f_view = 1; end
if ~isfield(opts,'show_pr'),     opts.show_pr = true; end
if ~isfield(opts,'show_roc'),    opts.show_roc = true; end
if ~isfield(opts,'show_curves'), opts.show_curves = false; end

% ---- Coerce to cell ----
Omega_true = coerce_cell_(Omega_true);
F = numel(Omega_true);
p = size(Omega_true{1},1);
assert(iscell(Omega_est) && numel(Omega_est)==F, 'Omega_est must be {F×1}');

% ---- Per-frequency metrics ----
tp=0; fp=0; fn=0;
prec=zeros(F,1); rec=zeros(F,1); F1=zeros(F,1);
rfe_mag=zeros(F,1);
fro_rel_ut=zeros(F,1); fro_rel_off=zeros(F,1);
kl_risk=zeros(F,1); stein_risk=zeros(F,1);

% optional: correlation |Gamma| vs |Omega_true|
haveGamma = isfield(opts,'Gamma') && ~isempty(opts.Gamma);
pear_r = NaN(F,1); spear_r = NaN(F,1);

for f=1:F
    T = Omega_true{f};  E = Omega_est{f};
    T = zerodiag_(T);   E = zerodiag_(E);

    % relative Frobenius (magnitudes, full)
    rfe_mag(f) = norm(abs(E)-abs(T),'fro') / max(norm(abs(T),'fro'), 1e-12);

    % ===== Off-diagonal Frobenius relative error =====
    E_ut = triu(E,1);   T_ut = triu(T,1);
    fro_rel_ut(f)  = norm(E_ut - T_ut, 'fro') / max(1e-12, norm(T_ut, 'fro'));
    E_off = E - diag(diag(E)); T_off = T - diag(diag(T));
    fro_rel_off(f) = norm(E_off - T_off, 'fro') / max(1e-12, norm(T_off, 'fro'));

    % ===== support sets (undirected; count upper triangle only) =====
    switch lower(opts.mode)
        case 'absolute'
            Strue = triu(abs(T) > opts.thr, 1);
            Sest  = triu(abs(E) > opts.thr, 1);
        case 'match_sparsity'
            % keep top-k estimated edges where k = #true edges
            Strue = triu(abs(T) > opts.thr, 1);
            k = nnz(Strue);
            Sest = topk_mask_ut(abs(E), k);   % ★ 只返回上三角掩码（无对称回填）
        otherwise
            error('Unknown opts.mode: %s', opts.mode);
    end

    % counts (upper triangle only)
    tp_f = nnz(Strue & Sest);
    fp_f = nnz(~Strue & Sest);
    fn_f = nnz(Strue & ~Sest);

    tp = tp + tp_f; fp = fp + fp_f; fn = fn + fn_f;

    prec(f) = tp_f / max(tp_f+fp_f,1);
    rec(f)  = tp_f / max(tp_f+fn_f,1);
    F1(f)   = 2*prec(f)*rec(f) / max(prec(f)+rec(f), eps);

    % ===== Stein / KL risks based on true covariance =====
    % KL(Σ_true || Σ_est) = 0.5 [ tr(Σ_est^{-1} Σ_true) − logdet(Σ_est^{-1} Σ_true) − p ]
    % Stein(Σ_true^{-1}, Σ_est^{-1}) = tr(Σ_true^{-1} Σ_est) − logdet(Σ_true^{-1} Σ_est) − p
    Sigma_true = symmetrize_(safe_inv_(T));
    Sigma_est  = symmetrize_(safe_inv_(E));

    A = Sigma_est \ Sigma_true;        % Σ_est^{-1} Σ_true
    kl_risk(f)    = 0.5 * ( trace(A) - logdet_spd(A) - p );

    B = Sigma_true \ Sigma_est;        % Σ_true^{-1} Σ_est
    stein_risk(f) = trace(B) - logdet_spd(B) - p;

    % ===== |Gamma| vs |Omega_true| correlation (upper triangle) =====
    if haveGamma
        if iscell(opts.Gamma), G = opts.Gamma{f}; else, G = opts.Gamma; end
        G_ut = upper_tri_vec(abs(G));
        T_utv = upper_tri_vec(abs(T));
        keep = isfinite(G_ut) & isfinite(T_utv);
        if any(keep) && numel(unique(G_ut(keep)))>1 && numel(unique(T_utv(keep)))>1
            pear_r(f)  = corr(G_ut(keep), T_utv(keep), 'type','Pearson');
            spear_r(f) = corr(G_ut(keep), T_utv(keep), 'type','Spearman');
        end
    end
end

% ---- micro-average ----
prec_all = tp / max(tp+fp,1);
rec_all  = tp / max(tp+fn,1);
F1_all   = 2*prec_all*rec_all / max(prec_all+rec_all, eps);

metrics = struct();
metrics.per_freq = struct( ...
    'precision',prec,'recall',rec,'F1',F1, ...
    'relFrob_mag',rfe_mag, ...
    'fro_rel_ut',fro_rel_ut,'fro_rel_off',fro_rel_off, ...
    'kl',kl_risk,'stein',stein_risk);
metrics.overall  = struct('precision',prec_all,'recall',rec_all,'F1',F1_all);

% ---- Aggregated PR & ROC over all freqs (upper-tri, truth=nonzero) ----
[pr_R, pr_P, aupr] = aggregate_pr_(Omega_true, Omega_est);
[roc_FPR, roc_TPR, auroc] = aggregate_roc_(Omega_true, Omega_est);
metrics.agg_pr  = struct('recall',pr_R,'precision',pr_P,'aupr',aupr);
metrics.agg_roc = struct('fpr',roc_FPR,'tpr',roc_TPR,'auroc',auroc);

fprintf(['[GT-compare] mode=%s | Overall: P=%.3f  R=%.3f  F1=%.3f | ' ...
         'meanRelFrobMag=%.3f | meanFroUT=%.3f | meanFroOff=%.3f\n'], ...
    opts.mode, prec_all, rec_all, F1_all, mean(rfe_mag), mean(fro_rel_ut), mean(fro_rel_off));
fprintf('[Curves ] AUPR=%.3f | AUROC=%.3f\n', aupr, auroc);
fprintf('[Risk   ] KL=%.3f | Stein=%.3f (mean over F)\n', mean(kl_risk), mean(stein_risk));

if haveGamma
    fprintf('[Corr   ] |Gamma|-|Omega_true| : Pearson=%.3f | Spearman=%.3f (means)\n', ...
        mean(pear_r,'omitnan'), mean(spear_r,'omitnan'));
    metrics.corr = struct('pearson_per_freq',pear_r,'spearman_per_freq',spear_r, ...
                          'pearson_mean',mean(pear_r,'omitnan'),'spearman_mean',mean(spear_r,'omitnan'));
end

% ---- Figure: per-frequency bars + heatmaps ----
figure('Name','GT vs Estimate (source domain)','NumberTitle','off'); clf;

subplot(2,3,1); bar(F1); grid on; ylim([0 1]); title('F1 per freq'); xlabel('f'); ylabel('F1');
subplot(2,3,2); plot(prec,'-o'); hold on; plot(rec,'-s'); grid on;
legend({'Precision','Recall'},'Location','best'); title('P/R per freq'); xlabel('f'); ylim([0 1]);
subplot(2,3,3); plot(rfe_mag,'-o'); grid on; title('Rel Frobenius (|.|) per freq'); xlabel('f'); ylabel('relFrob');

f_view = min(max(1,opts.f_view), F);
subplot(2,3,4);
imagesc(abs(Omega_true{f_view})); axis image; colorbar; title(sprintf('|\\Omega_{true}| (f=%d)',f_view));
subplot(2,3,5);
imagesc(abs(Omega_est{f_view}));  axis image; colorbar; title(sprintf('|\\Omega_{est}| (f=%d)',f_view));

% Scatter of edge magnitudes (upper triangle)
tv = upper_tri_vec(abs(Omega_true{f_view}));
ev = upper_tri_vec(abs(Omega_est{f_view}));
subplot(2,3,6); scatter(tv, ev, 6, 'filled'); grid on;
xlabel('|\\Omega_{true}|'); ylabel('|\\Omega_{est}|'); title('Edge magnitude scatter (UT)');

% ---- Aggregated PR/ROC plots ----
if opts.show_pr
    figure('Name','Aggregated PR','NumberTitle','off'); clf;
    plot(pr_R, pr_P, '-o','LineWidth',1.2,'MarkerSize',4); grid on; xlim([0 1]); ylim([0 1]);
    xlabel('Recall'); ylabel('Precision'); title(sprintf('Aggregated PR (AUPR=%.3f)', aupr));
end
if opts.show_roc
    figure('Name','Aggregated ROC','NumberTitle','off'); clf;
    plot(roc_FPR, roc_TPR, '-o','LineWidth',1.2,'MarkerSize',4); grid on; xlim([0 1]); ylim([0 1]);
    xlabel('FPR'); ylabel('TPR'); title(sprintf('Aggregated ROC (AUROC=%.3f)', auroc));
end

% ---- Optional: compare in whitened domain ----
if nargin >= 5 && ~isempty(D_src) && ~isempty(Gamma_tilde_star)
    G_true = cell(F,1);
    for f=1:F
        G_true{f} = D_src{f} * Omega_true{f} * D_src{f}; % whitened truth
    end
    rfe_w = zeros(F,1);
    for f=1:F
        rfe_w(f) = norm(abs(Gamma_tilde_star{f})-abs(G_true{f}),'fro') ...
                   / max(norm(abs(G_true{f}),'fro'),1e-12);
    end
    figure('Name','Whitened-domain relative error','NumberTitle','off'); clf;
    plot(rfe_w,'-o'); grid on; xlabel('f'); ylabel('relFrob'); title('Whitened-domain relFrob vs truth');
    metrics.whitened_relFrob = rfe_w;
end
end

% ================= helpers =================

function C = coerce_cell_(X)
if isa(X,'cell'), C = X(:); return; end
if isnumeric(X) && ndims(X)==3 && size(X,1)==size(X,2)
    F=size(X,3); C=cell(F,1); for f=1:F, C{f}=X(:,:,f); end; return;
end
if isnumeric(X) && ismatrix(X) && size(X,1)==size(X,2)
    C={X}; return;
end
error('coerce_cell_:unsupported','Expect cell{F,1} | p×p×F | single p×p');
end

function A = zerodiag_(A)
p=size(A,1); A(1:p+1:end)=0;
end

function X = symmetrize_(X)
X = (X+X')/2;
end

function M = topk_mask_ut(A, k)
% Return upper-triangle logical mask with top-k magnitudes (no symmetrization)
p=size(A,1); A(1:p+1:end)=0;              % zero diag
U = triu(A,1);                            % upper triangle only
[~,ord] = sort(U(:),'descend');           % sort within the p×p grid
ord = ord(1:min(k,numel(ord)));           % top-k positions in the same grid
M = false(p); M(ord) = true;              % set those positions
M = triu(M,1);                            % keep UT strictly (safety)
end

function v = upper_tri_vec(A)
p=size(A,1); A(1:p+1:end)=0; v = A(triu(true(p),1));
end

function L = logdet_spd(M)
% Stable logdet for (near) SPD matrices via Cholesky with fallback
M = (M+M')/2;
[Lc,flag] = chol(M,'lower');
if flag~=0
    e = eig(M); e = e(e>0);
    L = sum(log(e + eps));
else
    L = 2*sum(log(diag(Lc) + eps));
end
end

function Xinv = safe_inv_(X)
% Invert SPD/HPD with fallback regularization if needed
X = (X+X')/2;
[U,D] = eig(X);
d = real(diag(D));
d(d < 1e-12) = 1e-12;  % floor small/neg eigenvalues
Xinv = U * diag(1./d) * U';
Xinv = (Xinv+Xinv')/2;
end

function [R, P, AUPR] = aggregate_pr_(Otrue, Oest)
% PR over all freqs, truth = (|Omega_true| > 0) on UT; scores = |Omega_est|
truth = []; scores = [];
F = numel(Otrue); p = size(Otrue{1},1);
for f=1:F
    T=abs(Otrue{f}); E=abs(Oest{f});
    T(1:p+1:end)=0; E(1:p+1:end)=0;
    truth  = [truth;  upper_tri_vec(T) > 0];
    scores = [scores; upper_tri_vec(E)];
end
% sort by score desc
[s, ord] = sort(scores(:),'descend');
y = truth(:); y = y(ord);
tp = cumsum(y);
fp = cumsum(~y);
P = tp ./ max(1,(tp+fp));           % precision
R = tp / max(1,sum(y));             % recall
% remove NaNs at the very beginning if any
bad = isnan(P) | isnan(R);
P(bad)=[]; R(bad)=[];
% AUPR (trapz over recall)
AUPR = trapz(R, P);
end

function [FPR, TPR, AUROC] = aggregate_roc_(Otrue, Oest)
% ROC over all freqs, truth = (|Omega_true| > 0) on UT; scores = |Omega_est|
truth = []; scores = [];
F = numel(Otrue); p = size(Otrue{1},1);
for f=1:F
    T=abs(Otrue{f}); E=abs(Oest{f});
    T(1:p+1:end)=0; E(1:p+1:end)=0;
    truth  = [truth;  upper_tri_vec(T) > 0];
    scores = [scores; upper_tri_vec(E)];
end
[s, ord] = sort(scores(:),'descend');
y = truth(:); y = y(ord);
P = max(1,sum(y)); N = max(1,sum(~y));
tp = cumsum(y);      fn = P - tp;
fp = cumsum(~y);     tn = N - fp;
TPR = tp / P;  FPR = fp / N;
% prepend (0,0) and append (1,1) for a proper ROC polyline
FPR = [0; FPR; 1];
TPR = [0; TPR; 1];
% AUROC
AUROC = trapz(FPR, TPR);
end
