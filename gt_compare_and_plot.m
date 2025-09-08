function metrics = gt_compare_and_plot(Omega_true, Omega_est, opts, D_src, Gamma_tilde_star)
% GT_COMPARE_AND_PLOT - Compare estimated precision with ground truth.
% Inputs:
%   Omega_true : {F×1} or p×p×F   (ground truth, source domain)
%   Omega_est  : {F×1}            (estimate, source domain; e.g., Omega_src)
%   opts       : struct with fields:
%                  .mode   : 'absolute'|'match_sparsity' (default 'match_sparsity')
%                  .thr    : absolute threshold when mode='absolute' (default 1e-8)
%                  .f_view : frequency index for heatmaps (default 1)
%                  .show_pr: show aggregated PR curve (default true)
%   D_src (opt): {F×1} whitening matrices; if given且提供 Gamma_tilde_star，
%                会计算“白化域”真值用于和 Gamma_tilde_star 对比。
%   Gamma_tilde_star (opt): {F×1} estimated precision in whitened domain.
%
% Output:
%   metrics : struct with fields:
%       .per_freq.F1 / .precision / .recall / .relFrob
%       .overall.precision / .recall / .F1  (micro)
%       .agg_pr.thresholds / .precision / .recall  (if show_pr)

if nargin < 3, opts = struct(); end
if ~isfield(opts,'mode'),     opts.mode = 'match_sparsity'; end
if ~isfield(opts,'thr'),      opts.thr  = 1e-8; end
if ~isfield(opts,'f_view'),   opts.f_view = 1; end
if ~isfield(opts,'show_pr'),  opts.show_pr = true; end

% ---- Coerce to cell ----
Omega_true = coerce_cell_(Omega_true);
F = numel(Omega_true);
p = size(Omega_true{1},1);
assert(iscell(Omega_est) && numel(Omega_est)==F, 'Omega_est must be {F×1}');

% ---- Per-frequency metrics ----
tp=0; fp=0; fn=0;
prec=zeros(F,1); rec=zeros(F,1); F1=zeros(F,1); rfe=zeros(F,1);

for f=1:F
    T = Omega_true{f};  E = Omega_est{f};
    T(1:p+1:end)=0; E(1:p+1:end)=0;

    % relative Frobenius error (magnitudes)
    rfe(f) = norm(abs(E)-abs(T),'fro') / max(norm(abs(T),'fro'), 1e-12);

    % ---- support sets (undirected; count upper triangle only) ----
    switch lower(opts.mode)
        case 'absolute'
            Strue = triu(abs(T) > opts.thr, 1);
            Sest  = triu(abs(E) > opts.thr, 1);
        case 'match_sparsity'
            % keep top-k estimated edges where k = #true edges
            Strue = triu(abs(T) > opts.thr, 1);
            k = nnz(Strue);
            Sest = topk_mask_ut(abs(E), k);
        otherwise
    end

    % counts
    tp_f = nnz(Strue & Sest);
    fp_f = nnz(~Strue & Sest);
    fn_f = nnz(Strue & ~Sest);

    tp = tp + tp_f; fp = fp + fp_f; fn = fn + fn_f;

    prec(f) = tp_f / max(tp_f+fp_f,1);
    rec(f)  = tp_f / max(tp_f+fn_f,1);
    F1(f)   = 2*prec(f)*rec(f) / max(prec(f)+rec(f), eps);
end

% ---- micro-average ----
prec_all = tp / max(tp+fp,1);
rec_all  = tp / max(tp+fn,1);
F1_all   = 2*prec_all*rec_all / max(prec_all+rec_all, eps);

metrics = struct();
metrics.per_freq = struct('precision',prec,'recall',rec,'F1',F1,'relFrob',rfe);
metrics.overall  = struct('precision',prec_all,'recall',rec_all,'F1',F1_all);

fprintf('[GT-compare] mode=%s | Overall: P=%.3f  R=%.3f  F1=%.3f | meanRelFrob=%.3f\n', ...
    opts.mode, prec_all, rec_all, F1_all, mean(rfe));

% ---- Figure: per-frequency bars + heatmaps ----
figure('Name','GT vs Estimate (source domain)','NumberTitle','off'); clf;

subplot(2,3,1); bar(F1); grid on; ylim([0 1]); title('F1 per freq'); xlabel('f'); ylabel('F1');
subplot(2,3,2); plot(prec,'-o'); hold on; plot(rec,'-s'); grid on;
legend({'Precision','Recall'},'Location','best'); title('P/R per freq'); xlabel('f'); ylim([0 1]);
subplot(2,3,3); plot(rfe,'-o'); grid on; title('Rel Frobenius error per freq'); xlabel('f'); ylabel('relFrob');

f_view = min(max(1,opts.f_view), F);
subplot(2,3,4);
imagesc(abs(Omega_true{f_view})); axis image; colorbar; title(sprintf('|\\Omega_{true}| (f=%d)',f_view));
subplot(2,3,5);
imagesc(abs(Omega_est{f_view}));  axis image; colorbar; title(sprintf('|\\Omega_{est}| (f=%d)',f_view));

% Scatter of edge magnitudes (upper triangle)
tv = upper_tri_vec(abs(Omega_true{f_view}));
ev = upper_tri_vec(abs(Omega_est{f_view}));
subplot(2,3,6); scatter(tv, ev, 5, 'filled'); grid on;
xlabel('|\\Omega_{true}|'); ylabel('|\\Omega_{est}|'); title('Edge magnitude scatter (UT)');

% ---- Aggregated PR curve (sweep thresholds over all edges) ----
if opts.show_pr
    [pr_thr, prP, prR] = aggregate_pr_curve(Omega_true, Omega_est);
    figure('Name','Aggregated PR','NumberTitle','off'); clf;
    plot(prR, prP, '-o'); grid on; xlim([0 1]); ylim([0 1]);
    xlabel('Recall'); ylabel('Precision'); title('Aggregated Precision–Recall (all freqs)');
    metrics.agg_pr = struct('thresholds',pr_thr,'precision',prP,'recall',prR);
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

% ===== helpers =====
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

function m = topk_mask_ut(A, k)
% Return symmetric logical mask with top-k magnitudes in the upper triangle (off-diag)
p=size(A,1); A(1:p+1:end)=0; U = triu(A,1);
[~,idx] = sort(U(:),'descend'); idx = idx(1:min(k,numel(idx)));
M = false(p); M(idx) = true; M = M | M';
m = M;
end

function v = upper_tri_vec(A)
% Upper-triangle (off-diag) vector (column)
p=size(A,1); A(1:p+1:end)=0; v = A(triu(true(p),1));
end

function [thr, P, R] = aggregate_pr_curve(Otrue, Oest)
% Sweep thresholds over *estimated* magnitudes (all freqs) vs binary truth
F=numel(Otrue); p=size(Otrue{1},1);
% collect
truth=[]; scores=[];
for f=1:F
    T=abs(Otrue{f}); E=abs(Oest{f});
    T(1:p+1:end)=0; E(1:p+1:end)=0;
    truth  = [truth;  upper_tri_vec(T)>0]; % using tiny>0 as truth; change if needed
    scores = [scores; upper_tri_vec(E)];
end
% thresholds = unique scores (downsample for speed)
scores = scores(:);
us = unique(scores(~isnan(scores)));
if numel(us)>200, q=linspace(0,1,200); thr=quantile(us,q); else, thr=us; end
P=zeros(numel(thr),1); R=P;
for i=1:numel(thr)
    est = scores >= thr(i);
    tp = sum(truth & est);
    fp = sum(~truth & est);
    fn = sum(truth & ~est);
    P(i) = tp / max(tp+fp,1);
    R(i) = tp / max(tp+fn,1);
end
end
