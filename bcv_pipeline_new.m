function results = bcv_pipeline_new(data, methods, eval_cfg, outdir)
%BCV_PIPELINE_NEW  Minimal, modular evaluation pipeline following the
% "A/B/C" logic:
%   A) Partial-coherence side-by-side panels per frequency
%   B) Likelihood-driven soft-threshold scan to pick tau
%   C) SPD matrix distance table (KLS/LE/AIRM/JBLD/Bures–W)
%
% This file is self-contained (helpers at bottom) and does **not** depend
% on your simulators/estimators. You only need to prepare:
%   - data.emp_cov_cell : {F×1} cell of sensor covariances (n×n)
%   - data.L            : leadfield (n×p)
%   - data.T            : #samples per frequency (scalar or F×1)
%   - data.Omega_true   : (p×p×F) ground-truth precision  [optional]
%   - data.Sigma_true   : (p×p×F) ground-truth covariance [optional]
%   - data.bayes_truth  : noise info (see make_source_cov_proxy) [optional]
%
% and pass a list of methods:
%   methods = struct('name',{},'estimate',{},'Omega',{});
%     • If you already have the estimate, set .Omega = (p×p×F) or {F×1}.
%     • Or supply a function handle: Om = methods(k).estimate(data).
%
% eval_cfg fields (all optional):
%   - panel_clip : 0.99  (quantile to clip color scale)
%   - panel_caxis: []    (fixed [0 vmax])
%   - tau_scan   : linspace(0,0.3,31)
%   - dist_ref   : 'sigma_true' or 'omega_true'
%   - enable_prroc: false (kept for parity; not used by default here)
%   - adapter : struct with shrink params (spd_delta, spd_eps)
%
% Example:
%   data.emp_cov_cell = {randn(8,8)^2}; data.emp_cov_cell{1}=hermitizeC(data.emp_cov_cell{1});
%   data.L = randn(8,10); data.T = 500; 
%   methods(1).name='Dummy'; methods(1).Omega = repmat(eye(10),1,1,1);
%   R = bcv_pipeline_new(data, methods, struct('tau_scan',linspace(0,0.2,21)));
%
% Output "results" struct contains:
%   .cfg, .outdir, .Omega_true, .methods, .Omega_all, .bcv_eval
%   .bcv_eval.panel_files, .bcv_eval.lik_scan, .bcv_eval.dist_tables
%
% Author: you + ChatGPT (BCV-style template)
% ---------------------------------------------------------------------

if nargin < 3 || isempty(eval_cfg), eval_cfg = struct(); end
if nargin < 4 || isempty(outdir)
    ts = datestr(now,'yyyymmdd_HHMMSS');
    outdir = fullfile('results', ['exp_', ts]);
end
if ~exist(outdir,'dir'), mkdir(outdir); end

% ---- defaults ----
eval_cfg = set_default(eval_cfg, struct( ...
    'panel_clip', 0.99, ...
    'panel_caxis', [], ...
    'tau_scan', linspace(0,0.3,31), ...
    'dist_ref', 'sigma_true', ...
    'enable_prroc', false, ...
    'adapter', struct('spd_delta',0.03,'spd_eps',1e-6) ...
));

fprintf('================= bcv_pipeline_new =================\n');

% ---- gather basic shapes ----
F = numel(data.emp_cov_cell);
if isfield(data,'Omega_true') && ~isempty(data.Omega_true)
    OmTrue = coerce_precision_stack(data.Omega_true);
else
    OmTrue = [];
end
if isfield(data,'Sigma_true') && ~isempty(data.Sigma_true)
    SiTrue = coerce_precision_stack(data.Sigma_true);
else
    SiTrue = [];
end

% ========== (0) Produce/collect estimates ==========
Omega_all = struct();
for k = 1:numel(methods)
    nm = make_valid_name(methods(k).name);
    if isfield(methods(k),'Omega') && ~isempty(methods(k).Omega)
        Om = coerce_precision_stack(methods(k).Omega);
    elseif isfield(methods(k),'estimate') && ~isempty(methods(k).estimate)
        Om = methods(k).estimate(data);  % must return p×p×F or {F×1}
        Om = coerce_precision_stack(Om);
    else
        warning('Method %s has neither .Omega nor .estimate; skip.', nm);
        continue;
    end
    Omega_all.(nm) = Om;
end
method_names = fieldnames(Omega_all);

% ========== (A) Partial-coherence panels ==========
panel_files = plot_pcoh_panels(Omega_all, OmTrue, eval_cfg, outdir);

% ========== (B) Likelihood-driven tau scan ==========
Sii_proxy = make_source_cov_proxy(data, eval_cfg);
lik_scan = struct();
for k = 1:numel(method_names)
    name = method_names{k};
    OmEst = coerce_precision_stack(Omega_all.(name));
    [best_tau, tau_grid, ll_curve] = likelihood_scan_tau(OmEst, Sii_proxy, eval_cfg.tau_scan, eval_cfg.adapter);
    lik_scan.(name) = struct('best_tau',best_tau, 'tau_grid',tau_grid, 'll_curve',ll_curve);
    fprintf('  [lik-scan] %s: best tau=%.4g (grid in [%.3g, %.3g], %d pts)\n', ...
        name, best_tau, min(tau_grid), max(tau_grid), numel(tau_grid));
end
save(fullfile(outdir, 'likelihood_scan.mat'), 'lik_scan');

% ========== (C) SPD distance summary ==========
dist_summary = struct();
for k = 1:numel(method_names)
    name = method_names{k};
    OmEst = coerce_precision_stack(Omega_all.(name));
    % Prefer comparing covariances (numerically stabler)
    SigHat = zeros(size(OmEst));
    for f=1:size(OmEst,3), SigHat(:,:,f) = inv_psd(OmEst(:,:,f)); end
    if ~isempty(SiTrue) && strcmpi(eval_cfg.dist_ref,'sigma_true')
        Astack = SigHat; Bstack = SiTrue;
    elseif ~isempty(OmTrue) && strcmpi(eval_cfg.dist_ref,'omega_true')
        Astack = OmEst; Bstack = OmTrue;
    else
        warning('dist_ref=%s but reference missing; falling back to omega_true if available.', eval_cfg.dist_ref);
        if ~isempty(OmTrue), Astack=OmEst; Bstack=OmTrue; else Astack=SigHat; Bstack=SigHat; end
    end
    D = spd_distance_pack(Astack, Bstack);
    dist_summary.(name) = summarize_distances(D);
end
csv_path = fullfile(outdir,'spd_distance_summary.csv');
write_distance_csv(csv_path, dist_summary);
fprintf('  [dist] SPD distance summary saved: %s\n', csv_path);

% ========== pack results ==========
results = struct();
results.cfg = eval_cfg; results.outdir = outdir;
results.Omega_true = OmTrue; results.methods = methods;
results.Omega_all = Omega_all;
results.bcv_eval = struct('panel_files', {panel_files}, 'lik_scan', lik_scan, 'dist_tables', dist_summary);

save(fullfile(outdir, 'results_all.mat'), '-struct', 'results');
fprintf('\n[OK] Done. Artifacts saved under: %s\n', outdir);

end

% =====================================================================
%                              HELPERS
% =====================================================================
function panel_files = plot_pcoh_panels(Omega_all, OmTrue, eval_cfg, outdir)
% 1 row per frequency; columns: GroundTruth (if any) + each method.
% Color limits are clipped by eval_cfg.panel_clip for robust visualization.
    names = fieldnames(Omega_all);
    if ~isempty(OmTrue), [~,~,F] = size(OmTrue); else, F = size(Omega_all.(names{1}),3); end
    panel_files = cell(F,1);
    for f = 1:F
        M = numel(names) + (~isempty(OmTrue));
        Scells = cell(M,1); labels = cell(M,1); vmax = 0;
        pos = 1;
        if ~isempty(OmTrue)
            S = pcoh_scores(OmTrue(:,:,f));
            S(1:size(S,1)+1:end) = 0; v = abs(S(:)); v = v(~isnan(v));
            if ~isempty(v)
                q = quantile(v, eval_cfg.panel_clip); S = max(min(S, q), 0); vmax = max(vmax, q);
            end
            Scells{pos} = S; labels{pos} = 'Ground Truth'; pos = pos + 1;
        end
        for i=1:numel(names)
            K_stack = coerce_precision_stack(Omega_all.(names{i}));
            S = pcoh_scores(K_stack(:,:,f));
            S(1:size(S,1)+1:end) = 0; v = abs(S(:)); v = v(~isnan(v));
            if ~isempty(v)
                q = quantile(v, eval_cfg.panel_clip); S = max(min(S, q), 0); vmax = max(vmax, q);
            end
            Scells{pos} = S; labels{pos} = names{i}; pos = pos + 1;
        end
        if ~isempty(eval_cfg.panel_caxis), vmax = eval_cfg.panel_caxis(2); end
        figure('Color','w','Name',sprintf('PCoh panel (f=%d)', f),'Position',[100 100 260*M 220]);
        tiledlayout(1,M,'Padding','compact','TileSpacing','compact');
        for i=1:M
            nexttile; imagesc(Scells{i}); axis image off;
            title(strrep(labels{i},'_','\_'),'Interpreter','none','FontSize',10);
            if ~isempty(vmax), caxis([0 vmax]); end
        end
        colormap(parula); cb = colorbar; cb.Layout.Tile = 'east'; cb.Label.String = '|partial coherence|';
        file_i = fullfile(outdir, sprintf('pcoh_panel_f%02d.png', f));
        saveas(gcf, file_i); panel_files{f} = file_i; close(gcf);
    end
end

function Sii_proxy = make_source_cov_proxy(data, eval_cfg)
% GLS whitening + minimum-norm inverse to build a proxy of source-domain
% sample covariance per frequency. Falls back gracefully if pieces missing.
    F = numel(data.emp_cov_cell);
    L = getd(data,'L',[]); n = []; p = [];
    if ~isempty(L), [n,p] = size(L); end
    if isempty(L) || ~isnumeric(L) || any(size(L)<1)
        warning('No leadfield L provided; using identity proxy in source-space.');
        p = size(data.emp_cov_cell{1},1); L = eye(p); n = p;
    end
    Wcell = cell(F,1);
    bt = getd(data,'bayes_truth',struct());
    switch lower(getd(bt,'noise_type','scalar'))
        case 'matrix_per_freq'
            if isfield(bt,'Sigma_xixi_per_freq_true') && ~isempty(bt.Sigma_xixi_per_freq_true)
                for f=1:F, Wcell{f} = inv_chol(bt.Sigma_xixi_per_freq_true{f}); end
            else
                for f=1:F, Wcell{f} = eye(n); end
            end
        case 'matrix'
            if isfield(bt,'Sigma_xixi_true') && ~isempty(bt.Sigma_xixi_true)
                W = inv_chol(bt.Sigma_xixi_true); for f=1:F, Wcell{f} = W; end
            else
                for f=1:F, Wcell{f} = eye(n); end
            end
        otherwise
            sig2 = getd(bt,'sigma_xi2_true',1); W = eye(n)/sqrt(max(sig2,eps)); for f=1:F, Wcell{f} = W; end
    end
    Sii_proxy = zeros(p,p,F);
    lsq_ridge = getd(getd(eval_cfg,'adapter',struct()),'lsq_ridge',5e-3);
    cov_ridge = getd(getd(eval_cfg,'adapter',struct()),'cov_ridge',1e-6);
    for f=1:F
        Svv = hermitizeC(data.emp_cov_cell{f}); W = Wcell{f}; Lw = W * L; Svvw = hermitizeC(W * Svv * W');
        G = (Lw'*Lw + lsq_ridge*eye(p)) \ (Lw');
        Sii_proxy(:,:,f) = hermitizeC(G * Svvw * G') + cov_ridge*eye(p);
    end
end

function [best_tau, tau_grid, ll_curve] = likelihood_scan_tau(OmEst, Sii_proxy, tau_scan, adapter)
% For each relative tau, shrink off-diagonals in normalized magnitude units
% and sum source-domain log-likelihood across frequencies.
    [p,~,F] = size(OmEst); %#ok<ASGLU>
    up = triu(true(p),1); %#ok<NASGU>
    meds = zeros(F,1);
    for f=1:F
        S = pcoh_scores(OmEst(:,:,f)); v = abs(S(up));
        m = median(v); mad = median(abs(v - m)) / 0.6745;
        meds(f) = max(mad, 1e-6);
    end
    tau_grid = tau_scan(:); ll_curve = zeros(numel(tau_grid),1);
    for it = 1:numel(tau_grid)
        tau_rel = tau_grid(it); ll_sum = 0;
        for f=1:F
            K0 = hermitizeC(OmEst(:,:,f)); tau_eff = tau_rel * meds(f);
            K = shrink_offdiag(K0, tau_eff, adapter);
            Sf = hermitizeC(Sii_proxy(:,:,f));
            ll_sum = ll_sum + logdet_spd(K) - real(trace(Sf*K));
        end
        ll_curve(it) = ll_sum;
    end
    [~,idx] = max(ll_curve); best_tau = tau_grid(idx);
end

% --------------------- SPD distances & summaries ----------------------
function D = spd_distance_pack(Astack, Bstack)
    [~,~,F] = size(Astack);
    D = struct('KLS',zeros(F,1),'LE',zeros(F,1),'AIRM',zeros(F,1),'JBLD',zeros(F,1),'BW',zeros(F,1));
    for f=1:F
        A = spd_project_small(hermitizeC(Astack(:,:,f)));
        B = spd_project_small(hermitizeC(Bstack(:,:,f)));
        D.KLS(f)  = 0.5*( trace(B\A) + trace(A\B) - 2*size(A,1) );
        LA = logm(A); LB = logm(B); D.LE(f) = norm(LA - LB, 'fro');
        C = spd_project_small(A^(-1/2) * B * A^(-1/2)); D.AIRM(f) = norm(logm(C), 'fro');
        D.JBLD(f) = logdet_spd((A+B)/2) - 0.5*logdet_spd(A) - 0.5*logdet_spd(B);
        R = spd_project_small(A^(1/2) * B * A^(1/2)); D.BW(f) = sqrt( trace(A) + trace(B) - 2*trace(real(sqrtm(R))) );
    end
end

function S = summarize_distances(D)
    fn = fieldnames(D); S = struct();
    for i=1:numel(fn)
        x = D.(fn{i}); S.(fn{i}) = struct('mean',mean(x),'median',median(x),'std',std(x),'per_freq',x(:).');
    end
end

function write_distance_csv(path, dist_summary)
    methods = fieldnames(dist_summary); metrics = {'KLS','LE','AIRM','JBLD','BW'};
    fid = fopen(path,'w'); fprintf(fid, 'method,metric,mean,median,std\n');
    for i=1:numel(methods)
        m = methods{i}; S = dist_summary.(m);
        for j=1:numel(metrics)
            mt = metrics{j}; if isfield(S, mt)
                fprintf(fid, '%s,%s,%.6g,%.6g,%.6g\n', m, mt, S.(mt).mean, S.(mt).median, S.(mt).std);
            end
        end
    end
    fclose(fid);
end

% ----------------------- primitives & numerics ------------------------
function X3 = coerce_precision_stack(X)
    if iscell(X)
        F = numel(X); p = size(X{1},1); X3 = zeros(p,p,F,class(X{1}));
        for f = 1:F, X3(:,:,f) = hermitizeC(X{f}); end
    elseif isnumeric(X) && ndims(X)==3 && size(X,1)==size(X,2)
        X3 = X; for f=1:size(X3,3), X3(:,:,f) = hermitizeC(X3(:,:,f)); end
    else
        error('Unsupported precision format. Expect {F×1} or p×p×F.');
    end
end

function S = pcoh_scores(K)
    K = hermitizeC(K); d = real(diag(K)); d(d<=0) = eps; s = sqrt(d); D12 = s * s.';
    S = abs(K) ./ max(D12, eps); S(1:size(S,1)+1:end) = 0; S = hermitizeC(S);
end

function K = shrink_offdiag(K, tau, adapter)
    K = hermitizeC(K); p = size(K,1); d = real(diag(K)); d(d<=0)=eps; D12 = sqrt(d)*sqrt(d).';
    Mag = abs(K); Phs = complex(real(K), imag(K)) ./ max(Mag, eps);
    idx = triu(true(p),1); v = Mag(idx) - tau * D12(idx); v = max(v, 0);
    K(idx) = Phs(idx) .* v; K = hermitizeC(K);
    
    spd_delta = getd(adapter,'spd_delta',0.03); K = spd_project_small(K) + spd_delta*eye(p);
end

function S = inv_psd(K)
    K = hermitizeC(K); [V,D] = eig(K,'vector'); D = real(D);
    floorv = 1e-9*median(abs(D)+eps); if ~isfinite(floorv) || floorv<=0, floorv = 1e-9; end
    D(D < floorv) = floorv; S = V*diag(1./D)*V'; S = hermitizeC(S);
end

function s = logdet_spd(A)
    A = spd_project_small(A);
    try, C = chol(hermitizeC(A), 'lower'); s = 2*sum(log(diag(C)));
    catch, [V,D] = eig(hermitizeC(A), 'vector'); D = real(D); D(D<1e-12)=1e-12; s = sum(log(D)); end
end

function A = hermitizeC(A)
    A = (A + A')/2;  % Hermitian symmetrization
end

function K = spd_project_small(K)
    K = hermitizeC(K); [V,D] = eig(K,'vector'); D = real(D);
    floorv = max(1e-9, 1e-9*median(abs(D)+eps)); D(D<floorv) = floorv;
    K = V*diag(D)*V'; K = hermitizeC(K);
end

function W = inv_chol(S)
    S = hermitizeC(S); [V,D] = eig(S,'vector'); D = real(D); D(D<1e-12) = 1e-12; S = V*diag(D)*V';
    C = chol(S,'lower'); W = C \ eye(size(S));
end

function val = getd(s, name, def)
    if isstruct(s) && isfield(s,name) && ~isempty(s.(name)), val = s.(name); else, val = def; end
end

function val = set_default(S, D)
    val = S; f = fieldnames(D);
    for i=1:numel(f)
        if ~isfield(val,f{i}) || isempty(val.(f{i})), val.(f{i}) = D.(f{i}); end
    end
end

function nm2 = make_valid_name(nm)
    if exist('matlab.lang.makeValidName','file')
        nm2 = matlab.lang.makeValidName(nm);
    else
        nm2 = regexprep(nm,'[^a-zA-Z0-9_]','_');
        if isempty(nm2) || (~isletter(nm2(1)) && nm2(1)~='_'), nm2 = ['x_' nm2]; end
    end
end
