function [Sjj_hat, Sigma_post, Tjv, sigma2xi_used, info] = bcv_estep_scalar_noise(Svv_cell, L, sigma2xi, opts)
% BC-V style E-step with scalar sensor noise (σ^2 I)
% Inputs
%   Svv_cell : {F×1} sensor covs
%   L        : p×n leadfield
%   sigma2xi : scalar noise variance (if empty, auto from Svv)
%   opts.use_subspace    (bool, default true)
%   opts.subspace_energy (0-1,  default 0.99)
%   opts.do_scale_L      (bool, default true)   % 若缩放 L，会自动把 σ^2 做 1/c^2 补偿
%   opts.Sigma_prior     (cell|[]) prior covariances, default I
%
% Outputs
%   Sjj_hat   : {F×1} posterior source second moments
%   Sigma_post: {F×1} posterior covariances
%   Tjv       : {F×1} linear estimator
%   sigma2xi_used: final scalar σ^2 after scale compensation
%   info      : struct (Ur,Sr,Vr,scaleL,...)

if nargin<4, opts=struct(); end
Svv_cell = coerce_cov_cell_(Svv_cell);
F = numel(Svv_cell);
[p,n] = size(L);

use_sub = getf_(opts,'use_subspace',true);
eng     = getf_(opts,'subspace_energy',0.99);
scaleL  = getf_(opts,'do_scale_L',true);
Sigma_prior = getf_(opts,'Sigma_prior',[]);

if isempty(sigma2xi)
    % 一个稳健的初值：trace(Svv)/p 的 1%~2%
    s0 = 0;
    for f=1:F, s0 = s0 + trace(Svv_cell{f}); end
    s0 = s0 / (F*p);
    sigma2xi = max(1e-8, 1e-2*s0);
end

% 可选缩放 L，并对 σ^2 做等效补偿
scaleLvj = 1.0;
if scaleL
    scaleLvj = sqrt(trace(L*L')/p);
    if isfinite(scaleLvj) && scaleLvj>0
        L = L/scaleLvj;
        sigma2xi = sigma2xi/(scaleLvj^2);
    end
end

% 子空间（VARETA 风格）
Ur=[]; Sr=[]; Vr=[];
if use_sub
    [U,Sv,V] = svd(L,'econ');
    sing2 = diag(Sv).^2; cum = cumsum(sing2)/sum(sing2);
    r = find(cum>=eng,1,'first'); if isempty(r), r=size(Sv,1); end
    Ur=U(:,1:r); Sr=Sv(1:r,1:r); Vr=V(:,1:r);
end

Sjj_hat   = cell(F,1);
Sigma_post= cell(F,1);
Tjv       = cell(F,1);

I = eye(p);
for f=1:F
    Svv = 0.5*(Svv_cell{f}+Svv_cell{f}');
    if isempty(Sigma_prior)
        Om_prior = eye(n);   % Θ_prior = I
    else
        Sigp = 0.5*(Sigma_prior{f}+Sigma_prior{f}');
        Om_prior = inv_psd_(Sigp);
    end

    if ~isempty(Ur)      % 子空间加速：L'Σ^{-1}L = (1/σ^2) V_r S_r^2 V_r'
        A = Om_prior + (1/sigma2xi)*(Vr*(Sr*Sr)*Vr');
    else
        A = Om_prior + (1/sigma2xi)*(L'*L);
    end
    A = 0.5*(A+A');
    Sig_post = inv_psd_(A);

    T = Sig_post * (L'*(1/sigma2xi));     % = Σpost L' Σ^{-1}
    Sjj = 0.5*(Sig_post + T*Svv*T');      % = Σpost + T Svv T'

    Sjj_hat{f}    = Sjj;
    Sigma_post{f} = Sig_post;
    Tjv{f}        = T;
end

sigma2xi_used = sigma2xi;
info = struct('scaleL',scaleLvj,'Ur',Ur,'Sr',Sr,'Vr',Vr);

% ---- helpers
function C = coerce_cov_cell_(X)
    if iscell(X), C=X(:); return; end
    if isnumeric(X)&&ndims(X)==3&&size(X,1)==size(X,2)
        F=size(X,3); C=cell(F,1); for k=1:F, C{k}=X(:,:,k); end, return;
    end
    if isnumeric(X)&&ismatrix(X)&&size(X,1)==size(X,2), C={X}; return; end
    error('Svv must be {F×1} or p×p×F or single p×p');
end
function v = getf_(s,f,d), if isfield(s,f)&&~isempty(s.(f)), v=s.(f); else, v=d; end, end
function X = inv_psd_(X)
    X = 0.5*(X+X'); [U,S]=eig(X,'vector'); S = max(real(S),1e-10); X=U*diag(1./S)*U'; X=0.5*(X+X');
end
end
