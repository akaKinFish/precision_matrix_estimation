function warm = eloreta_warmstart_from_covs(Svv_cell, L, gamma_grid, opts)
% ELORETA_WARMSTART_FROM_COVS
% Inputs:
%   - Svv_cell : {F×1}, each p×p (complex Hermitian) sensor cross-spectrum
%   - L        : p×n leadfield (real or complex)
%   - gamma_grid : vector of ridge scales (e.g. logspace(-4,1,30))
%   - opts     : struct with fields:
%                  .maxit   (default 50)   eLORETA W-update max iters
%                  .tol     (1e-6)         relative tol for W update
%                  .verbose (false)
%
% Output struct 'warm' fields (all {F×1} cells):
%   .Omega_init{f}  : n×n  initial precision  ≈ inv_psd_robust(Sjj_e{f})
%   .Sjj_e{f}       : n×n  eLORETA reconstructed source covariance T*Svv*T'
%   .Tjv{f}         : n×p  eLORETA filter
%   .gamma_opt{f}   : scalar (chosen from gamma_grid by GCV)
%   .gcv_curve{f}   : vector of scanned GCV values

    if nargin < 4, opts = struct(); end
    maxit = getf(opts,'maxit',50);
    tol   = getf(opts,'tol',1e-6);
    verb  = getf(opts,'verbose',false);

    F = numel(Svv_cell);
    [p,n] = size(L);

    Tjv_cell   = cell(F,1);
    Sjj_cell   = cell(F,1);
    Om0_cell   = cell(F,1);
    gcv_cell   = cell(F,1);
    gopt_cell  = cell(F,1);

    for f = 1:F
        Svv = Svv_cell{f};
        % ---- run single-frequency eLORETA with GCV over gamma_grid
        [Tjv, Sjj, ~, gamma_opt, gcv] = eloreta_simple_(Svv, L, gamma_grid, maxit, tol, verb);

        % Hermitianize + PSD project (for numerical safety)
        Sjj = psd_project_(Sjj);

        % robust inverse as initial precision
        Om0 = inv_psd_robust_(Sjj, 1e-8, 1e-12);

        % store
        Tjv_cell{f}  = Tjv;
        Sjj_cell{f}  = Sjj;
        Om0_cell{f}  = Om0;
        gcv_cell{f}  = gcv;
        gopt_cell{f} = gamma_opt;
    end

    warm = struct('Omega_init',{Om0_cell}, ...
                  'Sjj_e',{Sjj_cell}, ...
                  'Tjv',{Tjv_cell}, ...
                  'gamma_opt',{gopt_cell}, ...
                  'gcv_curve',{gcv_cell});
end

% ---------- helpers ----------

function [Tjv, Sjj, W, gamma_opt, gcv] = eloreta_simple_(Svv, L, gamma_grid, maxit, tol, verb)
% Single-frequency eLORETA with GCV over gamma_grid (scalar orientation)

    [p,n] = size(L);
    gcv = zeros(numel(gamma_grid),1);

    best.T = []; best.W = []; best.gamma = NaN; best.score = Inf;

    for k = 1:numel(gamma_grid)
        gamma = gamma_grid(k);                % ridge scale
        % --- inner eLORETA loop (update W, then T)
        w = ones(n,1);                        % W = diag(w)
        for it = 1:maxit
            Winv = diag(1./w);
            A    = hermi_( L*Winv*L' );
            alpha= gamma * trace(A)/p;        % alpha normalized by mean eig
            M    = inv_psd_( A + alpha*eye(p) );

            w_old = w;
            % update each voxel weight: w_i = sqrt( l_i^H M l_i )
            for i=1:n
                li  = L(:,i);
                mii = real(li' * M * li);
                w(i)= sqrt(max(mii, eps));
            end
            if norm(w-w_old)/max(1,norm(w_old)) < tol, break; end
        end
        Winv = diag(1./w);
        T    = Winv * L' * M;                 % n×p: W^{-1} L^H (A+alpha I)^{-1}

        % --- GCV score
        Txiv = eye(p) - L*T;                  % I - L T
        num  = real(trace( hermi_(Txiv*Svv*Txiv') ))/p;
        den  = ( real(trace(Txiv))/p )^2 + eps;
        gcv(k) = num / den;

        if gcv(k) < best.score
            best.T = T; best.W = diag(w); best.gamma = gamma; best.score = gcv(k);
        end

        if verb && (mod(k,10)==1)
            fprintf('[eLORETA] gamma=%.3g, GCV=%.3g\n', gamma, gcv(k));
        end
    end

    Tjv       = best.T;
    gamma_opt = best.gamma;
    W         = best.W;
    Sjj       = Tjv * Svv * Tjv';
end

function A = hermi_(A), A = (A + A')/2; end
function X = inv_psd_(X)
    X = hermi_(X);
    [U,S] = eig(X,'vector'); S = max(real(S), eps); X = U*diag(1./S)*U'; X = hermi_(X);
end
function S = psd_project_(S)
    S = hermi_(S);
    [U,d] = eig(S,'vector'); d = max(real(d), 0); S = U*diag(d)*U'; S = hermi_(S);
end
function v = getf(s, f, d), if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end, end
function Om = inv_psd_robust_(A, eps_reg, min_ratio)
    A = hermi_(A);
    [V,D] = eig(A); d = real(diag(D));
    dmax  = max(d); floor_val = max(min_ratio * max(dmax, eps), 0);
    d(d<floor_val) = floor_val;
    if eps_reg > 0, d = (d + eps_reg * dmax) / (1 + eps_reg); end
    Om = V * diag(1./d) * V'; Om = hermi_(Om);
end
