function estep_out = module2_estep(estep_in, opts)
% MODULE2_ESTEP
% Wrapper for E-step with centralized numerical safeguards:
%   - Hermitian enforcement (inputs & outputs)
%   - Optional diagonal loading on outputs
%   - PSD enforcement via Cholesky with adaptive jitter
%
% This wrapper delegates the *core* E-step computation to
%   module2_estep_main(estep_in, core_opts)
% and then applies robust post-processing to the returned
% source second moments (Sjj_hat).
%
% Inputs:
%   estep_in : struct
%       .leadfield_matrix          (p×n)
%       .empirical_covariances     {F×1} cell of p×p
%       .source_prior_covariances  {F×1} cell of n×n
%       .noise_covariance          (p×p)
%       .frequencies               1×F
%   opts : struct (optional)
%       % numerical safety (defaults chosen to be robust)
%       .ensure_hermitian   (default: true)
%       .ensure_real_diag   (default: true)
%       .ensure_psd         (default: true)
%       .psd_tol            (default: 1e-10)
%       .diag_loading       (default: 1e-8)   % applied to outputs only
%       .max_jitter_tries   (default: 5)      % additional jitter attempts
%       .verbose            (default: false)
%       % passthrough to core (use only if your main accepts them)
%       .core_opts          (default: struct())
%
% Output:
%   estep_out : struct
%       .source_second_moments  {F×1} cell of n×n (Hermitian, PSD, loaded)
%       (other fields returned by module2_estep_main are preserved)
%
% Author: your team
% Date: 2025-09-09
% Version: 1.0

if nargin < 2, opts = struct(); end

% -------- defaults --------
d.ensure_hermitian = true;
d.ensure_real_diag = true;
d.ensure_psd       = true;
d.psd_tol          = 1e-10;
d.diag_loading     = 1e-8;   % unified loading place (outputs)
d.max_jitter_tries = 5;
d.verbose          = false;
d.core_opts        = struct();   % passed to *_main if needed

opts = set_defaults_(opts, d);

% -------- light input sanitation (Hermitian only, no loading) --------
if isfield(estep_in,'empirical_covariances') && iscell(estep_in.empirical_covariances)
    F = numel(estep_in.empirical_covariances);
    for f=1:F
        C = estep_in.empirical_covariances{f};
        if opts.ensure_hermitian
            C = 0.5*(C + C');
        end
        if opts.ensure_real_diag
            C(1:size(C,1)+1:end) = real(diag(C));
        end
        estep_in.empirical_covariances{f} = C;
    end
end

% -------- call core E-step (your existing algorithm) --------
% If your repo implements the logic in this very function, replace the line
% below with your original core. Otherwise, keep this delegation.
core_opts = opts.core_opts;
raw_out = module2_estep_main(estep_in, core_opts);

% -------- post-process outputs: central PSD & loading --------
estep_out = raw_out;
if ~isfield(raw_out,'source_second_moments') || ~iscell(raw_out.source_second_moments)
    error('module2_estep:missing_field', ...
        'module2_estep_main must return .source_second_moments as {F×1} cell.');
end

F = numel(raw_out.source_second_moments);
S_out = raw_out.source_second_moments;

for f = 1:F
    S = S_out{f};
    % unify hermitian/real-diag
    if opts.ensure_hermitian
        S = 0.5*(S + S');
    end
    if opts.ensure_real_diag
        S(1:size(S,1)+1:end) = real(diag(S));
    end

    if opts.ensure_psd
        % apply base loading
        if opts.diag_loading > 0
            S = S + opts.diag_loading * eye(size(S,1), 'like', S);
        end
        % check PSD via Cholesky; adaptively jitter if needed
        [R,flag] = chol(0.5*(S+S'), 'lower'); %#ok<ASGLU>
        jitter = opts.diag_loading;
        tries = 0;
        while flag ~= 0 && tries < opts.max_jitter_tries
            jitter = max(jitter*10, eps);
            S = 0.5*(S+S') + jitter * eye(size(S,1), 'like', S);
            [R,flag] = chol(S, 'lower'); %#ok<ASGLU>
            tries = tries + 1;
        end
        if flag ~= 0
            error('module2_estep:psd_enforce_failed', ...
                'Failed to enforce PSD after %d jitter tries.', opts.max_jitter_tries);
        end
        % optional tiny negative cleanup by eig thresholding
        if opts.psd_tol > 0
            [V,D] = eig(0.5*(S+S'));
            ev = real(diag(D));
            ev(ev < opts.psd_tol) = opts.psd_tol;
            S = V*diag(ev)*V';
            S = 0.5*(S+S');
        end
    end

    S_out{f} = S;
end

estep_out.source_second_moments = S_out;

if opts.verbose
    fprintf('[module2_estep] PSD/Loading applied: diag_loading=%.2e, tries<=%d\n', ...
        opts.diag_loading, opts.max_jitter_tries);
end

end

% ------------------------------- helpers ----------------------------------
function S = set_defaults_(S, D)
    f = fieldnames(D);
    for i=1:numel(f)
        k = f{i};
        if ~isfield(S,k) || isempty(S.(k)), S.(k) = D.(k); end
    end
end
