function [L_out, info] = normalize_graph_laplacian(L_in, mode, opts)
% NORMALIZE_GRAPH_LAPLACIAN
% Spectral or symmetric normalization for a (real, symmetric) Laplacian.
%
%   [L_out, info] = normalize_graph_laplacian(L_in, mode, opts)
%
% Inputs:
%   L_in : n×n real symmetric (approximately) Laplacian (unnormalized)
%   mode : 'spectral' (default) | 'sym'
%          'spectral' : L / lambda_max(L)
%          'sym'      : I - D^(-1/2) A D^(-1/2), where L = D - A (A>=0, diag=0)
%   opts : struct (optional)
%       .clip_negative_eigs : logical, default true (clip tiny negatives)
%       .eps_eig            : double, default 0
%
% Outputs:
%   L_out : normalized Laplacian (symmetric PSD)
%   info  : struct with eigen stats & symmetry error
%
% Notes:
%   - We symmetrize input with 0.5*(L+L'), and zero its tiny imag part if any.
%   - For 'sym' mode we reconstruct A=D-L and apply normalized Laplacian formula.
%
% Author: your team

    if nargin < 2 || isempty(mode), mode = 'spectral'; end
    if nargin < 3, opts = struct(); end
    opts = set_defaults_(opts, struct('clip_negative_eigs', true, 'eps_eig', 0));

    L = real(0.5*(L_in + L_in'));
    n = size(L,1);

    info = struct();
    info.symerr = norm(L - L','fro');
    ev = eig(L);
    info.min_eig_before = min(real(ev));
    info.max_eig_before = max(real(ev));

    switch lower(mode)
        case 'spectral'
            % L_bar = L / lambda_max(L)
            lam_max = max(real(ev));
            if lam_max <= 0
                error('normalize_graph_laplacian:spectral','Non-positive max eigenvalue.');
            end
            L_out = L / lam_max;

        case 'sym'
            % L = D - A  =>  A = D - L
            d  = diag(L);
            A  = diag(d) - L;
            A  = max(A, 0);        % no negative weights
            A(1:n+1:end) = 0;
            A  = 0.5*(A + A');     % sym
            deg = sum(A,2);
            Dm12 = diag(1 ./ max(sqrt(deg), eps));
            L_out = eye(n) - Dm12 * A * Dm12;

        otherwise
            error('normalize_graph_laplacian:mode','Unsupported mode: %s', mode);
    end

    % Optional PSD projection (clip tiny negatives)
    if opts.clip_negative_eigs
        [V,D] = eig(0.5*(L_out+L_out'));
        dvals = real(diag(D));
        dvals = max(dvals, opts.eps_eig);
        L_out = V*diag(dvals)*V';
    else
        L_out = 0.5*(L_out + L_out');
    end

    ev2 = eig(L_out);
    info.min_eig_after = min(real(ev2));
    info.max_eig_after = max(real(ev2));
    info.symerr_after  = norm(L_out - L_out','fro');
end

function S = set_defaults_(S, D)
    f = fieldnames(D);
    for i = 1:numel(f)
        k = f{i};
        if ~isfield(S,k) || isempty(S.(k)), S.(k) = D.(k); end
    end
end
