function [invA, logdetA, is_spd] = spd_inv_and_logdet(A)
%SPD_INV_AND_LOGDET  Robust inverse and log-determinant for SPD matrices.
%
% Strategy:
%   1) Enforce Hermitian symmetry.
%   2) Attempt Cholesky.
%   3) If it fails, add diagonal jitter and retry with increasing scale.
%
% Inputs:
%   A : (N x N) Hermitian matrix
%
% Outputs:
%   invA    : inverse of A (NaN if failure)
%   logdetA : log(det(A)) (-inf if failure)
%   is_spd  : true if Cholesky succeeds after jittering

    if ndims(A) ~= 2 || size(A, 1) ~= size(A, 2)
        error('spd_inv_and_logdet:SizeMismatch', 'Input must be a square matrix.');
    end

    A = 0.5 * (A + A');
    N = size(A, 1);
    A(1:N+1:end) = real(diag(A));

    [R, p] = chol(A);
    jitter = 0;

    if p > 0
        mean_diag = mean(real(diag(A)));
        if ~isfinite(mean_diag) || mean_diag <= 0
            mean_diag = 1.0;
        end
        jitter = 1e-6 * mean_diag;

        max_tries = 6;
        for k = 1:max_tries
            [R, p] = chol(A + jitter * eye(N, 'like', A));
            if p == 0
                break;
            end
            jitter = jitter * 10;
        end
    end

    if p > 0
        invA = NaN(N, N, 'like', A);
        logdetA = -inf;
        is_spd = false;
        return;
    end

    d = real(diag(R));
    logdetA = 2 * sum(log(d));

    I = eye(N, 'like', A);
    invA = R \ (R' \ I);
    is_spd = true;
end
