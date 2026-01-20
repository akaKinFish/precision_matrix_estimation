function Svv_out = sanitize_cross_spectrum(Svv_in, svv_cfg)
%SANITIZE_CROSS_SPECTRUM  Clean sensor cross-spectrum matrices.
%
% Behavior:
%   1) Enforce Hermitian symmetry: S = (S + S')/2
%   2) Enforce real diagonal
%   3) Optional diagonal jitter for stability
%
% Inputs:
%   Svv_in  : (Ne x Ne x F) complex tensor
%   svv_cfg : struct with optional field .jitter (scalar or length F)
%
% Output:
%   Svv_out : cleaned tensor

    if nargin < 2 || isempty(svv_cfg)
        svv_cfg = struct();
    end

    jitter = 0;
    if isfield(svv_cfg, 'jitter') && ~isempty(svv_cfg.jitter)
        jitter = svv_cfg.jitter;
    end

    is2d = (ndims(Svv_in) == 2);
    if is2d
        Svv_in = reshape(Svv_in, size(Svv_in, 1), size(Svv_in, 2), 1);
    end

    [Ne1, Ne2, F] = size(Svv_in);
    if Ne1 ~= Ne2
        error('sanitize_cross_spectrum:SizeMismatch', ...
            'Svv must be square per frequency.');
    end

    Svv_out = zeros(size(Svv_in), 'like', Svv_in);

    for t = 1:F
        S = Svv_in(:, :, t);

        % Enforce Hermitian symmetry.
        S = 0.5 * (S + S');

        % Enforce real diagonal.
        d = diag(S);
        S(1:Ne1+1:end) = real(d);

        % Optional diagonal jitter.
        if ~isempty(jitter)
            if numel(jitter) > 1
                jt = jitter(min(t, numel(jitter)));
            else
                jt = jitter;
            end
            if jt > 0
                S = S + jt * eye(Ne1, 'like', S);
            end
        end

        Svv_out(:, :, t) = S;
    end

    if is2d
        Svv_out = Svv_out(:, :, 1);
    end
end
