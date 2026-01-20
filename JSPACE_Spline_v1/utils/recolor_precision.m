function Omega_cell = recolor_precision(Gamma_cell, D_cell)
%RECOLOR_PRECISION  Convert whitened precision to original space.
%
% Formula: Omega = D * Gamma * D
%
% Inputs:
%   Gamma_cell : {F x 1}, whitened precision matrices
%   D_cell     : {F x 1}, diagonal whitening matrices
%
% Output:
%   Omega_cell : {F x 1}, recolored precision matrices

    F = numel(Gamma_cell);
    if numel(D_cell) ~= F
        error('recolor_precision:SizeMismatch', ...
            'Gamma_cell and D_cell must have the same length.');
    end

    Omega_cell = cell(F, 1);

    for t = 1:F
        G = Gamma_cell{t};
        D = D_cell{t};
        N = size(G, 1);

        if isvector(D)
            D = diag(D);
        end

        Om = D * G * D;
        Om = 0.5 * (Om + Om');
        Om(1:N+1:end) = real(diag(Om));

        Omega_cell{t} = Om;
    end
end
