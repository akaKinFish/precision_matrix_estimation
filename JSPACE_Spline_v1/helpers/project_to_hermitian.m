function Theta_cell = project_to_hermitian(Theta_cell)
%PROJECT_TO_HERMITIAN  Enforce Hermitian symmetry and real diagonals.
%
% T <- (T + T') / 2
% diag(T) <- real(diag(T))

    K = numel(Theta_cell);
    N = size(Theta_cell{1}, 1);

    for k = 1:K
        T = Theta_cell{k};
        T = 0.5 * (T + T');
        T(1:N+1:end) = real(diag(T));
        Theta_cell{k} = T;
    end
end
