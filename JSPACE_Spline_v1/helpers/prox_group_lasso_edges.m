function Theta_out = prox_group_lasso_edges(Theta_in, lambda1, alpha, w_edge)
%PROX_GROUP_LASSO_EDGES  Apply proximal operator for group lasso penalty.
%
% Optimization problem per edge (i,j):
%   min_{v}  0.5 ||v - z||_2^2 + (alpha * lambda1 * w_ij) ||v||_2
% Solution: Block Soft Thresholding.
%
% Inputs:
%   Theta_in : {K x 1} Current gradients step result
%   lambda1  : Group sparsity strength
%   alpha    : Current step size (from FISTA)
%   w_edge   : (N x N) Edge-specific weights
%
% Outputs:
%   Theta_out: {K x 1} Coefficient matrices after shrinkage

    Theta_out = Theta_in;

    if lambda1 <= 1e-12 || alpha <= 0
        return;
    end

    N = size(Theta_in{1}, 1);
    K = numel(Theta_in);

    if nargin < 4 || isempty(w_edge)
        w_edge = ones(N);
    end

    % Loop over upper triangle only
    for i = 1:N
        for j = i+1:N
            v = zeros(K, 1);
            for k = 1:K
                v(k) = Theta_in{k}(i, j);
            end

            tau = alpha * lambda1 * w_edge(i, j);
            if tau < 0
                tau = 0;
            end

            norm_v = norm(v, 2);
            if norm_v <= tau
                for k = 1:K
                    Theta_out{k}(i, j) = 0;
                    Theta_out{k}(j, i) = 0;
                end
            else
                scale = 1 - (tau / norm_v);
                v_new = scale * v;
                for k = 1:K
                    Theta_out{k}(i, j) = v_new(k);
                    Theta_out{k}(j, i) = conj(v_new(k));
                end
            end
        end
    end
end
