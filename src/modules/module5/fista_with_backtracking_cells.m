function [x, history, smoothf] = fista_with_backtracking_cells(f, grad_f, g, prox, x0, lambda, L0, eta, max_iter, tol, max_backtracking_iter)
% FISTA_WITH_BACKTRACKING_CELLS - FISTA with backtracking for cell-of-matrices variables.
%
% Inputs mirror the vector version but accept cells:
%   f        : smooth objective handle, f(x_cells) -> scalar
%   grad_f   : gradient handle, grad_f(x_cells) -> cell of same shape
%   g        : nonsmooth handle returning scalar (unscaled by lambda2)
%   prox     : proximal operator handle, prox(x_cells, tau) -> cell
%   x0       : cell initial iterate
%   lambda   : regularization scalar (multiplies g and prox tau)
%   L0       : initial Lipschitz estimate (default 1)
%   eta      : backtracking factor (default 2)
%   max_iter : max iterations (default 200)
%   tol      : relative objective tolerance (default 1e-6)
%   max_backtracking_iter : max backtracking steps (default 30)
%
% Outputs:
%   x        : final cell iterate
%   history  : vector of objective values f + lambda*g
%   smoothf  : final smooth part f(x)
%
% Notes:
%   - Uses quadratic upper bound with cell inner products.
%   - Keeps restart safeguard if objective increases.
%   - Assumes prox is pure (no gradient step inside).

    if nargin < 11 || isempty(max_backtracking_iter), max_backtracking_iter = 30; end
    if nargin < 10 || isempty(tol), tol = 1e-6; end
    if nargin < 9  || isempty(max_iter), max_iter = 200; end
    if nargin < 8  || isempty(eta), eta = 2; end
    if nargin < 7  || isempty(L0), L0 = 1; end

    x = x0;
    y = x0;
    t = 1;
    L = L0;
    history = zeros(max_iter, 1);
    f_old = f(x) + lambda * g(x);

    for k = 1:max_iter
        % Backtracking line search
        found_L = false;
        backtracking_iter = 0;
        fy = f(y);
        gy = g(y);
        grad_fy = grad_f(y);

        while ~found_L && backtracking_iter < max_backtracking_iter
            L_bar = L * eta^backtracking_iter;
            step = cell_sub(y, cell_scale(grad_fy, 1 / L_bar));
            x_new = prox(step, lambda / L_bar);

            fx_new = f(x_new);
            gx_new = g(x_new);

            delta = cell_sub(x_new, y);
            Q_L = fy + cell_inner(grad_fy, delta) + (L_bar / 2) * cell_norm_sq(delta) + lambda * gx_new;

            if fx_new + lambda * gx_new <= Q_L
                found_L = true;
                L = L_bar;
            else
                backtracking_iter = backtracking_iter + 1;
            end
        end

        f_new = fx_new + lambda * gx_new;

        % Restart if objective increased
        if f_new > f_old
            t = 1;
            y = x;
        else
            t_old = t;
            t = (1 + sqrt(1 + 4 * t^2)) / 2;
            y = cell_add(x_new, cell_scale(cell_sub(x_new, x), (t_old - 1) / t));
            x = x_new;
        end

        history(k) = f_new;

        if abs(f_new - f_old) / max(1, abs(f_old)) < tol
            history = history(1:k);
            smoothf = f(x);
            return;
        end
        f_old = f_new;
    end

    history = history(1:max_iter);
    smoothf = f(x);
end

% -------- cell helpers --------
function C = cell_add(A, B)
    F = numel(A); C = cell(F,1);
    for f = 1:F, C{f} = A{f} + B{f}; end
end

function C = cell_sub(A, B)
    F = numel(A); C = cell(F,1);
    for f = 1:F, C{f} = A{f} - B{f}; end
end

function C = cell_scale(A, s)
    F = numel(A); C = cell(F,1);
    for f = 1:F, C{f} = s * A{f}; end
end

function v = cell_inner(A, B)
    v = 0; F = numel(A);
    for f = 1:F
        v = v + real(sum(sum(conj(A{f}) .* B{f})));
    end
end

function v = cell_norm_sq(A)
    v = cell_inner(A, A);
end
