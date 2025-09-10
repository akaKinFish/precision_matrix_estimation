function Sigma = ensure_sigma_collection(C)
% ENSURE_SIGMA_COLLECTION  将 {F×1} 或 {1×F} cell 的 n×n 矩阵堆叠为 n×n×F 数组
% 并进行一次共轭对称化，保证严格 Hermitian。
%
% 用法：
%   Sigma = ensure_sigma_collection(Sjj_hat);
%
% 要求：
%   C 为 cell 数组，每个单元是同尺寸的 n×n 数值矩阵（可复数）。

    assert(iscell(C) && isvector(C), ...
        'Input must be a cell vector {F×1} or {1×F} containing n×n matrices.');

    F = numel(C);
    n = size(C{1}, 1);
    % 类型沿用第一个元素
    Sigma = zeros(n, n, F, 'like', C{1});

    for f = 1:F
        A = C{f};
        assert(ismatrix(A) && isequal(size(A), [n, n]), ...
            'Cell %d has size %s, expected %dx%d.', f, mat2str(size(A)), n, n);
        Sigma(:, :, f) = A;
    end

    % 共轭对称化（避免轻微数值不对称）
    Sigma = 0.5 * (Sigma + permute(conj(Sigma), [2 1 3]));
end