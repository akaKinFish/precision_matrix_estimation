function total_changes = count_sparsity_changes(prec_matrices)
% COUNT_SPARSITY_CHANGES - Count how many times edge patterns change across frequencies
%
% INPUT:
%   prec_matrices - Cell array of precision matrices {F x 1}
%
% OUTPUT:
%   total_changes - Total number of edge pattern changes across frequencies

F = length(prec_matrices);
if F < 2
    total_changes = 0;
    return;
end

n = size(prec_matrices{1}, 1);
total_changes = 0;

% Convert to binary sparsity patterns
sparsity_patterns = cell(F, 1);
for f = 1:F
    P = abs(prec_matrices{f});
    % Create binary matrix for upper triangular part (excluding diagonal)
    sparsity_patterns{f} = triu(P > 1e-10, 1);
end

% Count changes between consecutive frequencies
for f = 1:(F-1)
    pattern_f = sparsity_patterns{f};
    pattern_f1 = sparsity_patterns{f+1};
    
    % Count positions where patterns differ
    differences = xor(pattern_f, pattern_f1);
    changes_this_step = sum(differences(:));
    total_changes = total_changes + changes_this_step;
end

end