function cohort_stats = cohort_consistency_from_sjj(Sjj_list, subject_ids)
%COHORT_CONSISTENCY_FROM_SJJ Compute cohort consistency from Sjj estimates.
%   cohort_stats = cohort_consistency_from_sjj(Sjj_list, subject_ids)
%   Sjj_list: cell array, each item is Nr x Nr x Nw (or cell of Nr x Nr)
%   subject_ids: optional cell array with subject labels
%
%   Outputs include:
%     - mean_abs_sjj: cohort mean of abs(Sjj)
%     - edge_cv_sjj: edge-wise CV across subjects for abs(Sjj)
%     - mean_coherence: cohort mean coherence
%     - edge_cv_coherence: edge-wise CV across subjects for coherence
%     - similarity_sjj: pairwise subject correlation on vec_{i<j,w}(abs(Sjj))
%     - similarity_coherence: pairwise subject correlation on vec_{i<j,w}(coherence)

    if nargin < 2 || isempty(subject_ids)
        subject_ids = cellfun(@(k) sprintf('S%03d', k), num2cell(1:numel(Sjj_list)), 'UniformOutput', false);
    end

    S = numel(Sjj_list);
    if S == 0
        cohort_stats = struct();
        return;
    end

    first_tensor = to_tensor_(Sjj_list{1});
    [Nr, ~, Nw] = size(first_tensor);
    edge_mask = triu(true(Nr), 1);
    nEdges = nnz(edge_mask);

    abs_sjj_edges = zeros(S, nEdges, Nw);
    coh_edges = zeros(S, nEdges, Nw);
    mean_abs_sjj = zeros(Nr, Nr, Nw);
    mean_coherence = zeros(Nr, Nr, Nw);

    for s = 1:S
        T = to_tensor_(Sjj_list{s});
        if ~isequal(size(T), [Nr, Nr, Nw])
            error('All Sjj tensors must have the same size.');
        end

        absT = abs(T);
        C = coherence_from_sjj_(T);

        mean_abs_sjj = mean_abs_sjj + absT;
        mean_coherence = mean_coherence + C;

        for w = 1:Nw
            Aw = absT(:, :, w);
            Cw = C(:, :, w);
            abs_sjj_edges(s, :, w) = Aw(edge_mask);
            coh_edges(s, :, w) = Cw(edge_mask);
        end
    end

    mean_abs_sjj = mean_abs_sjj / S;
    mean_coherence = mean_coherence / S;

    edge_cv_sjj = squeeze(std(abs_sjj_edges, 0, 1) ./ (mean(abs_sjj_edges, 1) + eps));
    edge_cv_coherence = squeeze(std(coh_edges, 0, 1) ./ (mean(coh_edges, 1) + eps));

    vec_sjj = reshape(abs_sjj_edges, [S, nEdges * Nw]);
    vec_coh = reshape(coh_edges, [S, nEdges * Nw]);

    similarity_sjj = pairwise_corr_(vec_sjj);
    similarity_coherence = pairwise_corr_(vec_coh);

    cohort_stats = struct();
    cohort_stats.subject_ids = subject_ids(:);
    cohort_stats.n_subjects = S;
    cohort_stats.n_roi = Nr;
    cohort_stats.n_freq = Nw;
    cohort_stats.mean_abs_sjj = mean_abs_sjj;
    cohort_stats.edge_cv_sjj = edge_cv_sjj;
    cohort_stats.mean_coherence = mean_coherence;
    cohort_stats.edge_cv_coherence = edge_cv_coherence;
    cohort_stats.similarity_sjj = similarity_sjj;
    cohort_stats.similarity_coherence = similarity_coherence;

    cohort_stats.similarity_sjj_mean = mean(similarity_sjj(triu(true(S), 1)), 'omitnan');
    cohort_stats.similarity_sjj_std = std(similarity_sjj(triu(true(S), 1)), 0, 'omitnan');
    cohort_stats.similarity_coherence_mean = mean(similarity_coherence(triu(true(S), 1)), 'omitnan');
    cohort_stats.similarity_coherence_std = std(similarity_coherence(triu(true(S), 1)), 0, 'omitnan');
end

function T = to_tensor_(Sjj)
    if iscell(Sjj)
        T = cat(3, Sjj{:});
    else
        T = Sjj;
    end
end

function C = coherence_from_sjj_(T)
    [Nr, ~, Nw] = size(T);
    C = zeros(Nr, Nr, Nw);
    for w = 1:Nw
        Sw = T(:, :, w);
        p = real(diag(Sw));
        p = max(p, eps);
        denom = sqrt(p * p.');
        Cw = abs(Sw) ./ denom;
        Cw(1:Nr+1:end) = 1;
        C(:, :, w) = Cw;
    end
end

function R = pairwise_corr_(X)
    S = size(X, 1);
    R = eye(S);
    for i = 1:S
        xi = X(i, :);
        for j = i+1:S
            xj = X(j, :);
            v = isfinite(xi) & isfinite(xj);
            if nnz(v) < 3
                r = NaN;
            else
                c = corrcoef(xi(v), xj(v));
                r = c(1, 2);
            end
            R(i, j) = r;
            R(j, i) = r;
        end
    end
end
