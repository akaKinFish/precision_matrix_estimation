function report = summarize_jspace_cohort_analysis(analysis)
%SUMMARIZE_JSPACE_COHORT_ANALYSIS Summarize cohort consistency analysis result.
%   report = summarize_jspace_cohort_analysis(analysis)
%   analysis can be either:
%     - output struct from analyze_jspace_stoch_results_folder
%     - cohort_stats struct directly

    if nargin < 1 || isempty(analysis)
        error('Please provide analysis or cohort_stats input.');
    end

    if isfield(analysis, 'cohort_stats')
        s = analysis.cohort_stats;
    else
        s = analysis;
    end

    if ~isfield(s, 'n_subjects') || s.n_subjects < 2
        report = struct();
        report.n_subjects = getfield_with_default_(s, 'n_subjects', 0); %#ok<GFLD>
        report.overall = 'Not enough subjects for pairwise consistency analysis.';
        fprintf('[Report] %s\n', report.overall);
        return;
    end

    cv_sjj = s.edge_cv_sjj(:);
    cv_coh = s.edge_cv_coherence(:);

    report = struct();
    report.n_subjects = s.n_subjects;
    report.similarity_sjj_mean = s.similarity_sjj_mean;
    report.similarity_sjj_std = s.similarity_sjj_std;
    report.similarity_coherence_mean = s.similarity_coherence_mean;
    report.similarity_coherence_std = s.similarity_coherence_std;

    report.edge_cv_sjj_median = median(cv_sjj, 'omitnan');
    report.edge_cv_coherence_median = median(cv_coh, 'omitnan');
    report.edge_cv_sjj_p75 = prctile(cv_sjj, 75);
    report.edge_cv_coherence_p75 = prctile(cv_coh, 75);

    report.level_sjj = stability_level_(report.similarity_sjj_mean, report.edge_cv_sjj_median);
    report.level_coherence = stability_level_(report.similarity_coherence_mean, report.edge_cv_coherence_median);

    report.overall = sprintf('Sjj=%s | coherence=%s', report.level_sjj, report.level_coherence);

    fprintf('[Report] Subjects: %d\n', report.n_subjects);
    fprintf('[Report] Similarity(abs(Sjj)): mean=%.4f std=%.4f\n', report.similarity_sjj_mean, report.similarity_sjj_std);
    fprintf('[Report] Similarity(coherence): mean=%.4f std=%.4f\n', report.similarity_coherence_mean, report.similarity_coherence_std);
    fprintf('[Report] Median edge CV(abs(Sjj))=%.4f | P75=%.4f\n', report.edge_cv_sjj_median, report.edge_cv_sjj_p75);
    fprintf('[Report] Median edge CV(coherence)=%.4f | P75=%.4f\n', report.edge_cv_coherence_median, report.edge_cv_coherence_p75);
    fprintf('[Report] Overall: %s\n', report.overall);
end

function level = stability_level_(sim_mean, cv_median)
    if sim_mean >= 0.70 && cv_median <= 0.50
        level = 'stable';
    elseif sim_mean >= 0.50 && cv_median <= 0.80
        level = 'moderate';
    else
        level = 'unstable';
    end
end

function v = getfield_with_default_(s, f, default_v)
    if isfield(s, f)
        v = s.(f);
    else
        v = default_v;
    end
end
