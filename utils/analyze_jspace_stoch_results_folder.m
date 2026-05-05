function analysis = analyze_jspace_stoch_results_folder(results_dir, output_mat_path)
%ANALYZE_JSPACE_STOCH_RESULTS_FOLDER Analyze cohort consistency from result folder.
%   analysis = analyze_jspace_stoch_results_folder(results_dir, output_mat_path)
%   reads all *_jspace_stoch.mat files in results_dir, extracts Sjj_est,
%   computes cohort consistency metrics using cohort_consistency_from_sjj,
%   and optionally saves analysis to output_mat_path.
%
%   Input:
%     results_dir      - folder containing subject result .mat files
%     output_mat_path  - optional output .mat path; default is
%                        fullfile(results_dir, 'cohort_consistency_analysis.mat')
%
%   Output fields:
%     - cohort_stats
%     - loaded_files
%     - n_loaded
%     - n_failed
%     - failed_files

    if nargin < 1 || isempty(results_dir)
        error('Please provide results_dir.');
    end
    if nargin < 2 || isempty(output_mat_path)
        output_mat_path = fullfile(results_dir, 'cohort_consistency_analysis.mat');
    end

    if ~isfolder(results_dir)
        error('results_dir does not exist: %s', results_dir);
    end

    files = dir(fullfile(results_dir, '*_jspace_stoch.mat'));
    if isempty(files)
        warning('No *_jspace_stoch.mat files found in %s', results_dir);
    end

    Sjj_list = {};
    subject_ids = {};
    loaded_files = {};
    failed_files = {};

    for k = 1:numel(files)
        file_path = fullfile(files(k).folder, files(k).name);
        try
            d = load(file_path);
            if ~isfield(d, 'Sjj_est')
                failed_files{end + 1} = file_path; %#ok<AGROW>
                continue;
            end

            Sjj_list{end + 1} = d.Sjj_est; %#ok<AGROW>
            loaded_files{end + 1} = file_path; %#ok<AGROW>

            if isfield(d, 'subject_id') && ~isempty(d.subject_id)
                if isstring(d.subject_id)
                    subject_ids{end + 1} = char(d.subject_id); %#ok<AGROW>
                else
                    subject_ids{end + 1} = d.subject_id; %#ok<AGROW>
                end
            else
                [~, name] = fileparts(files(k).name);
                subject_ids{end + 1} = erase(name, '_jspace_stoch'); %#ok<AGROW>
            end
        catch
            failed_files{end + 1} = file_path; %#ok<AGROW>
        end
    end

    cohort_stats = cohort_consistency_from_sjj(Sjj_list, subject_ids);

    analysis = struct();
    analysis.cohort_stats = cohort_stats;
    analysis.loaded_files = loaded_files;
    analysis.n_loaded = numel(loaded_files);
    analysis.failed_files = failed_files;
    analysis.n_failed = numel(failed_files);

    save(output_mat_path, 'analysis');

    fprintf('[Analyze] Loaded %d files, failed %d files.\n', analysis.n_loaded, analysis.n_failed);
    if isfield(cohort_stats, 'n_subjects') && cohort_stats.n_subjects >= 2
        fprintf('[Analyze] Similarity | abs(Sjj): mean=%.4f std=%.4f\n', ...
            cohort_stats.similarity_sjj_mean, cohort_stats.similarity_sjj_std);
        fprintf('[Analyze] Similarity | coherence: mean=%.4f std=%.4f\n', ...
            cohort_stats.similarity_coherence_mean, cohort_stats.similarity_coherence_std);
    end
    fprintf('[Analyze] Analysis saved to %s\n', output_mat_path);
end
