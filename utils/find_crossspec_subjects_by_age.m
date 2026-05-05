function subjects = find_crossspec_subjects_by_age(site_dirs, age_range)
%FIND_CROSSSPEC_SUBJECTS_BY_AGE Collect subject MAT files and filter by age.
%   subjects = find_crossspec_subjects_by_age(site_dirs, age_range)
%   site_dirs: cellstr/string array of site folders, each containing
%              2pre/crossSpec/<subject_id>/<subject_id>.mat.
%   age_range: [min_age max_age], inclusive.
%
%   Output subjects is a struct array with fields:
%     - site_dir
%     - subject_id
%     - mat_path
%     - age

    if nargin < 2 || isempty(age_range)
        age_range = [-inf, inf];
    end
    if isstring(site_dirs)
        site_dirs = cellstr(site_dirs);
    end

    subjects = struct('site_dir', {}, 'subject_id', {}, 'mat_path', {}, 'age', {});

    for iSite = 1:numel(site_dirs)
        site_dir = site_dirs{iSite};
        crossspec_dir = fullfile(site_dir, '2pre', 'crossSpec');
        if ~isfolder(crossspec_dir)
            continue;
        end

        folders = dir(crossspec_dir);
        folders = folders([folders.isdir] & ~startsWith({folders.name}, '.'));

        for iSub = 1:numel(folders)
            subject_id = folders(iSub).name;
            mat_path = fullfile(crossspec_dir, subject_id, [subject_id, '.mat']);
            if ~isfile(mat_path)
                continue;
            end

            loaded = load(mat_path, 'data_struct');
            if ~isfield(loaded, 'data_struct') || ~isstruct(loaded.data_struct)
                continue;
            end
            if ~isfield(loaded.data_struct, 'age') || ~isfield(loaded.data_struct, 'CrossM')
                continue;
            end
            if isnumeric(loaded.data_struct.age)
                age = loaded.data_struct.age;
            else
                age = double(str2double(loaded.data_struct.age));
            end

            
            if isempty(age) || ~isfinite(age) || age < age_range(1) || age > age_range(2)
                continue;
            end

            subjects(end + 1) = struct( ...
                'site_dir', site_dir, ...
                'subject_id', subject_id, ...
                'mat_path', mat_path, ...
                'age', age); %#ok<AGROW>
        end
    end
end
