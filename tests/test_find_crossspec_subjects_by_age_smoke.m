function test_find_crossspec_subjects_by_age_smoke()
%TEST_FIND_CROSSSPEC_SUBJECTS_BY_AGE_SMOKE Fast smoke test for age filtering.

    root = tempname;
    mkdir(root);
    c = onCleanup(@() cleanup_tmp_(root));

    site1 = fullfile(root, 'Cuba2003');
    site2 = fullfile(root, 'CHBMP');
    make_subject_(site1, 'SUBJ_A', 21, true);
    make_subject_(site1, 'SUBJ_B', 17, true);
    make_subject_(site2, 'SUBJ_C', 30, true);
    make_subject_(site2, 'SUBJ_D', 45, true);
    make_subject_(site2, 'SUBJ_E', 25, false); % missing CrossM -> excluded

    subjects = find_crossspec_subjects_by_age({site1, site2}, [20 35]);

    ids = sort({subjects.subject_id});
    assert(isequal(ids, {'SUBJ_A', 'SUBJ_C'}), 'Age selection mismatch');
    assert(all([subjects.age] >= 20 & [subjects.age] <= 35), 'Age range check failed');

    clear c;
    cleanup_tmp_(root);
    fprintf('test_find_crossspec_subjects_by_age_smoke: PASS\n');
end

function make_subject_(site_dir, subject_id, age, with_crossm)
    subject_dir = fullfile(site_dir, '2pre', 'crossSpec', subject_id);
    if ~isfolder(subject_dir)
        mkdir(subject_dir);
    end

    data_struct = struct();
    data_struct.age = age;
    data_struct.freqrange = 1:5;
    if with_crossm
        data_struct.CrossM = repmat(eye(3), 1, 1, 5);
    end

    save(fullfile(subject_dir, [subject_id, '.mat']), 'data_struct');
end

function cleanup_tmp_(root)
    if isfolder(root)
        rmdir(root, 's');
    end
end
