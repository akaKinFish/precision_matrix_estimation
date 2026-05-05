function results = run_jspace_stoch_age_cohort(age_range, site_dirs, parameters_mat_path, out_dir, cfg_overrides)

    if nargin < 1 || isempty(age_range)
        age_range = [18 30];
    end
    if nargin < 2 || isempty(site_dirs)
        site_dirs = default_site_dirs_();
    end
    if nargin < 3 || isempty(parameters_mat_path)
        parameters_mat_path = "F:\hope_the_final\precision_matrix_estimation\data\parameters.mat";
    end
    if nargin < 4 || isempty(out_dir)
        out_dir = fullfile(pwd, 'result_v3', 'jspace_stoch_age_cohort');
    end
    if nargin < 5
        cfg_overrides = struct();
    end

    if ~isfolder(out_dir)
        mkdir(out_dir);
    end

    P = load(parameters_mat_path);
    Model = P.Compact_Model;
    L = Model.K;
    dwi_C = [];
    if isfield(P, 'Compact_Model') && isfield(P.Compact_Model, 'C')
        dwi_C = P.Compact_Model.C;
    end

    subjects = find_crossspec_subjects_by_age(site_dirs, age_range);
    subjects = filter_existing_subjects(subjects, out_dir);
    nSub = numel(subjects);

    % ---------- 预分配 results（parfor 必须这么做） ----------
    template = struct( ...
        'subject_id', '', ...
        'age', NaN, ...
        'site_dir', '', ...
        'site_name', '', ...
        'out_file', '', ...
        'status', '', ...
        'message', '');
    results = repmat(template, nSub, 1);

    % ---------- 预创建 site 输出文件夹 ----------
    site_names = unique(arrayfun(@(s) string(site_name_from_dir_(s.site_dir)), subjects));
    for i = 1:numel(site_names)
        out_site_dir = fullfile(out_dir, char(site_names(i)));
        if ~isfolder(out_site_dir)
            mkdir(out_site_dir);
        end
    end

    % ---------- outer parallel 时，inner parallel / GPU 关掉 ----------
    cfg_run = cfg_overrides;
    cfg_run.opt_use_parallel = false;
    cfg_run.use_gpu = false;

    % ---------- 把真正会生效的 stoch 参数补上 ----------
    if ~isfield(cfg_run, 'stoch') || isempty(cfg_run.stoch)
        cfg_run.stoch = struct();
    end
    cfg_run.stoch.n_per_band = 2;       
    cfg_run.obj_stoch_max_iter = 10;    % surrogate objective 阶段
    cfg_run.em_stoch_max_iter  = 30;    % full EM 阶段

    % ---------- 显式开 pool，别完全依赖自动创建 ----------
    pool = gcp('nocreate');
    if isempty(pool)
        parpool("Processes", 8);  % 先从 6 开始试
    end

    parfor k = 1:nSub
        s = subjects(k);

        one = template;
        one.subject_id = s.subject_id;
        one.age = s.age;
        one.site_dir = s.site_dir;
        one.site_name = site_name_from_dir_(s.site_dir);

        try
            loaded = load(s.mat_path, 'data_struct');
            ds = loaded.data_struct;
            Svv_cross = ds.CrossM;

            if isfield(ds, 'freqrange')
                freq = ds.freqrange(:)';
                freq = freq(1:min(numel(freq), size(Svv_cross,3)));
            else
                freq = 1:size(Svv_cross, 3);
            end

            [Omega_est, Sjj_est, outs_js] = run_jspace_real(Svv_cross, L, freq, dwi_C, cfg_run);

            out_site_dir = fullfile(out_dir, one.site_name);
            out_file = fullfile(out_site_dir, [s.subject_id, '_jspace_stoch.mat']);

            subject_id = s.subject_id; %#ok<NASGU>
            age        = s.age; %#ok<NASGU>
            site_dir   = s.site_dir; %#ok<NASGU>
            site_name  = one.site_name; %#ok<NASGU>

            parsave_subject(out_file, subject_id, age, site_dir, site_name, Omega_est, Sjj_est, outs_js);

            one.out_file = out_file;
            one.status = 'ok';
            one.message = '';

        catch ME
            one.out_file = '';
            one.status = 'failed';
            one.message = ME.message;
        end

        results(k) = one;
    end

    summary_file = fullfile(out_dir, 'cohort_run_summary.mat');
    save(summary_file, 'results', 'subjects', 'age_range', 'site_dirs', 'parameters_mat_path', 'cfg_overrides');
end

function site_name = site_name_from_dir_(site_dir)
    [~, site_name] = fileparts(char(site_dir));
end

function site_dirs = default_site_dirs_()
%DEFAULT_SITE_DIRS_ Default site folders provided by the user.
    site_dirs = {
        'D:\\result\\Cuba90', 'D:\\result\\Cuba2003', 'D:\\result\\Cuba2024', ...
        'D:\\result\\Germany', 'D:\\result\\Malaysia', 'D:\\result\\NewYork', ...
        'D:\\result\\Rusia', 'D:\\result\\Rusia1', 'D:\\result\\Switzerland', ...
        'G:\\BC-V_result_12_6_part1\\Barbados', 'G:\\BC-V_result_12_6_part1\\Bern', ...
        'G:\\BC-V_result_12_6_part1\\CHBMP', 'G:\\BC-V_result_12_6_part1\\Chengdu', ...
        'G:\\BC-V_result_12_6_part1\\Chongqing', 'G:\\BC-V_result_12_6_part1\\Colombia' ...
    };
end

function subjects = filter_existing_subjects(subjects, out_dir)
    % FILTER_EXISTING_SUBJECTS 检查输出目录，移除已存在结果的样本
    if isempty(subjects)
        return;
    end
    
    nSub = numel(subjects);
    is_done = false(nSub, 1);
    
    for i = 1:nSub
        s = subjects(i);
        % 这里的逻辑要与主循环中生成 out_file 的逻辑严格一致

        % 将路径按分隔符拆分，取最后一个非空的部分
        parts = strsplit(char(s.site_dir), filesep);
        parts = parts(~cellfun(@isempty, parts)); 
        site_name = parts{end};
        out_file = fullfile(out_dir, site_name, [s.subject_id, '_jspace_stoch.mat']);
        
        if isfile(out_file)
            is_done(i) = true;
        end
    end
    
    % 报告进度
    fprintf('检查完毕：共发现 %d 个样本，其中 %d 个已存在，跳过。剩余 %d 个待处理。\n', ...
        nSub, sum(is_done), nSub - sum(is_done));
    
    subjects = subjects(~is_done);
end