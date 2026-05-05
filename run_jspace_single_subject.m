function run_jspace_single_subject(input_mat, parameters_mat, out_path, cfg_overrides)
    % RUN_JSPACE_SINGLE_SUBJECT 针对单个 scalp .mat 文件运行 J-SPACE
    %
    % 参数:
    %   input_mat      : 输入的 scalp 数据路径 (包含 data_struct.CrossM)
    %   parameters_mat : 包含 Compact_Model (L 阵和 C 阵) 的路径
    %   out_path       : 结果保存的完整路径 (例如 'C:\res\sub01_res.mat')
    %   cfg_overrides  : 可选的配置覆盖结构体

    if nargin < 4, cfg_overrides = struct(); end

    %% 1. 环境准备
    fprintf('--- 正在处理单样本 ---\n输入: %s\n', input_mat);
    
    % 创建输出目录
    out_dir = fileparts(out_path);
    if ~isempty(out_dir) && ~isfolder(out_dir), mkdir(out_dir); end

    %% 2. 加载模型参数 (L 和 DWI Prior)
    P = load(parameters_mat);
    L = P.Compact_Model.K;
    dwi_C = [];
    if isfield(P.Compact_Model, 'C')
        dwi_C = P.Compact_Model.C;
    end

    %% 3. 加载被试数据
    loaded = load(input_mat, 'data_struct');
    ds = loaded.data_struct;
    Svv_cross = ds.CrossM;
    
    % 获取频率信息
    if isfield(ds, 'freqrange')
        freq = ds.freqrange(:)';
    else
        freq = 1:size(Svv_cross, 3);
    end

    %% 4. 配置算法参数
    % 这里复用你原有的 run_jspace_real 逻辑，但确保它是单机单核运行
    cfg_run = cfg_overrides;
    cfg_run.use_gpu = false;         % 单次运行建议先用 CPU 调试
    cfg_run.verbose = true;          % 开启详情打印
    cfg_run.opt_use_parallel = true; 

    %% 5. 执行算法
    try
        fprintf('开始运行 J-SPACE 算法...\n');
        % 注意：这里调用你已经定义好的 run_jspace_real 函数
        [Omega_est, Sjj_est, outs_js] = run_jspace_real(Svv_cross, L, freq, dwi_C, cfg_run);
        
        %% 6. 保存结果
        % 获取被试 ID（如果文件名就是 ID）
        [~, sub_id] = fileparts(input_mat);
        
        save(out_path, 'sub_id', 'Omega_est', 'Sjj_est', 'outs_js', 'input_mat', '-v7.3');
        fprintf('成功！结果已保存至: %s\n', out_path);
        
    catch ME
        fprintf('运行失败: %s\n', ME.message);
        rethrow(ME);
    end
end