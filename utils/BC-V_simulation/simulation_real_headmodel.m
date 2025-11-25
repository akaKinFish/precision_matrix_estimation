function [] = simulation_real_headmodel(sens_system, output_source, paths, varargin)
%% SIMULATION_REAL_HEADMODEL - 使用真实头模型进行MEEG源连接性仿真
%
% 用法:
%   simulation_real_headmodel(sens_system, output_source, paths)
%   simulation_real_headmodel(sens_system, output_source, paths, 'Name', Value, ...)
%
% 输入:
%   sens_system    - 传感器系统类型，固定为 'real_head'
%   output_source  - 输出目录路径
%   paths          - 结构体，包含以下字段：
%       .headmodel_mat - 头模型文件路径 (必需)
%                        应包含: HeadModel.Gain, GridOrient 或 VertNormals
%       .surf_mat      - 皮层表面文件路径 (必需)
%                        应包含: Sc.Vertices, Sc.Faces, VertNormals, (可选)Atlas
%       .channel_mat   - 通道文件路径 (可选，为将来扩展保留)
%
% 可选参数 (Name-Value pairs):
%   'Nseed'           - 种子点数量 (默认: 22)
%   'Nsamp'           - 样本数量 (默认: 6000)
%   'd0'              - patch 测地线半径 (默认: 1e-2，单位：米)
%   'seed_method'     - 种子点选择方法 (默认: 'maxdist')
%                       'random'  - 随机采样
%                       'maxdist' - 最大距离采样（推荐）
%                       'atlas'   - 基于脑区采样
%   'random_seed'     - 随机数种子 (默认: 2024)
%   'Nsim'            - 仿真次数 (默认: 3)
%   'db_source'       - 源噪声水平 (默认: 0.1)
%   'db_sens'         - 传感器噪声水平 (默认: 0.1)
%
% 示例:
%   paths.headmodel_mat = 'D:\result\Cuba2003\2pre\BC-V_Structure\F4XW6TBLHF3Q\leadfield\headmodel.mat';
%   paths.surf_mat = 'D:\result\Cuba2003\2pre\BC-V_Structure\F4XW6TBLHF3Q\surf\surf.mat';
%   paths.channel_mat = 'D:\result\Cuba2003\2pre\BC-V_Structure\F4XW6TBLHF3Q\channel\channel.mat';
%   
%   simulation_real_headmodel('real_head', './output', paths, ...
%       'Nseed', 22, 'Nsamp', 6000, 'd0', 1e-2, 'seed_method', 'maxdist');
%
% 作者: Augment Agent
% 日期: 2024
% 基于: Sim4_h_head_model_comparison/Simulation.m 和 Initialize_simulation.m

%% 参数解析
p = inputParser;
addRequired(p, 'sens_system', @(x) strcmp(x, 'real_head'));
addRequired(p, 'output_source', @ischar);
addRequired(p, 'paths', @isstruct);
addParameter(p, 'Nseed', 22, @isnumeric);
addParameter(p, 'Nsamp', 6000, @isnumeric);
addParameter(p, 'd0', 1e-2, @isnumeric);
addParameter(p, 'seed_method', 'maxdist', @ischar);
addParameter(p, 'random_seed', 2024, @isnumeric);
addParameter(p, 'Nsim', 3, @isnumeric);
addParameter(p, 'db_source', 0.1, @isnumeric);
addParameter(p, 'db_sens', 0.1, @isnumeric);
parse(p, sens_system, output_source, paths, varargin{:});
params = p.Results;

%% 设置随机种子
rng(params.random_seed);
fprintf('随机种子已设置为: %d\n', params.random_seed);

%% 验证路径
if ~isfield(paths, 'headmodel_mat') || ~exist(paths.headmodel_mat, 'file')
    error('必须提供有效的 headmodel_mat 路径');
end
if ~isfield(paths, 'surf_mat') || ~exist(paths.surf_mat, 'file')
    error('必须提供有效的 surf_mat 路径');
end

%% 加载头模型和导联场
fprintf('正在加载头模型: %s\n', paths.headmodel_mat);
headmodel_data = load(paths.headmodel_mat);

% 提取导联场矩阵
if isfield(headmodel_data, 'HeadModel') && isfield(headmodel_data.HeadModel, 'Gain')
    K = headmodel_data.HeadModel.Gain;
elseif isfield(headmodel_data, 'Gain')
    K = headmodel_data.Gain;
else
    error('头模型文件中未找到 Gain 矩阵');
end

% 提取方向信息（如果需要）
if isfield(headmodel_data, 'GridOrient')
    GridOrient = headmodel_data.GridOrient;
elseif isfield(headmodel_data, 'HeadModel') && isfield(headmodel_data.HeadModel, 'GridOrient')
    GridOrient = headmodel_data.HeadModel.GridOrient;
elseif isfield(headmodel_data, 'VertNormals')
    GridOrient = headmodel_data.VertNormals;
else
    warning('未找到 GridOrient 或 VertNormals');
    GridOrient = [];
end

fprintf('导联场矩阵大小: [%d × %d]\n', size(K, 1), size(K, 2));

%% 加载皮层表面
fprintf('正在加载皮层表面: %s\n', paths.surf_mat);
surf_data = load(paths.surf_mat);

% 提取顶点和面片
if isfield(surf_data, 'Sc')
    vertices = surf_data.Sc.Vertices;
    faces = surf_data.Sc.Faces;
elseif isfield(surf_data, 'Vertices') && isfield(surf_data, 'Faces')
    vertices = surf_data.Vertices;
    faces = surf_data.Faces;
else
    error('皮层表面文件中未找到 Vertices 和 Faces');
end

% 提取顶点法向量
if isfield(surf_data, 'VertNormals')
    VertNormals = surf_data.VertNormals;
elseif isfield(surf_data, 'Sc') && isfield(surf_data.Sc, 'VertNormals')
    VertNormals = surf_data.Sc.VertNormals;
else
    warning('未找到顶点法向量');
    VertNormals = [];
end

% 提取 Atlas（可选）
if isfield(surf_data, 'Atlas')
    Atlas = surf_data.Atlas;
elseif isfield(surf_data, 'Sc') && isfield(surf_data.Sc, 'Atlas')
    Atlas = surf_data.Sc.Atlas;
else
    Atlas = [];
end

fprintf('皮层表面: %d 个顶点, %d 个面片\n', size(vertices, 1), size(faces, 1));

%% 验证坐标单位并给出建议
vertex_range = max(vertices) - min(vertices);
max_range = max(vertex_range);
fprintf('顶点坐标范围: [%.4f, %.4f, %.4f]\n', vertex_range);
if max_range < 1
    fprintf('✓ 检测到坐标单位为米(m)，d0=%.4f 表示 %.1f 厘米\n', params.d0, params.d0*100);
elseif max_range > 10
    fprintf('⚠ 警告：坐标单位可能是毫米(mm)，建议 d0=%.1f\n', params.d0*1000);
else
    fprintf('⚠ 警告：无法确定坐标单位，请手动验证 d0 参数\n');
end

%% 加载通道信息（可选）
if isfield(paths, 'channel_mat') && exist(paths.channel_mat, 'file')
    fprintf('正在加载通道信息: %s\n', paths.channel_mat);
    channel_data = load(paths.channel_mat);
    if isfield(channel_data, 'Channel')
        elec_pos = zeros(length(channel_data.Channel), 3);
        for i = 1:length(channel_data.Channel)
            if isfield(channel_data.Channel(i), 'Loc') && ~isempty(channel_data.Channel(i).Loc)
                elec_pos(i, :) = channel_data.Channel(i).Loc(:, 1)';
            end
        end
    else
        warning('通道文件格式不支持，使用默认电极位置');
        elec_pos = zeros(size(K, 1), 3);
    end
else
    fprintf('未提供通道文件，使用默认电极位置\n');
    elec_pos = zeros(size(K, 1), 3);
end

%% 设置仿真参数
LeadFields = {K};  % 封装为cell数组以保持与原代码一致
Nv = size(vertices, 1);
Nsubj = 1;  % 单个受试者
Nsim = params.Nsim;
Nsamp = params.Nsamp;
Nseed = params.Nseed;
d0 = params.d0;
db_source = params.db_source;
db_sens = params.db_sens;

fprintf('\n仿真参数设置:\n');
fprintf('  传感器数量: %d\n', size(K, 1));
fprintf('  源点数量: %d\n', Nv);
fprintf('  种子点数量: %d\n', Nseed);
fprintf('  样本数量: %d\n', Nsamp);
fprintf('  测地线半径: %.4f\n', d0);
fprintf('  源噪声水平: %.2f\n', db_source);
fprintf('  传感器噪声水平: %.2f\n', db_sens);
fprintf('  仿真次数: %d\n', Nsim);

%% 生成种子点（data_tips）
fprintf('\n正在使用 %s 方法生成 %d 个种子点...\n', params.seed_method, Nseed);
data_tips = generate_seed_points(vertices, faces, Nseed, params.seed_method, Atlas);

fprintf('种子点生成完成！\n');
for i = 1:min(5, Nseed)
    fprintf('  种子点 %d: [%.4f, %.4f, %.4f]\n', i, ...
        data_tips(i).Position(1), data_tips(i).Position(2), data_tips(i).Position(3));
end
if Nseed > 5
    fprintf('  ... (共 %d 个种子点)\n', Nseed);
end

%% 将种子点坐标转换为顶点索引
Seeders = zeros(1, Nseed);
for cont = 1:Nseed
    coord_tmp = data_tips(cont).Position;
    vx = coord_tmp(1);
    vy = coord_tmp(2);
    vz = coord_tmp(3);
    Seeders(cont) = pickpoint(vx, vy, vz, vertices, 1E-5);
end
Seeders = sort(Seeders);

fprintf('种子点索引: ');
fprintf('%d ', Seeders(1:min(10, Nseed)));
if Nseed > 10
    fprintf('...');
end
fprintf('\n');

%% 设置 gen_hggm 选项
options.config = 2;
options.var = 2;
options.connections = [1 2; 2 3];
options.extensions = [ceil(Nseed/3); ceil(Nseed/3); Nseed - 2*ceil(Nseed/3)];

%% 生成 patch 索引
fprintf('\n正在生成皮层 patches...\n');
index_seed = cell(1, Nseed);
index_full = [];
for point = 1:Nseed
    Source = Seeders(point);
    [index, findex] = surfpatch(Source, vertices, faces, d0);
    index_seed{point} = index;
    index_full = [index_full; index];
    if mod(point, 5) == 0 || point == Nseed
        fprintf('  进度: %d/%d patches\n', point, Nseed);
    end
end
Nnoise = length(index_full);
fprintf('Patches 生成完成！总噪声源数量: %d\n', Nnoise);

%% 运行仿真
fprintf('\n开始仿真...\n');
J_sim = cell(1, Nsim);
Svv0_sim = cell(1, Nsim);
Svv_sim = cell(1, Nsim);
Theta_sim = cell(1, Nsim);
Seeders_sim = zeros(Nseed, Nsim);

process_waitbar1 = waitbar(0, 'Please wait...');

for sim = 1:Nsim
    waitbar((sim)/(Nsim), process_waitbar1, ...
        strcat('simulation # ', num2str(sim), '  to sens-system: ', sens_system));
    disp(['simulation # ', num2str(sim)])
    
    % 设置相关结构
    [S, Data, X] = gen_hggm2(Nsamp, sum(options.extensions), options);
    
    %% 模拟传感器时间序列数据
    V0 = cell(1, Nsubj);
    for cont_seed = 1:Nseed
        Seeders_sim(cont_seed, sim) = index_seed{cont_seed}(randi(length(index_seed{cont_seed})));
    end
    Seeders_sim(:, sim) = Seeders_sim(randperm(Nseed), sim);
    
    for cont4 = 1:Nsubj
        K_tmp = LeadFields{cont4};
        V0{cont4} = K_tmp(:, Seeders_sim(:, sim)) * Data;
    end
    
    %% 计算源方差
    J_tmp = sum(abs(Data).^2, 2);
    J = zeros(Nv, 1);
    J(Seeders_sim(:, sim), :) = J_tmp / Nsamp;
    
    %% 模拟生物噪声
    rs0 = randn(Nnoise, Nsamp, Nsubj) + 1i*randn(Nnoise, Nsamp, Nsubj);
    
    %% 投影生物噪声到传感器
    noisesources = cell(1, Nsubj);
    for cont4 = 1:Nsubj
        K_tmp = LeadFields{cont4};
        Ne = size(K_tmp, 1);
        rs_tmp = K_tmp(:, index_full) * squeeze(rs0(:, :, cont4));
        noisesources{cont4} = db_source * sum(abs(V0{cont4}(:)).^2)^(1/2) * rs_tmp / sum(abs(rs_tmp(:)).^2)^(1/2);
    end
    
    %% 模拟传感器噪声
    noisesensors = cell(1, Nsubj);
    for cont4 = 1:Nsubj
        K_tmp = LeadFields{cont4};
        Ne = size(K_tmp, 1);
        rs_tmp = randn(Ne, Nsamp) + 1i*randn(Ne, Nsamp);
        noisesensors{cont4} = db_sens * sum(abs(V0{cont4}(:)).^2)^(1/2) * rs_tmp / sum(abs(rs_tmp(:)).^2)^(1/2);
    end
    
    %% 模拟数据: K*J + 传感器噪声 + 生物噪声
    V = cell(1, Nsubj);
    for cont4 = 1:Nsubj
        V{cont4} = V0{cont4} + noisesources{cont4} + noisesensors{cont4};
    end
    
    %% 计算协方差
    Svv = cell(1, Nsubj);
    Svv0 = cell(1, Nsubj);
    for cont4 = 1:Nsubj
        Svv0{cont4} = squeeze(V0{cont4}) * squeeze(V0{cont4})' / Nsamp;
        Svv{cont4} = squeeze(V{cont4}) * squeeze(V{cont4})' / Nsamp;
    end
    
    %% 保存结果
    J_sim{:, sim} = J;
    Svv0_sim{:, sim} = Svv0;
    Svv_sim{:, sim} = Svv;
    Theta_sim{:, sim} = X;
end

delete(process_waitbar1);

%% 保存仿真结果
save_path = strcat(output_source, filesep, 'Pseudorand_Net');
save(save_path, 'J_sim', 'Seeders_sim', 'index_full', 'Svv0_sim', 'Svv_sim', ...
    'Theta_sim', 'vertices', 'faces', 'elec_pos', 'LeadFields', 'Nsubj', 'Nsamp', ...
    'sens_system', 'options');

fprintf('\n仿真完成！结果已保存到: %s.mat\n', save_path);

end

%% ========== 辅助函数 ==========

function data_tips = generate_seed_points(vertices, faces, Nseed, method, Atlas)
%% 生成种子点的三种方法
% 输入:
%   vertices - 顶点坐标 [Nv × 3]
%   faces    - 面片 [Nf × 3]
%   Nseed    - 种子点数量
%   method   - 选择方法: 'random', 'maxdist', 'atlas'
%   Atlas    - 脑区信息（可选）
% 输出:
%   data_tips - 结构体数组，每个元素包含 Position 字段

switch lower(method)
    case 'random'
        % 方法1: 随机均匀采样
        Nv = size(vertices, 1);
        seed_indices = randperm(Nv, Nseed);

    case 'maxdist'
        % 方法2: 最大距离采样（推荐）
        % 使用贪心算法选择彼此距离最大的点
        Nv = size(vertices, 1);
        seed_indices = zeros(1, Nseed);

        % 第一个点：选择最中心的点
        center = mean(vertices, 1);
        [~, seed_indices(1)] = min(sum((vertices - center).^2, 2));

        fprintf('  最大距离采样进度: ');
        for i = 2:Nseed
            % 计算所有点到已选点的最小距离
            min_dist = inf(Nv, 1);
            for j = 1:(i-1)
                dist = sum((vertices - vertices(seed_indices(j), :)).^2, 2);
                min_dist = min(min_dist, dist);
            end

            % 选择距离已选点最远的点
            [~, seed_indices(i)] = max(min_dist);

            if mod(i, 5) == 0 || i == Nseed
                fprintf('%d/%d ', i, Nseed);
            end
        end
        fprintf('\n');

    case 'atlas'
        % 方法3: 基于脑区采样
        if isempty(Atlas)
            warning('Atlas信息为空，改用maxdist方法');
            data_tips = generate_seed_points(vertices, faces, Nseed, 'maxdist', []);
            return;
        end

        % 获取所有脑区
        if isfield(Atlas, 'Scouts')
            scouts = Atlas.Scouts;
            Nregions = length(scouts);

            % 如果脑区数量少于种子点数量，需要从某些脑区选择多个点
            if Nregions < Nseed
                % 每个脑区选择的点数
                points_per_region = floor(Nseed / Nregions);
                extra_points = Nseed - points_per_region * Nregions;

                seed_indices = [];
                for i = 1:Nregions
                    region_vertices = scouts(i).Vertices;
                    n_select = points_per_region;
                    if i <= extra_points
                        n_select = n_select + 1;
                    end

                    % 从该脑区随机选择点
                    if length(region_vertices) >= n_select
                        selected = randperm(length(region_vertices), n_select);
                        seed_indices = [seed_indices, region_vertices(selected)];
                    else
                        seed_indices = [seed_indices, region_vertices];
                    end
                end
            else
                % 从前Nseed个脑区各选一个代表点
                seed_indices = zeros(1, Nseed);
                for i = 1:Nseed
                    region_vertices = scouts(i).Vertices;
                    % 选择脑区中心点
                    region_center = mean(vertices(region_vertices, :), 1);
                    [~, idx] = min(sum((vertices(region_vertices, :) - region_center).^2, 2));
                    seed_indices(i) = region_vertices(idx);
                end
            end
        else
            warning('Atlas格式不支持，改用maxdist方法');
            data_tips = generate_seed_points(vertices, faces, Nseed, 'maxdist', []);
            return;
        end

    otherwise
        error('未知的种子点选择方法: %s', method);
end

% 构建 data_tips 结构体
data_tips = struct('Position', cell(1, Nseed));
for i = 1:Nseed
    data_tips(i).Position = vertices(seed_indices(i), :);
end

end

