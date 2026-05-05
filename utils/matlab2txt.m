% 设置项目路径和输出文件路径
projectFolder = 'F:\hope_the_final\precision_matrix_estimation-develop\precision_matrix_estimation-develop';  
outputFile   = fullfile('C:\Users\kin\Desktop\code\all_code.txt');

% 递归查找所有 .m 文件
files = dir(fullfile(projectFolder, '**', '*.m'));

% === 关键过滤：排除 ** 文件夹下的文件 ===
is_bad = false(numel(files),1);
for k = 1:numel(files)
    % 统一路径分隔符，避免 Windows / Linux 差异
    folder_k = strrep(files(k).folder, '/', filesep);
    
    % 判断是否包含 \**\ 或 以 \** 结尾
    if contains(folder_k, [filesep 'JSPACE_Spline_v1' filesep]) || ...
       endsWith(folder_k, [filesep 'JSPACE_Spline_v1'])
        is_bad(k) = true;
    end
end
files = files(~is_bad);

fprintf('保留 %d 个 .m 文件（已排除 ** 文件夹）\n', numel(files));

% 打开输出文件
fid_out = fopen(outputFile, 'w', 'n', 'UTF-8');

for k = 1:numel(files)
    filePath = fullfile(files(k).folder, files(k).name);

    % 文件分隔标记
    fprintf(fid_out, '%%%% ==================================================\n');
    fprintf(fid_out, '%%%% File: %s\n', filePath);
    fprintf(fid_out, '%%%% ==================================================\n');

    % 读入并写出源码
    fid_in = fopen(filePath, 'r', 'n', 'UTF-8');
    tline = fgetl(fid_in);
    while ischar(tline)
        fprintf(fid_out, '%s\n', tline);
        tline = fgetl(fid_in);
    end
    fclose(fid_in);

    fprintf(fid_out, '\n\n');
end

fclose(fid_out);

disp('✅ 所有代码已导出到 all_code.txt（** 文件夹已排除）');