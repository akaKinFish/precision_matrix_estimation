function listFolderStructure(folder, indent)
    if nargin < 2
        indent = 0;  % 初始缩进
    end

    % 获取当前文件夹内容
    files = dir(folder);
    
    for i = 1:length(files)
        % 跳过 . 和 ..
        if ~isempty(regexp(files(i).name, '^\.{1,2}$'))
            continue;
        end
        
        % 缩进显示层级
        fprintf('%s%s\n', blanks(indent * 2), files(i).name);
        
        % 如果是文件夹，递归进入
        if files(i).isdir
            subfolder = fullfile(folder, files(i).name);
            listFolderStructure(subfolder, indent + 1);
        end
    end
end