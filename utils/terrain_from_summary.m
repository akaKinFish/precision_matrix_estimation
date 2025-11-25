function out = terrain_from_summary(csv_path, varargin)
% TERRAIN_FROM_SUMMARY  从 summary.csv 画 3D 地形图，找最优超参
% 用法：
%   terrain_from_summary('results/.../summary.csv');                % 综合名次
%   terrain_from_summary('summary.csv','Score','AUPR_mean');        % 只看 AUPR
%   terrain_from_summary('summary.csv','SaveDir','figs');           % 保存图片
%
% 参数（Name-Value）:
%   'Score'   : 'balanced'(默认) 或表中某列名（如 'LE_mean','AUPR_mean','F1_mean','Jaccard_mean',...）
%   'SaveDir' : 字符串；若提供则把每个切片的图保存为 PNG
%
% 说明：
%   - Z 轴(越低越好)：
%       * Score='balanced'：对 LE_mean 升序名次 + 对 AUPR/F1/Jaccard 降序名次 的总和
%       * Score=某“最大化”指标（AUROC/AUPR/F1/Jaccard/Prec/Rec）：绘图时用 1 - 指标
%       * Score=某“最小化”指标（如 LE_mean）：绘图时直接用该值
%
% 返回：
%   out 结构体，含 parsed 表、每个切片的最优点、全局最优点等

%% 读取
T = readtable(csv_path, 'TextType','string');

% 确保有 key 列
assert(any(strcmpi(T.Properties.VariableNames,'key')), ...
    'CSV 中需要名为 "key" 的列（例如 l1_0_12__l2f_2_1__l3r_0_5__sf_4）。');

% 尝试把常见指标列转为实数（有时 CSV 会带极小的虚部字符串）
numcols = ["LE_mean","AUPR_mean","F1_mean","Jaccard_mean","AUROC_mean", ...
    "AIRM_mean","JBLD_mean","BW_mean","KLS_mean", ...
    "Prec_mean","Rec_mean"];
for c = numcols
    if any(strcmp(T.Properties.VariableNames, c))
        T.(c) = to_real_numeric(T.(c));
    end
end

%% 解析 key -> 超参
P = arrayfun(@parse_key_to_params, T.key, 'uni', 1);
T.l1   = [P.l1]';
T.l2f  = [P.l2f]';
T.l3r  = [P.l3r]';
T.sf   = [P.sf]';

%% 解析参数
p = inputParser;
addParameter(p,'Score','balanced',@(s)ischar(s)||isstring(s));
addParameter(p,'SaveDir','',@(s)ischar(s)||isstring(s));
parse(p,varargin{:});
score_name = string(p.Results.Score);
save_dir   = string(p.Results.SaveDir);

if ~isempty(save_dir) && ~exist(save_dir,'dir')
    mkdir(save_dir);
end

% 决定绘图用的 Z 值（越低越好）
is_max_metric = @(name) any(strcmpi(name, ...
    {'AUPR_mean','AUROC_mean','F1_mean','Jaccard_mean','Prec_mean','Rec_mean'}));

switch lower(score_name)
    case "balanced"
        need = ["LE_mean","AUPR_mean","F1_mean","Jaccard_mean"];
        miss = need(~ismember(need, T.Properties.VariableNames));
        assert(isempty(miss), 'balanced 需要列: %s', strjoin(miss, ', '));
        % 名次（LE 越小越好；其他越大越好）
        rLE  = tiedrank(T.LE_mean, 1);           % 升序
        rAUPR= tiedrank(-T.AUPR_mean, 1);        % 降序 -> 对负号升序
        rF1  = tiedrank(-T.F1_mean, 1);
        rJ   = tiedrank(-T.Jaccard_mean, 1);
        Zval = rLE + rAUPR + rF1 + rJ;           % 越小越好
        zlab = "balanced rank (↓ is better)";
    otherwise
        % 单指标
        assert(any(strcmp(T.Properties.VariableNames, score_name)), ...
            '找不到列 "%s"。', score_name);
        v = T.(score_name);
        if is_max_metric(score_name)
            Zval = 1 - v;          % 最大化 -> 用 1-值，这样“越低越好”
            zlab = score_name + " (plotted as 1 - metric; ↓ is better)";
        else
            Zval = v;              % 最小化 -> 原值
            zlab = score_name + " (↓ is better)";
        end
end

T.Z = Zval;

%% 切片并绘图：每个 (λ3, σf) 一张
l3s = unique(T.l3r);
sfs = unique(T.sf);

best_per_slice = [];
figs = gobjects(numel(l3s)*numel(sfs),1);
fidx = 0;
zmin = nanmin(T.Z);
zmax = nanmax(T.Z);
for i = 1:numel(l3s)
    for j = 1:numel(sfs)
        sl = T(T.l3r==l3s(i) & T.sf==sfs(j), :);
        if isempty(sl), continue; end

        % 构造规则网格：X=λ2f, Y=λ1
        xVals = sort(unique(sl.l2f));
        yVals = sort(unique(sl.l1));
        [X,Y] = meshgrid(xVals, yVals);
        Z = nan(size(X));

        for ii = 1:numel(yVals)
            for jj = 1:numel(xVals)
                q = sl(sl.l1==yVals(ii) & sl.l2f==xVals(jj),:);
                if ~isempty(q)
                    Z(ii,jj) = q.Z(1);
                end
            end
        end

        % 找最优点
        [minZ, idxMin] = min(sl.Z);
        best = sl(idxMin,:);
        best_per_slice = [best_per_slice; best]; %#ok<AGROW>

        % 绘图
        fidx = fidx + 1;
        figs(fidx) = figure('Color','w');
        surf(X,Y,Z,'EdgeColor',[0.7 0.7 0.7]); hold on;
        zlim([zmin zmax]);
        caxis([zmin zmax]);   % colorbar 同一尺度
        scatter3(best.l2f, best.l1, best.Z, 80, 'filled','MarkerEdgeColor','k');
        % —— 这两行放在绘图循环的 surf()、scatter3(best...) 之后 ——
        [globalMin, gidx] = min(T.Z);
        global_best = T(gidx,:);

        isGlobal = best.l1==global_best.l1 & best.l2f==global_best.l2f & ...
            best.l3r==global_best.l3r & best.sf==global_best.sf;
        if isGlobal
            scatter3(best.l2f, best.l1, best.Z, 140, '^', 'filled', ...
                'MarkerEdgeColor','k'); % 红色三角表示“全局最优”
            text(best.l2f, best.l1, best.Z, '  GLOBAL BEST', ...
                'Color',[0.85 0 0], 'FontWeight','bold');
        end
        grid on; box on;
        xlabel('\lambda_2 factor');
        ylabel('\lambda_1');
        zlabel(zlab);
        title(sprintf('\\lambda_3 ratio = %.3g,  \\sigma_f = %.3g', l3s(i), sfs(j)));

        % 注释最佳点
        txt = sprintf('best: l1=%.3g, l2f=%.3g, Z=%.4g', best.l1, best.l2f, best.Z);
        text(best.l2f, best.l1, best.Z, ['  ' txt], 'FontSize', 10, 'Color','k');

        view(40,30); colormap(parula); colorbar; caxis([nanmin(Z(:)) nanmax(Z(:))]);

        if ~isempty(save_dir)
            fn = sprintf('terrain_l3r_%s_sf_%s_%s.png', num2str_nice(l3s(i)), num2str_nice(sfs(j)), lower(score_name));
            exportgraphics(figs(fidx), fullfile(save_dir, fn), 'Resolution', 200);
        end
    end
end

%% 全局最优（跨所有切片）
[globalMin, gidx] = min(T.Z);
global_best = T(gidx,:);

% 打印结果
fprintf('=== Global best (%s) ===\n', zlab);
fprintf('  l1=%.6g, l2f=%.6g, l3r=%.6g, sf=%.6g, Z=%.6g\n', ...
    global_best.l1, global_best.l2f, global_best.l3r, global_best.sf, global_best.Z);

% 每个切片最优
fprintf('=== Best per (l3r, sf) slice ===\n');
for k = 1:height(best_per_slice)
    r = best_per_slice(k,:);
    fprintf('  (l3r=%.6g, sf=%.6g) -> l1=%.6g, l2f=%.6g, Z=%.6g\n', ...
        r.l3r, r.sf, r.l1, r.l2f, r.Z);
end

% 可选导出 CSV
if ~isempty(save_dir)
    writetable(best_per_slice, fullfile(save_dir, sprintf('best_per_slice_%s.csv', lower(score_name))));
    writetable(global_best,     fullfile(save_dir, sprintf('global_best_%s.csv', lower(score_name))));
end

% 返回
out = struct('table',T,'best_per_slice',best_per_slice,'global_best',global_best,'figs',figs(~arrayfun(@isempty,figs)));

end % === 主函数 ===


%% ====== 工具函数 ======

function x = to_real_numeric(col)
% 把可能带小虚部的字符串列转成实数
if isnumeric(col)
    x = real(col);
    return;
end
if iscell(col), col = string(col); end
% 去掉 “ + a*i” 这样的虚部
col = regexprep(col, '\s*\+\s*[-+0-9.eE]+\s*i', '');
x = str2double(col);
end

function s = parse_key_to_params(keystr)
% 把 key 'l1_0_12__l2f_2_1__l3r_0_5__sf_4' 解析为数值
parts = split(string(keystr), "__");
m = struct('l1',NaN,'l2f',NaN,'l3r',NaN,'sf',NaN);
for k = 1:numel(parts)
    p = parts{k};
    if startsWith(p,"l1_")
        m.l1 = str2num_local(extractAfter(p,"l1_")); %#ok<ST2NM>
    elseif startsWith(p,"l2f_")
        m.l2f = str2num_local(extractAfter(p,"l2f_"));
    elseif startsWith(p,"l3r_")
        m.l3r = str2num_local(extractAfter(p,"l3r_"));
    elseif startsWith(p,"sf_")
        m.sf  = str2num_local(extractAfter(p,"sf_"));
    end
end
s = m;
end

function v = str2num_local(tok)
% 把 '0_12' -> 0.12；'4' -> 4
tok = strrep(string(tok), '_', '.');
v = str2double(tok);
end

function r = tiedrank(x, direction)
% 简单名次：direction=1 升序，-1 降序
if nargin<2, direction = 1; end
if direction<0, x = -x; end
[~,~,r] = unique(x);
% 把相同值的名次设为该组平均名次
r2 = zeros(size(x));
for v = unique(r).'
    idx = (r==v);
    r2(idx) = mean(find(idx));
end
r = r2;
end

function s = num2str_nice(v)
% 把 0.5 -> '0_5' 方便与 key 风格一致（仅用于文件名）
s = erase(string(strrep(num2str(v),'.','_')),' ');
end
