function judge_multi_metric(csvPath, varargin)
% 综合评估 summary.csv：多指标分组加权的分位秩汇总 + 可视化
% Name-Value:
%   'SaveDir'      : 导出图像目录（空则不导出）
%   'GroupWeights' : [wDist, wID, wStab]
%   'IDWeights'    : [AUPR, F1, AUROC, Prec, Rec]
%   'ScoreName'    : 输出打印/图上显示的分数名（默认 'balanced'）

%% -------- 参数 --------
p = inputParser;
addParameter(p,'SaveDir','',@ischar);
addParameter(p,'GroupWeights',[0.5 0.35 0.15],@(x)isnumeric(x)&&numel(x)==3);
addParameter(p,'IDWeights',[0.4 0.35 0.15 0.05 0.05],@(x)isnumeric(x)&&numel(x)==5);
addParameter(p,'ScoreName','balanced',@ischar);
parse(p,varargin{:});
savedir = strtrim(p.Results.SaveDir);
wGroup  = p.Results.GroupWeights(:).';
wID     = p.Results.IDWeights(:).';
sname   = p.Results.ScoreName;

if ~exist(csvPath,'file'), error('找不到 CSV：%s',csvPath); end
T = readtable(csvPath);

%% -------- 取超参列（兼容别名；若缺少则尝试从 key 解析） --------
l1  = getParam(T, {'l1','lambda1','Lambda1'});
l2f = getParam(T, {'l2f','lambda2_fac','lambda2fac','Lambda2_fac'});
l3r = getParam(T, {'l3r','lambda3_ratio','lambda3','Lambda3_ratio'});
sf  = getParam(T, {'sf','sigma_f','Sigma_f'});

%% -------- 可用的指标列（有就用） --------
distCols = intersect({'LE_mean','AIRM_mean','JBLD_mean','BW_mean','KLS_mean'}, T.Properties.VariableNames,'stable');
idCols   = intersect({'AUPR_mean','F1_mean','AUROC_mean','Prec_mean','Rec_mean'}, T.Properties.VariableNames,'stable');
stabCols = intersect({'Jaccard_mean'}, T.Properties.VariableNames,'stable');

if isempty(distCols)
    error('表中未找到任何距离列（LE/AIRM/JBLD/BW/KLS 的 *_mean）');
end
if isempty(idCols)
    warning('识别指标列为空（AUPR/F1/AUROC/Prec/Rec 的 *_mean）—将只用距离+稳定性。');
end
if isempty(stabCols)
    warning('稳定性列 Jaccard_mean 缺失—将只用距离+识别。');
end

% 组内权重（如果少列，自动重标化）
wDist = ones(1,numel(distCols)); wDist = wDist/sum(wDist);
if ~isempty(idCols)
    mapID = containers.Map( {'AUPR_mean','F1_mean','AUROC_mean','Prec_mean','Rec_mean'}, num2cell(wID) );
    wIDv = zeros(1,numel(idCols));
    for i=1:numel(idCols), wIDv(i) = mapID(idCols{i}); end
    if sum(wIDv)==0, wIDv = ones(size(wIDv)); end
    wIDv = wIDv/sum(wIDv);
else
    wIDv = [];
end

% 组间权重（缺组时重标化）
present = [~isempty(distCols), ~isempty(idCols), ~isempty(stabCols)];
wGroup = wGroup .* present;
wGroup = wGroup / sum(wGroup);

%% -------- 把各指标转成 0~1 分位秩（0=最好，1=最差）并汇总 --------
N = height(T);
rankMat = [];  % 用于 Pareto/相关分析（统一方向：越小越好）
labels  = {};

% 距离（越小越好）
for i=1:numel(distCols)
    labels{end+1} = distCols{i}; %#ok<AGROW>
    rankMat(:,end+1) = fracRank(T.(distCols{i}), 'ascend'); %#ok<AGROW>
end
distScore = sum(rankMat(:,end-numel(distCols)+1:end).*wDist, 2); % 组内加权平均

% 识别（越大越好 -> 取负号再做秩）
if ~isempty(idCols)
    R_id = zeros(N, numel(idCols));
    for i=1:numel(idCols)
        labels{end+1} = idCols{i}; %#ok<AGROW>
        R_id(:,i) = fracRank(T.(idCols{i}), 'descend'); % descend=>大好→小秩
    end
    idScore = sum(R_id.*wIDv, 2);
    rankMat = [rankMat, R_id];
else
    idScore = zeros(N,1);
end

% 稳定性（越大越好）
if ~isempty(stabCols)
    R_stab = fracRank(T.(stabCols{1}), 'descend');
    labels{end+1} = stabCols{1};
else
    R_stab = zeros(N,1);
end
% ---------------- 综合分 ----------------
S = wGroup(1)*distScore + wGroup(2)*idScore + wGroup(3)*R_stab;

%% -------- 打印全局最优 & 每切片最优 --------
[~,ix] = min(S);
fprintf('=== Global best (%s score, ↓更好) ===\n', sname);
fprintf('  l1=%.4g, l2f=%.4g, l3r=%.4g, sf=%.4g | %s=%.3f\n', l1(ix), l2f(ix), l3r(ix), sf(ix), sname, S(ix));

% 每 (l3r, sf) 切片
pairs = unique([l3r, sf],'rows');
fprintf('=== Best per (l3r, sf) slice ===\n');
for k=1:size(pairs,1)
    m = (l3r==pairs(k,1) & sf==pairs(k,2));
    [v,j] = min(S(m));
    idx = find(m); idx = idx(j);
    fprintf('  (l3r=%.4g, sf=%.4g) -> l1=%.4g, l2f=%.4g | %s=%.3f\n', pairs(k,1), pairs(k,2), l1(idx), l2f(idx), sname, v);
end

% 撞边检测
msg = boundaryMsg(l1,l2f,l3r,sf,ix);
fprintf('=== Boundary check ===\n%s\n', msg);

% Top-10
[~,ord] = sort(S,'ascend');
K = min(10,N);
fprintf('=== Top-%d breakdown ===\n', K);
for t=1:K
    i = ord(t);
    fprintf('#%d: l1=%.4g, l2f=%.4g, l3r=%.4g, sf=%.4g | %s=%.3f | Dist=%.3f, ID=%.3f, Jacc=%.3f\n', ...
        t, l1(i), l2f(i), l3r(i), sf(i), sname, S(i), distScore(i), idScore(i), R_stab(i));
end

%% -------- 可视化 --------
% 1) (l1,l2f)->S 的二维地形（每个 l3r,sf 切片）
uL3 = unique(l3r); uSF = unique(sf);
tl = tiledlayout(numel(uL3), numel(uSF), 'TileSpacing','compact','Padding','compact');
title(tl, sprintf('%s terrain over (\\lambda_1,\\lambda_2fac)', sname));
allZ=[]; AX=[];
for i=1:numel(uL3)
    for j=1:numel(uSF)
        ax = nexttile; AX(end+1)=ax; %#ok<AGROW>
        m  = (l3r==uL3(i) & sf==uSF(j));
        if nnz(m)<3, axis off; text(0.5,0.5,'N/A','HorizontalAlignment','center'); continue; end
        [xg,yg,zg] = grid2surf(l1(m), l2f(m), S(m));
        surf(ax,xg,yg,zg); view(35,35); shading interp; grid on;
        xlabel('\lambda_1'); ylabel('\lambda_2 factor'); zlabel([sname,' (↓)']);
        title(ax, sprintf('l3r=%.3g, sf=%.3g', uL3(i), uSF(j)));
        hold on; % 标注该切片最优
        [vmin, jj] = min(S(m)); idx = find(m); idx=idx(jj);
        plot3(l1(idx), l2f(idx), vmin, 'ko', 'MarkerFaceColor','y', 'MarkerSize',6);
        allZ = [allZ; zg(:)]; %#ok<AGROW>
    end
end
if ~isempty(allZ)
    zlim_all = [min(allZ,[],'omitnan'), max(allZ,[],'omitnan')];
    for k=1:numel(AX), if isgraphics(AX(k)), caxis(AX(k), zlim_all); end, end
    colormap turbo; colorbar;
end

% 2) 单维趋势（跨另一维取最优）：Z vs l2f / l1 / l3r
figure('Name','Best across \lambda_1: score vs \lambda_2 factor'); hold on; grid on;
[l2u,~,ic] = unique(l2f); zbest = accumBest(S, ic);
plot(sort(l2u), zbest, '-o'); xlabel('\lambda_2 factor'); ylabel([sname,' (best across \lambda_1)']);

figure('Name','Best across \lambda_2: score vs \lambda_1'); hold on; grid on;
[l1u,~,ic] = unique(l1); zbest = accumBest(S, ic);
plot(sort(l1u), zbest, '-o'); xlabel('\lambda_1'); ylabel([sname,' (best across \lambda_2)']);

figure('Name','Best across (\lambda_1,\lambda_2): score vs \lambda_3'); hold on; grid on;
[l3u,~,ic] = unique(l3r); zbest = accumBest(S, ic);
plot(sort(l3u), zbest, '-o'); xlabel('\lambda_3 ratio'); ylabel([sname,' (best across \lambda_1,\lambda_2)']);

% 3) 3D 散点
figure('Name','3D scatter colored by score');
scatter3(l1,l2f,l3r,36,S,'filled'); grid on; colorbar; colormap turbo;
xlabel('\lambda_1'); ylabel('\lambda_2 factor'); zlabel('\lambda_3 ratio');
title([sname,' (↓ better)']);

% 4) 指标相关性（Spearman）
M = []; nm = {};
for i=1:numel(distCols), M=[M, T.(distCols{i})]; nm{end+1}=distCols{i}; end %#ok<AGROW>
for i=1:numel(idCols),   M=[M, T.(idCols{i})];   nm{end+1}=idCols{i};   end %#ok<AGROW>
for i=1:numel(stabCols), M=[M, T.(stabCols{i})]; nm{end+1}=stabCols{i}; end %#ok<AGROW>
if ~isempty(M)
    R = corr(M,'Type','Spearman','Rows','pairwise');
    figure('Name','Spearman correlation of raw metrics');
    imagesc(R); axis image; colorbar; colormap turbo;
    xticks(1:numel(nm)); xticklabels(nm); xtickangle(45);
    yticks(1:numel(nm)); yticklabels(nm);
    title('Metric correlation (Spearman)');
end

% 保存
if ~isempty(savedir)
    if ~exist(savedir,'dir'), mkdir(savedir); end
    figs = findobj('Type','figure');
    for k=1:numel(figs)
        nm = get(figs(k),'Name'); if isempty(nm), nm=sprintf('fig_%d',k); end
        print(figs(k), fullfile(savedir,[matlab.lang.makeValidName(nm),'.png']), '-dpng','-r150');
    end
    fprintf('[SAVE] 图像已保存到 %s\n', savedir);
end
end

%% ========= 辅助 =========
function r = fracRank(x, mode)
% 0~1 分位秩（0=最佳），mode='ascend'(小好)/'descend'(大好)
x = x(:);
keep = isfinite(x);
r = nan(size(x));
if all(~keep), r(:)=0.5; return; end
xk = x(keep);
if strcmpi(mode,'ascend')
    s = tiedrank(xk);  % 小的排前
else
    s = tiedrank(-xk); % 大的排前
end
r(keep) = (s-1) / max(numel(s)-1,1);
r(~keep) = 0.5; % 缺失给中位秩
end

function z = accumBest(S, ic)
% 按组取最小 S（↓更好），并按组标签排序输出
z = accumarray(ic, S, [], @min);
z = z(sort(unique(ic)));
end

function [xg,yg,zg] = grid2surf(x,y,z)
xu = unique(x); yu = unique(y);
[xg,yg] = ndgrid(xu, yu);
zg = nan(size(xg));
for i=1:numel(xu)
    for j=1:numel(yu)
        m = (x==xu(i) & y==yu(j));
        if any(m), zg(i,j) = min(z(m)); end % 同格取最小分数
    end
end
end

function msg = boundaryMsg(l1,l2f,l3r,sf,ix)
fmt = @(v) sprintf('[%.4g, %.4g]', min(v), max(v));
hit = {};
if l1(ix)==min(l1)  || l1(ix)==max(l1),   hit{end+1}='λ1';  end %#ok<AGROW>
if l2f(ix)==min(l2f)|| l2f(ix)==max(l2f), hit{end+1}='λ2';  end %#ok<AGROW>
if l3r(ix)==min(l3r)|| l3r(ix)==max(l3r), hit{end+1}='λ3';  end %#ok<AGROW>
if sf(ix)==min(sf)  || sf(ix)==max(sf),   hit{end+1}='σ_f'; end %#ok<AGROW>
msg = sprintf('Ranges: λ1%s, λ2%s, λ3%s, σ_f%s\n', fmt(l1),fmt(l2f),fmt(l3r),fmt(sf));
if isempty(hit), msg = [msg, 'No boundary hit.']; else, msg=[msg, 'Boundary hit at: ', strjoin(hit,', ')]; end
end

function v = getParam(T, names)
for i=1:numel(names)
    if ismember(names{i}, T.Properties.VariableNames)
        v = T.(names{i}); return;
    end
end
% 尝试从 key 解析：l1_...__l2f_...__l3r_...__sf_...
if ismember('key',T.Properties.VariableNames)
    v = nan(height(T),1);
    pat = strjoin(regexprep(names,'_','\\_'), '|'); %#ok<NASGU>
    % 针对 l1/l2f/l3r/sf 四类分别解析
    if any(strcmpi(names, 'l1')|strcmpi(names,'lambda1')|strcmpi(names,'Lambda1'))
        v = parseFromKey(T.key,'l1');
        return;
    elseif any(contains(lower(names),'l2'))
        v = parseFromKey(T.key,'l2f'); return;
    elseif any(contains(lower(names),'l3'))
        v = parseFromKey(T.key,'l3r'); return;
    elseif any(contains(lower(names),'sf'))
        v = parseFromKey(T.key,'sf'); return;
    end
end
error('未找到列：%s，且无法从 key 解析', strjoin(names,', '));
end

function v = parseFromKey(keyCell, tag)
v = nan(numel(keyCell),1);
for i=1:numel(keyCell)
    s = string(keyCell{i});
    % 形如 l1_1_2__l2f_10__l3r_1_3__sf_4 → 把下划线恢复成小数点
    expr = sprintf('%s_([0-9_\\.]+)', tag); % 捕获 tag 后的数字串
    m = regexp(s, expr, 'tokens', 'once');
    if ~isempty(m)
        v(i) = str2double(strrep(m{1}, '_', '.'));
    end
end
end
