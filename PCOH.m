
%% ======================= 工具函数区 =======================
function PCOH(Omega_true_in, varargin)
% 画 ground-truth partial coherence（每个频点一幅热力图）
% 支持输入 {F×1} cell 或 n×n×F 数组
% 可选参数（Name,Value）：
%   'freq_labels' : 1×F 元胞（每幅子图标题用）
%   'clim'        : [lo hi] 颜色范围，默认 [0 1]
%   'colormap'    : 颜色图，默认 parula
%   'fig_title'   : 总标题，默认 'Ground-truth Partial Coherence'
%   'show_colorbar': true/false，默认 true

% ---- 解析参数 ----
p = inputParser;
addParameter(p,'freq_labels',[]);
addParameter(p,'clim',[0 1]);
addParameter(p,'colormap',parula(256));
addParameter(p,'fig_title','Ground-truth Partial Coherence');
addParameter(p,'show_colorbar',true);
parse(p,varargin{:});
opt = p.Results;

% ---- 统一成 cell{F,1} ----
[Om_cell, n, F] = coerce_prec_cell(Omega_true_in);

% ---- 计算每个频点的 partial coherence ----
PC = cell(F,1);
for f = 1:F
    Om = (Om_cell{f} + Om_cell{f}')/2;          % Hermitian
    d  = real(diag(Om));
    % 数值地板：避免除0/负噪
    d  = max(d, max(1e-12, 1e-12*max(d)));
    D  = sqrt(d) * sqrt(d).';
    pc = abs(-Om) ./ (D + eps);
    pc(1:n+1:end) = 0;                           % 对角置 0（自相关不画）
    PC{f} = pc;
end

% ---- 打印一些指标 ----
fprintf('\n[GT Partial Coherence] per-frequency stats:\n');
for f = 1:F
    pc = PC{f};
    ut = triu(true(n),1);
    vals = pc(ut);
    vals = vals(isfinite(vals));
    fprintf('  f=%d: max=%6.4f | median=%6.4f | 95%%=%6.4f\n', ...
        f, max(vals), median(vals), prctile(vals,95));
end

% ---- 作图 ----
ncol = min(F, 4); nrow = ceil(F / ncol);
fig = figure('Name', 'GT Partial Coherence', 'Color','w'); %#ok<NASGU>
tiledlayout(nrow, ncol, 'Padding','compact','TileSpacing','compact');

for f = 1:F
    nexttile;
    imagesc(PC{f});
    axis image off;
    colormap(opt.colormap);
    % if ~isempty(opt.clim), caxis(opt.clim); end
    title(get_label(opt,f),'FontWeight','bold');
    if opt.show_colorbar, colorbar; end
end
sgtitle(opt.fig_title, 'FontWeight','bold');

end

% ----------------- helpers -----------------
function [C, n, F] = coerce_prec_cell(Om)
    if iscell(Om)
        F = numel(Om); n = size(Om{1},1);
        C = Om(:);
    elseif isnumeric(Om) && ndims(Om)==3 && size(Om,1)==size(Om,2)
        [n,~,F] = size(Om);
        C = cell(F,1);
        for f=1:F, C{f} = Om(:,:,f); end
    else
        error('Expect {F×1} cell or n×n×F array for Omega_true.');
    end
end

function s = get_label(opt,f)
    if ~isempty(opt.freq_labels) && numel(opt.freq_labels)>=f
        s = opt.freq_labels{f};
    else
        s = sprintf('f=%d', f);
    end
end
