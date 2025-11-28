function recoloring_results = module8_recoloring(input_data, recoloring_params)
% MODULE8_RECOLORING - Recoloring wrapper with robust input normalization.
%
% Inputs:
%   input_data.whitened_precision_matrices : {F×1} cell | p×p | p×p×F
%   input_data.whitening_matrices          : {F×1} cell | p×p | p×p×F | p(×)1(diag)
%   （兼容别名：'Gamma_tilde','D','D_src'）
%
% Outputs:
%   recoloring_results.recolored_precision_matrices : {F×1} cell of Ω_f = D_f * Γ̃_f * D_f
%   recoloring_results.success : logical
%   recoloring_results.details : struct (F,p 及归一化信息)

    if nargin < 1, error('module8_recoloring:insufficient_input','input_data is required'); end
    if nargin < 2, recoloring_params = struct(); end

    % -------- Normalize & validate inputs to cell-of-matrices ----------
    try
        [in_data, norminfo] = normalize_recoloring_input_(input_data);
    catch ME
        error('module8_recoloring:bad_input','%s', ME.message);
    end

    % -------- Route to main or fallback core ----------
    try
        if exist('module8_recoloring_main','file') == 2
            % 把"已归一化"的 in_data 传给主实现
            recoloring_results = module8_recoloring_main(in_data, recoloring_params);
            % 兼容：补充 success 字段
            if ~isfield(recoloring_results,'success'), recoloring_results.success = true; end
        else
            % 兜底：本地实现 Ω = D * Γ̃ * D
            recoloring_results = recolor_core_(in_data);
        end
    catch ME_original
        % 包装更友好的错误
        try
            ME = MException('module8_recoloring:computation_failed', ...
                            'Module 8 recoloring computation failed.');
            ME = addCause(ME, ME_original);
            throw(ME);
        catch
            error('module8_recoloring:computation_failed', ...
                  'Module 8 recoloring computation failed: %s', ME_original.message);
        end
    end

    % -------- Ensure required outputs & attach details ----------
    if ~isfield(recoloring_results, 'recolored_precision_matrices')
        error('module8_recoloring:missing_output', ...
              'Expected recolored_precision_matrices field not found in results');
    end
    recoloring_results.success = true;
    % 附带一次归一化信息，便于排查
    recoloring_results.details = norminfo;
end

% =====================================================================
% =============== helpers =============================================
% =====================================================================

function [out, info] = normalize_recoloring_input_(S)
    % 取字段（兼容别名）
    G = pick_first_(S, {'whitened_precision_matrices','Gamma_tilde','Gamma_tilde_star','G'});
    D = pick_first_(S, {'whitening_matrices','D','D_src','whiteners'});

    if isempty(G)
        error('normalize_recoloring_input:missingGamma', ...
              'Missing field ''whitened_precision_matrices'' (or alias).');
    end
    if isempty(D)
        error('normalize_recoloring_input:missingD', ...
              'Missing field ''whitening_matrices'' (or alias).');
    end

    % 统一转 cell{F,1}
    Gc = coerce_to_cell_square_(G);
    Dc = coerce_to_cell_square_(D);

    Fg = numel(Gc); Fd = numel(Dc);
    if Fg ~= Fd
        % 若某一侧只有 1 个，允许广播
        if Fg == 1 && Fd > 1
            Gc = repmat(Gc, Fd, 1);
            Fg = Fd;
        elseif Fd == 1 && Fg > 1
            Dc = repmat(Dc, Fg, 1);
            Fd = Fg;
        else
            error('normalize_recoloring_input:freq_mismatch', ...
                  'Frequency count mismatch: |Gamma|=%d, |D|=%d', Fg, Fd);
        end
    end

    % 尺寸一致性 & Hermitian/实对角
    pG = size(Gc{1},1); pD = size(Dc{1},1);
    if pG ~= pD
        error('normalize_recoloring_input:size_mismatch', ...
              'Size mismatch: Gamma is %dx%d, D is %dx%d', pG,pG,pD,pD);
    end
    p = pG;
    for f=1:Fg
        Gc{f} = symfix_(Gc{f});
        Dc{f} = expand_diag_if_needed_(Dc{f}, p);
        Dc{f} = symfix_(Dc{f});
        % 数值下界以免奇异
        Dc{f} = spd_floor_(Dc{f}, 1e-12);
        Gc{f} = spd_floor_(Gc{f}, 0);       % Γ̃ 不必 SPD，但对称即可
    end

    out = struct();
    out.whitened_precision_matrices = Gc;
    out.whitening_matrices          = Dc;

    info = struct('F',Fg,'p',p,'notes','inputs normalized to cell; hermitian enforced');
end

function S = symfix_(S)
    S = 0.5*(S + S');
    % 对角取实数，兼容复 Hermitian
    S(1:size(S,1)+1:end) = real(diag(S));
end

function A = spd_floor_(A, eps_ld)
    % 仅将非常接近奇异的对角进行抬升，避免彻底改变矩阵
    % 这里不强制 SPD（Γ̃ 可半正定/不定）；只在 D 上使用 eps_ld>0。
    if eps_ld<=0, return; end
    [V,d] = eig(full(0.5*(A+A')),'vector');
    d = real(d);
    d(d < eps_ld) = eps_ld;
    A = V*diag(d)*V';
    A = 0.5*(A + A');
end

function C = coerce_to_cell_square_(X)
    if iscell(X)
        C = X(:);
        return;
    end
    if isnumeric(X)
        if ndims(X)==2 && size(X,1)==size(X,2)
            C = {X};
            return;
        elseif ndims(X)==3 && size(X,1)==size(X,2)
            F = size(X,3);
            C = cell(F,1);
            for f=1:F, C{f} = X(:,:,f); end
            return;
        end
    end
    error('coerce_to_cell_square:badType', ...
          'Expect cell{F,1} or p×p or p×p×F numeric.');
end

function D = expand_diag_if_needed_(D, p)
    % 允许把向量/对角提法转成完整方阵
    if isvector(D) && numel(D)==p
        D = diag(D(:));
    end
    if ~ismatrix(D) || any(size(D)~=[p p])
        error('expand_diag_if_needed:badSize','whitener must be p×p (or length-p diag).');
    end
end

function recol = recolor_core_(in_data)
    Gc = in_data.whitened_precision_matrices;
    Dc = in_data.whitening_matrices;
    F  = numel(Gc); p = size(Gc{1},1);
    Omega = cell(F,1);
    for f=1:F
        G = symfix_(Gc{f});
        D = symfix_(Dc{f});
        % 经典回色：Ω = D * Γ̃ * D
        Of = D * G * D;
        Of = symfix_( Of );
        % 为稳健可做轻微 SPD floor（仅正则到极小正特征）
        Of = spd_floor_(Of, 1e-12);
        Omega{f} = Of;
    end
    recol = struct();
    recol.recolored_precision_matrices = Omega;
    recol.success = true;
    recol.details = struct('F',F,'p',p,'impl','fallback_core');
end

function v = pick_first_(S, names)
    v = [];
    for k=1:numel(names)
        if isfield(S, names{k}) && ~isempty(S.(names{k}))
            v = S.(names{k});
            return;
        end
    end
end
