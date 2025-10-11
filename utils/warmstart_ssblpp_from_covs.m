function warm = warmstart_ssblpp_from_covs(Svv_cell, L, opts)
% WARMSTART_SSBLPP_FROM_COVS
% 用 sSSBLpp / sSSBLpp_ultralarge 生成 EM 的 warm-start:
%   warm.Omega_init{f}  = inv_psd_robust(Sjj_e{f})
%   warm.Sjj_e{f}       = Tjv * Svv_scaled * Tjv'
%   warm.Tjv{f}         = 传递算子 (n×p)
%
% Inputs
%   Svv_cell : {F×1}, 每项 p×p 传感器互谱（复厄米）
%   L        : p×n lead field
%   opts.m   : 有效样本数（缺省用  max(p, 20) ）
%   opts.W   : q×q 权重矩阵（若不给或太大则自动走 _ultralarge 版本）
%   opts.parcellation : {a×1}，每个元素是索引向量；若不给则每个节点单独成区
%   opts.variant : 'auto' | 'standard' | 'ultralarge'
%
% Output struct warm 同 eLORETA 版本，字段完全对齐：
%   .Omega_init{f}, .Sjj_e{f}, .Tjv{f}, .info{f}

    if nargin < 3, opts = struct(); end
    F = numel(Svv_cell);
    [p,n] = size(L);

    % 选项
    m   = getf(opts, 'm', max(p,20));            % 有效样本数
    W   = getf(opts, 'W', []);                   % 若为空→走 ultralarge
    parc= getf(opts, 'parcellation', default_parcellation(n));
    var = lower(getf(opts, 'variant', 'auto'));  % 'auto'|'standard'|'ultralarge'

    use_ul = strcmp(var,'ultralarge') || (strcmp(var,'auto') && (isempty(W) || numel(W)>5e6));

    Omega_init = cell(F,1);
    Sjj_e      = cell(F,1);
    Tjv_cell   = cell(F,1);
    info_cell  = cell(F,1);

    for f = 1:F
        Svv = hermi_(Svv_cell{f});   % 保证厄米
        if use_ul
            % --- 不带 W 的大规模版本 ---
            [~,~,Tjv,Svv_scaled,~,scaleJ,scaleLvj] = sSSBLpp_ultralarge(Svv, L, m, parc);
            used_variant = 'ultralarge';
        else
            % --- 标准版（带 W）---
            if isempty(W), W = speye(n); end
            [~,~,Tjv,Svv_scaled,~,scaleJ,scaleLvj] = sSSBLpp(Svv, L, m, W, parc);
            used_variant = 'standard';
        end

        % 源协方差估计（用 sSSBL 内部的缩放后 Svv）
        Sjj = hermi_( Tjv * Svv_scaled * Tjv' );

        % 数值稳健投影 + 求逆，得到初值 Ω^(0)
        Sjj  = psd_project_(Sjj);
        Om0  = inv_psd_robust_(Sjj, 1e-8, 1e-12);

        Omega_init{f} = Om0;
        Sjj_e{f}      = Sjj;
        Tjv_cell{f}   = Tjv;

        info_cell{f}  = struct('variant',used_variant, ...
                               'scaleJ',scaleJ,'scaleLvj',scaleLvj);
    end

    warm = struct('Omega_init',{Omega_init}, ...
                  'Sjj_e',{Sjj_e}, ...
                  'Tjv',{Tjv_cell}, ...
                  'info',{info_cell});
end

% --------- tiny utils (与工程里已有重名时可删去) ---------
function A = hermi_(A), A = (A + A')/2; end
function S = psd_project_(S)
    S = hermi_(S);
    [U,d] = eig(full(S),'vector'); d = max(real(d), 0);
    S = U*diag(d)*U'; S = hermi_(S);
end
function Om = inv_psd_robust_(A, eps_reg, min_ratio)
    A = hermi_(A);
    [V,D] = eig(full(A)); d = real(diag(D)); dmax = max(d);
    floor_val = max(min_ratio * max(dmax, eps), 0);
    d(d<floor_val) = floor_val;
    if eps_reg > 0, d = (d + eps_reg * dmax) / (1 + eps_reg); end
    Om = V * diag(1./d) * V'; Om = hermi_(Om);
end
function parc = default_parcellation(n)
    parc = arrayfun(@(i) i, 1:n, 'uni', 0);
end
function v = getf(s,f,def)
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = def; end
end
