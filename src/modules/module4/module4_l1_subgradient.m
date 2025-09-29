function [l1_subgradients, computation_stats] = module4_l1_subgradient(precision_matrices, params)
% MODULE4_L1_SUBGRADIENT  Compute L1 subgradients for complex Hermitian precision matrices (per-frequency).
%
% Inputs
%   precision_matrices : {F×1} cell, each p×p (complex allowed), intended Hermitian
%   params             : struct with fields (all optional)
%       .lambda2              (double, >=0)    default 0.01
%       .penalize_diagonal    (logical)        default false
%       .zero_threshold       (double)         default 1e-12
%       .force_hermitian      (logical)        default true
%       .subgradient_selection('zero'|'random'|'minimal'), default 'zero'
%       .verbose              (logical)        default false
%
% Outputs
%   l1_subgradients  : {F×1} cell, each p×p complex (Hermitian if force_hermitian=true)
%   computation_stats: struct with fields
%       .computation_time
%       .zero_entries_count     (Fx1)
%       .nonzero_entries_count  (Fx1)
%       .hermitian_violations   (Fx1)
%       .subgradient_norms      (Fx1)

% -------- Input checks --------
if nargin < 1
    error('module4_l1_subgradient:insufficient_input', 'precision_matrices is required');
end
if nargin < 2, params = struct(); end
if ~iscell(precision_matrices) || isempty(precision_matrices)
    error('module4_l1_subgradient:invalid_input', 'precision_matrices must be a non-empty cell array');
end

% -------- Defaults --------
defs = struct('lambda2', 0.01, ...
              'penalize_diagonal', false, ...
              'zero_threshold', 1e-12, ...
              'force_hermitian', true, ...
              'subgradient_selection', 'zero', ...
              'verbose', false);
fn = fieldnames(defs);
for k = 1:numel(fn)
    f = fn{k};
    if ~isfield(params, f), params.(f) = defs.(f); end
end
if params.lambda2 < 0
    error('module4_l1_subgradient:invalid_lambda2', 'lambda2 must be non-negative, got %.6f', params.lambda2);
end
valid_sel = {'zero','random','minimal'};
if ~ismember(params.subgradient_selection, valid_sel)
    error('module4_l1_subgradient:invalid_selection', ...
          'subgradient_selection must be one of: %s', strjoin(valid_sel, ', '));
end

% -------- Init --------
F = numel(precision_matrices);
p = size(precision_matrices{1}, 1);
l1_subgradients = cell(F,1);

computation_stats = struct();
computation_stats.zero_entries_count    = zeros(F,1);
computation_stats.nonzero_entries_count = zeros(F,1);
computation_stats.hermitian_violations  = zeros(F,1);
computation_stats.subgradient_norms     = zeros(F,1);

tic_main = tic;

% λ2 = 0 -> early return
if params.lambda2 == 0
    for f = 1:F
        l1_subgradients{f} = zeros(p, p);
    end
    computation_stats.computation_time = toc(tic_main);
    if params.verbose
        fprintf('  L1 subgradient: lambda2=0, returning zero subgradients\n');
    end
    return;
end

if params.verbose
    fprintf('Computing L1 subgradients (lambda2=%.6f)...\n', params.lambda2);
end

% -------- Per-frequency --------
for f = 1:F
    Gf = precision_matrices{f};
    if ~isnumeric(Gf) || ~isequal(size(Gf), [p, p])
        error('module4_l1_subgradient:invalid_matrix', ...
              'precision_matrices{%d} must be %dx%d numeric matrix', f, p, p);
    end

    % 只在上三角计算，再用 Hermitian 规则回填
    S = zeros(p,p);
    zero_cnt = 0; nonzero_cnt = 0;

    for i = 1:p
        for j = i:p   % 上三角含对角
            penalize = params.penalize_diagonal || (i ~= j);
            if ~penalize
                % 不惩罚的项（通常对角）
                S(i,j) = 0;
                continue;
            end

            gij = Gf(i,j);
            a   = abs(gij);

            if a > params.zero_threshold
                % 非零处：∂|z| = z/|z|
                S(i,j) = params.lambda2 * (gij / a);
                nonzero_cnt = nonzero_cnt + 1;
            else
                % 零邻域：|u| <= lambda2
                zero_cnt = zero_cnt + 1;
                switch params.subgradient_selection
                    case 'zero'
                        S(i,j) = 0;
                    case 'random'
                        ph = 2*pi*rand();
                        mg = params.lambda2 * rand();
                        S(i,j) = mg * exp(1i*ph);
                    case 'minimal'
                        S(i,j) = 0; % 最小范数选择
                end
            end
        end
    end

    % Hermitian 回填下三角
    if params.force_hermitian
        S_full = S;
        for i = 1:p
            for j = i+1:p
                S_full(j,i) = conj(S_full(i,j));
            end
        end
        % 最后做一次数值对称化，进一步抑制浮点误差
        S_full = 0.5*(S_full + S_full');
        herm_err = norm(S - S_full, 'fro');  % 仅作统计参考
        S = S_full;
    else
        herm_err = 0;
    end

    l1_subgradients{f} = S;
    computation_stats.zero_entries_count(f)    = zero_cnt;
    computation_stats.nonzero_entries_count(f) = nonzero_cnt;
    computation_stats.hermitian_violations(f)  = herm_err;
    computation_stats.subgradient_norms(f)     = norm(S, 'fro');

    if params.verbose && f <= 3
        fprintf('  Freq %d: %d zeros, %d nonzeros, ||subgrad||_F=%.6f\n', ...
                f, zero_cnt, nonzero_cnt, computation_stats.subgradient_norms(f));
    end
end

computation_stats.computation_time = toc(tic_main);
end
