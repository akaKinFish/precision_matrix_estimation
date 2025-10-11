function metrics = evaluate_all_methods(Omega_true, Omega_all, opts, Dsrc, Gamma_tilde_star)
% EVALUATE_ALL_METHODS  对多个方法统一评测（基于你的 gt_compare_and_plot）
%
% Omega_all: struct with fields matching method names, each {F×1}

if nargin < 3, opts = struct(); end
if ~isfield(opts,'f_view'), opts.f_view = 1; end

mnames = fieldnames(Omega_all);
metrics = struct();
for i = 1:numel(mnames)
    name = mnames{i};
    try
        fprintf('[EVAL] %s ...\n', name);
        % 传入白化矩阵与 Gamma_tilde_star 时，函数会额外计算“白化域”一致性指标
        metrics.(name) = gt_compare_and_plot(Omega_true, Omega_all.(name), opts, Dsrc, Gamma_tilde_star);
        saveas(gcf, sprintf('results/fig_metrics_%s.png', name));
    catch ME
        warning('Evaluation failed for %s: %s', name, ME.message);
    end
end
save('results/metrics_all.mat','-struct','metrics');
end
