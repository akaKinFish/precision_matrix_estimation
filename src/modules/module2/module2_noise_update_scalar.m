function sigma2 = module2_noise_update_scalar(Psixixi_cell)
% 频率平均迹/p 作为标量噪声估计
F = numel(Psixixi_cell);
if F==0, sigma2 = []; return; end
p = size(Psixixi_cell{1},1);
acc = 0;
for f=1:F
    acc = acc + trace(Psixixi_cell{f})/p;
end
sigma2 = real(acc / F);
sigma2 = max(sigma2, 1e-10);
end
