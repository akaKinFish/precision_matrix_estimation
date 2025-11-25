function [D_cell, S_tilde_cell] = module2_whitening_from_covset(Sjj_cell, whitening_opts)
% 将 Psijj → 白化后的 Sjj_tilde，并返回 D_src
if nargin<2, whitening_opts = struct(); end
pre = module1_preproc_from_covset(Sjj_cell, whitening_opts);
D_cell       = pre.D;
S_tilde_cell = pre.Sigma_tilde;
end
