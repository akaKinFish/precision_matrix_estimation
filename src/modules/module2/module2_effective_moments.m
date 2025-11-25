function [Psijj, T_xi_v, S_xi_xi, Sig_xi_post, Psixixi] = module2_effective_moments(Svv, L, Sigma_post, Tjv)
% BC-V E-step 有效矩：
% Psijj   = Σ_post + T Svv T'
% T_xi_v  = I - L T
% S_xi_xi = T_xi_v Svv T_xi_v'
% Σξξ,post= L Σ_post L'
% Psixixi = S_xi_xi + Σξξ,post

[p,~] = size(L);
I_p   = eye(p);

% Psijj
Psijj = Tjv * Svv * Tjv';
Psijj = (Psijj + Psijj')/2;
Psijj = Psijj + Sigma_post;

% Txi_v
T_xi_v = I_p - L * Tjv;
% 保持 Hermitian 可选：这里不强制

% S_xi_xi
S_xi_xi = module2_residual_empirical_covariance(T_xi_v, Svv, ...
    struct('enforce_hermitian',true,'regularization_factor',0, ...
           'min_eigenvalue_threshold',1e-12,'psd_only',true, ...
           'numerical_tolerance',1e-12,'verbose',false));

% Σξξ,post
Sig_xi_post = L * Sigma_post * L';
Sig_xi_post = (Sig_xi_post + Sig_xi_post')/2;

% Psixixi
Psixixi = S_xi_xi + Sig_xi_post;
end
