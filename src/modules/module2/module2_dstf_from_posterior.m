function T_jv = module2_dstf_from_posterior(L, Sigma_post, Sigma_xi_xi)
% T_jv = Σ_post * L' * Σξ^{-1}

T_jv = Sigma_post * (L' / Sigma_xi_xi);
if isreal(L) && isreal(Sigma_post) && isreal(Sigma_xi_xi)
    if max(abs(imag(T_jv(:)))) < 1e-13, T_jv = real(T_jv); end
end
end
