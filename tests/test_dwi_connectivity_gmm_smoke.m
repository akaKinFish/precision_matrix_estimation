function test_dwi_connectivity_gmm_smoke()
% Fast smoke test for DWI connectivity GMM masking.

    rng(0);
    p = 10;
    C = rand(p);
    C = triu(C, 1);
    C = C + C';
    C(1:p+1:end) = 0;
    C(1:3, 1:3) = C(1:3, 1:3) + 5;
    C(1:p+1:end) = 0;

    cfg = struct('dwi_gmm_log', false);
    [mask, info] = utils_build_dwi_connectivity_mask_gmm(C, cfg);

    assert(all(diag(mask)), 'Diagonal should be true.');
    assert(isequal(mask, mask.'), 'Mask should be symmetric.');
    off = mask & ~eye(p);
    assert(nnz(off) > 0, 'Expected at least one off-diagonal true.');
    assert(nnz(off) < p * (p - 1), 'Expected mask to be sparse.');
    assert(isfinite(info.threshold), 'Threshold should be finite.');
end
