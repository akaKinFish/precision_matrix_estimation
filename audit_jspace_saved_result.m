function audit = audit_jspace_saved_result(result_mat_path, alpha_band, eps_floor)
    if nargin < 2 || isempty(alpha_band), alpha_band = [8 12]; end
    if nargin < 3 || isempty(eps_floor), eps_floor = 1e-12; end

    d = load(result_mat_path);
    if ~isfield(d, 'Omega_est') || ~isfield(d, 'Sjj_est')
        error('Result file must contain Omega_est and Sjj_est.');
    end

    Omega = to_tensor_(d.Omega_est);
    Sjj   = to_tensor_(d.Sjj_est);

    Omega = sanitize_cross_tensor_(Omega, eps_floor);
    Sjj   = sanitize_cross_tensor_(Sjj,   eps_floor);

    [Nr, ~, F] = size(Omega);

    % optional freq
    if isfield(d, 'freq') && ~isempty(d.freq)
        freq = d.freq(:)';
        freq = freq(1:min(numel(freq), F));
    else
        freq = 1:F;
    end

    alpha_idx = find(freq >= alpha_band(1) & freq <= alpha_band(2));
    if isempty(alpha_idx), alpha_idx = 1:F; end

    % ----- build paired covariance from final Omega -----
    Sigma_from_Omega = zeros(Nr, Nr, F);
    Omega_from_Sjj   = zeros(Nr, Nr, F);

    relerr_Sjj_vs_invOm = zeros(F,1);
    inv_resid_Omega     = zeros(F,1);
    inv_resid_Sjj       = zeros(F,1);
    cond_Omega          = zeros(F,1);
    cond_Sjj            = zeros(F,1);
    min_eig_Omega       = zeros(F,1);
    min_eig_Sjj         = zeros(F,1);

    coh_max_Sjj         = zeros(F,1);
    coh_max_invOm       = zeros(F,1);
    pcor_max_Omega      = zeros(F,1);
    pcor_max_invSjj     = zeros(F,1);

    leadfrac_Sjj        = zeros(F,1);
    leadfrac_invOm      = zeros(F,1);

    for f = 1:F
        % ----- Omega side -----
        Om = (Omega(:,:,f) + Omega(:,:,f)')/2;
        [V,D] = eig(Om);
        ev = real(diag(D));
        ev(ev < 1e-8) = 1e-8;
        Om_spd = V * diag(ev) * V';
        Om_spd = (Om_spd + Om_spd')/2;

        Sigma_f = inv(Om_spd);
        Sigma_f = (Sigma_f + Sigma_f')/2;
        Sigma_from_Omega(:,:,f) = Sigma_f;

        inv_resid_Omega(f) = norm(Om_spd*Sigma_f - eye(Nr), 'fro') / sqrt(Nr);
        cond_Omega(f) = cond(Om_spd);
        min_eig_Omega(f) = min(real(eig(Om_spd)));

        % ----- Sjj side -----
        Sj = (Sjj(:,:,f) + Sjj(:,:,f)')/2;
        [V2,D2] = eig(Sj);
        ev2 = real(diag(D2));
        ev2(ev2 < 1e-8) = 1e-8;
        Sj_spd = V2 * diag(ev2) * V2';
        Sj_spd = (Sj_spd + Sj_spd')/2;

        Om_from_Sj = inv(Sj_spd);
        Om_from_Sj = (Om_from_Sj + Om_from_Sj')/2;
        Omega_from_Sjj(:,:,f) = Om_from_Sj;

        inv_resid_Sjj(f) = norm(Sj_spd*Om_from_Sj - eye(Nr), 'fro') / sqrt(Nr);
        cond_Sjj(f) = cond(Sj_spd);
        min_eig_Sjj(f) = min(real(eig(Sj_spd)));

        % ----- direct mismatch between saved Sjj and inv(final Omega) -----
        relerr_Sjj_vs_invOm(f) = norm(Sj_spd - Sigma_f, 'fro') / max(norm(Sigma_f, 'fro'), 1e-12);

        % ----- coherence and pcor summaries -----
        C1 = coherence_from_one_(Sj_spd, eps_floor);
        C2 = coherence_from_one_(Sigma_f, eps_floor);
        P1 = pcor_from_one_precision_(Om_spd, eps_floor);
        P2 = pcor_from_one_precision_(Om_from_Sj, eps_floor);

        coh_max_Sjj(f)    = max_offdiag_(C1);
        coh_max_invOm(f)  = max_offdiag_(C2);
        pcor_max_Omega(f) = max_offdiag_(P1);
        pcor_max_invSjj(f)= max_offdiag_(P2);

        % ----- leading eigenvalue fraction: detect strong common-mode structure -----
        e1 = sort(real(eig(Sj_spd)), 'descend');
        e2 = sort(real(eig(Sigma_f)), 'descend');
        leadfrac_Sjj(f)   = e1(1) / max(sum(e1), eps_floor);
        leadfrac_invOm(f) = e2(1) / max(sum(e2), eps_floor);
    end

    % ---- alpha-band summaries ----
    alpha = alpha_idx(:);

    audit = struct();
    audit.freq = freq;
    audit.alpha_idx = alpha_idx;
    audit.Sigma_from_Omega = Sigma_from_Omega;
    audit.Omega_from_Sjj   = Omega_from_Sjj;

    audit.relerr_Sjj_vs_invOm = relerr_Sjj_vs_invOm;
    audit.inv_resid_Omega = inv_resid_Omega;
    audit.inv_resid_Sjj   = inv_resid_Sjj;
    audit.cond_Omega      = cond_Omega;
    audit.cond_Sjj        = cond_Sjj;
    audit.min_eig_Omega   = min_eig_Omega;
    audit.min_eig_Sjj     = min_eig_Sjj;

    audit.coh_max_Sjj      = coh_max_Sjj;
    audit.coh_max_invOm    = coh_max_invOm;
    audit.pcor_max_Omega   = pcor_max_Omega;
    audit.pcor_max_invSjj  = pcor_max_invSjj;
    audit.leadfrac_Sjj     = leadfrac_Sjj;
    audit.leadfrac_invOm   = leadfrac_invOm;

    fprintf('\n===== JSPACE RESULT AUDIT =====\n');
    fprintf('median relerr( Sjj_est , inv(Omega_est) ) over alpha = %.3e\n', median(relerr_Sjj_vs_invOm(alpha)));
    fprintf('median inv residual of Omega_est over alpha           = %.3e\n', median(inv_resid_Omega(alpha)));
    fprintf('median inv residual of Sjj_est over alpha             = %.3e\n', median(inv_resid_Sjj(alpha)));

    fprintf('median cond(Omega_est) over alpha                     = %.3e\n', median(cond_Omega(alpha)));
    fprintf('median cond(Sjj_est) over alpha                       = %.3e\n', median(cond_Sjj(alpha)));

    fprintf('median max coherence from Sjj_est over alpha          = %.3e\n', median(coh_max_Sjj(alpha)));
    fprintf('median max coherence from inv(Omega_est) over alpha   = %.3e\n', median(coh_max_invOm(alpha)));

    fprintf('median max pcor from Omega_est over alpha             = %.3e\n', median(pcor_max_Omega(alpha)));
    fprintf('median max pcor from inv(Sjj_est) over alpha          = %.3e\n', median(pcor_max_invSjj(alpha)));

    fprintf('median leading-eig fraction of Sjj_est over alpha     = %.3e\n', median(leadfrac_Sjj(alpha)));
    fprintf('median leading-eig fraction of inv(Omega_est) alpha   = %.3e\n', median(leadfrac_invOm(alpha)));
end

function T = to_tensor_(X)
    if iscell(X), T = cat(3, X{:}); else, T = X; end
end

function T = sanitize_cross_tensor_(T, eps_floor)
    [N,~,F] = size(T);
    for f = 1:F
        Sf = (T(:,:,f) + T(:,:,f)')/2;
        Sf(1:N+1:end) = max(real(diag(Sf)), eps_floor);
        T(:,:,f) = Sf;
    end
end

function C = coherence_from_one_(S, eps_floor)
    N = size(S,1);
    p = max(real(diag(S)), eps_floor);
    C = (abs(S).^2) ./ (p * p.');
    C(1:N+1:end) = 1;
    C = min(max(C,0),1);
end

function P = pcor_from_one_precision_(Om, eps_floor)
    N = size(Om,1);
    d = max(real(diag(Om)), eps_floor);
    P = abs(-Om ./ sqrt(d*d.'));
    P(1:N+1:end) = 0;
end

function m = max_offdiag_(A)
    N = size(A,1);
    mask = triu(true(N),1);
    v = A(mask);
    v = v(isfinite(v));
    if isempty(v), m = NaN; else, m = max(v); end
end