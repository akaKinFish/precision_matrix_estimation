function figs = plot_jspace_audit_summary(audit, alpha_band)
%PLOT_JSPACE_AUDIT_SUMMARY Plot frequency-wise audit diagnostics.
%   figs = plot_jspace_audit_summary(audit, [8 12])

    if nargin < 2 || isempty(alpha_band), alpha_band = [8 12]; end

    f = audit.freq(:)';
    if isempty(f), f = 1:numel(audit.relerr_Sjj_vs_invOm); end

    figs = struct();

    % ---------- Figure 1: inverse-consistency ----------
    figs.inverse_consistency = figure('Name', 'Audit - Inverse Consistency', 'Position', [100 100 900 350]);
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    nexttile; hold on;
    patch([alpha_band(1) alpha_band(2) alpha_band(2) alpha_band(1)], [0 0 1 1], ...
        [1 0.92 0.92], 'EdgeColor','none', 'FaceAlpha',0.4);
    plot(f, audit.relerr_Sjj_vs_invOm, 'b-', 'LineWidth', 1.8);
    plot(f, audit.relerr_Sjj_vs_invOm, 'bo', 'MarkerSize', 4, 'MarkerFaceColor', 'b');
    xlabel('Frequency (Hz)'); ylabel('Relative error');
    title('relerr( Sjj\_est, inv(Omega\_est) )'); grid on; box on;

    nexttile; hold on;
    patch([alpha_band(1) alpha_band(2) alpha_band(2) alpha_band(1)], [0 0 1 1], ...
        [1 0.92 0.92], 'EdgeColor','none', 'FaceAlpha',0.4);
    semilogy(f, audit.cond_Omega, 'r-', 'LineWidth', 1.8);
    semilogy(f, audit.cond_Sjj,   'k-', 'LineWidth', 1.8);
    xlabel('Frequency (Hz)'); ylabel('Condition number');
    title('Numerical conditioning'); grid on; box on;
    legend({'alpha band','cond(Omega)','cond(Sjj)'}, 'Location','best');

    % ---------- Figure 2: coherence and pcor scales ----------
    figs.scale_compare = figure('Name', 'Audit - Coherence / PCOR Scale', 'Position', [120 120 900 350]);
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    nexttile; hold on;
    patch([alpha_band(1) alpha_band(2) alpha_band(2) alpha_band(1)], [0 0 1 1], ...
        [1 0.92 0.92], 'EdgeColor','none', 'FaceAlpha',0.4);
    plot(f, audit.coh_max_Sjj,   'b-', 'LineWidth', 1.8);
    plot(f, audit.coh_max_invOm, 'r-', 'LineWidth', 1.8);
    xlabel('Frequency (Hz)'); ylabel('Max linear coherence');
    title('Coherence scale'); grid on; box on;
    legend({'alpha band','from Sjj\_est','from inv(Omega\_est)'}, 'Location','best');

    nexttile; hold on;
    patch([alpha_band(1) alpha_band(2) alpha_band(2) alpha_band(1)], [0 0 1 1], ...
        [1 0.92 0.92], 'EdgeColor','none', 'FaceAlpha',0.4);
    plot(f, audit.pcor_max_Omega,  'b-', 'LineWidth', 1.8);
    plot(f, audit.pcor_max_invSjj, 'r-', 'LineWidth', 1.8);
    xlabel('Frequency (Hz)'); ylabel('Max PCOR');
    title('PCOR scale'); grid on; box on;
    legend({'alpha band','from Omega\_est','from inv(Sjj\_est)'}, 'Location','best');

    % ---------- Figure 3: common-mode strength ----------
    figs.common_mode = figure('Name', 'Audit - Common-Mode Strength', 'Position', [140 140 900 350]);
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    nexttile; hold on;
    patch([alpha_band(1) alpha_band(2) alpha_band(2) alpha_band(1)], [0 0 1 1], ...
        [1 0.92 0.92], 'EdgeColor','none', 'FaceAlpha',0.4);
    plot(f, audit.leadfrac_Sjj,   'b-', 'LineWidth', 1.8);
    plot(f, audit.leadfrac_invOm, 'r-', 'LineWidth', 1.8);
    xlabel('Frequency (Hz)'); ylabel('\lambda_{max} / tr(\Sigma)');
    title('Leading-mode fraction'); grid on; box on;
    legend({'alpha band','Sjj\_est','inv(Omega\_est)'}, 'Location','best');

    nexttile; hold on;
    patch([alpha_band(1) alpha_band(2) alpha_band(2) alpha_band(1)], [0 0 1 1], ...
        [1 0.92 0.92], 'EdgeColor','none', 'FaceAlpha',0.4);
    ratio_coh  = audit.coh_max_Sjj ./ max(audit.coh_max_invOm, 1e-12);
    ratio_pcor = audit.pcor_max_Omega ./ max(audit.pcor_max_invSjj, 1e-12);
    plot(f, ratio_coh,  'b-', 'LineWidth', 1.8);
    plot(f, ratio_pcor, 'r-', 'LineWidth', 1.8);
    xlabel('Frequency (Hz)'); ylabel('Scale ratio');
    title('State-vs-final scale ratios'); grid on; box on;
    legend({'alpha band','coh(Sjj)/coh(inv\Omega)','pcor(\Omega)/pcor(invSjj)'}, 'Location','best');
end
