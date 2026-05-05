function out_files = export_jspace_audit_views(result_mat_path, audit, out_dir)
%EXPORT_JSPACE_AUDIT_VIEWS Export paired result MAT files for audit visualization.
%   This function does NOT change JSPACE outputs. It creates two auxiliary
%   MAT files that can be fed into visualize_subject_jspace_result:
%
%   1) FINAL-PAIR view:
%        Omega_est = original final Omega_est
%        Sjj_est   = inv(Omega_est)  (from audit.Sigma_from_Omega)
%
%   2) STATE-PAIR view:
%        Omega_est = inv(Sjj_est)    (from audit.Omega_from_Sjj)
%        Sjj_est   = original saved Sjj_est
%
%   Usage:
%       audit = audit_jspace_saved_result(result_mat_path, scalp_mat_path, [8 12]);
%       out_files = export_jspace_audit_views(result_mat_path, audit, out_dir);
%
%   Then visualize with your existing function:
%       out1 = visualize_subject_jspace_result(out_files.final_pair, scalp_mat_path, opts);
%       out2 = visualize_subject_jspace_result(out_files.state_pair, scalp_mat_path, opts);

    if nargin < 3 || isempty(out_dir)
        out_dir = fileparts(result_mat_path);
        if isempty(out_dir), out_dir = pwd; end
    end
    if ~isfolder(out_dir), mkdir(out_dir); end

    d = load(result_mat_path);
    [~, base, ~] = fileparts(result_mat_path);

    % ---------- FINAL-PAIR ----------
    d_final = d;
    d_final.Sjj_est = audit.Sigma_from_Omega;
    if ~isfield(d_final, 'audit_view') || ~isstruct(d_final.audit_view)
        d_final.audit_view = struct();
    end
    d_final.audit_view.mode = 'final_pair';
    d_final.audit_view.description = 'Omega_est is original final precision; Sjj_est is inv(Omega_est).';

    final_file = fullfile(out_dir, [base '_AUDIT_FINAL_PAIR.mat']);
    save(final_file, '-struct', 'd_final', '-v7.3');

    % ---------- STATE-PAIR ----------
    d_state = d;
    d_state.Omega_est = audit.Omega_from_Sjj;
    if ~isfield(d_state, 'audit_view') || ~isstruct(d_state.audit_view)
        d_state.audit_view = struct();
    end
    d_state.audit_view.mode = 'state_pair';
    d_state.audit_view.description = 'Sjj_est is original saved covariance state; Omega_est is inv(Sjj_est).';

    state_file = fullfile(out_dir, [base '_AUDIT_STATE_PAIR.mat']);
    save(state_file, '-struct', 'd_state', '-v7.3');

    out_files = struct();
    out_files.final_pair = final_file;
    out_files.state_pair = state_file;

    fprintf('\n[export_jspace_audit_views] Wrote:\n');
    fprintf('  final_pair : %s\n', final_file);
    fprintf('  state_pair : %s\n', state_file);
end
