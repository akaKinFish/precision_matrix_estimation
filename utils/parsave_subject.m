function parsave_subject(out_file, subject_id, age, site_dir, site_name, Omega_est, Sjj_est, outs_js)
    save(out_file, 'subject_id', 'age', 'site_dir', 'site_name', ...
         'Omega_est', 'Sjj_est', 'outs_js', '-v7.3');
end