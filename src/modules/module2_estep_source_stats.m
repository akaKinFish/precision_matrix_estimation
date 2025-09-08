function out = module2_estep_source_stats(input_data, estep_params)
% MODULE2_ESTEP_SOURCE_STATS
% Compute source-domain sufficient statistics for M-step:
%   S_hat_jj(f) = T_jv(f) * Sigma_vv_emp(f) * T_jv(f)' + Sigma_jj_post(f)

    if nargin < 2 || isempty(estep_params), estep_params = struct(); end

    ensure_hermitian = getfield_def(estep_params, 'ensure_hermitian', true);
    ensure_real_diag = getfield_def(estep_params, 'ensure_real_diag', true);
    ensure_psd       = getfield_def(estep_params, 'ensure_psd', false);
    psd_tol          = getfield_def(estep_params, 'psd_tol', 1e-10);
    diag_loading     = getfield_def(estep_params, 'diag_loading', 0);

    % ---- [ADDED] Normalize empirical_covariances BEFORE calling E-step ----
    must_have(input_data, 'empirical_covariances');
    [Sigma_vv_emp, Fnorm] = ensure_emp_cov_cell(input_data.empirical_covariances, input_data);

    in = input_data;
    in.empirical_covariances = Sigma_vv_emp;
    if ~isfield(in,'frequencies') || isempty(in.frequencies)
        in.frequencies = 1:Fnorm;
    else
        if numel(in.frequencies) ~= Fnorm
            % 如果给了单个 p×p 但 frequencies 长度是 F>1，则复制到 F 份
            if iscell(Sigma_vv_emp) && numel(Sigma_vv_emp)==1 && numel(in.frequencies) > 1
                Sigma_vv_emp = repmat(Sigma_vv_emp, numel(in.frequencies), 1);
                in.empirical_covariances = Sigma_vv_emp;
                Fnorm = numel(in.frequencies);
            else
                error('module2_estep_source_stats:freq_mismatch', ...
                      'frequencies length (%d) must match F=%d.', numel(in.frequencies), Fnorm);
            end
        end
    end

    % ---- Call existing E-step wrapper ----
    estep_results = module2_estep(in, estep_params);

    % ---- Pull T_jv and Sigma_jj_post (support canonical + legacy names) ----
    if isfield(estep_results, 'transfer_functions')
        Tcell = estep_results.transfer_functions;
    elseif isfield(estep_results, 'source_transfer_functions')
        Tcell = estep_results.source_transfer_functions;
    else
        error('module2_estep_source_stats:missing_T','Missing transfer functions in estep_results.');
    end

    if isfield(estep_results, 'posterior_source_covariances')
        SigPost = estep_results.posterior_source_covariances;
    elseif isfield(estep_results, 'posterior_covariances')
        SigPost = estep_results.posterior_covariances;
    else
        error('module2_estep_source_stats:missing_SigPost','Missing posterior source covariances.');
    end

    assert(numel(Tcell)==Fnorm, 'transfer_functions length=%d, expected F=%d', numel(Tcell), Fnorm);
    assert(numel(SigPost)==Fnorm,'posterior_source_covariances length=%d, expected F=%d', numel(SigPost), Fnorm);

    % ---- Dimensions from the first frequency ----
    [n_from_T, p_from_T] = size(Tcell{1});

    % ---- Compute S_hat_jj(f) ----
    Sjj_hat = cell(Fnorm,1);
    diagnostics = struct('hermitian_error', zeros(Fnorm,1), ...
                         'min_eig',         zeros(Fnorm,1), ...
                         'added_loading',   zeros(Fnorm,1));

    for f = 1:Fnorm
        T  = Tcell{f};        % n x p
        Sv = Sigma_vv_emp{f}; % p x p
        Sp = SigPost{f};      % n x n

        % sanity checks
        [nT,pT] = size(T); [pS1,pS2] = size(Sv); [nS1,nS2] = size(Sp);
        assert(nT==n_from_T && pT==p_from_T, 'T dims mismatch at f=%d', f);
        assert(pS1==p_from_T && pS2==p_from_T,'Sigma_vv dims mismatch at f=%d', f);
        assert(nS1==n_from_T && nS2==n_from_T,'Sigma_post dims mismatch at f=%d', f);

        % core
        S = T * Sv * T' + Sp;

        if diag_loading > 0
            S = S + diag_loading * eye(n_from_T, 'like', S);
            diagnostics.added_loading(f) = diagnostics.added_loading(f) + diag_loading;
        end
        if ensure_hermitian, S = 0.5*(S+S'); end
        if ensure_real_diag, S(1:n_from_T+1:end) = real(diag(S)); end

        if ensure_psd
            m = min(real(eig((S+S')/2)));
            diagnostics.min_eig(f) = m;
            if m < psd_tol
                add = (psd_tol - m);
                S = S + add*eye(n_from_T, 'like', S);
                diagnostics.added_loading(f) = diagnostics.added_loading(f) + add;
            end
        else
            diagnostics.min_eig(f) = NaN;
        end

        diagnostics.hermitian_error(f) = norm(S-S','fro')/max(1,norm(S,'fro'));
        Sjj_hat{f} = S;
    end

    % ---- Pack outputs ----
    out = struct();
    out.estep_results         = estep_results;
    out.source_second_moments = Sjj_hat;
    out.diagnostics           = diagnostics;
end

% ===================== Local helpers =====================

function val = getfield_def(S, fname, default)
    if isstruct(S) && isfield(S, fname) && ~isempty(S.(fname)), val = S.(fname);
    else, val = default;
    end
end

function must_have(S, fieldname)
    assert(isstruct(S) && isfield(S, fieldname), 'Missing field: %s', fieldname);
end

function validate_cell_of_square(C, name)
    assert(iscell(C) && ~isempty(C), '%s must be a non-empty cell array', name);
    for k = 1:numel(C)
        A = C{k};
        assert(isnumeric(A) || islogical(A), '%s{%d} must be numeric', name, k);
        s = size(A);
        assert(numel(s) == 2 && s(1) == s(2), '%s{%d} must be square. Got %dx%d', name, k, s(1), s(2));
    end
end

function [C, F] = ensure_emp_cov_cell(val, input_data)
% Normalize empirical covariances into {F x 1} cell of p x p.

    % A) already cell
    if iscell(val)
        validate_cell_of_square(val, 'input_data.empirical_covariances');
        C = val(:); F = numel(C); return;
    end

    % B) numeric 3D stack p x p x F
    if isnumeric(val) && ndims(val)==3 && size(val,1)==size(val,2)
        F = size(val,3); C = cell(F,1);
        for f=1:F, C{f}=val(:,:,f); end
        return;
    end

    % C) single p x p
    if isnumeric(val) && ismatrix(val) && size(val,1)==size(val,2)
        if isfield(input_data,'frequencies') && ~isempty(input_data.frequencies)
            F = numel(input_data.frequencies);
            C = repmat({val}, F, 1);
        else
            C = {val}; F = 1;
        end
        return;
    end

    % D) struct wrapper（常见误用：把 sim 或含字段的包传进来）
    if isstruct(val)
        keys = {'Sigma_emp','emp_covariance','Sigma_vv_emp','Sigma_vv_observed','Sigma_vv_true'};
        for k=1:numel(keys)
            if isfield(val, keys{k}) && ~isempty(val.(keys{k}))
                [C,F] = ensure_emp_cov_cell(val.(keys{k}), input_data); return;
            end
        end
    end

    error('module2_estep_source_stats:emp_cov_format', ...
          ['empirical_covariances must be cell{F,1} or p×p×F (numeric) or single p×p (with/without frequencies), ', ...
           'or a struct wrapping Sigma_emp/emp_covariance/Sigma_vv_emp/...']);
end
