function pre = module1_preproc_from_covset(cov_cell, params)
% MODULE1_PREPROC_FROM_COVSET
% Reuse Module 1 preprocessing (whitening) **in the source domain** by
% packaging {F×1} source covariances as "simulation-like" inputs.

    % ---- Validate input covariance set ----
    assert(iscell(cov_cell) && ~isempty(cov_cell), 'cov_cell must be a non-empty cell array');
    F = numel(cov_cell);
    n = size(cov_cell{1}, 1);
    for f = 1:F
        A = cov_cell{f};
        assert(ismatrix(A) && size(A,1)==size(A,2), 'cov_cell{%d} must be square', f);
        assert(size(A,1)==n, 'cov_cell{%d} has inconsistent dimension', f);
    end

    % ---- Build "simulation-like" input for Module 1 ----
    input_data = struct();
    input_data.mode = 'simulation';
    input_data.sim_results = struct('Sigma_emp', {cov_cell}, 'n', n, 'F', F, 'T', 1);

    % ---- Parse params; map legacy names to current ones ----
    if nargin < 2 || isempty(params)
        nv = {};
    elseif isstruct(params)
        % map aliases
        if isfield(params,'min_power_threshold') && ~isfield(params,'min_power')
            params.min_power = params.min_power_threshold;   % alias -> official
            params = rmfield(params,'min_power_threshold');
        end
        if isfield(params,'diagonal_load') && ~isfield(params,'diagonal_loading')
            params.diagonal_loading = params.diagonal_load;  % alias -> official
            params = rmfield(params,'diagonal_load');
        end
        if isfield(params,'loading') && ~isfield(params,'loading_factor')
            params.loading_factor = params.loading;          % alias -> official
            params = rmfield(params,'loading');
        end
        nv = struct_to_nv(params);
    elseif iscell(params)
        % name-value already; best-effort map known aliases inline
        nv = params;
        nv = map_aliases_nv(nv);
    else
        error('params must be a struct or a name-value cell array');
    end

    % ---- Call Module 1 main; fallback to wrapper if needed ----
    try
        pre1 = module1_preprocessing_main(input_data, nv{:});
    catch
        % fallback through your public wrapper if main isn't visible
        pre1 = module1_preprocessing(input_data, cell2struct(nv(2:2:end), nv(1:2:end), 2));
    end

    % ---- Basic sanity checks and meta info ----
    must_have(pre1, 'Sigma_emp_loaded');
    must_have(pre1, 'D');
    must_have(pre1, 'Sigma_tilde');

    pre = pre1;
    pre.meta = struct('n', n, 'F', F, 'T', 1);
end

% ===================== Helpers =====================

function nv = struct_to_nv(S)
    f = fieldnames(S);
    nv = cell(1, 2*numel(f));
    for k = 1:numel(f)
        nv{2*k-1} = f{k};
        nv{2*k}   = S.(f{k});
    end
end

function nv2 = map_aliases_nv(nv)
% Map legacy names in a name-value cell to current Module 1 parameters.
    nv2 = nv;
    for k = 1:2:numel(nv)
        name = nv{k};
        if strcmpi(name,'min_power_threshold'), nv2{k} = 'min_power'; end
        if strcmpi(name,'diagonal_load'),      nv2{k} = 'diagonal_loading'; end
        if strcmpi(name,'loading'),            nv2{k} = 'loading_factor'; end
    end
end

function must_have(S, fieldname)
    assert(isstruct(S) && isfield(S, fieldname), 'Missing field: %s', fieldname);
end
