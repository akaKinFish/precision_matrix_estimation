function [Gamma_cells, results] = module5_stoch_fista_main(input_data, params)
% MODULE5_STOCH_FISTA_MAIN - Stochastic / Block FISTA for J-SPACE M-step (cell Hermitian SPD).
%
% Key design (as we agreed):
%   - Randomly sample frequency batch (band-stratified) each iter
%   - Optional neighbor-closure to stabilize cross-frequency smoothing
%   - Update only sampled frequencies; others frozen
%   - Cross-frequency smooth term uses full cached Gamma at current Y (exact, cheap)
%   - Fixed number of stochastic iterations (no early-stop by default)
%
% Modified based on feedback:
%   - Corrected local objective calculation for cross-freq terms (A1).
%   - Added switch to disable Nesterov momentum for stochastic stability (A2).
%   - Added L-parameter shrinking and clamping (A3).
%
% Compatible I/O with module5_fista_main:
%   input_data.whitened_covariances : {F x 1}
%   input_data.smoothing_kernel     : (F x F)
%   input_data.weight_matrix        : (p x p)
%   input_data.precision_matrices   : {F x 1} (warm start, optional)
%   input_data.active_mask          : {F x 1} (optional)
%
% params fields:
%   lambda1, lambda2, lambda3, weight_mode ('matrix'|'hadamard'), min_eig, alpha0, verbose
%   params.stoch.*:
%     .mode              : 'B' (band-stratified) | 'U' (uniform)
%     .freq              : 1xF frequency values (recommended)
%     .band_edges        : edges in Hz, e.g. [0 4 8 13 inf] (optional)
%     .band_groups       : cell of indices (optional, overrides band_edges)
%     .m_per_band        : #samples per band per iter (default 1)
%     .neighbor_closure  : 0/1/2 ... (default 1)
%     .max_iter          : stochastic iterations per call (default 40)
%     .seed              : RNG seed (default 0)
%     .use_backtracking  : true/false (default true)
%     .backtracking_factor : (default 2.0)
%     .max_backtracking  : (default 25)
%     .use_nesterov      : true/false (default false) [New A2]
%     .L_min             : Minimum L value (default 1e-3) [New A3]
%     .L_max             : Maximum L value (default 1e8) [New A3]
%     .L_shrink_on_success : Reduce L after successful step (default true) [New A3]
    % ---------------- Inputs ----------------
    Sigmas = input_data.whitened_covariances;
    Kernel = input_data.smoothing_kernel;
    W      = input_data.weight_matrix;
    F = numel(Sigmas);
    p = size(Sigmas{1}, 1);
    if isfield(input_data, 'precision_matrices') && ~isempty(input_data.precision_matrices)
        Gamma_curr = input_data.precision_matrices;
    else
        Gamma_curr = cell(F,1);
        for f = 1:F, Gamma_curr{f} = eye(p, 'like', Sigmas{1}); end
    end
    if isfield(input_data, 'active_mask') && ~isempty(input_data.active_mask)
        active_mask = input_data.active_mask;
    else
        active_mask = [];
    end
    if nargin < 2, params = struct(); end
    % ---------------- Defaults (core) ----------------
    params = set_default_(params, 'lambda1', 0);
    params = set_default_(params, 'lambda2', 0);
    params = set_default_(params, 'lambda3', 0);
    params = set_default_(params, 'weight_mode', 'hadamard'); 
    params = set_default_(params, 'min_eig', 1e-8);
    params = set_default_(params, 'alpha0', 0.1);
    params = set_default_(params, 'verbose', false);
    params = set_default_(params, 'cond_cap', 1e12);
    params = set_default_(params, 'rcond_min', 1e-12);
    params = set_default_(params, 'penalize_diagonal', false);
    params = set_default_(params, 'enforce_unit_diagonal', false);
    % ---------------- Defaults (stoch) ----------------
    if ~isfield(params, 'stoch'), params.stoch = struct(); end
    st = params.stoch;
    st = set_default_(st, 'mode', 'B');
    st = set_default_(st, 'm_per_band', 1);
    st = set_default_(st, 'neighbor_closure', 1);
    st = set_default_(st, 'max_iter', 40);
    st = set_default_(st, 'seed', 0);
    st = set_default_(st, 'use_backtracking', true);
    st = set_default_(st, 'backtracking_factor', 2.0);
    st = set_default_(st, 'max_backtracking', 25);
    
    % [B1] New defaults for stability
    st = set_default_(st, 'use_nesterov', false);      % A2: Default OFF for stochastic
    st = set_default_(st, 'L_min', 1e-3);             % A3: Minimum L
    st = set_default_(st, 'L_max', 1e8);              % A3: Maximum L (anti-explosion)
    st = set_default_(st, 'L_shrink_on_success', true); % A3: Allow L to decrease
    % band definition
    if isfield(st,'band_groups') && ~isempty(st.band_groups)
        band_groups = st.band_groups;
    else
        freq = [];
        if isfield(st,'freq') && ~isempty(st.freq), freq = st.freq; end
        if isempty(freq)
            % fallback: uniform groups (single band)
            band_groups = {1:F};
        else
            edges = [0 4 8 13 inf];
            if isfield(st,'band_edges') && ~isempty(st.band_edges)
                edges = st.band_edges(:).';
            end
            band_groups = build_band_groups_(freq, edges);
            if isempty(band_groups)
                band_groups = {1:F};
            end
        end
    end
    % local RNG stream
    stream = RandStream('twister','Seed', st.seed);
    % ---------------- Init: ensure Hermitian SPD ----------------
    for f = 1:F
        Gamma_curr{f} = utils_math.make_hermitian(Gamma_curr{f});
        [Gamma_curr{f}, ~] = utils_math.project_spd(Gamma_curr{f}, params.min_eig);
    end
    Y_curr = Gamma_curr;
    t_curr = 1;
    Ksym = (Kernel + Kernel')/2;
    % initial L for backtracking
    if params.alpha0 > 0
        L_curr = 1 / params.alpha0;
    else
        L_curr = 10;
    end
    % [A3] Initial Clamp
    L_curr = max(st.L_min, min(L_curr, st.L_max));
    % For histories
    batch_trace = cell(st.max_iter, 1);
    alpha_trace = zeros(st.max_iter, 1);
    bt_trace    = zeros(st.max_iter, 1);
    % ---------------- Main loop (fixed iters) ----------------
    for iter = 1:st.max_iter
        % ---- (A) sample batch indices ----
        batch = sample_batch_(band_groups, st.m_per_band, st.mode, F, stream);
        batch = neighbor_closure_(batch, st.neighbor_closure, F);
        batch = unique(batch(:).');
        
        % [A1] Create membership lookup for this batch (fast O(1) access)
        inB = false(1, F);
        inB(batch) = true;
        batch_trace{iter} = batch;
        % ---- (B) gradient at Y (ONLY for batch freqs) ----
        grads = compute_grad_batch_(Y_curr, Sigmas, Ksym, W, params, batch);
        % ---- (C) (optional) local backtracking on smooth part ----
        bt_count = 0;
        accepted = false;
        % Precompute f(Y) local (smooth only)
        % [A1] Updated call signature to include inB
        fY_local = smooth_obj_local_(Y_curr, Sigmas, Ksym, W, params, batch, inB);
        while ~accepted
            alpha = 1 / L_curr;
            
            % Candidate = Y - alpha * grad, then prox (ONLY batch updated)
            Gamma_cand = Gamma_curr; % start from current (frozen others)
            
            for ii = 1:numel(batch)
                f = batch(ii);
                Z = Y_curr{f} - alpha * grads{f};
                mask_f = [];
                if ~isempty(active_mask), mask_f = active_mask{f}; end
                
                % prox for L1-off + projections (Hermitian + SPD floor)
                [Gamma_cand{f}, ~] = module_proximal_operator_fista.compute( ...
                    Z, params.lambda2 * alpha, mask_f, params);
            end
            if ~st.use_backtracking
                accepted = true;
            else
                % majorization check: f(Gcand) <= f(Y) + <grad,Delta> + (L/2)||Delta||^2
                % [A1] Updated call signature to include inB
                f_new_local = smooth_obj_local_(Gamma_cand, Sigmas, Ksym, W, params, batch, inB);
                
                Delta = cell(F,1);
                for ii = 1:numel(batch)
                    f = batch(ii);
                    Delta{f} = Gamma_cand{f} - Y_curr{f};
                end
                
                diff_norm_sq = cell_norm_sq_batch_(Delta, batch);
                lin_term     = cell_inner_batch_(grads, Delta, batch);
                majorant     = fY_local + lin_term + (L_curr/2) * diff_norm_sq;
                if f_new_local <= majorant * (1 + 1e-12)
                    accepted = true;
                else
                    L_curr = L_curr * st.backtracking_factor;
                    bt_count = bt_count + 1;
                    if bt_count >= st.max_backtracking
                        % emergency accept with very small step
                        accepted = true;
                    end
                end
            end
        end
        
        % [A3] L-step logic: Shrink if successful immediately, and always clamp
        if st.use_backtracking
            if bt_count == 0 && st.L_shrink_on_success
                L_curr = L_curr / st.backtracking_factor;
            end
            % Clamp L to prevent underflow/overflow
            L_curr = min(max(L_curr, st.L_min), st.L_max);
        end
        % store step stats
        alpha_trace(iter) = 1 / L_curr;
        bt_trace(iter)    = bt_count;
        % ---- (D) Nesterov update (ONLY batch; freeze others) ----
        % [A2] Configurable Nesterov
        if st.use_nesterov
            t_next = (1 + sqrt(1 + 4 * t_curr^2)) / 2;
            Y_next = Y_curr;   % freeze by default
            for ii = 1:numel(batch)
                f = batch(ii);
                Y_next{f} = Gamma_cand{f} + ((t_curr - 1)/t_next) * (Gamma_cand{f} - Gamma_curr{f});
            end
            t_curr = t_next;
        else
            % Standard Proximal Gradient (no momentum) - safer for stochastic
            Y_next = Y_curr;
            for ii = 1:numel(batch)
                f = batch(ii);
                Y_next{f} = Gamma_cand{f};
            end
            t_curr = 1; % Reset momentum counter
        end
        Gamma_curr = Gamma_cand;
        Y_curr     = Y_next;
        if params.verbose && (iter == 1 || mod(iter,10)==0)
            fprintf('[StochFISTA] iter=%d/%d | |B|=%d | alpha=%.2e | bt=%d | Nes=%d\n', ...
                iter, st.max_iter, numel(batch), alpha_trace(iter), bt_count, st.use_nesterov);
        end
    end
    Gamma_cells = Gamma_curr;
    
    % [B3] Enhanced Results
    results = struct();
    results.batch_trace = batch_trace;
    results.alpha_trace = alpha_trace;
    results.bt_trace    = bt_trace;
    results.final_iter  = st.max_iter;
    results.seed        = st.seed;
    results.mode        = st.mode;
    results.m_per_band  = st.m_per_band;
    results.neighbor_closure = st.neighbor_closure;
    % Diagnostics
    results.L_final     = L_curr;
    results.alpha_final = 1/L_curr;
    results.used_nesterov = st.use_nesterov;
    results.L_min       = st.L_min;
    results.L_max       = st.L_max;
end
% =====================================================================
% Helpers
% =====================================================================
function s = set_default_(s, field, val)
    if ~isfield(s, field) || isempty(s.(field))
        s.(field) = val;
    end
end
function band_groups = build_band_groups_(freq, edges)
    % edges: e.g. [0 4 8 13 inf]
    freq = freq(:).';
    F = numel(freq);
    band_groups = {};
    for k = 1:(numel(edges)-1)
        lo = edges(k);
        hi = edges(k+1);
        idx = find(freq >= lo & freq < hi);
        if isempty(idx), continue; end
        band_groups{end+1} = idx; %#ok<AGROW>
    end
    % ensure coverage
    covered = false(1,F);
    for k=1:numel(band_groups)
        covered(band_groups{k}) = true;
    end
    rest = find(~covered);
    if ~isempty(rest)
        band_groups{end+1} = rest;
    end
end
function batch = sample_batch_(band_groups, m_per_band, mode, F, stream)
    if nargin < 4, F = 0; end %#ok<NASGU>
    if strcmpi(mode,'U')
        % uniform
        K = numel(band_groups);
        m = max(1, K*m_per_band);
        all_idx = 1:F;
        rp = randperm(stream, F, min(m, F));
        batch = all_idx(rp);
        return;
    end
    % band-stratified
    batch = [];
    for k = 1:numel(band_groups)
        idx = band_groups{k};
        nk = numel(idx);
        if nk <= 0, continue; end
        mk = min(m_per_band, nk);
        rp = randperm(stream, nk, mk);
        batch = [batch, idx(rp)]; %#ok<AGROW>
    end
end
function batch2 = neighbor_closure_(batch, radius, F)
    if radius <= 0
        batch2 = batch;
        return;
    end
    batch = unique(batch(:).');
    acc = batch;
    for r = 1:radius
        acc = [acc, acc-1, acc+1]; %#ok<AGROW>
    end
    acc = unique(acc);
    acc = acc(acc>=1 & acc<=F);
    batch2 = acc;
end
function grads = compute_grad_batch_(Gamma_cells, Sigma_cells, Ksym, W, params, batch)
    F = numel(Gamma_cells);
    p = size(Gamma_cells{1},1);
    grads = cell(F,1);
    
    lambda1 = params.lambda1;
    lambda3 = params.lambda3;
    mode    = params.weight_mode;
    
    degrees = sum(Ksym, 2);
    
    for ii = 1:numel(batch)
        f = batch(ii);
        G = Gamma_cells{f};
        S = Sigma_cells{f};
        
        % --- robust SPD ---
        G0 = utils_math.make_hermitian(G);
        [G_pd, ~] = utils_math.project_spd(G0, params.min_eig);
        
        rc = rcond(G_pd);
        if rc < params.rcond_min
            G_pd = project_spd_cond_clip_(G0, params.min_eig, params.cond_cap);
        end
        
        invG = inv_spd_chol_(G_pd);
        grad_fit = -invG + S;
        
        % --- cross-frequency smoothing gradient ---
        grad_smooth = zeros(p,p,'like',G);
        if lambda1 > 0
            % NeighborSum = sum_{f2} K(f,f2) * Gamma_{f2}
            NeighborSum = zeros(p,p,'like',G);
            rowK = Ksym(f,:);
            nz = find(rowK ~= 0);
            for kk = 1:numel(nz)
                f2 = nz(kk);
                NeighborSum = NeighborSum + rowK(f2) * Gamma_cells{f2};
            end
            DiffTerm = degrees(f) * G - NeighborSum; 
            if strcmp(mode,'matrix')
                grad_smooth = 2 * lambda1 * (W * DiffTerm);
            elseif strcmp(mode,'hadamard')
                grad_smooth = 2 * lambda1 * (W .* DiffTerm);
            else
                error('module5_stoch_fista_main:UnknownWeightMode','weight_mode must be matrix or hadamard');
            end
        end
        
        % --- single-frequency smoothing gradient ---
        grad_space = zeros(p,p,'like',G);
        if lambda3 > 0
            if strcmp(mode,'matrix')
                grad_space = 2 * lambda3 * (W * G);
            elseif strcmp(mode,'hadamard')
                grad_space = 2 * lambda3 * (W .* G);
            else
                error('module5_stoch_fista_main:UnknownWeightMode','weight_mode must be matrix or hadamard');
            end
        end
        
        grads{f} = utils_math.make_hermitian(grad_fit + grad_smooth + grad_space);
    end
end
% [A1] & [B2] Modified: Added inB for weight logic and nz for speed
function fval = smooth_obj_local_(Gamma_cells, Sigma_cells, Ksym, W, params, batch, inB)
    % Local smooth objective for backtracking.
    % CORRECTED LOGIC:
    % Term = lambda1 * 0.5 * K(f,f2) * ||Gf-Gf2||^2
    % If f2 is ALSO in batch: we count it here, and later when we visit f2.
    % So weight = 0.5 is correct (sum of two halves = 1).
    % If f2 is NOT in batch: we visit this edge only once (now).
    % So weight must be 1.0 to account for the full edge penalty in the global objective.
    
    lambda1 = params.lambda1;
    lambda3 = params.lambda3;
    mode    = params.weight_mode;
    
    fval = 0;
    
    % data terms
    for ii = 1:numel(batch)
        f = batch(ii);
        G = Gamma_cells{f};
        S = Sigma_cells{f};
        [ld, ok] = utils_math.safe_log_det(G);
        if ~ok
            fval = Inf;
            return;
        end
        fval = fval - ld + real(trace(S * G));
    end
    
    % cross-frequency smoothing term 
    if lambda1 > 0
        for ii = 1:numel(batch)
            f = batch(ii);
            rowK = Ksym(f,:);
            nz = find(rowK ~= 0); % [B2] Only iterate non-zero neighbors
            
            for kk = 1:numel(nz)
                f2 = nz(kk);
                kff = rowK(f2);
                
                % [A1] CRITICAL FIX: Determine weight based on batch membership
                if inB(f2)
                    weight = 0.5; % Edge will be visited twice (f->f2 and f2->f)
                else
                    weight = 1.0; % Edge visited only once (f->f2)
                end
                
                Gdiff = Gamma_cells{f} - Gamma_cells{f2};
                
                if strcmp(mode,'matrix')
                    term = real(trace(Gdiff' * (W * Gdiff)));
                else
                    term = real(sum(sum(conj(Gdiff) .* (W .* Gdiff))));
                end
                
                fval = fval + lambda1 * weight * kff * term;
            end
        end
    end
    
    % spatial term
    if lambda3 > 0
        for ii = 1:numel(batch)
            f = batch(ii);
            G = Gamma_cells{f};
            if strcmp(mode,'matrix')
                term = real(trace(G' * (W * G)));
            else
                term = real(sum(sum(conj(G) .* (W .* G))));
            end
            fval = fval + lambda3 * term;
        end
    end
end
function val = cell_norm_sq_batch_(Delta, batch)
    val = 0;
    for ii = 1:numel(batch)
        f = batch(ii);
        if isempty(Delta{f}), continue; end
        val = val + norm(Delta{f}, 'fro')^2;
    end
end
function val = cell_inner_batch_(A, B, batch)
    val = 0;
    for ii = 1:numel(batch)
        f = batch(ii);
        if isempty(A{f}) || isempty(B{f}), continue; end
        val = val + real(sum(sum(conj(A{f}) .* B{f})));
    end
end
function A_clip = project_spd_cond_clip_(A, min_eig, cond_cap)
    A = (A + A')/2;
    [V,d] = eig(A,'vector');
    d = real(d);
    d = max(d, min_eig);
    d = min(d, min_eig * cond_cap);
    A_clip = V * diag(d) * V';
    A_clip = (A_clip + A_clip')/2;
end
function invA = inv_spd_chol_(A)
    A = (A + A')/2;
    [R,flag] = chol(A);
    if flag ~= 0
        tau = 1e-6 * trace(A)/size(A,1);
        A2 = A + tau*eye(size(A),'like',A);
        [R,flag2] = chol(A2);
        if flag2 ~= 0
            invA = A \ eye(size(A,1),'like',A);
            return;
        end
        R = chol(A2);
    end
    I = eye(size(A,1),'like',A);
    invA = R \ (R' \ I);
end