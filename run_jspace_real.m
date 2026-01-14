function [Omega_est, Sjj_est, outs_js] = run_jspace_real(Svv_cross, L, freq, dwi_C, cfg_overrides)
% Minimal wrapper for J-SPACE on real data

cfg_js = struct();
cfg_js.max_em_iter   = 10;
cfg_js.verbose       = true;
cfg_js.use_gpu       = false;

cfg_js.m_samples     = 1000;
cfg_js.optimizer     = 'stoch_fista';
cfg_js.freq          = freq;

% Soft prior (optional)
if nargin >= 4 && ~isempty(dwi_C)
    cfg_js.dwi_C = dwi_C;
    cfg_js.dwi_weight_mode   = 'power';
    cfg_js.dwi_weight_alpha  = 1.0;
    cfg_js.dwi_weight_clip   = [0.25, 4];
    cfg_js.dwi_weight_normalize_median = true;
end

% Stochastic config (optional but recommended if you use _stoch solver)
cfg_js.stoch = struct();
cfg_js.stoch.mode = 'B';
cfg_js.stoch.seed = 0;
cfg_js.stoch.n_per_band = 2;
cfg_js.stoch.Nsfreq = cfg_js.stoch.n_per_band;
cfg_js.stoch.max_total  = 16;
cfg_js.stoch.band_edges = [0 4 8 13 20];
cfg_js.stoch.max_iter   = 40;
cfg_js.stoch.tol        = 0;
cfg_js.stoch.backtracking_beta = 0.5;
cfg_js.stoch.max_backtracking  = 20;
cfg_js.stoch.monotone          = true;
cfg_js.stoch.use_restart       = true;
cfg_js.stoch.rescale_kernel_mass = true;

cfg_js.opt_max_evals   = 80;     % 先加一点
cfg_js.opt_min_points  = 10;
cfg_js.opt_use_parallel = true;

% ---- eta bounds (looser lbs) ----
cfg_js.eta1_lb = 1e-4;  cfg_js.eta1_ub = 3e-2;
cfg_js.eta2_lb = 1e-4;  cfg_js.eta2_ub = 3e-1;
cfg_js.eta3_lb = 1e-4;  cfg_js.eta3_ub = 5e-1;

% ---- multi init points in log10 domain ----
lb = [cfg_js.eta1_lb cfg_js.eta2_lb cfg_js.eta3_lb];
ub = [cfg_js.eta1_ub cfg_js.eta2_ub cfg_js.eta3_ub];

init_eta = [
    1e-2  5e-2  1e-2
    5e-3  1e-1  2e-2
    1e-3  5e-2  5e-3
];

% LHS fill to min_points
need = max(0, cfg_js.opt_min_points - size(init_eta,1));
if need>0
    X = lhsdesign(need,3);
    init_eta = [init_eta; 10.^(log10(lb) + X.*(log10(ub)-log10(lb)))];
end

init_eta = min(max(init_eta, lb), ub);
cfg_js.opt_initial_points_log = log10(init_eta);

% 让 x0 是第一行（和你打印一致）
cfg_js.eta1_init = init_eta(1,1);
cfg_js.eta2_init = init_eta(1,2);
cfg_js.eta3_init = init_eta(1,3);


% Override cfg if provided
if nargin >= 5 && ~isempty(cfg_overrides)
    fn = fieldnames(cfg_overrides);
    for k=1:numel(fn)
        cfg_js.(fn{k}) = cfg_overrides.(fn{k});
    end
end

% --- sanity fixes (recommended) ---
Ne = size(Svv_cross,1);
Nw = size(Svv_cross,3);

for i = 1:Nw
    Si = Svv_cross(:,:,i);

    % Hermitianize
    Si = (Si + Si')/2;

    % force real diagonal
    Si(1:Ne+1:end) = real(diag(Si));

    Svv_cross(:,:,i) = Si;
end


% Call solver
[Omega_est, Sjj_est, outs_js] = solver_jspace_3d_opt_stoch(Svv_cross, L, [], cfg_js);

end
