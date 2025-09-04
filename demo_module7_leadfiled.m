% MODULE7_LEADFIELD_VIS_DEMO - Visual demo for Module 7 with leadfield (source & sensor views)
%
% This script generates a complete set of figures to *see* what Module 7 is
% producing in both source-space and sensor-space:
%   1) Geometry: 3D electrodes + sources + head spheres
%   2) Sensor covariances: heatmaps of Σ_vv_true and Σ_vv_observed across selected freqs
%   3) Source connectome: 3D network from Ω_source_true at a chosen frequency
%   4) Leadfield diagnostics: column norms, singular values, condition number
%   5) Forward topography: project a sparse source pattern to sensors (3D colored electrodes)
%   6) Frequency traces: trace(Σ_vv_true), trace(Σ_noise), achieved SNR_dB per frequency
%
% Requirements:
%   - module7_simulation_improved_complex.m (new API with leadfield)
%
% Notes:
%   - All comments are in English, per your request.
%   - The demo is robust to complex-valued covariances (uses magnitudes where needed).
%   - No dependency on Module 2; this is a standalone visualization demo.

clear; clc; close all;

fprintf('========================================\n');
fprintf('Module 7 Leadfield Visualization Demo\n');
fprintf('========================================\n\n');

%% ------------------------- Simulation Settings ---------------------------
cfg = struct();
cfg.n_nodes        = 32;            % = n_sources (source-space dimensionality)
cfg.n_sensors      = 19;            % EEG sensors (10-20), p can be changed to 32/64 etc.
cfg.n_freq         = 12;            % number of frequency bins
cfg.n_samples      = 200;           % samples per frequency (for empirical covariance)
cfg.leadfield_type = 'spherical3layer'; % 'simple' or 'spherical3layer'
cfg.source_space   = 'cortical';    % 'cortical' or 'volumetric'
cfg.sensor_snr_db  = 15;            % target sensor SNR (trace-based)
cfg.random_seed    = 20250904;      % deterministic demo
cfg.edge_density   = 0.25;          % baseline sparsity
cfg.sparsity_variation = 0.30;      % dynamic changes across frequency

% Run Module 7
[Omega_true, Sigma_xx_true, Sigma_emp, P] = module7_simulation_improved_complex( ...
    'n_nodes', cfg.n_nodes, ...
    'n_sensors', cfg.n_sensors, ...
    'n_freq', cfg.n_freq, ...
    'n_samples', cfg.n_samples, ...
    'generate_leadfield', true, ...
    'leadfield_type', cfg.leadfield_type, ...
    'source_space', cfg.source_space, ...
    'sensor_snr_db', cfg.sensor_snr_db, ...
    'edge_density', cfg.edge_density, ...
    'sparsity_variation', cfg.sparsity_variation, ...
    'random_seed', cfg.random_seed);

% Short-hands
L   = P.leadfield_matrix;         % p x n
Ele = P.electrode_positions;      % p x 3
Src = P.source_positions;         % n x 3
p   = size(L,1); n = size(L,2); F = cfg.n_freq;
head_radius = get_head_radius(P, 0.10); % meters by default in our new API; adjust if you used cm

% True & observed sensor covariances (cells over frequency)
Svv_true = P.Sigma_vv_true;               % 1xF cell, p x p
Svv_obs  = P.Sigma_vv_observed;           % 1xF cell, p x p
% Per-frequency noise covariances (if available, else derive)
if isfield(P,'Sigma_xixi_cell') && numel(P.Sigma_xixi_cell)==F
    Snoise_cell = P.Sigma_xixi_cell;
else
    Snoise_cell = cell(1,F);
    for f=1:F, Snoise_cell{f} = (Svv_obs{f}-Svv_true{f} + (Svv_obs{f}-Svv_true{f})')/2; end
end

fprintf('Leadfield: [%d sensors x %d sources], cond = %.2e\n', p, n, cond(L));
fprintf('Freq bins: %d | Samples per freq: %d | SNR target: %g dB\n\n', F, cfg.n_samples, cfg.sensor_snr_db);

%% ------------------------- Figure 1: Geometry ----------------------------
nice_figure('Geometry: electrodes, sources, spheres');
tiledlayout(1,2, 'Padding','compact','TileSpacing','compact');

% Left: electrodes & head spheres
nexttile;
plot_head_spheres(head_radius);
hold on;
plot3(Ele(:,1), Ele(:,2), Ele(:,3), 'ko','MarkerFaceColor',[0.2 0.6 1], 'MarkerSize',6);
axis equal vis3d; grid on;
title(sprintf('Sensors (p=%d) on scalp (R=%.2f)', p, head_radius));
xlabel('x'); ylabel('y'); zlabel('z'); view(35,20);

% Right: sources on cortex/volume
nexttile;
plot_head_spheres(head_radius, 0.1); % faint shells
hold on;
plot3(Src(:,1), Src(:,2), Src(:,3), 'ro','MarkerFaceColor',[1 0.3 0.3], 'MarkerSize',5);
axis equal vis3d; grid on;
title(sprintf('Sources (n=%d) in %s space', n, cfg.source_space));
xlabel('x'); ylabel('y'); zlabel('z'); view(-135,20);

%% -------- Figure 2: Sensor covariances heatmaps (true vs observed) -------
nice_figure('Sensor covariances: true vs observed');
idx = pick_freqs(F); % pick [1, mid, end]
tiledlayout(2, numel(idx), 'Padding','compact','TileSpacing','compact');

for k = 1:numel(idx)
    f = idx(k);
    nexttile;
    imagesc(abs(Svv_true{f})); axis image; colorbar;
    title(sprintf('\\Sigma_{vv,true} @ f=%d', f)); set(gca,'YDir','normal');
    nexttile;
    imagesc(abs(Svv_obs{f})); axis image; colorbar;
    title(sprintf('\\Sigma_{vv,obs} @ f=%d', f)); set(gca,'YDir','normal');
end
colormap(turbo);

%% -------- Figure 3: Source connectome (from Ω at a chosen frequency) -----
nice_figure('Source connectome (top-|Ω| edges)');
f_show = idx(min(2,numel(idx))); % show the mid if possible
Kedges = min(150, n*(n-1)/2);    % draw up to K strongest edges
plot_source_connectome3D(Src, Omega_true{f_show}, Kedges, head_radius);
title(sprintf('Top |Ω| edges @ f=%d (n=%d)', f_show, n));

%% -------- Figure 4: Leadfield diagnostics -------------------------------
nice_figure('Leadfield diagnostics');
tiledlayout(2,2, 'Padding','compact','TileSpacing','compact');

% Column norms
nexttile;
cn = sqrt(sum(L.^2,1));
stem(1:n, cn, 'filled'); grid on;
xlabel('source index'); ylabel('||L(:,j)||_2');
title('Column norms');

% Singular values
nexttile;
s = svd(L);
semilogy(s,'o-','LineWidth',1.2); grid on;
xlabel('index'); ylabel('singular value');
title(sprintf('Singular values (cond=%.2e)', cond(L)));

% Heatmap of |L|
nexttile;
imagesc(abs(L)); axis image; colorbar; colormap(turbo);
xlabel('source'); ylabel('sensor');
title('|L| heatmap');

% Histogram of |L(:)|
nexttile;
histogram(abs(L(:)), 30, 'FaceColor',[0.2 0.6 1]); grid on;
xlabel('|L_{ij}|'); ylabel('count');
title('Distribution of |L| entries');

%% -------- Figure 5: Forward topography (3D colored electrodes) ----------
% Create a sparse source pattern (few active sources with complex weights)
active_k = max(3, round(0.05*n));
act_idx = randperm(n, active_k);
s_vec = zeros(n,1);
s_vec(act_idx) = (randn(active_k,1) + 1i*randn(active_k,1))/sqrt(2);

v = L * s_vec;            % sensor projection
vals = abs(v);            % magnitude for coloring

nice_figure('Forward topography (3D colored electrodes)');
scatter3(Ele(:,1), Ele(:,2), Ele(:,3), 150, vals, 'filled'); hold on;
plot_head_spheres(head_radius, 0.05);
axis equal vis3d; grid on; colorbar; colormap(turbo);
title(sprintf('Topography from %d active sources (|v|)', active_k));
xlabel('x'); ylabel('y'); zlabel('z'); view(40,20);

%% -------- Figure 6: Frequency traces (trace & SNR) ----------------------
true_trace = zeros(1,F);
noise_trace = zeros(1,F);
snr_db_ach = zeros(1,F);
for f=1:F
    % traces
    true_trace(f)  = real(trace(Svv_true{f}));
    noise_trace(f) = real(trace(Snoise_cell{f}));
    snr_db_ach(f)  = 10*log10( true_trace(f) / max(noise_trace(f),eps) );
end

nice_figure('Frequency traces (trace and SNR)');
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

nexttile;
plot(1:F, true_trace, '-o','LineWidth',1.4); hold on;
plot(1:F, noise_trace,'-s','LineWidth',1.2);
grid on; legend('trace(\Sigma_{vv,true})','trace(\Sigma_{noise})','Location','best');
xlabel('frequency index'); ylabel('trace');
title('Trace vs frequency');

nexttile;
plot(1:F, snr_db_ach, '-^','LineWidth',1.4); yline(cfg.sensor_snr_db,'--r','Target');
grid on; xlabel('frequency index'); ylabel('SNR (dB)');
title(sprintf('Achieved SNR (target = %g dB, trace-based)', cfg.sensor_snr_db));

fprintf('\n🎉 Visualization complete. Explore the figures!\n');

%% ========================= Helper Functions ==============================

function r = get_head_radius(P, default_val)
% Try to read head radius from params; else return default
    r = default_val;
    if isfield(P, 'forward_model_params') && isstruct(P.forward_model_params)
        if isfield(P.forward_model_params, 'head_radius')
            r = P.forward_model_params.head_radius;
            return;
        end
        % sometimes nested in layers info
        if isfield(P.forward_model_params,'sphere_layers') ...
                && isfield(P.forward_model_params.sphere_layers,'scalp_radius')
            r = P.forward_model_params.sphere_layers.scalp_radius;
            return;
        end
    end
    if isfield(P, 'head_radius'), r = P.head_radius; end
end

function fidx = pick_freqs(F)
% Pick representative frequency indices: [1, round(F/2), F]
    if F<=2, fidx = 1:F; return; end
    fidx = unique([1, round((F+1)/2), F]);
end

function plot_head_spheres(R, alpha_val)
% Draw semi-transparent scalp/brain spheres for context
% R is scalp radius. If alpha_val provided, spheres are faint.
    if nargin<2, alpha_val = 0.15; end
    [xs,ys,zs] = sphere(48);
    surf(R*xs, R*ys, R*zs, 'FaceColor',[0.85 0.85 0.85], ...
        'EdgeColor','none', 'FaceAlpha', alpha_val); hold on;
    % brain shell (slightly smaller)
    Rb = 0.87*R;
    surf(Rb*xs, Rb*ys, Rb*zs, 'FaceColor',[0.95 0.85 0.85], ...
        'EdgeColor','none', 'FaceAlpha', 0.07);
end

function plot_source_connectome3D(src_pos, Omega, K, head_radius)
% Plot top-K edges by |Omega_ij| on top of head spheres and source points
    n = size(src_pos,1);
    W = abs(Omega); W(1:n+1:end) = 0; % zero diagonal
    % upper triangle vector
    mask = triu(true(n),1);
    vals = W(mask);
    [~,ord] = sort(vals,'descend');
    K = min(K, numel(ord));
    [ii, jj] = find(mask);
    ii = ii(ord(1:K)); jj = jj(ord(1:K));

    plot_head_spheres(head_radius, 0.05); hold on;
    plot3(src_pos(:,1), src_pos(:,2), src_pos(:,3), ...
        'ro','MarkerFaceColor',[1 0.3 0.3],'MarkerSize',4);

    % draw edges
    cmap = turbo(256);
    vmax = max(vals(:))+eps;
    for k=1:K
        i = ii(k); j = jj(k);
        % color by weight magnitude
        cidx = max(1, round(255 * W(i,j)/vmax)+1);
        plot3([src_pos(i,1) src_pos(j,1)], ...
              [src_pos(i,2) src_pos(j,2)], ...
              [src_pos(i,3) src_pos(j,3)], ...
              '-', 'Color', cmap(cidx,:), 'LineWidth', 1.2);
        hold on;
    end
    axis equal vis3d; grid on; view(35,20);
    xlabel('x'); ylabel('y'); zlabel('z');
end

function nice_figure(tstr)
% Create a figure with good defaults
    f = figure('Color','w','Name',tstr);
    set(f, 'Units','normalized','Position',[0.08 0.08 0.84 0.78]);
end
