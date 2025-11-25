function Results_3D_realhead(output_sourse, sim_idx)
% RESULTS_3D_REALHEAD 3D brain connectivity and power maps for real-head simulations.
%
% This function assumes that you have already run:
%   Main_realistic_em_penalty_test_real_headmodel(output_sourse)
% which in turn calls:
%   simulation_real_headmodel  -> ./data/Pseudorand_Net.mat
%   InverseSolvers             -> ./result/Solutions_higgs.mat & Solutions_my_adapter.mat
%
% It then produces six 3D brain figures (if data available):
%   1) Simulated PCoh (ground-truth precision) as a 3D connectivity graph
%   2) Simulated source power spatial distribution on cortex
%   3) HIGGS-lasso PCoh 3D connectivity graph
%   4) HIGGS-lasso source power spatial distribution
%   5) My-adapter PCoh 3D connectivity graph
%   6) My-adapter source power spatial distribution
%
% USAGE
%   Results_3D_realhead;                 % use current folder, sim #1
%   Results_3D_realhead(out_dir, 2);     % use out_dir, simulation index 2
%
% NOTES
%   - All colormaps use "scientific" continuous maps (turbo/parula/hot) for
%     better readability.
%   - Power maps are normalized to [0,1] before plotting.
%
% Author: Augment Agent
% Date  : 2025-11-14

if nargin < 1 || isempty(output_sourse)
    output_sourse = pwd;
end
if nargin < 2 || isempty(sim_idx)
    sim_idx = 1;
end

%% Load simulation substrate (real-head)
net_path = fullfile('.', 'data', 'Pseudorand_Net.mat');
if exist(net_path, 'file') ~= 2
    error('Pseudorand_Net.mat not found at %s. Run simulation_real_headmodel first.', net_path);
end
S = load(net_path, 'J_sim', 'Seeders_sim', 'Theta_sim', 'vertices', 'faces', 'elec_pos', 'sens_system');
J_sim        = S.J_sim;
Seeders_sim  = S.Seeders_sim;
Theta_sim    = S.Theta_sim;
vertices     = S.vertices;
faces        = S.faces;
elec_pos     = S.elec_pos;

[Nseed, Nsim] = size(Seeders_sim);
if sim_idx < 1 || sim_idx > Nsim
    error('sim_idx out of range. Must be between 1 and %d.', Nsim);
end

seed_idx = Seeders_sim(:, sim_idx);
Theta_true = Theta_sim{sim_idx};
J_true     = J_sim{sim_idx};

%% Load HIGGS solutions (for HIGGS-lasso)
higgs_path = fullfile('.', 'result', 'Solutions_higgs.mat');
Theta_higgs = [];
if exist(higgs_path, 'file') == 2
    Sh = load(higgs_path, 'sol_higgs');
    if isfield(Sh, 'sol_higgs') && numel(Sh.sol_higgs) >= 3 && size(Sh.sol_higgs, 2) >= sim_idx
        % sol_higgs{3,sim}(:,:,1) : precision for HIGGS-lasso on simulation sim
        Theta_higgs = Sh.sol_higgs{3, sim_idx}(:,:,1);
    else
        warning('Solutions_higgs.mat does not contain expected sol_higgs{3,%d}(:,:,1).', sim_idx);
    end
else
    warning('Solutions_higgs.mat not found. HIGGS-lasso plots will be skipped.');
end

%% Load my-adapter solutions (external Rayleigh version)
my_path = fullfile('.', 'result', 'Solutions_my_adapter.mat');
Theta_my = [];
if exist(my_path, 'file') == 2
    Sm = load(my_path, 'sol_my_adapter');
    if isfield(Sm, 'sol_my_adapter') && isfield(Sm.sol_my_adapter, 'Theta') && ndims(Sm.sol_my_adapter.Theta) == 3
        if size(Sm.sol_my_adapter.Theta, 3) >= sim_idx
            Theta_my = Sm.sol_my_adapter.Theta_unb(:,:,sim_idx);
        else
            warning('sol_my_adapter.Theta does not contain simulation %d.', sim_idx);
        end
    else
        warning('Solutions_my_adapter.mat found but does not contain sol_my_adapter.Theta.');
    end
else
    warning('Solutions_my_adapter.mat not found. My-adapter plots will be skipped.');
end

%% Compute global connectivity scale (shared across methods)
q = numel(seed_idx);
global_max_conn = 0;

Theta_list = {Theta_true, Theta_higgs, Theta_my};
for kk = 1:numel(Theta_list)
    Th = Theta_list{kk};
    if ~isempty(Th)
        % 用真正的 |partial coherence| 幅度
        W_loc = theta_to_abs_pcoh(Th);
        mloc = max(W_loc(:));
        if isfinite(mloc) && mloc > global_max_conn
            global_max_conn = mloc;
        end
    end
end



%% Ensure output directory exists
if exist(output_sourse, 'dir') ~= 7
    mkdir(output_sourse);
end

%% 1) Simulated PCoh: 3D connectivity
fig1 = figure('Name', sprintf('Simulated PCoh 3D (sim %d)', sim_idx), ...
              'Position', [100 100 800 700]);
plot_pcoh_3d_brain(vertices, faces, seed_idx, Theta_true, ...
    sprintf('Simulated PCoh (ground truth, sim %d)', sim_idx), global_max_conn);
saveas(fig1, fullfile(output_sourse, sprintf('realhead_sim%d_pcoh_true.fig', sim_idx)));

%% 2) Simulated source power: 3D map
J_true = real(J_true(:));
if max(J_true) > 0
    J_true = J_true / max(J_true);
end
fig2 = figure('Name', sprintf('Simulated power 3D (sim %d)', sim_idx), ...
              'Position', [950 100 800 700]);
plot_power_3d_brain(vertices, faces, J_true, elec_pos, ...
    sprintf('Simulated source power (sim %d)', sim_idx));
saveas(fig2, fullfile(output_sourse, sprintf('realhead_sim%d_power_true.fig', sim_idx)));

%% 3–4) HIGGS-lasso PCoh & power (if available)
if ~isempty(Theta_higgs)
    fig3 = figure('Name', sprintf('HIGGS-lasso PCoh 3D (sim %d)', sim_idx), ...
                  'Position', [100 450 800 700]);
    plot_pcoh_3d_brain(vertices, faces, seed_idx, Theta_higgs, ...
        sprintf('HIGGS-lasso PCoh (sim %d)', sim_idx), global_max_conn);
    saveas(fig3, fullfile(output_sourse, sprintf('realhead_sim%d_pcoh_higgs_lasso.fig', sim_idx)));

    J_higgs = estimate_power_from_theta(Theta_higgs, vertices, seed_idx);
    fig4 = figure('Name', sprintf('HIGGS-lasso power 3D (sim %d)', sim_idx), ...
                  'Position', [950 450 800 700]);
    plot_power_3d_brain(vertices, faces, J_higgs, elec_pos, ...
        sprintf('HIGGS-lasso source power (sim %d)', sim_idx));
    saveas(fig4, fullfile(output_sourse, sprintf('realhead_sim%d_power_higgs_lasso.fig', sim_idx)));
end

%% 5–6) My-adapter PCoh & power (if available)
if ~isempty(Theta_my)
    fig5 = figure('Name', sprintf('My-adapter PCoh 3D (sim %d)', sim_idx), ...
                  'Position', [100 800 800 700]);
    plot_pcoh_3d_brain(vertices, faces, seed_idx, Theta_my, ...
        sprintf('My-adapter PCoh (sim %d)', sim_idx), global_max_conn);
    saveas(fig5, fullfile(output_sourse, sprintf('realhead_sim%d_pcoh_my_adapter.fig', sim_idx)));

    J_my = estimate_power_from_theta(Theta_my, vertices, seed_idx);
    fig6 = figure('Name', sprintf('My-adapter power 3D (sim %d)', sim_idx), ...
                  'Position', [950 800 800 700]);
    plot_power_3d_brain(vertices, faces, J_my, elec_pos, ...
        sprintf('My-adapter source power (sim %d)', sim_idx));
    saveas(fig6, fullfile(output_sourse, sprintf('realhead_sim%d_power_my_adapter.fig', sim_idx)));
end

end

%% =======================================================================
function plot_pcoh_3d_brain(vertices, faces, seed_idx, Theta, title_str, global_max_conn)
% Plot a 3D brain connectivity graph from a precision matrix (PCoh proxy).
%
% Inputs:
%   vertices        : [Nv x 3] cortical coordinates
%   faces           : [Nf x 3] triangulation
%   seed_idx        : [q x 1] indices of active generators (rows/cols of Theta)
%   Theta           : [q x q] precision (complex allowed)
%   title_str       : plot title
%   global_max_conn : optional global max |Theta| used for shared scaling

if nargin < 6
    global_max_conn = [];
end

Nv = size(vertices,1);
q  = numel(seed_idx);
if size(Theta,1) ~= q || size(Theta,2) ~= q
    error('Theta size (%dx%d) does not match number of seeds (%d).', ...
          size(Theta,1), size(Theta,2), q);
end

% Symmetrize and build partial-coherence-like weight matrix
W = theta_to_abs_pcoh(Theta);

% Use global scale if provided, otherwise fall back to local max
if nargin >= 6 && ~isempty(global_max_conn) && isfinite(global_max_conn) && global_max_conn > 0
    W = W / global_max_conn;
else
    maxW = max(W(:));
    if isfinite(maxW) && maxW > 0
        W = W / maxW;
    else
        W = zeros(size(W));   % no reliable connectivity information
    end
end

% Select strongest edges to avoid clutter
mask_triu = triu(true(q),1);
[i_idx, j_idx, w_vals] = find(W .* mask_triu);

% Keep only finite, strictly positive weights
valid = isfinite(w_vals) & (w_vals > 0);
i_idx = i_idx(valid);
j_idx = j_idx(valid);
w_vals = w_vals(valid);

if isempty(w_vals)
    warning('No valid off-diagonal connectivity found; skipping edges.');
end

[~, order] = sort(w_vals, 'descend');
max_edges = min(40, numel(order));
order = order(1:max_edges);
i_idx = i_idx(order);
j_idx = j_idx(order);
w_vals = w_vals(order);

% Nonlinear mapping to enhance weak edges visually
gamma = 0.5;                % you can tune this between ~0.4 and 0.6 if needed
w_vals_disp = w_vals .^ gamma;

% Node positions (in brain coordinates)
seed_pos = vertices(seed_idx, :);

% Draw cortex surface (light gray)
patch('Faces', faces, 'Vertices', vertices, ...
      'FaceVertexCData', ones(Nv,1), ...
      'FaceColor', [0.85 0.85 0.85], ...
      'EdgeColor', [0.3 0.3 0.3], ...
      'FaceAlpha', 0.15);
hold on;
axis vis3d;
axis equal off;
view(90, 0);  % lateral view; adjust as needed
camlight headlight; lighting gouraud;

% Plot generator nodes
scatter3(seed_pos(:,1), seed_pos(:,2), seed_pos(:,3), ...
         40, [0 0 0], 'filled');

% Edge colors using hot colormap (shared scale)
cmap  = hot(256);
alpha = 0.3;                     % lift dark end so weak edges are visible
cmap  = alpha + (1-alpha)*cmap;  % simple lightening
Nc    = size(cmap,1);
for e = 1:numel(w_vals)
    i = i_idx(e); j = j_idx(e);
    w = w_vals_disp(e);          % use gamma-corrected weight for display
    ci = max(1, min(Nc, round(1 + w*(Nc-1))));
    col = cmap(ci, :);
    lw  = 1.0 + 4.0*w;           % thicker lines, scaled by strength
    plot3([seed_pos(i,1) seed_pos(j,1)], ...
          [seed_pos(i,2) seed_pos(j,2)], ...
          [seed_pos(i,3) seed_pos(j,3)], ...
          '-', 'Color', col, 'LineWidth', lw);
end

% A simple color scale bar for edge strength (0-1)
colormap(cmap);
caxis([0 1]);
cb = colorbar();
set(cb, 'Ticks', [0 0.5 1], 'TickLabels', {'weak','medium','strong'});
cb.Label.String = '|partial coherence| (norm., \gamma = 0.5)';

title(title_str, 'Interpreter', 'none');
end

%% =======================================================================
function plot_power_3d_brain(vertices, faces, J, elec_pos, title_str)
% Plot a 3D brain surface with source power map J on vertices.

Nv = size(vertices,1);
if numel(J) ~= Nv
    error('Length of J (%d) does not match number of vertices (%d).', numel(J), Nv);
end

J = real(J(:));
if max(J) > 0
    J = J / max(J);
end

cmap = get_scientific_colormap();

patch('Faces', faces, 'Vertices', vertices, ...
      'FaceVertexCData', J, ...
      'FaceColor', 'interp', ...
      'EdgeColor', [0.3 0.3 0.3], ...
      'FaceAlpha', 0.95);
hold on;
axis vis3d;
axis equal off;
view(90, 0);
camlight headlight; lighting gouraud;

% Sensor positions (if provided)
if nargin >= 4 && ~isempty(elec_pos)
    if size(elec_pos,2) == 3
        scatter3(elec_pos(:,1), elec_pos(:,2), elec_pos(:,3), ...
                 10, [0 0 0], 'filled');
    else
        scatter(elec_pos(:,1), elec_pos(:,2), 10, 'k', 'filled');
    end
end

colormap(cmap);
colorbar;
caxis([0 1]);

title(title_str, 'Interpreter', 'none');
end

%% =======================================================================
function J_est = estimate_power_from_theta(Theta, vertices, seed_idx)
% Estimate a seed-level power map from a precision matrix via SPD inverse.
%
% We compute a stabilized inverse of Theta to obtain a covariance estimate
% Sigma, then take its diagonal as power at each generator, and finally
% place those values onto the cortical vertices at positions seed_idx.

Nv = size(vertices,1);
q  = numel(seed_idx);
if size(Theta,1) ~= q || size(Theta,2) ~= q
    error('Theta size (%dx%d) does not match number of seeds (%d).', ...
          size(Theta,1), size(Theta,2), q);
end

ThetaH = 0.5*(Theta + Theta');
ThetaH = full(ThetaH);

% Eigen-decomposition with small eigenvalue floor for stability
[U, d] = eig(ThetaH, 'vector');
d = real(d);
if isempty(d)
    J_est = zeros(Nv,1);
    return;
end

maxd = max(d);
if maxd <= 0
    d_safe = ones(size(d));
else
    d_safe = max(d, 1e-8 * max(maxd,1));
end

Sigma_est = U * diag(1./d_safe) * U';
Sigma_est = 0.5*(Sigma_est + Sigma_est');

p_seed = real(diag(Sigma_est));
p_seed = max(p_seed, 0);
if max(p_seed) > 0
    p_seed = p_seed / max(p_seed);
end

J_est = zeros(Nv,1);
J_est(seed_idx) = p_seed;
end

%% =======================================================================
function cmap = get_scientific_colormap()
% Get a perceptually reasonable continuous colormap.

try
    % MATLAB R2020b+ has turbo
    cmap = turbo(256);
catch
    try
        cmap = parula(256);
    catch
        cmap = turbo(256);
    end
end

end
%% =======================================================================
function W = theta_to_abs_pcoh(Theta)
% Convert a precision matrix Theta into |partial coherence| matrix.
%
% Theta : [q x q] complex precision matrix (should be Hermitian p.d.)
% W     : [q x q] matrix, W_ij = |pcoh_ij| in [0,1], with zero diagonal.

q = size(Theta,1);
if size(Theta,2) ~= q
    error('Theta must be square.');
end

% 强制 Hermitian（数值上对称）
ThetaH = 0.5*(Theta + Theta');   % 这里的 ' 是共轭转置

% 对角线（理论上应为正实数）
d = real(diag(ThetaH));

% 防止 0 或负数（数值抖动）
pos = d > 0;
if ~any(pos)
    % 极端情况下，退回单位阵
    d_safe = ones(size(d));
else
    dmin = max(min(d(pos)), 1e-8 * max(d(pos)));
    d_safe = max(d, dmin);
end

% sqrt(Theta_ii * Theta_jj) 的外积
D = sqrt(d_safe * d_safe.');

% 计算复偏相干
P = -ThetaH ./ D;
P(1:q+1:end) = 0;    % 对角线强制为 0

% 幅度
W = abs(P);
W(~isfinite(W)) = 0;
end

