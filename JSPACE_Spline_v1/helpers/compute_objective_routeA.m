function [obj, minEig_all, ok_spd] = compute_objective_routeA( ...
    Theta_cell, Psi_tilde_cell, B, Gdiff, w_edge, W_spatial, mcfg)
%COMPUTE_OBJECTIVE_ROUTEA  Calculate total objective function value.
%
% J(Theta) = DataTerm + Penalty_Group + Penalty_Pspline + Penalty_Spatial
% DataTerm = sum_t m_t * ( trace(Psi_t * Gamma_t) - logdet(Gamma_t) )
%
% Inputs:
%   Theta_cell     : Current spline coefficients
%   Psi_tilde_cell : Whitened sufficient stats
%   B, Gdiff       : Basis and P-spline difference matrices
%   w_edge         : Edge weights for Group Lasso
%   W_spatial      : Spatial weights for optional L2
%   mcfg           : Struct with .lambda1, .lambda_ps, .lambda2, .m_samples
%
% Outputs:
%   obj        : Scalar objective value (inf if not SPD)
%   minEig_all : Minimum eigenvalue across freq
%   ok_spd     : Boolean flag

    if nargin < 7 || isempty(mcfg)
        mcfg = struct();
    end

    eps_pd = 1e-6;
    if isfield(mcfg, 'spd') && isstruct(mcfg.spd) && isfield(mcfg.spd, 'eps_pd')
        if ~isempty(mcfg.spd.eps_pd)
            eps_pd = mcfg.spd.eps_pd;
        end
    elseif isfield(mcfg, 'eps_pd') && ~isempty(mcfg.eps_pd)
        eps_pd = mcfg.eps_pd;
    end

    % 1. Assemble Gamma matrices
    [Gamma_cell, ~, logdet_cell, minEig_all, ok_spd] = ...
        assemble_Gamma_from_Theta(Theta_cell, B, eps_pd);

    if ~ok_spd
        obj = inf;
        return;
    end

    F = numel(Gamma_cell);
    N = size(Gamma_cell{1}, 1);
    K = numel(Theta_cell);

    % Weights and coefficients
    lambda1 = get_field_default(mcfg, 'lambda1', 0);
    lambda_ps = get_field_default(mcfg, 'lambda_ps', 0);
    lambda2 = get_field_default(mcfg, 'lambda2', 0);
    m_samples = get_field_default(mcfg, 'm_samples', 1);

    % --- Term 1: Data Fitting (GGM Negative Log-Likelihood) ---
    obj_data = 0;
    for t = 1:F
        m_t = m_samples;
        if numel(m_samples) > 1
            m_t = m_samples(t);
        end

        tr_val = real(trace(Psi_tilde_cell{t} * Gamma_cell{t}));
        ld_val = logdet_cell{t};
        obj_data = obj_data + m_t * (tr_val - ld_val);
    end

    % --- Term 2: Group Lasso on Edges (Penalty) ---
    obj_group = 0;
    if lambda1 > 0
        if isempty(w_edge)
            w_edge = ones(N);
        end
        for i = 1:N
            for j = i+1:N
                v = zeros(K, 1);
                for k = 1:K
                    v(k) = Theta_cell{k}(i, j);
                end
                obj_group = obj_group + w_edge(i, j) * norm(v, 2);
            end
        end
        obj_group = lambda1 * obj_group;
    end

    % --- Term 3: P-spline Smoothing (Quadratic Penalty) ---
    % Include diagonal to keep objective consistent with gradients.
    obj_ps = 0;
    if lambda_ps > 0
        for i = 1:N
            for j = i:N
                v = zeros(K, 1);
                for k = 1:K
                    v(k) = Theta_cell{k}(i, j);
                end
                val = real(v' * Gdiff * v);
                obj_ps = obj_ps + val;
            end
        end
        obj_ps = lambda_ps * obj_ps;
    end

    % --- Term 4: Spatial L2 (Optional) ---
    obj_sp = 0;
    if lambda2 > 0 && ~isempty(W_spatial)
        W2 = W_spatial .^ 2;
        for k = 1:K
            term_k = sum(sum(W2 .* (abs(Theta_cell{k}) .^ 2)));
            obj_sp = obj_sp + term_k;
        end
        obj_sp = lambda2 * obj_sp;
    end

    % Total Objective
    obj = real(obj_data + obj_group + obj_ps + obj_sp);
end

function val = get_field_default(s, field_name, default_val)
    if isstruct(s) && isfield(s, field_name) && ~isempty(s.(field_name))
        val = s.(field_name);
    else
        val = default_val;
    end
end
