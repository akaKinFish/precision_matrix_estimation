function [active_mask, stats] = module3_active_set(InputMatrices, params)
% MODULE3_ACTIVE_SET - Structural Active Set Selection (Edge + Node).
%
% Purpose:
%   Generates a binary mask to restrict optimization to "active" elements.
%   Implements the "Node Active Set" hypothesis: Edges connected to 
%   inactive nodes should be pruned, enforcing structural sparsity.
%
% Usage:
%   1. Initialization: Pass Whitened Covariance (Sigma_tilde)
%   2. Dynamic Update: Pass Current Precision (Gamma)
%
% Logic:
%   1. Compute Edge Proxies C_ij = |M_ij|.
%   2. Compute Node Proxies R_i = max_j(C_ij) (exclude diag).
%   3. Determine Threshold tau (based on quantile).
%   4. Edge Mask: C_ij > tau.
%   5. Node Mask: R_i > tau.
%   6. Combined Mask: Edge_Mask AND Node_Mask(i) AND Node_Mask(j).
%
% Inputs:
%   InputMatrices : {F x 1} Cell array (Sigma or Gamma)
%   params        : Struct with fields:
%       .quantile_level : (double) Keep top q% edges (e.g., 0.10).
%                         If negative, interpreted as absolute threshold value.
%       .strategy       : 'intersection' (Edge & Node, default) or 'edge_only'.
%       .force_diagonal : (bool) Always keep diagonal active? (default: true)
%
% Output:
%   active_mask : {F x 1} Logical matrices.
%   stats       : Struct with sparsity details.

    % ============================================================
    % 1. Setup
    % ============================================================
    if ~iscell(InputMatrices), InputMatrices = {InputMatrices}; end
    F = numel(InputMatrices);
    p = size(InputMatrices{1}, 1);
    
    if nargin < 2, params = struct(); end
    
    % Default Parameters
    if isfield(params, 'quantile_level')
        q_level = params.quantile_level;
    else
        q_level = 0.10;
    end
    use_node_rule = ~isfield(params, 'strategy') || strcmp(params.strategy, 'intersection');
    keep_diag = ~isfield(params, 'force_diagonal') || params.force_diagonal;
    
    active_mask = cell(F, 1);
    
    % Stats containers
    stats.num_active_edges = zeros(F, 1);
    stats.num_active_nodes = zeros(F, 1);
    stats.thresholds = zeros(F, 1);
    
    % ============================================================
    % 2. Processing
    % ============================================================
    
    % We can compute threshold globally (across all F) or locally.
    % Let's do Per-Frequency for flexibility in brain imaging.
    
    for f = 1:F
        M = abs(InputMatrices{f});
        
        % Remove diagonal from proxy calculation (we only care about connectivity)
        M_off = M;
        M_off(1:p+1:end) = 0;
        
        % --- A. Threshold Determination ---
        if q_level < 1 && q_level > 0
            % Quantile based: Keep top q%
            % We look at the distribution of off-diagonal elements
            vals = M_off(tril(true(p), -1)); % Lower triangle values
            tau = quantile(vals, 1 - q_level);
        else
            % Absolute threshold provided (e.g., 1e-4 for zero cleaning)
            tau = abs(q_level);
        end
        stats.thresholds(f) = tau;
        
        % --- B. Edge Active Set ---
        % "Is this connection strong?"
        Mask_Edge = (M_off >= tau);
        
        % --- C. Node Active Set (Your Idea) ---
        % "Is this node participating in the network?"
        if use_node_rule
            % Node Proxy: Max connection strength of this node
            % R_i = max( |M_ij| ) for all j != i
            Node_Proxy = max(M_off, [], 2); 
            
            % Active Nodes
            Mask_Node_Vec = (Node_Proxy >= tau);
            
            % Expand to Matrix: Mask_Node(i,j) = Node(i) AND Node(j)
            % An edge exists only if BOTH endpoints are active nodes
            Mask_Node_Mat = Mask_Node_Vec & Mask_Node_Vec';
            
            % --- D. Combined Rule ---
            Final_Mask_Off = Mask_Edge & Mask_Node_Mat;
            stats.num_active_nodes(f) = sum(Mask_Node_Vec);
        else
            Final_Mask_Off = Mask_Edge;
            stats.num_active_nodes(f) = p;
        end
        
        % --- E. Final Assembly ---
        if keep_diag
            Final_Mask = Final_Mask_Off | eye(p);
        else
            Final_Mask = Final_Mask_Off;
        end
        
        active_mask{f} = logical(Final_Mask);
        
        % Count edges (undirected)
        stats.num_active_edges(f) = sum(sum(Final_Mask & ~eye(p))) / 2;
    end
    
    stats.density = mean(stats.num_active_edges) / (p*(p-1)/2);

end