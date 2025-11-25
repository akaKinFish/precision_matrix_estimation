function [Omega_true, Sigma_true, emp_cov, L, T, bayes_truth, sim] = module7_simulation_bayesian(cfg)
% MODULE7_SIMULATION_BAYESIAN  Fully-Bayesian simulator (generator) for multi-frequency sparse precisions
% 新增支持集跨频策略：
%   cfg.support_strategy = 'static' | 'threshold' | 'markov'  (default: 'threshold')
%   若 'threshold'：cfg.edge_density 用作每频非零率目标（或 cfg.edge_density_per_freq）
%   若 'markov'   ：cfg.persist ∈ [0,1] 控制“支持沿频率保持不变”的概率（默认 0.85）

% ---------- defaults ----------
p   = getd(cfg,'p',10);
m   = getd(cfg,'m',3);
F   = getd(cfg,'F',3);
T   = getd(cfg,'T',200);
rho = getd(cfg,'edge_density',0.15);                 % 目标稀疏度（默认用于阈值策略）
rhoF = getd(cfg,'edge_density_per_freq', []);        % 可选：长度F的每频目标稀疏度
modeN = getd(cfg,'modeN','diff');
support_strategy = lower(getd(cfg,'support_strategy','threshold'));
persist = getd(cfg,'persist', 0.85);                 % markov 的保持概率
sigma_xi2 = getd(cfg,'sigma_xi2',0.2);
delta_eps = getd(cfg,'delta_eps',[0.1,1e-6]);
[a1,b1] = dealvec(getd(cfg,'gamma_a1b1',[1e-3,1e-3]));
[a2,b2] = dealvec(getd(cfg,'gamma_a2b2',[1e-3,1e-3]));
wantComplex = getd(cfg,'complex_samples', true);
noise_type = getd(cfg,'noise_type','scalar');   % 'scalar' | 'matrix' | 'matrix_per_freq'
Sigma_xixi = [];
Sigma_xixi_per_freq = [];

% leadfield
L = getd(cfg,'L',[]);
if isempty(L), L = randn(m,p); end

switch lower(noise_type)
    case 'scalar'
        sigma_xi2 = getd(cfg,'sigma_xi2',0.2);
    case 'matrix'
        if isfield(cfg,'Sigma_xixi') && ~isempty(cfg.Sigma_xixi)
            Sx = (cfg.Sigma_xixi + cfg.Sigma_xixi')/2;
        else
            % 随机生成一个 SPD 噪声协方差（可按需改更贴近真实）
            G = randn(m); Sx = (G*G')/m + 0.1*eye(m);
        end
        % 轻度投影，确保 SPD
        [U,D] = eig(Sx); d = max(real(diag(D)), 1e-10);
        Sigma_xixi = U*diag(d)*U';
    case 'matrix_per_freq'
        if isfield(cfg,'Sigma_xixi_per_freq') && ~isempty(cfg.Sigma_xixi_per_freq)
            Sigma_xixi_per_freq = cfg.Sigma_xixi_per_freq; % cell(F,1), each m×m SPD
        else
            Sigma_xixi_per_freq = cell(F,1);
            for f=1:F
                G = randn(m); Sx = (G*G')/m + 0.1*eye(m);
                [U,D] = eig(Sx); d = max(real(diag(D)), 1e-10);
                Sigma_xixi_per_freq{f} = U*diag(d)*U';
            end
        end
    otherwise
        error('Unknown noise_type=%s', noise_type);
end

% Sigma_xixi = sigma_xi2 * eye(m);



% ---------- edge index & "eligible"（可候选边） ----------
[iu, ju] = upper_edge_index(p);  d = numel(iu);
if strcmpi(getd(cfg,'support_mode','erdos'),'mask')
    mask = cfg.support_mask;        % p×p logical upper
    if ~isequal(size(mask), [p,p])
        error('support_mask must be p-by-p logical.');
    end
    eligible = mask(sub2ind([p,p], iu, ju));   % 长度 d
else
    eligible = true(d,1);                      % 所有边皆可候选
end

% ---------- transforms M (L1) and N (L2) ----------
M = speye(d*F);
switch lower(modeN)
    case 'diff'
        Df = first_diff(F);                    % (F-1)×F
        N  = kron(Df, speye(d));              % q = (F-1)*d
    case 'kernel'
        K = getd(cfg,'K',[]);
        if isempty(K), error('modeN="kernel" requires cfg.K'); end
        K = (K+K')/2; [U,S] = eig(K); s = max(real(diag(S)), 1e-12);
        Kmh = U*diag(1./sqrt(s))*U.';         % K^{-1/2}
        N  = kron(Kmh, speye(d));
    otherwise
        error('Unknown modeN: %s', modeN);
end

% ---------- hyperparameters (fixed or sampled) ----------
if isfield(cfg,'lambda1_star') && ~isempty(cfg.lambda1_star)
    lambda1 = cfg.lambda1_star;
else
    lambda1_sq = gamrnd(a1, 1/b1); lambda1 = sqrt(lambda1_sq);
end
if isfield(cfg,'lambda2_star') && ~isempty(cfg.lambda2_star)
    lambda2 = cfg.lambda2_star;
else
    lambda2 = gamrnd(a2, 1/b2);
end

% ---------- 先按 L1+L2 高斯层级生成连续权重 beta（跨频已“平滑相关”） ----------
dF = d*F;
tau2 = exprnd(2/(lambda1^2), dF, 1);          % τ^2 ~ Exp(λ1^2/2)
Prec = spdiags(1./tau2, 0, dF, dF);
Q = M'*Prec*M + lambda2*(N'*N); Q = (Q+Q')/2;
R = chol(full(Q),'upper'); z = randn(dF,1);
beta = R \ z;                                 % β ~ N(0, Q^{-1})
BETA = reshape(beta, d, F);                   % d × F

% ---------- 让“支持随频率变化但平滑” ----------
switch support_strategy
    case 'static'
        % 单套支持（如需限制可候选边，用 eligible）
        S = repmat(eligible, 1, F);

    case 'threshold'
        % 每个频点单独按目标稀疏度阈值截断（在 eligible 内部做分位数）
        S = false(d, F);
        if isempty(rhoF), rhoF = rho * ones(1,F); end
        for f = 1:F
            idx = eligible;
            x = abs(BETA(idx, f));
            if ~any(idx), continue; end
            % 目标非零率 rhoF(f) → 选取分位数阈值
            k = max(1, round((1 - rhoF(f)) * nnz(idx)));
            xs = sort(x, 'ascend');
            tau_f = xs(min(k, numel(xs)));
            S(idx, f) = abs(BETA(idx,f)) >= tau_f;
        end

    case 'markov'
        % 为每条边建立沿频率的一阶马尔可夫链，持久度 persist
        base_rho = rho; if isempty(rhoF), rhoF = rho * ones(1,F); end
        S = false(d, F);
        % 初始
        S(:,1) = eligible & (rand(d,1) < rhoF(1));
        % 递推
        for f = 2:F
            stay = rand(d,1) < persist;
            newon = rand(d,1) < rhoF(f);
            Sf = S(:,f-1);
            Sf(~stay) = newon(~stay);
            S(:,f) = eligible & Sf;
        end

    otherwise
        error('Unknown support_strategy: %s', support_strategy);
end

% 应用支持（让支持不同频点不同；底层 BETA 已平滑 → 支持变化也会相对“顺滑”）
BETA = BETA .* S;

% ---------- 组装 Ω(f) 并做 SPD 对角补偿 ----------
[delta, epsv] = deal(delta_eps(1), delta_eps(2));
Omega_true = zeros(p,p,F);
for f = 1:F
    b = BETA(:,f);
    Bf = zeros(p,p);
    Bf(sub2ind([p,p], iu, ju)) = b;
    Bf = Bf + Bf.'; Bf(1:p+1:end) = 0;
    abs_row = sum(abs(Bf),2);
    Df = diag(delta + abs_row + epsv);
    Of = Df + Bf;
    Omega_true(:,:,f) = Of;
end

% ---------- 生成 Σ_true / emp_cov ----------
Sigma_true = cell(F,1); emp_cov = cell(F,1);
for f = 1:F
    Sjj = inv(Omega_true(:,:,f));
    switch lower(noise_type)
        case 'scalar'
            Svv = L*Sjj*L' + sigma_xi2 * eye(m);
        case 'matrix'
            Svv = L*Sjj*L' + Sigma_xixi;
        case 'matrix_per_freq'
            Svv = L*Sjj*L' + Sigma_xixi_per_freq{f};
    end
    Svv = (Svv + Svv')/2;  % 数值对称

    Svv = (Svv+Svv')/2;
    C = chol(Svv,'lower');
    if wantComplex
        Z = (randn(m,T) + 1i*randn(m,T))/sqrt(2);
    else
        Z = randn(m,T);
    end
    V = C * Z;
    emp = (V*V')/T; emp = (emp+emp')/2;
    Sigma_true{f} = Svv; emp_cov{f} = emp;
end

% ---------- truths & meta ----------
bayes_truth = struct('lambda1_true',lambda1, ...
    'lambda2_true',lambda2, ...
    'eligible', find(eligible), ...
    'support_strategy', support_strategy, ...
    'support_mask', S, ...
    'M', M, 'N', N, ...
    'delta', delta, 'eps', epsv, ...
    'Sigma_xixi', Sigma_xixi);
bayes_truth.noise_type = noise_type;
if strcmpi(noise_type,'scalar')
    bayes_truth.sigma_xi2_true = sigma_xi2;
elseif strcmpi(noise_type,'matrix')
    bayes_truth.Sigma_xixi_true = Sigma_xixi;
else
    bayes_truth.Sigma_xixi_per_freq_true = Sigma_xixi_per_freq;
end

sim = struct('p',p,'m',m,'F',F,'T',T, ...
    'edge_density',rho,'modeN',modeN, ...
    'timestamp', datestr(now));
end

% ===== helpers =====
function val = getd(s, name, def)
if isfield(s,name) && ~isempty(s.(name)), val = s.(name); else, val = def; end
end
function [a,b] = dealvec(v), a=v(1); b=v(2); end
function [iu, ju] = upper_edge_index(p)
K = nchoosek(1:p,2); iu=K(:,1); ju=K(:,2);
end
function D = first_diff(F)
D = zeros(F-1,F); for r=1:F-1, D(r,r)=-1; D(r,r+1)=1; end
end
