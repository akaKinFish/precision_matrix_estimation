function [true_precision, true_covariance, emp_covariance, params] = module7_simulation_improved_complex(varargin)
% MODULE7_SIMULATION_IMPROVED_COMPLEX - Complex Hermitian simulator with optional EEG leadfield
%
% (内容略，同上一版；为节省篇幅仅保留关键实现。若你需要我再贴一次完整头部注释，也可以。)

%% ---- Parse inputs ----
p = inputParser;
addParameter(p,'n_nodes',8,@(x)isscalar(x)&&x>2);
addParameter(p,'n_freq',12,@(x)isscalar(x)&&x>1);
addParameter(p,'n_samples',100,@(x)isscalar(x)&&x>0);
addParameter(p,'graph_type','random',@(x)any(strcmpi(x,{'random','chain','hub'})));
addParameter(p,'edge_density',0.4,@(x)isscalar(x)&&x>0&&x<1);
addParameter(p,'sparsity_variation',0.3,@(x)isscalar(x)&&x>=0&&x<=1);
addParameter(p,'edge_activation_smoothness',0.8,@(x)isscalar(x)&&x>0);
addParameter(p,'n_basis',4,@(x)isscalar(x)&&x>0);
addParameter(p,'sigma_coef',0.5,@(x)isscalar(x)&&x>0);
addParameter(p,'complex_strength',1.0,@(x)isscalar(x)&&x>=0&&x<=2);
addParameter(p,'epsilon_reg',0.1,@(x)isscalar(x)&&x>0);
addParameter(p,'random_seed',[],@(x)isempty(x)||(isscalar(x)&&x>=0));
addParameter(p,'coefficient_complex_fraction',1.0,@(x)isscalar(x)&&x>=0&&x<=1);

% leadfield
addParameter(p,'generate_leadfield',false,@(x)islogical(x));
addParameter(p,'n_sensors',[],@(x)isempty(x)||(isscalar(x)&&x>0));
addParameter(p,'leadfield_type','simple',@(x)any(strcmpi(x,{'simple','spherical3layer'})));

% spherical model geometry
addParameter(p,'head_radius',0.10,@(x)isscalar(x)&&x>0);
addParameter(p,'conductivities',[0.33,0.0042,0.33],@(x)isvector(x)&&numel(x)==3&&all(x>0));
addParameter(p,'layer_radii',[],@(x)isempty(x)||(isvector(x)&&numel(x)==3&&all(x>0)));

% electrodes & sources
addParameter(p,'electrode_layout','1020-19',@(x)any(strcmpi(x,{'1020-19','auto-uniform'})));
addParameter(p,'source_space','cortical',@(x)any(strcmpi(x,{'cortical','volumetric'})));
addParameter(p,'cortical_radius_ratio',0.82,@(x)isscalar(x)&&x>0&&x<1);

% realism knobs
addParameter(p,'source_perturbation_std',0,@(x)isscalar(x)&&x>=0);
addParameter(p,'leadfield_jitter_std',0,@(x)isscalar(x)&&x>=0);

% sensor noise
addParameter(p,'sensor_noise_mode','snr',@(x)any(strcmpi(x,{'snr','variance'})));
addParameter(p,'sensor_snr_db',20,@(x)isscalar(x));
addParameter(p,'sensor_noise_variance',[],@(x)isempty(x)||(isscalar(x)&&x>=0));

parse(p,varargin{:});
params = p.Results;

if params.generate_leadfield && isempty(params.n_sensors)
    params.n_sensors = 19;
end
if ~isempty(params.random_seed), rng(params.random_seed,'twister'); end

if params.generate_leadfield
    fprintf('=== Complex Simulation (leadfield true) ===\n');
else
    fprintf('=== Complex Simulation (leadfield false) ===\n');
end
fprintf('Sources (n): %d, Freq: %d, Samples per freq: %d\n', ...
    params.n_nodes, params.n_freq, params.n_samples);
if params.generate_leadfield
    fprintf('Sensors (p): %d (10-20 layout default) | leadfield_type=%s\n', ...
        params.n_sensors, params.leadfield_type);
end

%% Step 1: base graph
base_edge_list = generateGraphStructure(params.n_nodes, params.graph_type, params.edge_density);
n_base_edges = size(base_edge_list,1);
fprintf('Base edges: %d\n', n_base_edges);

%% Step 2: variable edges
all_edges_all=[]; for i=2:params.n_nodes, for j=1:i-1, all_edges_all=[all_edges_all; i,j]; end, end %#ok<AGROW>
base_set=[base_edge_list; base_edge_list(:,[2 1])];
is_base=false(size(all_edges_all,1),1);
for e=1:size(all_edges_all,1)
    is_base(e)=any(all(base_set==all_edges_all(e,:),2));
end
candidate_vars=all_edges_all(~is_base,:);
min_variable_edges=max(3,round(n_base_edges*0.4));
n_variable=round(params.sparsity_variation*n_base_edges)+min_variable_edges;
n_variable=min(n_variable,size(candidate_vars,1));
if n_variable>0
    idx=randperm(size(candidate_vars,1),n_variable);
    variable_edge_list=candidate_vars(idx,:);
else
    variable_edge_list=[];
end
n_variable_edges=size(variable_edge_list,1);
fprintf('Variable edges: %d\n', n_variable_edges);

%% Step 3: basis & coefficients
freq_norm=linspace(0,1,params.n_freq)';
basis_functions=generateSmoothBasis(freq_norm, params.n_basis);
fprintf('Coefficient complex fraction: %.2f\n', params.coefficient_complex_fraction);

n_total_edges = n_base_edges + n_variable_edges;
n_total_coeff = n_total_edges + params.n_nodes;
coefficients = zeros(n_total_coeff, params.n_basis);
for c=1:n_total_coeff
    for b=1:params.n_basis
        if rand() < params.coefficient_complex_fraction
            % 实部/虚部分别缩放，更直观
            coefficients(c,b) = randn()*params.sigma_coef + 1i*randn()*params.sigma_coef*params.complex_strength;
        else
            coefficients(c,b) = randn()*params.sigma_coef;
        end
    end
end

%% Step 4: variable edge activations  [修复 n_freq 很小时的 randperm 溢出]
variable_activations = ones(params.n_freq, n_variable_edges);
if n_variable_edges>0
    for e=1:n_variable_edges
        t = (0:params.n_freq-1)'/max(1,params.n_freq-1);
        slow   = sin(2*pi*(0.9+0.4*rand())*t + 2*pi*rand());
        medium = 0.7*sin(2*pi*(2.5+1.5*rand())*t + 2*pi*rand());
        fast   = 0.5*sin(2*pi*(6+6*rand())*t   + 2*pi*rand());
        finer  = 0.3*sin(2*pi*(15+10*rand())*t + 2*pi*rand());
        prob = 1./(1+exp(-params.edge_activation_smoothness*1.5*(slow+medium+fast+finer)));
        % neighbor smoothing
        alpha=0.3; prob2=prob;
        for f=2:params.n_freq-1
            prob2(f)=(1-alpha)*prob(f)+alpha*0.5*(prob(f-1)+prob(f+1));
        end
        thr = 0.4+0.3*rand();
        act = double(prob2 + 0.03*randn(size(prob2)) > thr + 0.2*sin(2*pi*4*t+2*pi*rand()));

        % ---- 强制最小变化，安全处理小 F ----
        want_changes = max(1, round(0.25*params.n_freq));      % 至少 1
        have_changes = sum(abs(diff(act)));
        if have_changes < want_changes
            nslots = max(1, params.n_freq-1);                  % 可翻转槽位数
            need   = min(nslots, want_changes - have_changes); % 保证 need ≤ nslots
            if need>0
                flip_pos = randperm(nslots, need);
                for kk = flip_pos
                    act(kk+1) = 1 - act(kk);
                end
            end
        end

        % 若仍全相同，再随即翻 1 位（n_freq 可能为 2 的极端情形）
        if all(act==act(1)) && params.n_freq>1
            kk = randi(params.n_freq-1);
            act(kk+1) = 1 - act(kk);
        end
        variable_activations(:,e)=act;
    end
end

%% Step 5: 构造 Ω_f 与 Σ_xx,f
true_precision  = cell(params.n_freq,1);
true_covariance = cell(params.n_freq,1);
sparsity_changes=0; prev_pat=[]; matrices_with_complex=0; hermitian_count=0;

for f=1:params.n_freq
    if mod(f-1, max(1,round(params.n_freq/6)))==0
        fprintf('Processing frequency %d/%d\n', f, params.n_freq);
    end
    Lc=zeros(params.n_nodes); edge_idx=1;
    % base
    for e=1:n_base_edges
        i=base_edge_list(e,1); j=base_edge_list(e,2);
        if i<j, tmp=i; i=j; j=tmp; end
        val = sum(coefficients(edge_idx,:).*basis_functions(f,:));
        Lc(i,j)=val; edge_idx=edge_idx+1;
    end
    % variable
    for e=1:n_variable_edges
        i=variable_edge_list(e,1); j=variable_edge_list(e,2);
        if i<j, tmp=i; i=j; j=tmp; end
        val = sum(coefficients(edge_idx,:).*basis_functions(f,:));
        val = val * variable_activations(f,e);
        if abs(val)>0.05*params.sigma_coef, Lc(i,j)=val; end
        edge_idx=edge_idx+1;
    end
    % diagonals
    for i=1:params.n_nodes
        dval = real(sum(coefficients(edge_idx,:).*basis_functions(f,:)));
        Lc(i,i)=abs(dval)+0.1; edge_idx=edge_idx+1;
    end
    % dominance
    for i=1:params.n_nodes
        off = sum(abs(Lc(i,1:i-1)))+sum(abs(Lc(i+1:end,i)));
        if abs(Lc(i,i)) < off + params.epsilon_reg
            if abs(Lc(i,i))>1e-12
                Lc(i,i)=(off+params.epsilon_reg)*(Lc(i,i)/abs(Lc(i,i)));
            else
                Lc(i,i)=off+params.epsilon_reg;
            end
        end
    end
    Omega = Lc*Lc'; Omega=(Omega+Omega')/2;
    for i=1:params.n_nodes
        Omega(i,i)=real(Omega(i,i));
        if Omega(i,i)<=0, Omega(i,i)=params.epsilon_reg; end
    end
    true_precision{f}=Omega;
    true_covariance{f}=inv_psd(Omega);

    if any(abs(imag(Omega(:)))>1e-12), matrices_with_complex=matrices_with_complex+1; end
    if ishermitian(Omega), hermitian_count=hermitian_count+1; end

    cur_pat=abs(triu(Omega,1))>1e-10;
    if ~isempty(prev_pat) && ~isequal(cur_pat,prev_pat), sparsity_changes=sparsity_changes+1; end
    prev_pat=cur_pat;
end

%% Step 6: Leadfield & sensor sampling
emp_covariance=cell(params.n_freq,1);
if params.generate_leadfield
    [L_sensors,Epos,Spos]=build_leadfield_isolated(params);

    Sigma_vv_true=cell(params.n_freq,1);
    Sigma_vv_obs =cell(params.n_freq,1);
    Sigma_xixi_cell=cell(params.n_freq,1);
    snr_lin=10^(params.sensor_snr_db/10);

    for f=1:params.n_freq
        Sxx=true_covariance{f};
        Svv=L_sensors*Sxx*L_sensors'; Svv=(Svv+Svv')/2; Sigma_vv_true{f}=Svv;
        if strcmpi(params.sensor_noise_mode,'variance')
            alpha = ifelse(isempty(params.sensor_noise_variance), 0, params.sensor_noise_variance);
        else
            alpha = trace(Svv)/max(1e-12, params.n_sensors*snr_lin);
        end
        Sxi = alpha*eye(params.n_sensors); Sigma_xixi_cell{f}=Sxi;
        Sobs=Svv+Sxi; Sobs=(Sobs+Sobs')/2; Sigma_vv_obs{f}=Sobs;
        emp_covariance{f}=sample_covariance(Sobs, params.n_samples);
    end

    params.leadfield_matrix=L_sensors;
    params.electrode_positions=Epos;
    params.source_positions=Spos;
    params.Sigma_vv_true=Sigma_vv_true;
    params.Sigma_vv_observed=Sigma_vv_obs;
    params.Sigma_xixi_cell=Sigma_xixi_cell;
    if ~isempty(Sigma_xixi_cell), params.Sigma_xixi=Sigma_xixi_cell{1}; end
else
    for f=1:params.n_freq
        emp_covariance{f}=sample_covariance(true_covariance{f}, params.n_samples);
    end
    params.leadfield_matrix=[]; params.electrode_positions=[]; params.source_positions=[];
end

%% Step 7: finalize
params.base_edge_list=base_edge_list;
params.variable_edge_list=variable_edge_list;
params.n_base_edges=n_base_edges;
params.n_variable_edges=n_variable_edges;
params.basis_functions=basis_functions;
params.coefficients=coefficients;
params.sparsity_changes=sparsity_changes;
params.matrices_with_complex=matrices_with_complex;
params.hermitian_count=hermitian_count;

fprintf('=== Done. leadfield=%d | Hermitian: %d/%d | Var edges: %d | Sparsity changes: %d ===\n', ...
    params.generate_leadfield, hermitian_count, params.n_freq, n_variable_edges, sparsity_changes);
end

%% ===== helpers =====
function y=ifelse(cond,a,b), if cond, y=a; else, y=b; end, end

function Ainv = inv_psd(A)
try
    Ainv = inv(A);
catch
    [V,D]=eig((A+A')/2);
    d=real(diag(D)); d=max(d,1e-8*max(d));
    Ainv=V*diag(1./d)*V'; Ainv=(Ainv+Ainv')/2;
end
end

function C = sample_covariance(Sigma, n_samples)
d=size(Sigma,1);
if any(abs(imag(Sigma(:)))>1e-12)
    L=chol_regularized(Sigma);
    Z=(randn(d,n_samples)+1i*randn(d,n_samples))/sqrt(2);
    X=L*Z; C=(X*X')/n_samples;
else
    L=chol_regularized(real(Sigma));
    Z=randn(d,n_samples); X=L*Z; C=(X*X')/n_samples;
end
C=(C+C')/2;
end

function L = chol_regularized(Sigma)
try
    L=chol((Sigma+Sigma')/2,'lower');
catch
    [V,D]=eig((Sigma+Sigma')/2);
    d=real(diag(D)); d=max(d,1e-10*max(d));
    L=V*diag(sqrt(d));
end
end

function edge_list = generateGraphStructure(n_nodes, graph_type, edge_density)
switch lower(graph_type)
    case 'random'
        n_possible=n_nodes*(n_nodes-1)/2;
        n_edges=round(edge_density*n_possible);
        all_edges=[]; for i=2:n_nodes, for j=1:i-1, all_edges=[all_edges; i,j]; end, end %#ok<AGROW>
        if n_edges>0 && n_edges<=size(all_edges,1)
            perm=randperm(size(all_edges,1)); edge_list=all_edges(perm(1:n_edges),:);
        else, edge_list=[]; end
    case 'chain'
        edge_list=[]; for i=2:n_nodes, edge_list=[edge_list; i,i-1]; end %#ok<AGROW>
        n_chain=size(edge_list,1);
        n_possible=n_nodes*(n_nodes-1)/2; n_target=round(edge_density*n_possible);
        n_add=max(0,n_target-n_chain);
        if n_add>0
            others=[]; for i=3:n_nodes, for j=1:i-2, others=[others; i,j]; end, end %#ok<AGROW>
            n_add=min(n_add,size(others,1));
            if n_add>0, perm=randperm(size(others,1)); edge_list=[edge_list; others(perm(1:n_add),:)]; end
        end
    case 'hub'
        hub=ceil(n_nodes/2); edge_list=[];
        for i=1:n_nodes
            if i==hub, continue; end
            if i>hub, edge_list=[edge_list; i,hub]; else, edge_list=[edge_list; hub,i]; end %#ok<AGROW>
        end
        n_hub=size(edge_list,1);
        n_possible=n_nodes*(n_nodes-1)/2; n_target=round(edge_density*n_possible);
        n_add=max(0,n_target-n_hub);
        if n_add>0
            others=[]; for i=2:n_nodes, for j=1:i-1
                if ~(i==hub || j==hub), others=[others; i,j]; end %#ok<AGROW>
            end, end
            n_add=min(n_add,size(others,1));
            if n_add>0, perm=randperm(size(others,1)); edge_list=[edge_list; others(perm(1:n_add),:)]; end
        end
    otherwise, error('Unknown graph_type: %s', graph_type);
end
if isempty(edge_list), edge_list=[2,1]; end
end

function B = generateSmoothBasis(freq_norm, n_basis)
B=zeros(numel(freq_norm),n_basis);
for b=1:n_basis
    if b==1, B(:,b)=1;
    elseif b==2, B(:,b)=freq_norm;
    else
        k=floor((b-2)/2)+1;
        if mod(b,2)==1, B(:,b)=sin(2*pi*k*freq_norm);
        else,            B(:,b)=cos(2*pi*k*freq_norm);
        end
    end
    nb=norm(B(:,b)); if nb>0, B(:,b)=B(:,b)/nb; end
end
end

%% -------- Leadfield (RNG isolation) --------
function [L,Epos,Spos]=build_leadfield_isolated(params)
s=rng; %#ok<NASGU>
base_seed=uint32(0); if ~isempty(params.random_seed), base_seed=uint32(params.random_seed); end
geo_key = uint32(double(base_seed) + 131*uint32(params.n_sensors) + 137*uint32(params.n_nodes) + ...
                  139*uint32(round(1e6*params.head_radius)) + 149*uint32(round(1e6*params.conductivities(2))));
rng(geo_key,'twister');

switch lower(params.electrode_layout)
    case '1020-19', Epos=electrodes_1020_19(params.n_sensors, params.head_radius);
    otherwise,      Epos=uniform_points_on_sphere(params.n_sensors, params.head_radius);
end

switch lower(params.source_space)
    case 'cortical'
        r_c=params.cortical_radius_ratio*params.head_radius;
        Spos=cortical_sources(params.n_nodes, r_c, params.source_perturbation_std);
    otherwise
        r1=0.3*params.head_radius; r2=0.87*params.head_radius;
        Spos=volumetric_sources(params.n_nodes, r1, r2);
end

switch lower(params.leadfield_type)
    case 'simple'
        L=simple_leadfield(Epos,Spos);
        cn=sqrt(sum(L.^2,1)); cn(cn==0)=1; L=L./cn;
    case 'spherical3layer'
        if isempty(params.layer_radii)
            r_brain=0.87*params.head_radius; r_skull=0.92*params.head_radius; r_scalp=1.00*params.head_radius;
        else
            r_brain=params.layer_radii(1); r_skull=params.layer_radii(2); r_scalp=params.layer_radii(3);
            if ~(r_brain<r_skull && r_skull<r_scalp), error('layer_radii must satisfy r_brain < r_skull < r_scalp'); end
        end
        sig=params.conductivities;
        L=spherical3layer_leadfield(Epos,Spos,[r_brain,r_skull,r_scalp], sig, params.leadfield_jitter_std);
    otherwise
        error('Invalid leadfield_type: %s', params.leadfield_type);
end

rng(s);
end

function E=electrodes_1020_19(p,R)
tpl=[0,90;-22.5,67.5;22.5,67.5;-90,22.5;-45,45;0,45;45,45;90,22.5;-135,0;-45,0;0,0;45,0;135,0;-90,-45;-45,-45;0,-45;45,-45;112.5,-45;-22.5,-67.5;22.5,-67.5];
if p>size(tpl,1)
    add=uniform_points_on_sphere(p-size(tpl,1),R);
    E=[sph2cart_list(tpl,R); add];
else
    E=sph2cart_list(tpl(1:p,:),R);
end
end

function P=sph2cart_list(tpl,R)
P=zeros(size(tpl,1),3);
for i=1:size(tpl,1)
    theta=deg2rad(tpl(i,1)); phi=deg2rad(tpl(i,2));
    pol=pi/2 - phi;
    P(i,1)=R*sin(pol)*cos(theta);
    P(i,2)=R*sin(pol)*sin(theta);
    P(i,3)=R*cos(pol);
end
end

function P=uniform_points_on_sphere(m,R)
P=randn(m,3); P=P./vecnorm(P,2,2); P=R*P;
for it=1:60
    F=zeros(size(P));
    for i=1:m
        for j=1:m
            if i==j, continue; end
            d=P(i,:)-P(j,:); nd=norm(d); if nd>0, F(i,:)=F(i,:)+d/(nd^3); end
        end
    end
    P=P+0.1*F; P=R*(P./max(vecnorm(P,2,2),1e-12));
end
end

function S=cortical_sources(n,Rc,pert_std)
S=zeros(n,3); idx=(0:n-1)'; th=acos(1-2*(idx+0.5)/n); ph=pi*(1+sqrt(5))*(idx+0.5);
S(:,1)=Rc.*sin(th).*cos(ph); S(:,2)=Rc.*sin(th).*sin(ph); S(:,3)=Rc.*cos(th);
if pert_std>0, S=S+(pert_std*Rc)*randn(size(S)); S=Rc*(S./max(vecnorm(S,2,2),1e-12)); end
end

function S=volumetric_sources(n,rmin,rmax)
S=zeros(n,3); k=1;
while k<=n
    v=2*rand(1,3)-1; nv=norm(v);
    if nv<=1 && nv>0
        r=rmin+(rmax-rmin)*rand(); S(k,:)=v*(r/nv); k=k+1;
    end
end
end

function L=simple_leadfield(E,S)
p=size(E,1); n=size(S,1); L=zeros(p,n);
for j=1:n
    sj=S(j,:); rs=norm(sj)+eps;
    for i=1:p
        ei=E(i,:); re=norm(ei)+eps;
        d=norm(ei-sj)+eps; cosang=dot(sj,ei)/(rs*re);
        orient=(1+cosang)/2;
        L(i,j)=orient/(d);
    end
end
L=L/(4*pi);
end

function L=spherical3layer_leadfield(E,S,rads,sig,jitter_std)
r_brain=rads(1); r_skull=rads(2); r_scalp=rads(3);
sigma_scalp=sig(1); sigma_skull=sig(2); sigma_brain=sig(3);
p=size(E,1); n=size(S,1); L=zeros(p,n);
Ck=3;
for j=1:n
    sj=S(j,:); rs=norm(sj);
    if rs>=r_brain, sj=sj*(0.95*r_brain/rs); rs=norm(sj); end
    for i=1:p
        ei=E(i,:); re=norm(ei);
        if abs(re-r_scalp)>1e-6, ei=ei*(r_scalp/re); re=r_scalp; end
        d=norm(ei-sj)+1e-9; cosang=dot(sj,ei)/(rs*re); orient=(1+cosang)/2;
        base=1/(4*pi*sigma_scalp*d);
        dist_shape=1/(1+(d/r_scalp)^2);
        cond_factor=(sigma_brain/sigma_scalp) * 1/(1 + Ck*(1/max(sigma_skull,1e-6)));
        pot=base*cond_factor*dist_shape*orient;
        if jitter_std>0, pot=pot*(1+jitter_std*randn()); end
        L(i,j)=pot;
    end
end
end
