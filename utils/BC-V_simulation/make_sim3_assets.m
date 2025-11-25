function make_sim3_assets()
% 生成 small / large 两套资产到当前目录
gen_one('small',  30,  1500, 22);   % p=30,  n_src≈1.5k
gen_one('large', 30,  1500, 100);   % p=128, n_src≈4k
gen_one('extraLarge', 128,  8000, 100);   % p=128, n_src≈4k
gen_one('ultraLarge', 128,  4000, 4000);   % p=128, n_src≈4k

disp('Done. Generated: HeadModel_*.mat, LeadFields_*.mat, data_tips*.mat');
end

function gen_one(tag, p, nsrc, nseed)
r = [80 90 100];             % [brain, skull, scalp] 半径 (mm)
sig = [1 0.0125 1];          % [scalp, skull, brain] 电导率相对量
V  = fibonacci_sphere(nsrc, r(1));
F  = convhull(V(:,1), V(:,2), V(:,3));
cortex.vertices = V; cortex.faces = F;

E  = fibonacci_sphere(p, r(3));      % 电极在头皮上
coor = E;

% Leadfield：三层球近似（更“现实感”）；想要超快可换为 simple_leadfield(E,V)
L = spherical3layer_leadfield(E, V, r, sig, 0.0);

LeadFields = {L};
save(['HeadModel_' tag '.mat'],'cortex','coor');
save(['LeadFields_' tag '.mat'],'LeadFields');

% 远点采样做 22/24 个种子，等价手动 pickpoint
idx = fps_on_points(V, nseed);
for k=1:nseed, data_tips(k).Position = V(idx(k),:); end %#ok<AGROW>
if strcmp(tag,'small'), save('data_tips22_small.mat','data_tips');
elseif strcmp(tag,'large'),save('data_tips24_large.mat','data_tips'); 
elseif strcmp(tag,'extraLarge'), save('data_tips300_extraLarge', 'data_tips');
elseif strcmp(tag,'ultraLarge'), save('data_tips4000_ultraLarge', 'data_tips'); 
end
end

function P = fibonacci_sphere(n, rad)
i = (0:n-1)'; g = (1+sqrt(5))/2;
phi = acos(1 - 2*(i+0.5)/n); theta = 2*pi*(i+0.5)/g;
P = [rad*sin(phi).*cos(theta), rad*sin(phi).*sin(theta), rad*cos(phi)];
end

function idx = fps_on_points(P, k)
% farthest-point sampling（欧氏距离）在点集 P(n×3) 上选 k 个
n = size(P,1); idx = zeros(k,1); idx(1) = randi(n);
D = sum((P - P(idx(1),:)).^2,2);
for t=2:k
    [~,idx(t)] = max(D);
    D = min(D, sum((P - P(idx(t),:)).^2,2));
end
end

% ---- 物理近似 Leadfield（简单版）----
function L=simple_leadfield(E,S)
p=size(E,1); n=size(S,1); L=zeros(p,n);
for j=1:n, sj=S(j,:); rs=norm(sj)+eps;
  for i=1:p, ei=E(i,:); re=norm(ei)+eps;
    d=norm(ei-sj)+eps; cosang=dot(sj,ei)/(rs*re); orient=(1+cosang)/2;
    L(i,j)=orient/(d);
  end
end
L=L/(4*pi);
end

% ---- 三层球近似 Leadfield（默认使用）----
function L=spherical3layer_leadfield(E,S,rads,sig,jstd)
r_brain=rads(1); r_skull=rads(2); r_scalp=rads(3);
sigma_scalp=sig(1); sigma_skull=sig(2); sigma_brain=sig(3);
p=size(E,1); n=size(S,1); L=zeros(p,n); Ck=3;
for j=1:n
  sj=S(j,:); rs=norm(sj); if rs>=r_brain, sj=sj*(0.95*r_brain/rs); end
  for i=1:p
    ei=E(i,:); re=norm(ei); if abs(re-r_scalp)>1e-6, ei=ei*(r_scalp/re); end
    d=norm(ei-sj)+1e-9; cosang=dot(sj,ei)/(norm(sj)*norm(ei)); orient=(1+cosang)/2;
    base=1/(4*pi*sigma_scalp*d);
    dist_shape=1/(1+(d/r_scalp)^2);
    cond=(sigma_brain/sigma_scalp) * 1/(1 + Ck*(1/max(sigma_skull,1e-6)));
    pot=base*cond*dist_shape*orient; if jstd>0, pot=pot*(1+jstd*randn()); end
    L(i,j)=pot;
  end
end
end
