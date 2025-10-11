function plot_pcoh_panel(Omega_true, Omega_all, f_view)
% 画 f_view 频点的偏相干热图
if nargin<3, f_view = 1; end
pc = @(Om) pcorr_from_precision(Om);
Gt = as_cell_(Omega_true, f_view);

figure('Name',sprintf('PCOH compare (f=%d)', f_view), 'Position',[100 100 1200 350]);
names = fieldnames(Omega_all);
ncol  = numel(names) + 1;

subplot(1,ncol,1); imagesc(pc(Gt)); axis image; colorbar; title('GT PCOH'); colormap hot;
for k = 1:numel(names)
    Om = as_cell_(Omega_all.(names{k}), f_view);
    subplot(1,ncol,k+1); imagesc(pc(Om)); axis image; colorbar; title(names{k},'interpreter','none');
end
for ax = findall(gcf,'Type','axes')', caxis(ax,[0,0.45]); end

end

% === helpers ===
function P = pcorr_from_precision(Om)
Om = (Om + Om')/2;
d  = sqrt(max(real(diag(Om)), 1e-12));
D  = d*d.';
P  = abs(-Om ./ (D + eps));
P(1:size(Om,1)+1:end) = 0;
end

function A = as_cell_(X, f)
if iscell(X), A = X{min(f,numel(X))};
else, A = X(:,:,min(f,size(X,3))); end
end
