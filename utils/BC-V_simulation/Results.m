function [] = Results(output_sourse)

%% Run Results
%% Loading Simulation Substrate or Real Data
load('colormap2.mat')
load('./data/Pseudorand_Net.mat')
measures_label    = {'auc' 'sens' 'spec' 'prec' 'f1'};
methods_label     = {'higgs-lasso';'higgs-ridge';'higgs-naive';'eloreta-hglasso';'lcmv-hglasso';'my-adapter'};
Nmeasures         = length(measures_label);
%%
%% h_hggm
load('./result/Solutions_higgs.mat')
[measures]        = quality_measures(sol_higgs,Theta_sim);
measures(isnan(measures)) = 0.5;
measures_mean     = round(100*mean(measures,3));
measures_std      = round(100*std(measures,0,3));
%%
%%
process_waitbar = waitbar(0,'Please wait...');
for meth = 1:5
    for meas = 1:5
        waitbar((meth*meas)/(25),process_waitbar,strcat('Outputing results....'));
        measures_higgs{meth,meas}  = [num2str(measures_mean(meth,meas)) '+/-' num2str(measures_std(meth,meas))];
    end
end
delete(process_waitbar);

%% Compute my-adapter metrics across simulations (final only, with internal Rayleigh)
my_measures_row = cell(1, Nmeasures);  % strings like 'mean+/-std'
try
    if exist('my_test2_adapter_bcvE_v2','file') == 2 && exist('compute_spd_and_classification_metrics.m','file') == 2
        Nsim_local = ceil(size(sol_higgs,2));
        scores_accum = nan(Nsim_local, Nmeasures);
        Lvj0 = LeadFields{1};
        for sim = 1:Nsim_local
            try
                Svv_sim_this = Svv_sim{1,sim}{1};
                seeders_this = Seeders_sim(:,sim);
                L_sel = Lvj0(:, seeders_this);
                T = Nsamp;
                [Omega_my_cell, ~, ~, ~] = my_test2_adapter_bcvE_v2({Svv_sim_this}, L_sel, T, struct('verbose', false, 'rayleigh', struct('enable', true)));
                if iscell(Omega_my_cell) && ~isempty(Omega_my_cell)
                    Theta_my_sim = 0.5*(Omega_my_cell{1} + Omega_my_cell{1}');
                    Theta_true   = Theta_sim{sim};
                    args = struct('Theta_true',Theta_true,'Theta_ray',Theta_my_sim);
                    opts = struct('normalize_mode','maxabs','spd',struct('symmetrize',true,'project',true,'eps',1e-10), ...
                                  'alpha',struct('value',0,'use_bcv_bug',false), 'plot', struct('radar', false));
                    eval_out = compute_spd_and_classification_metrics(args, opts);
                    m = eval_out.metrics.ray;
                    scores_accum(sim,:) = [m.auc, m.sens, m.spec, m.prec, m.f1];
                end
            catch
                % skip this sim
            end
        end
        mu = round(100*nanmean(scores_accum,1));
        sd = round(100*nanstd(scores_accum,0,1));
        for j=1:Nmeasures
            my_measures_row{j} = [num2str(mu(j)) '+/-' num2str(sd(j))];
        end
    else
        % Fallback: leave empty if function not available
        for j=1:Nmeasures, my_measures_row{j} = 'N/A'; end
    end
catch ME
    warning(ME.identifier, '%s', ME.message);
    for j=1:Nmeasures, my_measures_row{j} = 'N/A'; end
end

%% Table with quality measures (append my-adapter as the 6th method)
Nmethods = numel(methods_label);
Table = cell(Nmethods+1, Nmeasures+1);
Table(2:end,1)    = methods_label;
Table(1,2:end)    = measures_label;
Table(2:6,2:end)  = measures_higgs;   % 5 HIGGS methods
Table(7,2:end)    = my_measures_row;  % my-adapter

save(strcat(output_sourse,filesep,'Table_sens_system_',sens_system,'.mat'));
disp(strcat('Saving Table ---->  Table with quality measures to  ---> ', output_sourse) );

%%
%% Plot likelihood


Nsim = ceil(size(sol_higgs,2));
figure_likelihood = figure; 

subplot(3,2,1);
llh = zeros(length(sol_higgs{4,1}{1}{1}),Nsim);

for sim = 1:Nsim
  
    llh(:,sim) = sol_higgs{4,sim}{1}{1};
end

plot(llh);
ylabel('likelihood')
xlabel('iterations')
title('higgs-lasso likelihood')


%%
%%

subplot(3,2,2); 
llh = zeros(length(sol_higgs{4,1}{1}{2}),Nsim);
for sim = 1:Nsim
   
    llh(:,sim) = sol_higgs{4,sim}{1}{2};  
end

plot(llh);
ylabel('likelihood')
xlabel('iterations')
title('higgs-ridge likelihood')


%%

subplot(3,2,3); 
llh = zeros(length(sol_higgs{4,1}{1}{3}),Nsim);
for sim = 1:Nsim
    
    llh(:,sim) = sol_higgs{4,sim}{1}{3};
end

plot(llh);
ylabel('likelihood')
xlabel('iterations')
title('higgs-naive likelihood')


%%


subplot(3,2,4); 
for sim = 1:Nsim
    [gcv_opt,idx_gamma]       = min(sol_higgs{4,sim}{3});
    plot(sol_higgs{4,1}{2},sol_higgs{4,sim}{3},...
        '-',sol_higgs{4,1}{2}(idx_gamma),...
        gcv_opt,'b*');
    hold on;
end

ylabel('gcv value')
xlabel('regularization parameter')
title('eloreta-hglasso gcv function')

%% Plot corticaL map 
cortex.vertices = vertices;
cortex.faces    = faces;
[qL,qR,qfull,indvL,indvR,indv,verticesL,verticesR,vertices,facesL,facesR,faces,elec_pos_trans] = split_hemispheres(cortex,elec_pos);
J = zeros(qfull,1);
J(index_full) = 1;
%%
if strcmp(sens_system,'large') == 1 || strcmp(sens_system,'small') == 1
    subplot(3,2,5); 
    patch('Faces',facesL,'Vertices',verticesL,'FaceVertexCData',J(indvL),'FaceColor','interp',...
    'EdgeColor',[0.313725501298904 0.313725501298904 0.313725501298904],'FaceAlpha',.95);
    if strcmp(sens_system,'large') == 1
        axis([-0.1 0.1 -0.005 0.1 -0.1 0.1]); axis off; view([-179.6 24.8]); % large
    elseif strcmp(sens_system,'small') == 1
        axis([-90 5 -110 110 -65 95]); axis off; view([-88.8 20]); % small
    end
    hold on 
    scatter3(elec_pos_trans(:,1),elec_pos_trans(:,2),elec_pos_trans(:,3),'Marker','.','MarkerFaceColor','k','MarkerEdgeColor','k')
    colormap(cmap);
%     set(gcf,'Color','k');
    caxis([0 1]);
    title('realistic-head model','color','k')
    %%
    subplot(3,2,6);
    patch('Faces',facesR,'Vertices',verticesR,'FaceVertexCData',J(indvR),'FaceColor','interp',...
    'EdgeColor',[0.313725501298904 0.313725501298904 0.313725501298904],'FaceAlpha',.95);
    if strcmp(sens_system,'large') == 1
        axis([-0.1 0.1 -0.1 0.005 -0.1 0.1]); axis off; view([7.60000000000002 27.2]); % large
    elseif strcmp(sens_system,'small') == 1
        axis([-5 90 -110 110 -65 95]); axis off; view([85.2 26.4]); % small
    end
    hold on 
    scatter3(elec_pos_trans(:,1),elec_pos_trans(:,2),elec_pos_trans(:,3),'Marker','.','MarkerFaceColor','k','MarkerEdgeColor','k')
    colormap(cmap);
%     set(gcf,'Color','k');
    caxis([0 1]);
    title('realistic-head model','color','k')
else
    subplot(3,2,5); 
    patch('Faces',faces,'Vertices',vertices,'FaceVertexCData',J,'FaceColor','interp','Marker','o','MarkerFaceColor','y','EdgeColor',[0.313725501298904 0.313725501298904 0.313725501298904],'FaceAlpha',.95);
    hold on 
    scatter(elec_pos(:,1),elec_pos(:,2),'Marker','.','MarkerFaceColor','k','MarkerEdgeColor','k')
    axis off;
    colormap(cmap);
%     set(gcf,'Color','k');
    caxis([0 1]);
    title('pseudo-head model','color','k')
end

saveas( figure_likelihood,strcat(output_sourse,filesep,'higgs likelihood_sens_system_',sens_system,'.fig'));
disp(strcat('Saving figure ---->  higgs likelihood to  ---> ', output_sourse) );
delete(figure_likelihood);


%%
figure_partial_coherences = figure;
load('colormap2')

% Compute my-adapter (final, with internal Rayleigh) on sim#1 / subj#1
Theta_my = [];
try
    if exist('my_test2_adapter_bcvE_v2','file') == 2
        % Extract data consistent with sol_higgs first panel
        Svv_11 = Svv_sim{1,1}{1};
        Lvj0   = LeadFields{1};
        seeders_1 = Seeders_sim(:,1);
        L_sel  = Lvj0(:, seeders_1);
        T      = Nsamp;
        [Omega_my_cell, ~, ~, ~] = my_test2_adapter_bcvE_v2({Svv_11}, L_sel, T, struct('verbose', false, 'rayleigh', struct('enable', true)));
        if iscell(Omega_my_cell) && ~isempty(Omega_my_cell)
            Theta_my = 0.5*(Omega_my_cell{1} + Omega_my_cell{1}');
        end
    else
        warning('my_test2_adapter_bcvE_v2.m not found on path. Skipping my-adapter panel.');
    end
catch ME
    warning('my_test2_adapter_bcvE_v2 failed: %s', ME.message);
end

%% Plot partial correlations (expanded to 2x4 to include my-adapter)
X  = Theta_sim{1};
X  = X - diag(diag(X));
X  = X/max(abs(X(:))+eps);
subplot(2,4,1); imagesc(abs(X));
ylabel('generators')
xlabel('generators')
title('simulated PCoh')
%%
X  = sol_higgs{3,1}(:,:,1);
X  = X - diag(diag(X));
X  = X/max(abs(X(:))+eps);
subplot(2,4,2); imagesc(abs(X));
ylabel('generators')
xlabel('generators')
title('higgs-lasso PCoh')
%%
X  = sol_higgs{3,1}(:,:,2);
X  = X - diag(diag(X));
X  = X/max(abs(X(:))+eps);
subplot(2,4,3); imagesc(abs(X));
ylabel('generators')
xlabel('generators')
title('higgs-ridge PCoh')
%%
X  = sol_higgs{3,1}(:,:,3);
X  = X - diag(diag(X));
X  = X/max(abs(X(:))+eps);
subplot(2,4,4); imagesc(abs(X));
ylabel('generators')
xlabel('generators')
title('higgs-naive PCoh')
%%
X  = sol_higgs{3,1}(:,:,4);
X  = X - diag(diag(X));
X  = X/max(abs(X(:))+eps);
subplot(2,4,5); imagesc(abs(X));
ylabel('generators')
xlabel('generators')
title('eloreta-hglasso PCoh')
%%
X  = sol_higgs{3,1}(:,:,5);
X  = X - diag(diag(X));
X  = X/max(abs(X(:))+eps);
subplot(2,4,6); imagesc(abs(X));
ylabel('generators')
xlabel('generators')
title('lcmv-hglasso PCoh')
%%
if ~isempty(Theta_my)
    X  = Theta_my;
    X  = X - diag(diag(X));
    X  = X/max(abs(X(:))+eps);
    subplot(2,4,7); imagesc(abs(X));
    ylabel('generators')
    xlabel('generators')
    title('my-adapter (final) PCoh')
end
%%
colormap(cmap);

saveas( figure_partial_coherences,strcat(output_sourse,filesep,'partial_coherence_maps_sens_system_',sens_system,'.fig'));
disp(strcat('Saving figure ---->  Partial Coherence Maps to  ---> ', output_sourse) );
delete(figure_partial_coherences);
end