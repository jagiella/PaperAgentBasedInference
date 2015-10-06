clear all
close all;
clc;

%% General option
% Path
addpath('/Users/janhasenauer/Documents/02_work/01_publications/01_manuscripts/01_ideas/2015 PLoS CB - Jagiella/Data');
addpath('/Users/janhasenauer/Documents/02_work/05_matlab/Toolboxes/export_fig');
filepath = '/Users/janhasenauer/Documents/02_work/01_publications/01_manuscripts/01_ideas/2015 PLoS CB - Jagiella/Data/';

% Text size
fs_label = 12;
ts_label = 9;

% Figure size & axis position
f_size_lp = [7,6];
f_axis_lp = [0.2 0.2 0.75 0.75];

f_size_sp = [20,20];

data_col = [0,0.7,0];
data_lw = 1;
data_ms = 8;

%
m_sigma = 2;
conf_level = 0.90;
n_grid = 100;

%% Datasets
% filename = 'TumorToyData2D_0.001merr_100pop_GCKI67ECM';
% filename = 'Tumor2dfull'; % -> no results
% filename = 'Tumor2dfull_all_trunc_P500_'; % -> no results
% filename = 'Tumor2dGCECM';

%Tumor2dGCKI67ECM28
filename = 'Tumor2dfull_P1000_';

switch filename
    case 'TumorToyData2D_0.001merr_100pop_GCKI67ECM'
        % Artifical data
        R.generations_max = 40;
        R.generations = [4,7,12,20,40];
        R.n_generations = length(R.generations);

        R.exp{1}.data_type = 'spheroid radius';
        R.exp{1}.xlabel = 'time [d]';
        R.exp{1}.ylabel = 'spheroid radius [{\mu}m]';
        R.exp{1}.xlim = [0,20];
        R.exp{1}.xtick = [0:5:20];
        R.exp{1}.ylim = [0,768.6];
        R.exp{1}.ytick = [0:250:750];
        R.exp{1}.filename = 'spheroid_radius';

        R.exp{2}.data_type = 'proliferating cells';
        R.exp{2}.xlabel = 'distance to boundary [{\mu}m]';
        R.exp{2}.ylabel = 'proliferating cells [%]';
        R.exp{2}.xlim = [0,250];
        R.exp{2}.xtick = [0:50:250];
        R.exp{2}.ylim = [0,60];
        R.exp{2}.ytick = [0:20:60];
        R.exp{2}.filename = 'fraction_of_proliferating_cells';

        R.exp{3}.data_type = 'extracellular matrix intensity';
        R.exp{3}.xlabel = 'distance to boundary [{\mu}m]';
        R.exp{3}.ylabel = 'ECM intensity [au]';
        R.exp{3}.xlim = [0,250];
        R.exp{3}.xtick = [0:50:250];
        R.exp{3}.ylim = [0,0.8];
        R.exp{3}.ytick = [0:0.2:0.8];
        R.exp{3}.filename = 'extracellular_matrix_intensity';

    case {'Tumor2dfull','Tumor2dGCECM','Tumor2dfull_all_trunc_P500_',...
            'Tumor2dfull_P1000_','Tumor2dfull_trunc_fixerr_minmax_P100_'}
        % Artifical data
        R.generations_max = 25;
        R.generations = [2,5,8,14,25];
        switch filename
            case 'Tumor2dfull_all_trunc_P500_'
                R.generations_max = 18;
                R.generations = [2,5,8,12,18];
            case 'Tumor2dfull_trunc_fixerr_minmax_P100_'
                R.generations_max = 39;
                R.generations = [2,5,12,24,39];
        end
        R.n_generations = length(R.generations);

        R.exp{1}.data_type = 'spheroid radius';
        R.exp{1}.xlabel = 'time [d]';
        R.exp{1}.ylabel = 'spheroid radius [{\mu}m]';
        R.exp{1}.xlim = [0,20];
        R.exp{1}.xtick = [0:5:20];
        R.exp{1}.ylim = [0,768.6];
        R.exp{1}.ytick = [0:250:750];
        R.exp{1}.filename = 'spheroid_radius';

        R.exp{2}.data_type = 'proliferating cells';
        R.exp{2}.xlabel = 'distance to boundary [{\mu}m]';
        R.exp{2}.ylabel = 'proliferating cells [%]';
        R.exp{2}.xlim = [0,250];
        R.exp{2}.xtick = [0:50:250];
        R.exp{2}.ylim = [0,100];
        R.exp{2}.ytick = [0:25:100];
        R.exp{2}.filename = 'fraction_of_proliferating_cells';

        R.exp{3}.data_type = 'extracellular matrix intensity';
        R.exp{3}.xlabel = 'distance to boundary [{\mu}m]';
        R.exp{3}.ylabel = 'ECM intensity [au]';
        R.exp{3}.xlim = [0,250];
        R.exp{3}.xtick = [0:50:250];
        R.exp{3}.ylim = [0,0.2];
        R.exp{3}.ytick = [0:0.05:0.2];
        R.exp{3}.filename = 'extracellular_matrix_intensity';
end

% Number of iterations
R.generations_max = 0;
while exist([filepath filename num2str(R.generations_max+1) '.mat'], 'file')
    R.generations_max = R.generations_max+1;
end
%R.generations = round(linspace(5,R.generations_max,5));

%% Definition of colormap
fh = figure(999);
cm = colormap(gray(R.n_generations)); cm = 0.85*cm(end:-1:1,:); cm(end,:) = [1,0,0];
% cm = colormap(jet(R.n_generations));
close(fh)

%% Load the relevant datasets
% Measurement data
warning off;
load([filename '.mat']); 
warning on;
D = options.data;
if ~isfield(D,'raw')
    for k = 1:length(D)
        D(k).raw = D(k).y(:)';
    end
end

% Simulation
% Loop: Generations
for i = 1:R.n_generations
    % Load data
    load([filepath filename num2str(R.generations(i)) '.mat' ],'raw_data','x','f','neval' );
    
    % Re-organize simulation results
    if i == 1
        R.n_samples = size(x,2);
        R.n_theta = size(x,1);
        R.n_experiments = size(raw_data,2);
    end
    
    % Determine index
    J = [];
    for j = 1:size(raw_data,1)
        if ~isempty(raw_data(j,1).x)
            J(end+1) = j;
        end
    end
    
    % Loop: Samples
    for j = 1:R.n_samples
        % Initialization
        if (i == 1) && (j == 1)
            R.theta = nan(R.n_theta,R.n_samples,R.n_generations);
        end
        % Assignment of simulation results
        R.theta(:,j,i) = x(:,j);
        % Loop: Measurements
        for k = 1:R.n_experiments
            % Assignment of simulation results
            try
            R.exp{k}.x(:,j,i) = raw_data(J(j),k).x(:);
            R.exp{k}.y(:,j,i) = raw_data(J(j),k).y(:);
            end
        end
    end
end

%% Visualization of measurement data and simulation
% plot_options.lineplot = 'lines';
plot_options.lineplot = 'percentile intervals';

% Loop: Measurements
for k = 1:R.n_experiments
    % Open subplot
    figure
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 f_size_lp],'renderer','zbuffer');

    axes('Position',f_axis_lp);
    set(gca,'fontsize',ts_label,'xtick',R.exp{k}.xtick,'ytick',R.exp{k}.ytick,'box','on');
    
    % Scaling
    switch R.exp{k}.data_type
        case 'proliferating cells'
            scaling = 100;
            thinning = 20;
        case 'extracellular matrix intensity'
            scaling = 1;
            thinning = 20;
        otherwise
            scaling = 1;
            thinning = 1;
    end

    % Simulation
    % Loop: Generations
    for i = 1:R.n_generations
        switch plot_options.lineplot
            case 'lines'
                % Loop: Samples
                for j = 1:R.n_samples
                    plot(R.exp{k}.x(:,j,i),R.exp{k}.y(:,j,i),'-','color',cm(R.generations(i),:)); hold on
                end
            case 'percentile intervals'
                R.exp{k}.y_perc(:,:,i) = prctile(R.exp{k}.y(:,:,i),100*[(1-conf_level)/2,0.5,1-(1-conf_level)/2],2);
                ind = find(R.exp{k}.x(:,1,i) <= R.exp{k}.xlim(2));
                fill([R.exp{k}.x(ind,1,i);R.exp{k}.x(ind(end:-1:1),1,i)],...
                     scaling*[R.exp{k}.y_perc(ind,1,i);R.exp{k}.y_perc(ind(end:-1:1),3,i)],...
                     'b','facecolor',cm(i,:),'edgecolor',cm(i,:)); hold on;
        end

        % Boundary
        plot(R.exp{k}.xlim([1,2,2,1,1]),R.exp{k}.ylim([1,1,2,2,1]),'k-','linewidth',0.5);

        % Measurement data
        ind = 1:thinning:size(D(k).raw,2);
        for l = ind(:)'
        	plot(D(k).x(l),scaling*nanmean(D(k).raw(:,l)',2),'.','color',data_col,'markersize',data_ms);
        	plot(D(k).x(l)*[1,1],scaling*(nanmean(D(k).raw(:,l)',2)+[-1,1]*m_sigma*nanstd(D(k).raw(:,l))),...
                    'color',data_col,'linewidth',data_lw);
        end
i
%         % Limits
%         set(gca,'xlim',R.exp{k}.xlim,'xtick',R.exp{k}.xtick,...
%                 'ylim',R.exp{k}.ylim,'ytick',R.exp{k}.ytick,...
%                 'fontsize',ts_label);

        % Labels
        xlabel(R.exp{k}.xlabel,'fontsize',fs_label);
        ylabel(R.exp{k}.ylabel,'fontsize',fs_label);

        % Save figure
        print('-depsc2','-r1200',['./figure/' filename '__' R.exp{k}.filename '_' num2str(R.generations(i)) '.eps']);
    end
end

%% Visualization of parameter samples
% Open figure
figure;
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 f_size_sp],'renderer','zbuffer');

% Loop: Generations
for i = 1:R.n_generations
    % Scatter plot
    % Loop: parameter on y-axis
    for jy = 1:R.n_theta
    	% Loop: parameter on x-axis
        for jx = 1:(jy-1)
            % Open subplot
            subplottight(R.n_theta,R.n_theta,R.n_theta*(jy-1)+jx);

            % Lower and upper ound for grid
            lbx = min(R.theta(jx,:,i));
            ubx = max(R.theta(jx,:,i));
            lbx = max(lbx-0.2*(ubx-lbx),options.lb(jx));
            ubx = min(ubx+0.2*(ubx-lbx),options.ub(jx));

            lby = min(R.theta(jy,:,i));
            uby = max(R.theta(jy,:,i));
            lby = max(lby-0.2*(uby-lby),options.lb(jy));
            uby = min(uby+0.2*(uby-lby),options.ub(jy));

            % Density estimation
            [Theta1,Theta2] = meshgrid(linspace(lbx,ubx,n_grid),linspace(lby,uby,n_grid));
            KDest = kde_simple(R.theta([jx,jy],:,i),[Theta1(:),Theta2(:)]');
            KDest = KDest/sum(KDest);
            
            % Lower bound
            [KDest_sorted,I] = sort(KDest,'ascend');
            [C,h] = contourf(Theta1,Theta2,reshape(KDest,size(Theta1)),...
                             KDest_sorted(find(cumsum(KDest_sorted) > (1-conf_level),1,'first'))*[1,1],...
                             'facecolor',cm(i,:),'edgecolor',cm(i,:)); hold on;
            if i > R.n_generations-1
                plot(C(1,:),C(2,:),'.','color',cm(i,:),'markersize',2*i);
            end
            
            % Boundary
            plot([options.lb(jx),options.ub(jx),options.ub(jx),options.lb(jx),options.lb(jx)],...
                 [options.lb(jy),options.lb(jy),options.ub(jy),options.ub(jy),options.lb(jy)],...
                 'k-','linewidth',0.5);
             
            % Limits
            xlim([options.lb(jx),options.ub(jx)]);
            ylim([options.lb(jy),options.ub(jy)]);

            % Label
            if jx == 1
                ylabel(options.legend{jy},'fontsize',fs_label);
            else
                set(gca,'yticklabel',[]);
            end
            if jy == R.n_theta
                xlabel(options.legend{jx},'fontsize',fs_label);
            else
                set(gca,'xticklabel',[]);
            end
        end
    end
    
    % Denisty plot
    % Loop: parameter on y-axis
    for j = 1:R.n_theta
        % Open subplot
        subplottight(R.n_theta,R.n_theta,R.n_theta*(j-1)+j);

        % Lower and upper ound for grid
        lbx = min(R.theta(j,:,i));
        ubx = max(R.theta(j,:,i));
        lbx = max(lbx-0.2*(ubx-lbx),options.lb(j));
        ubx = min(ubx+0.2*(ubx-lbx),options.ub(j));

        % Density estimation
        theta = linspace(lbx,ubx,n_grid);
        KDest = kde_simple(R.theta(j,:,i),theta);
        if i == R.n_generations
            fill(theta([1,1:end,end]),[0,KDest,0]/max(KDest),...
                         'b','facecolor',cm(i,:),'edgecolor',cm(i,:)); hold on;
        else
            plot(theta(1:end),KDest/max(KDest),'-','color',cm(i,:),'linewidth',2); hold on;
        end
                 
        % Limits
        xlim([options.lb(j),options.ub(j)]);
        ylim([0,1.05]);

        % Label
        if j == 1
            ylabel(options.legend{j},'fontsize',fs_label);
        else
            set(gca,'yticklabel',[]);
        end
        if j ~= R.n_theta
            set(gca,'xticklabel',[]);
        end
    end
    
    % Save figure
    print('-depsc2','-r1200',['./figure/' filename '__scatterplot_' num2str(R.generations(i)) '.eps']);
end

%% Visualization of acceptance rate
% % Open figure
% figure;
% set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 f_size_lp],'renderer','zbuffer');
% 
% %


%function makeFigures()
    
%     close all;
%     visualiseComparison({'../Data/Tumor2dGC', '../Data/Tumor2dGCKI67', '../Data/Tumor2dGCECM', '../Data/Tumor2dGCKI67ECM'});
%     
%     close all;
%     visualiseComparison({'../Data/TumorToyData2D_0.001merr_200pop_GCKI67ECM','../Data/TumorToyData2D_0.001merr_100pop_GCKI67ECM','../Data/TumorToyData2D_0.001merr_10pop_GCKI67ECM'});
%     
%     close all;
%     visualiseABCResultsNEW( '../Data/Tumor3dGCKI67ECM');
%     
%     close all;
%     visualiseComparison({'../Data/Tumor2dGCKI67ECM','../Data/Tumor3dGCKI67ECM'});
% 
%     i=1;
%     M = [-1.4 2 1.1 -0.2 -2.3 -3.3 -2.3];
%     for p=1:7
%         m=M(p);
%         
%         for e=1:4
%             s=e;
%             
%             X(:,i) = normrnd(m,s, [100,1]);
%             L(1,i) = p;
%             G(1,i) = e;
%             i=i+1; 
%         end
%     end
%     
%     figure(1)
%     boxplot(X, {L,G})
%     
%     figure(2)
%     [m,s, idTab] = identifiabilityTab(X, {L,G})
%     errorbar(m,s)
%     
%     matrix2latex(idTab, unique(G), unique(L))
%end
