clear all
close all;
clc;

% Comparison of two sampling results

%% General option
% Path
addpath('/Users/janhasenauer/Documents/02_work/01_publications/01_manuscripts/01_ideas/2015 PLoS CB - Jagiella/Data');
addpath('/Users/janhasenauer/Documents/02_work/05_matlab/own_tools/PESTO');
filepath = '/Users/janhasenauer/Documents/02_work/01_publications/01_manuscripts/01_ideas/2015 PLoS CB - Jagiella/Data/';

% Text size
fs_label = 9;
ts_label = 7;

conf_level = [0.99,0.95,0.8,0];

% Printing options
resolution = '-r300';
color_true = [1,0.5,0.5];
lw_true = 5;


%% Datasets
% figure_name = 'fig__ExpDesign__ArtData';
% figure_name = 'fig__ExpDesign__ExpData';
figure_name = 'fig__2D_3D_comparison';

switch figure_name
    %
    case 'fig__ExpDesign__ArtData'
        % experimental data, 2D model, 1 conditions, spheroid radius, 100 samples
        filename_load{1} = 'TumorToyData2D_0.001adderr_100pop_GC'; 
        R_col{1}.n_generations = 31;
        R_col{1}.generations = R_col{1}.n_generations;
        % experimental data, 2D model, 1 conditions, spheroid radius / proliferating fraction, 100 samples
        filename_load{2} = 'TumorToyData2D_0.001adderr_100pop_GCKI67'; 
        R_col{2}.n_generations = 31;
        R_col{2}.generations = R_col{2}.n_generations;
        % experimental data, 2D model, 1 conditions, spheroid radius / ECM, 100 samples
        filename_load{3} = 'TumorToyData2D_0.001adderr_100pop_GCECM'; 
        R_col{3}.n_generations = 46;
        R_col{3}.generations = R_col{3}.n_generations;
        % experimental data, 2D model, 1 conditions, spheroid radius / proliferating fraction / ECM, 100 samples
        filename_load{4} = 'TumorToyData2D_0.001adderr_100pop_GCKI67ECM'; 
        R_col{4}.n_generations = 46;
        R_col{4}.generations = R_col{4}.n_generations;
        % Options
        plot_true = true;
        % Legend
        leg = {'estimates obtained using growth curves',...
               'estimates obtained using growth curves and fractions of proliferating cells',...
               'estimates obtained using growth curves and ECM intensities',...
               'estimates obtained using growth curves, fractions of proliferating cells and ECM intensities',...
               'true parameter estimates'};
        % Figure size & axis position
        f_size_comparison = [18,12];
        f_axis_comparison = [0.35 0.3 0.61 0.67];
        f_axis_comparison_legend = [0.12 0.02 0.75 0.15];
    %
    case 'fig__ExpDesign__ExpData'
        % experimental data, 2D model, 1 conditions, spheroid radius, 100 samples
        filename_load{1} = 'Tumor2dGC'; 
        R_col{1}.n_generations = 18; % max: 18
        R_col{1}.generations = R_col{1}.n_generations;
        % experimental data, 2D model, 1 conditions, spheroid radius / proliferating fraction, 100 samples
        filename_load{2} = 'Tumor2dGCKI67'; 
        R_col{2}.n_generations = 26; % max: 26
        R_col{2}.generations = R_col{2}.n_generations;
        % experimental data, 2D model, 1 conditions, spheroid radius / ECM, 100 samples
        filename_load{3} = 'Tumor2dGCECM'; 
        R_col{3}.n_generations = 30; % max: 30
        R_col{3}.generations = R_col{3}.n_generations;
        % experimental data, 2D model, 1 conditions, spheroid radius / proliferating fraction / ECM, 100 samples
        filename_load{4} = 'Tumor2dGCKI67ECM'; 
        R_col{4}.n_generations = 28; % max: 28
        R_col{4}.generations = R_col{4}.n_generations;
        % Options
        plot_true = false;
        % Legend
        leg = {'estimates obtained using growth curves',...
               'estimates obtained using growth curves and fractions of proliferating cells',...
               'estimates obtained using growth curves and ECM intensities',...
               'estimates obtained using growth curves, fractions of proliferating cells and ECM intensities - {\bf{reference solution}}'};
        % Figure size & axis position
        f_size_comparison = [18,9];
        f_axis_comparison = [0.35 0.30 0.61 0.67];
        f_axis_comparison_legend = [0.12 0.04 0.75 0.15];
    %
    case 'fig__2D_3D_comparison'
        % experimental data, 2D model, 1 conditions, spheroid radius, 100 samples
        filename_load{1} = 'Tumor2dGCKI67ECM'; 
        R_col{1}.n_generations = 26; % max: 28
        R_col{1}.generations = R_col{1}.n_generations;
        % experimental data, 2D model, 1 conditions, spheroid radius / proliferating fraction, 100 samples
        filename_load{2} = 'Tumor3dGCKI67ECM'; 
        R_col{2}.n_generations = 26; % max: 26
        R_col{2}.generations = R_col{2}.n_generations;
        % Options
        plot_true = false;
        % Legend
        leg = {'estimates obtained using 2D model',...
               'estimates obtained using 3D model'};
        % Figure size & axis position
        f_size_comparison = [18,7];
        f_axis_comparison = [0.35 0.28 0.61 0.70];
        f_axis_comparison_legend = [0.12 0.04 0.55 0.08];
end

%% Color-coding
switch figure_name
    %
    case {'fig__ExpDesign__ExpData',...
          'fig__ExpDesign__ArtData'}
        cm{1}(1,:) = 0.9*min([1.0000,0.8000,0.8000]   +0.0,1);
        cm{1}(3,:) = 0.9*min([1.0000,0.5020,0.5020]   +0.0,1);
        cm{1}(2,:) = 0.5*(cm{1}(1,:)+cm{1}(3,:));
        cm{1}(3,:) = cm{1}(2,:);
        cm{1}(2,:) = 0.5*(cm{1}(1,:)+cm{1}(3,:));

        cm{2}(1,:) = 0.9*min([1.0000,1.0000,0.8000]   +0.0,1);
        cm{2}(3,:) = 0.9*min([1.0000,1.0000,0.5020]   +0.0,1);
        cm{2}(2,:) = 0.5*(cm{2}(1,:)+cm{2}(3,:));
        cm{2}(3,:) = cm{2}(2,:);
        cm{2}(2,:) = 0.5*(cm{2}(1,:)+cm{2}(3,:));

        cm{3}(1,:) = 0.9*min([0.8000,1.0000,0.8000]   +0.0,1);
        cm{3}(3,:) = 0.9*min([0.5020,1.0000,0.5020]   +0.0,1);
        cm{3}(2,:) = 0.5*(cm{3}(1,:)+cm{3}(3,:));
        cm{3}(3,:) = cm{3}(2,:);
        cm{3}(2,:) = 0.5*(cm{3}(1,:)+cm{3}(3,:));

        cm{4}(1,:) = [0.5020,0.8000,1.0000];
        cm{4}(3,:) = [0.0000,0.5020,1.0000];
        cm{4}(2,:) = 0.5*(cm{4}(1,:)+cm{4}(3,:));
        cm{4}(3,:) = cm{4}(2,:);
        cm{4}(2,:) = 0.5*(cm{4}(1,:)+cm{4}(3,:));
        
    %
    case 'fig__2D_3D_comparison'
        cm{2}(1,:) = [255,204,102]/255;
        cm{2}(3,:) = [255,128,0]/255;
        cm{2}(2,:) = 0.5*(cm{2}(1,:)+cm{2}(3,:));
        cm{2}(3,:) = cm{2}(2,:);
        cm{2}(2,:) = 0.5*(cm{2}(1,:)+cm{2}(3,:));

        cm{1} = cm{2}(:,[3,2,1]);
end


%% Load datasets
for m = 1:length(filename_load)
    
R = R_col{m};

% Measurement data
warning off;
load([filename_load{m} '.mat']); 
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
    load([filepath filename_load{m} num2str(i) '.mat' ],'raw_data','x','f','neval' );
    
    % Re-organize simulation results
    if i == 1
        R.n_samples = size(x,2);
        R.n_theta = size(x,1);
        R.n_experiments = size(raw_data,2);
    end
    
    % Assignment of function evaluations
    R.n_acc(i) = R.n_samples/neval;
    R.n_eval(i) = neval;
    R.f(:,i) = f;
    
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
            R.exp{k}.x(:,j,i) = raw_data(J(j),k).x(:);
            R.exp{k}.y(:,j,i) = raw_data(J(j),k).y(:);
        end
    end
end

% Assignment
R_col{m} = R;

% Parameter collection
parameters.name = options.legend;
parameters.number = length(options.legend);
parameters.S.par = R_col{m}.theta(:,:,end);
parameters = getParameterConfidenceIntervals(parameters,1-conf_level);
parameters.min = options.lb;
parameters.max = options.ub;

parameters_col{m} = parameters;

end

%% Visualization of the parameter uncertainties
lh = [];

% Open figure
figure;
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 f_size_comparison],'renderer','zbuffer');

axes('Position',f_axis_comparison);
set(gca,'fontsize',ts_label,'box','on');

% Width of bars for different confidence levels
conf_level_w = 0.5*((1-conf_level)/max(1-conf_level)).^(1/3);

lb = min(options.lb);
ub = max(options.ub);

% Background
switch figure_name
    %
    case 'fig__ExpDesign__ExpData'
        % Rewrite files
        R_col_all = R_col; clear R_col
        parameters_col_all = parameters_col; clear parameters_col
        
        R_ref = R_col_all{4};
        parameters_ref = parameters_col_all{4};
        
        R_col{1} = R_col_all{1};
        R_col{2} = R_col_all{2};
        R_col{3} = R_col_all{3};
        parameters_col{1} = parameters_col_all{1};
        parameters_col{2} = parameters_col_all{2};
        parameters_col{3} = parameters_col_all{3};
        
        % Visualization
        for i = 1:parameters_ref.number
            for j = 1:(length(conf_level)-1)
                ind = find(cm{4}(j,:) ~= 1);
                c = cm{4}(j,:)+0*([1,1,1]-cm{4}(j,:));
                clh = fill(10.^parameters_ref.CI.S(i,[1,2,2,1,1],j),...
                     length(R_col)*(parameters_ref.number-i+1)-1+length(R_col)/2*[-1,-1,1,1,-1],'b',...
                     'facecolor',c,...
                     'edgecolor',c,'linewidth',0.1); hold on;
                 if (i == 1) && (j == 1)
                     lh(4) = clh;
                 end
            end
            plot(10.^parameters_ref.CI.S(i,[1,2],end),...
                length(R_col)*(parameters_ref.number-i+1)-1+length(R_col)/2*[-1,1],...
                '-','color',0.5*[1,1,1],'linewidth',2);
       end
end

% Loop: Parameters
for i = 1:parameters_col{1}.number
    % Plot true parameters
    if plot_true
    	clh = plot(10.^options.true([i,i]),...
            [length(R_col)*(parameters_col{m}.number-i+1)+0.5,...
             length(R_col)*(parameters_col{m}.number-i)+0.5],'-','linewidth',lw_true,'color',color_true); hold on;
        if (i == 1)
            lh(length(conf_level)+1) = clh;
        end
    end

    % Confidence intervals
    for j = 1:(length(conf_level)-1)
        for m = 1:length(R_col)
            clh = fill(10.^parameters_col{m}.CI.S(i,[1,2,2,1,1],j),...
                 length(R_col)*(parameters_col{m}.number+1-i)-m+1+conf_level_w(j)*[-1,-1,1,1,-1],'b',...
                 'facecolor',cm{m}(j,:),...
                 'edgecolor','k','linewidth',0.2); hold on;
            if (i == 1) && (j == (length(conf_level)-1))
                lh(m) = clh;
            end
        end
    end
    for m = 1:length(R_col)
        plot(10.^parameters_col{m}.CI.S(i,[1,2],end),length(R_col)*(parameters_col{m}.number+1-i)-m+1+[-0.35,+0.35],'k-');
    end
    
    % Seperating line
    plot(10.^[lb,ub],length(R_col)*(parameters_col{m}.number-i)+0.5*[1,1],'--','color',0.8*[1,1,1],'linewidth',1);
end

% Limits
parameter_label = {'ECM division threshold, e^{dic}',...
                   'ECM degradation rate, k_{e}^{deg}',...
                   'ECM production rate, k_{e}^{pro}',...
                   'initial quiescent cell fraction, q_{init}',...
                   'initial population radius, L^{init}',...
                   'division depth, L^{div}',...
                   'division rate, k_{div}^{max}',...
                   };
for i = 1:length(parameter_label)
    text(10.^(-5.3),length(R_col)*(i-0.5)+0.5,parameter_label{i},'HorizontalAlignment','right','fontsize',fs_label)
end
set(gca,'xlim',10.^[lb,ub],'xscale','log',...
        'ylim',[0.5,length(R_col)*parameters_col{m}.number+0.5],'ytick',[],...
        'fontsize',ts_label);

% Labels & legend
xlabel('parameters','fontsize',fs_label);
if length(leg) >= 1
    legend(lh,leg,'fontsize',fs_label,'box','off','position',f_axis_comparison_legend)
end

% Save figure
print('-depsc2',resolution,['./figures/' figure_name '.eps']);
    
