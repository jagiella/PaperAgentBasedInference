clear all
close all;
clc;

%% General option
% Path
addpath('/Users/janhasenauer/Documents/02_work/01_publications/01_manuscripts/01_ideas/2015 PLoS CB - Jagiella/Data');
addpath('/Users/janhasenauer/Documents/02_work/05_matlab/own_tools/PESTO');
filepath = '/Users/janhasenauer/Documents/02_work/01_publications/01_manuscripts/01_ideas/2015 PLoS CB - Jagiella/Data/';

% Text size
fs_label = 9;
ts_label = 7;

% Figure size & axis position
f_size_lp = [4.2,4];
f_axis_lp = [0.27 0.24 0.67 0.71];

f_size_wlp = [4.8,4];
f_axis_wlp = [0.3612 0.24 0.5863 0.71];

f_size_sp = [13,13];

f_size_acc = f_size_lp;
f_axis_acc = f_axis_lp;

data_col = [0,0.7,0];
data_lw = 0.75;
data_ms = 9;
data_errorbar_w = 1/2;

true_parameter_color = [255,255,0]/255;
%[255,102,102]/255; 
%[1 0.2 0.2];
true_parameter_lw = 0.5;
true_parameter_lt = '-';

%
m_sigma = 2;
conf_level = 0.90;
n_grid = 100;

margins_scatter_plot = [0.11 0.11 0.03 0.03];

% Printing options
resolution = '-r300';

%% Data types
E_sr.data_type = 'radius';
E_sr.xlabel = 'time [d]';
E_sr.ylabel = 'radius [{\mu}m]';
E_sr.filename = 'spheroid_radius';

E_pc.data_type = 'proliferating cells';
E_pc.xlabel = 'distance to rim [{\mu}m]';
E_pc.ylabel = 'proliferating cells [%]';
E_pc.filename = 'fraction_of_proliferating_cells';

E_emc.data_type = 'extracellular matrix intensity';
E_emc.xlabel = 'distance to rim [{\mu}m]';
E_emc.ylabel = 'ECM intensity [au]';
E_emc.filename = 'extracellular_matrix_intensity';

%% Datasets
Dataset_ind = 3;

switch Dataset_ind
    case 1
        % experimental data, 2D model, 1 condition, all data, 100 samples
        filename_load = 'Tumor2dGCKI67ECM'; 
        filename_save = 'ExpData_2D_1condition_sr_pf_ecm_100samples'; 
    case 2
        % experimental data, 3D model, 1 condition, all data, 100 samples
        filename_load = 'Tumor3dGCKI67ECM'; 
        filename_save = 'ExpData_3D_1condition_sr_pf_ecm_100samples'; 
    case 3
        % experimental data, 3D model, 4 conditions, all data, 100 samples
        filename_load = 'Tumor2dfull_trunc_fixerr_minmax_P100_'; 
        filename_save = 'ExpData_3D_4conditions_sr_pf_ecm_100samples'; 
    case 4
        % artifical data, 2D model, 1 conditions, all data, 100 samples
        filename_load = 'TumorToyData2D_0.001merr_100pop_GCKI67ECM'; 
        filename_save = 'ArtData_2D_1condition_sr_pf_ecm_100samples'; 
        
    % Experimental design - Artificial data
    case 40
        % experimental data, 2D model, 1 conditions, spheroid radius, 100 samples
        filename_load = 'TumorToyData2D_0.001adderr_100pop_GC'; % 37 
        filename_save = 'ArtData_2D_1condition_sr_100samples__ExpDesign'; 
    case 41
        % experimental data, 2D model, 1 conditions, spheroid radius / proliferating fraction, 100 samples
        filename_load = 'TumorToyData2D_0.001adderr_100pop_GCKI67'; % 32 
        filename_save = 'ArtData_2D_1condition_sr_pf_100samples__ExpDesign'; 
    case 42
        % experimental data, 2D model, 1 conditions, spheroid radius / ECM, 100 samples
        filename_load = 'TumorToyData2D_0.001adderr_100pop_GCECM'; % 47
        filename_save = 'ArtData_2D_1condition_sr_ecm_100samples__ExpDesign'; 
    case 43
        % experimental data, 2D model, 1 conditions, spheroid radius / proliferating fraction / ECM, 100 samples
        filename_load = 'TumorToyData2D_0.001adderr_100pop_GCKI67ECM'; % 48
        filename_save = 'ArtData_2D_1condition_sr_pf_ecm_100samples__ExpDesign'; 

    % Experimental design - Experimental data
    case 50
        % experimental data, 2D model, 1 conditions, spheroid radius, 100 samples
        filename_load = 'Tumor2dGC'; 
        filename_save = 'ExpData_2D_1condition_sr_100samples__ExpDesign'; 
    case 51
        % experimental data, 2D model, 1 conditions, spheroid radius / proliferating fraction, 100 samples
        filename_load = 'Tumor2dGCKI67'; 
        filename_save = 'ExpData_2D_1condition_sr_pf_100samples__ExpDesign'; 
    case 52
        % experimental data, 2D model, 1 conditions, spheroid radius / ECM, 100 samples
        filename_load = 'Tumor2dGCECM'; 
        filename_save = 'ExpData_2D_1condition_sr_ecm_100samples__ExpDesign'; 
    case 53
        % experimental data, 2D model, 1 conditions, spheroid radius / proliferating fraction / ECM, 100 samples
        filename_load = 'Tumor2dGCKI67ECM'; 
        filename_save = 'ExpData_2D_1condition_sr_pf_ecm_100samples__ExpDesign'; 
end
% filename = 'TumorToyData2D_0.001merr_100pop_GCKI67ECM';
% filename = 'Tumor2dfull'; % -> no results
% filename = 'Tumor2dfull_all_trunc_P500_'; % -> no results
% filename = 'Tumor2dGCECM';
% filename = 'Tumor2dfull_P1000_';
% filename = 'Tumor2dfull_trunc_fixerr_minmax_P100_';

switch filename_save        
    case {'ExpData_3D_4conditions_sr_pf_ecm_100samples'}
        % Generations
        R.n_generations = 33; % max: 33
        R.generations = round(linspace(3,R.n_generations,5));
        final_color = [255,128,0]/255;
        final_color_sp = [1.0000    0.5510    0.0500];

        % Options
        plot_true_parameters = false;

        % Limites & Ticks
        xtick_gen = 0:5:R.n_generations;

        ylim_acc = [0,0.75];
        ytick_acc = 0:0.25:0.75;
        
        ylim_dist = [5e-1,1e2];
        ytick_dist = 10.^[0:2];

        ylim_fval = [1e4/2,2e7];
        ytick_fval = 10.^[4:7];
        
        E_sr.xlim = [0,20];
        E_sr.xtick = [0:5:20];
        E_sr.ylim = [0,500];
        E_sr.ytick = [0:100:500];
        E_sr.plot_data = true;

        E_pc.xlim = [0,250];
        E_pc.xtick = [0:50:250];
        E_pc.ylim = [0,100];
        E_pc.ytick = [0:25:100];
        E_pc.plot_data = true;
      
        E_emc.xlim = [0,250];
        E_emc.xtick = [0:50:250];
        E_emc.ylim = [0,1];
        E_emc.ytick = [0:0.25:1];
        E_emc.plot_data = true;
        
        % Datasets
        R.exp{1} = E_sr;
        R.exp{2} = E_sr;
        R.exp{3} = E_sr;
        R.exp{4} = E_sr;
        R.exp{5} = E_pc;
        R.exp{6} = E_emc;
        R.exp{7} = E_pc;
        R.exp{8} = E_emc;
        R.exp{9} = E_pc;
        R.exp{10} = E_emc;
        R.exp{11} = E_pc;
        R.exp{12} = E_emc;
        R.exp{13} = E_pc;
        R.exp{14} = E_emc;
        R.exp{15} = E_pc;
        R.exp{16} = E_emc;
    
    case {'ExpData_2D_1condition_sr_pf_ecm_100samples',...
          'ExpData_3D_1condition_sr_pf_ecm_100samples'}
        % Generations
        switch filename_save
            case 'ExpData_2D_1condition_sr_pf_ecm_100samples'
                R.n_generations = 26; % max: 28
                final_color = [0,128,255]/255;
                final_color_sp = [0.0500    0.5510    1.0000];
            case 'ExpData_3D_1condition_sr_pf_ecm_100samples'
                R.n_generations = 26; % max: 26
                final_color = [255,128,0]/255;
                final_color_sp = [1.0000    0.5510    0.0500];
        end
        R.generations = round(linspace(3,R.n_generations,5));

        % Options
        plot_true_parameters = false;

        % Limites & Ticks
        xtick_gen = 0:5:R.n_generations;

        ylim_acc = [0,0.5];
        ytick_acc = 0:0.1:0.5;
        
        ylim_dist = [1e0,1e3];
        ytick_dist = 10.^[0:3];

        ylim_fval = [1e2/2,1e5];
        ytick_fval = 10.^[2:5];

        E_sr.xlim = [0,20];
        E_sr.xtick = [0:5:20];
        E_sr.ylim = [0,500];
        E_sr.ytick = [0:100:500];
        E_sr.plot_data = true;

        E_pc.xlim = [0,250];
        E_pc.xtick = [0:50:250];
        E_pc.ylim = [0,100];
        E_pc.ytick = [0:25:100];
        E_pc.plot_data = true;
      
        E_emc.xlim = [0,250];
        E_emc.xtick = [0:50:250];
        E_emc.ylim = [0,0.2];
        E_emc.ytick = [0:0.05:0.2];
        E_emc.plot_data = true;
        
        % Datasets
        R.exp{1} = E_sr;
        R.exp{2} = E_pc;
        R.exp{3} = E_emc;
        
    case {'ArtData_2D_1condition_sr_pf_ecm_100samples'}
        % Generations
        switch filename_save
            case 'ArtData_2D_1condition_sr_pf_ecm_100samples'
                R.n_generations = 41; % max: 45
                final_color    = [0         0.5020    1.0000];
                final_color_sp = [0.0500    0.5510    1.0000];
        end
        R.generations = [2,8,16,26,41];
        
        % Options
        plot_true_parameters = true;
        
        % Limites & Ticks
        xtick_gen = 0:10:R.n_generations;
        
        ylim_acc = [0,0.6];
        ytick_acc = 0:0.1:0.6;

        ylim_dist = [1e-1,1e5];
        ytick_dist = 10.^[-1:5];
        
        ylim_fval = [1e2,1e6];
        ytick_fval = 10.^[2:6];

        E_sr.xlim = [0,20];
        E_sr.xtick = [0:5:20];
        E_sr.ylim = [0,776];
        E_sr.ytick = [0:250:776];
        E_sr.plot_data = true;

        E_pc.xlim = [0,300];
        E_pc.xtick = [0:100:300];
        E_pc.ylim = [0,100];
        E_pc.ytick = [0:25:100];
        E_pc.plot_data = true;
      
        E_emc.xlim = [0,300];
        E_emc.xtick = [0:100:300];
        E_emc.ylim = [0,0.7];
        E_emc.ytick = [0:0.1:0.7];
        E_emc.plot_data = true;
        
        % Datasets
        R.exp{1} = E_sr;
        R.exp{2} = E_pc;
        R.exp{3} = E_emc;
        
    case {'ArtData_2D_1condition_sr_100samples__ExpDesign',...
          'ArtData_2D_1condition_sr_pf_100samples__ExpDesign',...
          'ArtData_2D_1condition_sr_ecm_100samples__ExpDesign',...
          'ArtData_2D_1condition_sr_pf_ecm_100samples__ExpDesign'}
        % Generations
        switch filename_save
            case 'ArtData_2D_1condition_sr_100samples__ExpDesign'
                R.n_generations = 31; % max: 37
                final_color    = [1.0000    0.5020    0.5020];
                final_color_sp = [1.0000    0.5510    0.5510];
                final_color    = 0.9*min(final_color   +0.0,1);
                final_color_sp = 0.9*min(final_color_sp+0.0,1);
            case 'ArtData_2D_1condition_sr_pf_100samples__ExpDesign'
                R.n_generations = 31; % max: 32
                final_color    = [1.0000    1.0000    0.5020];
                final_color_sp = [1.0000    1.0000    0.5510];
                final_color    = 0.9*min(final_color   +0.0,1);
                final_color_sp = 0.9*min(final_color_sp+0.0,1);
            case 'ArtData_2D_1condition_sr_ecm_100samples__ExpDesign'
                R.n_generations = 46; % max: 47
                final_color    = [0.5020    1.0000    0.5020];
                final_color_sp = [0.5510    1.0000    0.5510];
                final_color    = 0.9*min(final_color   +0.0,1);
                final_color_sp = 0.9*min(final_color_sp+0.0,1);
            case 'ArtData_2D_1condition_sr_pf_ecm_100samples__ExpDesign'
                R.n_generations = 46; % max: 48
                final_color    = [0         0.5020    1.0000];
                final_color_sp = [0.0500    0.5510    1.0000];
        end
        R.generations = round(linspace(3,R.n_generations,5));
        
        % Options
        plot_true_parameters = false;
        
        % Limites & Ticks
        xtick_gen = 0:10:R.n_generations;
        
        ylim_acc = [0,0.6];
        ytick_acc = 0:0.1:0.6;

        ylim_dist = [1e-2,1e3];
        ytick_dist = 10.^[-2:3];
        
        ylim_fval = [1e2,1e6];
        ytick_fval = 10.^[2:6];

        E_sr.xlim = [0,20];
        E_sr.xtick = [0:5:20];
        E_sr.ylim = [0,776];
        E_sr.ytick = [0:250:776];

        E_pc.xlim = [0,300];
        E_pc.xtick = [0:100:300];
        E_pc.ylim = [0,100];
        E_pc.ytick = [0:25:100];
      
        E_emc.xlim = [0,300];
        E_emc.xtick = [0:100:300];
        E_emc.ylim = [0,0.7];
        E_emc.ytick = [0:0.1:0.7];
        
        switch filename_save
            case 'ArtData_2D_1condition_sr_100samples__ExpDesign'
                E_sr.plot_data = true;
                E_pc.plot_data = false;
                E_emc.plot_data = false;
            case 'ArtData_2D_1condition_sr_pf_100samples__ExpDesign'
                E_sr.plot_data = true;
                E_pc.plot_data = true;
                E_emc.plot_data = false;
            case 'ArtData_2D_1condition_sr_ecm_100samples__ExpDesign'
                E_sr.plot_data = true;
                E_pc.plot_data = false;
                E_emc.plot_data = true;
            case 'ArtData_2D_1condition_sr_pf_ecm_100samples__ExpDesign'
                E_sr.plot_data = true;
                E_pc.plot_data = true;
                E_emc.plot_data = true;
        end
        
        % Datasets
        R.exp{1} = E_sr;
        R.exp{2} = E_pc;
        R.exp{3} = E_emc;

    case {'ExpData_2D_1condition_sr_100samples__ExpDesign',...
          'ExpData_2D_1condition_sr_pf_100samples__ExpDesign',...
          'ExpData_2D_1condition_sr_ecm_100samples__ExpDesign',...
          'ExpData_2D_1condition_sr_pf_ecm_100samples__ExpDesign'}
        % Generations
        switch filename_save
            case 'ExpData_2D_1condition_sr_100samples__ExpDesign'
                R.n_generations = 18; % max: 18
                final_color    = [1.0000    0.5020    0.5020];
                final_color_sp = [1.0000    0.5510    0.5510];
                final_color    = 0.9*min(final_color   +0.0,1);
                final_color_sp = 0.9*min(final_color_sp+0.0,1);
            case 'ExpData_2D_1condition_sr_pf_100samples__ExpDesign'
                R.n_generations = 26; % max: 26
                final_color    = [1.0000    1.0000    0.5020];
                final_color_sp = [1.0000    1.0000    0.5510];
                final_color    = 0.9*min(final_color   +0.0,1);
                final_color_sp = 0.9*min(final_color_sp+0.0,1);
            case 'ExpData_2D_1condition_sr_ecm_100samples__ExpDesign'
                R.n_generations = 30; % max: 30
                final_color    = [0.5020    1.0000    0.5020];
                final_color_sp = [0.5510    1.0000    0.5510];
                final_color    = 0.9*min(final_color   +0.0,1);
                final_color_sp = 0.9*min(final_color_sp+0.0,1);
            case 'ExpData_2D_1condition_sr_pf_ecm_100samples__ExpDesign'
                R.n_generations = 28; % max: 28
                final_color    = [0         0.5020    1.0000];
                final_color_sp = [0.0500    0.5510    1.0000];
        end
        R.generations = round(linspace(3,R.n_generations,5));
        
        % Options
        plot_true_parameters = false;
        
        % Limites & Ticks
        xtick_gen = 0:10:R.n_generations;
        
        ylim_acc = [0,0.6];
        ytick_acc = 0:0.1:0.6;

        ylim_dist = [1e-2,1e3];
        ytick_dist = 10.^[-2:3];
        
        ylim_fval = [1e2,1e6];
        ytick_fval = 10.^[2:6];

        E_sr.xlim = [0,20];
        E_sr.xtick = [0:5:20];
        E_sr.ylim = [0,500];
        E_sr.ytick = [0:100:500];

        E_pc.xlim = [0,300];
        E_pc.xtick = [0:100:300];
        E_pc.ylim = [0,100];
        E_pc.ytick = [0:25:100];
      
        E_emc.xlim = [0,300];
        E_emc.xtick = [0:100:300];
        E_emc.ylim = [0,0.2];
        E_emc.ytick = [0:0.05:0.2];
        
        switch filename_save
            case 'ExpData_2D_1condition_sr_100samples__ExpDesign'
                E_sr.plot_data = true;
                E_pc.plot_data = false;
                E_emc.plot_data = false;
            case 'ExpData_2D_1condition_sr_pf_100samples__ExpDesign'
                E_sr.plot_data = true;
                E_pc.plot_data = true;
                E_emc.plot_data = false;
            case 'ExpData_2D_1condition_sr_ecm_100samples__ExpDesign'
                E_sr.plot_data = true;
                E_pc.plot_data = false;
                E_emc.plot_data = true;
            case 'ExpData_2D_1condition_sr_pf_ecm_100samples__ExpDesign'
                E_sr.plot_data = true;
                E_pc.plot_data = true;
                E_emc.plot_data = true;
        end
        
        % Datasets
        R.exp{1} = E_sr;
        R.exp{2} = E_pc;
        R.exp{3} = E_emc;
end

%% Definition of colormap
fh = figure(999);

cm = colormap(gray(R.n_generations));
cm = 0.85*cm(end:-1:1,:);
cm(end,:) = final_color;

cm_sp = colormap(gray(R.n_generations));
cm_sp = 0.85*cm_sp(end:-1:1,:);
cm_sp(end,:) = final_color_sp;

close(fh)

%% Load the relevant datasets
% Measurement data
warning off;
load([filename_load '.mat']); 
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
    load([filepath filename_load num2str(i) '.mat' ],'raw_data','x','f','neval' );
    
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

%% Visualization of measurement data and simulation
% plot_options.lineplot = 'lines';
plot_options.lineplot = 'percentile intervals';

% Loop: Measurements
for k = 1:R.n_experiments
    % Open subplot
    figure
    set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 f_size_lp],'renderer','zbuffer');

    axes('Position',f_axis_lp);
    set(gca,'fontsize',ts_label,'box','on');
    
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
    for i = R.generations
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
        if R.exp{k}.plot_data
            ind = 1:thinning:size(D(k).raw,2);
            w_errorbar = mean(diff(D(k).x(ind(:))))*data_errorbar_w;
            for l = ind(:)'
                if isfield(D(k),'s')
                    plot(D(k).x(l),scaling*D(k).y(l),'.','color',data_col,'markersize',data_ms);
                    plot(D(k).x(l)*[1,1],scaling*(D(k).y(l)+[-1,1]*m_sigma*D(k).s(l)),...
                            'color',data_col,'linewidth',data_lw);
                    plot(D(k).x(l)+w_errorbar*[-1,+1]/2,scaling*(D(k).y(l)+[1,1]*m_sigma*D(k).s(l)),...
                            'color',data_col,'linewidth',data_lw);    
                    plot(D(k).x(l)+w_errorbar*[-1,+1]/2,scaling*(D(k).y(l)-[1,1]*m_sigma*D(k).s(l)),...
                            'color',data_col,'linewidth',data_lw);    
                else
                    plot(D(k).x(l),scaling*nanmean(D(k).raw(:,l)',2),'.','color',data_col,'markersize',data_ms);
                    plot(D(k).x(l)*[1,1],scaling*(nanmean(D(k).raw(:,l)',2)+[-1,1]*m_sigma*nanstd(D(k).raw(:,l))),...
                            'color',data_col,'linewidth',data_lw);
                end
            end
        end
        
        % Limits
        set(gca,'xlim',R.exp{k}.xlim,'xtick',R.exp{k}.xtick,...
                'ylim',R.exp{k}.ylim,'ytick',R.exp{k}.ytick,...
                'fontsize',ts_label);

        % Labels
        xlabel(R.exp{k}.xlabel,'fontsize',fs_label);
        ylabel(R.exp{k}.ylabel,'fontsize',fs_label);

        % Save figure
        print('-depsc2',resolution,['./figures/' filename_save '__' num2str(k) '_' R.exp{k}.filename '_' num2str(i) '.eps']);
    end
end

%% Visualization of parameter samples
% Open figure
figure;
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 f_size_sp],...
         'renderer','zbuffer');

% Loop: Generations
for i = R.generations
    % Scatter plot
    % Loop: parameter on y-axis
    for jy = 1:R.n_theta
    	% Loop: parameter on x-axis
        for jx = 1:(jy-1)
            % Open subplot
            subplottight(R.n_theta,R.n_theta,R.n_theta*(jy-1)+jx,margins_scatter_plot);

            % Lower and upper bound for grid
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
                             'facecolor',cm_sp(i,:),'edgecolor',cm_sp(i,:),'linewidth',2); hold on;
            if i > R.generations(end-1)
                plot(C(1,:),C(2,:),'.','color',cm_sp(i,:),'markersize',11);
            end
            
            % Boundary
            plot([options.lb(jx),options.ub(jx),options.ub(jx),options.lb(jx),options.lb(jx)],...
                 [options.lb(jy),options.lb(jy),options.ub(jy),options.ub(jy),options.lb(jy)],...
                 'k-','linewidth',0.5);
             
            % Plot true parameter
            if plot_true_parameters
%                 plot(options.true(jx)*[1,1],[options.lb(jy),options.ub(jy)],true_parameter_lt,...
%                     'color',true_parameter_color,'linewidth',true_parameter_lw);
%                 plot([options.lb(jx),options.ub(jx)],options.true(jy)*[1,1],true_parameter_lt,...
%                     'color',true_parameter_color,'linewidth',true_parameter_lw);
                plot(options.true(jx),options.true(jy),'x',...
                    'color',0.7*true_parameter_color,'linewidth',1.5,'markersize',3.5);
                plot(options.true(jx),options.true(jy),'x',...
                    'color',true_parameter_color,'linewidth',1,'markersize',3.5);
            end
    
            % Limits
            set(gca,'xlim',[options.lb(jx),options.ub(jx)],...
                    'ylim',[options.lb(jy),options.ub(jy)],...
                    'fontsize',ts_label);

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
        subplottight(R.n_theta,R.n_theta,R.n_theta*(j-1)+j,margins_scatter_plot);

        % Lower and upper ound for grid
        lbx = min(R.theta(j,:,i));
        ubx = max(R.theta(j,:,i));
        lbx = max(lbx-0.2*(ubx-lbx),options.lb(j));
        ubx = min(ubx+0.2*(ubx-lbx),options.ub(j));

        % Density estimation
        theta = linspace(lbx,ubx,n_grid);
        KDest = kde_simple(R.theta(j,:,i),theta);
        if i == R.generations(end)
            fill(theta([1,1:end,end]),[0,KDest,0]/max(KDest),...
                         'b','facecolor',cm_sp(i,:),'edgecolor',cm_sp(i,:),'linewidth',1); hold on;
        else
            plot(theta(1:end),KDest/max(KDest),'-','color',cm_sp(i,:),'linewidth',1.5); hold on;
        end
                 
        % Plot true parameter
        if plot_true_parameters
%             plot(options.true(j)*[1,1],[0,1.1],true_parameter_lt,...
%                 'color',true_parameter_color,'linewidth',true_parameter_lw);
            plot(options.true(j),0.03,'x',...
                'color',0.7*true_parameter_color,'linewidth',1.5,'markersize',3.5);
            plot(options.true(j),0.03,'x',...
                'color',true_parameter_color,'linewidth',1,'markersize',3.5);
        end
        
        % Limits
        set(gca,'xlim',[options.lb(j),options.ub(j)],...
                'ylim',[0,1.1],...
                'fontsize',ts_label);

        % Label
        if j == 1
            ylabel(options.legend{j},'fontsize',fs_label);
        else
            set(gca,'yticklabel',[]);
        end
        if j == R.n_theta
            xlabel(options.legend{j},'fontsize',fs_label);
        else
            set(gca,'xticklabel',[]);
        end
    end
    
    % Save figure
    print('-depsc2',resolution,['./figures/' filename_save '__scatterplot_' num2str(i) '.eps']);
end

%% Visualization of acceptance rate
% Open figure
figure;
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 f_size_acc],'renderer','zbuffer');

axes('Position',f_axis_acc);
set(gca,'fontsize',ts_label,'xtick',R.exp{k}.xtick,'ytick',R.exp{k}.ytick,'box','on');

% Loop: Generations
for i = 1:R.n_generations
    % Line and circles
    if i < R.n_generations
        plot([i-1,i],100*R.n_samples./R.n_eval([i,i+1]),'-','color',cm(i,:)); hold on;
    end
    plot(i-1,100*R.n_samples./R.n_eval(i),'.','color',cm(i,:),'markersize',10);
    % Arrows
    if max(i == R.generations)
        text(i-1,100*(R.n_samples./R.n_eval(i)+0.05),'\downarrow','color',cm(i,:),...
            'fontsize',fs_label,'HorizontalAlignment','center');
    end
end

% Limits
set(gca,'xlim',[-0.5-(R.n_generations-1)/20,R.n_generations-0.5+(R.n_generations-1)/20],...
        'xtick',xtick_gen,...
        'ylim',100*ylim_acc,'ytick',100*ytick_acc,...
        'fontsize',ts_label);

% Labels
xlabel('generation','fontsize',fs_label);
ylabel('acceptance rate [%]','fontsize',fs_label);

% Save figure
print('-depsc2',resolution,['./figures/' filename_save '__acceptance_rate.eps']);

%% Visualization of cumulative number of function evaluations
% Open figure
figure;
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 f_size_wlp],'renderer','zbuffer');

axes('Position',f_axis_wlp);
set(gca,'fontsize',ts_label,'xtick',R.exp{k}.xtick,'ytick',R.exp{k}.ytick,'box','on');

% Loop: Generations
for i = 1:R.n_generations
    % Line and circles
    if i < R.n_generations
        plot([i-1,i],[sum(R.n_eval(1:i)),sum(R.n_eval(1:i+1))],'-','color',cm(i,:)); hold on;
    end
    plot(i-1,sum(R.n_eval(1:i)),'.','color',cm(i,:),'markersize',10);
    % Arrows
    if max(i == R.generations)
        text(i-1,3*sum(R.n_eval(1:i)),'\downarrow','color',cm(i,:),...
            'fontsize',fs_label,'HorizontalAlignment','center');
    end
end

% Limits
set(gca,'xlim',[-0.5-(R.n_generations-1)/20,R.n_generations-0.5+(R.n_generations-1)/20],...
        'xtick',xtick_gen,...
        'yscale','log',...
        'ylim',ylim_fval,'ytick',ytick_fval,...
        'fontsize',ts_label);

% Labels
xlabel('generation','fontsize',fs_label);
ylabel({'cumulative number of','function evaluations'},'fontsize',fs_label);

% Save figure
print('-depsc2',resolution,['./figures/' filename_save '__cumulative_fval.eps']);

%% Visualization of distance function
% Open figure
figure;
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 f_size_acc],'renderer','zbuffer');

axes('Position',f_axis_acc);
set(gca,'fontsize',ts_label,'xtick',R.exp{k}.xtick,'ytick',R.exp{k}.ytick,'box','on');

% Loop: Generations
for i = 1:R.n_generations
    % Line and circles
    if i < R.n_generations
        semilogy([i-1,i],[median(R.f(:,i)),median(R.f(:,i+1))],'-','color',cm(i,:)); hold on;
    end
    semilogy((i-1)*ones(R.n_samples,1),R.f(:,i),'.','color',cm(i,:),'markersize',4); 
    % Arrows
    if max(i == R.generations)
        text(i-1,3*max(R.f(:,i)),'\downarrow','color',cm(i,:),...
            'fontsize',fs_label,'HorizontalAlignment','center');
    end
end

% Limits
set(gca,'xlim',[-0.5-(R.n_generations-1)/20,R.n_generations-0.5+(R.n_generations-1)/20],...
        'xtick',xtick_gen,...
        'ylim',ylim_dist,'ytick',ytick_dist,...
        'fontsize',ts_label);

% Labels
xlabel('generation','fontsize',fs_label);
ylabel('distance measure','fontsize',fs_label);

% Save figure
print('-depsc2',resolution,['./figures/' filename_save '__distance_function.eps']);

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
