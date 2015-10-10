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

% Figure size & axis position
f_size_comparison = [18,6];
f_axis_comparison = [0.35 0.14 0.61 0.84];

% Printing options
resolution = '-r300';

%% Datasets
% experimental data, 2D model, 1 condition, all data, 100 samples
filename1_load = 'Tumor2dGCKI67ECM'; 
filename1_save = 'ExpData_2D_1condition_sr_pf_ecm_100samples'; 

% experimental data, 3D model, 1 condition, all data, 100 samples
filename2_load = 'Tumor3dGCKI67ECM'; 
filename2_save = 'ExpData_3D_1condition_sr_pf_ecm_100samples'; 

% Generations
R1.n_generations = 26; % max: 28
R2.n_generations = 26; % max: 26

R1.generations = R1.n_generations;
R2.generations = R2.n_generations;

% %% Definition of colormap
% fh = figure(999);
% cm = colormap(gray(R.n_generations)); cm = 0.85*cm(end:-1:1,:); cm(end,:) = [1,0,0];
% % cm = colormap(jet(R.n_generations));
% close(fh)

%% Load datasets 1
R = R1;
filename_load = filename1_load;

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

% Assignment
R1 = R;

%% Load datasets 1
R = R2;
filename_load = filename2_load;

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

% Assignment
R2 = R;

%% Visualization of the parameter uncertainties
options.conf_level = [0.99,0.95,0.8,0];

parameters1.name = options.legend;
parameters1.number = length(options.legend);
parameters1.S.par = R1.theta(:,:,end);
parameters1 = getParameterConfidenceIntervals(parameters1,1-options.conf_level);
parameters1.min = options.lb;
parameters1.max = options.ub;

parameters2.name = options.legend;
parameters2.number = length(options.legend);
parameters2.S.par = R2.theta(:,:,end);
parameters2 = getParameterConfidenceIntervals(parameters2,1-options.conf_level);
parameters2.min = options.lb;
parameters2.max = options.ub;

% Open figure
figure;
set(gcf, 'PaperUnits','centimeters', 'PaperPosition',[0 0 f_size_comparison],'renderer','zbuffer');

axes('Position',f_axis_comparison);
set(gca,'fontsize',ts_label,'box','on');

% Selection of color-coding
cm2(1,:) = [255,204,102]/255;
cm2(3,:) = [255,128,0]/255;
cm2(2,:) = 0.5*(cm2(1,:)+cm2(3,:));
cm2(3,:) = cm2(2,:);
cm2(2,:) = 0.5*(cm2(1,:)+cm2(3,:));

cm1 = cm2(:,[3,2,1]);

%cm1 = colormap(gray(length(options.conf_level))); cm1 = min(cm1(end:-1:1,:)+0,1); cm1(:,1) = 0.9;
%cm2 = colormap(gray(length(options.conf_level))); cm2 = min(cm2(end:-1:1,:)+0,1); cm2(:,3) = 0.9;

% Width of bars for different confidence levels
options.conf_level_w = 0.5*((1-options.conf_level)/max(1-options.conf_level)).^(1/3);

lb = min(options.lb);
ub = max(options.ub);

% Loop: Parameters
for i = 1:parameters1.number
    for j = 1:(length(options.conf_level)-1)
        % 2D
        fill(10.^parameters1.CI.S(i,[1,2,2,1,1],j),...
             2*(parameters1.number+1-i)+options.conf_level_w(j)*[-1,-1,1,1,-1],'b',...
             'facecolor',cm1(j,:),...
             'edgecolor','k','linewidth',0.2); hold on;
        % 3D
        fill(10.^parameters2.CI.S(i,[1,2,2,1,1],j),...
             2*(parameters2.number+1-i)-1+options.conf_level_w(j)*[-1,-1,1,1,-1],'b',...
             'facecolor',cm2(j,:),...
             'edgecolor','k','linewidth',0.2); hold on;
    end
    % 2D
    plot(10.^parameters1.CI.S(i,[1,2],end),2*(parameters1.number+1-i)+[-0.35,+0.35],'k-');
    % 3D
    plot(10.^parameters2.CI.S(i,[1,2],end),2*(parameters2.number+1-i)-1+[-0.35,+0.35],'k-');
    
%     % Regions which are ruled out
%     fill([lb,parameters2.min(i),parameters2.min(i),lb,lb],...
%          2*(parameters2.number+1-i)+[-1.5,-1.5,0.5,0.5,-1.5],'b',...
%          'facecolor',0.95*[1,1,1],...
%          'edgecolor',0*[1,1,1],'linewidth',0.5); hold on;
%     fill([parameters2.max(i),ub,ub,parameters2.max(i),parameters2.max(i)],...
%          2*(parameters2.number+1-i)+[-1.5,-1.5,0.5,0.5,-1.5],'b',...
%          'facecolor',0.95*[1,1,1],...
%          'edgecolor',0*[1,1,1],'linewidth',0.5); hold on;
    
    % Seperating line
    plot(10.^[lb,ub],2*(parameters1.number+1-i)-1.5*[1,1],'--','color',0.8*[1,1,1],'linewidth',1);
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
    text(10.^(-5.3),2*i-0.5,parameter_label{i},'HorizontalAlignment','right','fontsize',fs_label)
end
set(gca,'xlim',10.^[lb,ub],'xscale','log',...
        'ylim',[0.5,2*parameters1.number+0.5],'ytick',[],...
        'fontsize',ts_label);

% Labels
xlabel('parameters','fontsize',fs_label);


% Save figure
print('-depsc2',resolution,['./figures/figure__2D_3D_comparison.eps']);
    
