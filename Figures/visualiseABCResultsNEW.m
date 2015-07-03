function visualiseABCResultsNEW( prefix)
figure(1); clf(); hold off;
figure(2); clf(); hold off;
figure(3); clf(); hold off;
figure(4); clf(); hold off;
figure(5); clf(); hold off;
figure(6); clf(); hold off;


addpath('..'); % kde_simple(...)
addpath ../Hash; % DataHash(...)

%% GET GENERAL STUFF
FontSize = 10;
true=[];
truef=[];
load( [ prefix '.mat' ] ); 

h_sp = figure(4); clf()
     
%% LOAD RAW DATA
likelihood = zeros(100,1);

if( isfield( options, 'data'))
    for i=1:length(options.data) % each parameter
        if( isfield(options.data(i), 'raw'))
            n = size(options.data(i).raw,1);
            for j=1:n % each independent experiment
                likelihood(j) = likelihood(j) + 0.5 * sum( ( (options.data(i).y - options.data(i).raw(j,:)) ./ options.data(i).s ).^2)...
                    / length(options.data(i).y);
            end
        else
            n = 10;

            for j=1:n
                RND = normrnd(0,1);
                likelihood(j) = likelihood(j) + 0.5 *  RND^2;
                %likelihood(j) = likelihood(j) + 0.5 *  nansum( ( normrnd(0, options.data(i).s) ./ options.data(i).s ).^2) / length(options.data(i).y);
            end
        end
    end
    likelihood = likelihood(1:n)
end
    



%% PLOT
median_f=[];


%
%% PLOT RAW DATA
h_fi = figure(5); clf();


max_i=0;
while( exist([prefix num2str(max_i+1) '.mat'], 'file'))
    max_i=max_i+1;
end

color_map_size = round( max_i);
color_map = min( jet(color_map_size) +0.15, 1);
%color_map = min( jet(10) + 0.2, 1);
%color_map = colormapwhite( 5, color_map(1,:), color_map(10,:))
%color_map = jet(5)
c=1;
%IterationArray=[1 10 20 30 40];
%IterationArray = [1 4 7 10 14];
IterationArray = round(linspace(1,max_i,5));
for i=IterationArray(end)
    
    load( [ prefix num2str(i) '.mat' ], 'raw_data_x', 'raw_data_y', 'raw_data', 'x','f', 'neval' ); 

   %10.^mean(x,2)
    
    figure(10); clf
    n = size(raw_data,2);
    nx= ceil(n^0.5);
    ny= ceil(n/nx);

    for j=1:n
        subplot(nx,ny, j);
        raw_Y=[];
        raw_x=[];
        for k=1:size(raw_data,1)
            try
                %plot(raw_data(k,j).x, raw_data(k,j).y, 'color', color_map(i,:),'LineWidth',2); hold on;
                raw_Y(k,:) = raw_data(k,j).y;
                raw_x      = raw_data(k,j).x;
            catch
                warning('error reading raw data')
            end
        end
        plotCurveDensity( raw_x, raw_Y');
        hold on;
        try
            num_data_points = length(options.data(j).y);
            step_size = ceil(num_data_points/20);
            indices = 1:step_size:num_data_points;
            errorbar(options.data(j).x(indices,:), options.data(j).y(indices), options.data(j).s(indices), 'k-'); hold off;
        end

%         xlim([min(options.data(j).x), max(options.data(j).x)]);
%         idx = (options.data(j).s ~= 0);
%         ylim([min(options.data(j).y(idx)), max(options.data(j).y(idx) + options.data(j).s(idx))]);
%         xlabel(options.xlabel(j))
%         ylabel(options.ylabel(j))
    end
end
%return

 for i=IterationArray
    
    load( [ prefix num2str(i) '.mat' ], 'raw_data_x', 'raw_data_y', 'raw_data', 'x','f', 'neval' ); 

  
    %% PLOT RAW SIM DATA
    h_fi = figure(5);
    n = size(raw_data,2);
    nx= ceil(n^0.5);
    ny= ceil(n/nx);
    for j=1:n
        %subplot(1,n, j);
        subplot(nx,ny, j);
        for k=1:size(raw_data,1)
            try
                plot(raw_data(k,j).x, raw_data(k,j).y, 'color', color_map(i,:),'LineWidth',2); hold on;
            end
        end
        try
            num_data_points = length(options.data(j).y);
            step_size = ceil(num_data_points/20);
            indices = 1:step_size:num_data_points;
            errorbar(options.data(j).x(indices,:), options.data(j).y(indices), options.data(j).s(indices), 'k-'); hold off;
        end
        
        xlim([min(options.data(j).x), max(options.data(j).x)]);
        idx = (options.data(j).s ~= 0);
        ylim([min(options.data(j).y(idx)), max(options.data(j).y(idx) + options.data(j).s(idx))]);
        xlabel(options.xlabel(j))
        ylabel(options.ylabel(j))
        hold on;
    end

    
    if( i==IterationArray(end))
    end


    %saveas( gcf, [ prefix num2str(i) '_fits.png' ], 'png');

    %length(f);
    
    
    %% SCATTER MATRIX PLOT
     h_sp = figure(4); %clf()
     %scatterPlotMatrix(x, lb, ub, true, xr); %drawnow;
     scatter_plot_options.color   = color_map(i,:);
     scatter_plot_options.true    = 1;
     scatter_plot_options.contour = 0;
     scatter_plot_options.iterationIdx = find(i==IterationArray);
     scatterPlotMatrix(x, options.lb, options.ub, options.true, [], scatter_plot_options); %drawnow;
     
     
     
     %
     c=c+1;
end


data_it=[];
data_fe=[];
data_ar=[];
data_max_ob=[];
i = 1;
while exist( [ prefix num2str(i) '.mat' ], 'file' )
    
    load( [ prefix num2str(i) '.mat' ], 'raw_data_x', 'raw_data_y', 'raw_data', 'x','f', 'neval' ); 
     
    %% PLOT OBJ FUNC
    h_of = figure(2);
    median_f(i) = median(f);
    data_max_ob(i) = max(f);
    plot(size(median_f,2)*ones(size(f)),f,'o','color', color_map(i,:)); hold on;
    %plot(size(median_f,2)*ones(size(f)),f,'bo'); hold on;
    plot(median_f,'k-','LineWidth',2); 
    %plot(median_f,'color', color_map(5,:),'LineWidth',2); 
    %if( USE_RAW_DATA)
        median_true_f(i) = median(likelihood);
        plot(median_true_f, 'k--','LineWidth',2); 
        %plot(median_true_f, '-', 'color', color_map(end,:),'LineWidth',2); 
    %end
    set(gca,'yscale','log');
    legend('obj. fun. of sampled parameters', 'median obj. fun. of sampled parameters', 'median obj. fun. for true parameters)');
    xlabel('Iteration'); ylabel('Objective function');

    %% PLOT ACCEPTANCE RATE
%     h_ar = figure(3);
    data_ar(end+1) = length(f)/double(neval);
%     plot( size(median_f,2), length(f)/double(neval), 'bo-'); 
%     ylim([0 1]);
%     xlabel('Iteration', 'FontSize', FontSize ); ylabel('Acceptance rate');
%     hold on;

    %% PLOT #FUNC EVAL
%     h_fe = figure(6);
     data_it(end+1) = i;
     data_fe(end+1) = neval;
%     semilogy( size(median_f,2), neval, 'bo-'); 
%     %ylim([0 1]);
%     xlabel('Iteration', 'FontSize', FontSize ); ylabel('Number of Function Evaluations / Iteration');
%     hold on;

    i = i + 1;
end

%% PLOT ACCEPTANCE RATE
h_ar = figure(3);
%plot( data_it, data_ar, 'bo-');
for i=1:length(data_ar)
    plot( data_it(i), data_ar(i), 'o','color', color_map(i,:)); hold on;
end
plot( data_it, data_ar, 'k-','LineWidth',2);
ylim([0 1]);
xlabel('Iteration', 'FontSize', FontSize ); ylabel('Acceptance rate');
hold on;


%% PLOT #FUNC EVAL
h_fe = figure(6);
semilogy( data_it, data_fe, 'bo-'); 
%ylim([0 1]);
xlabel('Iteration', 'FontSize', FontSize ); ylabel('Number of Function Evaluations / Iteration');
hold on;

%% OBJ FUNC
% figure(h_of);
% for i=IterationArray
%     apos = [ i, data_max_ob(i)];
%     annotation('arrow', [i], [0]);
% end

[d,n] = size(x);
h_sp.Position = [0 0 1 1]*d*100
saveas( h_sp, [ prefix '_NEW_scatterPlotMatrix.png' ], 'png');
%plotEPS(h_sp,[ prefix '-scatterPlotMatrix.eps' ],800, 800);

%h_of = figure(2);
saveas( h_of, [ prefix '_NEW_objFunc.png' ], 'png');
%plotEPS(h_of,[ prefix '-objFunc.eps'],400, 300);

%h_ar = figure(3);
saveas( h_ar, [ prefix '_NEW_acceptanceRate.png' ], 'png');
%plotEPS(h_ar,[ prefix '-acceptanceRate.eps'],400, 300);

saveas( h_fi, [ prefix '_NEW_fits.png' ], 'png');
%plotEPS(h_fi,[ prefix '-fits.eps'],800, 600/3);

saveas( h_fe, [ prefix '_NEW_functionEvaluations.png' ], 'png');
%plotEPS(h_fi,'/Users/jagiella/Dropbox/Work/PUBLICATIONS/Paper-Helmholtz-ParameterEstimation/Figures/FitToyData1.eps',800, 800/3);



end



function cmap = bluewhitered( n)

cmap=[];

a=floor(n/2);
b=a+1;
for i=1:a
    I = (i-1)/double(n-b);
    cmap(i,:) = [I I 1];
end
for i=b:n
    I = 1 - (i-b)/(n-b);
    cmap(i,:) = [1 I I];
end

end