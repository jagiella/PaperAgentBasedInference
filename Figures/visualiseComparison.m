function visualiseComparison( prefixes)
    %close all;
    %clear all;
[substr, c1, c2] = diffSubStr( prefixes{:})

    AltLabels_x = {'Time (days)', 'Distance to Spheroid Border (\mu m)', 'Distance to Spheroid Border (\mu m)'};
    AltLabels_y = {'Spheroid Radius (\mu m)', 'Ki67-positive Cell Fraction', 'Extra-Cellular Matrix Density'};
    
    colmap = min( jet(10) + 0.1, 1);
    
    
    switch( size(prefixes,2) )
        case 2
            colmap = colmap( [1,10], :);
        case 3
            colmap = colmap( [1,4,10], :);
        case 4
            colmap = colmap( [1,4,7,10], :);
    end
    %colmap = colmap( [1,4,10], :);
    %colmap2 = colormapwhite( 10, colmap(1,:), colmap(10,:));
    
    colors='brgk';
    FontSize=20;
    
    size(prefixes)
    if( iscell(prefixes))
        figure(1); clf();
        h_fe = figure(2); clf();
        h_e = figure(3); clf();
        X = [];
        G = [];
        L = [];
        Lab = {};
        options=[];
        for i=1:size(prefixes,2)
            %subplot(2,2, i)
            %% find last
            it = [];
            fe_vs_it = [];
            e_vs_it = [];
            j=0;
            %exist([prefixes{i} num2str(j+1) '.mat'], 'file')
            while( exist([prefixes{i} num2str(j+1) '.mat'], 'file'))
                j=j+1;
                disp([prefixes{i} num2str(j) '.mat']);
                %% make computational cost comparision
                load( [prefixes{i} num2str(j) '.mat'], 'neval','epsilon');
                it(j) = j;
                fe_vs_it(j) = neval;
                e_vs_it(j)  = epsilon;
                %figure(2);
                %semilogy(j, neval, [colors(i) 'o']); hold on;
            end
            figure(2);
            %semilogy(it, fe_vs_it, [colors(i) 'o-']); hold on;
            semilogy(it, fe_vs_it, '-','LineWidth',2, 'color', colmap(i,:)); hold on;
            
            figure(3);
            semilogy(it, e_vs_it, '-','LineWidth',2, 'color', colmap(i,:)); hold on;
            
            disp([prefixes{i} num2str(j) '.mat']);
            
            load( [prefixes{i} num2str(j) '.mat'], 'x', 'raw_data')
            load( [prefixes{i} '.mat'], 'options')
            %prctile(x,[05 10 25 50 75 90 95],2)

            %boxplot( x', 'orientation', 'horizontal' ); hold on;
            multiple = 200 / size(x', 1);
            temp = [];
            for m=1:multiple
                % absolute error to true
                %temp = [temp; bsxfun(@minus, x', options.true)]; %
                % parameter values
                temp = [temp; x'];
            end
            %X = [X temp];
            X = [temp X];
            
            G = [G -i*ones(1,size(x,1))];
            L = [L 1:size(x,1)];
            %L = horzcat(L,  options.legend{1:size(x,1)});
            %size(G), size(X), size(L), size(Lab)
            
            
            %% final fit
            h_ff = figure(4);
            for j=1:length(options.data)
                %subplot(2,2, j);
                subplot(1,length(options.data), j);
                %for i=1:size(prefixes,2)

                % 
                xmax=0;
                for k=1:size(raw_data,1)
                    xmax = max( [xmax; raw_data(k,j).x(:)]);
                    plot( raw_data(k,j).x, raw_data(k,j).y, 'color', colmap(i,:)); hold on;
                    %plot( raw_data(k,j).x, raw_data(k,j).y, colors(i)); hold on;
                end
                num_data_points = length(options.data(j).y);
                step_size = ceil(num_data_points/20);
                indices = 1:step_size:num_data_points; % sparsen out
                %indices = options.data(j).x(indices,:) <= xmax; % smaller than data
                errorbar(options.data(j).x(indices,:), options.data(j).y(indices), options.data(j).s(indices), 'k-'); 
 
                %errorbar( options.data(j).x, options.data(j).y, options.data(j).s, 'k-');
                %xlim( xl(1), xmax);
                %xlim([min(options.data(j).x), max(options.data(j).x)]);
                idx = (options.data(j).s ~= 0);
                ylim([min(options.data(j).y(idx)), max(options.data(j).y(idx) + options.data(j).s(idx)*2)]);
                %xlabel(options.xlabel(j))
                %ylabel(options.ylabel(j))
                xlabel(AltLabels_x(j));
                ylabel(AltLabels_y(j));
                box on;
                %hold off;
            end
            
            

        end
        
        disp('Plot function evaluations')
        h_fe = figure(2);
        legend( findobj(h_fe,'Tag','Box'), prefixes, 'interpreter', 'none', 'Location', 'northoutside');
        xlabel('iteration');
        ylabel('function evaluations');
        
        disp('Plot objective function')
        h_e = figure(3);
        legend( findobj(h_e,'Tag','Box'), prefixes, 'interpreter', 'none', 'Location', 'northoutside');
        xlabel('iteration');
        ylabel('objective function threshold \epsilon');
        
        disp('Prepare Labels')
        for p=1:size(x,1) % parameters
            for m=1:size(prefixes,2) % models
                if( isfield(options, 'legend'))
                    if( m==size(prefixes,2))
                        Lab = horzcat(Lab, ['' options.legend{p} '']);
                    else
                        Lab = horzcat(Lab, ' ');
                    end
                else
                    AltLab = {'k_{div}^{max}', 'l^{crit}', 'l^{init}', 'q_{init}', 'k_{e}^{pro}', 'k_{e}^{deg}', 'e^{crit}'};
                    Lab = horzcat(Lab, ['' AltLab{p} '']);
                end
            end
        end
        
        identifiabilityTab( X, {L, G})
        
        h_bp = figure(1);
        gap_size =  0;
        boxplot( X, {L, G}, 'orientation', 'horizontal', 'factorgap',gap_size, 'factorseparator', 1, 'color', colmap(size(prefixes,2):-1:1,:), 'Labels', Lab);
        %boxplot( X, {L, G}, 'orientation', 'horizontal', 'factorgap',gap_size, 'factorseparator', 1, 'color', colors(size(prefixes,2):-1:1), 'Labels', Lab);
        legend(findobj(h_bp,'Tag','Box'), prefixes, 'interpreter', 'none', 'Location', 'northoutside');
        hold on;
        for p=1:size(x,1)
            %plot( [options.true(p), options.true(p)], size(prefixes,2) * [p-1, p] + 0.5 + gap_size, 'k')
        end
        xlabel('parameter value');
        ylabel('parameter');
        set(gca,'TickLabelInterpreter','tex');
        

        %% Independent Boxplots
        
        figure(7)
        groupedBoxplots( ceil(size(x,1)/2),2, X, L, options.legend, G, prefixes, colmap);
        
        %figure(5); clf();
        for p=1:size(x,1)
            h_5 = figure(5);
            h = subplot(ceil(size(x,1)/2),2, p);
            %h = subplot( size(x,1), 1, p);
            
            Xp = X(:,p:size(x,1):size(X,2));
            %plot( h, [options.true(p), options.true(p)], size(prefixes,2) * [0, 1] + 0.5, 'k')
            boxplot( h, Xp, 'symbol','+', 'orientation', 'horizontal', 'color', colmap); hold on;
            %boxplot( X(:,p:size(x,1):size(X,2)), 'orientation', 'horizontal', 'color', colors); hold on;
            plot( h, [options.true(p), options.true(p)], size(prefixes,2) * [0, 1] + 0.5, 'k','LineWidth',2)
            set(h,'yticklabel',[])
%             if(p==size(x,1))
%                 xlabel(h, 'parameter values');
%             end
            ylabel(h,options.legend{p});
            set(h,'TickLabelInterpreter','tex');
            %set(h,'TickLabelInterpreter','tex');
            
            xl = xlim(h);
            xl = [ min( options.true(p), xl(1) - 0.1*(xl(2)-xl(1)))   max( options.true(p), xl(2) + 0.1*(xl(2)-xl(1))) ];
            xlim(h,xl);
            %xlim([min( [Xp(:); options.true(p)]), max([Xp(:); options.true(p)])]);
            
            %if p==size(x,1)
            %    legend( subplot(ceil(size(x,1)/2),2, p+1), Lab);
            %end
            
            hold off;
        end
        

    end
    %h_bp = figure(1);
    %saveas( gcf, [prefixes{:} '.png'], 'png');

    LatexLabRow = {'$k_{div}^{max}$', '$l^{crit}$', '$l^{init}$', '$q_{init}$', '$k_{e}^{pro}$', '$k_{e}^{deg}$', '$e^{crit}$'};
    [m,s, idTab] = identifiabilityTab(X, {L,G})
    LatexLabCol = diffSubStr( prefixes{:})
    matrix2latex([ c1 'XXX' c2 'uncertainty.tex'], idTab, LatexLabCol, LatexLabRow)
    matrix2latex([ c1 'XXX' c2 'identifiability.tex'], idTab, LatexLabCol, LatexLabRow, 1)
    
    plotEPS( h_fe, [ c1 'XXX' c2 'functionEvaluations.eps'], 800/2, 300);
    plotEPS( h_fe, [ c1 'XXX' c2 'functionEvaluations.eps'], 800/2, 300);
    plotEPS( h_e,  [ c1 'XXX' c2 'objectiveFunction.eps'], 800/2, 300);
    plotEPS( h_5,  [ c1 'XXX' c2 'independentBoxplots.eps'], 800, 300);
    plotEPS( h_ff, [ c1 'XXX' c2 'fit.eps'], 800, 200);

end