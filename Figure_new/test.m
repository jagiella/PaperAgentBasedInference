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
                ylabel(options.legend{jx},'fontsize',fs_label);
            else
                set(gca,'yticklabel',[]);
            end
            if jy == R.n_theta
                xlabel(options.legend{jy},'fontsize',fs_label);
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
    print('-depsc2',['./figure/' filename '__scatterplot_' num2str(i) '.eps']);
end

