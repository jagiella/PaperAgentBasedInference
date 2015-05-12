function scatterPlotMatrix(x, lb, ub, true, xr, varargin)
    
    if( nargin < 6)
        options.color   = [0, 0, 1]; % blue
        options.true    = 1;
        options.contour = 1
    else
        options = varargin{1};
    end

    [d,n] = size(x)
    %[d,nr] = size(xr);
    %lb, ub, true
    
    colmap = min( jet(10) + 0.2, 1);
    %colmap = bluewhitered( 10);
    colmap2 = colormapwhite( 10, colmap(1,:), colmap(10,:));
    
    
    
    R = corrcoef(x');
    lab = {'k_{div}^{max}', 'l^{crit}', 'l^{init}', 'q_{init}', 'k_{e}^{pro}', 'k_{e}^{deg}', 'e^{crit}'};
    
    for i=1:d
        for j=1:d
            idx = i+(j-1)*d;
            subplottight(d,d,idx);
        
            %for j=(i+1):d
            if j>i
                %% SAMPLES & POSTERIOR
                %subplottight(d,d,i+(j-1)*d);

                %% plot kernel density estimate (KDE) of points
                if( options.contour)
                    [xt,yt] = meshgrid(linspace(lb(i),ub(i),100),linspace(lb(j),ub(j),100));
                    zt = reshape( kde_simple( x([i j], :), [xt(:),yt(:)]' ), size(xt));
                    v = [0.001 0.01 0.1 0.25 0.5, 0.75, 0.9, 0.95 0.99 1]*max(zt(:));
                    contour(xt,yt,zt, v);
                    hold on;
                    %clear xt yt zt;
                end

                %% profile
    %             [~, idx] = max(zt,[],1);
    %             IDX = idx+(0:100:9900);
    %             plot( xt(IDX), yt(IDX), 'r-');

                %% plot points
    %             if( ~isempty(xr))
    %                 plot( xr(i,:), xr(j,:), 'g.'); 
    %             end        
    %            hold on;
                %plot( x(i,:), x(j,:), 'b.'); hold on;
                
                %plot( x(i,:), x(j,:), 'Color', colmap(2*options.iterationIdx,:), 'LineStyle','none', 'Marker','.','markersize',20); hold on;
                plot( x(i,:), x(j,:), 'Color', options.color, 'LineStyle','none', 'Marker','.','markersize',20); hold on;
                %plot( 10.^x(i,:), 10.^x(j,:), 'Color', colmap(2*options.iterationIdx,:), 'LineStyle','none', 'Marker','.','markersize',20); hold on;
                %plot( x(i,:), x(j,:), 'Color', options.color, 'LineStyle','none', 'Marker','.'); hold on;
                
                %set(gca,'xscale','log','yscale','log');
                if( options.true)
                    if( ~isempty(true))
                        %plot(10.^true(i),10.^true(j), 'r*'); 
                        plot(true(i),true(j), 'r*'); 
                    end
                end

                %% set axis ranges
                xlim([lb(i),ub(i)]);
                ylim([lb(j),ub(j)]);
                %xlim([10.^lb(i),10.^ub(i)]);
                %ylim([10.^lb(j),10.^ub(j)]);
               % axis([lb(i) ub(i) lb(j) ub(j)]);

            end

            if i==j
                %subplottight(d,d,i+(i-1)*d);

%                 %% DIAGONAL
%                 % plot histogram
%                 nb = 20;
%                 xb = linspace(lb(i),ub(i),nb);
%                 yb = hist(x(i, :), xb); 
%                 bar(xb, yb/n*nb ); hold on;
% 
%                 % plot true
%                 if( ~isempty(true))
%                     true(i)
%                     plot( true(i), 0, 'r*');
%                 end
% 
%                 % plot KDE
%                 xt = linspace(lb(i), ub(i), 100);
%                 zt = kde_simple( x(i, :), xt );
%                 plot(xt, zt, 'r-');
%                 xlim([lb(i) ub(i)]);
% 
%                 if( i==d)
%                     xlabel(lab{i});
%                 end
% 
%                 clear xb yb xt zt;
%                  hold off;
                % plot KDE
                xt = linspace(lb(i), ub(i), 1000);
                zt = kde_simple( x(i, :), xt );

                fill([lb(i),xt,ub(i)],[0,zt,0],options.color)
                %fill([lb(i),xt,ub(i)],[0,zt,0],colmap(2*options.iterationIdx,:))
                xlim([lb(i) ub(i)]);

                if( i==d)
                    xlabel(lab{i});
                end

                clear xb yb xt zt;
                 hold off;
            end

            %% CORRELATION COEFFICIENT
            %for j=1:(i-1)
            if j<i
                %subplottight(d,d,i+(j-1)*d);
%                 colmap = colormap(jet(10));
%                 colmap = 1*min(colmap + 0.2,1);
%                fill([0 1 1 0], [0 0 1 1]', colmap( ceil(5*(1+R(i,j))), : ) ) 
%                fill([0 1 1 0], [0 0 1 1]', [1,1,1] ) 

%                cmap=jet(100);
%                fill([0 1 1 0], [0 0 1 1]', cmap( floor(50*(1+R(i,j))), : ) ) 
                fill([0 1 1 0], [0 0 1 1]', colmap2( ceil(5*(1+R(i,j))), : ) ) 
                %fill([0 1 1 0], [0 0 1 1]', bluewhitered_scalar(R(i,j))) 
                set(gca,'YTickLabel',[],'XTickLabel',[]);
                if options.iterationIdx == 5
                    text(0.5,0.5, num2str(R(i,j), '%.2f'), 'FontSize', 15, 'HorizontalAlignment','center'); hold on;
                end
            end
            
            if( i==1)
                ylabel(['log_{10}' lab{j}]);
            end
            if( j==d)
                xlabel(['log_{10}' lab{i}]);
            end
        
            %hold off;

            %xlabel(['i=' num2str(i) ', j=' num2str(j)]);
            if(  i ~= 1)
                set(gca,'YTickLabel',[]);            
            end
            if(  j ~= d)
                set(gca,'XTickLabel',[]);            
            end
        end
    end
    
    
    clear lab cmap R d n;
end

function c = bluewhitered_scalar( value)

if value < 0
    I = 1+value;
    c = [I I 1];
else
    I = 1-value;
    c = [1 I I];
end

end

function cmap = bluewhitered( n)

cmap=[];

a=n/2;
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

function cmap = c1whitec2( n, c1, c2)

cmap=[];

a=floor(n/2);
b=a+1;

for i=1:a % excludes white
    alpha = (i-1)/a;
    cmap(i,:) = alpha * [1 1 1] + (1-alpha) * c1;
end
for i=b:n % includes white
    alpha = (i-b)/(n-b);
    cmap(i,:) = (1-alpha) * [1 1 1] + alpha * c2;
end

end


function h = subplottight(n,m,i)
    [c,r] = ind2sub([m n], i);
    margins=[0.15 0.15 0.1 0.1];
    
    dx=(1-margins(1)-margins(3))/m;
    dy=(1-margins(2)-margins(4))/n;
    gap=0.01;
    
    ax = subplot('Position', [margins(1) + dx*(c-1), 1 - margins(4) - dy*r, dx-gap, dy-gap]);
    if(c>1)
        set(ax, 'YTickLabel', []);
    end
    if(r<n)
        set(ax, 'XTickLabel', []);
    end

    if(nargout > 0)
      h = ax;
    end
end