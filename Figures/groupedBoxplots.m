function groupedBoxplots( m,n, X, Groups, GroupLabels, SubGroups, SubGroupLabels, ColorMap)

nX1 = size(X,1); % nb of data sets
nX2 = size(X,2); % sample size


UniqueGroups = unique(Groups);
nG = length(UniqueGroups);

for p=1:nG
    h_subfig = subplot( m,n, p);
    
    Xp = X(:, Groups==UniqueGroups(p) );
    Gp = SubGroups( Groups==UniqueGroups(p) );
    
    
    boxplot( h_subfig, Xp, 'colorgroup', Gp, 'symbol','+', 'orientation', 'horizontal', 'colors', ColorMap); hold on;
    ylabel(  h_subfig, GroupLabels{p});
    set(     h_subfig,'TickLabelInterpreter','tex');
    
    if p==nG
        legend( findall(gcf,'tag','legend'), SubGroupLabels, 'interpreter', 'none', 'Location', 'southoutside');
    end
end



% for p=1:size(X,1)
%     h_fig = gcf;
%     h_subfig = subplot( m,n, p);
% 
%     Xp = X(:,p:size(x,1):size(X,2));
%     
%     boxplot( h_subfig, Xp, 'symbol','+', 'orientation', 'horizontal', 'color', colmap); hold on;
%     %plot( h, [options.true(p), options.true(p)], size(prefixes,2) * [0, 1] + 0.5, 'k','LineWidth',2)
%     set(h,'yticklabel',[])
%     ylabel( h_subfig,options.legend{p});
%     set( h_subfig,'TickLabelInterpreter','tex');
% 
%     xl = xlim(h);
%     xl = [ min( options.true(p), xl(1) - 0.1*(xl(2)-xl(1)))   max( options.true(p), xl(2) + 0.1*(xl(2)-xl(1))) ];
%     xlim(h,xl);
% 
%     hold off;
% end

end