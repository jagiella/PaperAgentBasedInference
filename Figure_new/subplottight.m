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