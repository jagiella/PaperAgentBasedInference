function plotCurveDensity(x, Y)
    
%     x  = 0:1:100; 
    
     y  = linspace(min(Y(:)),max(Y(:)),100);
%     dy = normrnd( ones(length(x),1000), ones(length(x),1000) );
%     Y = cumsum(dy);
    
%   figure(10); clf();
    
   
    size(Y)
    
    %subplot(2,1, 2)
    h2d = hist(Y', y) / 1000 * 256 * 1;% / size(Y,2)
    image(x,y,h2d)
    colormap(1-gray)
    set(gca,'YDir','normal')
    
    %subplot(2,1, 1)
    %hold on
    %errorbar(x, mean(Y,2), std(Y,[],2))
end