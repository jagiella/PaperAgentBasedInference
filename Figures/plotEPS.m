function plotEPS( h, filename, width, height)
set(h,'Position',[width height width height]);
set(h,'PaperPositionMode','auto');
print(h,'-depsc2',filename);
end