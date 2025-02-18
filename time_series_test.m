datestr(time(1),'yyyymmdd HH:MM:SS')

figure(1)
h = pcolor(time,range,Offline_Temp_Spatial_Avg');
set(h, 'EdgeColor', 'none');
caxis([1e-2 1e4]);
datetick('x','HH','keeplimits', 'keepticks');
colorbar('EastOutside');
colormap(jet)
set(gca,'Zscale', 'log')
set(gca,'Colorscale', 'log')
set(gca,'Zscale', 'linear')