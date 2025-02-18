figure(1)
caxis([1E-8 1E-4])
colormap(jet)



FigH = figure(1);
  %set(gca,'Fontsize',16,'Fontweight','b');  
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920 275*3]);
  name=strcat('NY_MPD_HSRL_jet'); 
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution 


  FigH = figure(3);
  %set(gca,'Fontsize',16,'Fontweight','b');  
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920 275*3]);
  name=strcat('NY_MPD_T'); 
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution 


  FigH = figure(2);
  %set(gca,'Fontsize',16,'Fontweight','b');  
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920 275*3]);
  name=strcat('NY_MPD_WV'); 
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution 