figure(101)
 semilogy(comb_duration, comb_ABC_MPD02(9,:))
hold on
 semilogy(comb_duration, comb_ABC_MPD03(9,:))
hold off
  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
  title('M2HATS Aerosol Backscatter Coefficient m^{-1} sr^{-1}')
 legend(['MPD02 @' num2str(y(9))], ['MPD03 @' num2str(y(9))])
 grid on
 ylim([10^-9, 10^-4])



 FigH = figure(101);
  FigH.Position = [100 100 1920 300]; % x, y, width, height in pixels
  name=char('backscatter_coeff'); 
  exportgraphics(FigH, [name, '.png'], 'Resolution', 150);