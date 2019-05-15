clear all
close all

filename = '/scr/sci/mhayman/DIAL/Processed_Data/RELAMPAGO/wv_dial01.181108.Python.nc';
%filename = '/net/ftp/pub/temp/users/mhayman/LAFE/wv_dial.170809.Python.nc';
%filename = '/scr/eldora1/wvdial_2_processed_data/wv_dial.170814.Python.nc';
%filename = '/net/ftp/pub/temp/users/mhayman/DIAL-PERDIGAO/WVDIAL1_WVDIAL_20170515T0000_20170518T0000_created_20180319__SondeEval.nc';
%filename = '/net/ftp/pub/temp/users/mhayman/DIAL-PERDIGAO/wv_dial.170517.Python.nc';
date = filename(end-15:end-10);
n = datenum(date, 'yymmdd');
%n = datenum('2017-08-09 00:00:00', 'yyyy-mm-dd HH:MM:SS');

ncid = netcdf.open(filename, 'NC_NOWRITE');
  %ncdisp(filename, '/', 'min') % use this to display all variables
  %ncdisp(filename, 'Absolute_Humidity') 
  

  variable{1} = 'Absolute_Humidity';
  variable{2} = 'Attenuated_Backscatter';
 % variable{2} = 'WV_Offline_Backscatter_Channel';
 % variable{3} = 'Denoised_Aerosol_Backscatter_Coefficient'; 
 % variable{4} = 'Denoised_Backscatter_Ratio';
  
  for i = 1:size(variable,2)
   var_units{i} = ncreadatt(filename, variable{i},'units');
   var_units{i} = erase(var_units{i}, '$');
   var_time{i} = ncread(filename,horzcat('time_',variable{i}));   
   var_range{i} = ncread(filename,horzcat('range_',variable{i}));   
   var{i} = ncread(filename, variable{i});
   var_variance{i} = ncread(filename, horzcat(variable{i},'_variance')); 
   var_mask{i} = ncread(filename,horzcat(variable{i},'_mask')); 
   x{i} = n+double(var_time{i}/3600/24);
   y{i} = double(var_range{i}');
  end  

netcdf.close(ncid);

scrsz = get(0,'ScreenSize');
date=datestr(n, 'dd mmm yyyy');
plot_size1 = [scrsz(4)/1.5 scrsz(4)/10 scrsz(3)/1.5 scrsz(4)/3];
font_size = 14;
xData =  linspace(fix(min(x{1})),  ceil(max(x{1})), 25);

% add masking
  for i = 1:size(variable,2)
     var{i}(var_mask{i} == 1) = nan;
  end

 for i = 1:size(variable,2) 
% plot variable 1
  figure1 = figure('Position',plot_size1);
  set(gcf,'renderer','zbuffer');
  h = pcolor(x{i}, y{i}/1000, real(log10(var{i})));
  set(h, 'EdgeColor', 'none'); 
  axis xy; colorbar('EastOutside'); %caxis([0 10]);
  title({[date ,' ',replace(variable{i}, '_', ' '), ' [',var_units{i},']']},...
       'fontweight','b','fontsize',font_size)%,'Interpreter', 'none')
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
  axis([n n+1 0 12])
  set(gca, 'XTick',  xData)
  datetick('x','HH','keeplimits', 'keepticks');
  colormap(jet)
%  shading interp
  set(gca,'Fontsize',font_size,'Fontweight','b');
 end
%caxis([0 2.5]);
 
  %plot combined backscatter and humidity
  figure10 = figure('visible', 'off','Position',[scrsz(4)/1.5 scrsz(4)/10 scrsz(3)/1.5 scrsz(4)/1.5]);
  set(figure10, 'visible', 'on', 'PaperUnits', 'points', 'PaperPosition', [0 0 1280 800]);
  subplot1=subplot(2,1,2,'Parent',figure10);
  box(subplot1,'on');
  set(gcf,'renderer','zbuffer');
  h = pcolor(x{1},y{1}/1000, var{1}); %real(log10(var{1})));
  set(h, 'EdgeColor', 'none');
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  set(gca, 'XTick',  xData)
  colorbar('EastOutside');
  axis([fix(min(x{1})) fix(min(x{1}))+1 0 6])
  caxis([1 20]);
  %caxis([1 5]);
  datetick('x','HH','keeplimits', 'keepticks');
  colormap(jet)
  %shading interp 
  hh =  title({[date ,' ',replace(variable{1}, '_', ' '), ' [',var_units{1},']']},...
       'fontweight','b','fontsize',font_size)%,'Interpreter', 'none')
  P_t = get(hh, 'Position');
  set(hh,'Position', [P_t(1) P_t(2)+0.2 P_t(3)])
  xlabel('Time (UTC)','fontweight','b','fontsize',font_size);
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size);
  set(gca,'Fontsize',font_size,'Fontweight','b');

set(gca,'Zscale', 'linear')
set(gca,'Colorscale', 'linear')
set(gca,'Zscale', 'linear')

  subplot1=subplot(2,1,1,'Parent',figure10);
  box(subplot1,'on'); %(number density in mol/cm3)(1e6 cm3/m3)/(N_A mol/mole)*(18g/mole)
  set(gcf,'renderer','zbuffer');
  h = pcolor(x{2},y{2}/1000,var{2});
  set(h, 'EdgeColor', 'none');
  set(gca, 'XTick',  xData)
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  colorbar('EastOutside');
  axis([fix(min(x{2})) fix(min(x{2}))+1 0 12])
  caxis([1e3 1e9]);
  datetick('x','HH','keeplimits', 'keepticks');
  colormap(jet)
  %shading interp
  hh =  title({[date ,' ',replace(variable{2}, '_', ' '), ' [',var_units{2},']']},...
       'fontweight','b','fontsize',font_size)%,'Interpreter', 'none')
  P_t = get(hh, 'Position');
  set(hh,'Position', [P_t(1) P_t(2)+0.2 P_t(3)])
  xlabel('Time (UTC)','fontweight','b','fontsize',font_size); 
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
  set(gca,'Fontsize',font_size,'Fontweight','b');
 
set(gca,'Zscale', 'log')
set(gca,'Colorscale', 'log')
set(gca,'Zscale', 'linear')
  

  cd('/Volumes/documents/WV_DIAL_data/plots/') % point to the directory where data is stored 
  %size = [scrsz(4)/1 scrsz(4)/1 scrsz(3)/0.35 scrsz(4)/2.05]; % use for long plots 
  size = [scrsz(4)/1.5 scrsz(4)/10 scrsz(3)/1.5 scrsz(4)/1.5];
  FigH = figure(3);
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', size);
  name=strcat(date, 'Python'); 
  print(FigH, name, '-dpng', '-r0') % set the resolution as 300 dpiFigH = figure(1);
  
  %cd('/Volumes/documents/WV_DIAL_data/processed_data/') % point to the directory where data is stored 
  %name=strcat(date,'Python');
  %N_avg = Zah';
  %N_error = absolute_humidity_variance';
  %time_new = x;
  %range = range';
  %save(name, 'N_avg', 'range', 'time_new', 'N_error')
