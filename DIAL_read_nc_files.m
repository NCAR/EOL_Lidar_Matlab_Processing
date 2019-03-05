clear all
close all

filename = '/scr/sci/mhayman/DIAL/Processed_Data/RELAMPAGO/wv_dial01.181101.Python.nc';
%filename = '/net/ftp/pub/temp/users/mhayman/LAFE/wv_dial.170809.Python.nc';
%filename = '/scr/eldora1/wvdial_2_processed_data/wv_dial.170814.Python.nc';
%filename = '/net/ftp/pub/temp/users/mhayman/DIAL-PERDIGAO/WVDIAL1_WVDIAL_20170515T0000_20170518T0000_created_20180319__SondeEval.nc';
%filename = '/net/ftp/pub/temp/users/mhayman/DIAL-PERDIGAO/wv_dial.170517.Python.nc';
n = datenum('2017-08-09 00:00:00', 'yyyy-mm-dd HH:MM:SS');
sub_grid = 0;

ncid = netcdf.open(filename, 'NC_NOWRITE');
  %ncdisp(filename, '/', 'min') % use this to display all variables
  ncdisp(filename, 'Absolute_Humidity') 
  
  variable = 'Denoised_Backscatter_Ratio';
  var_units = ncreadatt(filename, variable,'units');
  var_units = erase(var_units, '$');
  var_time = ncread(filename,horzcat('time_',variable));   
  var_range = ncread(filename,horzcat('range_',variable));   
  var = ncread(filename, variable);
  var_variance = ncread(filename, horzcat(variable,'_variance')); 
  var_mask = ncread(filename,horzcat(variable,'_mask'));
  
  variable2 = 'Backscatter_Ratio';
  var2_units = ncreadatt(filename, variable2,'units');
  var2_units = erase(var2_units, '$');
  var2_time = ncread(filename,horzcat('time_',variable2));   
  var2_range = ncread(filename,horzcat('range_',variable2));   
  var2 = ncread(filename, variable2);
  var2_variance = ncread(filename, horzcat(variable2,'_variance')); 
  var2_mask = ncread(filename,horzcat(variable2,'_mask'));
  
  variable3 = 'Absolute_Humidity';
  var3_units = ncreadatt(filename, variable3,'units');
  var3_units = erase(var3_units, '$');
  var3_time = ncread(filename,horzcat('time_',variable3));   
  var3_range = ncread(filename,horzcat('range_',variable3));   
  var3 = ncread(filename, variable3);
  var3_variance = ncread(filename, horzcat(variable3,'_variance')); 
  var3_mask = ncread(filename,horzcat(variable3,'_mask'));
  
  variable4 = 'Denoised_Aerosol_Backscatter_Coefficient';
  var4_units = ncreadatt(filename, variable4,'units');
  var4_units = erase(var4_units, '$');
  var4_time = ncread(filename,horzcat('time_',variable4));   
  var4_range = ncread(filename,horzcat('range_',variable4));   
  var4 = ncread(filename, variable4);
  var4_variance = ncread(filename, horzcat(variable4,'_variance')); 
  var4_mask = ncread(filename,horzcat(variable4,'_mask'));
  
  
  % quick way to generate next day file
  %time(time/3600/24<=1)=inf0;
  
  x = n+double(var_time/3600/24);
  y = double(var_range');
  x2 = n+double(var2_time/3600/24);
  y2 = double(var2_range');
  x3 = n+double(var3_time/3600/24);
  y3 = double(var3_range');
  x4 = n+double(var4_time/3600/24);
  y4 = double(var4_range');
  
  
  
%  if (size(absolute_humidity,2)== size(time,1))
%    x2 = x;
%  else
%    time_AH = ncread(filename,'time_Absolute_Humidity');
%    x2 = n+double(time_AH/3600/24);
%  end
%  if (size(absolute_humidity,1)== size(range,1))
%    y2 = y;
%  else
%    range_AH = ncread(filename,'range_Absolute_Humidity');
%    y2 = double(range_AH');
%  end
netcdf.close(ncid);

scrsz = get(0,'ScreenSize');
date=datestr(n, 'dd mmm yyyy');
plot_size1 = [scrsz(4)/1.5 scrsz(4)/10 scrsz(3)/1.5 scrsz(4)/3];
font_size = 14;
xData =  linspace(fix(min(x)),  ceil(max(x)), 25);

%Zba = real(log10(beta_a_backscat_DN));
%Zbm = real(log10(beta_m_backscat));
%Zcm = real(log10(counts_mol));
%Zcc = real(log10(counts_com));
%Zcon = real(log10(counts_on));
%Zv1 = real(var);
%Zv2 = real(var2);


% if sub_grid == 1
% % subgrid interpolation
%  x_int = interpn(x,3);
%  y_int = interpn(y,3);
%  y2_int = interpn(y2,3);
%  qc_mask_int = interpn(double(beta_a_backscat_mask),3 , 'nearest');
%  qc_mask_ah_int = interpn(double(absolute_humidity_mask),3 , 'nearest');
%  Zba_int = interpn(Zba, 3, 'cubic');
%  Zah_int = interpn(Zah, 3, 'cubic');
% end

% add masking
var(var_mask == 1) = nan;
var2(var2_mask == 1) = nan;
var3(var3_mask == 1) = nan;
var4(var4_mask == 1) = nan;
%if sub_grid == 1
%  Zba_int(qc_mask_int == 1) = nan;
%  Zah_int(qc_mask_ah_int == 1) = nan;
%end
  
% plot variable 1
  figure1 = figure('Position',plot_size1);
  set(gcf,'renderer','zbuffer');
  h = pcolor(x, y/1000, real(log10(var)));
  set(h, 'EdgeColor', 'none'); 
  axis xy; colorbar('EastOutside'); caxis([0 2.5]);
  title({[date ,' ',replace(variable, '_', ' '), ' [',var_units,']']},...
       'fontweight','b','fontsize',font_size)%,'Interpreter', 'none')
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
  axis([n n+1 0 12])
  set(gca, 'XTick',  xData)
  datetick('x','HH','keeplimits', 'keepticks');
  colormap(jet)
%  shading interp
  set(gca,'Fontsize',font_size,'Fontweight','b');

  % plot variable 2
  figure2 = figure('Position',plot_size1);
  set(gcf,'renderer','zbuffer');
  h = pcolor(x2, y2/1000, real(log10(var2)));
  set(h, 'EdgeColor', 'none'); 
  axis xy; colorbar('EastOutside');  caxis([0 3]);
  title({[date ,' ',replace(variable2, '_', ' '), ' [',var2_units,']']},...
       'fontweight','b','fontsize',font_size)%,'Interpreter', 'none')
  ylabel('range (km, AGL)','fontweight','b','fontsize',font_size); 
  axis([n n+1 0 12])
  set(gca, 'XTick',  xData)
  datetick('x','HH','keeplimits', 'keepticks');
  colormap(jet)
 %  shading interp
  set(gca,'Fontsize',font_size,'Fontweight','b');
  
  %cd('/Volumes/documents/WV_DIAL_data/plots/') % point to the directory where data is stored 
  %size = [scrsz(4)/1 scrsz(4)/1 scrsz(3)/0.35 scrsz(4)/2.05]; % use for long plots 
  %FigH = figure(1);
  %set(FigH, 'PaperUnits', 'points', 'PaperPosition', size);
  %name=strcat(date, 'Python'); 
  %print(FigH, name, '-dpng', '-r300') % set the resolution as 300 dpiFigH = figure(1);
  
  %cd('/Volumes/documents/WV_DIAL_data/processed_data/') % point to the directory where data is stored 
  %name=strcat(date,'Python');
  %N_avg = Zah';
  %N_error = absolute_humidity_variance';
  %time_new = x;
  %range = range';
  %save(name, 'N_avg', 'range', 'time_new', 'N_error')

  
%% plot the molecular counts
%  figure3 = figure('Position',plot_size1);
%  set(gcf,'renderer','zbuffer');
%  h = pcolor(x, y, Zcm);
%  set(h, 'EdgeColor', 'none'); 
%  axis xy; colorbar('EastOutside'); caxis([5 12]);
%  title({[date,' log10 Molecular Counts']},...
%       'fontweight','b','fontsize',font_size)
%  ylabel('range (m)','fontweight','b','fontsize',font_size); 
%  axis([n n+1 0 12000])
%  set(gca, 'XTick',  xData)
%  datetick('x','HH','keeplimits', 'keepticks');
%  colormap(jet)
% %  shading interp
%  set(gca,'Fontsize',font_size,'Fontweight','b');

% % plot the combined backscatter counts
%  figure4 = figure('Position',plot_size1);
%  set(gcf,'renderer','zbuffer');
%  h = pcolor(x, y, Zcc);
%  set(h, 'EdgeColor', 'none'); 
%  axis xy; colorbar('EastOutside'); caxis([5 12]);
%  title({[date,' log10 Combined Counts']},...
%       'fontweight','b','fontsize',font_size)
%  ylabel('range (m)','fontweight','b','fontsize',font_size); 
%  axis([n n+1 0 12000])
%  set(gca, 'XTick',  xData)
%  datetick('x','HH','keeplimits', 'keepticks');
%  colormap(jet)
% %  shading interp
%  set(gca,'Fontsize',font_size,'Fontweight','b');
  
%   % plot the denoised offline backscatter counts
%  figure5 = figure('Position',plot_size1);
%  set(gcf,'renderer','zbuffer');
%  h = pcolor(x2, y, Zcoff);
%  set(h, 'EdgeColor', 'none'); 
%  axis xy; colorbar('EastOutside'); caxis([0 5]);
%  title({[date,' log10 Denoised Offline Counts']},...
%       'fontweight','b','fontsize',font_size)
%  ylabel('range (m)','fontweight','b','fontsize',font_size); 
%  axis([n n+1 0 12000])
%  set(gca, 'XTick',  xData)
%  datetick('x','HH','keeplimits', 'keepticks');
%  colormap(jet)
% %  shading interp
%  set(gca,'Fontsize',font_size,'Fontweight','b');
  
%  % plot the denoised online backscatter counts
%  figure6 = figure('Position',plot_size1);
%  set(gcf,'renderer','zbuffer');
%  h = pcolor(x2, y, Zcon);
%  set(h, 'EdgeColor', 'none'); 
%  axis xy; colorbar('EastOutside'); caxis([0 5]);
%  title({[date,' log10 Denoised Online Counts']},...
%       'fontweight','b','fontsize',font_size)
%  ylabel('range (m)','fontweight','b','fontsize',font_size); 
%  axis([n n+1 0 12000])
%  set(gca, 'XTick',  xData)
%  datetick('x','HH','keeplimits', 'keepticks');
%  colormap(jet)
% %  shading interp
%  set(gca,'Fontsize',font_size,'Fontweight','b');
  
  %plot combined backscatter and humidity
  figure10 = figure('visible', 'off','Position',[scrsz(4)/1.5 scrsz(4)/10 scrsz(3)/1.5 scrsz(4)/1.5]);
  set(figure10, 'visible', 'on', 'PaperUnits', 'points', 'PaperPosition', [0 0 1280 800]);
  subplot1=subplot(2,1,2,'Parent',figure10);
  box(subplot1,'on');
  set(gcf,'renderer','zbuffer');
  h = pcolor(x4,y4/1000,real(log10(var4)));
  set(h, 'EdgeColor', 'none');
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  set(gca, 'XTick',  xData)
  colorbar('EastOutside');
  axis([fix(min(x)) fix(min(x))+1 0 12])
  caxis([-8 -4]);
  %caxis([1 5]);
  datetick('x','HH','keeplimits', 'keepticks');
  colormap(jet)
  %shading interp 
  hh =  title({[date ,' ',replace(variable4, '_', ' '), ' [',var4_units,']']},...
       'fontweight','b','fontsize',font_size)%,'Interpreter', 'none')
  P_t = get(hh, 'Position');
  set(hh,'Position', [P_t(1) P_t(2)+0.2 P_t(3)])
  xlabel('Time (UTC)','fontweight','b','fontsize',font_size);
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size);
  set(gca,'Fontsize',font_size,'Fontweight','b');
  
  %plot water vapor in g/m^3
  subplot1=subplot(2,1,1,'Parent',figure10);
  box(subplot1,'on'); %(number density in mol/cm3)(1e6 cm3/m3)/(N_A mol/mole)*(18g/mole)
  set(gcf,'renderer','zbuffer');
  h = pcolor(x3,y3/1000,var3);
  set(h, 'EdgeColor', 'none');
  set(gca, 'XTick',  xData)
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  colorbar('EastOutside');
  axis([fix(min(x)) fix(min(x))+1 0 6])
  caxis([0 20]);
  datetick('x','HH','keeplimits', 'keepticks');
  colormap(jet)
  %shading interp
  hh =  title({[date ,' ',replace(variable3, '_', ' '), ' [',var3_units,']']},...
       'fontweight','b','fontsize',font_size)%,'Interpreter', 'none')
  P_t = get(hh, 'Position');
  set(hh,'Position', [P_t(1) P_t(2)+0.2 P_t(3)])
  xlabel('Time (UTC)','fontweight','b','fontsize',font_size); 
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
  set(gca,'Fontsize',font_size,'Fontweight','b');
 

