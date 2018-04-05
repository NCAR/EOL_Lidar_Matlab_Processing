clear all
close all

%filename = '/net/ftp/pub/temp/users/mhayman/DIAL-LAFE/wv_dial.170808.Python.nc';
filename = '/net/ftp/pub/temp/users/mhayman/DIAL-PERDIGAO/WVDIAL1_WVDIAL_20170515T0000_20170518T0000_created_20180319__SondeEval.nc';
%filename = '/net/ftp/pub/temp/users/mhayman/DIAL-PERDIGAO/wv_dial.170517.Python.nc';
n = datenum('2017-05-15 00:00:00', 'yyyy-mm-dd HH:MM:SS');
sub_grid = 0;

ncid = netcdf.open(filename, 'NC_NOWRITE');
  ncdisp(filename) 
 % time = ncread(filename,'time');
 % range = ncread(filename,'range');   
  
 % beta_a_backscat = ncread(filename,'Aerosol_Backscatter_Coefficient');
 % beta_a_backscat_mask = ncread(filename,'Aerosol_Backscatter_Coefficient_mask');
 % beta_a_backscat_DN = ncread(filename,'Denoised_Aerosol_Backscatter_Coefficient');
 % beta_a_backscat_DN_mask = ncread(filename,'Denoised_Aerosol_Backscatter_Coefficient_mask');
 % beta_m_backscat = ncread(filename,'Molecular_Backscatter_Coefficient');
 % counts_mol = ncread(filename,'Molecular_Backscatter_Channel');
 % counts_mol_mask = ncread(filename,'Molecular_Backscatter_Channel_mask');
 % counts_com = ncread(filename,'Total_Backscatter_Channel');
 % counts_com_mask = ncread(filename,'Total_Backscatter_Channel_mask');
 % counts_on = ncread(filename,'Denoised_Online_Backscatter_Channel');
 % counts_on_mask = ncread(filename,'Denoised_Online_Backscatter_Channel_mask');
 % counts_off = ncread(filename,'Denoised_Offline_Backscatter_Channel');
 % counts_off_mask = ncread(filename,'Denoised_Offline_Backscatter_Channel_mask');
  
  time = ncread(filename,'time_Absolute_Humidity');   
  range = ncread(filename,'range_Absolute_Humidity');   
  absolute_humidity = ncread(filename,'Absolute_Humidity');
  absolute_humidity_variance = ncread(filename,'Absolute_Humidity_variance'); 
  absolute_humidity_mask = ncread(filename,'Absolute_Humidity_mask');
  
  time2 = ncread(filename,'time_Offline_Backscatter_Channel');   
  range2 = ncread(filename,'range_Offline_Backscatter_Channel');   
  counts_off = ncread(filename,'Offline_Backscatter_Channel');
  counts_off_variance = ncread(filename,'Offline_Backscatter_Channel_variance'); 
  counts_off_mask = ncread(filename,'Offline_Backscatter_Channel_mask');
  
  % quick way to generate next day file
  %time(time/3600/24<=1)=inf0;
  
  x = n+double(time/3600/24);
  y = double(range');
  x2 = n+double(time2/3600/24);
  y2 = double(range2');
  

  
  
  
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

%Zba = real(log10(beta_a_backscat_DN));
%Zbm = real(log10(beta_m_backscat));
%Zcm = real(log10(counts_mol));
%Zcc = real(log10(counts_com));
%Zcon = real(log10(counts_on));
Zcoff = real(log10(counts_off));
Zah = real(absolute_humidity);
xData =  linspace(fix(min(x)),  ceil(max(x)), 25);

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
%Zba(beta_a_backscat_DN_mask == 1) = nan;
Zah(absolute_humidity_mask == 1) = nan;
Zcoff(counts_off_mask == 1) = nan;
%Zcm(counts_mol_mask == 1) = nan;
%Zcc(counts_com_mask == 1) = nan;
%if sub_grid == 1
%  Zba_int(qc_mask_int == 1) = nan;
%  Zah_int(qc_mask_ah_int == 1) = nan;
%end
  
% % plot the beta_a_backscat
% figure1 = figure('Position',plot_size1);
%  set(gcf,'renderer','zbuffer');
%  if sub_grid == 1 
%     h = pcolor(x_int, y_int, Zba_int);
%  else
%     h = pcolor(x, y, Zba);
%  end
%  set(h, 'EdgeColor', 'none'); 
%  axis xy; colorbar('EastOutside'); caxis([-8 -3]);
%  title({[date,' log10 \beta_{a}backscatter']},...
%       'fontweight','b','fontsize',font_size)
%  ylabel('range (m)','fontweight','b','fontsize',font_size); 
%  axis([n n+1 0 12000])
%  set(gca, 'XTick',  xData)
%  datetick('x','HH','keeplimits', 'keepticks');
%  colormap(jet)
% %  shading interp
%  set(gca,'Fontsize',font_size,'Fontweight','b');

% plot the absoulte humdity
  figure2 = figure('Position',plot_size1);
  set(gcf,'renderer','zbuffer');
  %if sub_grid == 1 
  %   h = pcolor(x_int, y2_int, Zah_int);
  %else
     h = pcolor(x, y/1000, Zah);
  %end
  set(h, 'EdgeColor', 'none'); 
  axis xy; colorbar('EastOutside'); caxis([0 12]);
  title({[date,' Absolute Humidity']},...
       'fontweight','b','fontsize',font_size)
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
  axis([n n+3 0 6])
  set(gca, 'XTick',  xData)
  datetick('x','HH','keeplimits', 'keepticks');
  colormap(jet)
%  shading interp
  set(gca,'Fontsize',font_size,'Fontweight','b');

  cd('/Volumes/documents/WV_DIAL_data/plots/') % point to the directory where data is stored 
  size = [scrsz(4)/1 scrsz(4)/1 scrsz(3)/0.35 scrsz(4)/2.05]; % use for long plots 
  FigH = figure(1);
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', size);
  name=strcat(date, 'Python'); 
  print(FigH, name, '-dpng', '-r300') % set the resolution as 300 dpiFigH = figure(1);
  
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
  figure1 = figure('visible', 'off','Position',[scrsz(4)/1.5 scrsz(4)/10 scrsz(3)/1.5 scrsz(4)/1.5]);
  set(figure1, 'visible', 'on', 'PaperUnits', 'points', 'PaperPosition', [0 0 1280 800]);
  subplot1=subplot(2,1,2,'Parent',figure1);
  box(subplot1,'on');
  set(gcf,'renderer','zbuffer');
  h = pcolor(x2,y2/1000,Zcoff);
  set(h, 'EdgeColor', 'none');
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  set(gca, 'XTick',  xData)
  colorbar('EastOutside');
  axis([fix(min(x)) fix(min(x))+1 0 12])
 % caxis([-8 -4]);
  caxis([1 5]);
  datetick('x','HH','keeplimits', 'keepticks');
  colormap(jet)
  %shading interp 
  hh = title({[date,'  Denoised log10 \beta_{a} m^{-1}sr^{-1}']},'fontweight','b','fontsize',font_size);
  P_t = get(hh, 'Position');
  set(hh,'Position', [P_t(1) P_t(2)+0.2 P_t(3)])
  xlabel('Time (UTC)','fontweight','b','fontsize',font_size);
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size);
  set(gca,'Fontsize',font_size,'Fontweight','b');
  
  %plot water vapor in g/m^3
  subplot1=subplot(2,1,1,'Parent',figure1);
  box(subplot1,'on'); %(number density in mol/cm3)(1e6 cm3/m3)/(N_A mol/mole)*(18g/mole)
  set(gcf,'renderer','zbuffer');
  h = pcolor(x,y/1000,Zah);
  set(h, 'EdgeColor', 'none');
  set(gca, 'XTick',  xData)
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  colorbar('EastOutside');
  axis([fix(min(x)) fix(min(x))+1 0 6])
  caxis([0 12]);
  datetick('x','HH','keeplimits', 'keepticks');
  colormap(jet)
  %shading interp
  hh = title({[date,'  Water Vapor (g/m^{3})']},'fontweight','b','fontsize',font_size);
  P_t = get(hh, 'Position');
  set(hh,'Position', [P_t(1) P_t(2)+0.2 P_t(3)])
  xlabel('Time (UTC)','fontweight','b','fontsize',font_size); 
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
  set(gca,'Fontsize',font_size,'Fontweight','b');
 

