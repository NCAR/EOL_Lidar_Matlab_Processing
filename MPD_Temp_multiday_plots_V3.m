clear all; close all;
tic
dd=pwd;

 node = 'MPD03';
 date = '20230719';  
 days = 52; skip = 4;
%  date = '20230728';  
%  days = 15; skip = 3;  
  
lapse_rate = 0.0065; %standard atmosphere lapse rate
T_lapse = 0.0065; %standard atmosphere lapse rate
     
%serv_path = '/Volumes/documents/MPD';
serv_path = '/Volumes/fog1/rsfdata/MPD';
%serv_path = '/Volumes/eol/fog1/rsfdata/MPD';
%plot_path = '/Volumes/documents/MPD/Plots/';
plot_path = '/Users/spuler/Desktop/mpd/Plots/';
%flag.near = 0;  % read in the near range channel (0=off 1=on)
font_size = 36; % use this for 2018a version

%replot time vs range images at start of processing (0=off 1=on)
flag.replot = 1;
%save figures at end of processing (0=off 1=on)
flag.save_figs = 1;
%save data at end of processing (0=off 1=on)
flag.save_data = 0;
%plot NetCDF sonde data on top (0=off 1=on)
flag.plot_sonde_data = 0;
%decimate figures to the screen 2x size (1900 pixels x2)
flag.decimate = 1;
  
%C = importdata('NCAR_C_Map.mat');


if strcmp(node,'MPD01')==1
  % cd(strcat(serv_path,'/mpd_01_processed_data//Matlab_temp'))
   cd(strcat(serv_path,'/mpd_01_processed_data/Quickload/ReProcessing'))
   path_node = 'mpd01.';
 elseif strcmp(node,'MPD02')==1
   %cd(strcat(serv_path,'/mpd_02_processed_data/Quickload/ReProcessing'))
   cd(strcat(serv_path,'/mpd_02_processed_data/Quickload/FullProcessingBootstrap'))
   path_node = 'mpd02.';
 elseif strcmp(node,'MPD03')==1
   %cd(strcat(serv_path,'/mpd_03_processed_data/Quickload/ReProcessing'))
   cd(strcat(serv_path,'/mpd_03_processed_data/Quickload/FullProcessingBootstrap'))
   path_node = 'mpd03.';
 elseif strcmp(node,'MPD04')==1
   cd(strcat(serv_path,'/mpd_04_processed_data/Quickload/ReProcessing'))
   path_node = 'mpd04.';
 elseif strcmp(node,'MPD05')==1
   %cd(strcat(serv_path,'/mpd_05_processed_data/Quickload/ReProcessing'))
   cd(strcat(serv_path,'/mpd_05_processed_data/Quickload/FullProcessingBootstrap'))
   path_node = 'mpd05.';
end

i=1;


for i=1:days
  if i==1  
    if exist(strcat(path_node, date, '.Matlab.mat'))==2
      load(strcat(path_node, date, '.Matlab.mat'),'Retrievals')  
      % load(strcat(path_node, date, '.Matlab.mat'),'Retrievals','Data')  
      % reading in the 'Data' variable slows down the read considerably
      % save the file as save(strcat('mpd05.', date, '.Matlab.mat'), '-struct', 'Retreivals','Data');
      % then can load with the following line more quickly  
    end
    Temp_comb = Retrievals.Temperature.Value';
    Temp_comb_avg = Retrievals.Temperature.Smoothed';
    Temp_comb_var = Retrievals.Temperature.VarianceSm';
%    Alpha0_comb = Retrievals.PerturbOrders.Order0.Alpha.O2Online';
%    Alpha1_comb = Retrievals.PerturbOrders.Order1.Alpha.O2Online';
%    Alpha2_comb = Retrievals.PerturbOrders.Order2.Alpha.O2Online';
    WV_mask_comb = Retrievals.WaterVapor.Mask';
%    Alpha_comb = Retrievals.Alpha';
%    O2_online_comb = Data.Lidar.Interp.O2OnlineComb.Data';
%    BSR_comb = Retrievals.Python.BackRatio.Value';
%    BSR_comb = Retrievals.Python.MPD.BackRatio.Value';
%    BSR_duration = (datenum(date, 'yyyymmdd') + Retrievals.Python.MPD.BackRatio.TimeStamp/3600/24);
    duration =  (datenum(date, 'yyyymmdd') + Retrievals.Temperature.TimeStamp/3600/24)';
    range = Retrievals.Temperature.Range;
    % following four lines are to build a lapse rate temp field based on surface station data
%      T_surf = Data.Temperature;
%      WS_time = Data.TimeStamp/24;
%      T_surf_grid = interp1( WS_time, T_surf, Retrievals.Temperature.TimeStamp/3600/24, 'linear')';
%      T_lapse =  T_surf_grid-lapse_rate*range;
  else   
    date = datestr(addtodate(datenum(date, 'yyyymmdd'), 1, 'day'), 'yyyymmdd');
    if exist(strcat(path_node, date, '.Matlab.mat'))==2
     % load(strcat(path_node, date, '.Matlab.mat'),'Retrievals', 'Data')
      load(strcat(path_node, date, '.Matlab.mat'),'Retrievals')
    end
    % following four lines are to build a lapse rate temp field based on surface station data
%      T_surf = Data.Temperature;
%      WS_time = Data.TimeStamp/24;
%      T_surf_grid = interp1( WS_time, T_surf, Retrievals.Temperature.TimeStamp/3600/24, 'linear')';
%      Temp_lapse =  T_surf_grid-lapse_rate*range;
      
%    T_lapse = vertcat(T_lapse, Temp_lapse);
    Temp_comb = vertcat(Temp_comb, Retrievals.Temperature.Value');
    Temp_comb_avg = vertcat(Temp_comb_avg, Retrievals.Temperature.Smoothed');
    Temp_comb_var = vertcat(Temp_comb_var, Retrievals.Temperature.VarianceSm');
 %   Alpha0_comb = vertcat(Alpha0_comb, Retrievals.PerturbOrders.Order0.Alpha.O2Online');
 %   Alpha1_comb = vertcat(Alpha1_comb, Retrievals.PerturbOrders.Order1.Alpha.O2Online');
 %   Alpha2_comb = vertcat(Alpha2_comb, Retrievals.PerturbOrders.Order2.Alpha.O2Online');
    WV_mask_comb =  vertcat(WV_mask_comb, Retrievals.WaterVapor.Mask'); 
 %   Alpha_comb = vertcat(Alpha_comb, Retrievals.Alpha');
 %   BSR_comb = vertcat(BSR_comb, Retrievals.Python.MPD.BackRatio.Value');   
 %   BSR_duration = vertcat(BSR_duration, (datenum(date, 'yyyymmdd') + Retrievals.Python.MPD.BackRatio.TimeStamp/3600/24));
    duration = vertcat(duration, (datenum(date, 'yyyymmdd') + Retrievals.Temperature.TimeStamp/3600/24)');
  end
end


%Alpha0_var = movvar(Alpha0_comb, [2 2]);
% apply preliminary mask based on O_2 order absorption zeroth order variance 
Temp_comb_avg(Temp_comb_var > 5^2)= nan;
%     Temp_comb(Alpha0_var>1e-7)= nan;
%    Alpha_comb(Alpha0_var>1e-7)= nan;
 % Alpha0_comb(Alpha0_var>1e-7)= nan;
 % Alpha1_comb(Alpha0_var>1e-7)= nan;
 % Alpha2_comb(Alpha0_var>1e-7)= nan;
%  Temp_comb_avg(Alpha0_comb>7.5e-4 | Alpha0_comb<-2.5e-4)= nan;
%      Temp_comb(Alpha0_comb>7.5e-4 | Alpha0_comb<-2.5e-4)= nan;
%     Alpha_comb(Alpha0_comb>7.5e-4 | Alpha0_comb<-2.5e-4)= nan; 
 %  Alpha0_comb(Alpha0_comb>7.5e-4 | Alpha0_comb<-2.5e-4)= nan;
 %  Alpha1_comb(Alpha0_comb>7.5e-4 | Alpha0_comb<-2.5e-4)= nan;
 %  Alpha2_comb(Alpha0_comb>7.5e-4 | Alpha0_comb<-2.5e-4)= nan; 

    
  if days == 1
   xData =  linspace(fix(min(duration)),  ceil(max(duration)), 25);
 else
   xData =  linspace( fix(min(duration)),  ceil(max(duration)), round((ceil(max(duration))-fix(min(duration)))/skip)+1 );
  % xData =  linspace(fix(min(duration)),  round(max(duration)), 36);
  end

scrsz = [1  1  1920 1200];
Scrnsize=[scrsz(4)/2 scrsz(4)/10 scrsz(3)/1 scrsz(4)/2];


x = duration';
y = range./1000;
 
 
 %Z = Temp_comb_avg'-273.15; 
 Z = Temp_comb_avg'; 
 %Z(isnan(Alpha_comb)) = nan;
 %Z = Alpha_comb'; 
 figure('Position',Scrnsize)
 set(gcf,'renderer','zbuffer');
 h = pcolor(x,y,Z);
 set(h, 'EdgeColor', 'none');
 colorbar('EastOutside');
 axis([fix(min(x)) ceil(max(x)) 0 4])
 %caxis([1e-4 5e-4]);
 caxis([260 320]);
 caxis([273-10 273+35]);
 colormap(jet)
 %colormap(parula)
 %shading interp
 ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
 set(gca, 'XTick',  xData)
 set(gca,'TickDir','out');
 set(gca,'TickLength',[0.005; 0.0025]);
 if days == 1
   datetick('x','HH:MM','keeplimits', 'keepticks');
   xlabel('Time (UTC)','fontweight','b','fontsize',font_size);
   hh = title({[date, node, ' temperature (K)']},'fontweight','b','fontsize',font_size);
 else
   datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
   hh = title({[node, ' Temperature (K)']},'fontweight','b','fontsize',font_size);
   %hh = title({'DIAL Water Vapor (g m^{-3})'},'fontweight','b','fontsize',font_size);
 end
 set(gca,'Fontsize',font_size,'Fontweight','b');
 
 Z = Temp_comb'; 
 figure('Position',Scrnsize)
 set(gcf,'renderer','zbuffer');
 h = pcolor(x,y,Z);
 set(h, 'EdgeColor', 'none');
 colorbar('EastOutside');
 axis([fix(min(x)) ceil(max(x)) 0 4])
 %caxis([240 320]);
 caxis([-5+273.15 30+273.15]);
 caxis([240 320]);
 %caxis([-5 30]); % Celsius
 colormap(jet)
 %colormap(parula)
 %shading interp
 ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
 set(gca, 'XTick',  xData)
 set(gca,'TickDir','out');
 set(gca,'TickLength',[0.005; 0.0025]);
 if days == 1
   datetick('x','HH:MM','keeplimits', 'keepticks');
   xlabel('Time (UTC)','fontweight','b','fontsize',font_size);
   hh = title({[date, node, ' temperature (K)']},'fontweight','b','fontsize',font_size);
 else
   datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
   hh = title({[node, ' Temperature (K)']},'fontweight','b','fontsize',font_size);
   %hh = title({'DIAL Water Vapor (g m^{-3})'},'fontweight','b','fontsize',font_size);
 end
 set(gca,'Fontsize',font_size,'Fontweight','b');
 
 
%  Z = Alpha_comb';
%  figure('Position',Scrnsize)
%  set(gcf,'renderer','zbuffer');
%  h = pcolor(x,y,Z);
%  set(h, 'EdgeColor', 'none');
%  colorbar('EastOutside');
%  axis([fix(min(x)) ceil(max(x)) 0 6])
%  caxis([1e-6 3e-4]);
%  %colormap(jet)
%  colormap(parula)
%  %shading interp
%  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
%  set(gca, 'XTick',  xData)
%  set(gca,'TickDir','out');
%  set(gca,'TickLength',[0.005; 0.0025]);
%  if days == 1
%    datetick('x','HH:MM','keeplimits', 'keepticks');
%    xlabel('Time (UTC)','fontweight','b','fontsize',font_size);
%    hh = title({[date, node, ' Alpha (m^{-1})']},'fontweight','b','fontsize',font_size);
%  else
%    datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
%    hh = title({[node, ' Alpha (m^{-1})']},'fontweight','b','fontsize',font_size);
%    %hh = title({'DIAL Water Vapor (g m^{-3})'},'fontweight','b','fontsize',font_size);
%  end
%  set(gca,'Fontsize',font_size,'Fontweight','b');
%  
%  
%  Z = Alpha0_comb';
%  figure('Position',Scrnsize)
%  set(gcf,'renderer','zbuffer');
%  h = pcolor(x,y,Z);
%  set(h, 'EdgeColor', 'none');
%  colorbar('EastOutside');
%  axis([fix(min(x)) ceil(max(x)) 0 6])
%  caxis([1e-6 3e-4]);
%  %colormap(jet)
%  colormap(parula)
%  %shading interp
%  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
%  set(gca, 'XTick',  xData)
%  set(gca,'TickDir','out');
%  set(gca,'TickLength',[0.005; 0.0025]);
%  if days == 1
%    datetick('x','HH:MM','keeplimits', 'keepticks');
%    xlabel('Time (UTC)','fontweight','b','fontsize',font_size);
%    hh = title({[date, node, ' Alpha 0^{th} (m^{-1})']},'fontweight','b','fontsize',font_size);
%  else
%    datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
%    hh = title({[node, ' Alpha 0^{th} (m^{-1})']},'fontweight','b','fontsize',font_size);
%    %hh = title({'DIAL Water Vapor (g m^{-3})'},'fontweight','b','fontsize',font_size);
%  end
%  set(gca,'Fontsize',font_size,'Fontweight','b');
% 
% 
% Z = Alpha1_comb';  
%  figure('Position',Scrnsize)
%  set(gcf,'renderer','zbuffer');
%  h = pcolor(x,y,Z);
%  set(h, 'EdgeColor', 'none');
%  colorbar('EastOutside');
%  axis([fix(min(x)) ceil(max(x)) 0 6])
%  caxis([1e-6 4e-5]);
%  %colormap(jet)
%  colormap(parula)
%  %shading interp
%  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
%  set(gca, 'XTick',  xData)
%  set(gca,'TickDir','out');
%  set(gca,'TickLength',[0.005; 0.0025]);
%  if days == 1
%    datetick('x','HH:MM','keeplimits', 'keepticks');
%    xlabel('Time (UTC)','fontweight','b','fontsize',font_size);
%    hh = title({[date, node, ' Alpha 1^{st} (m^{-1})']},'fontweight','b','fontsize',font_size);
%  else
%    datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
%    hh = title({[node, ' Alpha 1^{st} (m^{-1})']},'fontweight','b','fontsize',font_size);
%    %hh = title({'DIAL Water Vapor (g m^{-3})'},'fontweight','b','fontsize',font_size);
%  end
%  set(gca,'Fontsize',font_size,'Fontweight','b');
% 
% 
% Z = Alpha2_comb';  
%  figure('Position',Scrnsize)
%  set(gcf,'renderer','zbuffer');
%  h = pcolor(x,y,Z);
%  set(h, 'EdgeColor', 'none');
%  colorbar('EastOutside');
%  axis([fix(min(x)) ceil(max(x)) 0 6])
%  caxis([1e-6 4e-6]);
%  %colormap(jet)
%  colormap(parula)
%  %shading interp
%  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
%  set(gca, 'XTick',  xData)
%  set(gca,'TickDir','out');
%  set(gca,'TickLength',[0.005; 0.0025]);
%  if days == 1
%    datetick('x','HH:MM','keeplimits', 'keepticks');
%    xlabel('Time (UTC)','fontweight','b','fontsize',font_size);
%    hh = title({[date, node, ' Alpha 2^{nd} (m^{-1})']},'fontweight','b','fontsize',font_size);
%  else
%    datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
%    hh = title({[node, ' Alpha 2^{nd} (m^{-1})']},'fontweight','b','fontsize',font_size);
%    %hh = title({'DIAL Water Vapor (g m^{-3})'},'fontweight','b','fontsize',font_size);
%  end
%  set(gca,'Fontsize',font_size,'Fontweight','b');

%  x = BSR_duration';
%  y = Retrievals.Python.MPD.BackRatio.Range./1000;  % note this is Matt's version and not mine
%  Z = BSR_comb';  
%  figure('Position',Scrnsize)
%  set(gcf,'renderer','zbuffer');
%  h = pcolor(x,y,Z);
%  set(h, 'EdgeColor', 'none');
%  colorbar('EastOutside');
%  axis([fix(min(x)) ceil(max(x)) 0 6])
%  caxis([1 1000]);
%  %caxis([1e-1 1e3]);
%  colormap(jet)
%  %shading interp
%  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
%  set(gca, 'XTick',  xData)
%  set(gca,'TickDir','out');
%  set(gca,'TickLength',[0.005; 0.0025]);
%  if days == 1
%    datetick('x','HH:MM','keeplimits', 'keepticks');
%    xlabel('Time (UTC)','fontweight','b','fontsize',font_size);
%    hh = title({[date, node, ' BSR']},'fontweight','b','fontsize',font_size);
%  else
%    datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
%    hh = title({[node, ' BSR']},'fontweight','b','fontsize',font_size);
%    %hh = title({'DIAL Water Vapor (g m^{-3})'},'fontweight','b','fontsize',font_size);
%  end
%  set(gca,'Fontsize',font_size,'Fontweight','b');
%   set(gca,'Zscale', 'log')
%   set(gca,'Colorscale', 'log')
%   set(gca,'Zscale', 'linear')
  
  
  


%% 
if flag.save_figs==1
  %cd('/Users/spuler/Desktop/WV_DIAL_data/plots/') % point to the directory where data is stored 
  %cd('/Volumes/documents/WV_DIAL_data/plots/') % point to the directory where data is stored 
  %cd('/Users/spuler/Desktop/mpd_05_processed_data/Plots/') % point to the directory where data is stored   
  %cd('/Users/lroot/Desktop/mpd/Plots/') % point to the directory where data is stored 
  cd(plot_path) % point to the directory where data is stored 

 
  FigH = figure(1);
  set(gca,'Fontsize',16,'Fontweight','b');  
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920 250]);       %1500 x 300
  name=strcat(date, node, ' Temp_Matlab_multi_avg'); 
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
 
  FigH = figure(2);
  set(gca,'Fontsize',16,'Fontweight','b');  
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920 250]);
  name=strcat(date, node, ' Temp_Matlab_multi'); 
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution    
   
%   FigH = figure(3);
%   set(gca,'Fontsize',16,'Fontweight','b');  
%   set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920 250]);   %1500 x 300
%   name=strcat(date, node, ' Alpha_Matlab_multi'); 
%   print(FigH, name, '-dpng', '-r0') % set at the screen resolution  
%   
%   FigH = figure(4);
%   set(gca,'Fontsize',16,'Fontweight','b');  
%   set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920 250]);   %1500 x 300
%   name=strcat(date, node, ' Alpha0_Matlab_multi'); 
%   print(FigH, name, '-dpng', '-r0') % set at the screen resolution  
%   
%   FigH = figure(5);
%   set(gca,'Fontsize',16,'Fontweight','b');  
%   set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920 250]);   %1500 x 300
%   name=strcat(date, node, ' Alpha1_Matlab_multi'); 
%   print(FigH, name, '-dpng', '-r0') % set at the screen resolution  
%   
%   FigH = figure(6);
%   set(gca,'Fontsize',16,'Fontweight','b');  
%   set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920 250]);   %1500 x 300
%   name=strcat(date, node, ' Alpha2_Matlab_multi'); 
%   print(FigH, name, '-dpng', '-r0') % set at the screen resolution  
  
  FigH = figure(3);
  set(gca,'Fontsize',16,'Fontweight','b');  
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920 250]);   %1500 x 300
  name=strcat(date, node, ' BSR_Matlab_multi'); 
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution
  
end


cd(dd);  % point back to original directory 
toc
