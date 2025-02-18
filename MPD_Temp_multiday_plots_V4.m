clear all; close all;
tic
dd=pwd;

 node = 'MPD04';
 date = '20241203';  
 days = 1; skip = 1; 
 days = 7; skip = 1; 

temp_range_offset = 0; % temp range offset for testing purposes only
flag.ReProcessing = 0; %use the ReProcessed folder to get the temp data
lapse_rate = 0.0065; %standard atmosphere lapse rate
T_lapse = 0.0065; %standard atmosphere lapse rate
     
%serv_path = '/Volumes/fog1/rsfdata/MPD';
serv_path = '/Volumes/smaug1/rsfdata/MPD';
%serv_path = '/Volumes/eol/sci/stillwel/MPDReprocess';
plot_path = '/Users/spuler/Desktop/mpd/Plots/';
font_size = 36; % use this for 2018a version


flag.replot = 1; %replot time vs range images at start of processing (0=off 1=on)
flag.save_figs = 1; %save figures at end of processing (0=off 1=on)
flag.save_data = 0; %save data at end of processing (0=off 1=on)
flag.plot_sonde_data = 0; %plot NetCDF sonde data on top (0=off 1=on)
flag.decimate = 1; %decimate figures to the screen 2x size (1900 pixels x2)
flag.PyHSRL = 0; %decimate figures to the screen 2x size (1900 pixels x2)

%C = importdata('NCAR_C_Map.mat');
addpath '/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing/matplotlib'

if strcmp(node,'MPD01')==1
   cd(strcat(serv_path,'/mpd_01_processed_data/Quickload/FullProcessingBootstrap'))
   if flag.ReProcessing == 1
     cd(strcat(serv_path,'/mpd_01_processed_data/Quickload/ReProcessing/FullProcessingBootstrap'))
   end
   path_node = 'mpd01.';
 elseif strcmp(node,'MPD02')==1
   cd(strcat(serv_path,'/mpd_02_processed_data/Quickload/FullProcessingBootstrap'))
   if flag.ReProcessing == 1
     cd(strcat(serv_path,'/mpd_02_processed_data/Quickload/ReProcessing/FullProcessingBootstrap'))
   end
   path_node = 'mpd02.';
 elseif strcmp(node,'MPD03')==1
    %cd(strcat(serv_path,'/mpd_03_processed_data/Quickload/FullProcessingBootstrap'))
    cd(strcat(serv_path,'/mpd_03_processed_data/Quickload/FullProcessingBootstrap/02_NewPressureAndHitran2020'))
    if flag.ReProcessing == 1
      cd(strcat(serv_path,...
           strcat('/mpd_03_processed_data/Quickload/ReProcessing',ReProc_ver,'/FullProcessingBootstrap')))  
    end
   path_node = 'mpd03.';
 elseif strcmp(node,'MPD04')==1
   cd(strcat(serv_path,'/mpd_04_processed_data/Quickload/FullProcessingBootstrap'))
   if flag.ReProcessing == 1
     cd(strcat(serv_path,'/mpd_04_processed_data/Quickload/ReProcessing/FullProcessingBootstrap'))
   end
   path_node = 'mpd04.';
 elseif strcmp(node,'MPD05')==1
   cd(strcat(serv_path,'/mpd_05_processed_data/Quickload/FullProcessingBootstrap'))
   if flag.ReProcessing == 1
     cd(strcat(serv_path,'/mpd_05_processed_data/Quickload/ReProcessing/FullProcessingBootstrap'))
   end
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
    Temp_duration =  (datenum(date, 'yyyymmdd') + Retrievals.Temperature.TimeStamp/3600/24);
    Temp_range = Retrievals.Temperature.Range - temp_range_offset;
    WV_comb_avg = Retrievals.WaterVapor.Smoothed';
    WV_mask_comb = Retrievals.WaterVapor.Mask';
    WV_duration =  (datenum(date, 'yyyymmdd') + Retrievals.WaterVapor.TimeStamp/3600/24);
    WV_range = Retrievals.WaterVapor.Range;

    if flag.PyHSRL == 1
        %BSR_comb = Retrievals.Python.BackRatio.Value';
        BSR_comb = Retrievals.Python.MPD.BackRatio.Value';
        BSR_duration = (datenum(date, 'yyyymmdd') + Retrievals.Python.MPD.BackRatio.TimeStamp/3600/24);
    else
        BSR_comb = Retrievals.HSRL.Smoothed';
        ABC_comb = Retrievals.HSRL.ABC';
        ABC_mask_comb =  Retrievals.HSRL.Mask';
        BSR_range = Retrievals.HSRL.Range';
        BSR_duration = (datenum(date, 'yyyymmdd') + Retrievals.HSRL.TimeStamp/3600/24);
        P_comb = Retrievals.HSRL.Press';
        T_surf = Retrievals.HSRL.Temp(1,:)';
    end
  else   
    date = datestr(addtodate(datenum(date, 'yyyymmdd'), 1, 'day'), 'yyyymmdd');
    if exist(strcat(path_node, date, '.Matlab.mat'))==2
      load(strcat(path_node, date, '.Matlab.mat'),'Retrievals')
    end
    Temp_comb = vertcat(Temp_comb, Retrievals.Temperature.Value');
    Temp_comb_avg = vertcat(Temp_comb_avg, Retrievals.Temperature.Smoothed');
    Temp_comb_var = vertcat(Temp_comb_var, Retrievals.Temperature.VarianceSm');
    WV_comb_avg = vertcat(WV_comb_avg, Retrievals.WaterVapor.Smoothed');
    WV_mask_comb =  vertcat(WV_mask_comb, Retrievals.WaterVapor.Mask'); 
 if flag.PyHSRL == 1
    BSR_comb = vertcat(BSR_comb, Retrievals.Python.MPD.BackRatio.Value');   
    BSR_duration = vertcat(BSR_duration, (datenum(date, 'yyyymmdd') + Retrievals.Python.MPD.BackRatio.TimeStamp/3600/24));
 else    
    BSR_comb = vertcat(BSR_comb, Retrievals.HSRL.Smoothed');  
    ABC_comb = vertcat(ABC_comb, Retrievals.HSRL.ABC');  
    ABC_mask_comb =  vertcat(ABC_mask_comb, Retrievals.HSRL.Mask');
    P_comb = vertcat(P_comb, Retrievals.HSRL.Press');
    T_surf = vertcat(T_surf, Retrievals.HSRL.Temp(1,:)');
    BSR_duration = vertcat(BSR_duration, (datenum(date, 'yyyymmdd') + Retrievals.HSRL.TimeStamp/3600/24));
 end
    Temp_duration = vertcat(Temp_duration, (datenum(date, 'yyyymmdd') + Retrievals.Temperature.TimeStamp/3600/24));
    WV_duration = vertcat(WV_duration, (datenum(date, 'yyyymmdd') + Retrievals.WaterVapor.TimeStamp/3600/24));
  end
end

% apply preliminary mask based on O_2 order absorption zeroth order variance 
Temp_comb_avg(Temp_comb_var > 5^2)= nan;
ABC_comb(ABC_mask_comb==1)= nan;


% put T_surf and P (BSR range and duration) on same grid 
% as Temp_comb (Temp range and duration) 

T_surf_grid = interp1(BSR_duration(~isnan(T_surf)), T_surf(~isnan(T_surf)), Temp_duration, 'linear');
  [range_in, duration_in] = meshgrid(BSR_range, BSR_duration);
  [range_out,duration_out] = meshgrid(Temp_range, Temp_duration);
P_grid = interp2(range_in, duration_in, P_comb, range_out, duration_out, 'linear');
  [range_in, duration_in] = meshgrid(WV_range, WV_duration);
WV_grid = interp2(range_in, duration_in, WV_comb_avg, range_out, duration_out, 'linear');
WV_surf = WV_grid(:,4);
P_surf = P_grid(:,1);

% calculate virtual potential temperature
  const.M_wv = 18.015;  %molar mass of water molecule
  const.M_air = 28.97; % molar mass gm/mol of air
  const.N_A = 6.022E23; % Avagadro number
  const.k_B = 1.3806488e-23; % (J/K)
  % water vapor number density 
  n_wv = WV_grid.*const.N_A./const.M_wv;
  n_wv_surf = WV_surf.*const.N_A./const.M_wv;
  % water vapor partial pressure (Pa to atm convert)
  e_wv =  n_wv.*const.k_B.*Temp_comb_avg./101300;
  e_wv_surf =  WV_surf.*const.k_B.*(T_surf_grid)./101300;
  % calculate virtual temperature
  comb_T_v = Temp_comb_avg./(1-((e_wv)./P_grid).*(1-const.M_wv/const.M_air));
  comb_T_v_surf = T_surf_grid./(1-((e_wv_surf)./P_surf).*(1-const.M_wv/const.M_air));
  % calculate virtual potential temperature
  P_0 = 0.987; % reference pressure in atm
  comb_T_p = comb_T_v.*(P_0./P_grid).^0.286; 
  comb_T_p_surf = comb_T_v_surf.*(P_0./P_surf).^0.286;

%   B = repmat(comb_T_p_surf, 1, size(comb_T_p,2));
%   B_test = B+2; % offset virtual potential temp at surface by 1.5deg
%   comb_T_p(comb_T_p>B_test) = nan;



    
  if days == 1
   xData =  linspace(fix(min(Temp_duration)),  ceil(max(Temp_duration)), 25);
 else
   xData =  linspace( fix(min(Temp_duration)),  ceil(max(Temp_duration)), round((ceil(max(Temp_duration))-fix(min(Temp_duration)))/skip)+1 );
  % xData =  linspace(fix(min(duration)),  round(max(duration)), 36);
  end

scrsz = [1  1  1920 1200];
Scrnsize=[scrsz(4)/2 scrsz(4)/10 scrsz(3)/1 scrsz(4)/2];


x = Temp_duration';
y = Temp_range./1000;
 
 
 Z = Temp_comb_avg'-273.15; 
 %Z = Temp_comb_avg'; 
 %Z(isnan(Alpha_comb)) = nan;
 %Z = Alpha_comb'; 
 figure('Position',Scrnsize)
 set(gcf,'renderer','zbuffer');
 h = pcolor(x,y,Z);
 set(h, 'EdgeColor', 'none');
 colorbar('EastOutside');
 axis([fix(min(x)) ceil(max(x)) 0 6])
 caxis([-15 20]);
 %caxis([250 310]);
 colormap(plasma)
%  colormap(jet)
%  shading interp
 ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
 set(gca, 'XTick',  xData)
 set(gca,'TickDir','out');
 set(gca,'TickLength',[0.005; 0.0025]);
 if days == 1
   datetick('x','HH:MM','keeplimits', 'keepticks');
   xlabel('Time (UTC)','fontweight','b','fontsize',font_size);
   hh = title({[date, node, ' T Standard (C)']},'fontweight','b','fontsize',font_size);
 else
   datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
%    hh = title({[node, ' Temperature (K)']},'fontweight','b','fontsize',font_size);
   hh = title({[node, ' Temperature (C)']},'fontweight','b','fontsize',font_size);
 end
 set(gca,'Fontsize',font_size,'Fontweight','b');
 grid on
   
%  Z = Temp_comb'-273.15;  
% Z = Temp_comb'; 
 Z =  comb_T_p'; 
 figure('Position',Scrnsize)
 set(gcf,'renderer','zbuffer');
 h = pcolor(x,y,Z);
 set(h, 'EdgeColor', 'none');
 colorbar('EastOutside');
 axis([fix(min(x)) ceil(max(x)) 0 6])
 %caxis([0 35]); % Celsius
 caxis([314 326]); 
 %colormap(jet)
 colormap(plasma)
 %shading interp
 ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
 set(gca, 'XTick',  xData)
 set(gca,'TickDir','out');
 set(gca,'TickLength',[0.005; 0.0025]);
 if days == 1
   datetick('x','HH:MM','keeplimits', 'keepticks');
   xlabel('Time (UTC)','fontweight','b','fontsize',font_size);
   hh = title({[date, node, 'Virt Potential Temperature (K)']},'fontweight','b','fontsize',font_size);
 else
   datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
   hh = title({[node, ' Virt Potential Temperature (K)']},'fontweight','b','fontsize',font_size);
%    hh = title({[node, ' Temperature (C)']},'fontweight','b','fontsize',font_size);
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

 x = BSR_duration';
 if flag.PyHSRL == 1
   y = Retrievals.Python.MPD.BackRatio.Range./1000;  % note this is Matt's version and not mine
 else
   y = Retrievals.HSRL.Range./1000;  % note this is Matt's version and not mine
 end
 Z = BSR_comb';  
 figure('Position',Scrnsize)
 set(gcf,'renderer','zbuffer');
 h = pcolor(x,y,Z);
 set(h, 'EdgeColor', 'none');
 colorbar('EastOutside');
 axis([fix(min(x)) ceil(max(x)) 0 6])
 caxis([1 10]);
 %caxis([1e-1 1e3]);
 %colormap(jet)
 colormap(flipud(viridis))
 colormap(viridis)
 %shading interp
 ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
 set(gca, 'XTick',  xData)
 set(gca,'TickDir','out');
 set(gca,'TickLength',[0.005; 0.0025]);
 if days == 1
   datetick('x','HH:MM','keeplimits', 'keepticks');
   xlabel('Time (UTC)','fontweight','b','fontsize',font_size);
   hh = title({[date, node, ' BSR']},'fontweight','b','fontsize',font_size);
 else
   datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
   hh = title({[node, ' BSR']},'fontweight','b','fontsize',font_size);
   %hh = title({'DIAL Water Vapor (g m^{-3})'},'fontweight','b','fontsize',font_size);
 end
 set(gca,'Fontsize',font_size,'Fontweight','b');
  set(gca,'Zscale', 'log')
  set(gca,'Colorscale', 'log')
  set(gca,'Zscale', 'linear')
  
  
 x = BSR_duration';
 if flag.PyHSRL == 1
   y = Retrievals.Python.MPD.BackRatio.Range./1000;  % note this is Matt's version and not mine
 else
   y = Retrievals.HSRL.Range./1000;  % note this is Matt's version and not mine
 end
 Z = ABC_comb';  
 figure('Position',Scrnsize)
 set(gcf,'renderer','zbuffer');
 h = pcolor(x,y,Z);
 set(h, 'EdgeColor', 'none');
 colorbar('EastOutside');
 axis([fix(min(x)) ceil(max(x)) 0 6])
 %caxis([1 10]);
 caxis([1e-8 1e-5]);
 %colormap(jet)
 colormap(flipud(viridis))
 colormap(viridis)
 %shading interp
 ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
 set(gca, 'XTick',  xData)
 set(gca,'TickDir','out');
 set(gca,'TickLength',[0.005; 0.0025]);
 if days == 1
   datetick('x','HH:MM','keeplimits', 'keepticks');
   xlabel('Time (UTC)','fontweight','b','fontsize',font_size);
   hh = title({[date, node, ' BSR']},'fontweight','b','fontsize',font_size);
 else
   datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
   hh = title({[node, ' Backscatter Coefficient [m^{-1} sr^{-1}]']},'fontweight','b','fontsize',font_size);
   %hh = title({'DIAL Water Vapor (g m^{-3})'},'fontweight','b','fontsize',font_size);
 end
 set(gca,'Fontsize',font_size,'Fontweight','b');
  set(gca,'Zscale', 'log')
  set(gca,'Colorscale', 'log')
  set(gca,'Zscale', 'linear')
  
  
  


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
  
  FigH = figure(4);
  set(gca,'Fontsize',16,'Fontweight','b');  
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920*1.5 250]);   %1500 x 300
  name=strcat(date, node, ' ABC_Matlab_multi'); 
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution
  
  
end


cd(dd);  % point back to original directory 
toc
