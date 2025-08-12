
clear all; close all
% day = datenum('23 Jul 2023'); % automated to pick the first day
UTC_convert = 0;  
flag.single_day = 0 % used to plot single days
flag.overlay_HRRR = 0;
if flag.single_day == 1
    UTC_convert = 7;  % convert to local time
end

offset1 = 0.5; % parcel method surface temp offset
offset2 = 1; % parcel method surface temp offset
offset3 = 1.5; % parcel method surface temp offset


skip = 1
node = 'MPD03';
flag.PTV = 1 %use the PTV data which has temperature
addpath '/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing/matplotlib/';

serv_path = '/Volumes/eol/sci/mhayman';
cd(strcat(serv_path,'/DIAL/Processed_Data/M2HATS/5min_release1.0/qc_masked/'))
%cd(strcat(serv_path,'/DIAL/Processed_Data/ECLIPSE/initial_test/'))

if flag.PTV == 1
%  serv_path = '/Volumes/eol/smaug1/rsfdata/MPD';
%  cd(strcat(serv_path,'/mpd_03_processed_data/PTV/temperature/M2HATS1.0'))
  serv_path = '/Volumes/eol/sci/mhayman';
  cd(strcat(serv_path,'/DIAL/Processed_Data/M2HATS/ptv1.1/'))
%   serv_path = '/Volumes/eol/scr-tmp/mhayman';
%   cd(strcat(serv_path,'/O2DIAL/PTV/gp_search'))
end
%cd(strcat(serv_path,'/mpd_03_processed_data/Python'))
d_read_data = pwd; % get the current path

plot_path = '/Users/spuler/Desktop';
cd(strcat(plot_path,'/mpd/Plots'))
d_save_data = pwd; %set the plot save path
flag.save_data = 0;  %save data at end of processing (0=off 1=on)
low_range_mask = 0;

cd(d_read_data);
[Pythonfilename, Pythondir] = uigetfile('*.*','Select the sonde file', 'MultiSelect', 'on');
jj=1;

 %variable{1} = 'Aerosol_Backscatter_Coefficient_variance';
  
 variable{1} = 'time';
 variable{2} = 'range';
 variable{3} = 'Temperature';
 variable{4} = 'Temperature_mask';
 variable{5} = 'Temperature_uncertainty';
 variable{6} = 'Absolute_Humidity'; 
 variable{7} = 'Absolute_Humidity_mask';
 variable{8} = 'Absolute_Humidity_variance'; 
 if flag.PTV == 1
   variable{8} = 'Absolute_Humidity_uncertainty'; 
 end
 variable{11} = 'Aerosol_Backscatter_Coefficient_variance';
 variable{9} = 'Aerosol_Backscatter_Coefficient';
 variable{10} = 'Aerosol_Backscatter_Coefficient_mask';
 variable{11} = 'Aerosol_Backscatter_Coefficient_variance';
 if flag.PTV == 1
   variable{11} = 'Aerosol_Backscatter_Coefficient_uncertainty';
   variable{12} = 'Pressure_Estimate';
   variable{13} = 'Pressure_Estimate_mask';
   variable{14} = 'Lapse_Rate';
   variable{15} = 'Lapse_Rate_mask';
   variable{16} = 'PBLH_Model';
   variable{17} =  'Virtual_Potential_Temperature';
   variable{18} =  'Virtual_Potential_Temperature_mask';
   variable{19} =  'Lifting_Level';
   variable{20} =  'Lifting_Level_mask';
   variable{21} =  'Surface_Temperature';
   variable{22} =  'Surface_Pressure';
   variable{23} =  'Surface_Absolute_Humidity';
   variable{24} = 'Backscatter_Photon_Counts_828'; 
 end

    

for jj = 1:size(Pythonfilename,2)
  filename = Pythonfilename{jj};
  %date = filename(end-15:end-10);
  %n = datenum(date, 'yymmdd');
  date = filename(end-15:end-8);
  n = datenum(date, 'yyyymmdd');
  
  ncid = netcdf.open(filename, 'NC_NOWRITE');
    %ncdisp(filename, '/', 'min') % use this to display all variables
    ncdisp(filename) % use this to display all variables
    time{jj} = ncread(filename,variable{1}); 
    alt{jj} = ncread(filename,variable{2});
    
if flag.PTV == 1
     T{jj}  = ncread(filename,variable{3});  
     T_mask{jj} = ncread(filename,variable{4}); 
     T_var{jj} = ncread(filename,variable{5}); 
     T{jj}(T_mask{jj} == 1) = nan;
     T_var{jj}(T_mask{jj} == 1) = nan; 
     P{jj} =  ncread(filename,variable{12});
     P_mask{jj} = ncread(filename,variable{13});
     L{jj} =  ncread(filename,variable{14});
     L_mask{jj} = ncread(filename,variable{15});
     PBLH{jj} =  ncread(filename,variable{16});
     T_vp{jj} =  ncread(filename,variable{17});
     T_vp_mask{jj} = ncread(filename,variable{18});
     Lift{jj} =  ncread(filename,variable{19});
     Lift_mask{jj} = ncread(filename,variable{20});
     T_surf{jj} =  ncread(filename,variable{21});
     P_surf{jj} =  ncread(filename,variable{22});
     AH_surf{jj} =  ncread(filename,variable{23});
end
    
    AH{jj}  = ncread(filename,variable{6});  
    AH_mask{jj} = ncread(filename,variable{7}); 
    AH_var{jj} = ncread(filename,variable{8}); 
    AH{jj}(AH_mask{jj} == 1) = nan;
    AH_var{jj}(AH_mask{jj} == 1) = nan; 
    AH{jj}(AH_var{jj} > 3) = nan;
    
    ABC{jj}  = ncread(filename,variable{9});   
    ABC_mask{jj} = ncread(filename,variable{10}); 
    ABC_var{jj} = ncread(filename,variable{11});
    ABC{jj}(ABC_mask{jj} == 1) = nan;
    ABC_var{jj}(ABC_mask{jj} == 1) = nan; 
   % ABC{jj}(ABC_var{jj} > 1e-11) = nan;
    ABC{jj}(ABC{jj} < 1e-12) = nan;
    Counts{jj} = ncread(filename,variable{24});   
    
    if flag.PTV == 1
      % P{jj}(P_mask{jj} == 1) = nan;
 %     L{jj}(L_mask{jj} == 1) = nan;
 %     T_vp{jj}(T_vp_mask{jj} == 1) = nan;
 %     Lift{jj}(Lift_mask{jj} == 1) = nan;
    end
    
  netcdf.close(ncid); 

  %convert from Unix time to date number (days since Jan 0 0000) 
  duration{jj} =  n+double(time{jj}/3600/24);
end

 % mask the lowest gates if desired
 % low_mask = repmat(alt{jj}, 1, size(AH{jj}, 2)); 
 % AH{jj}(low_mask < low_range_mask)= nan;
 % AH_var{jj}(low_mask < low_range_mask)= nan;
 
 for jj = 1:size(Pythonfilename,2)
  if jj == 1
      comb_duration = duration{jj};
      comb_AH = AH{jj};
      comb_AH_var = AH_var{jj};
      comb_ABC = ABC{jj};
      comb_ABC_var = ABC_var{jj};
      if flag.PTV == 1
        comb_T = T{jj};
        comb_T_var = T_var{jj};
        comb_P = P{jj};
        comb_L = L{jj};
        comb_PBLH = PBLH{jj};
        comb_T_vp = T_vp{jj};  
        comb_Lift = Lift{jj};  
        comb_T_surf = T_surf{jj}; 
        comb_P_surf = P_surf{jj}; 
        comb_AH_surf = AH_surf{jj}; 
        comb_Counts = Counts{jj}; 
      end
  else
      comb_duration = [comb_duration; duration{jj}];
      % find the maximum range to accumulate (in case it changes)
      % max_range = min(cellfun('size',alt,1))
      comb_AH = [comb_AH AH{jj}];
      comb_AH_var = [comb_AH_var AH_var{jj}];
      comb_ABC = [comb_ABC ABC{jj}];
      comb_ABC_var = [comb_ABC_var ABC_var{jj}];
      if flag.PTV == 1
        comb_T = [comb_T T{jj}];
        comb_T_var = [comb_T_var T_var{jj}];
        comb_P = [comb_P P{jj}];
        comb_L = [comb_L L{jj}];
        comb_PBLH = [comb_PBLH; PBLH{jj}];
        comb_T_vp = [comb_T_vp T_vp{jj}];
        comb_Lift = [comb_Lift; Lift{jj}];
        comb_T_surf = [comb_T_surf; T_surf{jj}];
        comb_P_surf = [comb_P_surf; P_surf{jj}];
        comb_AH_surf = [comb_AH_surf; AH_surf{jj}];
        comb_Counts = [comb_Counts Counts{jj}];
      end
  end
  
end

scrsz = get(0,'ScreenSize');
date=datestr(n, 'yyyy-mmm-dd');
plot_size1 = [scrsz(4)/1.5 scrsz(4)/10 scrsz(3)/1.5 scrsz(4)/3];
font_size = 16;

x = comb_duration-(UTC_convert/24);
y = (alt{1}./1000);
xData =  linspace( fix(min(x)),  ceil(max(x)), round((ceil(max(x))-fix(min(x)))/skip)+1 );
xData_m =  linspace( fix(min(x)),  ceil(max(x)), round((ceil(max(x))-fix(min(x)))/skip*24/2)+1 );


if flag.PTV == 1
% % plot the T
 % Z = real(comb_T-273.15);
  Z = real(comb_T);
  figure3 = figure('Position',plot_size1);
  set(gcf,'renderer','zbuffer');
  h = pcolor(x, y, Z);
  set(h, 'EdgeColor', 'none'); 
  axis xy; colorbar('EastOutside'); 
%  caxis([-15 30]);
  caxis([268 298]);
  axis([fix(min(x)) ceil(max(x)) 0 6]) 
%  shading interp
  set(gca, 'XTick',  xData)
  set(gca,'XMinorTick','on')
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  hh = title({[node, ' Temp (K)']},...
       'fontweight','b','fontsize',font_size);  
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
 % colormap(jet)
  colormap(plasma)
  set(gca,'Fontsize',font_size,'Fontweight','b');
  
%  % % plot the P or Potenial Temperature
% virtual potential temperature
  const.M_wv = 18.015;  %molar mass of water molecule
  const.M_air = 28.97; % molar mass gm/mol of air
  const.N_A = 6.022E23; % Avagadro number
  const.k_B = 1.3806488e-23; % (J/K)
  % water vapor number density 
  n_wv = comb_AH.*const.N_A./const.M_wv;
  n_wv_surf = comb_AH_surf.*const.N_A./const.M_wv;
  % water vapor partial pressure (Pa to atm convert)
  e_wv =  n_wv.*const.k_B.*comb_T./101300;
  e_wv_surf =  n_wv_surf.*const.k_B.*(comb_T_surf)./101300;
  % calculate virtual temperature
  comb_T_v = comb_T./(1-((e_wv)./comb_P).*(1-const.M_wv/const.M_air));
  comb_T_v_surf = comb_T_surf./(1-((e_wv_surf)./comb_P_surf).*(1-const.M_wv/const.M_air));
  % calculate virtual potential temperature
  P_0 = 0.987; % reference pressure in atm
  comb_T_p = comb_T_v.*(P_0./comb_P).^0.286; 
  comb_T_p_surf = comb_T_v_surf.*(P_0./comb_P_surf).^0.286;

  
   % add in the surface T 
   y(1)= 0;
   y(2)= 0.1;
   comb_T_p(1,:) = comb_T_p_surf; 
  
  B = repmat(comb_T_p_surf, 1, size(comb_T_p,1));

  VPT_surf = B'+offset1; % offset surface virtual potential temp 
  PBLH = comb_T_p; % set the PBLH equal to the VPT to set variable dim
  PBLH(comb_T_p>VPT_surf) = nan;
    PBLH(VPT_surf-comb_T_p>2) = nan;  % threshold of 2K difference
  PBLH_check_top = ~isnan(PBLH);
  Indices = arrayfun(@(x) find(PBLH_check_top(:,x), 1, 'last'), 1:size(PBLH, 2), 'UniformOutput', false);
  idx = cellfun(@isempty,Indices);
  Indices(idx) = {1};
  BLH_Indices = cell2mat(Indices);
  Thermo_BLH1 = y(BLH_Indices); 
  Thermo_BLH1(Thermo_BLH1==0) = nan; 
  
  VPT_surf = B'+offset2; % offset surface virtual potential temp 
  PBLH = comb_T_p;
  PBLH (PBLH >VPT_surf) = nan;
    PBLH(VPT_surf-comb_T_p>2) = nan;  % threshold of 2K difference
  PBLH_check_top = ~isnan(PBLH);
  Indices = arrayfun(@(x) find(PBLH_check_top(:,x), 1, 'last'), 1:size(PBLH, 2), 'UniformOutput', false);
  idx = cellfun(@isempty,Indices);
  Indices(idx) = {1};
  BLH_Indices = cell2mat(Indices);
  Thermo_BLH2 = y(BLH_Indices); 
  Thermo_BLH2(Thermo_BLH2==0) = nan; 
  
  VPT_surf = B'+offset3; % offset surface virtual potential temp  
  PBLH = comb_T_p;
  PBLH (PBLH >VPT_surf) = nan;
    PBLH(VPT_surf-comb_T_p>2) = nan;  % threshold of 2K difference
  PBLH_check_top = ~isnan(PBLH);
  Indices = arrayfun(@(x) find(PBLH_check_top(:,x), 1, 'last'), 1:size(PBLH, 2), 'UniformOutput', false);
  idx = cellfun(@isempty,Indices);
  Indices(idx) = {1};
  BLH_Indices = cell2mat(Indices);
  Thermo_BLH3 = y(BLH_Indices); 
  Thermo_BLH3(Thermo_BLH3==0) = nan; 
  
%   B_test = B'+offset4; % offset virtual potential temp at surface by 2deg
%   PBLH = comb_T_p;
%   PBLH (PBLH >B_test) = nan;
%   B_test_2 = ~isnan(PBLH);
%   Indices = arrayfun(@(x) find(B_test_2(:,x), 1, 'last'), 1:size(PBLH, 2), 'UniformOutput', false);
%   idx = cellfun(@isempty,Indices);
%   Indices(idx) = {1};
%   BLH_Indices = cell2mat(Indices);
%   Thermo_BLH4 = y(BLH_Indices); 
%   Thermo_BLH4(Thermo_BLH4==0) = nan; 
  
  
  % set this to just plot the first day selected for plots
  day = datenum( Pythonfilename{1}(end-15:end-8), 'yyyymmdd');

   offset_avg = mean([offset1 offset2 offset3]);
   Z = real((B'+ offset_avg)-comb_T_p);
   %Z = real(FY);
   %Z = real(comb_P);
   %Z = real(comb_T_vp);  % directly from Matt's processing   
 
   figure4 = figure('Position',plot_size1);
   set(gcf,'renderer','zbuffer');
   h = pcolor(x, y, Z);
   set(h, 'EdgeColor', 'none'); 
   axis xy; colorbar('EastOutside'); 
   caxis([-10 10]);
   %caxis([305 325]);
   axis([fix(min(x)) ceil(max(x)) 0 6]) 
   set(gca, 'XTick',  xData)
   set(gca,'XMinorTick','on')
   xAx = get(gca,'XAxis');
   xAx.MinorTickValues=xData_m;
    set(gca,'TickDir','out');
   set(gca,'TickLength',[0.005; 0.0025]);
   hh = title({[node, ' Virtual Potential Temperature Difference (Surface-Atmosphere) [K]']},...
        'fontweight','b','fontsize',font_size);  
   ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
   datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
 %  colormap(jet)
 %  colormap(plasma)
   colormap(redblue)
   set(gca,'Fontsize',font_size,'Fontweight','b'); 
   grid on 
   grid(gca, 'minor')
   if flag.single_day == 1 
       xlim([day day+1]) % just show one single day
     % xlim([day+5/24 day+20/24]) % to match Luke's MLH paper
     % ylim([0 4.5])
   end

  
  Thermo_BLH_avg = mean([Thermo_BLH1 Thermo_BLH2 Thermo_BLH3],2,'omitnan');
  % add the overlay of PBLH
  ax2 = axes(figure4);
  plot(x, Thermo_BLH_avg, 'ks')
  % hold on 
  %   plot(x, Thermo_BLH1, 'g+', 'MarkerSize', 1)
  %   plot(x, Thermo_BLH2, 'g+', 'MarkerSize', 1)
  %   plot(x, Thermo_BLH3, 'g+', 'MarkerSize', 1)
  % hold off
  colorbar('EastOutside'); 
  axis([fix(min(x)) ceil(max(x)) 0 6])
  ax2.Color = 'none';
  ax2.XAxisLocation = 'top';
  ax2.YAxisLocation = 'right';
  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
  set(gca,'Fontsize',font_size,'Fontweight','b');
  ax2.YColor = 'none';
  ax2.XColor = 'none';
  % Kludge to remove the colorbar font, but not rescale the axis
  c = colorbar;
  c.FontSize = 0.01;
  if flag.single_day == 1 
       xlim([day day+1]) % just show one single day
     % xlim([day+5/24 day+20/24]) % to match Luke's MLH paper
     % ylim([0 4.5])
  end
  
if flag.overlay_HRRR == 1 
  % add the overlay of PBLH
   ax2 = axes(figure4);
   plot(x, comb_PBLH/1000, 'ro')
   colorbar('EastOutside'); 
   axis([fix(min(x)) ceil(max(x)) 0 6])
   ax2.Color = 'none';
   ax2.XAxisLocation = 'top';
   ax2.YAxisLocation = 'right';
   datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
   set(gca,'Fontsize',font_size,'Fontweight','b');
   ax2.YColor = 'none';
   ax2.XColor = 'none';
   % Kludge to remove the colorbar font, but not rescale the axis
   c = colorbar;
   c.FontSize = 0.01;
   if flag.single_day == 1 
        xlim([day day+1]) % just show one single day
      % xlim([day+5/24 day+20/24]) % to match Luke's MLH paper
      % ylim([0 4.5])
   end
 end

   
%   % add the overlay of Lifting Level
%   ax3 = axes(figure4);
%   plot(comb_duration, comb_Lift/1000, 'go')
%   colorbar('EastOutside'); 
%   axis([fix(min(x)) ceil(max(x)) 0 6])
%   ax3.Color = 'none';
%   ax3.XAxisLocation = 'top';
%   ax3.YAxisLocation = 'right';
%   datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
%   set(gca,'Fontsize',font_size,'Fontweight','b');
%   ax3.YColor = 'none';
%   ax3.XColor = 'none';
%   % Kludge to remove the colorbar font, but not rescale the axis
%   c = colorbar;
%   c.FontSize = 0.01; 
   
   
   
%  
%   % % plot the Lapse Rate
%   Z = real(comb_L);
%   figure5 = figure('Position',plot_size1);
%   set(gcf,'renderer','zbuffer');
%   h = pcolor(x, y, Z);
%   set(h, 'EdgeColor', 'none'); 
%   axis xy; colorbar('EastOutside'); 
%   caxis([-15 15]);
%   axis([fix(min(x)) ceil(max(x)) 0 6]) 
%   set(gca, 'XTick',  xData)
%   set(gca,'TickDir','out');
%   set(gca,'TickLength',[0.005; 0.0025]);
%   hh = title({[node, ' Lapse Rate (K/km)']},...
%        'fontweight','b','fontsize',font_size);  
%   ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
%   datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
% %  colormap(jet)
%   set(gca,'Fontsize',font_size,'Fontweight','b');
  
end


% plot the AH
  Z = real(comb_AH);
  figure1 = figure('Position',plot_size1);
  set(gcf,'renderer','zbuffer');
  h = pcolor(x, y, Z);
  set(h, 'EdgeColor', 'none'); 
  axis xy; 
  colorbar('EastOutside'); 
  caxis([0 6]);
  axis([fix(min(x)) ceil(max(x)) 0 6]) 
%  shading interp
  set(gca, 'XTick',  xData)
  set(gca,'XMinorTick','on')
  xAx = get(gca,'XAxis');
  xAx.MinorTickValues=xData_m;
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  hh = title({[node, ' Absolute Humidity (g m^{-3})']},...
       'fontweight','b','fontsize',font_size);  
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
%  colormap(jet)
  colormap(CM_YlGnBu(64))
  set(gca,'Fontsize',font_size,'Fontweight','b');
  if flag.single_day == 1 
       xlim([day day+1]) % just show one single day
     % xlim([day+5/24 day+20/24]) % to match Luke's MLH paper
     % ylim([0 4.5])
   end
  
if flag.overlay_HRRR == 1 
  % add the overlay of PBLH
   ax2 = axes(figure1);
   plot(x, comb_PBLH/1000, 'ro')
   colorbar('EastOutside'); 
   axis([fix(min(x)) ceil(max(x)) 0 6])
   ax2.Color = 'none';
   ax2.XAxisLocation = 'top';
   ax2.YAxisLocation = 'right';
   datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
   set(gca,'Fontsize',font_size,'Fontweight','b');
   ax2.YColor = 'none';
   ax2.XColor = 'none';
   % Kludge to remove the colorbar font, but not rescale the axis
   c = colorbar;
   c.FontSize = 0.01;
   if flag.single_day == 1 
        xlim([day day+1]) % just show one single day
      % xlim([day+5/24 day+20/24]) % to match Luke's MLH paper
      % ylim([0 4.5])
   end
 end
  
 
 % plot the atmospheric backscatter coefficient 
  Z = real(comb_ABC);
  figure2 = figure('Position',plot_size1);
  set(gcf,'renderer','zbuffer');
  h = pcolor(x, y, Z);
  set(h, 'EdgeColor', 'none'); 
  axis xy; 
  colorbar('EastOutside'); 
  %caxis([0 12]);
  %ylim([-.1 6]);
  axis([fix(min(x)) ceil(max(x)) 0 6])
  %  shading interp
  set(gca, 'XTick',  xData)
  set(gca,'XMinorTick','on')
  xAx = get(gca,'XAxis');
  xAx.MinorTickValues=xData_m;
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  set(gca,'Zscale', 'log')
  set(gca,'Colorscale', 'log')
  set(gca,'Zscale', 'linear')
  hh = title({[node, ' Aerosol Backscatter Coefficient m^{-1} sr^{-1}']},...
       'fontweight','b','fontsize',font_size);     
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
 % datetick('x','HH:MM');
 % colormap(jet)
  colormap(viridis)
  set(gca,'Fontsize',font_size,'Fontweight','b');
  caxis([5e-9 1e-6]);
  if flag.single_day == 1 
       xlim([day day+1]) % just show one single day
     % xlim([day+5/24 day+20/24]) % to match Luke's MLH paper
     % ylim([0 4.5])
  end
  
if flag.overlay_HRRR == 1 
  % add the overlay of PBLH
   ax2 = axes(figure2);
   plot(x, comb_PBLH/1000, 'ro')
   colorbar('EastOutside'); 
   axis([fix(min(x)) ceil(max(x)) 0 6])
   ax2.Color = 'none';
   ax2.XAxisLocation = 'top';
   ax2.YAxisLocation = 'right';
   datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
   set(gca,'Fontsize',font_size,'Fontweight','b');
   ax2.YColor = 'none';
   ax2.XColor = 'none';
   % Kludge to remove the colorbar font, but not rescale the axis
   c = colorbar;
   c.FontSize = 0.01;
   if flag.single_day == 1 
        xlim([day day+1]) % just show one single day
      % xlim([day+5/24 day+20/24]) % to match Luke's MLH paper
      % ylim([0 4.5])
   end
 end

% plot the atmospheric backscatter coefficient 
  Z = real(comb_Counts);
  figure2 = figure('Position',plot_size1);
  set(gcf,'renderer','zbuffer');
  h = pcolor(x, y, Z);
  set(h, 'EdgeColor', 'none'); 
  axis xy; 
  colorbar('EastOutside'); 
  %caxis([0 12]);
  %ylim([-.1 6]);
  axis([fix(min(x)) ceil(max(x)) 0 6])
  %  shading interp
  set(gca, 'XTick',  xData)
  set(gca,'XMinorTick','on')
  xAx = get(gca,'XAxis');
  xAx.MinorTickValues=xData_m;
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  set(gca,'Zscale', 'log')
  set(gca,'Colorscale', 'log')
  set(gca,'Zscale', 'linear')
  hh = title({[node, ' Attenuated Backscatter ']},...
       'fontweight','b','fontsize',font_size);     
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
 % datetick('x','HH:MM');
 % colormap(jet)
  colormap(viridis)
  set(gca,'Fontsize',font_size,'Fontweight','b');
  caxis([1e3 1e6]);
  if flag.single_day == 1 
       xlim([day day+1]) % just show one single day
     % xlim([day+5/24 day+20/24]) % to match Luke's MLH paper
     % ylim([0 4.5])
  end


 
 cd(strcat(plot_path,'/mpd/Plots'))

if flag.single_day == 1
    plot_size = [1 1 1000 385];
else
    plot_size = [1 1 1920 250];
end

   
  FigH = figure(3);
  set(gca,'Fontsize',16,'Fontweight','b'); 
  set(FigH, 'PaperUnits', 'points', 'PaperPosition',  plot_size);  
  name=strcat(date, node, ' WV_Python_multi'); 
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
  
  FigH = figure(4);
  set(gca,'Fontsize',16,'Fontweight','b'); 
  set(FigH, 'PaperUnits', 'points', 'PaperPosition',  plot_size);
  name=strcat(date, node, ' Back_Coeff_Python_multi'); 
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution
  
  FigH = figure(2);
  set(gca,'Fontsize',16,'Fontweight','b'); 
  set(FigH, 'PaperUnits', 'points', 'PaperPosition',  plot_size);
  name=strcat(date, node, ' VPT_Diff_multi'); 
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution



 if flag.PTV == 1
  FigH = figure(1);
  set(gca,'Fontsize',16,'Fontweight','b'); 
  set(FigH, 'PaperUnits', 'points', 'PaperPosition',  plot_size);
  name=strcat(date, node, ' T_Python_multi'); 
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
  
 FigH = figure(2);
  set(gca,'Fontsize',16,'Fontweight','b'); 
  set(FigH, 'PaperUnits', 'points', 'PaperPosition',  plot_size);
  name=strcat(date, node, ' Potential_T_multi'); 
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution
  
  FigH = figure(3);
  set(gca,'Fontsize',16,'Fontweight','b'); 
  set(FigH, 'PaperUnits', 'points', 'PaperPosition',  plot_size);
  name=strcat(date, node, ' WV_Python_multi'); 
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
  
  FigH = figure(4);
  set(gca,'Fontsize',16,'Fontweight','b'); 
  set(FigH, 'PaperUnits', 'points', 'PaperPosition',  plot_size);
  name=strcat(date, node, ' Back_Coeff_Python_multi'); 
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution
  
  FigH = figure(5);
  set(gca,'Fontsize',16,'Fontweight','b'); 
  set(FigH, 'PaperUnits', 'points', 'PaperPosition',  plot_size);
  name=strcat(date, node, ' 828_Counts_multi'); 
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution

 
 end  
  

 if flag.save_data == 1
  
  cd(d_save_data);     
  range = alt{1}'; 
  N_avg_comb = (comb_AH./1e6.*6.022E23./18.015);
  duration = duration;
  % test1 = datestr(comb_duration(1:2), 'yyyy-mmm-dd HH:MM');
  % test2 = Thermo_BLH2(1:2);
  test3 = [comb_duration Thermo_BLH1 Thermo_BLH2 Thermo_BLH3];  
  figure(1)
  plot(test3(:,1), test3(:,2))
  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');


  
  
  %cd('/Users/spuler/Desktop/WV_DIAL_data') % point to the directory where data is stored 
  name=strcat(date, '_combined');
  if strcmp(node,'MPD01')==1
      MPD01.N_avg_comb = N_avg_comb;
     % MPD01.RB_comb = RB_comb;
      MPD01.range = range;
      MPD01.time = duration;
      save(name, 'MPD01')
  elseif strcmp(node,'MPD02')==1
      MPD02.N_avg_comb = N_avg_comb;
    %  MPD02.RB_comb = RB_comb;
      MPD02.range = range;
      MPD02.time = duration;
      save(name, 'MPD02')
  elseif strcmp(node,'MPD03')==1
      MPD03.N_avg_comb = N_avg_comb;
    %  MPD03.RB_comb = RB_comb;
      MPD03.range = range;
      MPD03.time = duration;
      save(name, 'MPD03')
  elseif strcmp(node,'MPD04')==1
      MPD04.N_avg_comb = N_avg_comb;
    %  MPD04.RB_comb = RB_comb;
      MPD04.range = range;
      MPD04.time = duration;
      save(name, 'MPD04')
   elseif strcmp(node,'MPD05')==1
      MPD05.N_avg_comb = N_avg_comb;
   %   MPD05.RB_comb = RB_comb;
      MPD05.range = range;
      MPD05.time =  duration;
      save(name, 'MPD05')
  end

 end
  