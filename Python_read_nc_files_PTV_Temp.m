
clear all; close all
% day = datenum('23 Jul 2023'); % automated to pick the first day


skip = 1
node = 'MPD04';
flag.PTV = 1 %use the PTV data which has temperature
offset = 2; % parcel method surface temp offset
gradient_limit = 0.2; % mask postive gradients
addpath '/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing/matplotlib/';

serv_path = '/Volumes/eol/sci/mhayman';
cd(strcat(serv_path,'/DIAL/Processed_Data/M2HATS/5min_release1.0/qc_masked/'))
%cd(strcat(serv_path,'/DIAL/Processed_Data/ECLIPSE/initial_test/'))


if flag.PTV == 1
%  serv_path = '/Volumes/eol/smaug1/rsfdata/MPD';
%  cd(strcat(serv_path,'/mpd_03_processed_data/PTV/temperature/M2HATS1.0'))
%  serv_path = '/Volumes/eol/sci/mhayman';
%  cd(strcat(serv_path,'/DIAL/Processed_Data/Marshall_2024/ptv1.0/'))
%  serv_path = '/Volumes/eol/scr-tmp/mhayman';
  %cd(strcat(serv_path,'/O2DIAL/PTV/gp_search/'))
  serv_path = '/Volumes/eol/sci/mhayman';
  cd(strcat(serv_path,'/DIAL/Processed_Data/BRIDGE_2025/ptv0.0/'))
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
 variable{6} = 'Absolute_Humidity_Standard'; 
 variable{7} = 'Absolute_Humidity_Standard_mask';
 variable{8} = 'Absolute_Humidity_Standard_uncertainty'; 
 variable{6} = 'Absolute_Humidity'; 
 variable{7} = 'Absolute_Humidity_mask';
 variable{8} = 'Absolute_Humidity_uncertainty'; 
 variable{9} = 'Aerosol_Backscatter_Coefficient';
 variable{10} = 'Aerosol_Backscatter_Coefficient_mask';
 variable{11} = 'Aerosol_Backscatter_Coefficient_variance';
 
 variable{9} = 'Backscatter_Photon_Counts_828'; 
 
 if flag.PTV == 1
   
   % variable{6} = 'Water_Vapor'; 
   % variable{8} = 'Water_Vapor_std';  
   variable{11} = 'Aerosol_Backscatter_Coefficient_uncertainty';
   % variable{9} = 'Backscatter_Ratio';
   % variable{10} = 'Backscatter_Ratio_mask';
   % variable{11} = 'Backscatter_Ratio_uncertainty';
   % variable{9} = 'Relative_Backscatter';
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
   variable{24} = 'Temperature_Model';
   variable{25} = 'Temperature_Standard'; 
   variable{26} = 'Absolute_Humidity_Standard'; 
   variable{27} = 'Absolute_Humidity_Standard_mask'; 
 end

    

for jj = 1:size(Pythonfilename,2)
  filename = Pythonfilename{jj};
  date = filename(end-15:end-8);
  %n = datenum(date, 'yymmdd');
  %date = filename(end-54:end-47);
  % %date = filename(end-15:end-8);
  % expression = '\d{8}T\d{4}'; % Matches 8 digits, 'T', then 4 digits
  % matches = regexp(filename, expression, 'match');
  % first_timestamp_str = matches{1};
  % date = first_timestamp_str(1:8); 


  n = datenum(date, 'yyyymmdd');
  
  ncid = netcdf.open(filename, 'NC_NOWRITE');
    %ncdisp(filename, '/', 'min') % use this to display all variables
    ncdisp(filename) % use this to display all variables
    time{jj} = ncread(filename,variable{1}); 
    alt{jj} = ncread(filename,variable{2});
    
if flag.PTV == 1
     T{jj}  = ncread(filename,variable{3});  
%     T_model{jj}  = ncread(filename,variable{24});
%     T_stand{jj}  = ncread(filename,variable{25});  
     T_mask{jj} = ncread(filename,variable{4}); 
 %    T_var{jj} = ncread(filename,variable{5}); 
     T{jj}(T_mask{jj} == 1) = nan;
     T_var{jj}(T_mask{jj} == 1) = nan; 
 %    P{jj} =  ncread(filename,variable{12});
 %    P_mask{jj} = ncread(filename,variable{13});
 %    L{jj} =  ncread(filename,variable{14});
 %    L_mask{jj} = ncread(filename,variable{15});
 %    PBLH{jj} =  ncread(filename,variable{16});
 %    T_vp{jj} =  ncread(filename,variable{17});
 %    T_vp_mask{jj} = ncread(filename,variable{18});
 %    Lift{jj} =  ncread(filename,variable{19});
 %    Lift_mask{jj} = ncread(filename,variable{20});
     T_surf{jj} =  ncread(filename,variable{21});
     P_surf{jj} =  ncread(filename,variable{22});
 %    AH_surf{jj} =  ncread(filename,variable{23});
   
end
    
    AH{jj}  = ncread(filename,variable{6});  
    AH_mask{jj} = ncread(filename,variable{7}); 
    AH_var{jj} = ncread(filename,variable{8}); 
    AH{jj}(AH_mask{jj} == 1) = nan;
    AH_var{jj}(AH_mask{jj} == 1) = nan; 

%    AH_stand{jj} =  ncread(filename,variable{26}); 
%    AH_stand_mask{jj} =  ncread(filename,variable{27}); 
%    AH_stand{jj}(AH_stand_mask{jj} == 1) = nan;
    % AH{jj}(AH_var{jj} > 3) = nan;
       
    ABC{jj}  = ncread(filename,variable{9});   
%     ABC_mask{jj} = ncread(filename,variable{10}); 
%     ABC_var{jj} = ncread(filename,variable{11});
  %   ABC{jj}(ABC_mask{jj} == 1) = nan;
   %  ABC_var{jj}(ABC_mask{jj} == 1) = nan; 
   % % ABC{jj}(ABC_var{jj} > 1e-11) = nan;
   %  ABC{jj}(ABC{jj} < 1e-12) = nan;
    
    if flag.PTV == 1
      % P{jj}(P_mask{jj} == 1) = nan;
%      L{jj}(L_mask{jj} == 1) = nan;
%      T_vp{jj}(T_vp_mask{jj} == 1) = nan;
%      Lift{jj}(Lift_mask{jj} == 1) = nan;
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
 %     comb_ABC_var = ABC_var{jj};
      if flag.PTV == 1
        comb_T = T{jj};
%        comb_T_model = T_model{jj};
 %       comb_T_stand = T_stand{jj};
        comb_T_var = T_var{jj};
 %       comb_P = P{jj};
 %       comb_L = L{jj};
 %       comb_PBLH = PBLH{jj};
 %       comb_T_vp = T_vp{jj};  
 %       comb_Lift = Lift{jj};  
        comb_T_surf = T_surf{jj}; 
        comb_P_surf = P_surf{jj}; 
 %       comb_AH_surf = AH_surf{jj};
  %      comb_AH_stand = AH_stand{jj}; 
      end
  else
      comb_duration = [comb_duration; duration{jj}];
      % find the maximum range to accumulate (in case it changes)
      % max_range = min(cellfun('size',alt,1))
      comb_AH = [comb_AH AH{jj}];
      comb_AH_var = [comb_AH_var AH_var{jj}];
      comb_ABC = [comb_ABC ABC{jj}];
 %     comb_ABC_var = [comb_ABC_var ABC_var{jj}];
      if flag.PTV == 1
        comb_T = [comb_T T{jj}];
  %      comb_T_model = [comb_T_model T_model{jj}];
  %      comb_T_stand = [comb_T_stand T_stand{jj}];        
  %      comb_AH_stand = [comb_AH_stand AH_stand{jj}];
  %      comb_T_var = [comb_T_var T_var{jj}];
  %      comb_P = [comb_P P{jj}];
  %      comb_L = [comb_L L{jj}];
  %      comb_PBLH = [comb_PBLH; PBLH{jj}];
  %      comb_T_vp = [comb_T_vp T_vp{jj}];
  %      comb_Lift = [comb_Lift; Lift{jj}];
        comb_T_surf = [comb_T_surf; T_surf{jj}];
        comb_P_surf = [comb_P_surf; P_surf{jj}];
  %      comb_AH_surf = [comb_AH_surf; AH_surf{jj}];
      end
  end
  
end

scrsz = get(0,'ScreenSize');
date=datestr(n, 'yyyy-mmm-dd');
plot_size1 = [scrsz(4)/1.5 scrsz(4)/10 scrsz(3)/1.5 scrsz(4)/3];
font_size = 16;

x = comb_duration;
y = (alt{1}./1000);
xData =  linspace( fix(min(x)),  ceil(max(x)), round((ceil(max(x))-fix(min(x)))/skip)+1 );
xData_m =  linspace( fix(min(x)),  ceil(max(x)), round((ceil(max(x))-fix(min(x)))/skip*24/2)+1 );


if flag.PTV == 1
% % plot the T
  Z = real(comb_T-273.15);
  figure3 = figure('Position',plot_size1);
  set(gcf,'renderer','zbuffer');
  h = pcolor(x, y, Z);
  set(h, 'EdgeColor', 'none'); 
  axis xy; colorbar('EastOutside'); 
  caxis([-20 30]);
  axis([fix(min(x)) ceil(max(x)) 0 6]) 
%  shading interp
  set(gca, 'XTick',  xData)
  set(gca,'XMinorTick','on')
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  hh = title({[node, ' Temp (C)']},...
       'fontweight','b','fontsize',font_size);  
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
 % colormap(jet)
  colormap(plasma)
  set(gca,'Fontsize',font_size,'Fontweight','b');
  
% %  % % plot the P or Potenial Temperature
% % virtual potential temperature
%   const.M_wv = 18.015;  %molar mass of water molecule
%   const.M_air = 28.97; % molar mass gm/mol of air
%   const.N_A = 6.022E23; % Avagadro number
%   const.k_B = 1.3806488e-23; % (J/K)
%   % water vapor number density 
%   n_wv = comb_AH.*const.N_A./const.M_wv;
%   n_wv_surf = comb_AH_surf.*const.N_A./const.M_wv;
%   % water vapor partial pressure (Pa to atm convert)
%   e_wv =  n_wv.*const.k_B.*comb_T./101300;
%   e_wv_surf =  n_wv_surf.*const.k_B.*(comb_T_surf)./101300;
%   % calculate virtual temperature
%   comb_T_v = comb_T./(1-((e_wv)./comb_P).*(1-const.M_wv/const.M_air));
%   comb_T_v_surf = comb_T_surf./(1-((e_wv_surf)./comb_P_surf).*(1-const.M_wv/const.M_air));
%   % calculate virtual potential temperature
%   P_0 = 0.987; % reference pressure in atm
%   comb_T_p = comb_T_v.*(P_0./comb_P).^0.286; 
%   comb_T_p_surf = comb_T_v_surf.*(P_0./comb_P_surf).^0.286;

  
%    % add in the surface T 
%    y(1)= 0;
%    y(2)= 0.1;
%    comb_T_p(1,:) = comb_T_p_surf; 
  
%   B = repmat(comb_T_p_surf, 1, size(comb_T_p,1));
%   B_test = B'+offset; % offset virtual potential temp at surface by 2deg
%   PBLH = comb_T_p;
%   PBLH (PBLH >B_test) = nan;
%   %[FX,FY] = gradient(comb_T_p);
%   %comb_T_p(FY>gradient_limit) = nan;
%   B_test_2 = ~isnan(PBLH);
%   Indices = arrayfun(@(x) find(B_test_2(:,x), 1, 'last'), 1:size(PBLH, 2), 'UniformOutput', false);
%   idx = cellfun(@isempty,Indices);
%   Indices(idx) = {1};
%   BLH_Indices = cell2mat(Indices);
%   Thermo_BLH = y(BLH_Indices); 
% % % set this to just plot the first day selected
% %   day = datenum( Pythonfilename{1}(end-15:end-8), 'yyyymmdd');

 

 %   %Z = comb_T-comb_T_model;
 %   Z = comb_T_model-273.115;
 % 
 %   figure4 = figure('Position',plot_size1);
 %   set(gcf,'renderer','zbuffer');
 %   h = pcolor(x, y, Z);
 %   set(h, 'EdgeColor', 'none'); 
 %   axis xy; colorbar('EastOutside'); 
 %   caxis([-20 20]);
 %   axis([fix(min(x)) ceil(max(x)) 0 6]) 
 %   set(gca, 'XTick',  xData)
 %   set(gca,'XMinorTick','on')
 %   xAx = get(gca,'XAxis');
 %   xAx.MinorTickValues=xData_m;
 %    set(gca,'TickDir','out');
 %   set(gca,'TickLength',[0.005; 0.0025]);
 %   hh = title({[node, ' Temp HRRR(C)']},...
 %        'fontweight','b','fontsize',font_size);  
 %   ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
 %   datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
 %   colormap(jet)
 % %  colormap(redblue)
 %   colormap(plasma)
 %   set(gca,'Fontsize',font_size,'Fontweight','b'); 
 %   grid on 
 %   grid(gca, 'minor')
% 
%     
%      % add the overlay of PBLH
%   ax2 = axes(figure4);
%   %plot(x, comb_PBLH/1000, 'ro')
%   plot(x, Thermo_BLH, 'go')
%   colorbar('EastOutside'); 
%   axis([fix(min(x)) ceil(max(x)) 0 6])
%   ax2.Color = 'none';
%   ax2.XAxisLocation = 'top';
%   ax2.YAxisLocation = 'right';
%   datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
%   set(gca,'Fontsize',font_size,'Fontweight','b');
%   ax2.YColor = 'none';
%   ax2.XColor = 'none';
%   % Kludge to remove the colorbar font, but not rescale the axis
%   c = colorbar;
%   c.FontSize = 0.01;

   
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
  
  %  Z = comb_T-comb_T_model;
  %  figure5 = figure('Position',plot_size1);
  %  set(gcf,'renderer','zbuffer');
  %  h = pcolor(x, y, Z);
  %  set(h, 'EdgeColor', 'none'); 
  %  axis xy; colorbar('EastOutside'); 
  %  caxis([-10 10]);
  %  axis([fix(min(x)) ceil(max(x)) 0 6]) 
  %  set(gca, 'XTick',  xData)
  %  set(gca,'XMinorTick','on')
  %  xAx = get(gca,'XAxis');
  %  xAx.MinorTickValues=xData_m;
  %   set(gca,'TickDir','out');
  %  set(gca,'TickLength',[0.005; 0.0025]);
  %  hh = title({[node, ' T_{PTV} - T_{HRRR} (K)']},...
  %       'fontweight','b','fontsize',font_size);  
  %  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
  %  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
  % colormap(redblue)
  %  set(gca,'Fontsize',font_size,'Fontweight','b'); 
  %  grid on 
  %  grid(gca, 'minor')


end


% plot the AH
  Z = real(comb_AH);
  figure1 = figure('Position',plot_size1);
  set(gcf,'renderer','zbuffer');
  h = pcolor(x, y, Z);
  set(h, 'EdgeColor', 'none'); 
  axis xy; 
  colorbar('EastOutside'); 
  caxis([0 25]);
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

  
%  % add the overlay of PBLH
%   ax2 = axes(figure1);
%   plot(x, comb_PBLH/1000, 'ro')
%   colorbar('EastOutside'); 
%   axis([fix(min(x)) ceil(max(x)) 0 6])
%   ax2.Color = 'none';
%   ax2.XAxisLocation = 'top';
%   ax2.YAxisLocation = 'right';
%   datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
%   set(gca,'Fontsize',font_size,'Fontweight','b');
%   ax2.YColor = 'none';
%   ax2.XColor = 'none';
%   % Kludge to remove the colorbar font, but not rescale the axis
%   c = colorbar;
%   c.FontSize = 0.01;

  
  

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

  caxis([1e2 1e6]);
    hh = title({[node, ' ', ' Relative Backscatter]']},'fontweight','b','fontsize',font_size);  
 % caxis([1e-8 1e-6]);
 %   hh = title({[node, ' Aerosol Backscatter Coefficient m^{-1} sr^{-1}']},'fontweight','b','fontsize',font_size);     
 % hh = title({[node, ' Backscatter Ratio']},...
 %     'fontweight','b','fontsize',font_size);     
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
 % datetick('x','HH:MM');
 % colormap(jet)
  colormap(viridis)
  set(gca,'Fontsize',font_size,'Fontweight','b');
  set(gca,'Zscale', 'log')
  set(gca,'Colorscale', 'log')
  set(gca,'Zscale', 'linear')

  
%   % add the overlay of PBLH
%   ax2 = axes(figure2);
%   plot(x, comb_PBLH/1000, 'ro')
%   colorbar('EastOutside'); 
%   axis([fix(min(x)) ceil(max(x)) 0 6])
%   ax2.Color = 'none';
%   ax2.XAxisLocation = 'top';
%   ax2.YAxisLocation = 'right';
%   datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
%   set(gca,'Fontsize',font_size,'Fontweight','b');
%   ax2.YColor = 'none';
%   ax2.XColor = 'none';
%   % Kludge to remove the colorbar font, but not rescale the axis
%   c = colorbar;
%   c.FontSize = 0.01;

 %  % plot the standard temp 
 %   Z = real(comb_T_std-273.15);
 %   figure7 = figure('Position',plot_size1);
 %   set(gcf,'renderer','zbuffer');
 %   h = pcolor(x, y, Z);
 %   set(h, 'EdgeColor', 'none'); 
 %   axis xy; colorbar('EastOutside'); 
 %   caxis([-20 20]);
 %   axis([fix(min(x)) ceil(max(x)) 0 6]) 
 %   set(gca, 'XTick',  xData)
 %   set(gca,'XMinorTick','on')
 %   xAx = get(gca,'XAxis');
 %   xAx.MinorTickValues=xData_m;
 %    set(gca,'TickDir','out');
 %   set(gca,'TickLength',[0.005; 0.0025]);
 %   hh = title({[node, ' T Standard (C)']},...
 %        'fontweight','b','fontsize',font_size);  
 %   ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
 %   datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
 %   colormap(jet)
 % %  colormap(redblue)
 %   colormap(plasma)
 %   set(gca,'Fontsize',font_size,'Fontweight','b'); 
 %   grid on 
 %   grid(gca, 'minor')
 % 
 % 
 %   Z = comb_T_std-comb_T_model;
 %   figure7 = figure('Position',plot_size1);
 %   set(gcf,'renderer','zbuffer');
 %   h = pcolor(x, y, Z);
 %   set(h, 'EdgeColor', 'none'); 
 %   axis xy; colorbar('EastOutside'); 
 %   caxis([-10 10]);
 %   axis([fix(min(x)) ceil(max(x)) 0 6]) 
 %   set(gca, 'XTick',  xData)
 %   set(gca,'XMinorTick','on')
 %   xAx = get(gca,'XAxis');
 %   xAx.MinorTickValues=xData_m;
 %    set(gca,'TickDir','out');
 %   set(gca,'TickLength',[0.005; 0.0025]);
 %   hh = title({[node, ' T_{standard} - T_{HRRR} (K)']},...
 %        'fontweight','b','fontsize',font_size);  
 %   ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
 %   datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
 %  colormap(redblue)
 %   set(gca,'Fontsize',font_size,'Fontweight','b'); 
 %   grid on 
 %   grid(gca, 'minor')
   

% figure(101)
%  semilogy(comb_duration, comb_ABC(8,:))
% hold on
%  semilogy(comb_duration, comb_ABC(9,:))
%  semilogy(comb_duration, comb_ABC(10,:))
% hold off
%   datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
%   title('MPD03 Aerosol Backscatter Coefficient m^{-1} sr^{-1}')
%  legend(num2str(y(8)), num2str(y(9)), num2str(y(10)))



   %MultiPlot
   
   figure10 = figure('Position',[1 1 1920 250*4]);
   font_size = 16;
   subplot1=subplot(3,1,1,'Parent',figure10,'YGrid','on', 'XGrid','on');
   box(subplot1,'on');
  % hold(subplot1,'all');
      Z = real(comb_ABC);
      set(gcf,'renderer','zbuffer');
      h = pcolor(x,y,Z);
      set(h, 'EdgeColor', 'none');
      set(gca,'TickDir','out');
      set(gca,'TickLength',[0.005; 0.0025]);  
      set(gca, 'XTick',  xData) 
      colorbar('EastOutside');
      axis([fix(min(x))  ceil(max(x)) 0 6])
      % caxis([1 10]);
      % hh = title({[node, ' ', ' Backscatter Ratio']},'fontweight','b','fontsize',font_size);  
      % caxis([1e-8 1e-3]);
      % hh = title({[node, ' ', ' Backscatter Coefficient [m^{-1}sr^{-1}]']},'fontweight','b','fontsize',font_size);  
      caxis([1e2 1e6]);
      hh = title({[node, ' ', ' Relative Backscatter ']},'fontweight','b','fontsize',font_size);   
      datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
      ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size);
      set(gca,'Fontsize',font_size,'Fontweight','b');
      %set(gca,'Zscale', 'log')
      set(gca,'Colorscale', 'log')
      set(gca,'Zscale', 'linear')
      colormap(subplot1, jet)
      grid on
  subplot2=subplot(3,1,2,'Parent',figure10,'YGrid','on', 'XGrid','on');
  box(subplot2,'on');
      Z = real(comb_AH);  
     % Z = real(comb_AH_stand); 
      set(gcf,'renderer','zbuffer');
      h = pcolor(x,y,Z);
      set(h, 'EdgeColor', 'none');
      colorbar('EastOutside');
      axis([fix(min(x)) ceil(max(x)) 0 6])
      caxis([0 25]);
      colormap(subplot2, CM_YlGnBu(64))
      ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
      set(gca, 'XTick',  xData)
      set(gca,'TickDir','out');
      set(gca,'TickLength',[0.005; 0.0025]);
      hh = title({[node, ' Water Vapor PTV (g m^{-3})']},'fontweight','b','fontsize',font_size);
      hh = title({[node, ' Water Vapor Standard (g m^{-3})']},'fontweight','b','fontsize',font_size);
      datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
      set(gca,'Fontsize',font_size,'Fontweight','b');
      grid on
  subplot3=subplot(3,1,3,'Parent',figure10,'YGrid','on', 'XGrid','on');
  box(subplot3,'on'); 
      Z = real(comb_T-273.15);
      set(gcf,'renderer','zbuffer');
      h = pcolor(x,y,Z);
      set(h, 'EdgeColor', 'none');
      set(gca,'TickDir','out');
      set(gca,'TickLength',[0.005; 0.0025]);
      set(gca, 'XTick',  xData)
      colorbar('EastOutside');
      axis([fix(min(x))  ceil(max(x)) 0 6])
      caxis([-20 30]);
      hh = title({[node, ' ', ' Temp PTV (C)']},'fontweight','b','fontsize',font_size);  
      datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
      ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size);
      set(gca,'Fontsize',font_size,'Fontweight','b');
      colormap(subplot3, plasma)
  % subplot4=subplot(4,1,4,'Parent',figure10,'YGrid','on', 'XGrid','on');
  % box(subplot4,'on'); 
  %     Z = real(comb_T_stand-273.15);    
  %     set(gcf,'renderer','zbuffer');
  %     h = pcolor(x,y,Z);
  %     set(h, 'EdgeColor', 'none');
  %     set(gca, 'XTick',  xData)
  %     set(gca,'TickDir','out');
  %     set(gca,'TickLength',[0.005; 0.0025]);
  %     colorbar('EastOutside');
  %     axis([fix(min(x)) ceil(max(x)) 0 6])
  %     caxis([-20 30]);
  %     hh = title({[node, ' Temp Standard (C.)']},'fontweight','b','fontsize',font_size);
  %     datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
  %     ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size);
  %     set(gca,'Fontsize',font_size,'Fontweight','b');
  %     colormap(subplot4, plasma)
  %    % set(gca,'Fontsize',font_size,'Fontweight','b');
  %    % set(gca,'Zscale', 'log')
  %    % set(gca,'Colorscale', 'log')
  %    % set(gca,'Zscale', 'linear')
   

 % bin = 3;    
 % range_plot = fix(alt{1}(bin))
 % 
 % 
 % figure(111)
 % plot(comb_duration, comb_AH_stand(bin,:))
 % hold on
 % plot(comb_duration, comb_AH(bin,:))
 % hold off            
 % ylabel('Absoultue Humidity (g m^{-3})');
 % title_string = sprintf('Absolute Humidity at %d m', range_plot);
 % title(title_string)
 % datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
 % legend ('standard AH', 'PTV AH')
 % grid on
 % grid minor
 % ylim([10 25])
 % 
 % figure(112)
 % plot(comb_duration, comb_T_stand(bin,:))
 % hold on
 % plot(comb_duration, comb_T(bin,:))
 % hold off            
 % ylabel('Temperature (K)');
 % title_string = sprintf('Temperature at %d m', range_plot);
 % title(title_string)
 % datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
 % legend ('standard T', 'PTV T')
 % grid on
 % grid minor
 % ylim([280 310])




 cd(strcat(plot_path,'/mpd/Plots'))

  FigH = figure(1);
  set(gca,'Fontsize',16,'Fontweight','b'); 
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 3000 250]);
  name=strcat(date, node, ' T_Python_multi'); 
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution
  
%   FigH = figure(2);
%   set(gca,'Fontsize',16,'Fontweight','b'); 
%   set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920 250]);
%   name=strcat(date, node, ' Back_Coeff_Python_multi'); 
%   print(FigH, name, '-dpng', '-r0') % set at the screen resolution
  
 if flag.PTV == 1
  FigH = figure(1);
  FigH.Units = 'pixels'; % Ensure units are pixels for direct mapping to your old 'PaperPosition' width/height
  FigH.Position = [100 100 1920 300]; % x, y, width, height in pixels
  name=char(strcat(node, "_", date, 'T_Python_multi')); 
  exportgraphics(FigH, [name, '.png'], 'Resolution', 150);
  
 % FigH = figure(2);
 %  set(gca,'Fontsize',16,'Fontweight','b'); 
 %  set(FigH, 'PaperUnits', 'points', 'PaperPosition',  [1 1 3000 250]);
 %  name=strcat(date, node, ' T_HRRR'); 
 %  print(FigH, name, '-dpng', '-r0') % set at the screen resolution
  
  % FigH = figure(3);
  % set(gca,'Fontsize',16,'Fontweight','b'); 
  % set(FigH, 'PaperUnits', 'points', 'PaperPosition',  [1 1 3000 250]);
  % name=strcat(date, node, ' T_global-T_HRRR'); 
  % print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
  
  % FigH = figure(4);
  % FigH.Units = 'pixels'; % Ensure units are pixels for direct mapping to your old 'PaperPosition' width/height
  % FigH.Position = [100 100 1920 300]; % x, y, width, height in pixels
  % name=char(strcat(node, "_", date, '_WV_Python_multi')); 
  % exportgraphics(FigH, [name, '.png'], 'Resolution', 150);


%   FigH = figure(5);
%   set(gca,'Fontsize',16,'Fontweight','b'); 
%   set(FigH, 'PaperUnits', 'points', 'PaperPosition',  [1 1 1920 250]);
%   name=strcat(date, node, ' Back_Coeff_Python_multi'); 
%   print(FigH, name, '-dpng', '-r0') % set at the screen resolution
  
%   FigH = figure(6);
%   set(gca,'Fontsize',16,'Fontweight','b'); 
%   set(FigH, 'PaperUnits', 'points', 'PaperPosition',  [1 1 1920 250]);
%   name=strcat(date, node, ' T_standard_multi'); 
%   print(FigH, name, '-dpng', '-r0') % set at the screen resolution
%   
    
  % FigH = figure(7);
  % set(gca,'Fontsize',16,'Fontweight','b'); 
  % set(FigH, 'PaperUnits', 'points', 'PaperPosition',  [1 1 3000 250]);
  % name=strcat(date, node, ' T_standard_minus_HRRR'); 
  % print(FigH, name, '-dpng', '-r0') % set at the screen resolution
  
   FigH = figure(4);
   FigH.Units = 'pixels'; % Ensure units are pixels for direct mapping to your old 'PaperPosition' width/height
   FigH.Position = [0 0 1920 275*4]; % x, y, width, height in pixels
   name=char(strcat(node, "_", date, '_all_comb')); 
   exportgraphics(FigH, [name, '.png'], 'Resolution', 150);

   % 
   % FigH = figure(111);
   % FigH.Units = 'pixels'; % Ensure units are pixels for direct mapping to your old 'PaperPosition' width/height
   % FigH.Position = [0 0 1920/1.5 275/1.5]; % x, y, width, height in pixels
   % name=char(strcat(node, "_", date, '_AH_timeseries')); 
   % exportgraphics(FigH, [name, '.png'], 'Resolution', 75);
   % 
   % FigH = figure(112);
   % FigH.Units = 'pixels'; % Ensure units are pixels for direct mapping to your old 'PaperPosition' width/height
   % FigH.Position = [0 0 1920/1.5 275/1.5]; % x, y, width, height in pixels
   % name=char(strcat(node, "_", date, '_temp_timeseries')); 
   % exportgraphics(FigH, [name, '.png'], 'Resolution', 75);


  
 end  
  

 if flag.save_data == 1
  
  cd(d_save_data);     
  range = alt{1}'; 
  N_avg_comb = (comb_AH./1e6.*6.022E23./18.015);
  duration = duration;


  
  
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
  