function[sonde_T_grid, MPD_T_grid, range_grid] = Sonde_read_M2HATS_temp_files_v2(jj, elevation, sondedir, sondefilename, Temp_comb_avg, Temp_comb_var, AH_comb_avg, AH_comb_var, T_lapse, duration, range_grid_size, range_grid_in, sonde_end_int, plot_path, comb_T_surf, comb_P_surf, comb_AH_surf, flag) 

%sondedir
filename = [sondedir sondefilename{jj}]; 
sonde_date = filename(end-17:end-10);
sonde_time = filename(end-8:end-3);
n = datenum([sonde_date sonde_time], 'yyyymmddHHMM');
datestr(n)

ncid = netcdf.open(filename, 'NC_NOWRITE');
  %ncdisp(filename, '/', 'min') % use this to display all variables
  %ncdisp(filename) % use this to display all variables
    
 variable{1} = 'time'; % time (seconds since 1970-1-1)
 variable{2} = 'time_offset'; % time (seconds since 1970-1-1)
 variable{3} = 'alt'; % altitude above mean sea level (m)
 variable{4} = 'pres'; % pressure (hPa or mbar)
 variable{5} = 'tdry'; %Dry Bulb Temperature (C)
 variable{6} = 'rh'; % relative humidity (%)
 variable{7} = 'lat'; % latitude
 variable{8} = 'lon'; % longitude 
% variable{9} = 'reference_alt'; 
  
 base_time  = ncread(filename,variable{1});   
 base_sonde_time = ncread(filename,variable{2}); 
 sonde_alt = ncread(filename,variable{3});
 sonde_P = ncread(filename,variable{4});  
 sonde_T = ncread(filename,variable{5});  
 sonde_RH = ncread(filename,variable{6}); 
 sonde_lat = ncread(filename,variable{7}); 
 sonde_lon = ncread(filename,variable{8}); 
% reference_alt = ncread(filename,variable{9});  
 
netcdf.close(ncid);

% only have the up portion of the profile
   sonde_range_stop = 10000; % set top of sonde at 10,000m
   for i=1:length(sonde_alt)
    if (sonde_range_stop<=sonde_alt(i)) == 1 
      sonde_top = i
      break
    else
      sonde_top = i;
    end
   end
  base_sonde_time =  base_sonde_time(1:sonde_top);
  sonde_alt = sonde_alt(1:sonde_top);
  sonde_P =  sonde_P(1:sonde_top);
  sonde_T =  sonde_T(1:sonde_top);
  sonde_RH = sonde_RH(1:sonde_top);

duration_sonde = n +  double(base_sonde_time)/24/60/60;
sonde_AGL = sonde_alt - elevation;


%figure(101)
%plot(T_sonde,alt)
%figure(102)
%plot(P_sonde,alt)
%figure(103)
%plot(sonde_RH,sonde_alt/1000)

%cd(d) % point back to original directory

%% convert RH to number density and AH
% vapor pressure of water
    a0 = 6.107799961;
    a1 = 4.436518521E-1; 
    a2 = 1.428945805E-2;
    a3 = 2.650648471E-4; 
    a4 = 3.031240396E-6;
    a5 = 2.034080948E-8;
    a6 = 6.136820929E-11;
e=((a0+sonde_T.*(a1+sonde_T.*(a2+sonde_T.*(a3+sonde_T.*(a4+sonde_T.*(a5+sonde_T.*a6))))))./1); %vapor pressure in hPa 

% constants
R = 8.31447215; %J mol^-1 K^-1
N_A= 6.0221415E23; %mol^-1

RH_surf=1;
T_surf=1;

% convert from RH to number density
sonde_N_H2O = 1.*(sonde_RH.*(RH_surf).*e./(R.*(sonde_T+273).*(T_surf))).*N_A*1e-6;  %cm^3
sonde_AH = sonde_N_H2O.*1e6./6.022E23.*18.015;
w_s = 621.9901.*(e./(sonde_P-e));
w = sonde_RH./100.*w_s;
sonde_MR=w;

%flag.plot_overlay = 1;

if flag.plot_overlay == 1
  figure(110)
  plot(sonde_T+273, sonde_AGL/1000)
  xlim([240 320])
  ylim([0 6])

%  figure(111)
%  plot(sonde_MR, sonde_AGL/1000)
%  xlim([0 20])
%  ylim([0 6])
  
  figure(113)
  scatter(duration_sonde, sonde_AGL/1000, 15, sonde_T+273, '+');
  %colormap(jet)
  colormap(parula)
  ylim([0 6])
  caxis([240 320])
  colorbar

  %if flag.MR == 1
  %    % sonde mixing ratio
  %  figure(1)
  %  hold on
  %  scatter(duration_sonde, sonde_AGL/1000, 10, sonde_MR, '+'); %make size 200 if only a single day
  %  ylim([0 6])
  %  colormap(jet)
  %else 
  % Sonde absolute humidity
  figure(1)  %overlay the sondes on the multiday on the next 4 lines
  hold on
   scatter(duration_sonde, sonde_AGL/1000, 15, sonde_T+273, '+');
  ylim([0 6])
  caxis([260 320])
  %colormap(jet)
 % colormap(parula)
  %end
  
  % how far has the sonde moved for the lower 4 km
  figure(104)
  sonde_lat(sonde_alt/1000>4)=NaN;
  sonde_lon(sonde_alt/1000>4)=NaN;
  plot(sonde_lat, sonde_lon)
  hold on
 % plot(36.31, -97.93, 'x') % MPD01 
 % plot(36.88, -97.07, 'x') % MPD02 
 % plot(36.82, -97.82, 'x') % MPD03 
 % plot(36.37, -97.073, 'x') % MPD04 
  hold off
  
  figure(105)
  plot(sonde_T+273.15, sonde_AGL/1000)
  %xlim([0 12])
  ylim([0 6])
  
  
end

% grid sonde data vs range 
%range_grid = 0:range_grid_size/1000:6;
range_grid = range_grid_in(1)/1000:range_grid_size/1000:6;
[sonde_AGL_km, index] = unique(sonde_AGL/1000); 
sonde_T_grid =interp1(sonde_AGL_km, sonde_T(index)+273.15, range_grid, 'linear');
sonde_AH_grid =interp1(sonde_AGL_km, sonde_AH(index), range_grid, 'linear');
sonde_RH_grid =interp1(sonde_AGL_km, sonde_RH(index), range_grid, 'linear');
sonde_P_grid =interp1(sonde_AGL_km, sonde_P(index), range_grid, 'linear');


% find the closes time index for the MPD water vapor
[minValue, closestIndex] = min(abs(min(duration_sonde)-duration))
[minValue, closestIndex_end] = min(abs(min(duration_sonde+sonde_end_int/24/60)-duration))
%MPD_AH = N_avg_comb(closestIndex,:).*1e6./6.022E23.*18.015;
%MPD_AH_var =  comb_AH_var(closestIndex,:);
%MPD_T_lapse = nanmedian(T_lapse(closestIndex:closestIndex_end,:),1)+273.15;
% MPD_T = median(Temp_comb(closestIndex:closestIndex_end,:), 'omitnan');
% MPD_T = median(Temp_comb_avg(closestIndex:closestIndex_end,:), 'omitnan');
% %MPD_T_var = var(Temp_comb(closestIndex:closestIndex_end,:),'includenan');
% MPD_T_var = median(Temp_comb_var(closestIndex:closestIndex_end,:),'includenan');
% MPD_T(isnan(MPD_T_var)) = nan;

% MPD_T = mean(Temp_comb_avg(closestIndex:closestIndex_end,1:round(6000/range_grid_size-1)),1, 'omitnan');
% MPD_T_var =  mean(Temp_comb_var(closestIndex:closestIndex_end,1:round(6000/range_grid_size-1)),1, 'omitnan');
 MPD_T = median(Temp_comb_avg(closestIndex:closestIndex_end,:),1, 'omitnan');
 MPD_T_var =  median(Temp_comb_var(closestIndex:closestIndex_end,:),1, 'omitnan');
 MPD_AH = median(AH_comb_avg(closestIndex:closestIndex_end,:),1, 'omitnan');
 MPD_AH_var =  median(AH_comb_var(closestIndex:closestIndex_end,:),1, 'omitnan');
 MPD_T_surf = median(comb_T_surf(closestIndex:closestIndex_end,:),1, 'omitnan');
 MPD_P_surf = median(comb_P_surf(closestIndex:closestIndex_end,:),1, 'omitnan');
 MPD_AH_surf = median(comb_AH_surf(closestIndex:closestIndex_end,:),1, 'omitnan');

 try
 %  MPD_T_lapse_grid = interp1(range_grid_in/1000, MPD_T_lapse, range_grid, 'linear');  
    MPD_T_grid = interp1(range_grid_in(~isnan(MPD_T))/1000, MPD_T(~isnan(MPD_T)), range_grid, 'linear');
    MPD_T_var_grid = interp1(range_grid_in(~isnan(MPD_T_var))/1000, MPD_T_var(~isnan(MPD_T_var)), range_grid, 'linear');
    MPD_AH_grid = interp1(range_grid_in(~isnan(MPD_AH))/1000, MPD_AH(~isnan(MPD_AH)), range_grid, 'linear');
    MPD_AH_var_grid = interp1(range_grid_in(~isnan(MPD_AH_var))/1000, MPD_AH_var(~isnan(MPD_AH_var)), range_grid, 'linear');
 catch
    MPD_T_grid = MPD_T;
    MPD_T_var_grid = MPD_T_var;
    MPD_AH_grid = MPD_AH;
    MPD_AH_var_grid = MPD_AH_var;
 end

if flag.plot_overlay == 1
  % overlay sonde vs MPD
  figure(115)
  plot(sonde_T_grid, range_grid)
  hold on
  plot(MPD_T_grid, range_grid, 'r')
  plot(MPD_T_grid, range_grid, 'ro')
  %plot the sonde T vs a standard lapse rate and surface station
%  plot(MPD_T_lapse_grid, range_grid, 'g+')
  eb(1) = errorbar(MPD_T_grid, range_grid, sqrt(MPD_T_var_grid), 'horizontal', 'LineStyle', 'none', 'HandleVisibility','off');
  set(eb, 'color', 'r', 'LineWidth', 1)
  hold off
  xlim([260 310])
  ylim([0 6])
  % grid(gca,'minor')
  grid on
  set(gca, 'YMinorTick','on', 'YMinorGrid','on')
  title(datestr(n))
  xlabel('Temperature (K)'); 
  ylabel('Range (km)'); 
  
  figure(116)
  plot(sonde_AH_grid, range_grid)
  hold on
  plot(MPD_AH_grid, range_grid, 'r')
  plot(MPD_AH_grid, range_grid, 'ro')
  eb(1) = errorbar(MPD_AH_grid, range_grid, sqrt(MPD_AH_var_grid), 'horizontal', 'LineStyle', 'none', 'HandleVisibility','off');
  set(eb, 'color', 'r', 'LineWidth', 1)
  hold off
  xlim([0 10])
  ylim([0 6])
  % grid(gca,'minor')
  grid on
  set(gca, 'YMinorTick','on', 'YMinorGrid','on')
  title(datestr(n))
  xlabel('Absolute Humidity (g/m^3)'); 
  ylabel('Range (km)'); 
  
  
  MPD_N_wv = MPD_AH_grid./1e6.*N_A./18.015;
  MPD_T_grid_C = MPD_T_grid-273.15;
  e_grid=((a0+MPD_T_grid_C.*(a1+MPD_T_grid_C.*(a2+MPD_T_grid_C.*(a3+MPD_T_grid_C.*(a4+MPD_T_grid_C.*(a5+MPD_T_grid_C.*a6))))))./1); %vapor pressure in hPa 
  MPD_RH_grid =  (MPD_N_wv.*R.*MPD_T_grid)./e_grid./N_A./1e-6;
 
  figure(117)
  plot(sonde_RH_grid, range_grid)
  hold on
  plot(MPD_RH_grid, range_grid, 'r')
  plot(MPD_RH_grid, range_grid, 'ro')
%   eb(1) = errorbar(MPD_AH_grid, range_grid, sqrt(MPD_AH_var_grid), 'horizontal', 'LineStyle', 'none', 'HandleVisibility','off');
%   set(eb, 'color', 'r', 'LineWidth', 1)
  hold off
  xlim([0 100])
  ylim([0 6])
  % grid(gca,'minor')
  grid on
  set(gca, 'YMinorTick','on', 'YMinorGrid','on')
  title(datestr(n))
  xlabel('Relative Humidity (%)'); 
  ylabel('Range (km)'); 
  
  
   
  
  figure(118)
  t = tiledlayout(1,1)
  ax1 = axes(t);
  plot(sonde_T_grid-273.15, range_grid, 'r')
  hold on
  plot(MPD_T_grid-273.15, range_grid, 'r')
  plot(MPD_T_grid-273.15, range_grid, 'ro')
  hold off
%  xlim([260 310])
  xlim([-40 40])
  ylim([-0.25 5])
  ax1.XAxisLocation = 'origin';
  ax1.YAxisLocation = 'left';
  ax1.XColor = 'r';
  xlabel('Temperature (C)'); 
  ylabel('Range AGL (km)'); 

   
  ax2 = axes(t);
  plot(ax2,sonde_RH_grid, range_grid, 'b')
  hold on
  plot(MPD_RH_grid, range_grid, 'b')
  plot(MPD_RH_grid, range_grid, 'bo')
  hold off
  xlim([0 100])
  ylim([-0.25 5])
  xlabel('Relative Humidity (%)'); 
  ylabel('Range AGL (km)'); 
  ax2.XAxisLocation = 'top';
  ax2.YAxisLocation = 'right';
  ax2.XColor = 'b';
  ax2.Color = 'none';
  title(datestr(n))
    
  ax3 = axes(t);
  plot(ax3,sonde_AH_grid, range_grid, 'g')
  hold on
  plot(MPD_AH_grid, range_grid, 'g')
  plot(MPD_AH_grid, range_grid, 'go')
  hold off
  xlim([0 10])
  ylim([-0.25 5])
  xlabel('Absolute Humidity (g/m^3)'); 
  ylabel('Range AGL (km)'); 
  ax3.XAxisLocation = 'bottom';
  ax3.YAxisLocation = 'left';
  ax3.XColor = 'g';
  ax3.Color = 'none';
  
  ax1.Box = 'off';
  ax2.Box = 'off';
  ax3.Box = 'off';
   
  
   % virtual potential temperature
   const.M_wv = 18.015;  % molar mass of water molecule
   const.M_air = 28.97; % molar mass gm/mol of air
   const.N_A = 6.022E23; % Avagadro number
   const.k_B = 1.3806488e-23; % (J/K)
   % water vapor number density 
   n_wv = MPD_AH_grid.*const.N_A./const.M_wv;
   n_wv_surf = MPD_AH_surf.*const.N_A./const.M_wv;
   n_wv_sonde = sonde_AH_grid.*const.N_A./const.M_wv;
   % water vapor partial pressure (convert from Pa to atm)
   e_wv =  n_wv.*const.k_B.*(MPD_T_grid)./101300;
   e_wv_surf =  n_wv_surf.*const.k_B.*(MPD_T_surf)./101300;
   sonde_e_wv = n_wv_sonde.*const.k_B.*(sonde_T_grid)./101300;
   % take sonde pressure in hPa and convert to atm
   MPD_P_grid = sonde_P_grid./1.01300e3;
   % calculate virtual temperature (assumes pressure in atm)
   MPD_T_v = MPD_T_grid./(1-((e_wv)./MPD_P_grid).*(1-const.M_wv/const.M_air));
   MPD_T_v_surf = MPD_T_surf./(1-((e_wv_surf)./MPD_P_surf).*(1-const.M_wv/const.M_air));
   sonde_T_v = sonde_T_grid./(1-((sonde_e_wv)./(sonde_P_grid./1.01300e3)).*(1-const.M_wv/const.M_air));
   % calculate virtual potential temperature
   P_0 = 0.987; % reference pressure in atm
   MPD_T_p = MPD_T_v.*(P_0./MPD_P_grid).^0.286; 
   MPD_T_p_surf = MPD_T_v_surf.*(P_0./MPD_P_surf).^0.286;
   sonde_T_p = sonde_T_v.*(P_0./(sonde_P_grid./1.01300e3)).^0.286; 
  
  if  isnan(sonde_T_p(1)) 
    surf_T =  sonde_T_p(2); 
  else
   surf_T =  sonde_T_p(1); 
  end
 
  figure(119)
  
  plot(sonde_T_p, range_grid, 'k')
  hold on
  plot(MPD_T_p, range_grid, 'r')
  plot(MPD_T_p_surf(1:5:end)*ones(size(sonde_T_grid(1:5:end))), range_grid(1:5:end), 'r.')
  plot(surf_T(1:5:end)*ones(size(sonde_T_grid(1:5:end))), range_grid(1:5:end), 'k.')
  hold off
  %xlim([0 100])
  ylim([0 5])
  % grid(gca,'minor')
  grid on
  set(gca, 'YMinorTick','on', 'YMinorGrid','on')
  title(datestr(n))
  xlabel('Virtual Potential Temp (K)'); 
  ylabel('Range (km)'); 
   
   
  
  Scrsize=[1 1 400 400];
  %cd('/Users/lroot/Desktop/mpd/Plots/')
  %cd('/Volumes/documents/mpd/Plots/')
  cd(plot_path)
  
   FigH = figure(115);
 %  set(gca,'Fontsize',30,'Fontweight','b'); % 
   set(FigH, 'PaperUnits', 'points', 'PaperPosition', Scrsize);
   name=strcat(sonde_date, '_', sonde_time, 'Sonde_temp_profile'); 
   print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
  
  FigH = figure(118);
 % set(gca,'Fontsize',30,'Fontweight','b'); % 
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', Scrsize);
  name=strcat(sonde_date, '_', sonde_time, 'Sonde_MPD_profile'); 
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
  
  FigH = figure(119);
 % set(gca,'Fontsize',30,'Fontweight','b'); % 
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', Scrsize);
  name=strcat(sonde_date, '_', sonde_time, 'Virtual_Potential_T_profile'); 
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
  
  
 % save(name, 'range_grid', 'sonde_AH_grid', 'MPD_AH_grid', 'MPD_AH_var_grid')
  
  
  
end



%pause

