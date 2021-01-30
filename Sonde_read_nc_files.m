function[sonde_AH_grid, MPD_AH_grid, range_grid] = Sonde_read_nc_files(jj, elevation, sondedir, sondefilename, N_avg_comb, duration, range_grid_size, range_grid_in, flag) 

%sondedir

filename = [sondedir sondefilename{jj}]; 
%filename = '/scr/sci/tammy/mpd/sgp/soundings/sgpsondewnpnC1.b1.20190429.023100.cdf';
%filename = '/Users/spuler/downloads/sgpsondewnpnC1.b1.20190501.083100.cdf';
%filename = '/Users/spuler/Desktop/mpd_03_processed_data/Sondes/Marshall_Field_Site_20201016_163106.nc';
date = filename(end-15:end-10);
n = datenum(date, 'yymmdd');

ncid = netcdf.open(filename, 'NC_NOWRITE');
  %ncdisp(filename, '/', 'min') % use this to display all variables
  %ncdisp(filename) % use this to display all variables
    
 variable{1} = 'base_time'; % time (seconds since 1970-1-1)
 variable{2} = 'time_offset'; % time (seconds since 1970-1-1)
 variable{3} = 'alt'; % altitude above mean sea level (m)
 variable{4} = 'pres'; % pressure (hPa or mbar)
 variable{5} = 'tdry'; %Dry Bulb Temperature (C)
 variable{6} = 'rh'; % relative humidity (%)
 variable{7} = 'lat'; % latitude
 variable{8} = 'lon'; % longitude 
  
 base_time  = ncread(filename,variable{1});   
 sonde_time = ncread(filename,variable{2});   
 sonde_alt = ncread(filename,variable{3});
 sonde_P = ncread(filename,variable{4});  
 sonde_T = ncread(filename,variable{5});  
 sonde_RH = ncread(filename,variable{6}); 
 sonde_lat = ncread(filename,variable{7}); 
 sonde_lon = ncread(filename,variable{8}); 
 
netcdf.close(ncid);

%convert sonde from Unix time to date number (days since Jan 0 0000) 
duration_sonde = datenum(datetime(int64(sonde_time) + int64(base_time), 'convertfrom', 'posixtime'));
%sonde_offset_min = 30;
%duration_sonde = duration_sonde + sonde_offset_min/24/60;
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
  plot(sonde_AH, sonde_AGL/1000)
  xlim([0 20])
  ylim([0 6])

  figure(111)
  plot(sonde_MR, sonde_AGL/1000)
  xlim([0 20])
  ylim([0 6])
  

  figure(113)
  scatter(duration_sonde, sonde_AGL/1000, 25, sonde_AH, '+');
  colormap(jet)
  ylim([0 6])
  caxis([0 16])
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
   scatter(duration_sonde, sonde_AGL/1000, 10, sonde_AH, '+');
  ylim([0 6])
  colormap(jet)
  %end
  
  % how far has the sonde moved for the lower 4 km
  figure(104)
  sonde_lat(sonde_alt/1000>4)=NaN;
  sonde_lon(sonde_alt/1000>4)=NaN;
  plot(sonde_lat, sonde_lon)
  hold on
  plot(36.31, -97.93, 'x') % MPD01 
  plot(36.88, -97.07, 'x') % MPD02 
  plot(36.82, -97.82, 'x') % MPD03 
  plot(36.37, -97.073, 'x') % MPD04 
  hold off

end

% grid sonde data vs range 
range_grid = 0:range_grid_size/1000:5.9;
[sonde_AGL_km, index] = unique(sonde_AGL/1000); 
sonde_AH_grid =interp1(sonde_AGL_km, sonde_AH(index), range_grid, 'nearest');
% find the closes time index for the MPD water vapor
[minValue, closestIndex] = min(abs(min(duration_sonde)-duration))
MPD_AH = N_avg_comb(closestIndex,:).*1e6./6.022E23.*18.015;
try
 MPD_AH_grid = interp1(range_grid_in(~isnan(MPD_AH))/1000, MPD_AH(~isnan(MPD_AH)), range_grid, 'nearest');
catch
 MPD_AH_grid = interp1(range_grid_in/1000, MPD_AH, range_grid, 'nearest');
end

if flag.plot_overlay == 1
  % overlay sonde vs MPD
  figure(115)
  plot(sonde_AH_grid, range_grid)
  hold on
  plot(MPD_AH_grid, range_grid, 'ro')
  hold off
  xlim([0 10])
  ylim([0 6])
  % grid(gca,'minor')
  grid on
  set(gca, 'YMinorTick','on', 'YMinorGrid','on')
  xlabel('Absolute humidity (g m^{-3})'); 
  ylabel('Range (km)'); 
 
  Scrsize=[1 1 800 800];
  cd('/Users/spuler/Desktop/mpd/Plots/')
  FigH = figure(115);
  set(gca,'Fontsize',30,'Fontweight','b'); % 
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', Scrsize);
  name=strcat(sonde_date, 'Sonde_profile'); 
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution
  
  
end



%pause

