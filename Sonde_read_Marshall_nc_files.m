function[sonde_AH_grid, MPD_AH_grid, range_grid] = Sonde_read_Marshall_nc_files(jj, elevation, sondedir, sondefilename, N_avg_comb, duration, range_grid_size, range_grid_in,  comb_AH_var,  sonde_end_int, flag) 

%sondedir

filename = [sondedir sondefilename{jj}]; 
%filename = '/scr/sci/tammy/mpd/sgp/soundings/sgpsondewnpnC1.b1.20190429.023100.cdf';
%filename = '/Users/spuler/downloads/sgpsondewnpnC1.b1.20190501.083100.cdf';
%filename = '/Users/spuler/Desktop/mpd_03_processed_data/Sondes/Marshall_Field_Site_20201016_163106.nc';
%filename = '/Users/spuler/Desktop/mpd_03_processed_data/Sondes/Marshall_Field_Site_20201007_183911.nc';

sonde_date = filename(end-15:end-10);
sonde_time = filename(end-8:end-3);
n = datenum([sonde_date sonde_time], 'yymmddHHMMSS');
datestr(n)

ncid = netcdf.open(filename, 'NC_NOWRITE');
  %ncdisp(filename, '/', 'min') % use this to display all variables
  ncdisp(filename) % use this to display all variables
    
 variable{1} = 'time'; % time (seconds since 1970-1-1)
 %variable{2} = 'time_offset'; % time (seconds since 1970-1-1)
 variable{3} = 'alt'; % altitude above mean sea level (m)
 variable{4} = 'pres'; % pressure (hPa or mbar)
 variable{5} = 'tdry'; %Dry Bulb Temperature (C)
 variable{6} = 'rh'; % relative humidity (%)
 variable{7} = 'lat'; % latitude
 variable{8} = 'lon'; % longitude 
 variable{9} = 'reference_alt'; 
  
 base_time  = ncread(filename,variable{1});   
 %sonde_time = ncread(filename,variable{2}); 
 sonde_alt = ncread(filename,variable{3});
 sonde_P = ncread(filename,variable{4});  
 sonde_T = ncread(filename,variable{5});  
 sonde_RH = ncread(filename,variable{6}); 
 sonde_lat = ncread(filename,variable{7}); 
 sonde_lon = ncread(filename,variable{8}); 
 reference_alt = ncread(filename,variable{9});  
 
netcdf.close(ncid);

%convert sonde from Unix time to date number (days since Jan 0 0000) 
%duration_sonde = datenum(datetime(int64(sonde_time) + int64(base_time), 'convertfrom', 'posixtime'));
%sonde_offset_min = 30;
duration_sonde = n +  base_time/24/60/60;
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
range_grid = 0:range_grid_size/1000:6;
[sonde_AGL_km, index] = unique(sonde_AGL/1000); 
sonde_AH_grid =interp1(sonde_AGL_km, sonde_AH(index), range_grid, 'nearest');
% find the closes time index for the MPD water vapor
[minValue, closestIndex] = min(abs(min(duration_sonde)-duration))
[minValue, closestIndex_end] = min(abs(min(duration_sonde+sonde_end_int/24/60)-duration))
%MPD_AH = N_avg_comb(closestIndex,:).*1e6./6.022E23.*18.015;
%MPD_AH_var =  comb_AH_var(closestIndex,:);
MPD_AH = mean(N_avg_comb(closestIndex:closestIndex_end,:),1, 'omitnan').*1e6./6.022E23.*18.015;
if flag.data_type == 0
  MPD_AH_var =  median(comb_AH_var(closestIndex:closestIndex_end,:),1, 'omitnan')./sqrt(sonde_end_int/10); %assumes 10 min in the Matlab data
else
  MPD_AH_var =  median(comb_AH_var(closestIndex:closestIndex_end,:),1, 'omitnan');
end
% remove isolated points
test = ~isnan(MPD_AH);
test2 = movmean(test, 5);
test3 = (test2>0.2);
test4 = (test==1 & test3==1);

try
 MPD_AH_grid = interp1(range_grid_in(test4)/1000, MPD_AH(test4), range_grid, 'linear');
 MPD_AH_var_grid = interp1(range_grid_in(test4)/1000, MPD_AH_var(test4), range_grid, 'linear');
% MPD_AH_grid = interp1(range_grid_in(~isnan(MPD_AH))/1000, MPD_AH(~isnan(MPD_AH)), range_grid, 'linear');
%MPD_AH_var_grid = interp1(range_grid_in(~isnan(MPD_AH_var))/1000, MPD_AH_var(~isnan(MPD_AH_var)), range_grid, 'linear');
end

% don't plot data above a large gap in information (don't use clouds)
if flag.data_type == 0
  range_index = find(test2<=0.2, 1, 'first')
  MPD_AH_grid(range_index:end) = nan;
  MPD_AH_var_grid(range_index:end) = nan;
end
% figure(501)
% plot(MPD_AH, range/1000, 'bo')
% hold on
% plot(MPD_AH_grid, range_grid, 'r+')
% plot(sonde_AH_grid, range_grid)
% ylim([0 9])
% xlim([0 5])



if flag.plot_overlay == 1
  % overlay sonde vs MPD
  figure(115)

  plot(sonde_AH_grid, range_grid)
  hold on

  plot(MPD_AH_grid, range_grid, 'ro')
  
  eb(1) = errorbar(MPD_AH_grid, range_grid, MPD_AH_var_grid, 'horizontal', 'LineStyle', 'none', 'HandleVisibility','off');
  set(eb, 'color', 'r', 'LineWidth', 1)
  
  %shadedErrorBar(range_grid, MPD_AH_grid, MPD_AH_var_grid); 
  %set(gca,'YDir','reverse');
  %camroll(90)
  
  hold off
  xlim([0 8])
  ylim([0 6])
  % grid(gca,'minor')
  grid on
  set(gca, 'YMinorTick','on', 'YMinorGrid','on')
  xlabel('Absolute humidity (g m^{-3})'); 
  ylabel('Range (km)'); 
  
   
  Scrsize=[1 1 800 800];
  cd('/Users/spuler/Desktop/mpd/Plots/')
 % cd('/Users/lroot/Desktop/mpd/Plots/') 
  FigH = figure(115);
  set(gca,'Fontsize',30,'Fontweight','b'); % 
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', Scrsize);
  name=strcat(sonde_date, 'Sonde_profile'); 
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
  
 % save(name, 'range_grid', 'sonde_AH_grid', 'MPD_AH_grid', 'MPD_AH_var_grid')
  
  
  
end



%pause

