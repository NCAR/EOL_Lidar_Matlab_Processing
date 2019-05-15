function[] = Sonde_read_nc_files(jj, elevation, sondedir, sondefilename) 

%sondedir

filename = [sondedir sondefilename{jj}] 
%filename = '/volumes/documents/WV_DIAL_data/SGP_sondes/sgpsondewnpnC1.b1.20190429.023100.cdf';
%filename = '/Users/spuler/downloads/sgpsondewnpnC1.b1.20190501.083100.cdf';
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
 
 base_time  = ncread(filename,variable{1});   

 sonde_time = ncread(filename,variable{2});   
 sonde_alt = ncread(filename,variable{3});
 sonde_P = ncread(filename,variable{4});  
 sonde_T = ncread(filename,variable{5});  
 sonde_RH = ncread(filename,variable{6});  
 
netcdf.close(ncid);

%convert sonde from Unix time to date number (days since Jan 0 0000) 
duration_sonde = datenum(datetime(int64(sonde_time) + int64(base_time), 'convertfrom', 'posixtime'));
sonde_AGL = sonde_alt - elevation;

%figure(1)
%plot(T_sonde,alt)
%figure(2)
%plot(P_sonde,alt)
%figure(3)
%plot(RH,alt)
    
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
 
figure(10)
plot(sonde_AH, sonde_AGL/1000)
xlim([0 20])
ylim([0 6])

figure(11)
scatter(duration_sonde, sonde_AGL/1000, 25, sonde_AH, '+');
colormap(jet)
ylim([0 6])
%clim([0 20])
colorbar

figure(1)
hold on
scatter(duration_sonde, sonde_AGL/1000, 10, sonde_AH, '+');
colormap(jet)

