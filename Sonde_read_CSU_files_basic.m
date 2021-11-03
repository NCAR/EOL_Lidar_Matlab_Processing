function[sonde_MR_grid, sonde_T_grid, sonde_P_grid, range_grid, sonde_date, sonde_time, duration_sonde] = Sonde_read_CSU_files_basic(jj, elevation, sondedir, sondefilename, step) 


filename = [sondedir sondefilename{jj}]; 
sonde_date = filename(end-16:end-9);
sonde_time = filename(end-7:end-4);
n = datenum([sonde_date sonde_time], 'yyyymmddHHMM');
datestr(n)

%% Initialize variables.
startRow = 28;
formatSpec = '%12s%10s%6s%6s%6s%6s%6s%6s%8s%10s%12s%s%[^\n\r]';
fileID = fopen(filename,'r');
textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false, 'EndOfLine', '\r\n');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'ReturnOnError', false);
fclose(fileID);
 

 sonde_t = str2double(dataArray{1,1});  % elapsed time (s)
 sonde_alt = str2double(dataArray{1,2});  % height above MSL (m)
 sonde_alt = str2double(dataArray{1,12});
 sonde_P = str2double(dataArray{1,3}); % P (hPa) (1 hPa = 1 mbar)
 sonde_T = str2double(dataArray{1,4}); % T (C)
 sonde_RH = str2double(dataArray{1,6}); % RH
 sonde_lat = str2double(dataArray{1,10});
 sonde_lon = str2double(dataArray{1,11});

 sonde_t = sonde_t(3:end,1);
 sonde_alt = sonde_alt(3:end,1);
 sonde_P = sonde_P(3:end,1);
 sonde_T = sonde_T(3:end,1);
 sonde_RH = sonde_RH(3:end,1);
 sonde_lat = sonde_lat(3:end,1);
 sonde_lon = sonde_lon(3:end,1);
%convert sonde from Unix time to date number (days since Jan 0 0000) 
%duration_sonde = datenum(datetime(int64(sonde_time) + int64(base_time), 'convertfrom', 'posixtime'));
%sonde_offset_min = 30;

duration_sonde = n + sonde_t/24/60/60;
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

%sonde_N_air = (sonde_P.*N_A)./((sonde_T+273.15).*R)./1e4; % number density of air in cm^3
%sonde_N_air = (sonde_P*9.869233E-4)./((sonde_T+273.15).*1.380649e-23)*101.325;
%sonde_vmr = sonde_N_H2O./sonde_N_air; %water vapor volume mixing ratio


% grid sonde data at 1km range increments to 8km
range_grid = 0:step:8;
[sonde_AGL_km, index] = unique(sonde_AGL/1000); 
sonde_AH_grid =interp1(sonde_AGL_km, sonde_AH(index), range_grid, 'linear');
sonde_MR_grid =interp1(sonde_AGL_km, sonde_MR(index), range_grid, 'linear');
sonde_T_grid =interp1(sonde_AGL_km, sonde_T(index), range_grid, 'linear');
sonde_P_grid =interp1(sonde_AGL_km, sonde_P(index), range_grid, 'linear');
  
end



% pause

