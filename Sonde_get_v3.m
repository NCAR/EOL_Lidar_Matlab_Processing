function[N_H2O, sonde_top, sonde_range, t, date, T_sonde, P_sonde, sonde_stop, Wind_speed] = Sonde_get_v3(sondedir, sondefilename, average) 

%d = pwd; % get the current path
%cd('/scr/rsf1/vanandel/HSRL/h2o_data/2012/sondes') % point to the directory where data is stored  

  sonde=strcat(sondedir, sondefilename);
      % read the date information from the sonde file name title
      sonde_year = strread(sondefilename(2:end), '%4c', 1);
      sonde_month = strread(sondefilename(6:end), '%2c', 1);
      sonde_day = strread(sondefilename(8:end), '%2c', 1);
      sonde_hour = strread(sondefilename(11:end), '%2c', 1);
      sonde_min = strread(sondefilename(13:end), '%2c', 1);
      sonde_sec = strread(sondefilename(15:end), '%2c', 1);
      header = importdata(sonde, ' ', 12)
      data = importdata(sonde, ' ', 14);
      P_sonde=data.data(:,5);  % pressure in mbar
      T_sonde=data.data(:,6);  % temperature in C
      RH = data.data(:,8); % relative humidity 
      Wind_speed=data.data(:,11); % wind speed (m/s)
      alt1=data.data(:,14); % altitude in meters
      alt=data.data(:,17); % gps altitude in meters

      
    % create a time stamp   
    t=datenum(1,1,1);
    t=addtodate(t, str2num(sonde_year)-1, 'year'); %
    t=addtodate(t, str2num(sonde_month)-1, 'month'); %
    t=addtodate(t, str2num(sonde_day)-1, 'day');
    t=addtodate(t, str2num(sonde_hour), 'hour'); % UTC
    t=addtodate(t, str2num(sonde_min), 'min');
    t=addtodate(t, str2num(sonde_sec), 'sec');
    
    % remove bad data points 
    %alt(alt==99999) = NaN;
    alt(alt==-999 | P_sonde==-999 | RH==-999 | T_sonde==-999 | alt<alt(1)) = NaN;
    RH = RH(~isnan(alt));
    T_sonde = T_sonde(~isnan(alt));
    P_sonde = P_sonde(~isnan(alt));
    Wind_speed = Wind_speed(~isnan(alt));
   % alt1 = alt1(~isnan(alt));
    alt = alt(~isnan(alt));
    

sonde_start=t;
sonde_stop=addtodate(t, average, 'min');  %25 min to reach 5 km

date=datestr(t, 'dd mmm yyyy, HH:MM')  % write sonde launch time to the screeen
 
%cd(d) % point back to original directory

%% convert RH to number density
% vapor pressure of water
    a0 = 6.107799961;
    a1 = 4.436518521E-1; 
    a2 = 1.428945805E-2;
    a3 = 2.650648471E-4; 
    a4 = 3.031240396E-6;
    a5 = 2.034080948E-8;
    a6 = 6.136820929E-11;
e=((a0+T_sonde.*(a1+T_sonde.*(a2+T_sonde.*(a3+T_sonde.*(a4+T_sonde.*(a5+T_sonde.*a6))))))./1); %vapor pressure in hPa 

% constants
R = 8.31447215; %J mol^-1 K^-1
N_A= 6.0221415E23; %mol^-1

RH_surf=1;
T_surf=1;

% convert from RH to number density
N_H2O = 1.*(RH.*(RH_surf).*e./(R.*(T_sonde+273).*(T_surf))).*N_A*1e-6;  %cm^3
 
% % get the spatial information
  ground = alt(1) %altidude of the sonde launch in meters
  sonde_range=(alt-ground)/1000; % convert sonde alt to range in km 
  sonde_range_stop = 12; % set top of sonde at 12 km 
  sonde_range = sonde_range(~isnan(N_H2O));
  T_sonde = T_sonde(~isnan(N_H2O));
  P_sonde = P_sonde(~isnan(N_H2O));
  Wind_speed = Wind_speed(~isnan(N_H2O));
  N_H2O = N_H2O(~isnan(N_H2O));

% 
 for i=1:length(sonde_range)
  if (sonde_range_stop<=sonde_range(i)) == 1 
    sonde_top = i
    break
  else
    sonde_top = i;
  end
 end


end

