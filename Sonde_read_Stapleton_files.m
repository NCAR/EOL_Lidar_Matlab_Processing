% clear all
% close all
%  download data from https://ruc.noaa.gov/raobs/

cd('/Users/spuler/Desktop/mpd/sondes/')
%folder = 'raob_soundings49530.cdf';
%folder = 'raob_soundings59417.cdf';
folder = 'raob_soundings18455.cdf';
ncid = netcdf.open(folder, 'NC_NOWRITE');
 % ncdisp(folder) % use this to display all variables

   R_Man = ncread(folder,'htMan'); %'Geopotential - Mandatory level'
   P_Man = ncread(folder,'prMan'); %'Pressure - Mandatory level'
   T_Man = ncread(folder,'tpMan'); %'Temperature (K) - Mandatory level'
   H_Man = ncread(folder,'tdMan'); %'Dew Point Depression - Mandatory level'
   
   R_sigT = ncread(folder,'htSigT'); %'Geopotential - Significant level wrt T'
   P_sigT = ncread(folder,'prSigT'); %'Pressure - Significant level wrt T'
   T_sigT = ncread(folder,'tpSigT'); %'Temperature (K) - Significant level wrt T' 
   H_sigT = ncread(folder,'tdSigT'); %'Dew Point Depression - Significant level wrt T'
  
   time = ncread(folder,'relTime');  %sounding release time
   Elev = ncread(folder,'staElev');  %'Station Elevation'
   %Type = ncread(folder,'sondTyp');  %'Instrument Type'

 netcdf.close(ncid);
 
 t = datetime(time,'ConvertFrom','posixtime');
 
 R_sigT(T_sigT > 9e36) = nan;
 T_sigT(T_sigT > 9e36) = nan;
 P_sigT(P_sigT > 9e36) = nan;
 
 R_Man(T_Man == 99999) = nan;
 R_Man(R_Man > 9e36) = nan;
 T_Man(T_Man == 99999) = nan;
 P_Man(P_Man == 99999) = nan;

 i=1;
for i=1:size(t,1) 
 % figure(i) 
%   plot(T_sigT(:,i),R_sigT(:,i)-Elev(i))
%   hold on
%   plot(T_Man(:,i),R_Man(:,i)-Elev(i))
%   hold off
%   title(datestr(t(i)))
%   ylim([0 10000])
  
  sonde_AGL = R_sigT(:,i)-Elev(i);
  sonde_T = T_sigT(:,i);
  sonde_duration = datenum(t(i))*ones(size(sonde_AGL));
 
  figure(1)
  hold on
  scatter(sonde_duration, sonde_AGL/1000, 50, sonde_T, 'o', 'filled');
  colormap(jet)
  %ylim([0 6])
  %xlim([datenum(t(1)) datenum(t(end))])
  caxis([240 320])
  colorbar
    
end
 

