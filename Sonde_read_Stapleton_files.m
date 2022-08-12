% clear all
% close all
%  download data from https://ruc.noaa.gov/raobs/

cd('/Users/spuler/Desktop/mpd/sondes/')
%folder = 'raob_soundings49530.cdf';
%folder = 'raob_soundings59417.cdf';
folder = 'raob_soundings18455.cdf';
folder = 'raob_soundings33531.cdf';
folder = 'raob_soundings6331.cdf';
folder = 'raob_soundings8014.cdf'; %GJT 72476 1-Jun to 10-Aug 2022
folder = 'raob_soundings26234.cdf'; %DNR 72469 1-Jun to 10-Aug 2022

ncid = netcdf.open(folder, 'NC_NOWRITE');
  ncdisp(folder) % use this to display all variables

   R_Man = ncread(folder,'htMan'); %'Geopotential - Mandatory level'
   P_Man = ncread(folder,'prMan'); %'Pressure - Mandatory level'
   T_Man = ncread(folder,'tpMan'); %'Temperature (K) - Mandatory level'
   H_Man = ncread(folder,'tdMan'); %'Dew Point Depression - Mandatory level'
   
   R_sigT = ncread(folder,'htSigT'); %'Geopotential - Significant level wrt T'
   P_sigT = ncread(folder,'prSigT'); %'Pressure - Significant level wrt T'
   T_sigT = ncread(folder,'tpSigT'); %'Temperature (K) - Significant level wrt T' 
   H_sigT = ncread(folder,'tdSigT'); %'Dew Point Depression - Significant level wrt T'
  
   synTime = ncread(folder,'synTime');  %Synoptic time
   time = ncread(folder,'relTime');  %sounding release time
   Elev = ncread(folder,'staElev');  %'Station Elevation'
   %Type = ncread(folder,'sondTyp');  %'Instrument Type'

 netcdf.close(ncid);
 
 t = datetime(time,'ConvertFrom','posixtime');
 t_synoptic = datetime(synTime,'ConvertFrom','posixtime');
 
 R_sigT(T_sigT > 9e36) = nan;
 T_sigT(T_sigT > 9e36) = nan;
 P_sigT(P_sigT > 9e36) = nan;
 
 R_Man(T_Man == 99999) = nan;
 R_Man(R_Man > 9e36) = nan;
 T_Man(T_Man == 99999) = nan;
 P_Man(P_Man == 99999) = nan;

% Offset (m) for testing
Offset = 750;

 
 
for i=1:size(t,1) 
  sonde_AGL_sigT = R_sigT(:,i)-Elev(i)+Offset;
  sonde_T_sigT = T_sigT(:,i);
  %sonde_duration_sigT = datenum(t(i))*ones(size(sonde_AGL_sigT));
    
  sonde_AGL_Man = R_Man(:,i)-Elev(i)+Offset;
  sonde_T_Man = T_Man(:,i);
  %sonde_duration_Man = datenum(t(i))*ones(size(sonde_AGL_Man));
  
   % grid sonde data vs range 
   range_grid = 0:37.5/1000:10;
   sonde_AGL =  vertcat(sonde_AGL_Man, sonde_AGL_sigT);
   sonde_T = vertcat(sonde_T_Man, sonde_T_sigT);
   [sonde_AGL_sort, ind] = sort(sonde_AGL);
   sonde_T_sort = sonde_T(ind);
   sonde_AGL_check = sonde_AGL_sort(~isnan(sonde_AGL_sort));
   [sonde_AGL_km, index] = unique(sonde_AGL_check/1000); 
   sonde_T_grid =interp1(sonde_AGL_km, sonde_T_sort(index), range_grid, 'linear');
   sonde_duration = datenum(t(i))*ones(size(range_grid));
   
   figure(1)
   hold on
   %scatter(sonde_duration_Man, sonde_AGL_Man/1000, 100, sonde_T_Man, 's', 'filled');
   %scatter(sonde_duration_sigT, sonde_AGL_sigT/1000, 100, sonde_T_sigT, 's', 'filled');
   scatter(sonde_duration, range_grid, 15, sonde_T_grid, 'o', 'filled');
   colormap(jet)
   ylim([0 6])
   %xlim([datenum(t(1)) datenum(t(end))])
   %caxis([240 320])
   colorbar
   
  figure(i+3) 
  hold on
  plot(sonde_T_grid, range_grid*1000,'g*')
  hold off
      
end

for i=1:(size(t,1)) 
  figure(i+3) 
  hold on
  plot(T_sigT(:,i),R_sigT(:,i)-Elev(i)+Offset, 'b+')
  plot(T_Man(:,i),R_Man(:,i)-Elev(i)+Offset, 'ro')
  hold off
  title(datestr(t(i)))
  ylim([0 10000])
end



cd('/Users/spuler/Desktop') % point to the directory where data is stor
FigH = figure(1);
set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920 250]);
name=strcat('Stapleton_sonde_overlay');
print(FigH, name, '-dpng', '-r0')  


  