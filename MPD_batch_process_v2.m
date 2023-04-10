% example batch process

start_date = '20230408';   
stop_date =  '20230410'; 

% start_date = '20220528'; %PRECIP   
% stop_date =  '20220810'; 

start_day = datenum(start_date,'yyyymmdd');
stop_day = datenum(stop_date,'yyyymmdd');
k=start_day;
flag.process = 1;
flag.plot = 0;

for n=1:1
tStart = tic   
    if flag.process == 1
      for k=start_day:stop_day
        file = datestr(k, 'yyyymmdd');  
        % save_quicklook, save_data, save_netCDF, save_catalog, channels, correction, node, daystr 
% %          MPD_process_NetCDF_function_v2(0,1,0,0,'WV','AP_OFF','MPD01',file)
%          MPD_process_NetCDF_function_v2(0,1,0,0,'WV','AP_ON','MPD02',file) 
%          MPD_process_NetCDF_function_v2(0,1,0,0,'WV','AP_ON','MPD03',file) 
%          MPD_process_NetCDF_function_v2(0,1,0,0,'WV','AP_ON','MPD04',file)  
% %          MPD_process_NetCDF_function_v2(0,1,0,0,'WV','AP_OFF','MPD05',file)    
%   MPD_process_NetCDF_function_v2(0,1,0,0,'O2','AP_ON','MPD01',file) %Smart switch testing started Feb 2023 
    MPD_process_NetCDF_function_v2(0,1,0,0,'O2','AP_ON','MPD03',file) % Temp began late Feb 2023           %
    MPD_process_NetCDF_function_v2(0,1,0,0,'O2','AP_ON','MPD05',file) 
      end
    end

  start_date = '20200929';
  stop_date =  '20201005';

    if flag.plot == 1
      % save_figs, save_data, near/low, afterpulse, node, daystr, daystr2, skip(tick marks at every day 1, other day, 2 etc)
     % MPD_multiday_plots_NetCDF_function(1,1,0,0,'MPD01', start_date, stop_date, 1) 
     % MPD_multiday_plots_NetCDF_function(1,1,0,0,'MPD02', start_date, stop_date, 1) 
     % MPD_multiday_plots_NetCDF_function(1,1,0,0,'MPD03', start_date, stop_date, 1) 
     % MPD_multiday_plots_NetCDF_function(1,1,0,0,'MPD04', start_date, stop_date, 1) 
     % MPD_multiday_plots_NetCDF_function(1,1,0,0,'MPD05', start_date, stop_date, 1) 
    end
tElapsed = toc(tStart) 
end