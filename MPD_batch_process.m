% example batch process 
start = '20200317';
stop = '20200323';

start_day = datenum(start,'yyyymmdd');
stop_day = datenum(stop,'yyyymmdd');
k=start_day;
flag.process = 1;
flag.plot = 1;

if flag.process == 1
  for k=start_day:stop_day
    file = datestr(k, 'yyyymmdd');  
    % save_quicklook, save_data, save_netCDF, save_catalog, near, afterpulse, node, daystr  
   % MPD_process_NetCDF_function(0,1,0,0,0,0,'MPD4',file) % MPD# 4 (high range)
   % MPD_process_NetCDF_function(0,1,0,0,1,0,'MPD4',file) % MPD# 4 (low range)
   % MPD_process_NetCDF_function(0,1,0,0,0,0,'MPD3',file) % MPD #3 (high gain)
   % MPD_process_NetCDF_function(0,1,0,0,1,0,'MPD3',file) % MPD #3 (low gain)
    MPD_process_NetCDF_function(0,1,0,0,0,1,'MPD3',file) % MPD #3 (high gain with afterpulse)
    MPD_process_NetCDF_function(0,1,0,0,1,1,'MPD3',file) % MPD #3 (low gain with afterpulse)
  end
end

start = '20200313';
stop = '20200323';

if flag.plot == 1
  %  MPD_multiday_plots_NetCDF_function(1,0,0,0,'MPD4', start, stop) % MPD# 4 (high range)
  %  MPD_multiday_plots_NetCDF_function(1,0,1,0,'MPD4', start, stop) % MPD# 4 (low range)
  %  MPD_multiday_plots_NetCDF_function(1,0,0,0,'MPD3', start, stop) % MPD #3 (high gain)
  %  MPD_multiday_plots_NetCDF_function(1,0,1,0,'MPD3', start, stop) % MPD #3 (low gain)
    MPD_multiday_plots_NetCDF_function(1,0,0,1,'MPD3', start, stop) % MPD #3 (high gain with afterpulse)
    MPD_multiday_plots_NetCDF_function(1,0,1,1,'MPD3', start, stop) % MPD #3 (low gain with afterpulse)
end