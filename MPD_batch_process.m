% example batch process
start_date = '20200331';
stop_date = '20200401';

start_day = datenum(start_date,'yyyymmdd');
stop_day = datenum(stop_date,'yyyymmdd');
k=start_day;
flag.process = 1;
flag.plot = 1;

for n=1:1
    if flag.process == 1
      for k=start_day:stop_day
        file = datestr(k, 'yyyymmdd');  
        % save_quicklook, save_data, save_netCDF, save_catalog, near, afterpulse, node, daystr
        MPD_process_NetCDF_function(0,1,0,0,0,0,'MPD4',file) % MPD #4 (high range) 
        MPD_process_NetCDF_function(0,1,0,0,1,0,'MPD4',file) % MPD #4(low range) 
        MPD_process_NetCDF_function(0,1,0,0,0,1,'MPD4',file) % MPD #4 (high range with afterpulse)
        MPD_process_NetCDF_function(0,1,0,0,1,1,'MPD4',file) % MPD #4 (low range with afterpulse)   
        MPD_process_NetCDF_function(0,1,0,0,0,0,'MPD3',file) % MPD #3 (high gain) 
        MPD_process_NetCDF_function(0,1,0,0,1,0,'MPD3',file) % MPD #3 (low gain)
        MPD_process_NetCDF_function(0,1,0,0,0,1,'MPD3',file) % MPD #3 (high gain with afterpulse)
        MPD_process_NetCDF_function(0,1,0,0,1,1,'MPD3',file) % MPD #3 (low gain with afterpulse)
      end
    end

    start_date = '20200319';
    stop_date = '20200401';

    if flag.plot == 1
      % save_figs, save_data, near, afterpulse, node, daystr, daystr2
      MPD_multiday_plots_NetCDF_function(1,0,0,0,'MPD4', start_date, stop_date) % MPD# 4 (high range) 
      MPD_multiday_plots_NetCDF_function(1,0,1,0,'MPD4', start_date, stop_date) % MPD# 4 (low range)
      MPD_multiday_plots_NetCDF_function(1,0,0,1,'MPD4', start_date, stop_date) % MPD# 4 (high range with afterpulse) 
      MPD_multiday_plots_NetCDF_function(1,0,1,1,'MPD4', start_date, stop_date) % MPD# 4 (low range with afterpulse)  
      MPD_multiday_plots_NetCDF_function(1,0,0,0,'MPD3', start_date, stop_date) % MPD%  #3 (high gain)
      MPD_multiday_plots_NetCDF_function(1,0,1,0,'MPD3', start_date, stop_date) % MPD #3 (low gain)
      MPD_multiday_plots_NetCDF_function(1,0,0,1,'MPD3', start_date, stop_date) % MPD #3 (high gain with afterpulse)
      MPD_multiday_plots_NetCDF_function(1,0,1,1,'MPD3', start_date, stop_date) % MPD #3 (low gain with afterpulse)
    end
end