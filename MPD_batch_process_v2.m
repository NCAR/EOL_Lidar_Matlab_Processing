% example batch process
start_date = '20200409';
stop_date = '20200413';

start_day = datenum(start_date,'yyyymmdd');
stop_day = datenum(stop_date,'yyyymmdd');
k=start_day;
flag.process = 1;
flag.plot = 1;

for n=1:1
tStart = tic   
    if flag.process == 1
      for k=start_day:stop_day
        file = datestr(k, 'yyyymmdd');  
        % save_quicklook, save_data, save_netCDF, save_catalog, channels, correction, node, daystr 
        MPD_process_NetCDF_function_v2(0,1,0,0,'ALL','AP_ON','MPD04',file) % MPD #4 
        MPD_process_NetCDF_function_v2(0,1,0,0,'ALL','AP_ON','MPD03',file) % MPD #3 
      end
    end

    start_date = '20200331';
    stop_date = '20200413';

    if flag.plot == 1
      % save_figs, save_data, near, afterpulse, node, daystr, daystr2
      MPD_multiday_plots_NetCDF_function(1,0,0,0,'MPD04', start_date, stop_date) % MPD# 4 (high range) 
      MPD_multiday_plots_NetCDF_function(1,0,1,0,'MPD04', start_date, stop_date) % MPD# 4 (low range)
      MPD_multiday_plots_NetCDF_function(1,0,0,1,'MPD04', start_date, stop_date) % MPD# 4 (high range with afterpulse) 
      MPD_multiday_plots_NetCDF_function(1,0,1,1,'MPD04', start_date, stop_date) % MPD# 4 (low range with afterpulse)  
      MPD_multiday_plots_NetCDF_function(1,0,0,0,'MPD03', start_date, stop_date) % MPD%  #3 (high gain)
      MPD_multiday_plots_NetCDF_function(1,0,1,0,'MPD03', start_date, stop_date) % MPD #3 (low gain)
      MPD_multiday_plots_NetCDF_function(1,0,0,1,'MPD03', start_date, stop_date) % MPD #3 (high gain with afterpulse)
      MPD_multiday_plots_NetCDF_function(1,0,1,1,'MPD03', start_date, stop_date) % MPD #3 (low gain with afterpulse)
    end
tElapsed = toc(tStart) 
end