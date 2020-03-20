% example batch process 
start = '20200301';
stop = '20200319';

start_day = datenum(start,'yyyymmdd');
stop_day = datenum(stop,'yyyymmdd');
k=start_day;


for k=start_day:stop_day
  file = datestr(k, 'yyyymmdd');  
  % save_quicklook, save_data, save_netCDF, save_catalog, near, afterpulse, node, daystr  
  DIAL_process_NetCDF_function(0,1,0,0,0,0,'MPD4',file) % MPD# 4 (high range)
  DIAL_process_NetCDF_function(0,1,0,0,1,0,'MPD4',file) % MPD# 4 (low range)
  DIAL_process_NetCDF_function(0,1,0,0,0,0,'MPD3',file) % MPD #3 (high gain)
  DIAL_process_NetCDF_function(0,1,0,0,1,0,'MPD3',file) % MPD #3 (low gain)
  DIAL_process_NetCDF_function(0,1,0,0,0,1,'MPD3',file) % MPD #3 (high gain with afterpulse)
  DIAL_process_NetCDF_function(0,1,0,0,1,1,'MPD3',file) % MPD #3 (low gain with afterpulse)
end
DIAL_multiday_plots_NetCDF_function(1,0,0,0,'MPD4', start, stop) % MPD# 4 (high range)
DIAL_multiday_plots_NetCDF_function(1,0,1,0,'MPD4', start, stop) % MPD# 4 (low range)
DIAL_multiday_plots_NetCDF_function(1,0,0,0,'MPD3', start, stop) % MPD #3 (high gain)
DIAL_multiday_plots_NetCDF_function(1,0,1,0,'MPD3', start, stop) % MPD #3 (low gain)
DIAL_multiday_plots_NetCDF_function(1,0,0,1,'MPD3', start, stop) % MPD #3 (high gain with afterpulse)
DIAL_multiday_plots_NetCDF_function(1,0,1,1,'MPD3', start, stop) % MPD #3 (low gain with afterpulse)

 
%  for k=1:1
%    % save_quicklook, save_data, save_netCDF, save_catalog, near, afterpulse, node, daystr  
%    DIAL_process_NetCDF_function(0,1,0,0,0,0,'MPD4','20200319') % MPD# 4 (high range)
%    DIAL_process_NetCDF_function(0,1,0,0,1,0,'MPD4','20200319') % MPD# 4 (low range)
%    DIAL_process_NetCDF_function(0,1,0,0,0,0,'MPD3','20200319') % MPD #3 (high gain)
%    DIAL_process_NetCDF_function(0,1,0,0,1,0,'MPD3','20200319') % MPD #3 (low gain)
%    DIAL_process_NetCDF_function(0,1,0,0,0,1,'MPD3','20200319') % MPD #3 (high gain with afterpulse)
%    DIAL_process_NetCDF_function(0,1,0,0,1,1,'MPD3','20200319') % MPD #3 (low gain with afterpulse)
%    % save_figs, save_data, near, afterpulse, node, daystr, daystr2
%    DIAL_multiday_plots_NetCDF_function(1,0,0,0,'MPD4','20200310','20200319') % MPD# 4 (high range)
%    DIAL_multiday_plots_NetCDF_function(1,0,1,0,'MPD4','20200310','20200319') % MPD# 4 (low range)
%    DIAL_multiday_plots_NetCDF_function(1,0,0,0,'MPD3','20200310','20200319') % MPD #3 (high gain)
%    DIAL_multiday_plots_NetCDF_function(1,0,1,0,'MPD3','20200310','20200319') % MPD #3 (low gain)
%    DIAL_multiday_plots_NetCDF_function(1,0,0,1,'MPD3','20200310','20200319') % MPD #3 (high gain with afterpulse)
%    DIAL_multiday_plots_NetCDF_function(1,0,1,1,'MPD3','20200310','20200319') % MPD #3 (low gain with afterpulse)
%  end