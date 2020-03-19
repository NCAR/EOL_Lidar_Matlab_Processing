% example batch process 10-Mar-2020 to 18-Mar-2020
for k=9:10
  % save_quicklook, save_data, save_netCDF, save_catalog, near, afterpulse, node, daystr  
  DIAL_process_NetCDF_function(0,1,0,0,0,0,'MPD4',strcat('2020031',num2str(k-1))) % MPD# 4 (high range)
  DIAL_process_NetCDF_function(0,1,0,0,1,0,'MPD4',strcat('2020031',num2str(k-1))) % MPD# 4 (low range)
  DIAL_process_NetCDF_function(0,1,0,0,0,0,'MPD3',strcat('2020031',num2str(k-1))) % MPD #3 (high gain)
  DIAL_process_NetCDF_function(0,1,0,0,1,0,'MPD3',strcat('2020031',num2str(k-1))) % MPD #3 (low gain)
  DIAL_process_NetCDF_function(0,1,0,0,0,1,'MPD3',strcat('2020031',num2str(k-1))) % MPD #3 (high gain with afterpulse)
  DIAL_process_NetCDF_function(0,1,0,0,1,1,'MPD3',strcat('2020031',num2str(k-1))) % MPD #3 (low gain with afterpulse)
end
 % save_figs, save_data, near, afterpulse, node, daystr, day_num
 DIAL_multiday_plots_NetCDF_function(1,0,0,0,'MPD4',20200310,10) % MPD# 4 (high range)
 DIAL_multiday_plots_NetCDF_function(1,0,1,0,'MPD4',20200310,10) % MPD# 4 (low range)
 DIAL_multiday_plots_NetCDF_function(1,0,0,0,'MPD3',20200310,10) % MPD #3 (high gain)
 DIAL_multiday_plots_NetCDF_function(1,0,1,0,'MPD3',20200310,10) % MPD #3 (low gain)
 DIAL_multiday_plots_NetCDF_function(1,0,0,1,'MPD3',20200310,10) % MPD #3 (high gain with afterpulse)
 DIAL_multiday_plots_NetCDF_function(1,0,1,1,'MPD3',20200310,10) % MPD #3 (low gain with afterpulse)