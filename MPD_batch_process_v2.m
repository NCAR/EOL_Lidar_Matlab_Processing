% example batch process

start_date = '20250214';   
stop_date =  '20250214'; 

if strcmp(getenv('HOSTNAME'),'eol-smaug.eol.ucar.edu') == 1  % when running on fog server
  prompt = 'Enter start date (e.g., 20230411)';
  start_date = input(prompt, 's');
  prompt = 'Enter stop date (e.g., 20230411)';
  stop_date = input(prompt, 's');
end

start_day = datenum(start_date,'yyyymmdd')
stop_day = datenum(stop_date,'yyyymmdd') 
k=start_day; 
flag.process = 1;
flag.plot = 0;

for n=1:1
tStart = tic   
    if flag.process == 1
      for k=start_day:stop_day
        file = datestr(k, 'yyyymmdd');  
        % save_quicklook, save_data, save_netCDF, save_catalog, channels, correction, node, daystr 
%          MPD_process_NetCDF_function_v2(0,1,0,0,'WV','AP_OFF','MPD01',file)
%          MPD_process_NetCDF_function_v2(0,1,0,0,'WV','AP_OFF','MPD02',file) 
%          MPD_process_NetCDF_function_v2(0,1,0,0,'WV','AP_ON','MPD03',file) 
%          MPD_process_NetCDF_function_v2(0,1,0,0,'WV','AP_ON','MPD04',file)  
%          MPD_process_NetCDF_function_v2(0,1,0,0,'WV','AP_OFF','MPD05',file)    
%          MPD_process_NetCDF_function_v2(0,1,0,0,'O2','AP_ON','MPD01',file)  
%          MPD_process_NetCDF_function_v2(0,1,0,0,'O2','AP_OFF','MPD02',file) 
%          MPD_process_NetCDF_function_v2(0,1,0,0,'O2','AP_OFF','MPD03',file)
%          MPD_process_NetCDF_function_v2(0,1,0,0,'O2','AP_OFF','MPD04',file) 
          MPD_process_NetCDF_function_v2(0,1,0,0,'O2','AP_OFF','MPD05',file) 
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