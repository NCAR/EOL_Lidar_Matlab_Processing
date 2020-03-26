function[] = MPD_process_NetCDF_function(save_quicklook, save_data, save_netCDF, save_catalog, near, afterpulse, node, daystr)
%clear all; 
%close all
%start_date = '20200325';
%save_quicklook=0; save_data=1; save_netCDF=0; save_catalog=0; near= 1; afterpulse=1; node='MPD3'; daystr=start_date; 

flag.save_quicklook = save_quicklook;  % save quicklook to local directory
flag.save_data = save_data;  % save files in matlab format
flag.save_netCDF = save_netCDF; % save files netCDF format
flag.save_catalog = save_catalog; % upload quicklook (and data) to field catalog

flag.mask_data = 1;  % mask applied to data based on error analysis threshold
flag.gradient_filter = 1;  % this is used to mask regions with 'high' backscatter gradients which tend to cause errors
flag.pileup = 1; % use pileup correction for detectors
flag.WS = 1; % use the surface weather station data to calcuate spectroscopy
flag.decimate = 0; % decimate all data to half the wv resoltuion
flag.int = 0; % interpolate nans in nanmoving_average
flag.mark_gaps = 1; % sets gaps in data to NaNs
flag.OF = 1; % correct for geometric overlap functions
flag.near = near; %process the near range channel (or low gain)
flag.afterpulse = afterpulse; % correct for afterpulsing (in progress on MPD#3 only)

flag.plot_data = 1;  % need to have this one to save the figs
flag.troubleshoot = 0; % shows extra plots used for troubleshooting
p_hour = 12; % hour to show troubleshooting profiles

ave_time.wv = 10.0; % averaging time (in minutes) for the water vapor 
ave_time.rb = 1.0; % averaging time (in minutes) for the relative backscatter
ave_time.gr = 0.5; % gridding time (in minutes) for the output files (HK data at 2 sec)

if strcmp(getenv('HOSTNAME'),'fog.eol.ucar.edu')
   serv_path = '/export/fog1/rsfdata/MPD/'; % when running on server
elseif strcmp(getenv('HOSTNAME'),'')
    serv_path = '/Users/spuler/Desktop/'; % when running on server   
else 
   serv_path = '/Volumes/eol/fog1/rsfdata/MPD/'; % 
end

if strcmp(node,'MPD1')==1
  write_data_folder = strcat(serv_path, 'wvdial_1_processed_data/Matlab');
  files = strcat(serv_path, 'wvdial_1_data/2020/', daystr);
  catalog = '/pub/incoming/catalog/operations';
elseif strcmp(node,'MPD2')==1
  write_data_folder = strcat(serv_path, 'wvdial_2_processed_data/Matlab');
  files = strcat(serv_path, 'wvdial_2_data/2020/', daystr);
 catalog = '/pub/incoming/catalog/operations';
elseif strcmp(node,'MPD3')==1
  write_data_folder = strcat(serv_path, 'wvdial_3_processed_data/Matlab');
  files = strcat(serv_path, 'wvdial_3_data/2020/', daystr);
  catalog = '/pub/incoming/catalog/operations';
  if flag.afterpulse == 1
    write_data_folder = strcat(serv_path, 'wvdial_3_processed_data/Matlab/afterpulse');
  end
elseif strcmp(node,'MPD4')==1
  write_data_folder = strcat(serv_path, 'wvdial_4_processed_data/Matlab');
  files = strcat(serv_path, 'wvdial_4_data/2020/', daystr);
 catalog = '/pub/incoming/catalog/operations';
elseif strcmp(node,'MPD5')==1
  write_data_folder = strcat(serv_path, 'wvdial_5_processed_data/Matlab');
  files = strcat(serv_path, 'wvdial_5_data/2020/', daystr);
 catalog = '/pub/incoming/catalog/operations';
end


folder = files;
date = textscan(folder(end-5:end), '%6f'); date=date{1};  % read date of file
MPD_read_calvals
folder_in=folder;
date_in = date;
MPD_Analysis_function_NetCDF_v4(folder, date, MCS, write_data_folder, flag, node, wavemeter_offset,...
    profiles2ave, P0, switch_ratio, ave_time, timing_range_correction, blank_range, p_hour, catalog, Afterpulse_File)%
   

end