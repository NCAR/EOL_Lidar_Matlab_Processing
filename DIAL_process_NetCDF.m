clear all; close all

node = 'MPD3';

%location to write files
    write_data_folder = uipickfiles('num',1,'out', 'char', 'prompt', ...
    'select folder to store data',  'FilterSpec', '/Users/spuler/Desktop/');
%    'select folder to store data',  'FilterSpec', '/Volumes/eol/fog1/rsfdata/MPD/');
%    'select folder to store data',  'FilterSpec', '/Volumes/documents/WV_DIAL_data/');

flag.save_quicklook = 0;  % save quicklook to local directory
flag.save_data = 1;  % save files in matlab format
flag.save_netCDF = 0; % save files netCDF format
flag.save_catalog = 0; % upload quicklook (and data) to field catalog

flag.mask_data = 1;  % mask applied to data based on error analysis threshold
flag.gradient_filter = 1;  % this is used to mask regions with 'high' backscatter gradients which tend to cause errors
flag.pileup = 1; % use pileup correction for detectors
flag.WS = 1; % use the surface weather station data to calcuate spectroscopy
flag.decimate = 0; % decimate all data to half the wv resoltuion
flag.int = 0; % interpolate nans in nanmoving_average
flag.mark_gaps = 1; % sets gaps in data to NaNs
flag.OF = 1; % correct for geometric overlap functions
flag.near = 0; %process the near range channel (or low gain)
flag.afterpulse = 1; % correct for afterpulsing (in progress on MPD#3 only)

flag.plot_data = 1;  % need to have this one to save the figs
flag.troubleshoot = 0; % shows extra plots used for troubleshooting
p_hour = 12; % hour to show troubleshooting profiles

ave_time.wv = 10.0; % averaging time (in minutes) for the water vapor 
ave_time.rb = 1.0; % averaging time (in minutes) for the relative backscatter
ave_time.gr = 0.5; % gridding time (in minutes) for the output files (HK data at 2 sec)

 if strcmp(getenv('HOSTNAME'),'fog.eol.ucar.edu')
     serv_path = '/export/fog1/rsfdata/MPD/' % when running on server
 elseif strcmp(getenv('HOSTNAME'),'')
    serv_path = '../' % running locally 
 else
   serv_path = '/Volumes/eol/fog1/rsfdata/MPD/' % 
 end

if strcmp(node,'MPD1')==1
  files = uipickfiles('prompt', 'select data files to process',  'FilterSpec', strcat(serv_path, 'wvdial_1_data/2020'));
  catalog = '/pub/incoming/catalog/operations';
elseif strcmp(node,'MPD2')==1
  files = uipickfiles('prompt', 'select data files to process', 'FilterSpec', strcat(serv_path, 'wvdial_2_data/2020'));
 catalog = '/pub/incoming/catalog/operations';
elseif strcmp(node,'MPD3')==1
  files = uipickfiles('prompt', 'select data files to process', 'FilterSpec', strcat(serv_path, 'wvdial_3_data/2020'));
 catalog = '/pub/incoming/catalog/operations';
elseif strcmp(node,'MPD4')==1
 files = uipickfiles('prompt', 'select data files to process', 'FilterSpec', strcat(serv_path, 'wvdial_4_data/2020'));
 catalog = '/pub/incoming/catalog/operations';
elseif strcmp(node,'MPD5')==1
 files = uipickfiles('prompt', 'select data files to process', 'FilterSpec', strcat(serv_path, 'wvdial_5_data/2020'));
 catalog = '/pub/incoming/catalog/operations';
end
j=1;


for j = 1:size(files,2)
    folder = (files{j});
    date = textscan(folder(end-5:end), '%6f'); date=date{1};  % read date of file
    if strcmp(node,'MPD1')==1
      MPD1_read_calvals % json format version of the above file
    elseif strcmp(node,'MPD2')==1
      MPD2_read_calvals % json format version of the above file
    elseif strcmp(node,'MPD3')==1
      MPD3_read_calvals % json format version of the above file
    elseif strcmp(node,'MPD4')==1
      MPD4_read_calvals % json format version of the above file
    elseif strcmp(node,'MPD5')==1
      MPD5_read_calvals % json format version of the above file
    end
    folder_in=folder;
    date_in = date;
    MPD_Analysis_function_NetCDF_v4(folder, date, MCS, write_data_folder, flag, node, wavemeter_offset,...
        profiles2ave, P0, switch_ratio, ave_time, timing_range_correction, blank_range, p_hour, catalog, Afterpulse_File)%
    %DIAL_Analysis_function_NetCDF(folder, date, MCS, write_data_folder, flag, node, ...
    %    profiles2ave, P0, switch_ratio, ave_time, timing_range_correction, blank_range, p_hour, catalog)%
    
    

end