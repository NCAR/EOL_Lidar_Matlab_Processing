clear all; close all

%location to write files
write_data_folder = uipickfiles('num',1,'out', 'char', 'prompt', ...
    'select folder to store data',  'FilterSpec', '/Volumes/documents/WV_DIAL_data/');
%write_data_folder = '/Volumes/documents/WV_DIAL_data/processed_data';
node = 'DIAL1';
    catalog = '/pub/incoming/catalog/perdigao';
%node = 'DIAL2'; %DIAL #2 HSRL data on NF channel, WV on FF channel
%    catalog = '/pub/incoming/catalog/operations';

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

flag.plot_data = 1;  % need to have this one to save the figs
flag.troubleshoot = 0; % shows extra plots used for troubleshooting
p_hour = 12; % hour to show troubleshooting profiles

ave_time.wv = 0.25; % averaging time (in minutes) for the water vapor 
ave_time.rb = 0.25; % averaging time (in minutes) for the relative backscatter
ave_time.gr = 0.25; % gridding time (in minutes) for the output files

if strcmp(node,'DIAL1')==1
    files = uipickfiles('prompt', 'select data files to process',  'FilterSpec', '/scr/eldora1/wvdial_1_data/');
else
  % flag.WS = 0; % use the surface weather station data to calcuate spectroscopy
  % flag.mark_gaps = 0; % sets gaps in data to NaNs
  files = uipickfiles('prompt', 'select data files to process', 'REFilter', 'FF', 'FilterSpec', '/scr/eldora1/wvdial_2_data/2018');
end
j=1;

for j = 1:size(files,2)
    folder = (files{j});
    date = textscan(folder(end-7:end), '%6f'); date=date{1};  % read date of file
    if strcmp(node,'DIAL1')==1
      %DIAL_settings_v2 % this loads the instrument settings which have changed over time
      read_dial1_calvals % json format version of the above file
    else
      %DIAL_2_settings_v1 % this loads the instrument settings which have changed over time
      read_dial2_calvals % json format version of the above file
    end
    folder_in=folder;
    DIAL_Analysis_v52_function(folder, MCS, write_data_folder, flag, node, ...
        profiles2ave, P0, switch_ratio, ave_time, timing_range_correction, blank_range, p_hour, catalog)%


end