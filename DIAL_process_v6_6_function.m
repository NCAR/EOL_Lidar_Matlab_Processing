function[] = DIAL_process_v6_6_function(files, save_quicklook, save_data, save_netCDF, save_catalog, node)
%clear all; 
close all

%location to write files
if strcmp(node,'NCAR')==1
  write_data_folder = '/scr/eldora1/wvdial_1_processed_data';
else
  write_data_folder = '/scr/eldora1/wvdial_2_processed_data';
end

flag.save_quicklook = save_quicklook; %1;  % save quicklook to local directory
flag.save_data = save_data; %1;  % save files in matlab format
flag.save_netCDF = save_netCDF; %1; % save files netCDF format
flag.save_catalog = save_catalog; %0; % upload quicklook (and data) to field catalog

flag.mask_data = 1;  % mask applied to data based on error analysis threshold
flag.gradient_filter = 1;  % this is used to mask regions with 'high' backscatter gradients which tend to cause errors
flag.pileup = 1; % use pileup correction for detectors
flag.WS = 1; % use the surface weather station data to calcuate spectroscopy
flag.decimate = 1; % decimate all data to half the wv resoltuion
flag.int = 0; % interpolate nans in nanmoving_average
flag.mark_gaps = 1; % sets gaps in data to NaNs

flag.plot_data = 0;  % need to have this one to save the figs
flag.troubleshoot = 0; % shows extra plots used for troubleshooting
p_hour = 12; % hour to show troubleshooting profiles

ave_time.wv = 5.0; % averaging time (in minutes) for the water vapor 
ave_time.rb = 0.5; % averaging time (in minutes) for the relative backscatter
ave_time.gr = 0.5; % gridding time (in minutes) for the output files

%files = uipickfiles('prompt', 'select data files to process',  'FilterSpec', '/scr/eldora1/h2o_data/');
j=1;

%for j = 1:size(files,2)
    %folder = (files{j});
    folder = files; 
    date = textscan(folder(end-7:end), '%6f'); date=date{1};  % read date of file
    if strcmp(node,'NCAR')==1
      DIAL_settings_v2 % this loads the instrument settings which have changed over time
    else
      MSU_DIAL_settings_v1 % this loads the instrument settings which have changed over time  
    end
    folder_in=folder;
    DIAL_Analysis_v52_function(folder, MCS, write_data_folder, flag, node, ...
        profiles2ave, P0, switch_ratio, ave_time, timing_range_correction, blank_range, p_hour)%
%end