function[] = MPD_process_NetCDF_function_v2(save_quicklook, save_data, save_netCDF, save_catalog, channels, correction, node, daystr)
%clear all; 
%close all
%start_date = '20200409';
%save_quicklook=0; save_data=1; save_netCDF=0; save_catalog=0; channels = 'WV'; correction = 'AP_ON'; node='MPD04'; daystr=start_date; 

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
%flag.near = near; %process the near range channel (or low gain)
%flag.afterpulse = afterpulse; % correct for afterpulsing (in progress on MPD 3 & 4 only)

flag.plot_data = 1;  % need to have this one to save the figs
flag.troubleshoot = 0; % shows extra plots used for troubleshooting
p_hour = 14.42; % hour to show troubleshooting profiles

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

nodeStr = extractAfter(node, 'MPD');
files = strcat(serv_path, 'mpd_', nodeStr, '_data/2020/', daystr);
catalog = '/pub/incoming/catalog/operations';

folder = files;
date = textscan(folder(end-5:end), '%6f'); date=date{1};  % read date of file
MPD_read_calvals % read in the calvals
folder_in=folder;
date_in = date;
read_time_in = 2; % set read data in time increments of seconds (default it 2sec) 
  
% read in all the data
 [data_on, data_off, data_near_on, data_near_off, MCS] = MPD_File_Retrieval_NetCDF_v5(flag, MCS, folder_in, read_time_in); %use to read binary data (bin number passed in) 

 
% process the main ch without afterpulse correction 
if strcmp(channels,'ALL') == 1 || strcmp(channels,'WV') == 1 
  write_data_folder = strcat(serv_path, 'mpd_', nodeStr, '_processed_data/Matlab'); 
  flag.near = 0; flag.afterpulse = 0; 
  MPD_Analysis_function_NetCDF_v5(data_on, data_off, folder, date, MCS, write_data_folder, flag, node, wavemeter_offset,...
        profiles2ave, P0, switch_ratio, ave_time, timing_range_correction, blank_range, p_hour, catalog, Afterpulse_File)%
end

% % experiment to combine the low and high range channels 
%     gate = round((MCS.bin_duration*1e-9*3e8/2)*100)/100;
% 
%     data_on_header = (data_on(:,1:9));
%     data_on_data =  (data_on(:,10:end));
%     data_on_background = mean(data_on_data(:,end-round(1200/gate):end),2)-0; % final ~1.2km background
%     data_on_back_sub = (bsxfun(@minus, data_on_data, data_on_background));
%     data_near_on_data =  (data_near_on(:,10:end));
%     data_near_on_background = mean(data_near_on_data(:,end-round(1200/gate):end),2)-0; % final ~1.2km background
%     data_near_on_back_sub = (bsxfun(@minus, data_near_on_data, data_near_on_background));
%     data_comb_on = [data_on_header, (data_on_back_sub + data_near_on_back_sub)];
% 
%     data_off_header = (data_off(:,1:9));
%     data_off_data =  (data_off(:,10:end));
%     data_off_background = mean(data_off_data(:,end-round(1200/gate):end),2)-0; % final ~1.2km background
%     data_off_back_sub = (bsxfun(@minus, data_off_data, data_off_background));
%     data_near_off_data =  (data_near_off(:,10:end));
%     data_near_off_background = mean(data_near_off_data(:,end-round(1200/gate):end),2)-0; % final ~1.2km background
%     data_near_off_back_sub = (bsxfun(@minus, data_near_off_data, data_near_off_background));
%     data_comb_off = [data_off_header, (data_off_back_sub + data_near_off_back_sub)];
% 
%     % process the combined low and high ch without afterpulse correction 
%     if strcmp(channels,'ALL') == 1 || strcmp(channels,'WV') == 1 
%       if strcmp(node,'MPD4') == 1 % MPD4 is using a low range channel
%         blank_range = 187.5; % low range 
%       end
%         write_data_folder = strcat(serv_path, 'mpd_', nodeStr, '_processed_data/Matlab'); 
%       flag.near = 0; flag.afterpulse = 0; 
%       MPD_Analysis_function_NetCDF_v5(data_comb_on, data_comb_off, folder, date, MCS, write_data_folder, flag, node, wavemeter_offset,...
%             profiles2ave, P0, switch_ratio, ave_time, timing_range_correction, blank_range, p_hour, catalog, Afterpulse_File)%
%     end


% ylim([0 3])
% grid on
% grid minor


% process the near/low ch
if strcmp(channels,'ALL') == 1 || strcmp(channels,'WVLow') == 1  
  flag.near = 1; flag.afterpulse = 0; 
  if strcmp(node,'MPD04') == 1 % MPD4 is using a low range channel
    blank_range = 187.5; % low range 
  end
  MPD_Analysis_function_NetCDF_v5(data_near_on, data_near_off, folder, date, MCS, write_data_folder, flag, node, wavemeter_offset,...
        profiles2ave, P0, switch_ratio, ave_time, timing_range_correction, blank_range, p_hour, catalog, Afterpulse_File)%
end

% process the main ch with afterpulse correction 
if (strcmp(channels,'ALL') == 1 || strcmp(channels,'WV') == 1) && strcmp(correction,'AP_ON') == 1 
  write_data_folder = strcat(serv_path, 'mpd_', nodeStr, '_processed_data/Matlab/afterpulse');   
  flag.near = 0; flag.afterpulse = 1; 
  MPD_Analysis_function_NetCDF_v5(data_on, data_off, folder, date, MCS, write_data_folder, flag, node, wavemeter_offset,...
        profiles2ave, P0, switch_ratio, ave_time, timing_range_correction, blank_range, p_hour, catalog, Afterpulse_File)%
end

% process the near/low ch with afterpulse correction 
if (strcmp(channels,'ALL') == 1 || strcmp(channels,'WVLow') == 1) && strcmp(correction,'AP_ON') == 1    
  flag.near = 1; flag.afterpulse = 1;
  if strcmp(node,'MPD04') == 1 % MPD4 is using a low range channel
    blank_range = 150; 
  end
  MPD_Analysis_function_NetCDF_v5(data_near_on, data_near_off, folder, date, MCS, write_data_folder, flag, node, wavemeter_offset,...
        profiles2ave, P0, switch_ratio, ave_time, timing_range_correction, blank_range, p_hour, catalog, Afterpulse_File)%
end
    
end