 function[] = MPD_process_NetCDF_function_v2(save_quicklook, save_data, save_netCDF, save_catalog, channels, correction, node, daystr)
%  clear all; 
%  close all
%  start_date = '20220812';
%  save_quicklook=0; save_data=1; save_netCDF=0; save_catalog=0; channels = 'O2'; correction = 'AP_ON'; node='MPD01'; daystr=start_date; 

flag.save_quicklook = save_quicklook;  % save quicklook to local directory
flag.save_data = save_data;  % save files in matlab format
flag.save_netCDF = save_netCDF; % save files netCDF format
flag.save_catalog = save_catalog; % upload quicklook (and data) to field catalog

flag.ap_quick = 0; 
flag.mask_data = 1;  % mask applied to data based on error analysis threshold
flag.gradient_filter = 1;  % this is used to mask regions with 'high' backscatter gradients which tend to cause errors
flag.pileup = 1; % use pileup correction for detectors
flag.WS = 1; % use the surface weather station data to calcuate spectroscopy
flag.decimate = 0; % decimate all data to half the wv resoltuion
flag.int = 0; % interpolate nans in nanmoving_average
flag.mark_gaps = 1; % sets gaps in data to NaNs
flag.OF = 1; % correct for geometric overlap functions
%flag.near = near; %process the near range channel (or low gain)

flag.plot_data = 1;  % need to have this one to save the figs
flag.troubleshoot = 0; % shows extra plots used for troubleshooting
p_hour = 20; % hour to show troubleshooting profiles
ave_time.wv = 10; %10.0; % averaging time (in minutes) for the water vapor and O2 
ave_time.rb = 1; %5.0; % averaging time (in minutes) for the relative backscatter
ave_time.gr = 1.0; % gridding time (in minutes) for the output files (native is 2 sec)
%ave_range.gr = 75; % grid data to this range (useful for SmartSwitch tests)

if strcmp(getenv('HOSTNAME'),'eol-smaug.eol.ucar.edu') == 1  % when running on fog server
   serv_path = '/export/smaug1/rsfdata/MPD/'; 
   cal_serv_path = '/export/smaug1/rsfdata/MPD/calibration/'
   cal_serv_path = '/home/rsfdata/Processing/Python/'
elseif strcmp(getenv('HOSTNAME'),'tikal.eol.ucar.edu') ==1  % when running on tikal server
    serv_path = '/export/smaug/rsfdata/MPD/'; 
    cal_serv_path = '/h/eol/spuler/' 
elseif strcmp(getenv('HOSTNAME'),'') % running locally 
    serv_path = '/Volumes/smaug1/rsfdata/MPD/'; 
%     serv_path = '/Volumes/eol/cirrus2/rsfdata/projects/mpd/'; 
    cal_serv_path = '../' 
    %cal_serv_path = '/Users/spuler/Documents/GitHub';  
else % running locally but use calibration from fog server
   serv_path = '/Volumes/smaug1/rsfdata/MPD/'; % 
   cal_serv_path = '/Volumes/eol/smaug1/rsfdata/MPD/calibration/' % 
end

nodeStr = extractAfter(node, 'MPD')
year_folder = cell2mat(textscan(daystr,'%4c', 1));
if nodeStr == "06"
 files = strcat(serv_path, 'adihsrl', '_data/', year_folder,'/', daystr);
else
 files = strcat(serv_path, 'mpd_', nodeStr, '_data/', year_folder,'/', daystr);
end

catalog = '/pub/incoming/catalog/operations';

folder = files;
date = textscan(folder(end-5:end), '%6f'); date=date{1};  % read date of file
folder_in=folder;
date_in = date;
serial_date = datenum(num2str(date_in),'yymmdd');
MPD_read_calvals % read in the calvals
read_time_in
profiles2ave.wv = 2*round(((ave_time.wv*60/read_time_in)+1)/2)   
profiles2ave.rb = 2*round(((ave_time.rb*60/read_time_in)+1)/2)

% read in all the data
 if  strcmp(channels,'WV') == 1 
   %[data_on, data_off, MCS] = MPD_File_Retrieval_NetCDF_v5(flag, MCS, folder_in, read_time_in); %use to read binary data (bin number passed in) 
   [data_wv_on, data_wv_off, MCS] = MPD_File_Retrieval_NetCDF_v5(flag, MCS, folder_in, read_time_in); %use to read binary data (bin number passed in) 
 end
 if  strcmp(channels,'O2') == 1 
   [data_wv_on, data_wv_off, data_O2_on_comb, data_O2_off_comb, data_O2_on_mol, data_O2_off_mol, MCS] = MPD_File_Retrieval_NetCDF_v6(flag, MCS, folder_in, read_time_in); %use to read binary data (bin number passed in) 
 end
 if  strcmp(channels,'HSRL') == 1 
   [data_mol, data_comb, MCS] = MPD_File_Retrieval_NetCDF_HSRL(flag, MCS, folder_in); %use to read binary data (bin number passed in) 
 end
 
% Pause here to create an afterpulse file if desired (run MPS_afterpulse_cals_v2.m)
 
% process the main ch without afterpulse correction 
%if strcmp(channels,'ALL') == 1 || strcmp(channels,'WV') == 1 
if (strcmp(channels,'ALL') == 1 || strcmp(channels,'WV') == 1) && strcmp(correction,'AP_OFF') == 1 
  write_data_folder = strcat(serv_path, 'mpd_', nodeStr, '_processed_data/Matlab'); 
  flag.near = 0; flag.afterpulse = 0; 
  MPD_Analysis_function_NetCDF_v5(data_wv_on, data_wv_off, folder, date, MCS, write_data_folder, flag, node, wavemeter_offset,...
        profiles2ave, P0, switch_ratio, ave_time, timing_range_correction, blank_range, p_hour, catalog, Afterpulse_File, MPD_elevation, cal_serv_path)%

 %[N_WV, N_WV_error] = MPD_WV_analysis_function_v2(data_wv_on, data_wv_off, folder_in, date_in, MCS, write_data_folder, flag, node, wavemeter_offset,...
 %        profiles2ave, switch_ratio, ave_time, timing_range_correction, blank_range, p_hour, catalog, Afterpulse_File, MPD_elevation);
end

% % process the oxygen channels
 if (strcmp(channels,'ALL') == 1 || strcmp(channels,'O2') == 1)  && strcmp(correction,'AP_OFF') == 1 
   write_data_folder = strcat(serv_path, 'mpd_', nodeStr, '_processed_data/Matlab'); 
   flag.near = 0; flag.afterpulse = 0; 
   gates2ave = 1; %number of gates to average
    data_on = data_O2_on_comb;
    data_off = data_O2_off_comb;
  [O2_online_comb, O2_offline_comb, range, RB_comb, time_comb, Surf_T, Surf_P, O2_on_wavelength, O2_off_wavelength] =  MPD_Analysis_function_O2_v1(data_O2_on_comb, data_O2_off_comb, folder, date, MCS, write_data_folder, flag, node, wavemeter_offset,...
         profiles2ave, switch_ratio, ave_time, timing_range_correction, blank_range, p_hour, gates2ave, Afterpulse_File, cal_serv_path);%
  [O2_online_mol, O2_offline_mol, range, RB_mol,  time_mol] =  MPD_Analysis_function_O2_v1(data_O2_on_mol, data_O2_off_mol, folder, date, MCS, write_data_folder, flag, node, wavemeter_offset,...
         profiles2ave, switch_ratio, ave_time, timing_range_correction, blank_range, p_hour, gates2ave, Afterpulse_File, cal_serv_path);%
  [T, P, BSR, RD, HSRLMolecular_scan_wavelength, const, beta_m_profile] = Process_HSRL_K_data(O2_online_comb, O2_offline_comb,...
         O2_online_mol,O2_offline_mol, time_comb, range, Surf_T, Surf_P, flag, node, daystr, Receiver_Scan_File, write_data_folder,cal_serv_path, receiver_scale_factor);
  [N_WV, N_WV_error] = MPD_WV_analysis_function_v1(data_wv_on, data_wv_off, folder_in, date_in, MCS, write_data_folder, flag, node, wavemeter_offset,...
         profiles2ave, T, P, switch_ratio, ave_time, timing_range_correction, blank_range, p_hour, catalog, Afterpulse_File, MPD_elevation, cal_serv_path);
  O2_absorption(const, T, P, O2_online_comb, O2_offline_comb, ...
         time_comb, range, BSR, RD, HSRLMolecular_scan_wavelength, N_WV, beta_m_profile, O2_on_wavelength, node, daystr, write_data_folder, flag);
 end 

% % process the oxygen channels
 if (strcmp(channels,'HSRL') == 1)  && strcmp(correction,'AP_OFF') == 1 
   write_data_folder = strcat(serv_path, 'mpd_', nodeStr, '_processed_data/Matlab'); 
   flag.near = 0; flag.afterpulse = 0; 
   gates2ave = 1; %number of gates to average
  [HSRL_comb, ~, range, RB_comb, time_comb, Surf_T, Surf_P, O2_on_wavelength, O2_off_wavelength] =  MPD_Analysis_function_O2_v1(data_comb, data_comb, folder, date, MCS, write_data_folder, flag, node, wavemeter_offset,...
         profiles2ave, switch_ratio, ave_time, timing_range_correction, blank_range, p_hour, gates2ave, Afterpulse_File, cal_serv_path);%
  [HSRL_mol, ~, range, RB_mol,  time_mol] =  MPD_Analysis_function_O2_v1(data_mol, data_mol, folder, date, MCS, write_data_folder, flag, node, wavemeter_offset,...
         profiles2ave, switch_ratio, ave_time, timing_range_correction, blank_range, p_hour, gates2ave, Afterpulse_File, cal_serv_path);%
  [T, P, BSR, RD, HSRLMolecular_scan_wavelength, const, beta_m_profile] = Process_HSRL_K_data(HSRL_comb, HSRL_comb,...
         HSRL_mol,HSRL_mol, time_comb, range, Surf_T, Surf_P, flag, node, daystr, Receiver_Scan_File, write_data_folder,cal_serv_path, receiver_scale_factor);
 end 


 
 
%% change plots when runnning locally to focus on lower ranges 
%ylim([0 3])
%grid on
%grid minor


% process the near/low ch
%if strcmp(channels,'ALL') == 1 || strcmp(channels,'WVLow') == 1  
%  flag.near = 1; flag.afterpulse = 0; 
%  if strcmp(node,'MPD04') == 1 % MPD04 is using a low range channel
%    blank_range = 187.5; % low range 
%  end
%  MPD_Analysis_function_NetCDF_v5(data_near_on, data_near_off, folder, date, MCS, write_data_folder, flag, node, wavemeter_offset,...
%        profiles2ave, P0, switch_ratio, ave_time, timing_range_correction, blank_range, p_hour, catalog, Afterpulse_File)%
%end


% process the main ch with afterpulse correction 
if (strcmp(channels,'ALL') == 1 || strcmp(channels,'WV') == 1) && strcmp(correction,'AP_ON') == 1 
  write_data_folder = strcat(serv_path, 'mpd_', nodeStr, '_processed_data/Matlab/afterpulse');  
  flag.near = 0; flag.afterpulse = 1; 
%   if ((strcmp(node,'MPD04') == 1) && (serial_date >= 737902)) || ... 
%      ((strcmp(node,'MPD03') == 1) && (serial_date >= 737916)) % MPD04 and 03 use combined low range channels
%      blank_range = 150; % low range 
%   end
  MPD_Analysis_function_NetCDF_v5(data_wv_on, data_wv_off, folder, date, MCS, write_data_folder, flag, node, wavemeter_offset,...
        profiles2ave, P0, switch_ratio, ave_time, timing_range_correction, blank_range, p_hour, catalog, Afterpulse_File, MPD_elevation, cal_serv_path)%
end


% process the main O2 with afterpulse correction 
if (strcmp(channels,'ALL') == 1 || strcmp(channels,'O2') == 1) && strcmp(correction,'AP_ON') == 1 
  write_data_folder = strcat(serv_path, 'mpd_', nodeStr, '_processed_data/Matlab/afterpulse');  
  flag.near = 0; flag.afterpulse = 1; 
  gates2ave = 1; %number of gates to average
  data_on = data_O2_on_comb;
  data_off = data_O2_off_comb;
  flag.molecular = 0;
  [O2_online_comb, O2_offline_comb, range, RB_comb, time_comb, Surf_T, Surf_P, O2_on_wavelength, O2_off_wavelength] =  MPD_Analysis_function_O2_v1(data_O2_on_comb, data_O2_off_comb, folder, date, MCS, write_data_folder, flag, node, wavemeter_offset,...
         profiles2ave, switch_ratio, ave_time, timing_range_correction, blank_range, p_hour, gates2ave, Afterpulse_File,cal_serv_path);%
  flag.molecular = 1;
  [O2_online_mol, O2_offline_mol, range, RB_mol,  time_mol] =  MPD_Analysis_function_O2_v1(data_O2_on_mol, data_O2_off_mol, folder, date, MCS, write_data_folder, flag, node, wavemeter_offset,...
         profiles2ave, switch_ratio, ave_time, timing_range_correction, blank_range, p_hour, gates2ave, Afterpulse_File,cal_serv_path);%
  [T, P, BSR, RD, HSRLMolecular_scan_wavelength, const, beta_m_profile] = Process_HSRL_K_data(O2_online_comb, O2_offline_comb,...
         O2_online_mol,O2_offline_mol, time_comb, range, Surf_T, Surf_P, flag, node, daystr, Receiver_Scan_File, write_data_folder,cal_serv_path, receiver_scale_factor);
  [N_WV, N_WV_error] = MPD_WV_analysis_function_v1(data_wv_on, data_wv_off, folder_in, date_in, MCS, write_data_folder, flag, node, wavemeter_offset,...
         profiles2ave, T, P, switch_ratio, ave_time, timing_range_correction, blank_range, p_hour, catalog, Afterpulse_File, MPD_elevation, cal_serv_path);
  O2_absorption(const, T, P, O2_online_comb, O2_offline_comb, ...
         time_comb, range, BSR, RD, HSRLMolecular_scan_wavelength, N_WV, beta_m_profile, O2_on_wavelength, node, daystr, write_data_folder, flag);
end

% process the near/low ch with afterpulse correction 
%if (strcmp(channels,'ALL') == 1 || strcmp(channels,'WVLow') == 1) && strcmp(correction,'AP_ON') == 1    
%  flag.near = 1; flag.afterpulse = 1;
%  if strcmp(node,'MPD04') == 1 % MPD04 is using a low range channel
%    blank_range = 150; 
%  end
%  MPD_Analysis_function_NetCDF_v5(data_near_on, data_near_off, folder, date, MCS, write_data_folder, flag, node, wavemeter_offset,...
%        profiles2ave, P0, switch_ratio, ave_time, timing_range_correction, blank_range, p_hour, catalog, Afterpulse_File)%
%end
    
end