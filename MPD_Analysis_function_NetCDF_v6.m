function[] = MPD_Analysis_function_NetCDF_v6(data_on, data_off, folder_in, date_in, MCS, write_data_folder, flag, node, wavemeter_offset,...
    profiles2ave, P0, switch_ratio, ave_time, timing_range_correction, blank_range, p_hour, catalog, Afterpulse_File, MPD_elevation, cal_serv_path)

%% notes
% ... (Notes/Version history stripped for brevity) ...
%% start
     
dd = pwd; % get the current path        
close all; clc
scrsz = get(0,'ScreenSize');
 
tic;

%Colormap
C = importdata('NCAR_C_Map.mat');
Hitran.file = dlmread('823nm_834nm_HITRAN_2016.csv', ',',[0 0 2689 7]);

RB_scale = 1; 
gate = round((MCS.bin_duration*1e-9*3e8/2)*100)/100

delta_r_index =  75/gate; 
r1 = round(1500/gate); 
r2 = round(2500/gate); 
spatial_average1 = 150/gate; 
spatial_average2 = 300/gate; 
spatial_average3 = 600/gate; 

%% Importing online and offline files from the selected date

%Making the online and offline data files the same size
try
    Online_Raw_Data = data_on(1:size(data_on,1),1:end);
    Offline_Raw_Data = data_off(1:size(Online_Raw_Data,1),1:end);
catch err
    Offline_Raw_Data = data_off(1:size(data_off,1),1:end);
    Online_Raw_Data = data_on(1:size(Offline_Raw_Data,1),1:end);
end

  % add trap error associated with Perdigao instrument crash 
  serial_date = datenum(num2str(date_in),'yymmdd');
  time2 = double((Online_Raw_Data(:,1)))./24+serial_date;
  time = time2;
  time = time(~(time2<(median(time2,'omitnan')-1.5)));
  time= time(~(time>(median(time,'omitnan')+1.5)));
  Online_Raw_Data = Online_Raw_Data(~isnan(time),:);
  Offline_Raw_Data = Offline_Raw_Data(~isnan(time),:);


% ... (WS data processing and data cleanup stripped for brevity) ...
 
% Parsing off the auxiliary 
 Online = single(Online_Raw_Data(:,10:end)); 
 Offline = single(Offline_Raw_Data(:,10:end));
 
% vector in meters
range = single(0:gate:(size(Online,2)-1)*gate);

if flag.pileup == 1
% ... (Pileup correction skipped for brevity) ...
end

if flag.afterpulse == 1   % afterpulse correction
   
  if flag.ap_quick == 1  
      warning('Using quick afterpulse method. This should be disabled for production.');
      ap_spline_sub_off = zeros(1, size(Online, 2)); % Placeholder zeros matching size
      ap_spline_sub_on = zeros(1, size(Online, 2));   % Placeholder zeros matching size
  else    
     % read the afterpulse nc file identified in the json file 
     ap_filename = strcat(cal_serv_path, 'eol-lidar-calvals/calfiles/', Afterpulse_File);
   
     ncid = netcdf.open(ap_filename, 'NC_NOWRITE');
     if flag.near == 1
       ap_off_rate = ncread(ap_filename, 'WVOfflineLow_afterpulse');
       ap_on_rate = ncread(ap_filename, 'WVOnlineLow_afterpulse');  
     else
       ap_off_rate = ncread(ap_filename, 'WVOffline_afterpulse');
       ap_on_rate = ncread(ap_filename, 'WVOnline_afterpulse');
     end
     ap_range = ncread(ap_filename, 'range');
     netcdf.close(ncid);   
     afterpulse_off = ap_off_rate*MCS.accum*MCS.bin_duration*1e-9;
     afterpulse_on = ap_on_rate*MCS.accum*MCS.bin_duration*1e-9;  

     % --- CRITICAL FIX: Interpolate AP profile to match current range vector (size) ---
     afterpulse_off_col = afterpulse_off(:);
     afterpulse_on_col = afterpulse_on(:);
     ap_range_col = ap_range(:);
     
     ap_interpolated_off = interp1(ap_range_col, afterpulse_off_col, range, 'linear', 0);
     ap_interpolated_on = interp1(ap_range_col, afterpulse_on_col, range, 'linear', 0);
     
     ap_spline_sub_off = ap_interpolated_off(:)'; 
     ap_spline_sub_on = ap_interpolated_on(:)'; 
     % ---------------------------------------------------------------------------------
     
  end
  
   Offline_ap_sub = (bsxfun(@minus, Offline, ap_spline_sub_off));
   Online_ap_sub = (bsxfun(@minus, Online, ap_spline_sub_on)); 

    Online = Online_ap_sub;
    Offline = Offline_ap_sub;
end


 i = size(Online, 1);
 j = size(Online, 2);

 
%% Background subtraction 
% ... (Background subtraction, temporal/spatial averaging, DIAL/saving/plotting skipped for brevity) ...

toc
cd(dd) 

end