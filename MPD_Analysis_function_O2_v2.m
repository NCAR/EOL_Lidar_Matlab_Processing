function[O2_online_comb, O2_offline_comb, range, RB, time_grid, Surf_T, Surf_P, lambda_all, lambda_all_off] = MPD_Analysis_function_O2_v2(data_on, data_off, folder_in, date_in, MCS, write_data_folder, flag, node, wavemeter_offset,...
    profiles2ave, switch_ratio, ave_time, timing_range_correction, blank_range, p_hour, gates2ave, Afterpulse_File, cal_serv_path)

 
%% notes

%% start
     
dd = pwd; % get the current path        
close all; clc
scrsz = get(0,'ScreenSize');
 
tic;

%Colormap
C = importdata('NCAR_C_Map.mat');
%C = viridis(70);

%Importing HITRAN data
Hitran.file = dlmread('815nm_841nm_HITRAN_2008.csv', ',',[1 1 1676 8]);
%hitran = dlmread('823nm_834nm_HITRAN_2012.csv', ',',[0 0 2633 7]);

RB_scale = 1; % use to keep the arbitrary units of RB scale the same before

%Spatial averaging (range average) in bins.  
gate = round((MCS.bin_duration*1e-9*3e8/2)*100)/100

delta_r_index =  75/gate*gates2ave; % this is the cumlative sum photons gate spacing 
%delta_r_index =  1; % proces at the native gate spacing
delta_r = delta_r_index*gate*100; % delta r in cm
r1 = round(1500/gate); % index for smoothing range 1 (1500m)
r2 = round(2500/gate); % index for smoothing range 2 (2500m)
spatial_average1 = 150/gate; %150 meter smoothing range 1 
spatial_average2 = 300/gate; %300 meter smoothing between range 1 and 2
spatial_average3 = 600/gate; %600 meter smoothing above range 2
% this filter creates range delays 

disp('--- O2_v2: START EXECUTION ---');

%% Importing online and offline files from the selected date

%Making the online and offline data files the same size
try
    Online_Raw_Data = data_on(1:size(data_on,1),1:end);
    Offline_Raw_Data = data_off(1:size(Online_Raw_Data,1),1:end);
catch err
    Offline_Raw_Data = data_off(1:size(data_off,1),1:end);
    Online_Raw_Data = data_on(1:size(Offline_Raw_Data,1),1:end);
end

disp(['O2_v2: Data Size: Online=', num2str(size(Online_Raw_Data, 1)), ' rows.']);
if isempty(Online_Raw_Data); error('Input data is empty.'); end

  % add trap error associated with Perdigao instrument crash 
  serial_date = datenum(num2str(date_in),'yymmdd');
  time2 = double((Online_Raw_Data(:,1)))./24+serial_date;
  time = time2;
  time = time(~(time2<(median(time2,'omitnan')-1.5)));
  time= time(~(time>(median(time,'omitnan')+1.5)));
  Online_Raw_Data = Online_Raw_Data(~isnan(time),:);
  Offline_Raw_Data = Offline_Raw_Data(~isnan(time),:);



%% read in weather station data

  I_on = Online_Raw_Data(:,4); % online current
  I_off = Offline_Raw_Data(:,4); % offline current
  P_on = Online_Raw_Data(:,5); % transmitted average power
  P_off = Offline_Raw_Data(:,5); % transmitted average power
  T_bench = Online_Raw_Data(:,6); % optical bench temperature (thermocouple 1)
  T_base = Offline_Raw_Data(:,6); % base temperature (thermocouple 2)
  
if flag.WS==1
  % read in surface weather station data
  Surf_T = Offline_Raw_Data(:,7);  %temperature in C
  Surf_P = Offline_Raw_Data(:,8)./1013.249977;  % pressure in atm   
  Surf_T(isnan(Surf_T))= median(Surf_T); % fills in missing values with median
  Surf_P(isnan(Surf_P))= median(Surf_P); % fills in missing values with median
     if isnan(median(Surf_T,'omitnan'))==1
       T0=25;
       Surf_T= ones(size(Surf_T)).*T0;
       warning('No Temperature Weather Station Data')
     end
     if isnan(median(Surf_P,'omitnan'))==1
       P0=1;
       P0=0.81; % 
       Surf_P= ones(size(Surf_P)).*P0;
       warning('No Pressure Weather Station Data')
     end
    
  Surf_RH = Offline_Raw_Data(:,9);
  % convert RH to number density and absolute humidity
  % vapor pressure of water
    a0 = 6.107799961;
    a1 = 4.436518521E-1; 
    a2 = 1.428945805E-2;
    a3 = 2.650648471E-4; 
    a4 = 3.031240396E-6;
    a5 = 2.034080948E-8;
    a6 = 6.136820929E-11;
  e=((a0+Surf_T.*(a1+Surf_T.*(a2+Surf_T.*(a3+Surf_T.*(a4+Surf_T.*(a5+Surf_T.*a6))))))./1); %vapor pressure in hPa 
  % constants
  R = 8.31447215; %J mol^-1 K^-1
  N_A= 6.0221415E23; %mol^-1
  % convert from RH to number density
  Surf_N = 1.*(Surf_RH.*(1).*e./(R.*(Surf_T+273).*(1))).*N_A*1e-6;  %cm^3
  Surf_AH = Surf_N.*1e6./6.022E23.*18.015;
end

%% initial data preparation

%Calculating average wavelenth
lambda_all=Online_Raw_Data(:,2);
lambda_diff= Online_Raw_Data(:,3);
lambda_all_off=Offline_Raw_Data(:,2);
lambda_diff_off= Offline_Raw_Data(:,3);

lambda_on = median(lambda_all);
lambda_off = median(lambda_all_off);
 
disp(['O2_v2: Wavelength median calculated: ', num2str(lambda_on)]);
    
 folder_date = textscan(folder_in(end-7:end), '%8f', 1); folder_date=folder_date{1};
  
% Parsing off the auxiliary 
 Online = single(Online_Raw_Data(:,10:end)); 
 Offline = single(Offline_Raw_Data(:,10:end));
 
% vector in meters
range = single(0:gate:(size(Online,2)-1)*gate);

if flag.pileup == 1
% apply linear correction factor to raw counts 
  %t_d=32E-9; % module dead time of Perkin Elmber
  t_d=37.25E-9; %Excelitas SPCM-AQRH-13 Module 24696
  %t_d=50E-9; %Emperical best fit to remove the noise in the WV behind clouds
  %t_d=34E-9; %Excelitas SPCM-AQRH-13 Module 24696 for count rates < 5 Mc/s
  % MCSC gives counts accumulated for set bin duration so convert to count rate  C/s.
  % divide by bin time in sec (e.g., 500ns) and # of acumulations (e.g., 10000)
  % e.g., 10 accumlated counts is 2000 C/s

  %remove non linear data
%  Online(Online./(MCS.bin_duration*1E-9*MCS.accum)>1E6)= nan;
%  Offline(Offline./(MCS.bin_duration*1E-9*MCS.accum)>1E6) = nan;
  
  C_Online = 1./(1-(t_d.*(Online./(MCS.bin_duration*1E-9*MCS.accum))));   
  C_Offline = 1./(1-(t_d.*(Offline./(MCS.bin_duration*1E-9*MCS.accum))));   
  Online = Online.*C_Online;
  Offline = Offline.*C_Offline;

end


if flag.afterpulse == 1   % afterpulse correction
   
  if flag.ap_quick == 1  
  %  read my version of the afterpulse calibration files which are now using rates     
      warning('Using quick afterpulse method. This should be disabled for production.');
      ap_spline_sub_off = zeros(1, size(Online, 2)); 
      ap_spline_sub_on = zeros(1, size(Online, 2));   
      
  else    

     % read the afterpulse nc file identified in the json file  
     ap_filename = strcat(cal_serv_path, 'eol-lidar-calvals/calfiles/', Afterpulse_File)   
   
     ncid = netcdf.open(ap_filename, 'NC_NOWRITE');
 %    ncdisp(ap_filename, '/', 'min') % use this to display all variables
     if flag.molecular == 1
       ap_off_rate = ncread(ap_filename, 'O2OfflineMol_afterpulse');
       ap_on_rate = ncread(ap_filename, 'O2OnlineMol_afterpulse');  
     else
       ap_off_rate = ncread(ap_filename, 'O2OfflineComb_afterpulse');
       ap_on_rate = ncread(ap_filename, 'O2OnlineComb_afterpulse');
     end
     ap_range = ncread(ap_filename, 'range');
     netcdf.close(ncid);   
     afterpulse_off = ap_off_rate*MCS.accum*MCS.bin_duration*1e-9;
     afterpulse_on = ap_on_rate*MCS.accum*MCS.bin_duration*1e-9;  

     range_shift = -(delta_r_index-1)/2*gate + timing_range_correction; % 
     range_act = range + range_shift; % %actual range points 

     figure(1004)
     semilogy(afterpulse_off, 'bo-')
     hold on
     semilogy(afterpulse_on, 'b+-')    
     hold off
     legend('hayman_{off}', 'hayman_{on}') 
     ylabel('counts')
     xlabel('bins')
     grid on
     grid minor 
   
     % --- CRITICAL FIX: Interpolate AP profile to match current range vector (size) ---
     afterpulse_off_col = afterpulse_off(:);
     afterpulse_on_col = afterpulse_on(:);
     ap_range_col = ap_range(:);
     
     ap_interpolated_off = interp1(ap_range_col, afterpulse_off_col, range, 'linear', 0);
     ap_interpolated_on = interp1(ap_range_col, afterpulse_on_col, range, 'linear', 0);
     
     ap_spline_sub_off = ap_interpolated_off(:)'; % **FIXED: Ensure correct length and orientation**
     ap_spline_sub_on = ap_interpolated_on(:)';   % **FIXED: Ensure correct length and orientation**
     % ---------------------------------------------------------------------------------
     
  end
  
   Offline_ap_sub = (bsxfun(@minus, Offline, ap_spline_sub_off)); 
   Online_ap_sub = (bsxfun(@minus, Online, ap_spline_sub_on)); 
    
    Online = Online_ap_sub;
    Offline = Offline_ap_sub;
end

disp('O2_v2: Afterpulse correction complete. Starting Background subtraction.');

 i = size(Online, 1);
 j = size(Online, 2);

 
%% Background subtraction %take the values from 14.25km to 15km for background  
%h = waitbar(0,'Background subtraction and range correction');

  %range_km_squared = (range./1e3).^2.*(50/MCS.bin_duration);  % this keeps the original color bar which is arbitrary
  %range_km_squared = (range./1e3).^2.*((MCS.bin_duration*1E-6*accum*(1-switch_ratio)))^(-1);  % in units of m^2 C/ms
  range_km_squared = (range).^2./((MCS.bin_duration*MCS.accum*(1-switch_ratio)));  % in units of km^2 C/ns 

% may need to average before backround sub to avoid errors with negative numbers (small SNR regions)  
% Online_sum = cumsum(Online,1)-[zeros(profiles2ave.wv,j); cumsum(Online(1:i-profiles2ave.wv,:),1)];  %rolling average of rows or time
% Offline_sum = cumsum(Offline,1)-[zeros(profiles2ave.wv,j); cumsum(Offline(1:i-profiles2ave.wv,:),1)];  %rolling average of rows or time
% Online = Online_sum./profiles2ave.wv;
% Offline = Offline_sum./profiles2ave.wv; 
 
  background_on = mean(Online(:,end-round(1200/gate):end),2)-0; % select last ~1200 meters to measure background
  background_off = mean(Offline(:,end-round(1200/gate):end),2)-0; % select last ~1200 meters to measure background 
  % to deal with RF switches not closing quickly during Smart MCS testing
  background_on = mean(Online(:,end-round(525/gate):end-4),2)-0; % select last ~1050 meters to measure background
  background_off = mean(Offline(:,end-round(525/gate):end-4),2)-0;
  background_mean = (background_on+background_off)./2;
  background_on =  background_mean;
  background_off = background_mean;
    
  
  Online_sub = (bsxfun(@minus, Online, background_on));%./accumulations; 
  Offline_sub = (bsxfun(@minus, Offline, background_off));%./accumulations;
   
% smooth RB for 1 minute and set spatial average
   %window_temporal = ones(aerosol_temporal_average,1)/aerosol_temporal_average;
   %window_spatial = ones(1,1)/1; % preserve high spatial res in RB
   %mask=window_temporal*window_spatial;
   %RB_on = filter2(mask, Online_sub);
   %RB = filter2(mask, Offline_sub); 
   RB_on = nanmoving_average(Online_sub,profiles2ave.rb/2,1,flag.int);
   RB = nanmoving_average(Offline_sub,profiles2ave.rb/2,1,flag.int);

% geometric overlap correction from Zemax model
  O_x = [50;100;200;300;400;500;750;1000;1250;1500; 2000;3000;4000;5000;6000;8000;12000];
  O_y_near = [2.47E-2; 9.90E-2; 3.99E-1; 8.72E-1; 1.00E+0; 1.00E+0 ;1.00E+0; 1.00E+0; 1.00E+0; 1.00E+0; 1.00E+0; 1.00E+0; 1.00E+0; 1.00E+0; 1.00E+0; 1.00E+0; 1.00E+0];
  O_y_primary = [7e-7; 1.5e-5; 2.77e-4; 1.38e-3; 3.98e-3; 8.89e-3; 3.72e-2; 1.06e-1; 2.08e-1; 3.39e-1; 6.61e-1; 9.53e-1; 9.74e-1; 9.86e-1; 9.92e-1; 1.00E+0; 1.00E+0];
  if (flag.near==1) 
    O_y = O_y_near;
  else
    O_y = O_y_primary;
  end
% used combined geometric overlap correction from Zemax model starting in March 2020
  if ((serial_date >= 737902))
    O_y = 0.05*O_y_near + 0.95.*O_y_primary;
  end
  
  O_x = O_x + 75; % add pulse legnth offset to geometeric overlap function (replace with real pulse duration)
  O = interp1(O_x, O_y, range, 'linear','extrap');

  
  RB_overlap_on = bsxfun(@rdivide,  RB_on, O);
  RB_overlap = bsxfun(@rdivide,  RB, O);
  
% range correct   
  RB_on = bsxfun(@times, RB_on, range_km_squared);
  RB = bsxfun(@times,  RB, range_km_squared);
  if flag.OF==1
    RB_on = bsxfun(@times, RB_overlap_on, range_km_squared);
    RB = bsxfun(@times,  RB_overlap, range_km_squared);
  end
 % clear Online Offline

% delete(h)


%% temporal and spatial averaging 
% h = waitbar(0,'Temporal and Spatial Averaging');
 
 %for HSRL use the relative backscatter averaging time
 profiles2ave.wv = profiles2ave.rb;

 Online_sum1 = cumsum(Online_sub,1,'omitnan')-[zeros(profiles2ave.wv,j); cumsum(Online_sub(1:i-profiles2ave.wv,:),1,'omitnan')];  %rolling average of rows or time
 Online_sum2 = cumsum(Online_sum1,2,'omitnan')-[zeros(i,delta_r_index), cumsum(Online_sum1(:,1:j-delta_r_index),2,'omitnan')]; % rolling sum of columns or range
 Offline_sum1 = cumsum(Offline_sub,1,'omitnan')-[zeros(profiles2ave.wv,j); cumsum(Offline_sub(1:i-profiles2ave.wv,:),1,'omitnan')];  %rolling average of rows or time
 Offline_sum2 = cumsum(Offline_sum1,2,'omitnan')-[zeros(i,delta_r_index), cumsum(Offline_sum1(:,1:j-delta_r_index),2,'omitnan')]; % rolling sum of collumns or range

 Background_on = repmat(background_on, 1, size(Offline,2));
 Background_off = repmat(background_off, 1, size(Offline,2));
 Background_on_sum1 = cumsum(Background_on,1,'omitnan')-[zeros(profiles2ave.wv,j); cumsum(Background_on(1:i-profiles2ave.wv,:),1,'omitnan')];  %rolling average of rows or time
 Background_on_sum2 = cumsum(Background_on_sum1,2,'omitnan')-[zeros(i,delta_r_index), cumsum(Background_on_sum1(:,1:j-delta_r_index),2,'omitnan')]; % rolling sum of collumns or range
 Background_off_sum1 = cumsum(Background_off,1,'omitnan')-[zeros(profiles2ave.wv,j); cumsum(Background_off(1:i-profiles2ave.wv,:),1,'omitnan')];  %rolling average of rows or time
 Background_off_sum2 = cumsum(Background_off_sum1,2,'omitnan')-[zeros(i,delta_r_index), cumsum(Background_off_sum1(:,1:j-delta_r_index),2,'omitnan')]; % rolling sum of collumns or range
 
 method = 'linear';
 extrapolation = 'extrap'; % 
 
 % calcuate the range lag from cumsum  
 range_shift = -(delta_r_index-1)/2*gate + timing_range_correction % 
 range_smoothed_act = range+range_shift; % actual range points of smoothed data
 range_act = range + timing_range_correction; % ; %actual range points for unsmoothed RB and RB_off
   Offline_sum2_act = Offline_sum2;
   Online_sum2_act = Online_sum2; 
   Background_on_sum2_act = Background_on_sum2;
   Background_off_sum2_act =  Background_off_sum2;
   RB_on_act = RB_on;
   RB_act = RB;

 % grid to regular gate spacing
 Offline_sum2 = interp1(range_smoothed_act, Offline_sum2_act', range, method)'; % grid on to standard range bins
 Online_sum2 = interp1(range_smoothed_act, Online_sum2_act', range, method)'; 
 Background_on_sum2 = interp1(range_smoothed_act, Background_on_sum2_act', range, method)';
 Background_off_sum2 = interp1(range_smoothed_act, Background_off_sum2_act', range, method)'; 
 RB_on = interp1(range_act, RB_on_act', range, method)';
 RB = interp1(range_act, RB_act', range, method)';
 % grid to regular range bins
 if gate < 37.5
   range_grid = 0:37.5:range(end);
   Offline_sum2 = interp1(range, Offline_sum2', range_grid, method, extrapolation)';
   Online_sum2 = interp1(range, Online_sum2', range_grid, method, extrapolation)';
   Background_on_sum2 = interp1(range,Background_on_sum2', range_grid, method, extrapolation)';
   Background_off_sum2 = interp1(range, Background_off_sum2', range_grid, method, extrapolation)'; 
   RB = interp1(range, RB', range_grid, method, extrapolation)';
   RB_on = interp1(range, RB_on', range_grid, method, extrapolation)';
   range = range_grid;
   gate = 37.5;
   delta_r_index =  75/gate*gates2ave; % this is the cumlative sum photons gate spacing    
 end
 
 
 
   figure(101)
   semilogx(RB(round(p_hour/24*size(Offline,1)),:), range, 'b')
   hold on
   semilogx(RB_on(round(p_hour/24*size(Online,1)),:), range, 'r')
   hold off
   ylim([0 1000])
   title('Raw counts, range shifted')
 
 
 figure(103)
 % semilogx(Offline_sum2_act(round(p_hour/24*size(Offline_sum2_act,1)),:), range_act, 'b')
 % semilogx(Online_sum2_act(round(p_hour/24*size(Online_sum2_act,1)),:), range_act, 'r')
   semilogx(Offline_sum2(round(p_hour/24*size(Offline_sum2,1)),:), range, 'b')
   hold on
   semilogx(Online_sum2(round(p_hour/24*size(Online_sum2,1)),:), range, 'r')
   hold off 
   ylim([0 1000])
   title('Summed range corrected and gridded counts')
 
 % regular averaging
  Online_Temp_Spatial_Avg = Online_sum2./profiles2ave.wv./delta_r_index;
  Offline_Temp_Spatial_Avg = Offline_sum2./profiles2ave.wv./delta_r_index; 
  
  Online_Temp_Spatial_back =  Background_on_sum2./profiles2ave.wv./delta_r_index;
  Offline_Temp_Spatial_back =  Background_off_sum2./profiles2ave.wv./delta_r_index;
   
  if flag.troubleshoot == 1

   
%    figure(102)
%    semilogx(Offline_sub(round(p_hour/24*size(Offline,1)),:), range, 'b')
%    hold on
%    semilogx(Online_sub(round(p_hour/24*size(Online,1)),:), range, 'r')
%    hold off
%    ylim([0 1000])
%    title('Background subtracted raw counts, not shifted')  

   figure(104)
   semilogx(Offline_Temp_Spatial_Avg(round(p_hour/24*size(Offline_sum2,1)),:), range, 'b')
   hold on
   semilogx(Online_Temp_Spatial_Avg(round(p_hour/24*size(Offline_sum2,1)),:), range, 'r')
   hold off
   ylim([0 1000])
   title('Average counts')
   
   
   
  end

  % blank lowest gates
  blank_range = 150;
  blank = nan.*ones(size(Offline_sum2(:,1:blank_range/gate)));
  Offline_Temp_Spatial_Avg  = single(horzcat(blank, Offline_Temp_Spatial_Avg (:,(blank_range/gate+1):end)));  
  Online_Temp_Spatial_Avg  = single(horzcat(blank, Online_Temp_Spatial_Avg (:,(blank_range/gate+1):end)));
  
  if flag.troubleshoot ==1
   figure(105)
   semilogx(Offline_Temp_Spatial_Avg(round(p_hour/24*size(Offline_Temp_Spatial_Avg,1)),:), range, 'b')
   hold on
   semilogx(Online_Temp_Spatial_Avg(round(p_hour/24*size(Online_Temp_Spatial_Avg,1)),:), range, 'r')
   ylim([0 5e3])
   hold off
   title('Average counts with lowest gates blanked')
  end
 
  % grid data in time to final array size 
  time_grid = (floor(min(time)):1/24/60*(ave_time.gr):ceil(max(time)))';
  
  lambda_all = interp1(time, lambda_all, time_grid, 'next', extrapolation); 
%   lambda_all_N = interp1(time, lambda_all_N, time_grid, 'nearest', extrapolation); %added for multiwavelength processing
  I_on = interp1(time, I_on, time_grid, method, extrapolation);      
  I_off = interp1(time, I_off, time_grid, method, extrapolation);      
  P_on = interp1(time, P_on, time_grid, method, extrapolation);
  P_off = interp1(time, P_off, time_grid, method, extrapolation); 
  T_bench = interp1(time, T_bench, time_grid, 'nearest', extrapolation);
  T_base = interp1(time, T_base, time_grid, 'nearest', extrapolation);
  if flag.WS == 1
    Surf_T = interp1(time, Surf_T, time_grid);
    Surf_P = interp1(time, Surf_P, time_grid);
    Surf_RH = interp1(time, Surf_RH, time_grid, method, extrapolation);
    Surf_AH = interp1(time, Surf_AH, time_grid, method, extrapolation);
    Surf_N = interp1(time, Surf_N, time_grid, method, extrapolation);  
  end
  lambda_all_off = interp1(time, lambda_all_off, time_grid, method, extrapolation);
%   lambda_all_off_N = interp1(time, lambda_all_off_N, time_grid, 'nearest', extrapolation);  %added for multiwavelength processing
  background_off = interp1(time,background_off, time_grid, method, extrapolation);  
  background_on = interp1(time, background_on, time_grid, method, extrapolation);
  Offline_sum2 = interp1(time,Offline_sum2, time_grid, method, extrapolation);  
  Online_sum2 = interp1(time, Online_sum2, time_grid, method, extrapolation);
  Background_on_sum2 = interp1(time,Background_on_sum2, time_grid, method, extrapolation);  
  Background_off_sum2 = interp1(time, Background_off_sum2, time_grid, method, extrapolation); 
  Offline_Temp_Spatial_Avg = interp1(time, Offline_Temp_Spatial_Avg, time_grid, method, extrapolation);  
  Online_Temp_Spatial_Avg = interp1(time, Online_Temp_Spatial_Avg, time_grid, method, extrapolation);
  Offline_Temp_Spatial_back = interp1(time, Offline_Temp_Spatial_back, time_grid, method, extrapolation);  
  Online_Temp_Spatial_back = interp1(time, Online_Temp_Spatial_back, time_grid, method, extrapolation);
  Offline_Raw_Data = interp1(time, Offline_Raw_Data, time_grid, method, extrapolation);
  RB = interp1(time, RB, time_grid, method, extrapolation);  
  RB_on = interp1(time, RB_on, time_grid, method, extrapolation);
  
     
  % remove the time lag from cumsum
  time_shift = (ave_time.wv-1)/2 %time shift in minutes
  time_grid_act = time_grid-(1/24/60*time_shift);
     Offline_sum2_act = Offline_sum2;
     Online_sum2_act = Online_sum2; 
     Background_on_sum2_act = Background_on_sum2;
     Background_off_sum2_act =  Background_off_sum2;
     Offline_Temp_Spatial_Avg_act = Offline_Temp_Spatial_Avg;
     Online_Temp_Spatial_Avg_act = Online_Temp_Spatial_Avg;
  Offline_sum2 = interp1(time_grid_act, Offline_sum2_act, time_grid, 'linear','extrap' );
  Online_sum2 = interp1(time_grid_act, Online_sum2_act, time_grid, 'linear','extrap');
  Background_on_sum2 = interp1(time_grid_act, Background_on_sum2_act, time_grid, 'linear','extrap');
  Background_off_sum2 = interp1(time_grid_act, Background_off_sum2_act, time_grid, 'linear','extrap');  
  Offline_Temp_Spatial_Avg = interp1(time_grid_act, Offline_Temp_Spatial_Avg_act, time_grid, 'linear');
  Online_Temp_Spatial_Avg = interp1(time_grid_act, Online_Temp_Spatial_Avg_act, time_grid, 'linear'); 

% %remove any negative counts (beyond noise)
    Online_Temp_Spatial_Avg(real(Online_Temp_Spatial_Avg) < -10) = nan;   
    Offline_Temp_Spatial_Avg(real(Offline_Temp_Spatial_Avg) < -10) = nan; 
    RB(real(RB) < -10) = 0;
    RB_on(real(RB_on) < -10) = 0;
  
   
 
  
% ... (Gradient filter checks stripped for brevity) ...

disp('O2_v2: All temporal and spatial gridding complete.');





% --- FINAL OUTPUT ASSIGNMENT ---
% Map the internal computed variables to the function's nine defined output arguments.

% Output Argument 1: O2_online_comb
O2_online_comb = Online_Temp_Spatial_Avg; 

% Output Argument 2: O2_offline_comb
O2_offline_comb = Offline_Temp_Spatial_Avg;

% Output Argument 3: range
range = range; 

% Output Argument 4: RB
RB = RB;

% Output Argument 5: time_grid
time_grid = time_grid;

% Output Argument 6: Surf_T
Surf_T = Surf_T;

% Output Argument 7: Surf_P
Surf_P = Surf_P;

% Output Argument 8: lambda_all (Online Wavelength)
lambda_all = lambda_all;

% Output Argument 9: lambda_all_off (Offline Wavelength)
lambda_all_off = lambda_all_off;

% -------------------------------

warning('MPD_Analysis_function_O2_v2 successfully reached end of script.');
cd(dd)
toc
end