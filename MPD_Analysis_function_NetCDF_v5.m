function[] = MPD_Analysis_function_NetCDF_v5(data_on, data_off, folder_in, date_in, MCS, write_data_folder, flag, node, wavemeter_offset,...
    profiles2ave, P0, switch_ratio, ave_time, timing_range_correction, blank_range, p_hour, catalog, Afterpulse_File)

%% notes
%Amin Nehrir (original author)
        % last modification....8/25/2011
%addtional modifications by Spuler in 2012
%additional modification by Spuler in Oct-Dec 2013 
        % added ability to read binary data 
        % added ability to read near and far range channels
        % vectorized filtering and added bsnfun asf needed to speed up
        % added matrix decimation to speed up
%additional modification by Spuler in Apr 2014
        % added data masking based on signal/background ratio
        % removed vestigal code
        % rewrote time vector reconstruction to speed up
        % as of April 2014, processing a 24 hour day takes about 3 min
%additional modification by Spuler in Jun 2014     
        % modifed code to work with 75m range bins (500 ns accumulations)
        % vs 7.5m range bins with 50 ns accumulations 
        % calculate offline SNR - not just a S/B ratio 
        % as of June 2014, processing a 24 hour day takes about 36 sec
%additional modification by Spuler in Jul-Sept 2014     
        % corrected HITRAN T_ref=298K and P_ref=1atm (incorrectly using surface T and P)
        % read in the 2008 HITRAN pressure shifts (was using a single pressure shift hardcoded for the 828 line) 
        % corrected Voigt width 'y' term to use all lines (bug in code applied first line width to all lines)
        % as of Aug 2014, processing a 24 hour day takes about 30-60 sec  (10min vs 1min averaging)       
%additional modification by Spuler in Dec 2014 to Feb 2015     
        % added error driven thresholds 
        % switched to 7 kHz processing to get better background subtraction
        % processing a 24 hour day takes about 24-60 sec  (10min vs 1min averaging) 
        % Note about 7 column header (added 23-Feb-2015)
        % header has been : Julian day, wavelength, diode current, power transmitted, T(C), p(mbar), humidity
        % need to, wire diode current (has not been connected) and add: bin duration, # bins, # accumulations 
%additional modification by Spuler in Mar-Apr 2015     
        % removed routine to add time gaps for missing data and replaced with 30 sec interpolated gridding 
        % repalced imagesc plot calls with pcolor
%additional modification by Spuler in Jan 2016     
        % removed filter2 calls as these shift data backward 1/2 grid point when called with even number of grid averages
        % replaced with nanmoving_ave which has no lag or shift
        % corrected for timing of MCS to laser pulse center 
        % corrected the lag induced by the cumsum routine with a linear interpolation
        % masked out the lowest 4 gates, after summing as these are contaminated by the outgoing laser pulse
        % grid data right after summing/averaging (no more decimation step)
        % compute DIAL equation (N) and error equation (N_error) at native grid spacing (75m range)
        % still need to correct for cumsum rotine timelag  
        % Note about 7 column header (added 28-Mar-2016)
        % Online header: (1) Julian day, (2) online wavelength, (3) online diode current, (4) ave power transmitted,...
        % (5) MCS bin duration, (6) MCS # bins, (7) MCS # accumulations 
        % Offline header: (1) Julian day, (2) offline wavelength, (3) offline diode current, (4) bench temp, ...
        % (5) surface temp, C, (6) surface pressure, hPa, (7) surface relative humidity 
% additional modification by Spuler in Apr 2016
        % corrected for cumsum routing timelag
        % correct for 1/2 bin range lag within number density calculation
        % added toggle on/off extrapolation within nanmoving_average
        % added some troubleshooting graphs and ability to toggle on/off
        % added the ability to decimate all data to half wv averaged resolution
        % added housekeeping plots and save data (laser aver power, seed current, WS data)
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

delta_r_index =  75/gate; % this is the cumlative sum photons gate spacing 
%delta_r_index =  1; % proces at the native gate spacing
delta_r = delta_r_index*gate*100; % delta r in cm
r1 = round(1500/gate); % index for smoothing range 1 (1500m)
r2 = round(2500/gate); % index for smoothing range 2 (2500m)
spatial_average1 = 150/gate; %150 meter smoothing range 1 
spatial_average2 = 300/gate; %300 meter smoothing between range 1 and 2
spatial_average3 = 600/gate; %600 meter smoothing above range 2
% this filter creates range delays 

%lambda offset for testing purposes only -- set to zero
%lambda_offset = 0;  

%% Importing online and offline files from the selected date

%read_time_in = 2; % read data in time increments of seconds (default it 2sec) 
%if flag.near == 1
%  [dummy,dummy,near_on, near_off, folder_in] = MPD_File_Retrieval_NetCDF_v4(flag, MCS.bins, folder_in, read_time_in); %use to read binary data (bin number passed in) 
%   % this is a way to process the near range data instead of the normal
%   data_on = near_on;
%   data_off = near_off;
%   if strcmp(node,'MPD4') == 1 % MPD4 is using a low range channel
%     blank_range = 187.5; % low range 
%   end
%else
%  % [data_on,data_off,dummy, dummy,folder_in] = File_Retrieval_NetCDF_v4(flag, MCS.bins, folder_in); %use to read binary data (bin number passed in) 
%   [data_on,data_off,folder_in] = MPD_File_Retrieval_NetCDF_v3(flag, MCS.bins, folder_in, read_time_in); %use to read binary data (bin number passed in) 
%  %[data_on,data_off,folder_in] = File_Retrieval_NetCDF_v2_noPow(MCS.bins, folder_in); %use to read binary data (bin number passed in) 
%end

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
  %time2 = time2-0.25; % this was a fix for computer switch 20-Oct to 26-Oct 2017
  %time2(time2>292.5)= time2(time2>292.5)-0.25; % 19-Oct-2017 fix
  time = time2;
  time = time(~(time2<(nanmedian(time2)-1.5)));
  time= time(~(time>(nanmedian(time)+1.5)));
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

% analyze system stability 
time_step = diff(Online_Raw_Data(:,1)); % median time step in seconds
%time_step(time_step >= 4*median(time_step))= nan; % remove the outliers
%if flag.troubleshoot == 1
  figure(11)
  edges= 0:1:100;
  %histogram(time_step*3600*24,100, 'Normalization', 'probability')
  histogram(time_step*3600*24,edges, 'Normalization', 'probability')
  axis([0 100 0 1])
  xlabel('data aquistion time, seconds')
  ylabel('probability')
  set(gca, 'YScale', 'log')
%end

%Calculating average wavelenth
lambda_all=Online_Raw_Data(:,2);
lambda_diff= Online_Raw_Data(:,3);
lambda_all_off=Offline_Raw_Data(:,2);
lambda_diff_off= Offline_Raw_Data(:,3);

ii=1;
for ii=1:size(lambda_all,1)
    if ii>1&& lambda_all(ii)<828
        lambda_all(ii)=lambda_all(ii-1);
        lambda_diff(ii)=lambda_diff(ii-1);
    end
    if ii>1&& lambda_all_off(ii)<828
        lambda_all_off(ii)=lambda_all_off(ii-1);
        lambda_diff_off(ii)=lambda_diff_off(ii-1);
    end
end

lambda_on = median(lambda_all);
lambda_off = median(lambda_all_off);
 
% if nanstd(lambda_all) >= 5e-4    
%    h = msgbox('Online wavelength not stable during time period', 'Warning','warn');
% end
    % check for multiple wavelengths
    edges_on=828.1800:.00001:828.2200;
    mult_fac = 1;
    lambda_on_set = lambda_all-lambda_diff;
    [value,edges]=histcounts(round(mult_fac*lambda_on_set,3)/mult_fac,edges_on); % bin rounded wavelengths
    lambda_on_N = edges(value>=1000)  % wavelength values with occurance > 10
    %lambda_F = value(value~=0);  % frequency of occurance
    lambda_all_N = round(mult_fac*lambda_on_set,3)/mult_fac; 
    figure(1234)
    plot(lambda_all_N)
    hold on
    plot(lambda_all)
    hold off

    edges_off=828.280:.00005:828.320;
    mult_fac = 1;
    lambda_off_set = lambda_all_off-lambda_diff_off;
    [value,edges]=histcounts(round(mult_fac*lambda_off_set,3)/mult_fac,edges_off); % bin rounded wavelengths
    lambda_off_N = edges(value>=1000)  % wavelength values with occurance > 10
    % select the most common offline values associated with the online 
      value_sort = [value; edges(1:end-1)]';
      values_sorted = sortrows(value_sort, 1);
      lambda_off_N = sort(values_sorted(end-size(lambda_on_N,2)+1:end))
    %lambda_off_F = value(value~=0);  % frequency of occurance
    lambda_all_off_N=round(mult_fac*lambda_off_set,3)/mult_fac;
    figure(5678)
    plot(lambda_all_off_N)
    hold on
    plot(lambda_all_off)
    hold off

  if isempty(lambda_off_N) == 1
      lambda_off_N = lambda_off 
  end
    
 folder_date = textscan(folder_in(end-7:end), '%8f', 1); folder_date=folder_date{1}; % change to read ingore the NF and FF 
 %folder_CH = textscan(folder_in(end-1:end), '%s'); folder_CH=folder_CH{1}; % change to read ingore the NF and FF 
 %folder_CH=folder_CH{1:1};
  
% Parsing off the auxiliary 
 Online = single(Online_Raw_Data(:,10:end)); 
 Offline = single(Offline_Raw_Data(:,10:end));
 
% vector in meters
range = single(0:gate:(size(Online,2)-1)*gate);
%range = range+timing_range_correction;  %make range correction from trigger timing
%time = (Online_Raw_Data(:,1)); 
%time = time2;
%clear Online_Raw_Data  Offline_Raw_Data 
 
  if flag.afterpulse == 1   % afterpulse correction
   
%  read my vesrion of the afterpulse calibration files      
%      afterpulse_filename =  sscanf(Afterpulse_File, '%c', 25);   
%      if flag.near == 1
%        afterpulse_filename = strcat(afterpulse_filename, '_near');  
%      end
%      load (afterpulse_filename, 'ap_spline_sub_off', 'ap_spline_sub_on')  
      
    % read the afterpulse nc file identified in the json file 
    if strcmp(getenv('HOSTNAME'),'fog.eol.ucar.edu')
     serv_path = '/home/rsfdata/Processing/'; % when running on server
    elseif strcmp(getenv('HOSTNAME'),'')
      serv_path = '../'; % running locally 
    else
     serv_path = '/Volumes/eol/fog1/rsfdata/MPD/calibration/'; % 
    end
    ap_filename = strcat(serv_path, 'eol-lidar-calvals/calfiles/', Afterpulse_File)   
   
    ncid = netcdf.open(ap_filename, 'NC_NOWRITE');
    %ncdisp(ap_filename, '/', 'min') % use this to display all variables
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

    range_shift = -(delta_r_index-1)/2*gate + timing_range_correction; % 
    range_act = range + range_shift; % %actual range points 

%      figure(1004)
%      semilogy(ap_range, afterpulse_off, 'bo-')
%      hold on
%      semilogy(ap_range, afterpulse_on, 'b+-')    
% %     semilogy(ap_range, afterpulse_off, 'ro-')
% %     semilogy(ap_range, afterpulse_on, 'r+-')  
%      semilogy(range_act, ap_spline_sub_off, 'ro-')
%      semilogy(range_act, ap_spline_sub_on, 'r+-')
%      hold off
%      legend('hayman_{off}', 'hayman_{on}', 'spuler_{off}', 'spuler_{on}') 
% %    legend('high-gain_{off}', 'high-gain_{on}', 'low-gain_{off}', 'low-gain_{on}') 
%     ylim([1e-2 1e5])
%     xlim([-200 4000])
%     ylabel('counts')
%     xlabel('range (m)')
%     grid on
%     grid minor 
   
   %grid to the current range and substitude nc file for mat file
   ap_spline_sub_off = spline(ap_range, afterpulse_off, range_act);
   ap_spline_sub_on = spline(ap_range, afterpulse_on, range_act);
   
   Offline_ap_sub = (bsxfun(@minus, Offline, ap_spline_sub_off));
   Online_ap_sub = (bsxfun(@minus, Online, ap_spline_sub_on)); 

%     figure(1002)
%     semilogy(range, Online(100,:), 'b')
%     hold on
%     semilogy(range, Online_ap_sub(100,:), 'r')
%     hold off
    
    Online = Online_ap_sub;
    Offline = Offline_ap_sub;
  end

if flag.pileup == 1
% apply linear correction factor to raw counts 
  %t_d=32E-9; % module dead time of Perkin Elmber
  t_d=37.25E-9; %Excelitas SPCM-AQRH-13 Module 24696
  %t_d=50E-9; %Emperical best fit to remove the noise in the WV behind clouds
  %t_d=34E-9; %Excelitas SPCM-AQRH-13 Module 24696 for count rates < 5 Mc/s
  % MCSC gives counts accumulated for set bin duration so convert to count rate  C/s.
  % divide by bin time in sec (e.g., 500ns) and # of acumulations (e.g., 10000)
  % e.g., 10 accumlated counts is 2000 C/s
  C_Online = 1./(1-(t_d.*(Online./(MCS.bin_duration*1E-9*MCS.accum))));   
  C_Offline = 1./(1-(t_d.*(Offline./(MCS.bin_duration*1E-9*MCS.accum))));   
  Online = Online.*C_Online;
  Offline = Offline.*C_Offline;
end

% clear Online_Raw_Data Offline_Raw_Data data_on data_off

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
  %background_mean = (background_on+background_off)./2;
  
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
  if (flag.near==1) && (strcmp(node,'MPD4') == 1)
    O_y = [2.47E-2; 9.90E-2; 3.99E-1; 8.72E-1; 1.00E+0; 1.00E+0 ;1.00E+0; 1.00E+0; 1.00E+0; 1.00E+0; 1.00E+0; 1.00E+0; 1.00E+0; 1.00E+0; 1.00E+0; 1.00E+0; 1.00E+0];
  else
    O_y = [7e-7; 1.5e-5; 2.77e-4; 1.38e-3; 3.98e-3; 8.89e-3; 3.72e-2; 1.06e-1; 2.08e-1; 3.39e-1; 6.61e-1; 9.53e-1; 9.74e-1; 9.86e-1; 9.92e-1; 1.00E+0; 1.00E+0];
  end
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
   
  if flag.troubleshoot == 1

   
   figure(102)
   semilogx(Offline_sub(round(p_hour/24*size(Offline,1)),:), range, 'b')
   hold on
   semilogx(Online_sub(round(p_hour/24*size(Online,1)),:), range, 'r')
   hold off
   ylim([0 1000])
   title('Background subtracted raw counts, not shifted')  

   figure(104)
   semilogx(Offline_Temp_Spatial_Avg(round(p_hour/24*size(Offline_sum2,1)),:), range, 'b')
   hold on
   semilogx(Online_Temp_Spatial_Avg(round(p_hour/24*size(Offline_sum2,1)),:), range, 'r')
   hold off
   ylim([0 1000])
   title('Average counts')
  end

  % blank lowest gates
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
  lambda_all_N = interp1(time, lambda_all_N, time_grid, 'nearest', extrapolation); %added for multiwavelength processing
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
  lambda_all_off_N = interp1(time, lambda_all_off_N, time_grid, 'nearest', extrapolation);  %added for multiwavelength processing
  background_off = interp1(time,background_off, time_grid, method, extrapolation);  
  background_on = interp1(time, background_on, time_grid, method, extrapolation);
  Offline_sum2 = interp1(time,Offline_sum2, time_grid, method, extrapolation);  
  Online_sum2 = interp1(time, Online_sum2, time_grid, method, extrapolation);
  Background_on_sum2 = interp1(time,Background_on_sum2, time_grid, method, extrapolation);  
  Background_off_sum2 = interp1(time, Background_off_sum2, time_grid, method, extrapolation); 
  Offline_Temp_Spatial_Avg = interp1(time, Offline_Temp_Spatial_Avg, time_grid, method, extrapolation);  
  Online_Temp_Spatial_Avg = interp1(time, Online_Temp_Spatial_Avg, time_grid, method, extrapolation);
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
    Online_Temp_Spatial_Avg(real(Online_Temp_Spatial_Avg) < -10) = 0;   
    Offline_Temp_Spatial_Avg(real(Offline_Temp_Spatial_Avg) < -10) = 0;  
    RB(real(RB) < -10) = 0;
    RB_on(real(RB_on) < -10) = 0;
  
% clear Online_Raw_Data Online Offline Online_sub Offline_sub data_on data_off C_Online C_Offline ... 
%      Offline_sum1 Online_sum1 Background_on_sum1 Background_off_sum1 Background_off Background_on 
    
  if flag.gradient_filter == 1
    [FX,FY] = gradient(Offline_Temp_Spatial_Avg);
    Offline_Temp_Spatial_Avg(FX<-(2000*(MCS.bin_duration*MCS.accum)/500/14200)) = nan; % remove falling edge of clouds
    Offline_Temp_Spatial_Avg(FX> (2000*(MCS.bin_duration*MCS.accum)/500/14200)) = nan; % remove leading edge of clouds   
    
    if flag.troubleshoot == 1
      imagesc(time_grid,range./1e3, FX');
      axis xy; colorbar('EastOutside'); caxis([-1000 1000]);   
      title({[' Vertical Gradient']},...
      'fontweight','b','fontsize',16)
      ylabel('Altitude (km)','fontweight','b','fontsize',20); 
      datetick('x','HH','keeplimits');
      axis([fix(min(time_grid)) fix(min(time_grid))+1 0 7])
      colormap(C)
    end
  end



 % delete(h)

 
%% Spectral Line Fitting
% h = waitbar(0,'Line Fiting');

%Voigt profile calculation
WNmin = 1/(828+4)*1e7;  % 832.0 nm converted to wavenumber
WNmax = 1/(828-4)*1e7;  % 824.0 nm converted to wavenumber

%Find lines from WNmin to WNmax to calculate voigt profile
line_indices = Hitran.file(1:size(Hitran.file,1),1)>WNmin & Hitran.file(1:size(Hitran.file,1),1)<WNmax;

line = double(Hitran.file(line_indices, 1:size(Hitran.file,2)));

%Calculate temperature and pressure profile
if flag.WS == 1
    T0 = nanmedian(Surf_T)+273.15
    P0 = nanmedian(Surf_P)
   % T0 = median(Surf_T)+273.15
   % P0 = median(Surf_P)
else
  T0 = 273+30; % surface temperature
end

  T = T0-0.0065.*range; % set to match the sounding

% pressure in atmospheres
  %P0  = 0.83; % surface pressure in Boulder
  %P0  = 0.929; % surface pressure in Ellis KS
%  P = P0.*(T0./T).^-5.5;   % set this value to match sounding
  P = P0.*(T0./T).^-5.2199;   % set this value to match sounding
  %p_s = 7.75e-3.*exp(-range./3000).*P; % estimate of the partial pressure of wv for self broadening

Hitran.T00 = 296;              % HITRAN reference temperature [K]
Hitran.P00 = 1;                % HITRAN reference pressure [atm]
Hitran.nu0_0 = line(:,1);      % absorption line center wavenumber from HITRAN [cm^-1]
Hitran.S0 = line(:,2);         % initial linestrength from HITRAN [cm^-1/(mol*cm^-2)]   
Hitran.gammal0 = line(:,4);    % air-broadened halfwidth at T_ref and P_ref from HITRAN [cm^-1/atm]
Hitran.gamma_s = line(:,5);    % self-broadened halfwidth at T_ref and P_ref from HITRAN [cm^-1/atm]
Hitran.E = line(:,6);          % ground state transition energy from HITRAN [cm^-1]  
Hitran.alpha = line(:,7);      % linewidth temperature dependence factor from HITRAN
Hitran.delta = line(:,8);      % pressure shift from HiTRAN [cm^-1 atm^-1]

% new code to handle multiple wavelength changes during a single day
for l=1:length(lambda_on_N)

%offset = 0.00017;  % This is the absolute accuracy (nm) of the Bristol 671A at 828 nm  
offset = -wavemeter_offset*1e9; % This is the measured correction (nm) for each unit 

Hitran.nu_on = 1/(lambda_on_N(l)+offset)*1e7;
Hitran.nu_off = 1/(lambda_off_N(l)+offset)*1e7;
%Hitran.nu_on = 1/lambda_on_N(l)*1e7;
%Hitran.nu_off = 1/lambda_off_N(l)*1e7;

%sigma_on_total = zeros(size(Online_Temp_Spatial_Avg,1),size(Online_Temp_Spatial_Avg,2));


for i = 1:size(Online_Temp_Spatial_Avg,2); % calculate the absorption cross section at each range
    j=rem(i,25);
    if j==0 
 %     waitbar(i/(size(Online_Temp_Spatial_Avg,2)),h);
    end
  
    %nu0 = nu0_0-.0428.*P(i);    %Pressure shift
    %nu0 = nu0_0+delta.*P(i);     %calculate the pressure shifts for selected lines as function of range
    Hitran.nu0 = Hitran.nu0_0+Hitran.delta.*(P(i)./Hitran.P00); % unclear if it should be Pi/P00
    %gammal = gammal0.*P(i).*((T00./T(i)).^alpha);    %Calculate Lorentz lineweidth at P(i) and T(i)
    Hitran.gammal = Hitran.gammal0.*(P(i)./Hitran.P00).*((Hitran.T00./T(i)).^Hitran.alpha);    %Calculate Lorentz lineweidth at P(i) and T(i)
    % revise pressure broadened halfwidth to include correction for self broadening term 
    %gammal = ((P(i)-p_s(i)).*gammal0 + p_s(i).*gamma_s).*((T00./T(i)).^alpha);    %Calculate Lorentz lineweidth at P(i) and T(i)
   
    const.m = 18.015E-3./6.022E23; % mass of a single water molecule
    const.k_B = 1.3806488e-23; % (J/K)
    const.c = 299792458; % (m/s) (exact)
   
    Hitran.gammad = (Hitran.nu0).*((2.0.*const.k_B.*T(i).*log(2.0))./(const.m.*const.c^2)).^(0.5);  %Calculate HWHM Doppler linewidth at T(i)
    
    % term 1 in the Voigt profile
    y = (Hitran.gammal./Hitran.gammad).*((log(2.0)).^(0.5));
    
    % term 2 in the Voigt profile
    x_on = ((Hitran.nu_on-Hitran.nu0)./Hitran.gammad).*(log(2.0)).^(0.5);
    x_off = ((Hitran.nu_off-Hitran.nu0)./Hitran.gammad).*(log(2.0)).^(0.5);
    
    %setting up Voigt convolution
    t = (-(size(line,1))/2:1:size(line,1)/2-1); %set up the integration spectral step size
    t = repmat(t',1,length(x_on))';
    x_on = repmat(x_on,1,length(t));
    x_off = repmat(x_off,1,length(t));
    y = repmat(y,1,length(t));
    f_on = (exp(-t.^2.0))./(y.^2.0+(x_on-t).^2.0);  % combined Voigt term 1 and 2
    f_off = (exp(-t.^2.0))./(y.^2.0+(x_off-t).^2.0);
    %Voigt integration over all of the lines at the on and offline locations
    z_on = trapz(t(1,:),f_on,2);
    z_off =  trapz(t(1,:),f_off,2);
    integral_on = z_on;
    integral_off = z_off;
    %Calculate linestrength at temperature T
    S = Hitran.S0.*((Hitran.T00./T(i)).^(1.5)).*exp(1.439.*Hitran.E.*((1./Hitran.T00)-(1./T(i))));
    %term3= (1-exp(-1.4388.*nu0_0./T(i)))./(1-exp(-1.4388.*nu0_0./T00));
    %S = S0.*((T00./T(i)).^(1.5)).*exp(1.4388.*E.*((1./T00)-(1./T(i)))).*term3;
      
    %Calculate Doppler line-center cross section [cm^2]
    % K = (1./gammad).*(((log(2.0))./pi).^(0.5));
 
    %Calculate the Voigt profile
    K_on = (y(:,1)./pi).*integral_on;
    K_off = (y(:,1)./pi).*integral_off;
    
    %far_wing_on = nu_on./nu0.*tanh((1.986e-25.*nu_on)./(2.*1.38E-23.*T(i)))./tanh((1.986e-25.*nu0)./(2.*1.38E-23.*T(i)));
    %far_wing_off =nu_off./nu0.*tanh((1.986e-25.*nu_off)./(2.*1.38E-23.*T(i)))./tanh((1.986e-25.*nu0)./(2.*1.38E-23.*T(i)));
    
    %Calculate the Voigt profile absorption cross section [cm^2]
    sigmav_on = S.*(1./Hitran.gammad).*(((log(2.0))./pi).^(0.5)).*K_on; %.*far_wing_on;
    sigmav_off = S.*(1./Hitran.gammad).*(((log(2.0))./pi).^(0.5)).*K_off; %.*far_wing_off;
    
    %Sum contributions from all of the surrounding lines
    sigma_on_total(i) = sum(sigmav_on);
    sigma_off_total(i) = sum(sigmav_off);
    
end

sigma_on_total = repmat(sigma_on_total,size(Online_Temp_Spatial_Avg,1),1);
sigma_off_total = repmat(sigma_off_total,size(Online_Temp_Spatial_Avg,1),1);

% new lines to handle multiple wavelengths during the day
S_on_N(:,:,l) = sigma_on_total;
S_off_N(:,:,l) = sigma_off_total;

clear sigma_on_total sigma_off_total sigmav_on sigmav_off

end
% delete(h)
 
% combine the multiwavelenth into a single cross section matrix
sigma_on_total=zeros(size(S_on_N,1),size(S_on_N,2));
sigma_off_total=zeros(size(S_on_N,1),size(S_on_N,2));

for l=1:length(lambda_on_N)
  for i = 1:size(Online_Temp_Spatial_Avg,1)
    if round(single(lambda_all_N(i)),3)==round(single(lambda_on_N(l)),3)
     sigma_on_total(i,:)= S_on_N(i,:,l); 
    end
    if round(single(lambda_all_off_N(i)),3)==round(single(lambda_off_N(l)),3)
     sigma_off_total(i,:)= S_off_N(i,:,l); 
    end   
  end
end

%% DIAL Equation to calculate Number Density and error

 Inside = (Online_Temp_Spatial_Avg.*(circshift(Offline_Temp_Spatial_Avg, [0, -1])))./...
     ((circshift(Online_Temp_Spatial_Avg, [0, -1])).*Offline_Temp_Spatial_Avg);
 del_cross = single(1./(2.*(sigma_on_total-sigma_off_total).*gate*100));
 N =  (del_cross.*log(Inside)); 
 N(N == inf) = nan; 
  
 % error calculation   
 N_error = (1/2./(sigma_on_total-sigma_off_total)./(gate*100)...    
    .*sqrt(...
    (Online_sum2+Background_on_sum2)./Online_sum2.^2 + ...
    (circshift(Online_sum2, [0, -1])+Background_on_sum2)./circshift(Online_sum2, [0, -1]).^2 + ...
    (Offline_sum2+Background_off_sum2)./Offline_sum2.^2 + ...
    (circshift(Offline_sum2, [0, -1])+Background_off_sum2)./circshift(Offline_sum2, [0, -1]).^2));

 N_1_error = real(N_error./sqrt(spatial_average1));
 N_2_error = real(N_error./sqrt(spatial_average2));
 N_3_error = real(N_error./sqrt(spatial_average3));
 %Combine different averaging windows for different altitudes; 
 N_error=[N_1_error(:,1:r1)';N_2_error(:,r1+1:r2)';N_3_error(:,r2+1:size(N,2))']';
 %smooth again at the smallest resolution to avoid boundaries
 N_error = nanmoving_average(N_error,spatial_average1/2,2,flag.int);
 
  
  %% Smoothing Number Density for different range zones

 if flag.gradient_filter == 1
    N(N>1E18) = nan; % use this to remove high water vapor errors
 end
  
 % insert the surface water vapor number density 
 %Surf_N = Surf_AH./1e6.*6.022E23./18.015;
 %N(:,1) = Surf_N;
 %N(:,blank_range/gate) = Surf_N;
  
  N_avg = N;
  N_1_avg = nanmoving_average(N_avg,spatial_average1/2,2,flag.int);
  N_2_avg = nanmoving_average(N_avg,spatial_average2/2,2,flag.int);
  N_3_avg = nanmoving_average(N_avg,spatial_average3/2,2,flag.int);
  
  %Combine different averaging windows for different altitudes; 
  N_avg=[...
      N_1_avg(:,1:r1)';...
      N_2_avg(:,r1+1:r2)';...
      N_3_avg(:,r2+1:size(N,2))'    ]';
  clear N_1_avg N_2_avg N_3_avg dummy
  %smooth again at the smallest resolution to avoid boundaries
  N_avg = nanmoving_average(N_avg,spatial_average1/2,2,flag.int);


  
%% Mask the Number density data based on the error, correct for range center, and add WS data at lowest gate 

 N_masked = N_avg;
 N_masked(N_avg < 0) = nan; % remove non-pysical (negative) wv regions
 N_masked(abs(N_error./N_avg) > 3.00) = nan; % remove high error regions
% check_high_counts = Offline_Raw_Data(:,9:end)./(MCS.bin_duration*1e-9*MCS.accum);
% N_masked(check_high_counts > 5E6) = nan; % remove raw counts above linear count threshold (5MC/s)

%if strcmp(folder_CH,'NF') == 1
 %  N_masked(abs(N_error./N_avg) > 1.00) = nan; % remove high error regions
 %  N_masked(Offline_Raw_Data(:,8:end)./(MCS.bin_duration*1e-9*MCS.accum) > 1.25E6) = nan; % remove raw counts above linear count threshold (5MC/s)
 %end
 % blank first 300 meters when pulse is going out
 %blank = nan.*ones(size(N_masked(:,1:300/gate)));
 %N_masked = single(horzcat(blank, N_masked(:,(300/gate+1):end)));  


 % calcuate the range lag for number density (to center in range bin)
  range_shift_N = gate/2 % 
  range_act_N = range+range_shift_N; % actual range points of data
  N_act = N;
  N_avg_act = N_avg;
  N_masked_act = N_masked;
  N_error_act = N_error;
  % grid to regular (75 m) gate spacing
  N = interp1(range_act_N, N_act', range, method, extrapolation)'; % grid on to standard range bins
  N_avg = interp1(range_act_N, N_avg_act', range, method, extrapolation)'; 
  N_masked = interp1(range_act_N, N_masked_act', range, method, extrapolation)';
  N_error = interp1(range_act_N, N_error_act', range, method, extrapolation)'; 

  if flag.WS == 1 % use the weather station to fill in the bottom gates 
    N_avg(:,1) = Surf_N;  % gate 1, 0 meter
    N_masked(:,1) = Surf_N;  % gate 1, 0 meter
    % N(:,2) = Surf_N;  % gate 2, 75 meter
    % N_error(:,1) = N_error(:,5);  % gate 1, 0 meter
    % N_error(:,2) = N_error(:,5);  % gate 2, 75 meter 
    % N(:,3) = Surf_N;  % gate 3, 150 meters
    % N(:,4) = Surf_N;  % gate 4, 225 meters
  end

  
  if flag.troubleshoot == 1
    figure(502)
    plot(N(round(p_hour/24*size(N,1)),:).*1e6./6.022E23.*18.015, range, 'k')
    ylim([0 4e3])
    hold on
    plot(N_avg(round(p_hour/24*size(N_avg,1)),:).*1e6./6.022E23.*18.015, range, 'bo')
    plot(N_masked(round(p_hour/24*size(N_masked,1)),:).*1e6./6.022E23.*18.015, range, 'r+')
    ylim([0 4e3])
    hold off
    title('Absolute Humidty (black), averaged (blue points), avergaged after mask (red +)')
  end
    

    
% if you want to mask the data use this
if flag.mask_data == 1
  N_avg = N_masked;
end

  
%% mark gaps in data with NaNs (this should be improved as it uses a zero change in wavelength as an indication that the data has a gap)
if flag.mark_gaps ==1
  %lambda_blank=lambda_all;
  %lambda_blank(diff(lambda_all,1)==0) = NaN;  
  lambda_blank = background_off;
  lambda_blank(diff(background_off,1)==0) = NaN; 
  r_blank = bsxfun(@times, ones(length(lambda_all),length(range)), lambda_blank);
  lambda_all_off(isnan(lambda_blank)) = NaN; 
  T_bench(isnan(lambda_blank)) = NaN; 
  T_base(isnan(lambda_blank)) = NaN; 
  I_on(isnan(lambda_blank)) = NaN; 
  I_off(isnan(lambda_blank)) = NaN;     
  P_on(isnan(lambda_blank)) = NaN;
  P_off(isnan(lambda_blank)) = NaN;
  if flag.WS ==1
    Surf_T(isnan(lambda_blank)) = NaN; 
    Surf_P(isnan(lambda_blank)) = NaN; 
    Surf_RH(isnan(lambda_blank)) = NaN; 
    Surf_AH(isnan(lambda_blank)) = NaN; 
    Surf_N(isnan(lambda_blank)) = NaN; 
  end
  background_on(isnan(lambda_blank)) = NaN; 
  background_off(isnan(lambda_blank)) = NaN; 
  Offline_Temp_Spatial_Avg(isnan(r_blank)) = NaN; 
  Online_Temp_Spatial_Avg(isnan(r_blank)) = NaN;
  sigma_on_total(isnan(r_blank)) = NaN; 
  sigma_off_total(isnan(r_blank)) = NaN;
  RB(isnan(r_blank)) = NaN; 
  RB_on(isnan(r_blank)) = NaN;
  N_masked(isnan(r_blank)) = NaN; 
  N_avg(isnan(r_blank)) = NaN;  
  N_error(isnan(r_blank)) = NaN;   
end

  year = strread(folder_in(end-7:end-4), '%4f', 1); 
 % year = 2000+year;
  time_new = time_grid;%+datenum(year,1,0);
  date_plot = datestr(nanmean(time_new), 'dd mmm yyyy') % this was changed
 % date_plot = datestr(mean(time_new), 'dd mmm yyyy') % this was changed
  
 % OD is - ln(I/I.o), since offline is not the same as online it needs to
 % scaled by the first few good gates -- choose 300 m to 450 m
 scale_factor = nanmean(Online_Temp_Spatial_Avg(:,300/gate:525/gate),2)./nanmean(Offline_Temp_Spatial_Avg(:,300/gate:450/gate),2);
% scale_factor = mean(Online_Temp_Spatial_Avg(:,300/gate:525/gate),2)./mean(Offline_Temp_Spatial_Avg(:,300/gate:450/gate),2);
 scale = bsxfun(@times, Offline_Temp_Spatial_Avg, scale_factor);
 OD = -(log(Online_Temp_Spatial_Avg./scale)); % calculate column optical depth
  
  

% decimate data in time to final array size 
if flag.decimate == 1 
    decimate_time = ave_time.wv/ave_time.gr/2; %ave_time.wv/ave_time.gr;
    decimate_range = spatial_average1/2; %spatial_average1;
    decimate_range = 1; % keep native gate spacing 
    % average RB data before decimating
      RB = nanmoving_average(RB,decimate_time/2,1,flag.int);
      RB_on = nanmoving_average(RB_on,decimate_time/2,1,flag.int);
    % then decimate
    RB_on =  RB_on(1:decimate_time:end, 1:decimate_range:end); 
    RB =  RB(1:decimate_time:end, 1:decimate_range:end); 
    background_off = background_off(1:decimate_time:end, 1:decimate_range:end);  
    background_on = background_on(1:decimate_time:end, 1:decimate_range:end); 
    Online_Temp_Spatial_Avg = Online_Temp_Spatial_Avg(1:decimate_time:end, 1:decimate_range:end); 
    Offline_Temp_Spatial_Avg = Offline_Temp_Spatial_Avg(1:decimate_time:end, 1:decimate_range:end); 
    OD = OD(1:decimate_time:end, 1:decimate_range:end); 
    N_error = N_error(1:decimate_time:end, 1:decimate_range:end); 
    N_avg = N_avg(1:decimate_time:end, 1:decimate_range:end); 
    N_masked = N_masked(1:decimate_time:end, 1:decimate_range:end); 
    lambda_all = lambda_all(1:decimate_time:end);
    if flag.mark_gaps ==1
     lambda_blank = lambda_blank(1:decimate_time:end);
    end
    lambda_all_off = lambda_all_off(1:decimate_time:end);
    P=P(1:decimate_range:end); 
    T=T(1:decimate_range:end); 
    time_new = time_new(1:decimate_time:end);
    range= range(1:decimate_range:end);
    I_on = I_on(1:decimate_time:end);    
    I_off = I_off(1:decimate_time:end);
    P_on = P_on(1:decimate_time:end);
    P_off = P_off(1:decimate_time:end);
    T_bench = T_bench(1:decimate_time:end);
    T_base = T_base(1:decimate_time:end);
    if flag.WS ==1
        Surf_T=Surf_T(1:decimate_time:end);
        Surf_P=Surf_P(1:decimate_time:end);
        Surf_RH=Surf_RH(1:decimate_time:end);
        Surf_AH=Surf_AH(1:decimate_time:end);
        Surf_N=Surf_N(1:decimate_time:end);
    end
end
    
  background_off = background_off/(MCS.bin_duration*1e-9*MCS.accum*(1-switch_ratio)); % change background to counts/sec
  background_on = background_on/(MCS.bin_duration*1e-9*MCS.accum*switch_ratio); % change background to counts/sec
  

%% save data
  
 if flag.save_data == 1
  %cd('/Users/spuler/Desktop/WV_DIAL_data') % point to the directory where data is stored 
  %cd('/Volumes/documents/WV_DIAL_data/processed_data') % point to the directory where data is stored 
   cd(write_data_folder)
   if flag.near == 1
     name=strcat(date_plot,'_near');   
   else
    name=strcat(date_plot);
   end
   
   if flag.WS == 1
     save(name, 'N_avg', 'RB', 'range', 'time_new', 'T', 'P', 'OD', 'background_off', 'background_on', 'profiles2ave', 'N_error',...
     'Surf_T', 'Surf_P', 'Surf_RH', 'Surf_AH', 'I_on', 'I_off', 'P_on', 'P_off', 'T_bench', 'T_base', 'lambda_all', 'lambda_all_off', 'gate')
   else
     save(name, 'N_avg', 'RB', 'range', 'time_new', 'T', 'P', 'OD', 'background_off', 'background_on', 'profiles2ave', 'N_error', ...
     'I_on', 'I_off', 'P_on', 'P_off', 'T_bench', 'T_base', 'lambda_all', 'lambda_all_off', 'gate')
   end
 end
 
 if flag.save_netCDF == 1  % save the data as an nc file
    % convert NaN fill values (and Inf) to -1 flag
    time_new(isnan(time_new)==1) = -1;  % NaNs converted to -1 flag
    N_avg(isnan(N_avg)==1) = -1; 
    N_error(isnan(N_error)==1) = -1; 
    N_avg(isinf(N_avg)==1) = -1;  
    N_error(isinf(N_error)==1) = -1;  
    Offline_Temp_Spatial_Avg(isinf(Offline_Temp_Spatial_Avg)==1) = -1;  
    Online_Temp_Spatial_Avg(isinf(Online_Temp_Spatial_Avg)==1) = -1; 
    time_unix = (time_new-datenum(1970,1,1))*86400; % convert to unix time
   
    cdf_name = strcat('wv_dial.', datestr(date_plot, 'yymmdd'));
    ncid = netcdf.create([cdf_name '.nc'],'CLOBBER');       
    % define the dimensions and variables
    % netcdf.reDef(ncid);
    dimid1 = netcdf.defDim(ncid,'time',netcdf.getConstant('NC_UNLIMITED'));
    %dimid1 = netcdf.defDim(ncid,'time', length(time_new));
    dimid2 = netcdf.defDim(ncid, 'range', length(range));
    dimid3 = netcdf.defDim(ncid, 'lambda', length(lambda_on));
    
    myvarID1 = netcdf.defVar(ncid,'time','double',dimid1);
      netcdf.putAtt(ncid, myvarID1, 'units', 'days since January 0, 0000')
    myvarID2 = netcdf.defVar(ncid,'range','float',dimid2);
      netcdf.putAtt(ncid, myvarID2, 'units', 'meters')
    myvarID3 = netcdf.defVar(ncid,'N_avg','float',[dimid2 dimid1]);
      netcdf.putAtt(ncid, myvarID3, 'long_name', 'water_vapor_number_density')
      netcdf.putAtt(ncid, myvarID3, 'units', 'molecules/cm^3')
      netcdf.putAtt(ncid, myvarID3, 'FillValue', '-1')
    myvarID4 = netcdf.defVar(ncid,'N_error','float',[dimid2 dimid1]);
      netcdf.putAtt(ncid, myvarID4, 'long_name', 'water_vapor_number_density_error')
      netcdf.putAtt(ncid, myvarID4, 'units', 'molecules/cm^3')
      netcdf.putAtt(ncid, myvarID4, 'FillValue', '-1')
    myvarID5 = netcdf.defVar(ncid,'P','float', dimid2);
      netcdf.putAtt(ncid, myvarID5, 'long_name', 'pressure')
      netcdf.putAtt(ncid, myvarID5, 'units', 'atm')
    myvarID6 = netcdf.defVar(ncid,'T','float', dimid2);
      netcdf.putAtt(ncid, myvarID6, 'long_name', 'temperature')
      netcdf.putAtt(ncid, myvarID6, 'units', 'degK')
   myvarID7 = netcdf.defVar(ncid,'RB','float',[dimid2 dimid1]);
      netcdf.putAtt(ncid, myvarID7, 'long_name', 'relative_backscatter')
      netcdf.putAtt(ncid, myvarID7, 'units', 'arbitrary_units')
      netcdf.putAtt(ncid, myvarID7, 'FillValue', '-1')
   myvarID8 = netcdf.defVar(ncid,'time_unix','double',dimid1);
      netcdf.putAtt(ncid, myvarID8, 'long_name', 'unix time')
      netcdf.putAtt(ncid, myvarID8, 'units', 'seconds since 00:00:00 UTC, January 1, 1970')
      netcdf.putAtt(ncid, myvarID8, 'FillValue', '-1')
   myvarID9 = netcdf.defVar(ncid,'Offline_Temp_Spatial_Avg','float',[dimid2 dimid1]);
      netcdf.putAtt(ncid, myvarID9, 'long_name', 'offline counts')
      netcdf.putAtt(ncid, myvarID9, 'units', 'counts')
      netcdf.putAtt(ncid, myvarID9, 'FillValue', '-1')
    myvarID10 = netcdf.defVar(ncid,'Online_Temp_Spatial_Avg','float',[dimid2 dimid1]);
      netcdf.putAtt(ncid, myvarID10, 'long_name', 'online counts')
      netcdf.putAtt(ncid, myvarID10, 'units', 'counts')
      netcdf.putAtt(ncid, myvarID10, 'FillValue', '-1')
    myvarID11 = netcdf.defVar(ncid,'lambda_on','float', dimid3);
      netcdf.putAtt(ncid, myvarID11, 'long_name', 'online wavelength')
      netcdf.putAtt(ncid, myvarID11, 'units', 'nm')
      netcdf.putAtt(ncid, myvarID11, 'FillValue', '-1')
    myvarID12 = netcdf.defVar(ncid,'lambda_off','float', dimid3);
      netcdf.putAtt(ncid, myvarID11, 'long_name', 'offline wavelength')
      netcdf.putAtt(ncid, myvarID11, 'units', 'nm')
      netcdf.putAtt(ncid, myvarID11, 'FillValue', '-1')
      
    netcdf.endDef(ncid)  

    % save the variables to the file
    
    netcdf.putVar(ncid,myvarID1,0,length(time_new),time_new);   
    netcdf.putVar(ncid,myvarID2,range);   
    netcdf.putVar(ncid,myvarID3,[0,0],size(N_avg'),N_avg');  
    netcdf.putVar(ncid,myvarID4,[0,0],size(N_error'),N_error');
    netcdf.putVar(ncid,myvarID5,P);   
    netcdf.putVar(ncid,myvarID6,T);
    netcdf.putVar(ncid,myvarID7,[0,0],size(RB'),RB');
    netcdf.putVar(ncid,myvarID8,time_unix);
    netcdf.putVar(ncid,myvarID9,[0,0],size(Offline_Temp_Spatial_Avg'),Offline_Temp_Spatial_Avg');  
    netcdf.putVar(ncid,myvarID10,[0,0],size(Online_Temp_Spatial_Avg'),Online_Temp_Spatial_Avg');
    netcdf.putVar(ncid,myvarID11,lambda_on); 
    netcdf.putVar(ncid,myvarID12,lambda_off); 
    
    
  netcdf.close(ncid);

  % original data save to ascii format for Dave Turner   
  %   name=strcat(date, 'wv_number_density.txt');
  %   wv_number_density = [time_new, double(N_avg)]; 
  %   save(name, 'wv_number_density', '-ascii', '-double'); 
  %   name=strcat(date, 'wv_number_density_error.txt');  
  %   wv_number_density_error = [time_new, double(N_error)]; 
  %   save(name, 'wv_number_density_error', '-ascii', '-double');
  
 end
 
 %% plot data

plot_size1 = [scrsz(4)/1.5 scrsz(4)/10 scrsz(3)/1.5 scrsz(4)/3];
plot_size2 = [1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]; 

x = (time_new)';
y = (range./1e3);
font_size = 14;
xData =  linspace(fix(min(time_new)),  ceil(max(time_new)), 25);

  %plot relative backscatter
  %figure1 = figure('Position',[scrsz(4)/1.5 scrsz(4)/10 scrsz(3)/1.5 scrsz(4)/1.5]);
  figure1 = figure('visible', 'off','Position',[scrsz(4)/1.5 scrsz(4)/10 scrsz(3)/1.5 scrsz(4)/1.5]);
  set(figure1, 'visible', 'off', 'PaperUnits', 'points', 'PaperPosition', [0 0 1280 800]);
  if flag.plot_data == 1 
    set(figure1, 'visible', 'on');
  end
  subplot1=subplot(2,1,1,'Parent',figure1);
  box(subplot1,'on');
  set(gcf,'renderer','zbuffer');
  %Z = double(log10((real(RB')./RB_scale)));
  Z = double(real(RB')./RB_scale);
  %Z = real(double(log10(Offline_Temp_Spatial_Avg_act')));
  %Z(isnan(Z)) = -1;
  h = pcolor(x,y,Z);
  set(h, 'EdgeColor', 'none');
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  set(gca, 'XTick',  xData)
  colorbar('EastOutside');
  axis([fix(min(time_new)) fix(min(time_new))+1 0 12])
  caxis([1e1 1e6]);
  datetick('x','HH','keeplimits', 'keepticks');
  colormap(C)
  %shading interp 
  hh = title({[date_plot,'  Relative Backscatter (C/ns km^2)']},'fontweight','b','fontsize',font_size);
  P_t = get(hh, 'Position');
  set(hh,'Position', [P_t(1) P_t(2)+0.2 P_t(3)])
  xlabel('Time (UTC)','fontweight','b','fontsize',font_size);
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size);
  set(gca,'Fontsize',font_size,'Fontweight','b');
  set(gca,'Zscale', 'log')
  set(gca,'Colorscale', 'log')
  set(gca,'Zscale', 'linear')


  %plot water vapor in g/m^3
  subplot1=subplot(2,1,2,'Parent',figure1);
  box(subplot1,'on'); %(number density in mol/cm3)(1e6 cm3/m3)/(N_A mol/mole)*(18g/mole)
  set(gcf,'renderer','zbuffer');
  Z = double(real(N_avg'.*1e6./6.022E23.*18.015));
  %Z = real(double(log10(Online_Temp_Spatial_Avg_act')));
  %Z(isnan(Z)) = -1;
  h = pcolor(x,y,Z);
  set(h, 'EdgeColor', 'none');
  set(gca, 'XTick',  xData)
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  colorbar('EastOutside');
  axis([fix(min(time_new)) fix(min(time_new))+1 0 6])
  caxis([0 6]);
  datetick('x','HH','keeplimits', 'keepticks');
  colormap(C)
  %shading interp
  hh = title({[date_plot,'  Water Vapor (g/m^{3})']},'fontweight','b','fontsize',font_size);
  P_t = get(hh, 'Position');
  set(hh,'Position', [P_t(1) P_t(2)+0.2 P_t(3)])
  xlabel('Time (UTC)','fontweight','b','fontsize',font_size); 
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
  set(gca,'Fontsize',font_size,'Fontweight','b');
 
 %link the x axis for all 2 subplots
 %ax(1)=subplot(2,1,1);
 %ax(2)=subplot(2,1,2);
 %linkaxes([ax(1),ax(2)],'xy')

 if flag.save_quicklook == 1
  cd(write_data_folder) % point to the directory where data is stored 
  date_save=datestr(nanmean(time_new), 'yyyymmdd');
%  date_save=datestr(mean(time_new), 'yyyymmdd');
% save the image as a PNG to the local data folder 
  %name1=strcat('lidar.NCAR-WV-DIAL.', date, '0000.', folder_CH, '.png'); 
  %name=strcat('lidar.',node,'-WV-DIAL.', date_save, '0000.OL.png'); 
  print(figure1, name, '-dpng', '-r0') % set the resolution as 300 dpi
  if flag.save_catalog == 1 % upload figure to the field catalog
    test=ftp('catalog.eol.ucar.edu', 'anonymous', 'spuler@ucar.edu')
    %cd(test,'/pub/incoming/catalog/frappe');
    %cd(test,'/pub/incoming/catalog/pecan');
    %cd(test,'/pub/incoming/catalog/operations');
    %cd(test,'/pub/incoming/catalog/perdigao');
    cd(test, catalog);
    mput(test, name);
    dir(test)
    close(test);
  end
end

if flag.plot_data == 1 
  % plot the relative error
  figure10 = figure('Position',plot_size1);
  set(gcf,'renderer','zbuffer');
  Z = real(double(abs(N_error./N_avg)*100))';
  %Z(isnan(Z)) = -1;
  h = pcolor(x,y,Z);
  set(h, 'EdgeColor', 'none'); 
  axis xy; colorbar('EastOutside'); caxis([0 300]);
  title({[date_plot,' relative error']},...
     'fontweight','b','fontsize',font_size)
  ylabel('Altitude (km)','fontweight','b','fontsize',font_size); 
  axis([fix(min(time_new)) fix(min(time_new))+1 0 6])
  datetick('x','HH','keeplimits');
  colormap(C)
  set(gca,'Fontsize',font_size,'Fontweight','b');

  % plot the masked data
  figure11 = figure('Position',plot_size1);
  set(gcf,'renderer','zbuffer');
  Z = real(double(N_masked.*1e6./6.022E23.*18.015))';
  %Z(isnan(Z)) = -1;
  h = pcolor(x,y,Z);
  set(h, 'EdgeColor', 'none'); 
  axis xy; colorbar('EastOutside'); caxis([0 18]);
  title({[date_plot,' masked wv']},...
     'fontweight','b','fontsize',font_size)
  ylabel('Altitude (km)','fontweight','b','fontsize',font_size); 
  axis([fix(min(time_new)) fix(min(time_new))+1 0 6])
  datetick('x','HH','keeplimits');
  colormap(C)
  set(gca,'Fontsize',font_size,'Fontweight','b');

  % plot the calculated number density
  figure12 = figure('Position',plot_size1);
  set(gcf,'renderer','zbuffer');
  Z = double(real(log10(N_avg)'));
  %Z(isnan(Z)) = -1;
  h = pcolor(x,y,Z);
  set(h, 'EdgeColor', 'none'); 
  axis xy; colorbar('EastOutside'); caxis([15 19]);  % model assume 2E17 at ground and 2E15 at 8km 
  title({[date_plot,' number density']},...
     'fontweight','b','fontsize',font_size)
  ylabel('Altitude (km)','fontweight','b','fontsize',font_size); 
  axis([fix(min(time_new)) fix(min(time_new))+1 0 6])
  datetick('x','HH','keeplimits');
  colormap(C)
  set(gca,'Fontsize',font_size,'Fontweight','b');

  % plot the number density error
  figure13 = figure('Position',plot_size1);
  set(gcf,'renderer','zbuffer');
  Z = double(real(log10(N_error)'));
  %Z(isnan(Z)) = -1;
  h = pcolor(x,y,Z);
  set(h, 'EdgeColor', 'none');
  %shading flat;
  axis xy; colorbar('EastOutside'); caxis([15 19]);  % model assumes there will be ~2E17 at ground and 2E15 at 8km 
  title({[date_plot,' number density error']},...
     'fontweight','b','fontsize', font_size)
  ylabel('Altitude (km)','fontweight','b','fontsize',font_size); 
  axis([fix(min(time_new)) fix(min(time_new))+1 0 6])
  datetick('x','HH','keeplimits');
  colormap(C)
  set(gca,'Fontsize',font_size,'Fontweight','b');


  % plot the online wavelength
  figure20 = figure('Position', plot_size2);
  plot(time_new, lambda_all); 
  datetick('x','HH','keeplimits');
  title({[date_plot,'  Online wavelength']},...
     'fontweight','b','fontsize',20)
  ylim([lambda_on-.001 lambda_on+.001])
 
  % plot the offline wavelength
  figure21 = figure('Position', plot_size2);
  plot(time_new, lambda_all_off);
  datetick('x','HH','keeplimits');
  title({[date_plot,'  Offline wavelength']},...
   'fontweight','b','fontsize',20)
  ylim([lambda_off-.001 lambda_off+.001])

 figure(22)
 OD(OD == -Inf) = NaN;
 OD(OD == Inf) = NaN;
 plot(nanmean(OD,1), range)
 %plot(nanmean(OD(5:650,:)),)
 xlim([0 2])

%receiver_scale_factor=0.75;

figure(122)
profile_start = round(5.2/24.*size(N_avg,1));  
profile_end = round(5.3/24.*size(N_avg,1));  
semilogx(N_avg(profile_start:profile_end,:), range, 'r') % molecular
hold on
semilogx(nanmean(N_avg(profile_start:profile_end,:)), range, 'b') % combined 
hold off
ylim([0 9e3]);
%xlim([1e-2 1e6]);


 if flag.WS ==1
  figure(201)
   plot(Surf_T)
   hold on
   plot(T_bench)
   plot(T_base)
   hold off
   title('internal bench and external surface temp, C')
   grid on
  
  figure(202)
   plot(Surf_P)
   title('surface pressure, atm')
  
  figure(203)
   plot(Surf_RH)
   title('surface relative humidity, atm')
   figure(204)
   plot(Surf_AH)
   title('surface absolute humidity, atm')
 %  figure(206)
 %  plotyy([time_new', time_new'], [(P_on)', (P_off)'], [time_new', time_new',time_new'], [T_bench', T_base', Surf_T']);
 %  title('on and offline seed current (mA) and rel transmitted power (left), inside and outside temp, C (right)')
 end
   figure(205)
   plot(time_new, P_on, 'r') 
   hold on
   plot(time_new, P_off, 'k')
   legend('P_{on}', 'P_{off}') 
   hold off 
   figure(206)
   plot(time_new, P_on./mean(P_off), 'r') 
   hold on
   plot(time_new, P_off./mean(P_off), 'k')
   legend('P_{on}', 'P_{off}') 
   hold off 
end
 
toc
cd(dd) % point back to original directory


%% use this for troubleshooting raw data

%troubleshoot = 0;
if flag.troubleshoot == 1;

 %figure(3)
 %semilogx(N_avg(fix(end*14/24),:), range(1,1:end-1)./1e3, 'r')
 %ylim([0 6])
 %title('DIAL profile at about 14:00 MDT')
 %xlabel('Water Vapor Number Density (cm^{-3})');
 %ylabel('range(km)'); 

%column integrated water vapor...just testing
% pw=sum(N_avg(:,750/7.5:3500/7.5),2)./6.022e23*225*100;

% find high gradient errors that occur near cloud boundaries   
[FX,FY] = gradient(N_avg);
figure('Position',plot_size1);
imagesc(time_new,range./1e3, FX');
axis xy; colorbar('EastOutside'); caxis([-5e16 5e16]);  % model assumes there will be ~2E17 at ground and 2E15 at 8km 
title({[date,' Vertical Gradient']},...
     'fontweight','b','fontsize',16)
ylabel('Altitude (km)','fontweight','b','fontsize',20); 
datetick('x','HH','keeplimits');
axis([fix(min(time_new)) fix(min(time_new))+1 0 7])
colormap(C)  

figure('Position',plot_size1);
imagesc(time_new,range./1e3, FY');
axis xy; colorbar('EastOutside'); caxis([-5e16 5e16]);  % model assumes there will be ~2E17 at ground and 2E15 at 8km 
title({[date_plot,' Horizontal Gradient']},...
     'fontweight','b','fontsize',16)
ylabel('Altitude (km)','fontweight','b','fontsize',20); 
datetick('x','HH','keeplimits');
axis([fix(min(time_new)) fix(min(time_new))+1 0 7])
colormap(C)  
    
 %plot column OD
 figure('Position',[scrsz(4)/2 scrsz(4)/10 scrsz(3)/1.5 scrsz(4)/2])
 imagesc(time_new,range./1e3,real(OD')); 
 set(gca,'Fontsize',30,'Fontweight','b');
 axis xy;
 colorbar('EastOutside');
 axis([fix(min(time_new)) fix(min(time_new))+1 0 12])
 caxis([-0.1 2]);
 datetick('x','HH','keeplimits');
 colormap(C)
 title({[date_plot,'  Column Optical Depth']},...
      'fontweight','b','fontsize',30)
 xlabel('Time (UTC)','fontweight','b','fontsize',30); 
 ylabel('Height (km, AGL)','fontweight','b','fontsize',30); 
 set(gca,'XMinorTick','off');
 set(gca,'YMinorTick','off');
 set(gca,'Fontsize',30,'Fontweight','b');
    
    
 % plot the offline background
 figure('Position',[scrsz(4)/2 scrsz(4)/10 scrsz(3)/1.5 scrsz(4)/2])
 semilogy(time_new, background_off, time_new, background_on); %
 datetick('x','HH','keeplimits');
  title({[date_plot,'  Offline background C/s']},...
     'fontweight','b','fontsize',30)
 ylim([1e2  1e7])
 hold on
  semilogy(time_new, 5e6, 'green')   % Add a horizontal line to show the 5Mc/s linearity limit
 hold off

 
 % plot and save the cloud base (added 11-Sept-2014)
 figure('Position',[scrsz(4)/2 scrsz(4)/10 scrsz(3)/1.5 scrsz(4)/2])
 cloud_base = (diff(RB) > 25);
 for i=1:size(cloud_base,1)
   if isempty(find(cloud_base(i,:), 1, 'first'))==1
     cloud_base_idx(i,:) = NaN;
   else
     cloud_base_idx(i,:) = find(cloud_base(i,:), 1, 'first')*gate;
   end
 end
 plot(time_new(1:end-1), cloud_base_idx);
 datetick('x','HH','keeplimits');
 title({[date_plot,'  Cloud Base']},'fontweight','b','fontsize',30)
 ylim([0  16000])
 ylabel('Height (m, AGL)','fontweight','b','fontsize',30); 
 xlabel('UTC','fontweight','b','fontsize',30); 
% name=strcat(date, 'cloudbase.txt');
% cloud_base = [time_new(1:end-1), cloud_base_idx];
% cloud_base(isnan(cloud_base)==1) = -1;  % Dave Turner wants the NaNs converted to -1 flag
% save(name, 'cloud_base', '-ascii', '-double');
    

  
end


%figure(100)
%semilogy(mean(Offline_Raw_Data(1:300,8:566),1))
%grid on




