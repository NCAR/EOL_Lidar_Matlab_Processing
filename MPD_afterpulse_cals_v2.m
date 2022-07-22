% Generate an Afterpulse file

   ap.date = '20-Jul-2022'; 
   ap.filename = 'MPD01_afterpulse_20220720';  
   ap.SPCM = 'MPD01 SPCM s/n xxxxx (0.x%AP) gate on';  
   afterpulse_start = 27488; afterpulse_stop = 28478; 
   
   ap.date = '20-Jul-2022'; 
   ap.filename = 'MPD05_afterpulse_20220720';  
   ap.SPCM = 'MPD05 SPCM s/n xxxxx (0.x%AP) gate on';  
   afterpulse_start = 27714; afterpulse_stop = 28721; 
   
   ap.date = '22-Jul-2022'; 
   ap.filename = 'MPD04_afterpulse_20220722';  
   ap.SPCM = 'MPD04 SPCM s/n xxxxx (0.x%AP) gate on';  
   afterpulse_start = 1101; afterpulse_stop = 2905; 
   
   %Spatial averaging (range average) in bins.  
   gate = round((MCS.bin_duration*1e-9*3e8/2)*100)/100;
      
   afterpulse_num = (afterpulse_stop-afterpulse_start)+1; 

   ap_accum_counts_off = sum(data_wv_off(afterpulse_start:afterpulse_stop,10:end))./afterpulse_num;  %accum counts per bin
   ap_accum_counts_on = sum(data_wv_on(afterpulse_start:afterpulse_stop,10:end))./afterpulse_num; %accum counts per bin

   % apply linear correction factor to raw counts 
    t_d=37.25E-9; %Excelitas SPCM-AQRH-13 Module 24696
    % MCSC gives counts accumulated for set bin duration so convert to count rate  C/s.
    % divide by bin time in sec (e.g., 500ns) and # of acumulations (e.g., 10000)
    % e.g., 10 accumlated counts is 2000 C/s
    C_Offline = 1./(1-(t_d.*(ap_accum_counts_off./(MCS.bin_duration*1E-9*MCS.accum))));   
    C_Online = 1./(1-(t_d.*(ap_accum_counts_on./(MCS.bin_duration*1E-9*MCS.accum))));   
    ap_accum_counts_off = ap_accum_counts_off.*C_Offline;
    ap_accum_counts_on = ap_accum_counts_on.*C_Online;
   
   ap_count_rate_off = ap_accum_counts_off/MCS.accum/(MCS.bin_duration*1e-9); % counts per second
   ap_count_rate_on = ap_accum_counts_on/MCS.accum/(MCS.bin_duration*1e-9); % counts per second
   
   
   % vector in meters
   range = single(0:gate:(size(ap_accum_counts_off,2)-1)*gate);   
   ap_range = range; 
   
   ap_spline_off = spline(range, ap_count_rate_off,range);
   ap_spline_on = spline(range, ap_count_rate_on,range);
   back_off = mean(ap_spline_off(:,end-round(1200/gate):end),2)-0; % select last ~1200 meters to measure background
   back_on = mean(ap_spline_on(:,end-round(1200/gate):end),2)-0; % select last ~1200 meters to measure background
   back_off_array = back_off.*(ones(size(ap_spline_off)));
   back_on_array = back_on.*(ones(size(ap_spline_on)));
   afterpulse_sub_off = ap_count_rate_off-back_off_array; 
   afterpulse_sub_on =  ap_count_rate_on-back_on_array; 
   ap_spline_sub_off = ap_spline_off-back_off_array;
   ap_spline_sub_on= ap_spline_on - back_on_array;
     
%    figure(1000)
%    semilogy(range, ap_spline_off, 'bx')
%    hold on
%    semilogy(range, ap_spline_on, 'ro')
%    semilogy(range, back_off_array , 'b--')
%    semilogy(range, back_on_array, 'r--')
%    hold off
%    legend('afterpulse_{off}', 'afterpulse_{on}') 
%    title([ap.date,' ', ap.SPCM])
%    ylim([1 1e8])
%    xlim([0 5000])
%    ylabel('phonton rate (counts/sec)')
%    xlabel('range (m)')
%    grid on

   figure(1001)
   semilogy(range, afterpulse_sub_off, 'bx')
   hold on
   semilogy(range, afterpulse_sub_on, 'ro')
   semilogy(range, ap_spline_sub_off, 'b--')
   semilogy(range, ap_spline_sub_on, 'r--')
   hold off
   legend('afterpulse_{off}', 'afterpulse_{on}') 
   title([ap.date,' ', ap.SPCM])
   ylim([1 1e8])
   xlim([0 5000])
   ylabel('phonton rate (counts/sec)')
   xlabel('range (m)')
   grid on

   
save(ap.filename, 'ap_spline_sub_off', 'ap_spline_sub_on', 'ap_range');
   

   
   

  
  