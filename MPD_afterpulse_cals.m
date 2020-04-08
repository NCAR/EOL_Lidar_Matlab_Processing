%     ap.date = '10-Mar-2020';     
%     ap.SPCM_near = 'low gain channel SPCM s/n 41551 (0.9%AP) with gate on';  
%     ap.filename = 'MPD03_afterpulse_20200310';   
%     afterpulse_start = 32458; afterpulse_stop = 34413;

%    ap.date = '16-Mar-2020';
%    ap.SPCM_near = 'low gain channel SPCM s/n 41551 (0.9%AP) with gate on';  
%    ap.filename = 'MPD03_afterpulse_20200316';    
%    afterpulse_start = 32363; afterpulse_stop = 33068;
% 
%    ap.date = '24-Mar-2020'; 
%    ap.SPCM_near = 'low gain channel SPCM s/n 41551 (0.9%AP) gate off';  
%    ap.filename = 'MPD03_afterpulse_20200324';       
%    afterpulse_start = 30257; afterpulse_stop = 31334; 

%    ap.date = '31-Mar-2020'; 
%    ap.SPCM_near = 'low gain channel SPCM s/n 23375 (0.9%AP) gate off';  
%    ap.filename = 'MPD03_afterpulse_20200331';     
%    afterpulse_start = 30288; afterpulse_stop = 31181; 

   ap.date = '07-Apr-2020'; 
   ap.SPCM_near = 'low gain channel SPCM s/n 36274 (0.1%AP) gate off';  
   ap.filename = 'MPD03_afterpulse_20200407';   
   afterpulse_start = 28493; afterpulse_stop = 29463; 
   
   %Spatial averaging (range average) in bins.  
   gate = round((MCS.bin_duration*1e-9*3e8/2)*100)/100;
   
   afterpulse_num = (afterpulse_stop-afterpulse_start)+1; 

   ap_accum_counts_off = sum(data_off(afterpulse_start:afterpulse_stop,10:end))./afterpulse_num;  %accum counts per bin
   ap_accum_counts_on = sum(data_on(afterpulse_start:afterpulse_stop,10:end))./afterpulse_num; %accum counts per bin
   ap_count_rate_off = ap_accum_counts_off/MCS.accum/(MCS.bin_duration*1e-9); % counts per second
   ap_count_rate_on = ap_accum_counts_on/MCS.accum/(MCS.bin_duration*1e-9); % counts per second


   % vector in meters
   range = single(0:gate:(size(ap_accum_counts_off,2)-1)*gate);   
 

   ap_spline_off = spline(range,ap_count_rate_off,range);
   ap_spline_on = spline(range,ap_count_rate_on,range);
   back_off = mean(ap_spline_off(:,end-round(1200/gate):end),2)-0; % select last ~1200 meters to measure background
   back_on = mean(ap_spline_on(:,end-round(1200/gate):end),2)-0; % select last ~1200 meters to measure background
   afterpulse_sub_off = (bsxfun(@minus, ap_count_rate_off, back_off)); 
   afterpulse_sub_on = (bsxfun(@minus, ap_count_rate_off, back_on));
   ap_spline_sub_off = (bsxfun(@minus, ap_spline_off, back_off)); 
   ap_spline_sub_on= (bsxfun(@minus, ap_spline_on, back_on));
   
   figure(1001)
   semilogy(range, afterpulse_sub_off, 'bo')
   hold on
   semilogy(range, afterpulse_sub_on, 'ro')
   %semilogy(range, ap_spline_sub_off, 'b--')
   %semilogy(range, ap_spline_sub_on, 'r--')
   hold off
   legend('afterpulse_{off}', 'afterpulse_{on}') 
   title([ap.date,' ', ap.SPCM_near])
   ylim([1 1e8])
   xlim([0 5000])
   ylabel('phonton rate (counts/sec)')
   xlabel('range (m)')
   grid on
   
   if flag.near == 1
       
       ap_accum_counts_near_off = sum(data_near_off(afterpulse_start:afterpulse_stop,10:end))./afterpulse_num;  %accum counts per bin
       ap_accum_counts_near_on = sum(data_near_on(afterpulse_start:afterpulse_stop,10:end))./afterpulse_num; %accum counts per bin
       ap_count_rate_near_off = ap_accum_counts_near_off/MCS.accum/(MCS.bin_duration*1e-9); % counts per second
       ap_count_rate_near_on = ap_accum_counts_near_on/MCS.accum/(MCS.bin_duration*1e-9); % counts per second

       ap_spline_near_off = spline(range,ap_count_rate_near_off,range);
       ap_spline_near_on = spline(range,ap_count_rate_near_on,range);
       back_near_off = mean(ap_spline_near_off(:,end-round(1200/gate):end),2)-0; % select last ~1200 meters to measure background
       back_near_on = mean(ap_spline_near_on(:,end-round(1200/gate):end),2)-0; % select last ~1200 meters to measure background
       afterpulse_sub_near_off = (bsxfun(@minus, ap_count_rate_near_off, back_near_off)); 
       afterpulse_sub_near_on = (bsxfun(@minus, ap_count_rate_near_off, back_near_on));
       ap_spline_sub_near_off = (bsxfun(@minus, ap_spline_near_off, back_near_off)); 
       ap_spline_sub_near_on= (bsxfun(@minus, ap_spline_near_on, back_near_on));

       figure(1001)
       hold on
       semilogy(range, afterpulse_sub_near_off, 'b+')
       semilogy(range, afterpulse_sub_near_on, 'r+')
       %semilogy(range, ap_spline_sub_near_off, 'b--')
       %semilogy(range, ap_spline_sub_near_on, 'r--')
       hold off
       legend('afterpulse_{off}', 'afterpulse_{on}', 'afterpulse-near_{off}', 'afterpulse-near_{on}') 
       ylim([1 1e8])
       xlim([0 5000])
       ylabel('phonton rate (counts/sec)')
       xlabel('range (m)')
       grid on

   end
   
   
   save(ap.filename, 'ap_spline_sub_off', 'ap_spline_sub_on', 'ap_spline_sub_near_off', 'ap_spline_sub_near_on'); 

  
  