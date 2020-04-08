  % 10-Mar-2020 data; 
   % afterpulse_start = 32458; afterpulse_stop = 34413;
   % 16-Mar-2020 data; 
    ap_filename = 'MPD03_afterpulse_20200316'
    afterpulse_start = 32363; afterpulse_stop = 33068;
   % 24-Mar-2020 data; 
   % ap_filename = 'MPD03_afterpulse_20200324'
   % afterpulse_start = 30257; afterpulse_stop = 31334; 
   % 31-Mar-2020 data; 
   % ap_filename = 'MPD03_afterpulse_20200331';
   % afterpulse_start = 30288; afterpulse_stop = 31181; 
   % 07-Apr-2020 
   % ap_filename = 'MPD03_afterpulse_20200407';
   % afterpulse_start = 28493; afterpulse_stop = 29463; 
   
   %Spatial averaging (range average) in bins.  
   gate = round((MCS.bin_duration*1e-9*3e8/2)*100)/100;
   
   afterpulse_num = (afterpulse_stop-afterpulse_start)+1; 
   afterpulse_off = sum(data_off(afterpulse_start:afterpulse_stop,10:end))./afterpulse_num;
   afterpulse_on = sum(data_on(afterpulse_start:afterpulse_stop,10:end))./afterpulse_num;

   % vector in meters
   range = single(0:gate:(size(afterpulse_off,2)-1)*gate);   
   skip = round(10/gate);  % don't fit outgoing pulse, cut ranges < ~200 m

   ap_spline_off = spline(range(skip:end),afterpulse_off(skip:end),range);
   ap_spline_on = spline(range(skip:end),afterpulse_on(skip:end),range);
   back_off = mean(ap_spline_off(:,end-round(1200/gate):end),2)-0; % select last ~1200 meters to measure background
   back_on = mean(ap_spline_on(:,end-round(1200/gate):end),2)-0; % select last ~1200 meters to measure background
   afterpulse_sub_off = (bsxfun(@minus, afterpulse_off, back_off)); 
   afterpulse_sub_on = (bsxfun(@minus, afterpulse_on, back_on));
   ap_spline_sub_off = (bsxfun(@minus, ap_spline_off, back_off)); 
   ap_spline_sub_on= (bsxfun(@minus, ap_spline_on, back_on));
   
   figure(1001)
   semilogy(range(skip:end), afterpulse_sub_off(skip:end), 'bo')
   hold on
   semilogy(range(skip:end), afterpulse_sub_on(skip:end), 'ro')
   %semilogy(range, ap_spline_sub_off, 'b--')
   %semilogy(range, ap_spline_sub_on, 'r--')
   hold off
   legend('afterpulse_{off}', 'afterpulse_{on}') 
   %ylim([5 100])
   xlim([0 5000])
   ylabel('counts')
   xlabel('range (m)')
   grid on
   
   if flag.near == 1
       
       afterpulse_near_off = sum(data_near_off(afterpulse_start:afterpulse_stop,10:end))./afterpulse_num;
       afterpulse_near_on = sum(data_near_on(afterpulse_start:afterpulse_stop,10:end))./afterpulse_num;

       ap_spline_near_off = spline(range(skip:end),afterpulse_near_off(skip:end),range);
       ap_spline_near_on = spline(range(skip:end),afterpulse_near_on(skip:end),range);
       back_near_off = mean(ap_spline_near_off(:,end-round(1200/gate):end),2)-0; % select last ~1200 meters to measure background
       back_near_on = mean(ap_spline_near_on(:,end-round(1200/gate):end),2)-0; % select last ~1200 meters to measure background
       afterpulse_sub_near_off = (bsxfun(@minus, afterpulse_near_off, back_near_off)); 
       afterpulse_sub_near_on = (bsxfun(@minus, afterpulse_near_on, back_near_on));
       ap_spline_sub_near_off = (bsxfun(@minus, ap_spline_near_off, back_near_off)); 
       ap_spline_sub_near_on= (bsxfun(@minus, ap_spline_near_on, back_near_on));

       figure(1001)
       hold on
       semilogy(range(skip:end), afterpulse_sub_near_off(skip:end), 'b+')
       semilogy(range(skip:end), afterpulse_sub_near_on(skip:end), 'r+')
       %semilogy(range, ap_spline_sub_near_off, 'b--')
       %semilogy(range, ap_spline_sub_near_on, 'r--')
       hold off
       legend('afterpulse_{off}', 'afterpulse_{on}', 'afterpulse-near_{off}', 'afterpulse-near_{on}') 
       %ylim([5 100])
       xlim([0 5000])
       ylabel('counts')
       xlabel('range (m)')
       grid on

   end
   
   
   save(ap_filename, 'ap_spline_sub_off', 'ap_spline_sub_on', 'ap_spline_sub_near_off', 'ap_spline_sub_near_on'); 

  
  