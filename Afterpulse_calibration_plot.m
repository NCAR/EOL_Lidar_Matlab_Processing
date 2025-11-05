
header = 10; % not sure if this is needed
% show all the data summed to find calibration start stop periods
WV_off_sum_counts=sum(data_wv_off(:,header :end),2);
figure(1)
semilogy(WV_off_sum_counts)

AP_start = 14896;
AP_stop = 15370;

AP_start = 17342;
AP_stop = 18101;

profiles =  AP_stop-AP_start;
gap = 500; 

gate=150*MCS.bin_duration/1000;
range=0:gate:(size(data_wv_off,2)-header)*gate;

AP_off=sum(data_wv_off(AP_start:AP_stop,header:end));
AP_on=sum(data_wv_on(AP_start:AP_stop,header:end));
Off_raw=sum(data_wv_off(AP_start-gap-profiles:AP_start-gap,header:end));
On_raw=sum(data_wv_on(AP_start-gap-profiles:AP_start-gap,header:end));

figure(2)
semilogy(range,AP_off-median(AP_off),'r', 'DisplayName', 'AP_{off-sub}'); 
hold on
semilogy(range, AP_off,'r--', 'DisplayName', 'AP_{off}'); 
semilogy(range,Off_raw-median(Off_raw),'k', 'DisplayName', 'N_{off-sub}'); 
semilogy(range,Off_raw,'k--', 'DisplayName', 'N_{off}');
semilogy(range,On_raw-median(On_raw),'b', 'DisplayName', 'N_{on-sub}'); 
semilogy(range,On_raw,'b--', 'DisplayName', 'N_{on}');
hold off
grid on
xlabel('range (m)')
ylabel('background subtracted counts/bin')
xlim([0 6000])
ylim([1 1e5])
legend('Location', 'northeast','FontSize', 12);
%title('29 Oct MPD04 Afterpulse Cal, 1\mus pulse');
title('04 Nov MPD04 Afterpulse Cal, 1\mus pulse');
