driver_current = [90;100;110;120;130;140;150;160;170;180;190;200;210]; % driver current in mA
P_in_TSOA = [11.15;19.77;29.06;40.1;52;64;77.4;88.5;100.7;112.7;124.6;126.5;126.5]; % power in (uW)
P_out_TSOA = [14.4;16.5;18.9;20.6;21.9;22.7;23.8;24.38;24.46;24.56;24.62;24.63;24.64]; % power out (mW)
PRF= 8000; % pulse rep frequency
duty_cycle = 0.5; % online only
seed_pulse_length=1e-6; %seed pulse length 
TSOA_pulse_length=625e-9; %TSOA pulse length 
SOA_type = 'SOA-352-DBUT-SM-FC/APC';
SOA_serial = '226712';
SOA_max_current = '268 mA';
SOA_max_cw_ouput = '30 mW'; 
SOA_driver_max_current_peak = '250 mA'; 
SOA_driver_max_current_ave = '150 mA'; 

figure(1)
plot(driver_current, (P_in_TSOA*1e-6)./(PRF*duty_cycle)./seed_pulse_length*1000, 'bo')
xlabel('SOA pulse drive current (mA)')
ylabel('equivalent cw power (mW) after isolator')
ylim([0 40])
xlim([0 225])
legend([SOA_type, ' Serial #:', SOA_serial, ' Max current: ', SOA_max_current])


figure(2)
plot(driver_current, (P_out_TSOA/duty_cycle), 'ro')
xlabel('SOA pulse drive current (mA)')
ylabel('equivalent cw power (mW) after isolator')
ylim([0 80])
xlim([0 225])
legend([SOA_type, ' Serial #:', SOA_serial, ' Max current: ', SOA_max_current])