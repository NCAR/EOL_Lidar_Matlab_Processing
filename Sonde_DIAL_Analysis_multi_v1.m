clear all
close all

shift = -7.5; % number of profiles shift comparison data (15 min is reasonable, ~time to reach 4 km) 
    % note a zero time shift for a 0:00 sonde will have problems b/c of averageing edge effects
average = 15; % time in minutes to collect DIAL profiles 
error_threshold =  100.0; % DIAL error threshold
save_figs = 1; 

%cd('/scr/eldora1/h2o_data/2015/sondes/') % point to the directory where data is stored 
cd('/Users/spuler/Desktop/Perdigao/sondes/') % point to the directory where data is stored 
[sondefilename, sondedir] = uigetfile('*.*','Select the sonde file', 'MultiSelect', 'on');
j=1;

for j = 1:size(sondefilename,2)
    cd('/Users/spuler/Documents/GitHub/Matlab_DIAL_processing/')
    [N_H2O, sonde_top, sonde_range, t, date, T_sonde, P_sonde, sonde_stop, Wind_speed] = Sonde_get_v3(sondedir, sondefilename{j}, average); 
    Sonde_DIAL_comparison_funct_v4(N_H2O, sonde_top, sonde_range, t, date, T_sonde, P_sonde, sonde_stop, shift, error_threshold, Wind_speed, save_figs)
    %Sonde_DIAL_comparison_funct_Python(N_H2O, sonde_top, sonde_range, t, date, T_sonde, P_sonde, sonde_stop, shift, error_threshold, Wind_speed, save_figs)
end