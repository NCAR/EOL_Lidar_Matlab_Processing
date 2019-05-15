
elevation= 310.0; %MPD05 was at 310m elevation at SGP

%d=pwd;
%cd('/Volumes/documents/WV_DIAL_data/SGP_sondes/') % point to the directory where data is stored
cd('/scr/sci/tammy/mpd/sgp/soundings/')
[sondefilename, sondedir] = uigetfile('*.*','Select the sonde file', 'MultiSelect', 'on');
jj=1;
%cd= d;
for jj = 1:size(sondefilename,2)
    cd('/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_processing/')
    Sonde_read_nc_files(jj, elevation, sondedir, sondefilename); 
    %Sonde_DIAL_comparison_funct_v6(N_H2O, sonde_top, sonde_range, t, date, T_sonde, P_sonde, sonde_stop, shift, error_threshold, Wind_speed, save_figs, ID_sonde);
    %Sonde_DIAL_comparison_funct_Python(N_H2O, sonde_top, sonde_range, t, date, T_sonde, P_sonde, sonde_stop, shift, error_threshold, Wind_speed, save_figs)
end