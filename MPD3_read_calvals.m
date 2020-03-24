addpath('./jsonlab')
if strcmp(getenv('HOSTNAME'),'fog.eol.ucar.edu')
   serv_path = '/export/fog1/rsfdata/MPD/calibration/'; % when running on server
elseif strcmp(getenv('HOSTNAME'),'')
    serv_path = '../'; % running locally 
else
   serv_path = '/Volumes/eol/fog1/rsfdata/MPD/calibration/'; % 
end
dat=loadjson([strcat(serv_path, 'eol-lidar-calvals/calvals/dial3_calvals.json')],'SimplifyCell',1);

%t_date = '24-Mar-2020' % used for testing
wavemeter_offset = double(0); 
t_date = datetime(num2str(date),'InputFormat','yyMMdd')

for i=1:size(dat.Afterpulse_File,2)
if (t_date >= datetime(dat.Afterpulse_File(i).date,'InputFormat','d-MM-yyyy H:m')) == 1
    Afterpulse_File = dat.Afterpulse_File(i).value;
end
end

for i=1:size(dat.MCS_bins,2)
if (t_date >= datetime(dat.MCS_bins(i).date,'InputFormat','d-MM-yyyy H:m')) == 1
    MCS.bins = dat.MCS_bins(i).value;
end
end

for i=1:size(dat.Bin_Width,2)
if (t_date >= datetime(dat.Bin_Width(i).date,'InputFormat','d-MM-yyyy H:m')) == 1
    MCS.bin_duration = dat.Bin_Width(i).value*1e9; % convert to ns
end
end

for i=1:size(dat.MCS_accum,2)
if (t_date >= datetime(dat.MCS_accum(i).date,'InputFormat','d-MM-yyyy H:m')) == 1
    MCS.accum = dat.MCS_accum(i).value; % convert to ns
end
end

for i=1:size(dat.MCS_accum_delay,2)
if (t_date >= datetime(dat.MCS_accum_delay(i).date,'InputFormat','d-MM-yyyy H:m')) == 1
    MCS.accum_delay = dat.MCS_accum_delay(i).value*1e9; % convert to ns
end
end

for i=1:size(dat.MCS_Delay,2)
if (t_date >= datetime(dat.MCS_Delay(i).date,'InputFormat','d-MM-yyyy H:m')) == 1
    MCS.delay = dat.MCS_Delay(i).value*1e9; % convert to ns
end
end

for i=1:size(dat.Default_P,2)
if (t_date >= datetime(dat.Default_P(i).date,'InputFormat','d-MM-yyyy H:m')) == 1
    P0 = dat.Default_P(i).value;
end
end

for i=1:size(dat.Molecular_Gain_Matlab,2)
if (t_date >= datetime(dat.Molecular_Gain_Matlab(i).date,'InputFormat','d-MM-yyyy H:m')) == 1
    receiver_scale_factor = dat.Molecular_Gain_Matlab(i).value;
    diff_geo_on = dat.Molecular_Gain_Matlab(i).diff_geo;
end
end

for i=1:size(dat.switch_ratio,2)
if (t_date >= datetime(dat.switch_ratio(i).date,'InputFormat','d-MM-yyyy H:m')) == 1
    switch_ratio = dat.switch_ratio(i).value;
end
end

for i=1:size(dat.Location,2)
if (t_date >= datetime(dat.Location(i).date,'InputFormat','d-MM-yyyy H:m')) == 1
    location = dat.Location(i).location;
end
end

for i=1:size(dat.Wavemeter_offset,2)
if (t_date >= datetime(dat.Wavemeter_offset(i).date,'InputFormat','d-MM-yyyy H:m')) == 1
    wavemeter_offset = dat.Wavemeter_offset(i).value;
end
end

for i=1:size(dat.Laser_Pulse_Width,2)
if (t_date >= datetime(dat.Laser_Pulse_Width(i).date,'InputFormat','d-MM-yyyy H:m')) == 1
    laser_pulse_width = dat.Laser_Pulse_Width(i).value*1e9; % convert to ns;
end
end

for i=1:size(dat.Laser_Pulse_Delay,2)
if (t_date >= datetime(dat.Laser_Pulse_Delay(i).date,'InputFormat','d-MM-yyyy H:m')) == 1
    laser_pulse_delay = dat.Laser_Pulse_Delay(i).value*1e9; % convert to ns;
end
end

location %write the location to the screen 
wavemeter_offset %write the calibration offset to the screen 


%calcuate the accumuation time per MCS dwell 
time_per_column = MCS.accum*((MCS.bins*MCS.bin_duration)+MCS.accum_delay)/1e9; % acummulation time in seconds
profiles2ave.wv = 2*round(((ave_time.wv*60/time_per_column)+1)/2);  % 7kHz, 10k accum data rate is ~1.4s  
profiles2ave.rb = 2*round(((ave_time.rb*60/time_per_column)+1)/2); % 7kHz, 10k accum data rate is ~1.4s 
load('diff_geo_cor_170810.mat');
%timing_range_correction = ((1.25+1/2)-0.5/2)*150;  % changed hardware timing to start after pulse through
%timing_range_correction = (1.25-0.2+0.25/2-1.0/2)*150  % Delay of MCS - delay of TOSA trigger + MCS bin duration/2 - pulse duration/2
timing_range_correction = (MCS.delay-laser_pulse_delay+MCS.bin_duration/2-laser_pulse_width/2)/1000*150
blank_range = 300; % new pulse generator shifts gate timing so less outgoing pulse contamination