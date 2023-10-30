addpath('./jsonlab')
% if strcmp(getenv('HOSTNAME'),'fog.eol.ucar.edu') == 1
%    cal_serv_path = '/export/fog1/rsfdata/MPD/calibration/' % when running on server
% elseif strcmp(getenv('HOSTNAME'),'tikal.eol.ucar.edu') ==1 
%    cal_serv_path = '/h/eol/spuler/' % when running on tikal server
% elseif strcmp(getenv('HOSTNAME'),'') == 1
%     cal_serv_path = '../' % running locally 
% else
%    cal_serv_path = '/Volumes/eol/fog1/rsfdata/MPD/calibration/' % 
% end

%  blank_range = 525; % This was the original blank range in 2019
  blank_range = 300; %  This was the blank range prior to adding the WFOV receiver in summer of 2020
  blank_range = 150; %  For testing turn off the blanking

if strcmp(node,'MPD01')==1
  dat=loadjson([strcat(cal_serv_path, 'eol-lidar-calvals/calvals/dial1_calvals.json')],'SimplifyCell',1);
  if (serial_date >= 737988) % day the WFOV channel was installed
       blank_range = 37.5*7;
  end
elseif strcmp(node,'MPD02')==1
  dat=loadjson([strcat(cal_serv_path, 'eol-lidar-calvals/calvals/dial2_calvals.json')],'SimplifyCell',1);
  if (serial_date >= 738024) % day the WFOV channel was installed
       blank_range = 37.5*7  % 1µs pulse requires 6 range bins
  end
elseif strcmp(node,'MPD03')==1
  dat=loadjson([strcat(cal_serv_path, 'eol-lidar-calvals/calvals/dial3_calvals.json')],'SimplifyCell',1);
  if (serial_date >= 737916) % day the WFOV channel was installed
       blank_range = 37.5*7  % 1µs pulse, 8 range bins, based on M2HATS
  end
elseif strcmp(node,'MPD04')==1
  dat=loadjson([strcat(cal_serv_path, 'eol-lidar-calvals/calvals/dial4_calvals.json')],'SimplifyCell',1);
  if (serial_date >= 737902) % day the WFOV channel was installed
       blank_range = 37.5*7;
  end
elseif strcmp(node,'MPD05')==1
  dat=loadjson([strcat(cal_serv_path, 'eol-lidar-calvals/calvals/dial5_calvals.json')],'SimplifyCell',1);
  if (serial_date >= 738014) % day the WFOV channel was installed
       blank_range = 37.5*7; % 1µs pulse requires 6 range bins
  end
end


%t_date = '04-Jan-2022' % used for testing
wavemeter_offset = double(0); 
t_date = datetime(num2str(date),'InputFormat','yyMMdd')

if isfield(dat, 'Afterpulse_File') == 1
    for i=1:size(dat.Afterpulse_File,2)
        if size(dat.Afterpulse_File,2)== 1
            Afterpulse_File = dat.Afterpulse_File.value;
        elseif i >= 2 && strcmp(dat.Afterpulse_File(i).change_type, 'gradual') == 1
             index_n0 = t_date - datetime(dat.Afterpulse_File(i-1).date,'InputFormat','d-MMM-yyyy H:m');
             index_n1 = t_date - datetime(dat.Afterpulse_File(i).date,'InputFormat','d-MMM-yyyy H:m');
             if abs(index_n0) < abs(index_n1) == 1
                Afterpulse_File = dat.Afterpulse_File(i-1).value; 
             else
               Afterpulse_File = dat.Afterpulse_File(i).value;
             end
        elseif (t_date >= datetime(dat.Afterpulse_File(i).date,'InputFormat','d-MMM-yyyy H:m')) == 1
          Afterpulse_File = dat.Afterpulse_File(i).value;
        end
    end
    if exist('Afterpulse_File')==0  % if there is no afterpulse cal file then make a 
        Afterpulse_File = 'No Afterpulse Cal'
    end
else
     Afterpulse_File = 'No Afterpulse Cal' 
end
Afterpulse_File
  

for i=1:size(dat.MCS_bins,2)
if (t_date >= datetime(dat.MCS_bins(i).date,'InputFormat','d-MMM-yyyy H:m')) == 1
    MCS.bins = dat.MCS_bins(i).value;
end
end

for i=1:size(dat.Bin_Width,2)
if (t_date >= datetime(dat.Bin_Width(i).date,'InputFormat','d-MMM-yyyy H:m')) == 1
    MCS.bin_duration = dat.Bin_Width(i).value*1e9; % convert to ns
end
end

for i=1:size(dat.MCS_accum,2)
if (t_date >= datetime(dat.MCS_accum(i).date,'InputFormat','d-MMM-yyyy H:m')) == 1
    MCS.accum = dat.MCS_accum(i).value; % convert to ns
end
end

for i=1:size(dat.MCS_accum_delay,2)
if (t_date >= datetime(dat.MCS_accum_delay(i).date,'InputFormat','d-MMM-yyyy H:m')) == 1
    MCS.accum_delay = dat.MCS_accum_delay(i).value*1e9; % convert to ns
end
end

for i=1:size(dat.MCS_Delay,2)
if (t_date >= datetime(dat.MCS_Delay(i).date,'InputFormat','d-MMM-yyyy H:m')) == 1
    MCS.delay = dat.MCS_Delay(i).value*1e9; % convert to ns
end
end

for i=1:size(dat.Default_P,2)
if (t_date >= datetime(dat.Default_P(i).date,'InputFormat','d-MMM-yyyy H:m')) == 1
    P0 = dat.Default_P(i).value;
end
end

for i=1:size(dat.Molecular_Gain_Matlab,2)
if (t_date >= datetime(dat.Molecular_Gain_Matlab(i).date,'InputFormat','d-MMM-yyyy H:m')) == 1
    receiver_scale_factor = dat.Molecular_Gain_Matlab(i).value;
    diff_geo_on = dat.Molecular_Gain_Matlab(i).diff_geo;
end
end

for i=1:size(dat.switch_ratio,2)
if (t_date >= datetime(dat.switch_ratio(i).date,'InputFormat','d-MMM-yyyy H:m')) == 1
    switch_ratio = dat.switch_ratio(i).value;
end
end

for i=1:size(dat.Location,2)
if (t_date >= datetime(dat.Location(i).date,'InputFormat','d-MMM-yyyy H:m')) == 1
    MPD_location = dat.Location(i).location;
    MPD_elevation = dat.Location(i).elevation;
end
end


for i=1:size(dat.Wavemeter_offset,2)
if (t_date >= datetime(dat.Wavemeter_offset(i).date,'InputFormat','d-MMM-yyyy H:m')) == 1
    wavemeter_offset = dat.Wavemeter_offset(i).value;
end
end

for i=1:size(dat.Laser_Pulse_Width,2)
if (t_date >= datetime(dat.Laser_Pulse_Width(i).date,'InputFormat','d-MMM-yyyy H:m')) == 1
    laser_pulse_width = dat.Laser_Pulse_Width(i).value*1e9; % convert to ns;
end
end

for i=1:size(dat.Laser_Pulse_Delay,2)
if (t_date >= datetime(dat.Laser_Pulse_Delay(i).date,'InputFormat','d-MMM-yyyy H:m')) == 1
    laser_pulse_delay = dat.Laser_Pulse_Delay(i).value*1e9; % convert to ns;
end
end

for i=1:size(dat.Receiver_Scan_File,2)
if (t_date >= datetime(dat.Receiver_Scan_File(i).date,'InputFormat','d-MMM-yyyy H:m')) == 1
    Receiver_Scan_File = dat.Receiver_Scan_File(i).value; 
end
end


      
  

MPD_location %write the location to the screen 
MPD_elevation %write the location to the screen 
wavemeter_offset %write the calibration offset to the screen 

%calcuate the accumuation time per MCS dwell 
%time_per_column = MCS.accum*((MCS.bins*MCS.bin_duration)+MCS.accum_delay)/1e9 % acummulation time in seconds
% this should have the PRF right?  getting the wrong values
%time_per_column = 2.0 %override for now
%profiles2ave.wv = 2*round(((ave_time.wv*60/time_per_column)+1)/2)   
%profiles2ave.rb = 2*round(((ave_time.rb*60/time_per_column)+1)/2) 
load('diff_geo_cor_170810.mat');
%timing_range_correction = ((1.25+1/2)-0.5/2)*150;  % changed hardware timing to start after pulse through
%timing_range_correction = (1.25-0.2+0.5/2-1.0/2)*150;% Delay of MCS - delay of TOSA trigger + MCS bin duration/2 - pulse duration/2
timing_range_correction = (MCS.delay-laser_pulse_delay+MCS.bin_duration/2-laser_pulse_width/2)/1000*150 % where variables are in ns and correction is in m

    
