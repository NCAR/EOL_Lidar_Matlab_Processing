%% load_system_calvals.m: Loads and Computes System Calibration Constants

function calvals = load_system_calvals(node, date_str, cal_serv_path)
% NOTE: Relies on external function 'loadjson' from the 'jsonlab' library.

calvals.t_date = datetime(date_str, 'InputFormat', 'yyyyMMdd');
calvals.serial_date = datenum(date_str, 'yyyyMMdd');

% --- 1. Load Calibration JSON Data ---
calvals.config_file = '';
if strcmp(node, 'MPD01') == 1
    calvals.config_file = 'dial1_calvals.json';
elseif strcmp(node, 'MPD02') == 1
    calvals.config_file = 'dial2_calvals.json';
% [Image of MPD Lidar System Diagram with DIAL, HSRL, and O2 channels]
elseif strcmp(node, 'MPD03') == 1
    calvals.config_file = 'dial3_calvals.json';
elseif strcmp(node, 'MPD04') == 1
    calvals.config_file = 'dial4_calvals.json';
elseif strcmp(node, 'MPD05') == 1
    calvals.config_file = 'dial5_calvals.json';
elseif strcmp(node, 'MPD06') == 1
    calvals.config_file = 'dial6_calvals.json';
else
    error('Invalid node specified.');
end

% Construct the full path to the JSON file using the provided cal_serv_path
full_json_path = fullfile(cal_serv_path, 'eol-lidar-calvals', 'calvals', calvals.config_file);

% *** Call the required JSON loading library ***
try
    dat = loadjson(full_json_path, 'SimplifyCell', 1);
catch ME
    error('load_system_calvals:JSONLoadError', 'Failed to load JSON file at %s. Ensure path is correct and JSON file exists. Error: %s', full_json_path, ME.message);
end

% --- Helper to Find Date-Dependent Parameter (Logic copied from original MPD_read_calvals.m) ---
function value = get_date_dependent_val(dat_struct, current_date)
    value = []; 
    for j = 1:size(dat_struct, 2)
        % Using the same date format as the original code
        item_date = datetime(dat_struct(j).date, 'InputFormat', 'd-MMM-yyyy H:m');
        if current_date >= item_date
            value = dat_struct(j).value;
        end
    end
    if isempty(value) && size(dat_struct, 2) > 0
        value = dat_struct(1).value; 
    end
end

% --- 2. Extract and Set Base Calibration Values ---
calvals.P0 = get_date_dependent_val(dat.Default_P, calvals.t_date);
calvals.switch_ratio = get_date_dependent_val(dat.switch_ratio, calvals.t_date);
calvals.wavemeter_offset = get_date_dependent_val(dat.Wavemeter_offset, calvals.t_date);
calvals.laser_pulse_width = get_date_dependent_val(dat.Laser_Pulse_Width, calvals.t_date) * 1e9; % ns
calvals.laser_pulse_delay = get_date_dependent_val(dat.Laser_Pulse_Delay, calvals.t_date) * 1e9; % ns
calvals.receiver_scale_factor = get_date_dependent_val(dat.Molecular_Gain_Matlab, calvals.t_date); 
calvals.Afterpulse_File = get_date_dependent_val(dat.Afterpulse_File, calvals.t_date);
calvals.Receiver_Scan_File = get_date_dependent_val(dat.Receiver_Scan_File, calvals.t_date);
% Location values
location_data = dat.Location;
for i=1:size(location_data, 2)
    if calvals.t_date >= datetime(location_data(i).date,'InputFormat','d-MMM-yyyy H:m')
        calvals.MPD_location = location_data(i).location;
        calvals.MPD_elevation = location_data(i).elevation;
    end
end

% --- 3. MCS Parameters ---
MCS.bins = get_date_dependent_val(dat.MCS_bins, calvals.t_date);
MCS.bin_duration = get_date_dependent_val(dat.Bin_Width, calvals.t_date) * 1e9; % ns
MCS.accum = get_date_dependent_val(dat.MCS_accum, calvals.t_date);
MCS.accum_delay = get_date_dependent_val(dat.MCS_accum_delay, calvals.t_date) * 1e9; % ns
MCS_delay = get_date_dependent_val(dat.MCS_Delay, calvals.t_date) * 1e9; % ns
calvals.MCS = MCS;


% --- 4. Compute Calculated Constants ---
calvals.timing_range_correction = ...
    (MCS_delay - calvals.laser_pulse_delay + calvals.MCS.bin_duration/2 - calvals.laser_pulse_width/2) / 1000 * 150;

% Simplified Blank Range Logic (Matches logic from original MPD_read_calvals.m)
calvals.blank_range = 150; 
if strcmp(node,'MPD01')==1 && calvals.serial_date >= datenum(2023,11,29)
    calvals.blank_range = 37.5*4;
elseif (strcmp(node,'MPD02')==1 || strcmp(node,'MPD03')==1 || strcmp(node,'MPD04')==1 || strcmp(node,'MPD05')==1 || strcmp(node,'MPD06')==1)
    if calvals.serial_date >= datenum(2020,01,01) 
         calvals.blank_range = 37.5*7;
    end
end

disp(['Loaded CalVals for ', node, ' on ', date_str]);

end