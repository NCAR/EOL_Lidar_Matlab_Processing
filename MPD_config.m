%% MPD_config.m: Centralized Configuration for MPD Processing

function config = MPD_config(start_date_in, stop_date_in)

% --- 1. Processing Dates (Can be overridden by user input) ---
if nargin < 2
    config.dates.start_str = '20251106';   
    config.dates.stop_str =  '20251107'; 
else
    config.dates.start_str = start_date_in;
    config.dates.stop_str = stop_date_in;
end
config.dates.start_day = datenum(config.dates.start_str,'yyyymmdd');
config.dates.stop_day = datenum(config.dates.stop_str,'yyyymmdd');


% --- 2. Operational Flags (Control what to do) ---
config.flags.process_data = 1;   % Set this to 0 when you ONLY want plots
config.flags.run_plots = 0;      % Set this to 1 when you ONLY want plots
config.flags.save_quicklook  = 0;
config.flags.save_data       = 1;   % Save daily .mat file (used by plots)
config.flags.save_netCDF     = 0;
config.flags.save_catalog    = 0;
config.flags.dev_mode        = 1;  % Set to 1 for fast local processing (using copied files)

% --- 3. Processing Parameters (Global to the pipeline) ---
config.processing.read_time_in     = 4;     
config.processing.ave_time_wv      = 10;    
config.processing.ave_time_rb      = 2.0;   
config.processing.ave_time_gr      = 1.0;   
config.processing.p_hour           = 20;    
config.processing.flag_mask_data   = 1;
config.processing.flag_gradient_filter = 1;
config.processing.flag_pileup      = 1; 
config.processing.flag_WS          = 1; 
config.processing.flag_OF          = 1; 
config.processing.flag_ap_quick    = 0;
config.processing.flag_mark_gaps   = 1; % sets gaps in data to NaNs
config.processing.flag_decimate    = 0; % decimate all data to half the wv resoltuion

% --- 4. System/Channel Matrix (Define the jobs to run) ---
config.systems_to_process = {
    % The job originally requested:
    struct('node', 'MPD04', 'channels', 'O2', 'correction', 'AP_OFF', 'save_data', 1, 'save_netCDF', 0);
};

% --- 6. Path Definitions (for upload/catalog) ---
config.paths.catalog = '/pub/incoming/catalog/operations';
if config.flags.dev_mode == 1
    % Set the base path to your local copy
    config.paths.raw_data_base = '/Users/Spuler/Desktop/MPD_Dev_Data/'; 
else
    % Set the base path to the server (production)
    config.paths.raw_data_base = '/Volumes/smaug1/rsfdata/MPD/';
end

% --- 5. Multi-Day Plot Parameters ---
config.multiday_plots = {
    struct('start_date', '20251106', 'stop_date', '20251107', ...
           'node', 'MPD04', 'near', 0, 'afterpulse', 0, 'skip', 1); 
};

end