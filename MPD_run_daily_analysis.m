%% MPD_run_daily_analysis.m: Unified Daily Lidar Data Processing Pipeline

function [] = MPD_run_daily_analysis(channels, correction, node, daystr, global_config_struct)


% --- 0. Initialization & Configuration Extraction ---

% Extract key processing times for calculation
read_time_in = global_config_struct.processing.read_time_in; 
ave_time.wv = global_config_struct.processing.ave_time_wv; 
ave_time.rb = global_config_struct.processing.ave_time_rb; 
ave_time.gr = global_config_struct.processing.ave_time_gr; % CRITICAL FIX: Time grid interval structure

% Calculate the actual number of profiles for averaging (profiles2ave)
profiles2ave.wv = 2*round(((ave_time.wv*60/read_time_in)+1)/2);
profiles2ave.rb = 2*round(((ave_time.rb*60/read_time_in)+1)/2);

% Extract remaining job-specific parameters (p_hour)
p_hour = global_config_struct.processing.p_hour;
catalog = global_config_struct.paths.catalog;

flag.save_data = global_config_struct.job.save_data; 
flag.save_netCDF  = global_config_struct.job.save_netCDF;
flag = global_config_struct.flags;
flag.afterpulse = strcmp(correction, 'AP_ON');


% --- 1. Setup, Path, and Calibration ---
serv_path     = global_config_struct.paths.serv_path; 
cal_serv_path = global_config_struct.paths.cal_path; 

fprintf('>> Processing using path: %s\n', serv_path);

% Load Calibration & System Constants
calvals = load_system_calvals(node, daystr, cal_serv_path);
MCS = calvals.MCS;

% Prepare write folders
nodeStr = extractAfter(node, 'MPD');
write_data_folder = fullfile(serv_path, ['mpd_', nodeStr, '_processed_data'], 'Matlab');
if strcmp(correction, 'AP_ON') == 1
    write_data_folder = fullfile(write_data_folder, 'afterpulse');
end


% Define input data folder path
year_folder = daystr(1:4);
if strcmp(nodeStr, '06')
    raw_data_folder = fullfile(serv_path, 'adihsrl_data', year_folder, daystr);
else
    raw_data_folder = fullfile(serv_path, ['mpd_', nodeStr, '_data'], year_folder, daystr);
end


% --- 2. Data Ingestion (Using the unified retrieval service) ---

try
    [data_raw, MCS] = MPD_read_and_grid_data(node, raw_data_folder, channels, ...
                                             global_config_struct.processing, MCS);
catch ME
    warning(['Skipping day ', daystr, '. Data ingestion failed: ', ME.message]);
    return;
end

% === CRITICAL STABILITY CHECK: Did the ingestion service return valid data? ===
% If data_raw is not a structure, it means the underlying file reading failed.
if ~isstruct(data_raw) || isempty(fieldnames(data_raw))
    warning(['Skipping day ', daystr, '. Data ingestion failed (returned empty). Cannot proceed with analysis.']);
    return;
end
% ===========================


% Extract merged data arrays
data_wv_on = data_raw.online_merged;
data_wv_off = data_raw.offline_merged;
data_O2_on_comb = data_raw.O2online_comb_merged;
data_O2_off_comb = data_raw.O2offline_comb_merged;
data_O2_on_mol = data_raw.O2online_mol_merged;
data_O2_off_mol = data_raw.O2offline_mol_merged;


% --- 3. Analysis Orchestration (The Core Logic) ---

% A. WV DIAL (Stand-alone or Pre-O2/HSRL)
if (strcmp(channels, 'WV') == 1 || strcmp(channels, 'ALL') == 1) && strcmp(correction, 'AP_OFF') == 1
    MPD_Analysis_function_NetCDF_v6(data_wv_on, data_wv_off, raw_data_folder, daystr, MCS, write_data_folder, flag, node, calvals.wavemeter_offset,...
        profiles2ave, calvals.P0, calvals.switch_ratio, ave_time, calvals.timing_range_correction, calvals.blank_range, p_hour, catalog, calvals.Afterpulse_File, calvals.MPD_elevation, cal_serv_path);
    disp(['    --> Finished WV processing (No AP)']);
end


% B. O2 Processing 
if (strcmp(channels, 'O2') == 1 || strcmp(channels, 'ALL') == 1) 
    
    % Initialize all outputs for the O2 chain to prevent "Output argument not assigned" errors on crash
    O2_online_comb = []; O2_offline_comb = []; range = []; RB_comb = []; time_comb = [];
    Surf_T = []; Surf_P = []; O2_on_wavelength = []; O2_online_mol = []; O2_offline_mol = [];
    T = []; P = []; BSR = []; RD = []; HSRLMolecular_scan_wavelength = []; const = []; beta_m_profile = [];

    try
        % 3a. O2 Count Processing - Combined Channel (*** CALLS NEW V2 ***)
        flag.molecular = 0;
        [O2_online_comb, O2_offline_comb, range, RB_comb, time_comb, Surf_T, Surf_P, O2_on_wavelength, ~] = ...
            MPD_Analysis_function_O2_v2(data_O2_on_comb, data_O2_off_comb, raw_data_folder, daystr, MCS, write_data_folder, flag, node, calvals.wavemeter_offset,...
                                        profiles2ave, calvals.switch_ratio, ave_time, calvals.timing_range_correction, calvals.blank_range, p_hour, 1, calvals.Afterpulse_File, cal_serv_path);
        
        % 3b. O2 Count Processing - Molecular Channel (*** CALLS NEW V2 ***)
        flag.molecular = 1; 
        [O2_online_mol, O2_offline_mol, ~, ~, time_mol] = ... % Only need the Mol data, reusing time, range, etc.
            MPD_Analysis_function_O2_v2(data_O2_on_mol, data_O2_off_mol, raw_data_folder, daystr, MCS, write_data_folder, flag, node, calvals.wavemeter_offset,...
                                        profiles2ave, calvals.switch_ratio, ave_time, calvals.timing_range_correction, calvals.blank_range, p_hour, 1, calvals.Afterpulse_File, cal_serv_path);
        
        % 3c. HSRL/K-Ratio Processing (Must call new v2 version)
        [T, P, BSR, RD, HSRLMolecular_scan_wavelength, const, beta_m_profile] = ...
            Process_HSRL_K_data(O2_online_comb, O2_offline_comb, O2_online_mol, O2_offline_mol, ...
                                time_comb, range, Surf_T, Surf_P, flag, node, daystr, ...
                                calvals.Receiver_Scan_File, write_data_folder, cal_serv_path, calvals.receiver_scale_factor);
        
        % 3d. WV Re-Analysis (Must call new v2 version)
        [N_WV, N_WV_error] = MPD_WV_analysis_function_v2(data_wv_on, data_wv_off, raw_data_folder, daystr, MCS, write_data_folder, flag, node, calvals.wavemeter_offset,...
                                                        profiles2ave, T, P, calvals.switch_ratio, ave_time, calvals.timing_range_correction, calvals.blank_range, p_hour, catalog, calvals.Afterpulse_File, calvals.MPD_elevation, cal_serv_path);

        % 3e. O2 Absorption Calculation (Must call new v2 version)
        O2_absorption_v2(const, T, P, O2_online_comb, O2_offline_comb, ...
                      time_comb, range, BSR, RD, HSRLMolecular_scan_wavelength, N_WV, beta_m_profile, ...
                      O2_on_wavelength, node, daystr, write_data_folder, flag);
                          
        disp(['    --> Finished O2/HSRL processing (Correction: ', correction, ')']);

% Temporary diagnostic
if ~isscalar(flag.mark_gaps)
    fprintf('DEBUG: flag.mark_gaps is size %s\n', mat2str(size(flag.mark_gaps)));
end

    catch ME
        warning(['Core O2/HSRL Processing Chain Failed for ', daystr, '. Skipping rest of chain. Error: ', ME.message]);
    end % CORRECT END OF TRY-CATCH BLOCK
end % CORRECT END OF OUTER IF BLOCK (O2 check)

% C. WV DIAL (AP_ON)
if (strcmp(channels, 'WV') == 1 || strcmp(channels, 'ALL') == 1) && strcmp(correction, 'AP_ON') == 1
    flag.afterpulse = 1;
    MPD_Analysis_function_NetCDF_v6(data_wv_on, data_wv_off, raw_data_folder, daystr, MCS, write_data_folder, flag, node, calvals.wavemeter_offset,...
        profiles2ave, calvals.P0, calvals.switch_ratio, ave_time, calvals.timing_range_correction, calvals.blank_range, p_hour, catalog, calvals.Afterpulse_File, calvals.MPD_elevation, cal_serv_path);
    disp(['    --> Finished WV processing (AP ON)']);
end % CORRECT END OF WV AP_ON IF BLOCK

disp(['Daily processing for ', node, ' on ', daystr, ' complete.']);





end 