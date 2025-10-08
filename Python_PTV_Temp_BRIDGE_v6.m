clear all; close all

% --- 1. CONFIGURATION AND PATHS ---
node = 'MPD04';
% Use fullfile for robust path definition
serv_path = '/Volumes/eol/sci/mhayman';
plot_path = '/Users/spuler/Desktop';

data_dir = fullfile(serv_path, 'DIAL', 'Processed_Data', 'BRIDGE_2025', 'ptv0.2');
plot_dir = fullfile(plot_path, 'mpd', 'Plots');
addpath '/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing/matplotlib/';
addpath '/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing'

flag.save_data = 0;  %save data at end of processing (0=off 1=on)
low_range_mask = 0;

% --- PLOTTING EXCLUSION FLAGS (0=Exclude/Off, 1=Include/On) ---
flag.plot_multi_temp = 1;      % Toggle 5-Panel Temperature Figure (Figure 1)
flag.plot_multi_wv = 1;        % Toggle 7-Panel Water Vapor Figure (Figure 2)

% Flags for REDUNDANT SINGLE PANELS (Defaulted to 0/OFF)
flag.plot_single_temp_panels = 0; 
flag.plot_single_wv_panels = 0;   

% Flags for essential/less redundant panels
flag.plot_aerosol_counts = 1;  
flag.plot_uncertainty = 0;     
flag.plot_histograms_2d = 0;   % Toggle 2D Density Maps (T_Range, AH_Range)

flag.plot_1d_histograms = 0;   % Toggle 1D Histograms (T_Diff_Histogram, AH_Diff_Histogram)
% ------------------------------------

% --- DATA READING EXCLUSION FLAGS (0=Exclude/Off, 1=Include/On) ---
flag.read_temp_std = 1;        % Include reading Temperature_Standard and related fields
flag.read_wv_multi = 1;        % Include reading Absolute_Humidity_MultiPulse fields
flag.read_uncertainty = 0;     % Include reading ALL uncertainty/variance fields (T_var, AH_var, etc.)
% ------------------------------------------------------------------

% --- PLOTTING RESOLUTION SETTING ---
TARGET_PIXEL_WIDTH = 5000;     % Target output width (pixels) for visual quality (used to calculate N_DECIMATE)
flag.decimate_hist_data = 0;   % 0 = Use FULL RESOLUTION data for 2D Histograms (High Accuracy, Slower)
% -----------------------------------


disp('--- Script Start ---');
tic_total = tic; % Start total timing

cd(data_dir); % Move to data directory to select files

% --- 2. FILE SELECTION ---
[Pythonfilename, Pythondir] = uigetfile('*.*','Select the sonde file', 'MultiSelect', 'on');
num_files = size(Pythonfilename, 2);

if num_files == 0
    error('No files selected. Script terminated.');
end

% --- 3. VARIABLE DEFINITION ---
variables = cell(1, 41);
variables([1:12, 13:15, 16:18, 26:28, 30:32, 40:41]) = {
    'time', 'range', 'Temperature_PTV', 'Temperature_PTV_mask', 'Temperature_PTV_uncertainty', ...
    'Absolute_Humidity_Standard', 'Absolute_Humidity_Standard_mask', 'Absolute_Humidity_Standard_uncertainty', ...
    'Aerosol_Backscatter_Coefficient_PTV', 'Aerosol_Backscatter_Coefficient_PTV_mask', ...
    'Aerosol_Backscatter_Coefficient_PTV_uncertainty', 'Backscatter_Photon_Counts_828', ...
    'Temperature_Standard', 'Temperature_Standard_mask', 'Temperature_Standard_uncertainty', ... 
    'Absolute_Humidity_PTV', 'Absolute_Humidity_PTV_mask', 'Absolute_Humidity_PTV_uncertainty', ...
    'Absolute_Humidity_MultiPulse', 'Absolute_Humidity_MultiPulse_mask', 'Absolute_Humidity_MultiPulse_uncertainty', ...
    'Surface_Temperature', 'Surface_Pressure', 'Surface_Absolute_Humidity', ...
    'Absolute_Humidity_ERA5', 'Temperature_ERA5'
};

% --- 4. DATA READING AND MASKING LOOP (ORIGINAL LOGIC PRESERVED) ---
disp('Starting data reading and masking...');
tic_read = tic;

for jj = 1:num_files
    filename = Pythonfilename{jj};
    date_str = filename(end-15:end-8);
    n = datenum(date_str, 'yyyymmdd');
    
    ncid = netcdf.open(filename, 'NC_NOWRITE');
    

    % Read data (Base variables - always read)
    time{jj} = ncread(filename,variables{1}); 
    alt{jj} = ncread(filename,variables{2});
        
    % Get dimensions for fallback NaNs
    num_ranges_file = size(alt{jj}, 1);
    num_timesteps_file = size(time{jj}, 1);
 
    % Temperature PTV (Always read, but uncertainty is conditional)
    T{jj}  = ncread(filename,variables{3});  
    T_mask{jj} = ncread(filename,variables{4}); 
    T_model{jj}  = ncread(filename,variables{41}); 

    % Conditional: PTV Temperature Uncertainty
    if flag.read_uncertainty % Include reading T_var if flag is 1
        T_var{jj} = ncread(filename,variables{5}); 
    else
        T_var{jj} = nan(num_ranges_file, num_timesteps_file);
    end
    
    % Conditional: Temperature Standard (NEW)
    if flag.read_temp_std % Include reading T_Std if flag is 1
        T_Std{jj}  = ncread(filename,variables{13});  
        T_Std_mask{jj} = ncread(filename,variables{14}); 
        if flag.read_uncertainty
             T_Std_var{jj} = ncread(filename,variables{15}); 
        else
             T_Std_var{jj} = nan(num_ranges_file, num_timesteps_file);
        end
    else
        % Fill with NaNs if excluding (flag.read_temp_std is 0)
        T_Std{jj} = nan(num_ranges_file, num_timesteps_file);
        T_Std_mask{jj} = nan(num_ranges_file, num_timesteps_file);
        T_Std_var{jj} = nan(num_ranges_file, num_timesteps_file);
    end
    
    % Surface Data (Always read)
    T_surf{jj} =  ncread(filename,variables{30});
    P_surf{jj} =  ncread(filename,variables{31});
    AH_surf{jj} =  ncread(filename,variables{32});
    
    % Absolute Humidity Standard (Always read, but uncertainty is conditional)
    AH{jj}  = ncread(filename,variables{6});  
    AH_mask{jj} = ncread(filename,variables{7}); 
    AH_PTV{jj}  = ncread(filename,variables{16});  
    AH_PTV_mask{jj} = ncread(filename,variables{17}); 
    AH_model{jj}  = ncread(filename,variables{40});  
    
    % Conditional: AH Standard and PTV Uncertainty
    if flag.read_uncertainty % Include reading uncertainty if flag is 1
        AH_var{jj} = ncread(filename,variables{8}); 
        AH_PTV_var{jj} = ncread(filename,variables{18});
    else
        AH_var{jj} = nan(num_ranges_file, num_timesteps_file);
        AH_PTV_var{jj} = nan(num_ranges_file, num_timesteps_file);
    end

    % Conditional: Absolute Humidity MultiPulse
    if flag.read_wv_multi % Include reading MultiPulse data if flag is 1
        try
          AH_MultiPulse{jj}  = ncread(filename,variables{26});  
          AH_MultiPulse_mask{jj} = ncread(filename,variables{27}); 
          if flag.read_uncertainty
              AH_MultiPulse_var{jj} = ncread(filename,variables{28});
          else
              AH_MultiPulse_var{jj} = nan(num_ranges_file, num_timesteps_file);
          end
        catch
           warning(['MultiPulse data missing in file: ', filename]);
           AH_MultiPulse{jj} = nan(num_ranges_file, num_timesteps_file);
           AH_MultiPulse_mask{jj} = nan(num_ranges_file, num_timesteps_file);
           AH_MultiPulse_var{jj} = nan(num_ranges_file, num_timesteps_file);
        end
    else
        % Fill with NaNs if excluding
        AH_MultiPulse{jj} = nan(num_ranges_file, num_timesteps_file);
        AH_MultiPulse_mask{jj} = nan(num_ranges_file, num_timesteps_file);
        AH_MultiPulse_var{jj} = nan(num_ranges_file, num_timesteps_file);
    end
    
    % ABC/Counts (Uncertainty is conditional)
    ABC{jj}  = ncread(filename,variables{9});   
    ABC_mask{jj} = ncread(filename,variables{10}); 
    Counts{jj} = ncread(filename,variables{12});   
    
    if flag.read_uncertainty % Include reading uncertainty if flag is 1
        ABC_var{jj} = ncread(filename,variables{11});
    else
        ABC_var{jj} = nan(num_ranges_file, num_timesteps_file);
    end


    % --- MASKING ---
    T{jj}(T_mask{jj} == 1) = nan;
    T_var{jj}(T_mask{jj} == 1) = nan; 
    T_Std{jj}(T_Std_mask{jj} == 1) = nan;
    T_Std_var{jj}(T_Std_mask{jj} == 1) = nan; 
    mask_ah = (AH_mask{jj} == 1) | (AH_var{jj} > 25);
    AH{jj}(mask_ah) = nan;
    AH_var{jj}(AH_mask{jj} == 1) = nan; 
    mask_ah1 = (AH_PTV_mask{jj} == 1) | (AH_PTV_var{jj} > 25);
    AH_PTV{jj}(mask_ah1) = nan;
    AH_PTV_var{jj}(AH_PTV_mask{jj} == 1) = nan; 
    mask_ah2 = (AH_MultiPulse_mask{jj} == 1) | (AH_MultiPulse_var{jj} > 25);
    AH_MultiPulse{jj}(mask_ah2) = nan;
    AH_MultiPulse_var{jj}(AH_MultiPulse_mask{jj} == 1) = nan; 
    ABC{jj}(ABC_mask{jj} == 1) = nan;

    netcdf.close(ncid); 
    
    % Convert from Unix time to date number
    duration{jj} = n + double(time{jj}/86400); % 86400 = 3600*24
end

disp(['Data reading and masking complete. Time elapsed: ', num2str(toc(tic_read)), ' seconds.']);

% --- 5. PREALLOCATION AND COMBINATION (ORIGINAL LOGIC PRESERVED) ---
disp('Starting data combination...');
tic_combine = tic;

total_timesteps = sum(cellfun('size', duration, 1)); 
num_ranges = size(AH{1}, 1); 

% -----------------------------------------------------------
% --- FULL RESOLUTION COPIES (For Histograms if NOT decimating) ---
% -----------------------------------------------------------
comb_AH_full = nan(num_ranges, total_timesteps);
comb_AH_var_full = nan(num_ranges, total_timesteps);
comb_AH_PTV_full = nan(num_ranges, total_timesteps);
comb_AH_PTV_var_full = nan(num_ranges, total_timesteps);
comb_AH_MultiPulse_full = nan(num_ranges, total_timesteps);
comb_AH_model_full = nan(num_ranges, total_timesteps);
comb_ABC_full = nan(num_ranges, total_timesteps);
comb_T_full = nan(num_ranges, total_timesteps);
comb_T_Std_full = nan(num_ranges, total_timesteps);
comb_T_model_full = nan(num_ranges, total_timesteps);


% -----------------------------------------------------------------------
% --- DECIMATED/BASE COPIES (These are the main variables used for pcolor) ---
% -----------------------------------------------------------------------
comb_AH = nan(num_ranges, total_timesteps);
comb_AH_var = nan(num_ranges, total_timesteps);
comb_AH_PTV = nan(num_ranges, total_timesteps);
comb_AH_PTV_var = nan(num_ranges, total_timesteps);
comb_AH_MultiPulse = nan(num_ranges, total_timesteps);
comb_AH_MultiPulse_var = nan(num_ranges, total_timesteps);
comb_AH_model = nan(num_ranges, total_timesteps);
comb_ABC = nan(num_ranges, total_timesteps);
comb_ABC_var = nan(num_ranges, total_timesteps);
comb_T = nan(num_ranges, total_timesteps);
comb_T_var = nan(num_ranges, total_timesteps);
comb_T_Std = nan(num_ranges, total_timesteps); 
comb_T_Std_var = nan(num_ranges, total_timesteps); 
comb_T_model = nan(num_ranges, total_timesteps);
comb_Counts = nan(num_ranges, total_timesteps);

% Preallocate 1D arrays
comb_duration = nan(total_timesteps, 1);
comb_T_surf = nan(total_timesteps, 1); 
comb_P_surf = nan(total_timesteps, 1); 
comb_AH_surf = nan(total_timesteps, 1); 

start_col = 1;

for jj = 1:num_files
    current_timesteps = size(AH{jj}, 2); 
    end_col = start_col + current_timesteps - 1;

    % Fill 1D data
    comb_duration(start_col:end_col) = duration{jj};
    comb_T_surf(start_col:end_col) = T_surf{jj};
    comb_P_surf(start_col:end_col) = P_surf{jj};
    comb_AH_surf(start_col:end_col) = AH_surf{jj};
    
    % Fill 2D data:
    
    % Full Resolution Copies (Filled unconditionally for histogram option)
    comb_AH_full(:, start_col:end_col) = AH{jj};
    comb_AH_var_full(:, start_col:end_col) = AH_var{jj};
    comb_AH_PTV_full(:, start_col:end_col) = AH_PTV{jj};
    comb_AH_PTV_var_full(:, start_col:end_col) = AH_PTV_var{jj};
    comb_AH_MultiPulse_full(:, start_col:end_col) = AH_MultiPulse{jj};
    comb_AH_model_full(:, start_col:end_col) = AH_model{jj};
    comb_T_full(:, start_col:end_col) = T{jj};
    comb_T_Std_full(:, start_col:end_col) = T_Std{jj}; 
    comb_T_model_full(:, start_col:end_col) = T_model{jj};
    comb_ABC_full(:, start_col:end_col) = ABC{jj};
    
    % Decimated Copies (Base variables)
    comb_AH(:, start_col:end_col) = AH{jj};
    comb_AH_var(:, start_col:end_col) = AH_var{jj};
    comb_AH_PTV(:, start_col:end_col) = AH_PTV{jj};
    comb_AH_PTV_var(:, start_col:end_col) = AH_PTV_var{jj};
    comb_AH_MultiPulse(:, start_col:end_col) = AH_MultiPulse{jj};
    comb_AH_MultiPulse_var(:, start_col:end_col) = AH_MultiPulse_var{jj};
    comb_AH_model(:, start_col:end_col) = AH_model{jj};
    comb_ABC(:, start_col:end_col) = ABC{jj};
    comb_ABC_var(:, start_col:end_col) = ABC_var{jj};
    comb_T(:, start_col:end_col) = T{jj};
    comb_T_var(:, start_col:end_col) = T_var{jj};
    comb_T_Std(:, start_col:end_col) = T_Std{jj}; 
    comb_T_Std_var(:, start_col:end_col) = T_Std_var{jj}; 
    comb_T_model(:, start_col:end_col) = T_model{jj};
    comb_Counts(:, start_col:end_col) = Counts{jj}; 
    
    start_col = end_col + 1;
end

disp(['Data combination complete. Time elapsed: ', num2str(toc(tic_combine)), ' seconds.']);

% --- NaN GAP INSERTION FOR VISUAL BREAKS (ORIGINAL LOGIC PRESERVED) ---
disp('Starting gap insertion...');
tic_gap = tic;

% Calculate the time difference (in days) between consecutive profiles
dt = diff(comb_duration);

% Gap threshold: If the time jump is greater than 2 days, assume it's a data gap
gap_threshold = 2.0; 

% Find indices immediately preceding a gap
gap_indices = find(dt > gap_threshold);

% Combine all variables for gap insertion
data_arrays_2d_all = {'comb_AH', 'comb_AH_var', 'comb_AH_PTV', 'comb_AH_PTV_var', 'comb_AH_MultiPulse', 'comb_AH_MultiPulse_var', 'comb_AH_model', 'comb_ABC', 'comb_ABC_var', 'comb_T', 'comb_T_var', 'comb_T_Std', 'comb_T_Std_var', 'comb_T_model', 'comb_Counts', ...
                  'comb_AH_full', 'comb_AH_var_full', 'comb_AH_PTV_full', 'comb_AH_PTV_var_full', 'comb_AH_MultiPulse_full', 'comb_AH_model_full', 'comb_ABC_full', 'comb_T_full', 'comb_T_Std_full', 'comb_T_model_full'};

data_arrays_1d = {'comb_duration', 'comb_T_surf', 'comb_P_surf', 'comb_AH_surf'};

if ~isempty(gap_indices)
    disp(['Found ', num2str(length(gap_indices)), ' time gaps exceeding ', num2str(gap_threshold), ' days. Inserting 2 NaNs...']);
    
    % Loop through the detected gaps in reverse order to maintain correct indexing
    for i = length(gap_indices):-1:1
        idx = gap_indices(i); % Index *before* the gap starts
        
        gap_start_time = comb_duration(idx);
        gap_end_time = comb_duration(idx + 1);
        
        % Define two insertion time points to span the gap visually:
        nan_time_1 = gap_start_time + (gap_end_time - gap_start_time) * 0.001;
        nan_time_2 = gap_end_time - (gap_end_time - gap_start_time) * 0.001;

        % --- Insert two NaN rows into 1D arrays (duration and surface data) ---
        for arr = data_arrays_1d
            array_name = arr{1};
            current_array = eval(array_name);
            
            % Set insertion values (Time for duration, NaN for surface data)
            nan_val_1 = NaN;
            nan_val_2 = NaN;
            if strcmp(array_name, 'comb_duration')
                nan_val_1 = nan_time_1; 
                nan_val_2 = nan_time_2; 
            end
            
            % Insert the two new rows/values at index idx + 1 and idx + 2
            current_array = [current_array(1:idx); nan_val_1; nan_val_2; current_array(idx+1:end)];
            eval([array_name ' = current_array;']);
        end
        
        % --- Insert two NaN columns into 2D arrays ---
        nan_cols = nan(num_ranges, 2);
        for arr = data_arrays_2d_all
            array_name = arr{1};
            % Check if variable exists before trying to evaluate/modify
            if exist(array_name, 'var')
                current_array = eval(array_name);
                % Insert the two NaN columns at index idx + 1
                current_array = [current_array(:, 1:idx), nan_cols, current_array(:, idx+1:end)];
                eval([array_name ' = current_array;']);
            end
        end
    end
else
    disp('No significant time gaps found.');
end

disp(['Gap insertion complete. Time elapsed: ', num2str(toc(tic_gap)), ' seconds.']);
 
% --- APPLY SMART DECIMATION HERE (ORIGINAL LOGIC PRESERVED) ---
disp('Starting smart decimation...');
tic_decimate = tic;

% 1. DYNAMICALLY CALCULATE N_DECIMATE
N_time_actual = length(comb_duration);
% Calculate factor, ensuring minimum decimation of 1
N_DECIMATE_DYNAMIC = max(1, floor(N_time_actual / TARGET_PIXEL_WIDTH));
disp(['Calculated dynamic decimation factor N = ', num2str(N_DECIMATE_DYNAMIC), '.']);

% List of DECIMATED 2D matrices (only comb_* variables are used for pcolor plots)
decimate_vars = {'comb_AH', 'comb_AH_var', 'comb_AH_PTV', 'comb_AH_PTV_var', 'comb_AH_MultiPulse', 'comb_AH_MultiPulse_var', 'comb_AH_model', 'comb_ABC', 'comb_ABC_var', 'comb_T', 'comb_T_var', 'comb_T_Std', 'comb_T_Std_var', 'comb_T_model', 'comb_Counts'};

% Decimate each 2D matrix
for i = 1:length(decimate_vars)
    var_name = decimate_vars{i};
    current_matrix = eval(var_name);
    
    if ~isempty(current_matrix)
        % Overwrite the full-resolution variable with the decimated version
        decimated_matrix = decimate_matrix(current_matrix, N_DECIMATE_DYNAMIC);
        eval([var_name ' = decimated_matrix;']);
    end
end

% Decimate the time axis (x-axis) to match the new DECIMATED data
N_total = length(comb_duration);
comb_duration_decimated = comb_duration(1:N_DECIMATE_DYNAMIC:N_total);

% Overwrite the old x-axis variable with the decimated version
x = comb_duration_decimated;

disp(['Decimation complete. Time elapsed: ', num2str(toc(tic_decimate)), ' seconds.']);


% --- 6. PLOT SETUP AND CALLS (MODIFIED FOR SCROLLING) ---
disp('Starting plotting and saving...');
tic_plot_save = tic;

scrsz = get(0,'ScreenSize');
date_str = datestr(comb_duration(1), 'yyyy-mmm-dd'); % Use first day of combined data
plot_size_wide = [scrsz(4)/1.5 scrsz(4)/10 scrsz(3)/1.5 scrsz(4)/3]; % Used for single plots
plot_size_5panel = [100 100 1200 800]; % Custom size for 5 vertical panels
plot_size_7panel = [100 100 1200 1100]; % Custom size for 7 vertical panels

font_size = 16; % ORIGINAL FONT SIZE FOR ALL SINGLE PANEL PLOTS
font_size_small = 8; % REDUCED FONT SIZE FOR MULTI-PANEL PLOTS
y = (alt{1}./1000); % Assumes 'alt' is consistent across files

% --- DEFINE GLOBAL CONSTANTS FOR SCROLLING ---
GLOBAL_X_DATA = x; % Decimated time array
GLOBAL_Y_DATA = y; % Altitude data
DAYS_TO_DISPLAY = 7; % Scrolling window size (1 week)
MIN_DATE = floor(min(GLOBAL_X_DATA));
MAX_DATE = ceil(max(GLOBAL_X_DATA));
current_start_date = MIN_DATE; % Initial view start date

% Define the smooth_gray_cmap (Differential plots)
N = 64;
N_half = N / 2;
blue = [0 0 1];
red = [1 0 0];
light_gray = [0.9 0.9 0.9];
R1 = linspace(blue(1), light_gray(1), N_half)'; G1 = linspace(blue(2), light_gray(2), N_half)'; B1 = linspace(blue(3), light_gray(3), N_half)';
cmap_half1 = [R1, G1, B1];
R2 = linspace(light_gray(1), red(1), N_half)'; G2 = linspace(light_gray(2), red(2), N_half)'; B2 = linspace(light_gray(3), red(3), N_half)';
cmap_half2 = [R2, G2, B2];
smooth_gray_cmap = [cmap_half1; cmap_half2];

% Define PLASMAP_CMAP as a matrix, necessary for stable subplot colormaps
plasma_cmap = plasma(64); % 

% Define common parameters for plots
cmap_AH = CM_YlGnBu(64); % 
caxis_AH = [0 25];
caxis_diff_AH = [-10 10];
caxis_diff_T = [-5 5];
caxis_uncertainty = [-5 5]; 
caxis_T_abs = [265 305];
caxis_ABC = [5e-9 1e-6];
caxis_Counts = [1e2 1e6];

% === HISTOGRAM/DENSITY MAP PARAMETER DEFINITIONS (FIXED LOCATION) ===
AH_diff_bin_limits = [-10 10];
T_diff_bin_limits = [-10 10];
T_hist_limits = [-15 15];
AH_hist_limits = [-10 10];
T_diff_bin_width = 0.2;
AH_diff_bin_width = 0.2;
cmap_density = flipud(gray(256));
% ====================================================================

figure_idx = 1;
figure_list = {}; 
suffix_list = {}; 
SCROLLABLE_FIGURES = []; % Track figures that are scrollable

% ---------------------------------------------
% --- FIGURE 1: 5-PANEL TEMPERATURE COMPARISON (SCROLLABLE) ---
% ---------------------------------------------
if flag.plot_multi_temp
    
    current_end_date = min(current_start_date + DAYS_TO_DISPLAY, MAX_DATE); 

    % --- 6.1. SETUP INTERACTIVE FIGURE AND CONTROLS ---
    figure_T = figure(figure_idx); % Figure 1
    set(figure_T, 'Position', plot_size_5panel, 'renderer', 'zbuffer', 'Name', [node, ': Temperature Comparison (Scrollable)']);
    
    % Store all required data structures for the callback function
    setappdata(figure_T, 'node', node);
    setappdata(figure_T, 'x', GLOBAL_X_DATA);
    setappdata(figure_T, 'y', GLOBAL_Y_DATA);
    
    % Store 2D data matrices (comb_T* variables)
    comp_T_data = {comb_T_model, comb_T, comb_T_Std, comb_T - comb_T_model, comb_T_Std - comb_T_model};
    setappdata(figure_T, 'comp_T_data', comp_T_data);
    setappdata(figure_T, 'flag_read_std', flag.read_temp_std);
    setappdata(figure_T, 'plot_params', struct('caxis_T_abs', caxis_T_abs, 'caxis_diff_T', caxis_diff_T, 'plasma_cmap', plasma_cmap, 'smooth_gray_cmap', smooth_gray_cmap, 'font_size_small', font_size_small));
    setappdata(figure_T, 'DAYS_TO_DISPLAY', DAYS_TO_DISPLAY);
    setappdata(figure_T, 'MIN_DATE', MIN_DATE);
    setappdata(figure_T, 'MAX_DATE', MAX_DATE);
    setappdata(figure_T, 'current_start_date', current_start_date);
    
    % --- 6.2. INITIAL PLOT CALL ---
    plot_temperature_comparison(figure_T, current_start_date, current_end_date);
    
    SCROLLABLE_FIGURES = [SCROLLABLE_FIGURES figure_idx];
    figure_list{end+1} = figure_idx;
    suffix_list{end+1} = 'T_5Panel_Comparison';
    figure_idx=figure_idx+1; % Increment figure index
end


% ---------------------------------------------------
% --- FIGURE 2: 7-PANEL WATER VAPOR COMPARISON (NOT YET SCROLLABLE, REMAINS BATCH PLOT) ---
% ---------------------------------------------------
if flag.plot_multi_wv
    hf = figure(figure_idx); % Figure 2 (or next available index)
    set(hf, 'Position', plot_size_7panel, 'renderer', 'zbuffer'); % Use custom 7-panel size

    % --- Dynamic Tick Calculation (Original logic) ---
    TARGET_TICKS = 8; 
    time_span_days = ceil(max(x) - min(x));
    days_per_tick = max(1, ceil(time_span_days / TARGET_TICKS));
    if days_per_tick <= 2, days_per_tick = 1; 
    elseif days_per_tick <= 4, days_per_tick = 3;
    elseif days_per_tick <= 10, days_per_tick = 7;
    elseif days_per_tick <= 20, days_per_tick = 14;
    elseif days_per_tick <= 60, days_per_tick = 30;
    else days_per_tick = ceil(days_per_tick / 30) * 30; end
    start_date_tick = floor(min(x));
    end_date_tick = ceil(max(x)); 
    xData = start_date_tick : days_per_tick : end_date_tick;
    xData_m = start_date_tick : days_per_tick/2 : end_date_tick;
    xData(xData > end_date_tick) = [];
    xData_m(xData_m > end_date_tick) = [];
    
    % --- Plotting logic (Original Figure 2 logic) ---
    
    num_panels = 7;
    total_gap = 0.05 * (num_panels - 1);
    total_height = 1 - 0.15 - total_gap;
    panel_height = total_height / num_panels;

    wv_plots = {
        comb_AH_model, 'ERA5 Model (g m^{-3})', caxis_AH, cmap_AH;
        comb_AH, 'Standard (g m^{-3})', caxis_AH, cmap_AH;
        comb_AH_PTV, 'PTV (g m^{-3})', caxis_AH, cmap_AH;
        comb_AH_MultiPulse, 'MultiPulse (g m^{-3})', caxis_AH, cmap_AH;
        comb_AH - comb_AH_model, 'Standard - Model Difference (g m^{-3})', caxis_diff_AH, smooth_gray_cmap;
        comb_AH_PTV - comb_AH_model, 'PTV - Model Difference (g m^{-3})', caxis_diff_AH, smooth_gray_cmap;
        comb_AH_MultiPulse - comb_AH_model, 'MultiPulse - Model Difference (g m^{-3})', caxis_diff_AH, smooth_gray_cmap;
    };
    
    if ~flag.read_wv_multi
        wv_plots([4 7], :) = []; 
        num_panels = size(wv_plots, 1);
        total_gap = 0.05 * (num_panels - 1);
        total_height = 1 - 0.15 - total_gap;
        panel_height = total_height / num_panels;
    end

    sgtitle({[node, ': Absolute Humidity Comparison']}, 'fontweight','b','fontsize',font_size);
    
    SUBPLOT_LEFT = 0.1;
    SUBPLOT_WIDTH = 0.78; 
    
    for i = 1:num_panels
        bottom_pos = 0.1 + (num_panels - i) * (panel_height + 0.05);

        h_ax = subplot('Position', [SUBPLOT_LEFT, bottom_pos, SUBPLOT_WIDTH, panel_height]);
        
        data_plot = real(wv_plots{i, 1});
        plot_title = wv_plots{i, 2};
        caxis_val = wv_plots{i, 3};
        cmap_val = wv_plots{i, 4}; 
        
        h = pcolor(h_ax, x, y, data_plot); 
        set(h, 'EdgeColor', 'none'); 
        axis xy; 
        
        h_cb = colorbar('peer', h_ax, 'EastOutside');
        
        is_abs_plot = contains(plot_title, '(g m^{-3})') && ~contains(plot_title, 'Difference');
        
        if is_abs_plot
            ylabel(h_cb, 'Absolute Humidity (g m^{-3})', 'fontweight', 'b', 'fontsize', font_size_small - 2); 
        else
            ylabel(h_cb, 'AH Difference (g m^{-3})', 'fontweight', 'b', 'fontsize', font_size_small - 2); 
        end

        axis([floor(min(x)), ceil(max(x)), 0, 6]);
        caxis(caxis_val);
        h_ax.Colormap = cmap_val;
        
        set(gca,'Fontsize',font_size_small,'Fontweight','b'); 
        title({plot_title}, 'fontweight','b','fontsize',font_size_small);  
        ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size_small); 
        
        set(gca, 'XTick', xData);
        set(gca,'XMinorTick','on');
        xAx = get(gca,'XAxis');
        xAx.MinorTickValues=xData_m;
        datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
        
        if i ~= num_panels 
            xlabel('');
        else
            xlabel('Time (UTC)','fontweight','b','fontsize',font_size_small); 
        end
        
        set(gca,'TickDir','out');
        set(gca,'TickLength',[0.005; 0.0025]);
    end
    figure_list{end+1} = figure_idx;
    suffix_list{end+1} = 'WV_7Panel_Comparison';
    figure_idx=figure_idx+1; 
end

% ... (Original logic for single-panel plots, Aerosol/Counts, Uncertainty, and Histograms follows)
% NOTE: The batch plots use the original dynamic tick calculation logic defined above Figure 2.

% AEROSOL/COUNTS (Conditional)
if flag.plot_aerosol_counts
    % Figure X
    create_pcolor_plot(x, y, comb_ABC, 'Aerosol Backscatter Coefficient, PTV (m^{-1} sr^{-1})', caxis_ABC, 'viridis', node, plot_size_wide, font_size, xData, xData_m, figure_idx, 1);
    figure_list{end+1} = figure_idx; suffix_list{end+1} = 'Back_Coeff_PTV_multi'; figure_idx=figure_idx+1;
    
    % Figure X+1
    create_pcolor_plot(x, y, comb_Counts, 'Attenuated Backscatter, 828 nm', caxis_Counts, 'jet', node, plot_size_wide, font_size, xData, xData_m, figure_idx, 1); 
    figure_list{end+1} = figure_idx; suffix_list{end+1} = '828_Counts_multi'; figure_idx=figure_idx+1;
end


% --- PLOT SAVING LOOP (MODIFIED TO SKIP SCROLLABLE FIGURES) ---
cd(plot_dir);
plot_size = [100 100 1920 225];
plot_size_square = [100 100 750 750]; 
plot_size_5panel = [100 100 1900 800]; 
plot_size_7panel = [100 100 1900 1000]; 
num_plots = length(figure_list);

for k_idx = 1:num_plots
    k = figure_list{k_idx}; 
    suffix = suffix_list{k_idx}; 
    
    % Skip scrollable figures from batch saving
    if ismember(k, SCROLLABLE_FIGURES)
        disp(['Skipping interactive figure ', num2str(k), ' (', suffix, ') from batch save.']);
        continue; 
    end
    
    FigH = figure(k);
    drawnow;
    FigH.Units = 'pixels';
    
    % Determine size based on the specific figure created
    if strcmp(suffix, 'T_5Panel_Comparison')
        FigH.Position = plot_size_5panel;
    elseif strcmp(suffix, 'WV_7Panel_Comparison')
        FigH.Position = plot_size_7panel;
    % Check if it's a square plot (histogram/density maps)
    elseif contains(suffix, 'Histogram_Range') || contains(suffix, 'Histogram')
        FigH.Position = plot_size_square;
    else
        % Default wide size for single pcolor plots
        FigH.Position = plot_size; 
    end
    
    name = char(strcat(node, "_", date_str, '_', suffix));
    exportgraphics(FigH, [name, '.png'], 'Resolution', 150);
    
    close(FigH); 
end

disp(['Plotting and saving complete. Time elapsed: ', num2str(toc(tic_plot_save)), ' seconds.']);
disp(['Total script runtime: ', num2str(toc(tic_total)), ' seconds.']);
disp('--- Script End ---');

% ------------------------------------------------------------------
% --- LOCAL FUNCTIONS (SCROLLING IMPLEMENTATION) ---
% ------------------------------------------------------------------

function scroll_plot_T(hObject, eventdata, direction)
    % Callback function for the scroll buttons
    
    hf = gcf; % Get the current figure handle
    DAYS_TO_DISPLAY = getappdata(hf, 'DAYS_TO_DISPLAY');
    MIN_DATE = getappdata(hf, 'MIN_DATE');
    MAX_DATE = getappdata(hf, 'MAX_DATE');
    current_start_date = getappdata(hf, 'current_start_date');

    step = DAYS_TO_DISPLAY;
    
    if strcmp(direction, 'next')
        new_start_date = current_start_date + step;
        % Stop scrolling once the window passes the max date range
        if new_start_date >= MAX_DATE
             % Check if the last full window can be shown, otherwise show the last available data span
             last_possible_start = MAX_DATE - DAYS_TO_DISPLAY;
             new_start_date = max(current_start_date, last_possible_start);
             if new_start_date < MIN_DATE, new_start_date = MIN_DATE; end % Don't go before min date
        end
    elseif strcmp(direction, 'prev')
        new_start_date = current_start_date - step;
        % Stop scrolling before the min date
        if new_start_date < MIN_DATE
            new_start_date = MIN_DATE; 
        end
    end
    
    new_end_date = min(new_start_date + DAYS_TO_DISPLAY, MAX_DATE);
    
    % Update the stored start date
    setappdata(hf, 'current_start_date', new_start_date);
    
    % Redraw the plot with the new limits
    plot_temperature_comparison(hf, new_start_date, new_end_date);
end


function plot_temperature_comparison(hf, start_date, end_date)
    % Function to plot the 5-panel temperature comparison given time limits.
    
    % Retrieve stored data
    node = getappdata(hf, 'node');
    x = getappdata(hf, 'x');
    y = getappdata(hf, 'y');
    comp_T_data = getappdata(hf, 'comp_T_data');
    flag_read_std = getappdata(hf, 'flag_read_std');
    params = getappdata(hf, 'plot_params');
    
    % Use new dynamic tick calculation for the current visible span
    [xData, xData_m] = calculate_dynamic_ticks(start_date, end_date);
    
    % --- Define Plotting Scheme ---
    temp_plots = {
        comp_T_data{1}, 'ERA5 Model (K)', params.caxis_T_abs, params.plasma_cmap; 
        comp_T_data{2}, 'PTV (K)', params.caxis_T_abs, params.plasma_cmap;           
        comp_T_data{3}, 'Standard (K)', params.caxis_T_abs, params.plasma_cmap;    
        comp_T_data{4}, 'PTV - Model Difference (K)', params.caxis_diff_T, params.smooth_gray_cmap; 
        comp_T_data{5}, 'Standard - Model Difference (K)', params.caxis_diff_T, params.smooth_gray_cmap; 
    };
    
    % Filter out Standard plots if data was excluded
    if ~flag_read_std
        temp_plots([3 5], :) = []; % Remove Standard and Standard Diff rows
    end
    num_panels = size(temp_plots, 1);

    % Clear previous axes and set supertitle
    clf(hf); 
    sgtitle({[node, ': Temperature Comparison (', datestr(start_date, 'dd-mmm-yy'), ' to ', datestr(end_date, 'dd-mmm-yy'), ')']}, 'fontweight','b','fontsize',params.font_size_small * 2); 

    SUBPLOT_LEFT = 0.1;
    SUBPLOT_WIDTH = 0.78;
    total_gap = 0.05 * (num_panels - 1); 
    total_height = 1 - 0.15 - total_gap; 
    panel_height = total_height / num_panels;

    for i = 1:num_panels
        bottom_pos = 0.1 + (num_panels - i) * (panel_height + 0.05);
        
        % 1. Create subplot using explicit position
        h_ax = subplot('Position', [SUBPLOT_LEFT, bottom_pos, SUBPLOT_WIDTH, panel_height]);
        
        data_plot = real(temp_plots{i, 1});
        plot_title = temp_plots{i, 2};
        caxis_val = temp_plots{i, 3};
        cmap_val = temp_plots{i, 4};
        
        % 2. Pcolor plot (Core logic from create_pcolor_plot)
        h = pcolor(h_ax, x, y, data_plot); 
        set(h, 'EdgeColor', 'none'); 
        axis xy; 
        
        % Add colorbar and label it
        h_cb = colorbar('peer', h_ax, 'EastOutside');
        is_abs_plot = contains(plot_title, '(K)');
        if is_abs_plot
            ylabel(h_cb, 'Temperature (K)', 'fontweight', 'b', 'fontsize', params.font_size_small - 2); 
        else
            ylabel(h_cb, 'Temp. Difference (K)', 'fontweight', 'b', 'fontsize', params.font_size_small - 2); 
        end

        % 3. Set common axis limits and colormap
        axis([start_date, end_date, 0, 6]); % **CRITICAL: Use the scroll limits**
        caxis(caxis_val);
        h_ax.Colormap = cmap_val;
        
        % Apply SMALLER FONT SIZE
        set(gca,'Fontsize',params.font_size_small,'Fontweight','b'); 
        
        % 4. Titles and Labels
        title({plot_title}, 'fontweight','b','fontsize',params.font_size_small);  
        ylabel('Height (km, AGL)','fontweight','b','fontsize',params.font_size_small); 
        
        % 5. X-axis formatting: Use the dynamically calculated ticks
        set(gca, 'XTick', xData);
        set(gca,'XMinorTick','on');
        xAx = get(gca,'XAxis');
        xAx.MinorTickValues=xData_m;
        datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
        
        if i ~= num_panels 
            xlabel('');
        else
            % Add X-label to the bottom panel
            xlabel('Time (UTC)','fontweight','b','fontsize',params.font_size_small); 
        end
        
        set(gca,'TickDir','out');
        set(gca,'TickLength',[0.005; 0.0025]);
    end
    
    % Redraw buttons on top of the cleared figure
    uicontrol('Style', 'pushbutton', 'String', '<< Previous Week', ...
              'Units', 'normalized', 'Position', [0.05, 0.01, 0.15, 0.04], ...
              'Callback', {@scroll_plot_T, 'prev'});
    uicontrol('Style', 'pushbutton', 'String', 'Next Week >>', ...
              'Units', 'normalized', 'Position', [0.80, 0.01, 0.15, 0.04], ...
              'Callback', {@scroll_plot_T, 'next'});
end


function [xData, xData_m] = calculate_dynamic_ticks(start_date, end_date)
    % Helper function to generate sensible ticks for the current plot window.
    
    time_span_days = end_date - start_date;
    TARGET_TICKS = 7; 
    
    days_per_tick = max(1, floor(time_span_days / TARGET_TICKS));
    
    % Ensure days_per_tick is a readable interval
    if days_per_tick <= 2
        days_per_tick = 1; % Daily
    elseif days_per_tick <= 4
        days_per_tick = 3; % Every 3 days
    elseif days_per_tick <= 10
        days_per_tick = 7; % Weekly
    elseif days_per_tick <= 20
        days_per_tick = 14; % Bi-weekly
    elseif days_per_tick <= 60
        days_per_tick = 30; % Monthly
    else
        days_per_tick = ceil(days_per_tick / 30) * 30; % Round up to the nearest month
    end
    
    % Generate major ticks starting at the beginning of the first day visible
    xData = floor(start_date) : days_per_tick : ceil(end_date);
    % Generate minor ticks at half the major interval
    xData_m = floor(start_date) : days_per_tick/2 : ceil(end_date);

    % Remove ticks that fall outside the current display range
    xData(xData < start_date | xData > end_date) = [];
    xData_m(xData_m < start_date | xData_m > end_date) = [];
    
    % Add the start and end dates as ticks if they aren't already close to an existing tick
    xData = unique([start_date, xData, end_date]);
    xData_m = unique([start_date, xData_m, end_date]);
    
    % Ensure no duplicates and sort
    xData = sort(unique(xData));
    xData_m = sort(unique(xData_m));
end

% ------------------------------------------------------------------
% --- ORIGINAL HELPER FUNCTIONS (Preserved in full) ---
% ------------------------------------------------------------------

function hf = create_pcolor_plot(x, y, data, plot_title, caxis_val, cmap, node, plot_size1, font_size, xData, xData_m, figure_num, log_scale)
% Function to create and format a pcolor plot with common settings.
% ... (Original help text removed for brevity)

    Z = real(data);
    
    % Create figure
    hf = figure(figure_num);
    set(hf, 'Position', plot_size1, 'renderer', 'zbuffer');
    
    % Pcolor plot
    h = pcolor(x, y, Z);
    set(h, 'EdgeColor', 'none'); 
    axis xy; 
    colorbar('EastOutside'); 
    
    % Set common axis limits
    axis([fix(min(x)) ceil(max(x)) 0 6]);
    caxis(caxis_val);
    
    % Handle log scale for backscatter/counts plots
    if log_scale
        set(gca, 'Colorscale', 'log');
    else
        set(gca, 'Colorscale', 'linear'); % Ensure it's linear if not log
    end
    
    % Axis formatting
    set(gca, 'XTick',  xData)
    set(gca,'XMinorTick','on')
    xAx = get(gca,'XAxis');
    % The xData_m setting is only needed if xAx is defined. 
    % We'll only apply it if xData_m is actually passed and needed for minor ticks
    if nargin > 10 && ~isempty(xData_m)
        xAx.MinorTickValues=xData_m;
    end
    set(gca,'TickDir','out');
    set(gca,'TickLength',[0.005; 0.0025]);
    
    % Titles and Labels
    plot_title = {sprintf('%s %s', node, plot_title)};
    title(plot_title, 'fontweight','b','fontsize',font_size);  
    ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
    
    % Time axis format
    datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
    
    % Colormap and Font
    colormap(cmap);
    set(gca,'Fontsize',font_size,'Fontweight','b');
end

function create_2d_density_map(diff_data, y_range, title_str, xlabel_str, x_bin_limits, y_bin_limits, x_bin_width, fig_num, node, plot_size1, font_size, cmap)
    figure(fig_num);
    set(gcf, 'Position', plot_size1);
    
    % Prepare data (Flatten and Filter NaN)
    num_timesteps = size(diff_data, 2);
    alt_matrix = repmat(y_range, 1, num_timesteps);
    alt_vector = alt_matrix(:);
    diff_vector = diff_data(:);
    valid_idx = ~isnan(diff_vector) & ~isnan(alt_vector);
    alt_valid = alt_vector(valid_idx);
    diff_valid = diff_vector(valid_idx);

    % Define bin parameters
    Alt_bin_width = 0.12; % Consistent vertical resolution
    
    % Create the 2D histogram
    h_hist2 = histogram2(diff_valid, alt_valid, ...
        'BinWidth', [x_bin_width Alt_bin_width], ... 
        'XBinLimits', x_bin_limits, ... 
        'YBinLimits', y_bin_limits, ...   
        'Normalization', 'probability', ...  
        'DisplayStyle', 'tile');            

    % --- Formatting ---
    title({[node, ' ', title_str]}, 'fontweight', 'b', 'fontsize', font_size, 'Interpreter', 'latex'); % <-- ADDED
    xlabel(xlabel_str, 'fontweight', 'b', 'fontsize', font_size, 'Interpreter', 'latex');              % <-- ADDED
    ylabel('Height (km, AGL)', 'fontweight', 'b', 'fontsize', font_size);
    grid on;
    set(gca, 'Layer', 'top'); 

    xlim(x_bin_limits); 
    ylim(y_bin_limits); 

    colormap(cmap); 
    colorbar('EastOutside'); 
    set(gca, 'Fontsize', font_size, 'Fontweight', 'b');
    box on;
end

function create_1d_histogram(data_vector, title_str, xlabel_str, bin_limits, fig_num, node, font_size)
    figure(fig_num);
    
    % Get screen size (Assuming scrsz is defined globally or passed)
    scrsz = get(0,'ScreenSize');
    set(gcf, 'Position', [scrsz(4)/1.5 scrsz(4)/10 scrsz(3)/3 scrsz(4)/3]);

    % Remove NaN values
    data_vector = data_vector(~isnan(data_vector));
    
    % Plot the histogram
    h_hist = histogram(data_vector, 'BinLimits', bin_limits, 'NumBins', 100, 'Normalization', 'probability');
    
    % --- Formatting ---
    title({[node, ' ', title_str]}, 'fontweight', 'b', 'fontsize', font_size);
    xlabel(xlabel_str, 'fontweight', 'b', 'fontsize', font_size);
    ylabel('Probability Density', 'fontweight', 'b', 'fontsize', font_size);
    grid on;
    box on;
    set(gca, 'Fontsize', font_size, 'Fontweight', 'b');
    
    % Calculate and display stats
    data_in_limits = data_vector(data_vector >= bin_limits(1) & data_vector <= bin_limits(2));
    bias = mean(data_in_limits);
    std_dev = std(data_in_limits);
    ylim_max = max(h_hist.Values); 
    
    unit = strsplit(xlabel_str, '(');
    unit = strtrim(unit{2}(1:end-1));
    text_str = sprintf('Bias: %.2f %s\nStd Dev: %.2f %s', bias, unit, std_dev, unit);
    text(max(xlim)*0.6, ylim_max*0.9, text_str, ...
         'Fontsize', font_size-2, 'FontWeight', 'b', 'BackgroundColor', 'w');
end

function output_matrix = decimate_matrix(input_matrix, N_decimate)
    % Decimates the time (column) dimension of a 2D matrix by averaging N_decimate points.
    
    [N_range, N_time] = size(input_matrix);
    
    % Pad the end of the matrix with NaNs if the size isn't evenly divisible
    N_remainder = mod(N_time, N_decimate);
    if N_remainder ~= 0
        N_pad = N_decimate - N_remainder;
        padding = nan(N_range, N_pad);
        input_matrix = [input_matrix, padding];
        N_time = size(input_matrix, 2);
    end
    
    % Reshape and average: (N_range) x (N_time/N_decimate)
    output_matrix = mean(reshape(input_matrix, N_range, N_decimate, N_time/N_decimate), 2);
    output_matrix = squeeze(output_matrix);
end