
clear all; close all

% --- 1. CONFIGURATION AND PATHS ---
node = 'MPD04';
% Use fullfile for robust path definition
serv_path = '/Volumes/eol/sci/mhayman';
plot_path = '/Users/spuler/Desktop';

data_dir = fullfile(serv_path, 'DIAL', 'Processed_Data', 'BRIDGE_2025', 'ptv0.3');
plot_dir = fullfile(plot_path, 'mpd', 'Plots');
addpath '/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing/matplotlib/';
addpath '/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing'

flag.save_data = 0;  %save data at end of processing (0=off 1=on)
low_range_mask = 0;

% --- PLOTTING EXCLUSION FLAGS (0=Exclude/Off, 1=Include/On) ---
flag.plot_multi_temp = 0;      % Toggle multi panel T WV abs (Figure 1)
flag.plot_multi_wv = 0;        % Toggle multi panel T WV diff (Figure 2)

% Flags for REDUNDANT SINGLE PANELS (Defaulted to 0/OFF)
flag.plot_single_temp_panels = 0; 
flag.plot_single_wv_panels = 0;   

% Flags for essential/less redundant panels
flag.plot_aerosol_counts = 0;  
flag.plot_uncertainty = 0;     
flag.plot_histograms_2d = 1;   % Toggle 2D Density Maps (T_Range, AH_Range)
flag.plot_ah_scatter = 0;      % Toggle AH vs ERA5 Scatter Plot
flag.plot_temp_scatter = 0;      % Toggle AH vs ERA5 Scatter Plot
flag.plot_ah_density_map = 1;  % Toggle AH vs ERA5 Density Map
flag.plot_t_density_map = 1;   % Toggle Temp vs ERA5 Density Map

flag.plot_1d_histograms = 0;   % Toggle 1D Histograms (T_Diff_Histogram, AH_Diff_Histogram)
% ------------------------------------

% --- DATA READING EXCLUSION FLAGS (0=Exclude/Off, 1=Include/On) ---
flag.read_temp_std = 1;        % Include reading Temperature_Standard and related fields
flag.read_wv_multi = 1;        % Include reading Absolute_Humidity_MultiPulse fields
flag.read_uncertainty = 0;     % Include reading ALL uncertainty/variance fields (T_var, AH_var, etc.)
% ------------------------------------------------------------------

% --- PLOTTING RESOLUTION SETTING ---
TARGET_PIXEL_WIDTH = 10000;     % Target output width (pixels) for visual quality (used to calculate N_DECIMATE)
flag.decimate_hist_data = 1;   % 0 = Use FULL RESOLUTION data for 2D Histograms (High Accuracy, Slower)
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

% --- 4. DATA READING AND MASKING LOOP (Conditional Reading) ---
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

% --- 5. PREALLOCATION AND COMBINATION ---
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

% --- NaN GAP INSERTION FOR VISUAL BREAKS ---
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
 
% --- APPLY SMART DECIMATION HERE ---
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


% --- 6. PLOT SETUP AND CALLS (Conditional Plotting Implemented) ---
disp('Starting plotting and saving...');
tic_plot_save = tic;

scrsz = get(0,'ScreenSize');
date_str = datestr(comb_duration(1), 'yyyy-mmm-dd'); % Use first day of combined data
plot_size_wide = [scrsz(4)/1.5 scrsz(4)/10 scrsz(3)/1.5 scrsz(4)/3]; % Used for single plots
plot_size_5panel = [100 100 1200 800]; % Custom size for 5 vertical panels (Figure 1 uses this)
plot_size_7panel = [100 100 1200 1100]; % Custom size for 7 vertical panels (Figure 2 uses this)

font_size = 16; % ORIGINAL FONT SIZE FOR ALL SINGLE PANEL PLOTS
font_size_small = 8; % REDUCED FONT SIZE FOR MULTI-PANEL PLOTS
y = (alt{1}./1000); % Assumes 'alt' is consistent across files

% --- DYNAMIC TICK MARK CALCULATION (FINAL FIX) ---
disp('Calculating dynamic tick interval...');
% Target 8 major ticks for the main time span
TARGET_TICKS = 8; 
time_span_days = ceil(max(comb_duration) - min(comb_duration));
days_per_tick = max(1, ceil(time_span_days / TARGET_TICKS));

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

% --- CRITICAL FIX: Use the calculated step size to generate an exact date series.
start_date = floor(min(x));
end_date = ceil(max(x)); 

% Generate major ticks starting at the beginning of the first day and stepping by days_per_tick
xData = start_date : days_per_tick : end_date;
% Generate minor ticks at half the major interval
xData_m = start_date : days_per_tick/2 : end_date;

% Remove ticks that fall outside the display range for cleanliness
xData(xData > end_date) = [];
xData_m(xData_m > end_date) = [];

disp(['Dynamic Major Tick Interval: ', num2str(days_per_tick), ' days.']);
% ------------------------------------------

% % Define the smooth_gray_cmap (Differential plots)
% N = 64;
% N_half = N / 2;
% blue = [0 0 1];
% red = [1 0 0];
% light_gray = [0.9 0.9 0.9];
% R1 = linspace(blue(1), light_gray(1), N_half)'; G1 = linspace(blue(2), light_gray(2), N_half)'; B1 = linspace(blue(3), light_gray(3), N_half)';
% cmap_half1 = [R1, G1, B1];
% R2 = linspace(light_gray(1), red(1), N_half)'; G2 = linspace(light_gray(2), red(2), N_half)'; B2 = linspace(light_gray(3), red(3), N_half)';
% cmap_half2 = [R2, G2, B2];
% smooth_gray_cmap = [cmap_half1; cmap_half2];


% % Call the Local Function to Generate the Colormap (N_colors, max_value, max_gray_value)---
smooth_gray_cmap_T  = create_diverging_cmap(64, 5, 1);
smooth_gray_cmap_AH = create_diverging_cmap(64, 5, 1);


% Define PLASMAP_CMAP as a matrix, necessary for stable subplot colormaps
plasma_cmap = plasma(64); % Assuming 'plasma' is available/defined

% Define common parameters for plots
cmap_AH = CM_YlGnBu(64); % CM_YlGnBu is a custom function
caxis_AH = [0 25];
caxis_diff_AH = [-5 5];
caxis_diff_T = [-5 5];
caxis_uncertainty = [-5 5]; 
caxis_T_abs = [265 305];
caxis_ABC = [5e-9 1e-6];
caxis_Counts = [1e2 1e6];

% === HISTOGRAM/DENSITY MAP PARAMETER DEFINITIONS (FIXED LOCATION) ===
AH_diff_bin_limits = [-10 10];
T_diff_bin_limits = [-10 10]; %
DENSITY_CLIM = [1e-4 1e-2]; % for ERA5Global Density Plot CLim ---
AH_density_limits = [0 25];  % Absolute Humidity (g/m3)
T_density_limits = [265 305]; % Temperature (K)
T_hist_limits = [-15 15];
AH_hist_limits = [-10 10];

T_diff_bin_width = 0.2;
AH_diff_bin_width = 0.2;

cmap_density = flipud(magma(256));
% ====================================================================

figure_idx = 1;
figure_list = {}; % Dynamic list to store figure numbers that were actually created
suffix_list = {}; % Dynamic list to store suffixes for figures that were actually created

% --- FIX: Normalized Plotting Coordinates for Subplot Alignment ---
SUBPLOT_LEFT = 0.1;
SUBPLOT_WIDTH = 0.78; 

% ---------------------------------------------
% --- FIGURE 1: 7-PANEL ABSOLUTE VALUES (AH & T) ---
% ---------------------------------------------
if flag.plot_multi_temp % Keeping the flag name, but changing its function
    hf = figure(figure_idx); % Figure 1
    % Ensure this figure uses the 7-panel size
    set(hf, 'Position', plot_size_7panel, 'renderer', 'zbuffer');
    
    num_panels_base = 7;
    
    % Initial list of all 7 plots (Absolute Values, Top to Bottom)
    all_abs_plots = {
        % ABSOLUTE HUMIDITY
        comb_AH_model, 'Absolute Humidity ERA5 Model (g m^{-3})', caxis_AH, cmap_AH;
        comb_AH, 'Absolute Humidity Standard (g m^{-3})', caxis_AH, cmap_AH;
        comb_AH_PTV, 'Absolute Humidity PTV (g m^{-3})', caxis_AH, cmap_AH;
        comb_AH_MultiPulse, 'Absolute Humidity MultiPulse (g m^{-3})', caxis_AH, cmap_AH;
        % TEMPERATURE
        comb_T_model, 'Temperature ERA5 Model (K)', caxis_T_abs, plasma_cmap;
        comb_T_Std, 'Temperature Standard (K)', caxis_T_abs, plasma_cmap;
        comb_T, 'Temperature PTV (K)', caxis_T_abs, plasma_cmap;
    };
    
    % --- CONDITIONAL PANEL REMOVAL ---
    plots_to_remove = [];
    if ~flag.read_wv_multi
        plots_to_remove = [plots_to_remove, 4]; % Remove MultiPulse (Row 4)
    end
    if ~flag.read_temp_std
        plots_to_remove = [plots_to_remove, 6]; % Remove T_Standard (Row 6)
    end

    % Apply removal
    if ~isempty(plots_to_remove)
        all_abs_plots(plots_to_remove, :) = [];
    end
    
    num_panels = size(all_abs_plots, 1);
    Vertical_gap = 0.035;
    Total_Vertical_Margin = 0.1;

    % Recalculate panel geometry based on actual number of panels
    total_gap = Vertical_gap * (num_panels - 1);
    total_height = 1 - Total_Vertical_Margin - total_gap;
    panel_height = total_height / num_panels;

    % Loop through all_abs_plots (formerly temp_plots)
    for i = 1:num_panels
        % Calculate normalized bottom position for panel i
        bottom_pos = Total_Vertical_Margin/2 + (num_panels - i) * (panel_height + Vertical_gap);
        
        % 1. Create subplot using explicit position
        h_ax = subplot('Position', [SUBPLOT_LEFT, bottom_pos, SUBPLOT_WIDTH, panel_height]);
        
        data_plot = real(all_abs_plots{i, 1});
        plot_title = all_abs_plots{i, 2};
        caxis_val = all_abs_plots{i, 3};
        cmap_val = all_abs_plots{i, 4};
        
        % 2. Pcolor plot (Core logic from create_pcolor_plot)
        h = pcolor(h_ax, x, y, data_plot); 
        set(h, 'EdgeColor', 'none'); 
        axis xy; 
        
        % Add colorbar
        h_cb = colorbar('peer', h_ax, 'EastOutside');
        
        % Set dynamic colorbar label based on contents
        if contains(plot_title, 'Humidity')
            ylabel(h_cb, 'Absolute Humidity (g m^{-3})', 'fontweight', 'b', 'fontsize', font_size_small - 2); 
        else % Temperature plots
            ylabel(h_cb, 'Temperature (K)', 'fontweight', 'b', 'fontsize', font_size_small - 2); 
        end

        % 3. Set common axis limits and colormap
        axis([floor(min(x)), ceil(max(x)), 0, 6]);
        caxis(caxis_val);
        h_ax.Colormap = cmap_val;
        
        % Apply SMALLER FONT SIZE
        set(gca,'Fontsize',font_size_small,'Fontweight','b'); 
        
        % 4. Titles and Labels
        title({plot_title}, 'fontweight','b','fontsize',font_size_small);  
        ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size_small); 
        
        % 5. X-axis formatting: Plot labels on ALL panels for guaranteed width
        set(gca, 'XTick', xData);
        set(gca,'XMinorTick','on');
        xAx = get(gca,'XAxis');
        xAx.MinorTickValues=xData_m;
        datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
        
        % Only show X-label on the bottom panel
        if i ~= num_panels 
            xlabel('');
            %set(gca, 'XTickLabel', []); % <-- Original line (keep commented out as before)
        else
            % Add X-label to the bottom panel
            xlabel('Time (UTC)','fontweight','b','fontsize',font_size_small); 
        end
        
        set(gca,'TickDir','out');
        set(gca,'TickLength',[0.005; 0.0025]);
    end
    figure_list{end+1} = figure_idx;
    suffix_list{end+1} = 'T_WV_7Panel_Absolute'; % NEW SUFFIX
    figure_idx=figure_idx+1; % Increment figure index
end


% ---------------------------------------------------
% --- FIGURE 2: 7-PANEL DIFFERENCES (AH & T) ---
% ---------------------------------------------------
if flag.plot_multi_wv % Keeping the flag name, but changing its function
    hf = figure(figure_idx); % Figure 2 (or next available index)
    % Ensure this figure uses the 7-panel size
    set(hf, 'Position', plot_size_7panel, 'renderer', 'zbuffer'); % Use custom 7-panel size

    num_panels_base = 7;
    
    % Initial list of all 7 plots (Model and Differences, Top to Bottom)
    all_diff_plots = {
        % ABSOLUTE HUMIDITY
        comb_AH_model, 'Absolute Humidity ERA5 Model (g m^{-3})', caxis_AH, cmap_AH; % Model Plot
        comb_AH - comb_AH_model, 'AH Standard - Model Difference (g m^{-3})', caxis_diff_AH, smooth_gray_cmap_AH;
        comb_AH_PTV - comb_AH_model, 'AH PTV - Model Difference (g m^{-3})', caxis_diff_AH, smooth_gray_cmap_AH;
        comb_AH_MultiPulse - comb_AH_model, 'AH MultiPulse - Model Difference (g m^{-3})', caxis_diff_AH, smooth_gray_cmap_AH;
        % TEMPERATURE
        comb_T_model, 'Temperature ERA5 Model (K)', caxis_T_abs, plasma_cmap; % Model Plot
        comb_T_Std - comb_T_model, 'T Standard - Model Difference (K)', caxis_diff_T, smooth_gray_cmap_T;
        comb_T - comb_T_model, 'T PTV - Model Difference (K)', caxis_diff_T, smooth_gray_cmap_T;
    };
    
    % --- CONDITIONAL PANEL REMOVAL ---
    plots_to_remove = [];
    if ~flag.read_wv_multi
        plots_to_remove = [plots_to_remove, 4]; % Remove MultiPulse Diff (Row 4)
    end
    if ~flag.read_temp_std
        plots_to_remove = [plots_to_remove, 6]; % Remove T_Standard Diff (Row 6)
    end
    
    % Apply removal
    if ~isempty(plots_to_remove)
        all_diff_plots(plots_to_remove, :) = [];
    end

    num_panels = size(all_diff_plots, 1);
    Vertical_gap = 0.035;
    Total_Vertical_Margin = 0.1;

    % Recalculate panel geometry based on actual number of panels
    total_gap = Vertical_gap * (num_panels - 1);
    total_height = 1 - Total_Vertical_Margin - total_gap;
    panel_height = total_height / num_panels;

    
    for i = 1:num_panels
        % Calculate normalized bottom position for panel i
        bottom_pos = Total_Vertical_Margin/2  + (num_panels - i) * (panel_height + Vertical_gap);

        % 1. Create subplot using explicit position
        h_ax = subplot('Position', [SUBPLOT_LEFT, bottom_pos, SUBPLOT_WIDTH, panel_height]);
        
        data_plot = real(all_diff_plots{i, 1});
        plot_title = all_diff_plots{i, 2};
        caxis_val = all_diff_plots{i, 3};
        cmap_val = all_diff_plots{i, 4}; 
        
        % 2. Pcolor plot
        h = pcolor(h_ax, x, y, data_plot); 
        set(h, 'EdgeColor', 'none'); 
        axis xy; 
        
        % 3. Add colorbar
        h_cb = colorbar('peer', h_ax, 'EastOutside');
        
        % Check if this is an absolute plot (model) or difference plot
        is_diff_plot = contains(plot_title, 'Difference');
        
        if is_diff_plot
            if contains(plot_title, 'AH')
                ylabel(h_cb, 'AH Difference (g m^{-3})', 'fontweight', 'b', 'fontsize', font_size_small - 2); 
            else % Temperature Difference
                ylabel(h_cb, 'Temp. Difference (K)', 'fontweight', 'b', 'fontsize', font_size_small - 2); 
            end
        else % Absolute Plot (Model)
            if contains(plot_title, 'Humidity')
                ylabel(h_cb, 'Absolute Humidity (g m^{-3})', 'fontweight', 'b', 'fontsize', font_size_small - 2); 
            else % Temperature
                ylabel(h_cb, 'Temperature (K)', 'fontweight', 'b', 'fontsize', font_size_small - 2); 
            end
        end


        % 4. Set common axis limits and colormap
        axis([floor(min(x)), ceil(max(x)), 0, 6]); 
        caxis(caxis_val);
        h_ax.Colormap = cmap_val;
        
        set(gca,'Fontsize',font_size_small,'Fontweight','b'); 
        title({plot_title}, 'fontweight','b','fontsize',font_size_small);  
        ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size_small); 
        
        % 5. X-axis formatting: Plot labels on ALL panels for guaranteed width
        set(gca, 'XTick', xData);
        set(gca,'XMinorTick','on');
        xAx = get(gca,'XAxis');
        xAx.MinorTickValues=xData_m;
        datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
        
        % Only show X-label on the bottom panel
        if i ~= num_panels 
            xlabel('');
            %set(gca, 'XTickLabel', []); % <-- Original line (keep commented out as before)
        else
            % Add X-label to the bottom panel
            xlabel('Time (UTC)','fontweight','b','fontsize',font_size_small); 
        end
        
        set(gca,'TickDir','out');
        set(gca,'TickLength',[0.005; 0.0025]);
    end
    figure_list{end+1} = figure_idx;
    suffix_list{end+1} = 'T_WV_7Panel_Difference'; % NEW SUFFIX
    figure_idx=figure_idx+1; % Increment figure index
end

% ---------------------------------------------
% --- FIGURE 3: 4-PANEL CORE PARAMETERS (Counts, ABC, AH, T) ---
% ---------------------------------------------
if flag.plot_aerosol_counts % Using this flag to control the new plot
    
    % Define custom plot size for 4 vertical panels
    plot_size_4panel = [100 100 1200 650]; 
    hf = figure(figure_idx); % Figure 3 (or next available index)
    set(hf, 'Position', plot_size_4panel, 'renderer', 'zbuffer');
    
    num_panels = 4;
    
    % --- ADD DEFINITION FOR GEOMETRY VARIABLES (Use 0.03/0.10 for consistency) ---
    Vertical_gap = 0.06;
    Total_Vertical_Margin = 0.10;
    
    % Define the height and gap of each panel for the 4-panel configuration
    total_gap = Vertical_gap * (num_panels - 1);
    total_height = 1 - Total_Vertical_Margin - total_gap;
    panel_height = total_height / num_panels;

    core_plots = {
        comb_Counts, 'Attenuated Backscatter, 828 nm (Log Scale)', caxis_Counts, 'jet', 1;       % Data, Title, Clim, CMap, LogScale
        comb_ABC, 'Aerosol Backscatter Coefficient, PTV (m^{-1} sr^{-1})', caxis_ABC, 'viridis', 1;
        comb_AH_PTV, 'Absolute Humidity PTV (g m^{-3})', caxis_AH, cmap_AH, 0;
        comb_T, 'Temperature PTV (K)', caxis_T_abs, plasma_cmap, 0;
    };


%    sgtitle({[node, ': Core Lidar Parameters']}, 'fontweight','b','fontsize',font_size);
    
    for i = 1:num_panels
        % Calculate normalized bottom position: moves from top (high index i) to bottom (low index i)
        bottom_pos = Total_Vertical_Margin/2 + (num_panels - i) * (panel_height + Vertical_gap);

        % 1. Create subplot using explicit position
        h_ax = subplot('Position', [SUBPLOT_LEFT, bottom_pos, SUBPLOT_WIDTH, panel_height]);
        
        data_plot = real(core_plots{i, 1});
        plot_title = core_plots{i, 2};
        caxis_val = core_plots{i, 3};
        cmap_val = core_plots{i, 4};
        log_scale = core_plots{i, 5};
        
        % 2. Pcolor plot
        h = pcolor(h_ax, x, y, data_plot); 
        set(h, 'EdgeColor', 'none'); 
        axis xy; 
        
        % 3. Apply Log Scale if required (Counts and ABC)
        if log_scale
            set(gca, 'Colorscale', 'log');
        else
            set(gca, 'Colorscale', 'linear');
        end
        
        % 4. Add colorbar
        h_cb = colorbar('peer', h_ax, 'EastOutside');
        
        % Set dynamic colorbar label
        if contains(plot_title, 'Counts')
            ylabel(h_cb, 'Counts', 'fontweight', 'b', 'fontsize', font_size_small - 2); 
        elseif contains(plot_title, 'Aerosol')
            ylabel(h_cb, 'ABC (m^{-1} sr^{-1})', 'fontweight', 'b', 'fontsize', font_size_small - 2); 
        elseif contains(plot_title, 'Humidity')
            ylabel(h_cb, 'Absolute Humidity (g m^{-3})', 'fontweight', 'b', 'fontsize', font_size_small - 2); 
        elseif contains(plot_title, 'Temperature')
            ylabel(h_cb, 'Temperature (K)', 'fontweight', 'b', 'fontsize', font_size_small - 2); 
        end

        % 5. Set common axis limits and colormap
        axis([floor(min(x)), ceil(max(x)), 0, 6]);
        caxis(caxis_val);
        
        % *** FIX: Use the FUNCTIONAL FORM for colormap assignment to handle 'jet' and 'viridis' strings ***
        colormap(h_ax, cmap_val); 
        % h_ax.Colormap = cmap_val; % <-- Original failing line
        % ---------------------------------------------------------------------------------------------------
        
        % Apply SMALLER FONT SIZE
        set(gca,'Fontsize',font_size_small,'Fontweight','b'); 
        
        % 6. Titles and Labels
        title({plot_title}, 'fontweight','b','fontsize',font_size_small);  
        ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size_small); 
        
        % 7. X-axis formatting
        set(gca, 'XTick', xData);
        set(gca,'XMinorTick','on');
        xAx = get(gca,'XAxis');
        xAx.MinorTickValues=xData_m;
        datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
        
        % Only show X-label on the bottom panel
        if i ~= num_panels 
            xlabel('');
            %set(gca, 'XTickLabel', []); % Explicitly hide tick labels
        else
            % Add X-label to the bottom panel
            xlabel('Time (UTC)','fontweight','b','fontsize',font_size_small); 
        end
        
        set(gca,'TickDir','out');
        set(gca,'TickLength',[0.005; 0.0025]);
    end
    figure_list{end+1} = figure_idx;
    suffix_list{end+1} = 'Core_4Panel_Parameters';
    figure_idx=figure_idx+1; % Increment figure index
end



% ---------------------------------------------------
% --- SINGLE-PANEL PLOTS (CONDITIONAL) ---
% ---------------------------------------------------

% Single Temperature and WV panels (REDUNDANT with multi-plots, only run if requested)
if flag.plot_single_temp_panels || flag.plot_single_wv_panels
    warning('Redundant single plots enabled. Use flag.plot_single_*_panels=0 for efficiency.');
end

% Single Temperature Panels 
if flag.plot_single_temp_panels
    % T_model
    create_pcolor_plot(x, y, comb_T_model, 'Temp, ERA5 Model (K)', caxis_T_abs, plasma_cmap, node, plot_size_wide, font_size, xData, xData_m, figure_idx, 0); 
    figure_list{end+1} = figure_idx; suffix_list{end+1} = 'T_Model_multi'; figure_idx=figure_idx+1;
    
    % T_PTV
    create_pcolor_plot(x, y, comb_T, 'Temp, PTV (K)', caxis_T_abs, plasma_cmap, node, plot_size_wide, font_size, xData, xData_m, figure_idx, 0); 
    figure_list{end+1} = figure_idx; suffix_list{end+1} = 'T_PTV_multi'; figure_idx=figure_idx+1;
    
    if flag.read_temp_std
        % T_Std
        create_pcolor_plot(x, y, comb_T_Std, 'Temp, Standard (K)', caxis_T_abs, plasma_cmap, node, plot_size_wide, font_size, xData, xData_m, figure_idx, 0); 
        figure_list{end+1} = figure_idx; suffix_list{end+1} = 'T_Standard_multi'; figure_idx=figure_idx+1;

        % T_Std Difference
        create_pcolor_plot(x, y, comb_T_Std-comb_T_model, 'Temperature, Standard-ERA5 (K)', caxis_diff_T, smooth_gray_cmap_T, node, plot_size_wide, font_size, xData, xData_m, figure_idx, 0);
        figure_list{end+1} = figure_idx; suffix_list{end+1} = 'T_Diff_Standard_Model'; figure_idx=figure_idx+1;
    end
    
    % T_PTV Difference
    create_pcolor_plot(x, y, comb_T-comb_T_model, 'Temperature, PTV-ERA5 (K)', caxis_diff_T, smooth_gray_cmap_T, node, plot_size_wide, font_size, xData, xData_m, figure_idx, 0);
    figure_list{end+1} = figure_idx; suffix_list{end+1} = 'T_Diff_PTV_Model'; figure_idx=figure_idx+1;
end

% Single Water Vapor Panels 
if flag.plot_single_wv_panels
    % Absolute AH plots
    create_pcolor_plot(x, y, comb_AH, 'Absolute Humidity Standard (g m^{-3})', caxis_AH, cmap_AH, node, plot_size_wide, font_size, xData, xData_m, figure_idx, 0); 
    figure_list{end+1} = figure_idx; suffix_list{end+1} = 'WV_Standard_multi'; figure_idx=figure_idx+1;

    create_pcolor_plot(x, y, comb_AH_PTV, 'Absolute Humidity PTV (g m^{-3})', caxis_AH, cmap_AH, node, plot_size_wide, font_size, xData, xData_m, figure_idx, 0); 
    figure_list{end+1} = figure_idx; suffix_list{end+1} = 'WV_PTV_multi'; figure_idx=figure_idx+1;

    if flag.read_wv_multi
        create_pcolor_plot(x, y, comb_AH_MultiPulse, 'Absolute Humidity MultiPulse (g m^{-3})', caxis_AH, cmap_AH, node, plot_size_wide, font_size, xData, xData_m, figure_idx, 0); 
        figure_list{end+1} = figure_idx; suffix_list{end+1} = 'WV_MultiPulse_multi'; figure_idx=figure_idx+1;
    end
    
    create_pcolor_plot(x, y, comb_AH_model, 'Absolute Humidity ERA5 Model (g m^{-3})', caxis_AH, cmap_AH, node, plot_size_wide, font_size, xData, xData_m, figure_idx, 0); 
    figure_list{end+1} = figure_idx; suffix_list{end+1} = 'WV_Model_multi'; figure_idx=figure_idx+1;

    % Difference AH plots
    create_pcolor_plot(x, y, comb_AH-comb_AH_model, 'Absolute Humidity Standard-ERA5 (g m^{-3})', caxis_diff_AH, smooth_gray_cmap_AH, node, plot_size_wide, font_size, xData, xData_m, figure_idx, 0); 
    figure_list{end+1} = figure_idx; suffix_list{end+1} = 'WV_Diff_Standard_Model'; figure_idx=figure_idx+1;

    create_pcolor_plot(x, y, comb_AH_PTV-comb_AH_model, 'Absolute Humidity PTV-ERA5 (g m^{-3})', caxis_diff_AH, smooth_gray_cmap_AH, node, plot_size_wide, font_size, xData, xData_m, figure_idx, 0); 
    figure_list{end+1} = figure_idx; suffix_list{end+1} = 'WV_Diff_PTV_Model'; figure_idx=figure_idx+1;

    if flag.read_wv_multi
        create_pcolor_plot(x, y, comb_AH_MultiPulse-comb_AH_model, 'Absolute Humidity, MultiPulse-ERA5 (g m^{-3})', caxis_diff_AH, smooth_gray_cmap_AH, node, plot_size_wide, font_size, xData, xData_m, figure_idx, 0); 
        figure_list{end+1} = figure_idx; suffix_list{end+1} = 'WV_Diff_MultiPulse_Model'; figure_idx=figure_idx+1;
    end
end


% AEROSOL/COUNTS (Conditional)
if flag.plot_aerosol_counts
    % Figure X
    create_pcolor_plot(x, y, comb_ABC, 'Aerosol Backscatter Coefficient, PTV (m^{-1} sr^{-1})', caxis_ABC, 'viridis', node, plot_size_wide, font_size, xData, xData_m, figure_idx, 1);
    figure_list{end+1} = figure_idx; suffix_list{end+1} = 'Back_Coeff_PTV_multi'; figure_idx=figure_idx+1;
    
    % Figure X+1
    create_pcolor_plot(x, y, comb_Counts, 'Attenuated Backscatter, 828 nm', caxis_Counts, 'jet', node, plot_size_wide, font_size, xData, xData_m, figure_idx, 1); 
    figure_list{end+1} = figure_idx; suffix_list{end+1} = '828_Counts_multi'; figure_idx=figure_idx+1;
end

% UNCERTAINTY (Conditional)
if flag.plot_uncertainty && flag.read_uncertainty
    % Figure X
    create_pcolor_plot(x, y, comb_AH_var, 'Absolute Humidity Standard Uncertainty (g m^{-3})', caxis_uncertainty, smooth_gray_cmap_AH, node, plot_size_wide, font_size, xData, xData_m, figure_idx, 0); 
    figure_list{end+1} = figure_idx; suffix_list{end+1} = 'WV_Standard_Uncertainty'; figure_idx=figure_idx+1;
    
    % Figure X+1
    create_pcolor_plot(x, y, comb_AH_PTV_var, 'Absolute Humidity PTV Uncertainty (g m^{-3})', caxis_uncertainty, smooth_gray_cmap_AH, node, plot_size_wide, font_size, xData, xData_m, figure_idx, 0); 
    figure_list{end+1} = figure_idx; suffix_list{end+1} = 'WV_PTV_Uncertainty'; figure_idx=figure_idx+1;
    
    if flag.read_wv_multi
        % Figure X+2
        create_pcolor_plot(x, y, comb_AH_MultiPulse_var, 'Absolute Humidity MultiPulse Uncertainty (g m^{-3})', caxis_uncertainty, smooth_gray_cmap_AH, node, plot_size_wide, font_size, xData, xData_m, figure_idx, 0); 
        figure_list{end+1} = figure_idx; suffix_list{end+1} = 'WV_MultiPulse_Uncertainty'; figure_idx=figure_idx+1;
    end
    
    if flag.read_temp_std
        % Figure X+3
        create_pcolor_plot(x, y, comb_T_Std_var, 'Temperature Standard Uncertainty (K)', caxis_uncertainty, smooth_gray_cmap_T, node, plot_size_wide, font_size, xData, xData_m, figure_idx, 0); 
        figure_list{end+1} = figure_idx; suffix_list{end+1} = 'T_Standard_Uncertainty'; figure_idx=figure_idx+1; 
    end
elseif flag.plot_uncertainty && ~flag.read_uncertainty
    warning('Uncertainty plots requested, but flag.read_uncertainty is 0. Skipping uncertainty plots.');
end


% ---------------------------------------------
% --- AH VS. ERA5 SCATTER PLOT (CONDITIONAL) ---
% ---------------------------------------------
if flag.plot_ah_scatter
    % Figure X
    create_ah_scatter_plot(comb_AH_model, comb_AH, comb_AH_PTV, comb_AH_MultiPulse, y, ...
        figure_idx, node, plot_size_wide, font_size, flag.read_wv_multi);
    
    figure_list{end+1} = figure_idx; 
    suffix_list{end+1} = 'AH_vs_ERA5_Scatter'; 
    figure_idx=figure_idx+1; % Increment figure index
end

% ---------------------------------------------
% --- TEMP VS. ERA5 SCATTER PLOT (CONDITIONAL) ---
% ---------------------------------------------
if flag.plot_temp_scatter
    % Figure X+1
    create_temp_scatter_plot(comb_T_model, comb_T, comb_T_Std, ...
        figure_idx, node, plot_size_wide, font_size, flag.read_temp_std);
    
    figure_list{end+1} = figure_idx; 
    suffix_list{end+1} = 'T_vs_ERA5_Scatter'; 
    figure_idx=figure_idx+1; % Increment figure index
end

% ---------------------------------------------
% --- PRODUCT VS. MODEL DENSITY MAPS (CONDITIONAL) ---
% ---------------------------------------------
% AH Density Maps
if flag.plot_ah_density_map
    % AH PTV vs Model
    create_product_vs_model_density(comb_AH_PTV_full, comb_AH_model_full, ...
        'AH PTV vs ERA5 Model Density Map', 'ERA5 Model AH (g m^{-3})', 'Lidar PTV AH (g m^{-3})', ...
        AH_density_limits, DENSITY_CLIM, figure_idx, node, plot_size_wide, font_size, cmap_AH); % 
    figure_list{end+1} = figure_idx; 
    suffix_list{end+1} = 'AH_PTV_vs_ERA5_Density'; 
    figure_idx=figure_idx+1;

    % AH Standard vs Model
    create_product_vs_model_density(comb_AH_full, comb_AH_model_full, ...
        'AH Standard vs ERA5 Model Density Map', 'ERA5 Model AH (g m^{-3})', 'Lidar Standard AH (g m^{-3})', ...
        AH_density_limits, DENSITY_CLIM, figure_idx, node, plot_size_wide, font_size, cmap_AH); %
    figure_list{end+1} = figure_idx; 
    suffix_list{end+1} = 'AH_Std_vs_ERA5_Density'; 
    figure_idx=figure_idx+1;
    
    if flag.read_wv_multi
        % AH MultiPulse vs Model
        create_product_vs_model_density(comb_AH_MultiPulse_full, comb_AH_model_full, ...
            'AH MultiPulse vs ERA5 Model Density Map', 'ERA5 Model AH (g m^{-3})', 'Lidar MultiPulse AH (g m^{-3})', ...
            AH_density_limits, DENSITY_CLIM, figure_idx, node, plot_size_wide, font_size, cmap_AH); 
        figure_list{end+1} = figure_idx; 
        suffix_list{end+1} = 'AH_MP_vs_ERA5_Density'; 
        figure_idx=figure_idx+1;
    end
end

% Temperature Density Maps
if flag.plot_t_density_map
    % T PTV vs Model
    create_product_vs_model_density(comb_T_full, comb_T_model_full, ...
        'T PTV vs ERA5 Model Density Map', 'ERA5 Model T (K)', 'Lidar PTV T (K)', ...
        T_density_limits, DENSITY_CLIM, figure_idx, node, plot_size_wide, font_size, cmap_density); % <-- USES cmap_density (Magma)
    figure_list{end+1} = figure_idx; 
    suffix_list{end+1} = 'T_PTV_vs_ERA5_Density'; 
    figure_idx=figure_idx+1;
    
    if flag.read_temp_std
        % T Standard vs Model
        create_product_vs_model_density(comb_T_Std_full, comb_T_model_full, ...
            'T Standard vs ERA5 Model Density Map', 'ERA5 Model T (K)', 'Lidar Standard T (K)', ...
            T_density_limits, DENSITY_CLIM, figure_idx, node, plot_size_wide, font_size, cmap_density); % <-- USES cmap_density (Magma)
        figure_list{end+1} = figure_idx; 
        suffix_list{end+1} = 'T_Std_vs_ERA5_Density'; 
        figure_idx=figure_idx+1;
    end
end




% HISTOGRAMS and 2D DENSITY MAPS (Conditional)
if flag.plot_histograms_2d || flag.plot_1d_histograms
    
    % 1D HISTOGRAMS
    if flag.plot_1d_histograms
        % Figure X
        % Note: Histograms use FULL resolution data for better statistical accuracy
        create_1d_histogram(comb_T_full - comb_T_model_full, 'T PTV - ERA5 Model Difference', 'Temperature Difference (K)', T_hist_limits, figure_idx, node, font_size); 
        figure_list{end+1} = figure_idx; suffix_list{end+1} = 'T_Diff_Histogram'; figure_idx=figure_idx+1;
        
        % Figure X+1
        create_1d_histogram(comb_AH_full - comb_AH_model_full, 'AH Standard - ERA5 Model Difference', 'Absolute Humidity Difference (g m\textsuperscript{-3})', AH_hist_limits, figure_idx, node, font_size); 
        figure_list{end+1} = figure_idx; suffix_list{end+1} = 'AH_Diff_Histogram'; figure_idx=figure_idx+1;
    end

    % 2D DENSITY MAPS
    if flag.plot_histograms_2d
        y_range = y;
        
        % Define data based on decimation flag
        if flag.decimate_hist_data
            % Use decimated data for plotting (faster)
            T_PTV_diff_hist = comb_T - comb_T_model;
            T_Std_diff_hist = comb_T_Std - comb_T_model;
            AH_Std_diff_hist = comb_AH - comb_AH_model;
            AH_PTV_diff_hist = comb_AH_PTV - comb_AH_model;
            AH_MultiPulse_diff_hist = comb_AH_MultiPulse - comb_AH_model;
        else
            % Use full resolution data for high accuracy (slower)
            T_PTV_diff_hist = comb_T_full - comb_T_model_full;
            T_Std_diff_hist = comb_T_Std_full - comb_T_model_full;
            AH_Std_diff_hist = comb_AH_full - comb_AH_model_full;
            AH_PTV_diff_hist = comb_AH_PTV_full - comb_AH_model_full;
            AH_MultiPulse_diff_hist = comb_AH_MultiPulse_full - comb_AH_model_full;
        end

        
        % Figure X+2
        create_2d_density_map(T_PTV_diff_hist, y_range, ...
            'T PTV - ERA5 Difference vs. Range (Bin Width: 0.2 K)', 'Temperature Difference (K)', ...
            T_diff_bin_limits, [0 6], T_diff_bin_width, DENSITY_CLIM, figure_idx, node, plot_size_wide, font_size, cmap_density); 
        figure_list{end+1} = figure_idx; suffix_list{end+1} = 'T_Diff_PTV_Histogram_Range'; figure_idx=figure_idx+1;
        
        if flag.read_temp_std
            % Figure X+3
            create_2d_density_map(T_Std_diff_hist, y_range, ... 
                'T Standard - ERA5 Difference vs. Range (Bin Width: 0.2 K)', 'Temperature Difference (K)', ...
                T_diff_bin_limits, [0 6], T_diff_bin_width, DENSITY_CLIM, figure_idx, node, plot_size_wide, font_size, cmap_density); 
            figure_list{end+1} = figure_idx; suffix_list{end+1} = 'T_Diff_Standard_Histogram_Range'; figure_idx=figure_idx+1;
        end

        % Figure X+4
        create_2d_density_map(AH_Std_diff_hist, y_range, ...
            'AH Standard - ERA5 Difference vs. Range (Bin Width: 0.2 g/m\textsuperscript{3})', 'Absolute Humidity Difference (g/m\textsuperscript{3})', ...
            AH_diff_bin_limits, [0 6], AH_diff_bin_width, DENSITY_CLIM, figure_idx, node, plot_size_wide, font_size, cmap_AH); 
        figure_list{end+1} = figure_idx; suffix_list{end+1} = 'AH_Diff_Standard_Histogram_Range'; figure_idx=figure_idx+1;
        
        % Figure X+5
        create_2d_density_map(AH_PTV_diff_hist, y_range, ...
            'AH PTV - ERA5 Difference vs. Range (Bin Width: 0.2 g/m\textsuperscript{3})', 'Absolute Humidity Difference (g/m\textsuperscript{3})', ...
            AH_diff_bin_limits, [0 6], AH_diff_bin_width, DENSITY_CLIM, figure_idx, node, plot_size_wide, font_size, cmap_AH); 
        figure_list{end+1} = figure_idx; suffix_list{end+1} = 'AH_PTV_Diff_Histogram_Range'; figure_idx=figure_idx+1;

        if flag.read_wv_multi
            % Figure X+6
            create_2d_density_map(AH_MultiPulse_diff_hist, y_range, ...
                'AH MultiPulse - ERA5 Difference vs. Range (Bin Width: 0.2 g/m\textsuperscript{3})', 'Absolute Humidity Difference (g/m\textsuperscript{3})', ...
                AH_diff_bin_limits, [0 6], AH_diff_bin_width, DENSITY_CLIM, figure_idx, node, plot_size_wide, font_size, cmap_AH); 
            figure_list{end+1} = figure_idx; suffix_list{end+1} = 'AH_MultiPulse_Diff_Histogram_Range'; figure_idx=figure_idx+1;
        end
    end
end


% --- PLOT SAVING LOOP (Iterating over created figures only) ---
cd(plot_dir);
plot_size = [100 100 1920 225];
plot_size_square = [100 100 750 750]; 
plot_size_4panel = [100 100 1900 650]; 
plot_size_5panel = [100 100 1900 800]; 
plot_size_7panel = [100 100 1900 890]; 
num_plots = length(figure_list);

for k_idx = 1:num_plots
    k = figure_list{k_idx}; % Get the actual figure number
    suffix = suffix_list{k_idx}; % Get the corresponding suffix
    
    FigH = figure(k);
    drawnow;
    FigH.Units = 'pixels';
    
    % Determine size based on the specific figure created
    if strcmp(suffix, 'T_WV_7Panel_Absolute')
        FigH.Position = plot_size_7panel;
    elseif strcmp(suffix, 'T_WV_7Panel_Difference')
        FigH.Position = plot_size_7panel;
    elseif strcmp(suffix, 'Core_4Panel_Parameters')
        FigH.Position = plot_size_4panel;
%    elseif contains(suffix, 'Histogram_Range') || contains(suffix, 'Histogram') || strcmp(suffix, 'AH_vs_ERA5_Scatter') || strcmp(suffix, 'T_vs_ERA5_Scatter')l% 
    elseif contains(suffix, 'Histogram') || contains(suffix, 'Density') || contains(suffix, 'Scatter')
       FigH.Position = plot_size_square;
    else
        % Default wide size for single pcolor plots
        FigH.Position = plot_size; 
    end
    
    name = char(strcat(node, "_", date_str, '_', suffix));
    exportgraphics(FigH, [name, '.png'], 'Resolution', 150);
    
    % --- CRITICAL EFFICIENCY IMPROVEMENT ---
   % close(FigH); % Close the figure immediately to free up memory
end

disp(['Plotting and saving complete. Time elapsed: ', num2str(toc(tic_plot_save)), ' seconds.']);
disp(['Total script runtime: ', num2str(toc(tic_total)), ' seconds.']);
disp('--- Script End ---');


% ylim([0 3])
% colormap(viridis)
%clim([0 20])



% ------------------------------------------------------------------
% --- LOCAL FUNCTIONS ---
% ------------------------------------------------------------------

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


function cmap = create_diverging_cmap(N, T_max, T_gray)
% CREATE_DIVERGING_CMAP Generates a custom diverging colormap (Blue-Gray-Red)
% with a fixed gray center that corresponds to a user-defined data range.
%
%   N:       Total number of colors in the map (e.g., 256).
%   T_max:   The maximum absolute data value (e.g., 5 for +/- 5K).
%   T_gray:  The absolute data value defining the extent of the gray center 
%            (e.g., 2 for +/- 2K).

% --- Define Key Colors ---
blue = [0 0 1];
red = [1 0 0];
light_gray = [0.9 0.9 0.9];

% --- Determine Indices for the Gray Band ---
% Calculate the ratio of the gray band range to the total range
ratio_gray = T_gray / T_max; 
if ratio_gray >= 1
    error('T_gray must be less than T_max.');
end

% N_gray is the total number of indices that will be dedicated to light_gray.
N_gray = round(N * ratio_gray);
if mod(N_gray, 2) ~= 0
    N_gray = N_gray + 1; % Ensure N_gray is an even number
end

N_half = N / 2;
N_gray_half = N_gray / 2;

% N_transition is the number of indices that will handle the smooth fade
N_transition = N_half - N_gray_half;

% --- Build the First Half (Blue -> Transition -> Light Gray) ---
% Part A: Blue to Gray Transition
R1_trans = linspace(blue(1), light_gray(1), N_transition)';
G1_trans = linspace(blue(2), light_gray(2), N_transition)';
B1_trans = linspace(blue(3), light_gray(3), N_transition)';
cmap_half1_trans = [R1_trans, G1_trans, B1_trans];

% Part B: Solid Gray Center (Lower Half)
cmap_half1_gray = repmat(light_gray, N_gray_half, 1);

% Concatenate the first half: Transition + Solid Gray
cmap_half1 = [cmap_half1_trans; cmap_half1_gray];

% --- Build the Second Half (Light Gray -> Transition -> Red) ---
% Part A: Solid Gray Center (Upper Half)
cmap_half2_gray = repmat(light_gray, N_gray_half, 1);

% Part B: Gray to Red Transition
R2_trans = linspace(light_gray(1), red(1), N_transition)';
G2_trans = linspace(light_gray(2), red(2), N_transition)';
B2_trans = linspace(light_gray(3), red(3), N_transition)';
cmap_half2_trans = [R2_trans, G2_trans, B2_trans];

% Concatenate the second half: Solid Gray + Transition
cmap_half2 = [cmap_half2_gray; cmap_half2_trans];

% --- Final Colormap ---
cmap = [cmap_half1; cmap_half2];

end
% ==========================================================