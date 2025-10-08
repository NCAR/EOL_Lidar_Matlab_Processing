clear all; close all

skip = 1;
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

cd(data_dir); % Move to data directory to select files

% --- 2. FILE SELECTION ---
[Pythonfilename, Pythondir] = uigetfile('*.*','Select the sonde file', 'MultiSelect', 'on');
num_files = size(Pythonfilename, 2);

if num_files == 0
    error('No files selected. Script terminated.');
end

% --- 3. VARIABLE DEFINITION ---
% Using a cell array for clearer variable mapping (indices must match original)
variables = cell(1, 41);
variables([1:12, 16:18, 26:28, 30:32, 40:41]) = {
    'time', 'range', 'Temperature_PTV', 'Temperature_PTV_mask', 'Temperature_PTV_uncertainty', ...
    'Absolute_Humidity_Standard', 'Absolute_Humidity_Standard_mask', 'Absolute_Humidity_Standard_uncertainty', ...
    'Aerosol_Backscatter_Coefficient_PTV', 'Aerosol_Backscatter_Coefficient_PTV_mask', ...
    'Aerosol_Backscatter_Coefficient_PTV_uncertainty', 'Backscatter_Photon_Counts_828', ...
    'Absolute_Humidity_PTV', 'Absolute_Humidity_PTV_mask', 'Absolute_Humidity_PTV_uncertainty', ...
    'Absolute_Humidity_MultiPulse', 'Absolute_Humidity_MultiPulse_mask', 'Absolute_Humidity_MultiPulse_uncertainty', ...
    'Surface_Temperature', 'Surface_Pressure', 'Surface_Absolute_Humidity', ...
    'Absolute_Humidity_ERA5', 'Temperature_ERA5'
};

% --- 4. DATA READING AND MASKING LOOP ---
for jj = 1:num_files
    filename = Pythonfilename{jj};
    date_str = filename(end-15:end-8);
    n = datenum(date_str, 'yyyymmdd');
    
    ncid = netcdf.open(filename, 'NC_NOWRITE');
    ncdisp(filename) % use this to display all variables
    

    % Read data
    time{jj} = ncread(filename,variables{1}); 
    alt{jj} = ncread(filename,variables{2});
        
    % Get dimensions for fallback NaNs
    num_ranges_file = size(alt{jj}, 1);
    num_timesteps_file = size(time{jj}, 1);
 
    T{jj}  = ncread(filename,variables{3});  
    T_mask{jj} = ncread(filename,variables{4}); 
    T_var{jj} = ncread(filename,variables{5}); 
    T_model{jj}  = ncread(filename,variables{41}); 
    T_surf{jj} =  ncread(filename,variables{30});
    P_surf{jj} =  ncread(filename,variables{31});
    AH_surf{jj} =  ncread(filename,variables{32});
    AH{jj}  = ncread(filename,variables{6});  
    AH_mask{jj} = ncread(filename,variables{7}); 
    AH_var{jj} = ncread(filename,variables{8}); 
    AH_PTV{jj}  = ncread(filename,variables{16});  
    AH_PTV_mask{jj} = ncread(filename,variables{17}); 
    AH_PTV_var{jj} = ncread(filename,variables{18});
    try
      AH_MultiPulse{jj}  = ncread(filename,variables{26});  
      AH_MultiPulse_mask{jj} = ncread(filename,variables{27}); 
      AH_MultiPulse_var{jj} = ncread(filename,variables{28});
    catch
       warning(['MultiPulse data missing in file: ', filename]);
       % If variables are missing, fill array with NaNs of the correct size
       AH_MultiPulse{jj} = nan(num_ranges_file, num_timesteps_file);
       AH_MultiPulse_mask{jj} = nan(num_ranges_file, num_timesteps_file);
       AH_MultiPulse_var{jj} = nan(num_ranges_file, num_timesteps_file);
    end
    AH_model{jj}  = ncread(filename,variables{40});  
    ABC{jj}  = ncread(filename,variables{9});   
    ABC_mask{jj} = ncread(filename,variables{10}); 
    ABC_var{jj} = ncread(filename,variables{11});
    Counts{jj} = ncread(filename,variables{12});   


    % Masking (Vectorized for clarity)
    T{jj}(T_mask{jj} == 1) = nan;
    T_var{jj}(T_mask{jj} == 1) = nan; 

    mask_ah = (AH_mask{jj} == 1) | (AH_var{jj} > 25);
    AH{jj}(mask_ah) = nan;
    AH_var{jj}(AH_mask{jj} == 1) = nan; % Only mask uncertainty based on mask

    mask_ah1 = (AH_PTV_mask{jj} == 1) | (AH_PTV_var{jj} > 5);
    AH_PTV{jj}(mask_ah1) = nan;
    AH_PTV_var{jj}(AH_PTV_mask{jj} == 1) = nan; 

    mask_ah2 = (AH_MultiPulse_mask{jj} == 1) | (AH_MultiPulse_var{jj} > 5);
    AH_MultiPulse{jj}(mask_ah2) = nan;
    AH_MultiPulse_var{jj}(AH_MultiPulse_mask{jj} == 1) = nan; 
    
    ABC{jj}(ABC_mask{jj} == 1) = nan;

    netcdf.close(ncid); 
    
    % Convert from Unix time to date number
    duration{jj} = n + double(time{jj}/86400); % 86400 = 3600*24
end

% --- 5. PREALLOCATION AND COMBINATION (Your optimized section) ---
total_timesteps = sum(cellfun('size', duration, 1)); 
num_ranges = size(AH{1}, 1); 

% Preallocate 2D arrays (Range x Time)
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
comb_T_model = nan(num_ranges, total_timesteps);
comb_Counts = nan(num_ranges, total_timesteps);
% Preallocate 1D arrays (Time x 1)
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
    
    % Fill 2D data (Columns are time)
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
    comb_T_model(:, start_col:end_col) = T_model{jj};
    comb_Counts(:, start_col:end_col) = Counts{jj};
    
    start_col = end_col + 1;
end



% --- NaN GAP INSERTION FOR VISUAL BREAKS (Final Fix: Insert Two NaN Columns) ---

% Calculate the time difference (in days) between consecutive profiles
dt = diff(comb_duration);

% Gap threshold: If the time jump is greater than 2 days, assume it's a data gap
gap_threshold = 2.0; 

% Find indices immediately preceding a gap
gap_indices = find(dt > gap_threshold);

% Data arrays to modify (defined earlier)
data_arrays_2d = {'comb_AH', 'comb_AH_var', 'comb_AH_PTV', 'comb_AH_PTV_var', 'comb_AH_MultiPulse', 'comb_AH_MultiPulse_var', 'comb_AH_model', 'comb_ABC', 'comb_ABC_var', 'comb_T', 'comb_T_var', 'comb_T_model', 'comb_Counts'};
data_arrays_1d = {'comb_duration', 'comb_T_surf', 'comb_P_surf', 'comb_AH_surf'};

if ~isempty(gap_indices)
    disp(['Found ', num2str(length(gap_indices)), ' time gaps exceeding ', num2str(gap_threshold), ' days. Inserting 2 NaNs...']);
    
    % Loop through the detected gaps in reverse order to maintain correct indexing
    for i = length(gap_indices):-1:1
        idx = gap_indices(i); % Index *before* the gap starts
        
        % The new time points will be slightly offset from the data point *before* the gap, 
        % and slightly offset from the data point *after* the gap.
        
        gap_start_time = comb_duration(idx);
        gap_end_time = comb_duration(idx + 1);
        
        % Define two insertion time points to span the gap visually:
        % 1. Time slightly after the last valid point (to terminate smearing)
        nan_time_1 = gap_start_time + (gap_end_time - gap_start_time) * 0.001;
        % 2. Time slightly before the next valid point (to start the blank period)
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
        for arr = data_arrays_2d
            array_name = arr{1};
            current_array = eval(array_name);
            
            % Insert the two NaN columns at index idx + 1
            current_array = [current_array(:, 1:idx), nan_cols, current_array(:, idx+1:end)];
            eval([array_name ' = current_array;']);
        end
    end
else
    disp('No significant time gaps found.');
end
 




% --- 6. PLOT SETUP AND CALLS ---
scrsz = get(0,'ScreenSize');
date_str = datestr(comb_duration(1), 'yyyy-mmm-dd'); % Use first day of combined data
plot_size1 = [scrsz(4)/1.5 scrsz(4)/10 scrsz(3)/1.5 scrsz(4)/3];
font_size = 16;
x = comb_duration;
y = (alt{1}./1000); % Assumes 'alt' is consistent across files
xData =  linspace( fix(min(x)),  ceil(max(x)), round((ceil(max(x))-fix(min(x)))/skip)+1 );
xData_m =  linspace( fix(min(x)),  ceil(max(x)), round((ceil(max(x))-fix(min(x)))/skip*24/2)+1 );

% Define the smooth_gray_cmap (Differential plots)
N = 64;
N_half = N / 2;
blue = [0 0 1];
red = [1 0 0];
gray = [0.9 0.9 0.9];
R1 = linspace(blue(1), gray(1), N_half)'; G1 = linspace(blue(2), gray(2), N_half)'; B1 = linspace(blue(3), gray(3), N_half)';
cmap_half1 = [R1, G1, B1];
R2 = linspace(gray(1), red(1), N_half)'; G2 = linspace(gray(2), red(2), N_half)'; B2 = linspace(gray(3), red(3), N_half)';
cmap_half2 = [R2, G2, B2];
smooth_gray_cmap = [cmap_half1; cmap_half2];

% Define common parameters for plots
cmap_AH = CM_YlGnBu(64); % Assuming CM_YlGnBu is a custom function
caxis_AH = [0 25];
caxis_diff_AH = [-10 10];
caxis_diff_T = [-5 5];
caxis_uncertainty = [-5 5]; 
caxis_T_abs = [268 298];
caxis_ABC = [5e-9 1e-6];
caxis_Counts = [1e2 1e6];

% --- Plot Calls using the function (function 'create_pcolor_plot' not included here) ---
create_pcolor_plot(x, y, comb_T, 'Temp, PTV (K)', caxis_T_abs, 'plasma', node, plot_size1, font_size, xData, xData_m, 1, 0);
create_pcolor_plot(x, y, comb_ABC, 'Aerosol Backscatter Coefficient, PTV (m^{-1} sr^{-1})', caxis_ABC, 'viridis', node, plot_size1, font_size, xData, xData_m, 2, 1);
create_pcolor_plot(x, y, comb_Counts, 'Attenuated Backscatter, 828 nm', caxis_Counts, 'jet', node, plot_size1, font_size, xData, xData_m, 3, 1);
create_pcolor_plot(x, y, comb_AH, 'Absolute Humidity Standard (g m^{-3})', caxis_AH, cmap_AH, node, plot_size1, font_size, xData, xData_m, 4, 0);
create_pcolor_plot(x, y, comb_AH_PTV, 'Absolute Humidity PTV (g m^{-3})', caxis_AH, cmap_AH, node, plot_size1, font_size, xData, xData_m, 5, 0);
create_pcolor_plot(x, y, comb_AH_MultiPulse, 'Absolute Humidity MultiPulse (g m^{-3})', caxis_AH, cmap_AH, node, plot_size1, font_size, xData, xData_m, 6, 0);
create_pcolor_plot(x, y, comb_AH_model, 'Absolute Humidity ERA5 Model (g m^{-3})', caxis_AH, cmap_AH, node, plot_size1, font_size, xData, xData_m, 7, 0);
create_pcolor_plot(x, y, comb_AH-comb_AH_model, 'Absolute Humidity Standard-ERA5 (g m^{-3})', caxis_diff_AH, smooth_gray_cmap, node, plot_size1, font_size, xData, xData_m, 8, 0);
create_pcolor_plot(x, y, comb_AH_PTV-comb_AH_model, 'Absolute Humidity PTV-ERA5 (g m^{-3})', caxis_diff_AH, smooth_gray_cmap, node, plot_size1, font_size, xData, xData_m, 9, 0);
create_pcolor_plot(x, y, comb_AH_MultiPulse-comb_AH_model, 'Absolute Humidity, MultiPulse-ERA5 (g m^{-3})', caxis_diff_AH, smooth_gray_cmap, node, plot_size1, font_size, xData, xData_m, 10, 0);
create_pcolor_plot(x, y, comb_T-comb_T_model, 'Temperature, PTV-ERA5 (K)', caxis_diff_T, smooth_gray_cmap, node, plot_size1, font_size, xData, xData_m, 11, 0);
create_pcolor_plot(x, y, comb_AH_var, 'Absolute Humidity Standard Uncertainty (g m^{-3})', caxis_uncertainty, smooth_gray_cmap, node, plot_size1, font_size, xData, xData_m, 12, 0);
create_pcolor_plot(x, y, comb_AH_PTV_var, 'Absolute Humidity PTV Uncertainty (g m^{-3})', caxis_uncertainty, smooth_gray_cmap, node, plot_size1, font_size, xData, xData_m, 13, 0);
create_pcolor_plot(x, y, comb_AH_MultiPulse_var, 'Absolute Humidity MultiPulse Uncertainty (g m^{-3})', caxis_uncertainty, smooth_gray_cmap, node, plot_size1, font_size, xData, xData_m, 14, 0);

% Define parameters for density maps
AH_diff_bin_limits = [-10 10]; % Used for X-axis display and binning for AH
T_diff_bin_limits = [-10 10];
T_hist_limits = [-15 15];
AH_hist_limits = [-10 10];

T_diff_bin_width = 0.2;
AH_diff_bin_width = 0.2;

cmap_density = flipud(magma(256));
y_range = y;
plot_size_square = [100 100 750 750]; 

% --- 1D HISTOGRAMS (Figures 15, 17) ---
create_1d_histogram(comb_T - comb_T_model, 'T PTV - ERA5 Model Difference', 'Temperature Difference (K)', T_hist_limits, 15, node, font_size);
create_1d_histogram(comb_AH - comb_AH_model, 'AH Standard - ERA5 Model Difference', 'Absolute Humidity Difference (g m\textsuperscript{-3})', AH_hist_limits, 17, node, font_size);

% --- 2D DENSITY MAPS (Figures 16, 18, 19, 20) ---
create_2d_density_map(comb_T - comb_T_model, y_range, ...
    'T PTV - ERA5 Difference vs. Range (Bin Width: 0.2 K)', 'Temperature Difference (K)', ...
    T_diff_bin_limits, [0 6], T_diff_bin_width, 16, node, plot_size1, font_size, cmap_density);

create_2d_density_map(comb_AH - comb_AH_model, y_range, ...
    'AH Standard - ERA5 Difference vs. Range (Bin Width: 0.2 g/m\textsuperscript{3})', 'Absolute Humidity Difference (g/m\textsuperscript{3})', ...
    AH_diff_bin_limits, [0 6], AH_diff_bin_width, 18, node, plot_size1, font_size, cmap_density);

create_2d_density_map(comb_AH_PTV - comb_AH_model, y_range, ...
    'AH PTV - ERA5 Difference vs. Range (Bin Width: 0.2 g/m\textsuperscript{3})', 'Absolute Humidity Difference (g/m\textsuperscript{3})', ...
    AH_diff_bin_limits, [0 6], AH_diff_bin_width, 19, node, plot_size1, font_size, cmap_density);

create_2d_density_map(comb_AH_MultiPulse - comb_AH_model, y_range, ...
    'AH MultiPulse - ERA5 Difference vs. Range (Bin Width: 0.2 g/m\textsuperscript{3})', 'Absolute Humidity Difference (g/m\textsuperscript{3})', ...
    AH_diff_bin_limits, [0 6], AH_diff_bin_width, 20, node, plot_size1, font_size, cmap_density);


% --- PLOT SAVING LOOP (Major consolidation) ---
cd(plot_dir);
% plot_size = [100 100 1920*2 225*2];
plot_size = [100 100 1920 225];
plot_size_square = [100 100 750 750]; 
num_plots = 20;

filename_suffixes = {
    'T_PTV_multi', 'Back_Coeff_PTV_multi', '828_Counts_multi', ...
    'WV_Standard_multi', 'WV_PTV_multi', 'WV_MultiPulse_multi', ...
    'WV_Model_multi', 'WV_Diff_Standard_Model', 'WV_Diff_PTV_Model', ...
    'WV_Diff_MultiPulse_Model', 'T_Diff_PTV_Model', ...
    'WV_Standard_Uncertainty', 'WV_PTV_Uncertainty', 'WV_MultiPulse_Uncertainty',...
    'T_Diff_Histogram', 'T_Diff_Histogram_Range',  'AH_Diff_Histogram', 'AH_Diff_Histogram_Range', ...
    'AH_PTV_Diff_Histogram_Range', 'AH_MultiPulse_Diff_Histogram_Range'    
};



for k = 1:num_plots
    FigH = figure(k);
    drawnow;
    FigH.Units = 'pixels';
    if k >= 15
        FigH.Position = plot_size_square;
    else
        FigH.Position = plot_size;
    end
    name = char(strcat(node, "_", date_str, '_', filename_suffixes{k}));
    exportgraphics(FigH, [name, '.png'], 'Resolution', 150);
end