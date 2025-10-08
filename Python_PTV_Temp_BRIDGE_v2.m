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
    
    % ncid is only needed for netcdf.close, ncread handles file opening/closing internally
    % unless you are reading a large number of variables from a single file, 
    % in which case keeping it open can be slightly faster. We'll stick to 
    % your method to ensure compatibility.
    ncid = netcdf.open(filename, 'NC_NOWRITE');
    
    % Read data
    time{jj} = ncread(filename,variables{1}); 
    alt{jj} = ncread(filename,variables{2});
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

    AH_MultiPulse{jj}  = ncread(filename,variables{26});  
    AH_MultiPulse_mask{jj} = ncread(filename,variables{27}); 
    AH_MultiPulse_var{jj} = ncread(filename,variables{28});

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


% --- 15. TEMPERATURE DIFFERENCE HISTOGRAM ---
% 1. Calculate the point-by-point difference (T_observed - T_model)
diff_T = comb_T - comb_T_model;
% 2. Flatten the 2D difference matrix into a single vector
diff_vector = diff_T(:);
% 3. Remove NaN values
diff_vector = diff_vector(~isnan(diff_vector));
% Define the bin limits for the histogram 
bin_limits = [-15 15]; 
figure(15);
set(gcf, 'Position', [scrsz(4)/1.5 scrsz(4)/10 scrsz(3)/3 scrsz(4)/3]);
% Plot the histogram using 'BinLimits' to restrict the data considered.
% 'NumBins' is optional; 100 is a good default count for continuous data.
h_hist = histogram(diff_vector, 'BinLimits', bin_limits, 'NumBins', 100, 'Normalization', 'probability');
% --- Histogram Formatting ---
title({[node, ' T PTV - ERA5 Model Difference']}, ...
       'fontweight', 'b', 'fontsize', font_size);
xlabel('Temperature Difference (K)', 'fontweight', 'b', 'fontsize', font_size);
ylabel('Probability Density', 'fontweight', 'b', 'fontsize', font_size);
grid on;
box on;
set(gca, 'Fontsize', font_size, 'Fontweight', 'b');

% Calculate and display mean and standard deviation (only of the data *inside* the limits)
% NOTE: The mean/std here calculates over ALL data points, including those outside the limits.
% To calculate the stats only for the plotted data:
data_in_limits = diff_vector(diff_vector >= bin_limits(1) & diff_vector <= bin_limits(2));

mean_diff = mean(data_in_limits);
std_diff = std(data_in_limits);

% Find the maximum height of the plotted bars
ylim_max = max(h_hist.Values); 

text_str = sprintf('Bias: %.2f K\nStd Dev: %.2f K', mean_diff, std_diff);
text(max(xlim)*0.6, ylim_max*0.9, text_str, ...
     'Fontsize', font_size-2, 'FontWeight', 'b', 'BackgroundColor', 'w');



% --- 16. T DIFFERENCE VS. RANGE DENSITY MAP (Revised for 0.2 K Bin Width) ---

% 1. Calculate the difference (T_observed - T_model)
diff_T = comb_T - comb_T_model;

% 2. Define the range/altitude vector (Y)
y_range = y; 

% 3. Create the 2D Histogram (Density Map)
figure(16);
set(gcf, 'Position', plot_size1);

% Reshape the data for histogram2 
num_timesteps = size(comb_T, 2);
alt_matrix = repmat(y_range, 1, num_timesteps);

% Flatten the altitude and difference matrices for histogram2
alt_vector = alt_matrix(:);
diff_vector = diff_T(:);

% Filter out NaN values 
valid_idx = ~isnan(diff_vector) & ~isnan(alt_vector);
alt_valid = alt_vector(valid_idx);
diff_valid = diff_vector(valid_idx);

% Define bin parameters
T_diff_bin_width = 0.2;  % Your requested horizontal width (0.2 K)
Alt_bin_width = 0.12;   % Keeps roughly 50 vertical bins over 6 km (6/50=0.12)

% Create the 2D histogram
% We use 'BinWidth' to ensure the T difference resolution is exactly 0.2 K.
h_hist2 = histogram2(diff_valid, alt_valid, ...
    'BinWidth', [T_diff_bin_width Alt_bin_width], ... % Enforce 0.2 K width
    'Normalization', 'probability', ...  
    'DisplayStyle', 'tile');            

% --- Histogram2 Formatting ---
title({[node, ' T PTV - ERA5 Difference vs. Range (Bin Width: 0.2 K)']}, ...
       'fontweight', 'b', 'fontsize', font_size);
xlabel('Temperature Difference (K)', 'fontweight', 'b', 'fontsize', font_size); 
ylabel('Height (km, AGL)', 'fontweight', 'b', 'fontsize', font_size); 
grid on;
set(gca, 'Layer', 'top'); 

% Use BinLimits to match the displayed range with the data range
h_hist2.XBinLimits = [-10 10]; 
h_hist2.YBinLimits = [0 6]; 

% Use the figure axis limits to confirm the view
xlim([-10 10]); 
ylim([0 6]); 

% Use a sequential colormap to show density
colormap(flipud(magma(256))); 
colorbar('EastOutside'); 

set(gca, 'Fontsize', font_size, 'Fontweight', 'b');
box on;


% --- 17 WV DIFFERENCE HISTOGRAM ---
% 1. Calculate the point-by-point difference 
diff_AH = comb_AH - comb_AH_model;
% 2. Flatten the 2D difference matrix into a single vector
diff_vector = diff_AH(:);
% 3. Remove NaN values
diff_vector = diff_vector(~isnan(diff_vector));
% Define the bin limits for the histogram 
bin_limits = [-10 10]; 
figure(17);
set(gcf, 'Position', [scrsz(4)/1.5 scrsz(4)/10 scrsz(3)/3 scrsz(4)/3]);
% Plot the histogram using 'BinLimits' to restrict the data considered.
% 'NumBins' is optional; 100 is a good default count for continuous data.
h_hist = histogram(diff_vector, 'BinLimits', bin_limits, 'NumBins', 100, 'Normalization', 'probability');
% --- Histogram Formatting ---
title({[node, ' WV PTV - ERA5 Model Difference']}, ...
       'fontweight', 'b', 'fontsize', font_size);
xlabel('WV Difference (g m)', 'fontweight', 'b', 'fontsize', font_size);
ylabel('Probability Density', 'fontweight', 'b', 'fontsize', font_size);
grid on;
box on;
set(gca, 'Fontsize', font_size, 'Fontweight', 'b');
% Calculate and display mean and standard deviation (only of the data *inside* the limits)
% NOTE: The mean/std here calculates over ALL data points, including those outside the limits.
% To calculate the stats only for the plotted data:
data_in_limits = diff_vector(diff_vector >= bin_limits(1) & diff_vector <= bin_limits(2));
mean_diff = mean(data_in_limits);
std_diff = std(data_in_limits);
% Find the maximum height of the plotted bars
ylim_max = max(h_hist.Values); 
text_str = sprintf('Bias: %.2f g/m3\nStd Dev: %.2f g/m3', mean_diff, std_diff);
text(max(xlim)*0.6, ylim_max*0.9, text_str, ...
     'Fontsize', font_size-2, 'FontWeight', 'b', 'BackgroundColor', 'w');



% --- 18. WV DIFFERENCE VS. RANGE DENSITY MAP  
% 1. Calculate the difference 
diff_AH = comb_AH - comb_AH_model;
AH_diff_bin_limits = [-20 20]; 
% 2. Define the range/altitude vector (Y)
y_range = y; 
% 3. Create the 2D Histogram (Density Map)
figure(18);
set(gcf, 'Position', plot_size1);
% Reshape the data for histogram2 
num_timesteps = size(comb_AH, 2);
alt_matrix = repmat(y_range, 1, num_timesteps);
% Flatten the altitude and difference matrices for histogram2
alt_vector = alt_matrix(:);
diff_vector = diff_AH(:);
% Filter out NaN values 
valid_idx = ~isnan(diff_vector) & ~isnan(alt_vector);
alt_valid = alt_vector(valid_idx);
diff_valid = diff_vector(valid_idx);
% Define bin parameters
AH_diff_bin_width = 0.2;  % horizontal width (0.2 g/m^3)
Alt_bin_width = 0.12;   % Keeps roughly 50 vertical bins over 6 km (6/50=0.12)
% Create the 2D histogram
% We use 'BinWidth' to ensure the WV difference resolution is exactly 0.2 g/m^3
h_hist2 = histogram2(diff_valid, alt_valid, ...
    'BinWidth', [AH_diff_bin_width Alt_bin_width], ... % Enforce 0.2 g/m^3 width
    'XBinLimits', AH_diff_bin_limits, ... % <--- FIX 1: RESTRICTS X-DATA FOR BINNING
    'Normalization', 'probability', ...  
    'DisplayStyle', 'tile');            
% --- Histogram2 Formatting ---
title({[node, ' AH Standard - ERA5 Difference vs. Range (Bin Width: 0.2 g/m^3)']}, ...
       'fontweight', 'b', 'fontsize', font_size);
xlabel('Absolute Humidity Difference (g/m^3)', 'fontweight', 'b', 'fontsize', font_size); 
ylabel('Height (km, AGL)', 'fontweight', 'b', 'fontsize', font_size); 
grid on;
set(gca, 'Layer', 'top'); 
% Use BinLimits to match the displayed range with the data range
h_hist2.XBinLimits = [-10 10]; 
h_hist2.YBinLimits = [0 6]; 
% Use the figure axis limits to confirm the view
xlim([-10 10]); 
ylim([0 6]); 
% Use a sequential colormap to show density
colormap(flipud(magma(256))); 
colorbar('EastOutside'); 
set(gca, 'Fontsize', font_size, 'Fontweight', 'b');
box on;

% --- 19. WV DIFFERENCE VS. RANGE DENSITY MAP  
% 1. Calculate the difference 
diff_AH_PTV = comb_AH_PTV - comb_AH_model;
% 2. Define the range/altitude vector (Y)
y_range = y; 
% 3. Create the 2D Histogram (Density Map
figure(19);
set(gcf, 'Position', plot_size1);
% Reshape the data for histogram2 
num_timesteps = size(comb_AH_PTV, 2);
alt_matrix = repmat(y_range, 1, num_timesteps);
% Flatten the altitude and difference matrices for histogram2
alt_vector = alt_matrix(:);
diff_vector = diff_AH_PTV(:);
% Filter out NaN values 
valid_idx = ~isnan(diff_vector) & ~isnan(alt_vector);
alt_valid = alt_vector(valid_idx);
diff_valid = diff_vector(valid_idx);
% Define bin parameters
AH_diff_bin_width = 0.2;  % horizontal width (0.2 g/m^3)
Alt_bin_width = 0.12;   % Keeps roughly 50 vertical bins over 6 km (6/50=0.12)
% Create the 2D histogram
% We use 'BinWidth' to ensure the WV difference resolution is exactly 0.2 g/m^3
h_hist2 = histogram2(diff_valid, alt_valid, ...
    'BinWidth', [AH_diff_bin_width Alt_bin_width], ... % Enforce 0.2 g/m^3 width
    'XBinLimits', AH_diff_bin_limits, ... % <--- FIX 1: RESTRICTS X-DATA FOR BINNING
    'Normalization', 'probability', ...  
    'DisplayStyle', 'tile');            
% --- Histogram2 Formatting ---
title({[node, ' AH PTV - ERA5 Difference vs. Range (Bin Width: 0.2 g/m^3)']}, ...
       'fontweight', 'b', 'fontsize', font_size);
xlabel('Absolute Humidity Difference (g/m^3)', 'fontweight', 'b', 'fontsize', font_size); 
ylabel('Height (km, AGL)', 'fontweight', 'b', 'fontsize', font_size); 
grid on;
set(gca, 'Layer', 'top'); 
% Use BinLimits to match the displayed range with the data range
h_hist2.XBinLimits = [-10 10]; 
h_hist2.YBinLimits = [0 6]; 
% Use the figure axis limits to confirm the view
xlim([-10 10]); 
ylim([0 6]); 
% Use a sequential colormap to show density
colormap(flipud(magma(256))); 
colorbar('EastOutside'); 
set(gca, 'Fontsize', font_size, 'Fontweight', 'b');
box on;

% --- 20. WV DIFFERENCE VS. RANGE DENSITY MAP  
% 1. Calculate the difference 
diff_AH_MultiPulse = comb_AH_MultiPulse - comb_AH_model;
% 2. Define the range/altitude vector (Y)
y_range = y; 
% 3. Create the 2D Histogram (Density Map)
figure(20);
set(gcf, 'Position', plot_size1);
% Reshape the data for histogram2 
num_timesteps = size(comb_AH_MultiPulse, 2);
alt_matrix = repmat(y_range, 1, num_timesteps);
% Flatten the altitude and difference matrices for histogram2
alt_vector = alt_matrix(:);
diff_vector = diff_AH_MultiPulse(:);
% Filter out NaN values 
valid_idx = ~isnan(diff_vector) & ~isnan(alt_vector);
alt_valid = alt_vector(valid_idx);
diff_valid = diff_vector(valid_idx);
% Define bin parameters
AH_diff_bin_width = 0.2;  % horizontal width (0.2 g/m^3)
Alt_bin_width = 0.12;   % Keeps roughly 50 vertical bins over 6 km (6/50=0.12)
% Create the 2D histogram
% We use 'BinWidth' to ensure the WV difference resolution is exactly 0.2 g/m^3
h_hist2 = histogram2(diff_valid, alt_valid, ...
    'BinWidth', [AH_diff_bin_width Alt_bin_width], ... % Enforce 0.2 g/m^3 width
    'XBinLimits', AH_diff_bin_limits, ... % <--- FIX 1: RESTRICTS X-DATA FOR BINNING
    'Normalization', 'probability', ...  
    'DisplayStyle', 'tile');            
% --- Histogram2 Formatting ---
title({[node, ' AH MultiPulse - ERA5 Difference vs. Range (Bin Width: 0.2 g/m^3)']}, ...
       'fontweight', 'b', 'fontsize', font_size);
xlabel('Absolute Humidity Difference (g/m^3)', 'fontweight', 'b', 'fontsize', font_size); 
ylabel('Height (km, AGL)', 'fontweight', 'b', 'fontsize', font_size); 
grid on;
set(gca, 'Layer', 'top'); 
% Use BinLimits to match the displayed range with the data range
h_hist2.XBinLimits = [-10 10]; 
h_hist2.YBinLimits = [0 6]; 
% Use the figure axis limits to confirm the view
xlim([-10 10]); 
ylim([0 6]); 
% Use a sequential colormap to show density
colormap(flipud(magma(256))); 
colorbar('EastOutside'); 
set(gca, 'Fontsize', font_size, 'Fontweight', 'b');
box on;

% --- 7. PLOT SAVING LOOP (Major consolidation) ---
cd(plot_dir);
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