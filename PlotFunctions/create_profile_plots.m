function create_profile_plots(full_time_axis, y_alt, AH, AH_PTV, AH_MultiPulse, AH_model, ...
                              T, T_Std, T_model, AH_surf, T_surf, ...
                              target_time, node, figure_idx)
% CREATE_PROFILE_PLOTS Generates line plots of vertical profiles at a specific time,
% showing the median and a 3-sigma shaded area over a 1-hour window.

    % --- 1. DEFINE TIME WINDOW ---
    % Window size is set back to 1 hour
    window_hr = 1;
    TIME_WINDOW_HALF = 0.5 / 24; % 30 minutes in datenum units
    
    start_time = target_time - TIME_WINDOW_HALF;
    end_time = target_time + TIME_WINDOW_HALF;

    % Find the closest single index for the warning message fallback
    [~, idx_closest_single] = min(abs(full_time_axis - target_time)); % Uses full_time_axis
    
    disp(' ');
    disp('--- PROFILE PLOT DIAGNOSTICS ---');
    fprintf('Profile Window Size: %d Hours\n', window_hr);
    fprintf('Target Time: %s\n', datestr(target_time, 'dd-mmm-yyyy HH:MM:SS'));
    fprintf('Search Start: %s (datenum: %.10f)\n', datestr(start_time, 'dd-mmm-yyyy HH:MM:SS'), start_time);
    fprintf('Search End:   %s (datenum: %.10f)\n', datestr(end_time, 'dd-mmm-yyyy HH:MM:SS'), end_time);
    
    % --- Robustness Adjustment: Use a small tolerance for floating point comparison ---
    TOLERANCE = 1e-10; % Equivalent to microseconds of time tolerance (very small)

    % Find indices for the window using tolerance
    time_indices = find(full_time_axis >= start_time - TOLERANCE & full_time_axis <= end_time + TOLERANCE); % Uses full_time_axis
    
    num_profiles_found = length(time_indices);
    fprintf('Profiles Found in Window: %d\n', num_profiles_found);
    
    % --- Determine which indices to use for statistics ---
    if num_profiles_found < 2 % Need at least two profiles to calculate variance (std)
        if num_profiles_found == 0
            % Fallback: Find the single closest profile index
            idx_to_use = idx_closest_single;
            time_indices = idx_to_use; % Use only the single closest point
            warning(['No data found in the %d-hour window around: ', datestr(target_time)], window_hr);
        else
            % One profile found, but need more for shaded area
            idx_to_use = time_indices(1);
            warning(['Only ONE profile found in the %d-hour window. Cannot calculate variability (shaded area will be zero).'], window_hr);
        end
        
        fprintf('Falling back to single closest profile at: %s\n', datestr(full_time_axis(idx_to_use), 'dd-mmm-yyyy HH:MM:SS'));
    else
        % Use the found indices for the window
        idx_to_use = time_indices;
    end
    
    % --- 2. EXTRACT WINDOW DATA & MODEL DATA ---
    % NOTE: This extracts only a few columns, preventing memory overload.
    AH_Std_window = AH(:, idx_to_use);
    AH_PTV_window = AH_PTV(:, idx_to_use);
    AH_MultiPulse_window = AH_MultiPulse(:, idx_to_use);
    T_PTV_window = T(:, idx_to_use);
    T_Std_window = T_Std(:, idx_to_use);
    
    % For model data, just take the first profile's time index if multiple were found
    if length(idx_to_use) > 1
         model_idx = idx_to_use(1);
    else
         model_idx = idx_to_use;
    end
    
    AH_Model_profile = AH_model(:, model_idx);
    T_Model_profile = T_model(:, model_idx);
    
    AH_Surf_window = AH_surf(idx_to_use);
    T_Surf_window = T_surf(idx_to_use);

    % --- 3. CALCULATE MEDIAN AND 1-SIGMA (STANDARD DEVIATION) ---
    [AH_Std_median, AH_Std_std] = calc_nan_stats(AH_Std_window);
    [AH_PTV_median, AH_PTV_std] = calc_nan_stats(AH_PTV_window);
    [AH_MultiPulse_median, AH_MultiPulse_std] = calc_nan_stats(AH_MultiPulse_window);
    
    [T_PTV_median, T_PTV_std] = calc_nan_stats(T_PTV_window);
    [T_Std_median, T_Std_std] = calc_nan_stats(T_Std_window);
    
    % Surface Data Median/Std for the time window
    AH_Surf_valid = AH_Surf_window(isfinite(AH_Surf_window));
    T_Surf_valid = T_Surf_window(isfinite(T_Surf_window));
    
    AH_Surf_median = median(AH_Surf_valid);
    AH_Surf_std = std(AH_Surf_valid);
    T_Surf_median = median(T_Surf_valid);
    T_Surf_std = std(T_Surf_valid);
    
    
    % Handle case where only one profile was used (std is 0 or NaN)
    if num_profiles_found < 2
        AH_Std_std(:) = 0; AH_PTV_std(:) = 0; AH_MultiPulse_std(:) = 0;
        T_PTV_std(:) = 0; T_Std_std(:) = 0;
        AH_Surf_std = 0; T_Surf_std = 0;
        
        time_str_range = ['Single Profile: ', datestr(full_time_axis(idx_closest_single), 'dd-mmm-yyyy HH:MM:SS UTC')];
    else
        time_str_range = [datestr(full_time_axis(idx_to_use(1)), 'dd-mmm-yyyy HH:MM') ' to ' datestr(full_time_axis(idx_to_use(end)), 'HH:MM UTC')];
    end
    disp('--- END DIAGNOSTICS ---');
    disp(' ');
    
    % --- 4. AH PROFILE PLOT (Median + Shaded Area) ---
    hf1 = figure(figure_idx);
    set(hf1, 'Position', [50 50 500 700]);
    hold on;
    
    % Define common colors and transparency
    color_Std = [0 0 1]; color_PTV = [1 0 0]; color_MP = [0 0.6 0];
    alpha_val = 0.15; % Transparency for the shaded area

    % --- AH Standard (Blue) ---
    shaded_area(y_alt, AH_Std_median, AH_Std_std, color_Std, alpha_val, 'AH Std $\pm 3\sigma$');
    plot(AH_Std_median, y_alt, '-', 'Color', color_Std, 'LineWidth', 2, 'DisplayName', 'AH Std Median');
    
    % --- AH PTV (Red) ---
    shaded_area(y_alt, AH_PTV_median, AH_PTV_std, color_PTV, alpha_val, 'AH PTV $\pm 3\sigma$');
    plot(AH_PTV_median, y_alt, '-', 'Color', color_PTV, 'LineWidth', 2, 'DisplayName', 'AH PTV Median');
    
    % --- AH MultiPulse (Green) ---
    shaded_area(y_alt, AH_MultiPulse_median, AH_MultiPulse_std, color_MP, alpha_val, 'AH MultiPulse $\pm 3\sigma$');
    plot(AH_MultiPulse_median, y_alt, '-', 'Color', color_MP, 'LineWidth', 2, 'DisplayName', 'AH MultiPulse Median');
    
    % --- AH ERA5 Model (Black Dashed) ---
    plot(AH_Model_profile, y_alt, 'k--', 'LineWidth', 2, 'DisplayName', 'AH ERA5 Model');

    % --- AH Surface Station (Magenta) ---
    % Plot median as the marker center
    plot(AH_Surf_median, 0, 'ms', 'MarkerSize', 10, 'MarkerFaceColor', 'm', 'DisplayName', 'AH Surface Median');
    % Plot 1-sigma as error bars (horizontal)
    errorbar(AH_Surf_median, 0, AH_Surf_std, 'horizontal', 'Color', 'm', 'LineWidth', 1.5, 'CapSize', 10, 'HandleVisibility', 'off');

    grid on; box on;
    title({[node, ' Absolute Humidity Profile (', num2str(window_hr), '-Hour Window)'], ['Window: ', time_str_range]}, 'fontweight', 'b');
    xlabel('Absolute Humidity (g m^{-3})', 'fontweight', 'b');
    ylabel('Height (km, AGL)', 'fontweight', 'b');
    legend('show', 'Location', 'northeast', 'interpreter', 'latex');
    ylim([0 ceil(max(y_alt))]);
    set(gca, 'Fontsize', 12, 'Fontweight', 'b');
    hold off;
    
    figure_idx = figure_idx + 1;
    
    % --- 5. TEMPERATURE PROFILE PLOT (Median + Shaded Area) ---
    hf2 = figure(figure_idx);
    set(hf2, 'Position', [550 50 500 700]);
    hold on;
    
    % --- T PTV (Red) ---
    shaded_area(y_alt, T_PTV_median, T_PTV_std, color_PTV, alpha_val, 'T PTV $\pm 3\sigma$');
    plot(T_PTV_median, y_alt, '-', 'Color', color_PTV, 'LineWidth', 2, 'DisplayName', 'T PTV Median');
    
    % --- T Standard (Blue) ---
    shaded_area(y_alt, T_Std_median, T_Std_std, color_Std, alpha_val, 'T Std $\pm 3\sigma$');
    plot(T_Std_median, y_alt, '-', 'Color', color_Std, 'LineWidth', 2, 'DisplayName', 'T Std Median');

    % --- T ERA5 Model (Black Dashed) ---
    plot(T_Model_profile, y_alt, 'k--', 'LineWidth', 2, 'DisplayName', 'T ERA5 Model');

    % --- T Surface Station (Magenta) ---
    % Plot median as the marker center
    plot(T_Surf_median, 0, 'ms', 'MarkerSize', 10, 'MarkerFaceColor', 'm', 'DisplayName', 'T Surface Median');
    % Plot 1-sigma as error bars (horizontal)
    errorbar(T_Surf_median, 0, T_Surf_std, 'horizontal', 'Color', 'm', 'LineWidth', 1.5, 'CapSize', 10, 'HandleVisibility', 'off');

    grid on; box on;
    title({[node, ' Temperature Profile (', num2str(window_hr), '-Hour Window)'], ['Window: ', time_str_range]}, 'fontweight', 'b');
    xlabel('Temperature (K)', 'fontweight', 'b');
    ylabel('Height (km, AGL)', 'fontweight', 'b');
    legend('show', 'Location', 'northeast', 'interpreter', 'latex');
    ylim([0 ceil(max(y_alt))]);
    set(gca, 'Fontsize', 12, 'Fontweight', 'b');
    hold off;
end