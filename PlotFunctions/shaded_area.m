function shaded_area(y_alt, median_data, std_data, color_val, alpha_val, display_name)
% SHADED_AREA Plots a shaded region representing the median +/- 3 * standard deviation.
    
    % Calculate upper and lower bounds using 3 * std_data
    upper_bound = median_data + 3 * std_data;
    lower_bound = median_data - 3 * std_data;
    
    % Create coordinates for the patch: [lower_rev, upper] vs [y, y_rev]
    x_patch = [lower_bound; flipud(upper_bound)];
    y_patch = [y_alt; flipud(y_alt)];
    
    % Remove NaNs and corresponding height values
    valid_indices = isfinite(x_patch) & isfinite(y_patch);
    x_patch_clean = x_patch(valid_indices);
    y_patch_clean = y_patch(valid_indices);

    % Plot the filled area
    h_fill = fill(x_patch_clean, y_patch_clean, color_val);
    set(h_fill, 'EdgeColor', 'none', 'FaceAlpha', alpha_val, 'DisplayName', display_name);
end