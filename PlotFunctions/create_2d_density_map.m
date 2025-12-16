function create_2d_density_map(diff_data, y_range, title_str, xlabel_str, x_bin_limits, y_bin_limits, x_bin_width, caxis_val, fig_num, node, plot_size1, font_size, cmap)
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

    set(gca, 'Colorscale', 'log'); % Apply log scale to the color axis
    clim(caxis_val);               % Apply the fixed color limits

    % --- Formatting ---
   % title({[node, ' ', title_str]}, 'fontweight', 'b', 'fontsize', font_size, 'Interpreter', 'tex'); 
   % xlabel(xlabel_str, 'fontweight', 'b', 'fontsize', font_size, 'Interpreter', 'latex');       
    title({[node, ' ', title_str]}, 'fontweight', 'b', 'fontsize', font_size); 
    xlabel(xlabel_str, 'fontweight', 'b', 'fontsize', font_size);    
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