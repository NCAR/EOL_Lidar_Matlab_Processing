function [hf] = create_7panel_difference_plot(flag, node, comb_AH, comb_AH_PTV, comb_AH_MultiPulse, comb_AH_model, ...
    comb_T, comb_T_Std, comb_T_model, x, y, avail, caxis_AH, cmap_AH, caxis_diff_AH, smooth_gray_cmap_AH, ...
    caxis_T_abs, caxis_diff_T, smooth_gray_cmap_T, plasma_cmap, MAX_ALT_KM, font_size_small, plot_size_7panel, xData, xData_m)
% CREATE_7PANEL_DIFFERENCE_PLOT generates Figure 2, the time-height plots for 
% ERA5 Model and the differences between Lidar products and the ERA5 model.

    hf = figure(); 
    set(hf, 'Position', plot_size_7panel, 'renderer', 'zbuffer');
    
    % FIX: Normalized Plotting Coordinates for Subplot Alignment
    SUBPLOT_LEFT = 0.1;
    SUBPLOT_WIDTH = 0.78; 

    % Initial list of all 7 plots (Model and Differences, Top to Bottom)
    all_diff_plots = {
        % ABSOLUTE HUMIDITY
        comb_AH_model, 'Absolute Humidity ERA5 Model (g m^{-3})', caxis_AH, cmap_AH, 'Model', false; % Model Plot
        comb_AH - comb_AH_model, 'AH Standard - Model Difference (g m^{-3})', caxis_diff_AH, smooth_gray_cmap_AH, 'AH_Std', true;
        comb_AH_PTV - comb_AH_model, 'AH PTV - Model Difference (g m^{-3})', caxis_diff_AH, smooth_gray_cmap_AH, 'AH_PTV', true;
        comb_AH_MultiPulse - comb_AH_model, 'AH MultiPulse - Model Difference (g m^{-3})', caxis_diff_AH, smooth_gray_cmap_AH, 'AH_MultiPulse', true;
        % TEMPERATURE
        comb_T_model, 'Temperature ERA5 Model (K)', caxis_T_abs, plasma_cmap, 'Model', false; % Model Plot
        comb_T_Std - comb_T_model, 'T Standard - Model Difference (K)', caxis_diff_T, smooth_gray_cmap_T, 'T_Std', true;
        comb_T - comb_T_model, 'T PTV - Model Difference (K)', caxis_diff_T, smooth_gray_cmap_T, 'T_PTV', true;
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
        product_tag = all_diff_plots{i, 5};
        is_diff_plot = all_diff_plots{i, 6}; % New flag to control colorbar label
        
        % 2. Pcolor plot
        h = pcolor(h_ax, x, y, data_plot); 
        set(h, 'EdgeColor', 'none'); 
        axis xy; 
        
        % 3. Add colorbar
        h_cb = colorbar('peer', h_ax, 'EastOutside');
        
        % Set dynamic colorbar label
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
        else
            % Add X-label to the bottom panel
            xlabel('Time (UTC)','fontweight','b','fontsize',font_size_small); 
        end
        
        % 6. Add Availability Text Annotation
        avail_val = NaN; % Default to NaN
        if strcmp(product_tag, 'AH_Std')
            avail_val = avail.AH_Std;
        elseif strcmp(product_tag, 'AH_PTV')
            avail_val = avail.AH_PTV;
        elseif strcmp(product_tag, 'AH_MultiPulse')
            avail_val = avail.AH_MultiPulse;
        elseif strcmp(product_tag, 'T_Std')
            avail_val = avail.T_Std;
        elseif strcmp(product_tag, 'T_PTV')
            avail_val = avail.T_PTV;
        end
        
        if isfinite(avail_val)
            % Format the text string
            avail_str = {['<', num2str(MAX_ALT_KM), ' km Avail: '], [num2str(avail_val, '%.1f'), '%']};
            % Use the 'text' function in normalized axes units (relative to the subplot)
            text(0.99, 0.9, avail_str, 'Units', 'normalized', ...
                 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
                 'FontSize', font_size_small - 2, 'FontWeight', 'bold', 'Color', 'k');
        end
        
        set(gca,'TickDir','out');
        set(gca,'TickLength',[0.005; 0.0025]);
    end
end