function [hf] = create_4panel_core_plot(node, comb_Counts, comb_ABC, comb_AH_PTV, comb_T, x, y, avail, ...
    caxis_Counts, caxis_ABC, caxis_AH, cmap_AH, caxis_T_abs, plasma_cmap, ...
    MAX_ALT_KM, font_size_small, plot_size_4panel, xData, xData_m)
% CREATE_4PANEL_CORE_PLOT generates Figure 3, the time-height plots for 
% the core lidar and derived parameters (Counts, ABC, AH_PTV, T_PTV).

    hf = figure(); 
    set(hf, 'Position', plot_size_4panel, 'renderer', 'zbuffer');
    
    % FIX: Normalized Plotting Coordinates for Subplot Alignment
    SUBPLOT_LEFT = 0.1;
    SUBPLOT_WIDTH = 0.78; 
    
    num_panels = 4;
    
    % --- DEFINE GEOMETRY VARIABLES ---
    Vertical_gap = 0.06;
    Total_Vertical_Margin = 0.10;
    
    % Recalculate panel geometry
    total_gap = Vertical_gap * (num_panels - 1);
    total_height = 1 - Total_Vertical_Margin - total_gap;
    panel_height = total_height / num_panels;

    core_plots = {
        comb_Counts, 'Attenuated Backscatter, 828 nm (Log Scale)', caxis_Counts, 'jet', 1, 'Counts';
        comb_ABC, 'Aerosol Backscatter Coefficient, PTV (m^{-1} sr^{-1})', caxis_ABC, 'viridis', 1, 'ABC';
        comb_AH_PTV, 'Absolute Humidity PTV (g m^{-3})', caxis_AH, cmap_AH, 0, 'AH_PTV';
        comb_T, 'Temperature PTV (K)', caxis_T_abs, plasma_cmap, 0, 'T_PTV';
    };
    
    for i = 1:num_panels
        % Calculate normalized bottom position: moves from top to bottom
        bottom_pos = Total_Vertical_Margin/2 + (num_panels - i) * (panel_height + Vertical_gap);

        % 1. Create subplot using explicit position
        h_ax = subplot('Position', [SUBPLOT_LEFT, bottom_pos, SUBPLOT_WIDTH, panel_height]);
        
        data_plot = real(core_plots{i, 1});
        plot_title = core_plots{i, 2};
        caxis_val = core_plots{i, 3};
        cmap_val = core_plots{i, 4};
        log_scale = core_plots{i, 5};
        product_tag = core_plots{i, 6};
        
        % 2. Pcolor plot
        h = pcolor(h_ax, x, y, data_plot); 
        set(h, 'EdgeColor', 'none'); 
        axis xy; 
        
        % 3. Apply Log Scale if required
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
        
        colormap(h_ax, cmap_val); 
        
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
        else
            % Add X-label to the bottom panel
            xlabel('Time (UTC)','fontweight','b','fontsize',font_size_small); 
        end
        
        % 8. Add Availability Text Annotation (Only for AH and T)
        avail_val = NaN; % Default to NaN
        if strcmp(product_tag, 'AH_PTV')
            avail_val = avail.AH_PTV;
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