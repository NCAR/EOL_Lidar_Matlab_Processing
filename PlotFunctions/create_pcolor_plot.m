function hf = create_pcolor_plot(x, y, data, plot_title, caxis_val, cmap, node, plot_size1, font_size, xData, xData_m, figure_num, log_scale)
% Function to create and format a pcolor plot with common settings.
%
% Inputs:
%   x, y:          X and Y axis data (time and altitude).
%   data:          The 2D data matrix to plot (Z).
%   plot_title:    String for the plot title (e.g., 'Temp, PTV (K)').
%   clim_val:      [min max] for caxis/clim.
%   cmap:          Colormap to use (e.g., 'plasma', viridis, or custom map).
%   node:          String for the node name (e.g., 'MPD04').
%   plot_size1:    Figure position array.
%   font_size:     Font size for text.
%   xData:         Main XTick values.
%   xData_m:       Minor XTick values.
%   figure_num:    Figure number.
%   log_scale:     Boolean (1 for log colorscale, 0 otherwise).

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