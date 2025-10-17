function create_temp_scatter_plot(t_model, t_ptv, t_standard, fig_num, node, plot_size1, font_size, flag_read_temp_std)
% Function to create a scatter plot comparing temperature products to the ERA5 model.
%
% Inputs:
%   t_model:           ERA5 Temperature (X-axis data).
%   t_ptv:             PTV Temperature product (Y-axis data 1).
%   t_standard:        Standard Temperature product (Y-axis data 2, optional).
%   fig_num:           Figure number.
%   node:              String for the node name.
%   plot_size1:        Figure position array.
%   font_size:         Font size for text.
%   flag_read_temp_std: Boolean (1 if Standard data is included).

    hf = figure(fig_num);
    set(hf, 'Position', plot_size1, 'renderer', 'zbuffer');
    
    % --- 1. Prepare Data ---
    % Flatten 2D matrices into 1D vectors for scatter plotting
    t_model_flat = t_model(:);
    t_ptv_flat = t_ptv(:);
    t_standard_flat = t_standard(:);
    
    % --- 2. Filter out NaNs for required datasets ---
    
    % PTV data (always plotted)
    valid_idx_ptv = ~isnan(t_model_flat) & ~isnan(t_ptv_flat);
    x_ptv = t_model_flat(valid_idx_ptv);
    y_ptv = t_ptv_flat(valid_idx_ptv);
    
    % Prepare plot for multiple series on the same axes
    hold on;
    
    % --- 3. Plot Scatter Series ---
    
    % Plot 1: Temperature PTV (Always plot first for visual layering)
    s1 = scatter(x_ptv, y_ptv, 10, 'r', 's', 'filled', 'DisplayName', 'T PTV'); % Red squares
    s1.MarkerFaceAlpha = 0.5; % Set transparency

    % Plot 2: Temperature Standard (Conditional)
    if flag_read_temp_std
        valid_idx_std = ~isnan(t_model_flat) & ~isnan(t_standard_flat);
        x_std = t_model_flat(valid_idx_std);
        y_std = t_standard_flat(valid_idx_std);
        
        s2 = scatter(x_std, y_std, 10, 'b', 'o', 'filled', 'DisplayName', 'T Standard'); % Blue circles
        s2.MarkerFaceAlpha = 0.5;
    end
    
    % --- 4. Plot 1:1 Line (Reference) ---
    % Calculate plot limits based on the data to ensure the 1:1 line covers the range
    max_val = ceil(max([t_model_flat(:); t_ptv_flat(:)]));
    min_val = floor(min([t_model_flat(:); t_ptv_flat(:)]));
    
    % Use fixed bounds for cleaner temperature comparison (260 K to 310 K is a good atmospheric range)
    plot_min = max(260, min_val);
    plot_max = min(310, max_val);
    
    plot([plot_min plot_max], [plot_min plot_max], 'k--', 'LineWidth', 2, 'DisplayName', '1:1 Line');

    % --- 5. Formatting ---
    title({sprintf('%s Temperature Products vs. ERA5 Model', node)}, 'fontweight','b','fontsize',font_size);  
    xlabel('ERA5 Model Temperature (K)', 'fontweight','b','fontsize',font_size); 
    ylabel('Lidar Temperature (K)', 'fontweight','b','fontsize',font_size); 
    
    axis equal; % Ensure axes have equal scaling
    xlim([plot_min plot_max]); 
    ylim([plot_min plot_max]); 
    
    legend('Location', 'SouthEast', 'FontSize', font_size-4);
    grid on;
    box on;
    set(gca,'Fontsize',font_size,'Fontweight','b');
    
    hold off;
end