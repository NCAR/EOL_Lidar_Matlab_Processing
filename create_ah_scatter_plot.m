function create_ah_scatter_plot(ah_model, ah_standard, ah_ptv, ah_multipulse, y_range, fig_num, node, plot_size1, font_size, flag_wv_multi)
% Function to create a scatter plot comparing three AH products to the ERA5 model.
%
% Inputs:
%   ah_model:      ERA5 Absolute Humidity (X-axis data).
%   ah_standard:   Standard AH product (Y-axis data 1).
%   ah_ptv:        PTV AH product (Y-axis data 2).
%   ah_multipulse: MultiPulse AH product (Y-axis data 3, optional).
%   y_range:       Range/altitude data (used for altitude filtering, not plotting).
%   fig_num:       Figure number.
%   node:          String for the node name.
%   plot_size1:    Figure position array.
%   font_size:     Font size for text.
%   flag_wv_multi: Boolean (1 if MultiPulse data is included).

    hf = figure(fig_num);
    set(hf, 'Position', plot_size1, 'renderer', 'zbuffer');
    
    % --- 1. Prepare Data ---
    % Flatten 2D matrices into 1D vectors for scatter plotting
    ah_model_flat = ah_model(:);
    ah_standard_flat = ah_standard(:);
    ah_ptv_flat = ah_ptv(:);
    ah_multipulse_flat = ah_multipulse(:);
    
    % --- 2. Filter out NaNs for all main datasets ---
    % Create a master valid index that includes all three required AH sets
    valid_idx_std = ~isnan(ah_model_flat) & ~isnan(ah_standard_flat);
    valid_idx_ptv = ~isnan(ah_model_flat) & ~isnan(ah_ptv_flat);
    
    % Filter the data
    x_std = ah_model_flat(valid_idx_std);
    y_std = ah_standard_flat(valid_idx_std);

    x_ptv = ah_model_flat(valid_idx_ptv);
    y_ptv = ah_ptv_flat(valid_idx_ptv);
    
    % Prepare plot for multiple series on the same axes
    hold on;
    
    % --- 3. Plot Scatter Series ---
    % Plot 1: AH Standard
    s1 = scatter(x_std, y_std, 10, 'b', 'o', 'filled', 'DisplayName', 'AH Standard'); % Blue circles
    s1.MarkerFaceAlpha = 0.5; % Set transparency

    % Plot 2: AH PTV
    s2 = scatter(x_ptv, y_ptv, 10, 'r', 's', 'filled', 'DisplayName', 'AH PTV');      % Red squares
    s2.MarkerFaceAlpha = 0.5;

    % Plot 3: AH MultiPulse (Conditional)
    if flag_wv_multi
        valid_idx_mp = ~isnan(ah_model_flat) & ~isnan(ah_multipulse_flat);
        x_mp = ah_model_flat(valid_idx_mp);
        y_mp = ah_multipulse_flat(valid_idx_mp);
        s3 = scatter(x_mp, y_mp, 10, 'g', '^', 'filled', 'DisplayName', 'AH MultiPulse'); % Green triangles
        s3.MarkerFaceAlpha = 0.5;
    end
    
    % --- 4. Plot 1:1 Line (Reference) ---
    max_val = ceil(max([x_std; y_std; x_ptv; y_ptv])); 
    min_val = floor(min([x_std; y_std; x_ptv; y_ptv])); 
    plot([min_val max_val], [min_val max_val], 'k--', 'LineWidth', 2, 'DisplayName', '1:1 Line');

    % --- 5. Formatting ---
    title({sprintf('%s Absolute Humidity Products vs. ERA5 Model', node)}, 'fontweight','b','fontsize',font_size);  
    xlabel('ERA5 Model Absolute Humidity (g m^{-3})', 'fontweight','b','fontsize',font_size); 
    ylabel('Lidar Absolute Humidity (g m^{-3})', 'fontweight','b','fontsize',font_size); 
    
    axis equal; % Ensure axes have equal scaling
    xlim([0 25]); % Fixed limits for AH
    ylim([0 25]); % Fixed limits for AH
    
    legend('Location', 'SouthEast', 'FontSize', font_size-4);
    grid on;
    box on;
    set(gca,'Fontsize',font_size,'Fontweight','b');
    
    hold off;
end