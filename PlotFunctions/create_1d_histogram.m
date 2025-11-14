function create_1d_histogram(data_vector, title_str, xlabel_str, bin_limits, fig_num, node, font_size)
    figure(fig_num);
    
    % Get screen size (Assuming scrsz is defined globally or passed)
    scrsz = get(0,'ScreenSize');
    set(gcf, 'Position', [scrsz(4)/1.5 scrsz(4)/10 scrsz(3)/3 scrsz(4)/3]);

    % Remove NaN values
    data_vector = data_vector(~isnan(data_vector));
    
    % Plot the histogram
    h_hist = histogram(data_vector, 'BinLimits', bin_limits, 'NumBins', 100, 'Normalization', 'probability');
    
    % --- Formatting ---
    title({[node, ' ', title_str]}, 'fontweight', 'b', 'fontsize', font_size);
    xlabel(xlabel_str, 'fontweight', 'b', 'fontsize', font_size);
    ylabel('Probability Density', 'fontweight', 'b', 'fontsize', font_size);
    grid on;
    box on;
    set(gca, 'Fontsize', font_size, 'Fontweight', 'b');
    
    % Calculate and display stats
    data_in_limits = data_vector(data_vector >= bin_limits(1) & data_vector <= bin_limits(2));
    bias = mean(data_in_limits);
    std_dev = std(data_in_limits);
    ylim_max = max(h_hist.Values); 
    
    unit = strsplit(xlabel_str, '(');
    unit = strtrim(unit{2}(1:end-1));
    text_str = sprintf('Bias: %.2f %s\nStd Dev: %.2f %s', bias, unit, std_dev, unit);
    text(max(xlim)*0.6, ylim_max*0.9, text_str, ...
         'Fontsize', font_size-2, 'FontWeight', 'b', 'BackgroundColor', 'w');
end
