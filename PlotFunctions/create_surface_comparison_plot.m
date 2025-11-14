function create_surface_comparison_plot(x, AH_surf_full, AH_MP_dec, AH_model_dec, node, font_size, xData, xData_m, figure_idx, N_decimate)
% CREATE_SURFACE_COMPARISON_PLOT Generates a time-series plot comparing the AH surface 
% observation with the lowest valid bins of MultiPulse and ERA5 data, using a 
% manual median filter to smooth data anomalies.

    % --- 1. CONFIGURATION ---
    FILTER_ORDER = 20; 
    
    hf = figure(figure_idx);
    set(hf, 'Position', [100 100 1000 350], 'renderer', 'zbuffer'); % Custom wide plot size

    % --- 1D Surface Data Decimation (CRITICAL FIX) ---
    % Input AH_surf_full is the full-resolution data.
    % We need to decimate it to match the size of x (N_time_decimated).
    AH_surf_decimated = AH_surf_full(:); % Ensure column vector for decimation
    
    % Pad with NaNs if necessary (using logic derived from decimate_matrix)
    N_total = length(AH_surf_decimated);
    N_remainder = mod(N_total, N_decimate);
    if N_remainder ~= 0
        N_pad = N_decimate - N_remainder;
        padding = nan(N_pad, 1);
        AH_surf_decimated = [AH_surf_decimated; padding];
        N_total = length(AH_surf_decimated);
    end
    
    % Calculate the average over the N_decimate bins
    AH_surf_decimated = mean(reshape(AH_surf_decimated, N_decimate, N_total/N_decimate), 1);
    
    % Extract the already decimated lowest bin data
    AH_MP_lowest = AH_MP_dec(5, :); 
    AH_Model_lowest = AH_model_dec(1, :); 
    
    % Ensure data is a row vector for filtering consistency
    AH_surf_decimated = AH_surf_decimated(:)'; 
    AH_MP_lowest = AH_MP_lowest(:)'; 
    AH_Model_lowest = AH_Model_lowest(:)'; 
    
    % --- 2. APPLY MANUAL MEDIAN FILTER (Toolbox-independent) ---
    
    AH_surf_filt = manual_medfilt1(AH_surf_decimated, FILTER_ORDER); % Apply filter to decimated data
    AH_MP_lowest_filt = manual_medfilt1(AH_MP_lowest, FILTER_ORDER);
    AH_Model_lowest_filt = manual_medfilt1(AH_Model_lowest, FILTER_ORDER);


    % --- 3. PLOT THE FILTERED TIME SERIES ---
    
    color_MP = [0 0.6 0];  % Define the MultiPulse color for consistency with profile plots 

    plot(x, AH_surf_filt, 'k-', 'LineWidth', 2, 'DisplayName', 'Surface Station AH (Filtered)'); hold on;
    plot(x, AH_MP_lowest_filt, 'Color', color_MP, 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', 'MultiPulse Lowest Bin (Filtered)');
    plot(x, AH_Model_lowest_filt, 'r:', 'LineWidth', 1.5, 'DisplayName', 'ERA5 Model Lowest Bin (Filtered)');

    grid on; box on;

    title({[node, ': AH Lowest Bin vs. Surface Station Time Series (Median Filtered)'], ...
           ['Comparison Period: ', datestr(min(x), 'dd-mmm-yy'), ' to ', datestr(max(x), 'dd-mmm-yy')]}, ...
           'fontweight', 'b', 'fontsize', font_size);

    xlabel('Time (UTC)', 'fontweight', 'b', 'fontsize', font_size);
    ylabel('Absolute Humidity (g m^{-3})', 'fontweight', 'b', 'fontsize', font_size);

    % Apply X-axis formatting
    set(gca, 'XTick', xData);
    set(gca, 'XMinorTick', 'on');
    xAx = get(gca, 'XAxis');
    xAx.MinorTickValues = xData_m;
    datetick('x', 'dd-mmm-yy', 'keeplimits', 'keepticks');

    legend('show', 'Location', 'best', 'fontsize', font_size-4);
    set(gca, 'Fontsize', font_size-2, 'Fontweight', 'b', 'TickDir', 'out');

    % Set Y-axis limits based on AH max/min of the filtered data
    min_data = min([AH_surf_filt, AH_MP_lowest_filt, AH_Model_lowest_filt], [], 'all', 'omitnan');
    max_data = max([AH_surf_filt, AH_MP_lowest_filt, AH_Model_lowest_filt], [], 'all', 'omitnan');
    
    if isfinite(min_data) && isfinite(max_data)
        ylim([min_data - 1, max_data + 1]);
    end

    hold off;
end