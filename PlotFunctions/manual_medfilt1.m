function y_filtered = manual_medfilt1(y_input, order)
% MANUAL_MEDFILT1 Implements a toolbox-independent 1D median filter.
% Handles edges by truncating the window and omits NaN values.
    
    N = length(y_input);
    y_filtered = nan(1, N);
    
    % Ensure filter order is odd for a symmetrical window
    if mod(order, 2) == 0
        order = order - 1; 
    end
    
    half_window = floor(order / 2);
    
    for i = 1:N
        % Define the window boundaries, truncating at the edges
        start_idx = max(1, i - half_window);
        end_idx = min(N, i + half_window);
        
        % Extract data window
        window_data = y_input(start_idx:end_idx);
        
        % Calculate the median of the valid (non-NaN) data in the window
        valid_data = window_data(isfinite(window_data));
        
        if ~isempty(valid_data)
            y_filtered(i) = median(valid_data);
        end
    end
end