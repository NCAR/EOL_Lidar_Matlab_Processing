function [med, sd] = calc_nan_stats(data)
% Calculates median and standard deviation while ignoring NaNs,
% acting as a fallback for nanmedian and nanstd.
    
    % Get the number of rows (range gates)
    N_rows = size(data, 1); 
    med = nan(N_rows, 1);
    sd = nan(N_rows, 1);

    for r = 1:N_rows
        % Extract the time series for the current range gate
        row_data = data(r, :); 
        
        % Remove NaNs
        valid_data = row_data(isfinite(row_data));
        
        if ~isempty(valid_data)
            med(r) = median(valid_data);
            sd(r) = std(valid_data);
        end
    end
end