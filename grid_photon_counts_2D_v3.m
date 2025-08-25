function [totalCounts_inGridInterval, griddedFluxRate_countsPerSec] = grid_photon_counts_2D_V3(dataArray, time_grid, fixedGridInterval_seconds)
% grid_photon_counts_2D_V3 Resamples raw photon data by calculating a weighted-average flux.
%
%   [totalCounts_inGridInterval, griddedFluxRate_countsPerSec] = grid_photon_counts_2D_V3(dataArray, time_grid, fixedGridInterval_seconds)
%
%   This function addresses the issue of varying original accumulation times by:
%   1. Calculating the flux rate for each raw measurement based on its specific duration.
%   2. Gridding these flux rates using a time-weighted average to a new time base.
%   3. Converting the gridded flux back to total counts for the new grid intervals.
%
%   Inputs:
%     dataArray: N x M array with time stamps in fractional hours (column 1)
%                and photon counts (columns 2:end).
%     time_grid: A vector of time points (in fractional hours) defining the
%                new fixed time interval grid.
%     fixedGridInterval_seconds: The fixed time interval in seconds for the
%                                'time_grid'.
%
%   Outputs:
%     totalCounts_inGridInterval: P x (M-1) array of total photon counts within
%                                 each fixed time interval.
%     griddedFluxRate_countsPerSec: P x (M-1) array of the average photon flux
%                                   rate in counts/second for each interval.

% --- 1. Pre-processing: Remove duplicates and calculate time intervals ---
[~, ia, ~] = unique(dataArray(:,1), 'rows');
uA = dataArray(ia, :);

time_hours = uA(:, 1);
photonCountsData = uA(:, 2:end);

% Calculate the duration (in seconds) for each measurement.
deltaTime_hours = diff(time_hours);
% Handle the first point by assuming its interval is the same as the second.
deltaTime_hours = [deltaTime_hours(1); deltaTime_hours];

deltaTime_seconds = deltaTime_hours * 3600;
% Ensure deltaTime_seconds does not have non-positive values.
min_valid_dt = min(deltaTime_seconds(deltaTime_seconds > 0));
if isempty(min_valid_dt) || min_valid_dt == 0
    min_valid_dt = 1e-9;
end
deltaTime_seconds(deltaTime_seconds <= 0) = min_valid_dt;

% --- 2. Calculate the raw flux rate for each measurement ---
% This normalizes each measurement to a "per second" rate.
fluxRate_countsPerSec = photonCountsData ./ deltaTime_seconds;

% --- 3. Grid the flux data using a time-weighted average ---
time_grid_edges = [time_grid; time_grid(end) + (time_grid(2) - time_grid(1))];
bin_indices = discretize(time_hours, time_grid_edges);

num_channels = size(photonCountsData, 2);
num_grid_points = length(time_grid);
griddedFluxRate_countsPerSec = zeros(num_grid_points, num_channels);

% Loop through each channel and calculate the weighted average flux rate.
for colIdx = 1:num_channels
    currentFlux = fluxRate_countsPerSec(:, colIdx);
    
    % Use 'accumarray' to sum the weighted flux and the total weights for each bin.
    % The third argument [num_grid_points, 1] ensures the output array has the
    % correct size, filling empty bins with the default value (0).
    sum_of_flux_x_weight = accumarray(bin_indices, currentFlux .* deltaTime_seconds, [num_grid_points, 1], @sum, 0);
    sum_of_weights = accumarray(bin_indices, deltaTime_seconds, [num_grid_points, 1], @sum, 0);
    
    % Avoid division by zero for empty bins.
    zero_weights = (sum_of_weights == 0);
    sum_of_weights(zero_weights) = 1; % Temporarily set to 1 to avoid NaN
    
    % Calculate the weighted average flux for the current channel.
    griddedFluxRate_countsPerSec(:, colIdx) = sum_of_flux_x_weight ./ sum_of_weights;
    
    % Set the result for empty bins back to zero.
    griddedFluxRate_countsPerSec(zero_weights, colIdx) = 0;
end

% --- 4. Convert the gridded flux back to total counts ---
% Total Counts = Average Flux Rate (counts/sec) * Time Interval (sec).
totalCounts_inGridInterval = griddedFluxRate_countsPerSec * fixedGridInterval_seconds;

end