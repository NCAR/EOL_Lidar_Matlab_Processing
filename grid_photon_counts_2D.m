function [totalCounts_inGridInterval, griddedFluxRate_countsPerSec] = grid_photon_counts_2D(dataArray, time_grid, fixedGridInterval_seconds)
% CONVERTPHOTONDATA Processes raw photon count data, calculates flux, grids it, and converts back to total counts.
%
%   [totalCounts_inGridInterval, griddedFluxRate_countsPerSec] = CONVERTPHOTONDATA(dataArray, time_grid, fixedGridInterval_seconds)
%
%   Inputs:
%     dataArray: A N x M array where N is the number of data points and M is
%                the number of channels + 1. The first column (dataArray(:,1))
%                must be time stamps in fractional hours. The remaining columns
%                (dataArray(:,2:end)) are photon counts.
%     time_grid: A vector of time points (in fractional hours) defining the
%                new fixed time interval grid to which the flux data will be
%                interpolated.
%     fixedGridInterval_seconds: The fixed time interval in seconds corresponding
%                                to the 'time_grid' (e.g., 2 seconds), used for
%                                converting the gridded flux back to total counts.
%
%   Outputs:
%     totalCounts_inGridInterval: An P x (M-1) array (where P is length of time_grid)
%                                 containing the total photon counts within the
%                                 specified fixedGridInterval_seconds for each channel,
%                                 after gridding.
%     griddedFluxRate_countsPerSec: An P x (M-1) array (where P is length of time_grid)
%                                   containing the photon flux rate in counts/second
%                                   for each channel, after gridding.
%
%   Example Usage:
%   % Assume 'myRawData' is your 43065x561 array (like MCSsample.online)
%   % and 'myTimeGrid' is your desired gridded time vector (like time_grid)
%   % and 'avgIntervalSec' is your fixed interval (like ave_time, e.g., 2 seconds).
%   %
%   % [final_counts, final_flux] = convertPhotonData(myRawData, myTimeGrid, avgIntervalSec);

% --- 1. Remove duplicate time points ---
[~, ia, ~] = unique(dataArray(:,1),'rows'); % Get indices of unique time points
uA = dataArray(ia,:); % Apply those indices to get data with unique time points

% --- 2. Calculate initial flux rate (counts/sec) ---

% Extract the time stamp column (first column) from the unique-timestamped data
time_hours = uA(:, 1);

% Calculate the time differences between consecutive measurements in hours
% This gives the duration over which each row's counts were accumulated.
deltaTime_hours = diff(time_hours);

% Handle the first row: Assume the first interval is the same as the second.
% This ensures `deltaTime_seconds` has the same number of rows as the data.
if ~isempty(deltaTime_hours)
    deltaTime_hours = [deltaTime_hours(1); deltaTime_hours];
else
    % Handle case where there's only one row (no difference possible).
    % In this scenario, flux conversion might not be meaningful.
    deltaTime_hours = zeros(size(time_hours));
    warning('Input dataArray (after unique processing) has only one row. Flux conversion may not be meaningful.');
end

% Convert time differences from hours to seconds (1 hour = 3600 seconds)
deltaTime_seconds = deltaTime_hours * 3600;

% Ensure deltaTime_seconds does not have non-positive values (zero or negative)
% to prevent division by zero or invalid results. Replace with a very small positive number.
min_valid_dt = min(deltaTime_seconds(deltaTime_seconds > 0));
if isempty(min_valid_dt) || min_valid_dt == 0
    min_valid_dt = 1e-9; % Fallback if all intervals are zero or only one data point
end
deltaTime_seconds(deltaTime_seconds <= 0) = min_valid_dt; % Replace non-positive intervals

% Extract the photon counts data (all columns except the first one)
photonCountsData = uA(:, 2:end);

% Initialize an array to store the flux rates
fluxRate_countsPerSec = zeros(size(photonCountsData));

% Calculate flux rate (counts/second) for each data point
% Divide each photon count by its corresponding time interval in seconds using element-wise division (./).
for colIdx = 1:size(photonCountsData, 2)
    fluxRate_countsPerSec(:, colIdx) = photonCountsData(:, colIdx) ./ deltaTime_seconds;
end

% --- 3. Grid the flux rate to a fixed time interval ---
% Uses 'nearest' neighbor interpolation and 'extrap' to handle points outside the original range.
griddedFluxRate_countsPerSec = interp1(time_hours, fluxRate_countsPerSec, time_grid, 'nearest', 'extrap');

% --- 4. Convert gridded flux rate back to total counts in the grid interval ---
% Total Counts = Flux Rate (counts/sec) * Time Interval (sec)
totalCounts_inGridInterval = griddedFluxRate_countsPerSec * fixedGridInterval_seconds;

end
