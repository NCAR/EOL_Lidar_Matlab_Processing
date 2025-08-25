function [totalCounts_inGridInterval, griddedFluxRate_countsPerSec] = grid_photon_counts_2D_v2(dataArray, time_grid, fixedGridInterval_seconds)
% grid_photon_counts_2D_revised Processes raw photon count data by binning it to a fixed time interval grid.
%
%   [totalCounts_inGridInterval, griddedFluxRate_countsPerSec] = grid_photon_counts_2D_revised(dataArray, time_grid, fixedGridInterval_seconds)
%
%   Inputs:
%     dataArray: A N x M array where N is the number of data points and M is
%                the number of channels + 1. The first column (dataArray(:,1))
%                must be time stamps in fractional hours. The remaining columns
%                (dataArray(:,2:end)) are photon counts.
%     time_grid: A vector of time points (in fractional hours) defining the
%                new fixed time interval grid.
%     fixedGridInterval_seconds: The fixed time interval in seconds corresponding
%                                to the 'time_grid' (e.g., 2 seconds).
%
%   Outputs:
%     totalCounts_inGridInterval: An P x (M-1) array (where P is length of time_grid)
%                                 containing the total photon counts within the
%                                 specified fixedGridInterval_seconds for each channel.
%     griddedFluxRate_countsPerSec: An P x (M-1) array (where P is length of time_grid)
%                                   containing the photon flux rate in counts/second
%                                   for each channel.

% --- 1. Remove duplicate time points ---
[~, ia, ~] = unique(dataArray(:,1),'rows');
uA = dataArray(ia,:);
time_hours = uA(:, 1);
photonCountsData = uA(:, 2:end);

% --- 2. Discretize timestamps into the new time grid ---
% 'discretize' assigns each original time stamp to a bin on the new grid.
% The time_grid is a vector of bin edges, so we need to add one extra edge.
time_grid_edges = [time_grid; time_grid(end) + (time_grid(2) - time_grid(1))];
bin_indices = discretize(time_hours, time_grid_edges);

% --- 3. Bin and sum the photon counts for each channel ---
num_channels = size(photonCountsData, 2);
num_grid_points = length(time_grid);
totalCounts_inGridInterval = zeros(num_grid_points, num_channels);

% Use 'accumarray' to efficiently sum all counts that fall into each bin.
% This is a vectorized and robust method that avoids loops.
for colIdx = 1:num_channels
    % The third argument [num_grid_points, 1] ensures the output array has the
    % correct size, even if some bins are empty.
    totalCounts_inGridInterval(:, colIdx) = accumarray(bin_indices, photonCountsData(:, colIdx), [num_grid_points, 1], @sum, 0);
end

% --- 4. Calculate the gridded flux rate ---
% The flux is simply the total counts in each grid interval divided by the interval duration.
griddedFluxRate_countsPerSec = totalCounts_inGridInterval / fixedGridInterval_seconds;
end