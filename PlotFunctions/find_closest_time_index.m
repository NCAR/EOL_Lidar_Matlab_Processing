function [idx_to_plot] = find_closest_time_index(time_axis_data, target_time)
% FIND_CLOSEST_TIME_INDEX Finds the index of the time point closest to the target_time.
    [~, idx_to_plot] = min(abs(time_axis_data - target_time));
end
