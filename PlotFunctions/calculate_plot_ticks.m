function [xData, xData_m] = calculate_plot_ticks(x)
% CALCULATE_PLOT_TICKS Calculates dynamic major and minor tick marks for the time axis.
    
    % Target 8 major ticks for the main time span
    TARGET_TICKS = 8; 
    time_span_days = ceil(max(x) - min(x));
    days_per_tick = max(1, ceil(time_span_days / TARGET_TICKS));

    % Ensure days_per_tick is a readable interval
    if days_per_tick <= 2
        days_per_tick = 1; % Daily
    elseif days_per_tick <= 4
        days_per_tick = 3; % Every 3 days
    elseif days_per_tick <= 10
        days_per_tick = 7; % Weekly
    elseif days_per_tick <= 20
        days_per_tick = 14; % Bi-weekly
    elseif days_per_tick <= 60
        days_per_tick = 30; % Monthly
    else
        days_per_tick = ceil(days_per_tick / 30) * 30; % Round up to the nearest month
    end
    
    % Use the calculated step size to generate an exact date series.
    start_date = floor(min(x));
    end_date = ceil(max(x)); 
    
    % Generate major ticks
    xData = start_date : days_per_tick : end_date;
    % Generate minor ticks at half the major interval
    xData_m = start_date : days_per_tick/2 : end_date;

    % Remove ticks that fall outside the display range
    xData(xData > end_date) = [];
    xData_m(xData_m > end_date) = [];
    
    disp(['Dynamic Major Tick Interval: ', num2str(days_per_tick), ' days.']);
end