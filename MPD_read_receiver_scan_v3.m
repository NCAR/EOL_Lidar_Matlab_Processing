clear %all
close all
d=pwd;
%cd('/Volumes/smaug1/rsfdata/MPD/mpd_03_data/2025/20250825')
cd('/Volumes/Macintosh HD/Users/spuler/Downloads/EKO_Scans/20250825/')

% --- USER-DEFINED GRID RESOLUTION ---
time_interval_sec = 1; % Grid resolution in seconds.

% Define the four types of files
file_types = {'ReceiverScanEtalon', 'ReceiverScanLaser', 'ReceiverScanMCS', 'ReceiverScanWavemeter'};
% Initialize a master structure to store all combined data
all_combined_data = struct();
all_combined_timestamps = struct();
for i = 1:length(file_types)
    current_type = file_types{i};
    
    % Get a list of files for the current type
    nc_files = dir([current_type, '_*.nc']);

    % --- NEW DIAGNOSTIC OUTPUT ---
    num_files = length(nc_files);
    fprintf('Found %d files for type: %s\n', num_files, current_type);
    
    % Skip if no files are found for this type
    if isempty(nc_files)
        continue;
    end
    
    % Pre-allocate for timestamps
    timestamps = NaT(length(nc_files), 1);
    file_names = {nc_files.name};
    
    % Parse timestamps from filenames
    for j = 1:length(file_names)
        filename = file_names{j};
        
        % --- Use regexp for robust filename parsing ---
        pattern = '(\d{8})_(\d{6})\.nc';
        tokens = regexp(filename, pattern, 'tokens');
        
        if isempty(tokens)
            warning('Could not parse date and time from filename: %s. Skipping.', filename);
            continue;
        end
        
        date_str = tokens{1}{1};
        time_str = tokens{1}{2};
        
        temp_timestamp = datetime(date_str, 'InputFormat', 'yyyyMMdd') + ...
                        hours(str2double(time_str(1:2))) + ...
                        minutes(str2double(time_str(3:4))) + ...
                        seconds(str2double(time_str(5:6)));
        
        timestamps(j) = temp_timestamp;

    end
    
    % Sort files by time
    [sorted_timestamps, sorted_indices] = sort(timestamps);
    sorted_files = file_names(sorted_indices);
    
    % Store the sorted timestamps in the master structure
    all_combined_timestamps.(current_type) = sorted_timestamps;
    % Read and concatenate data
    temp_combined_data = struct();
    
    if ~isempty(sorted_files)
        
        % Get info for the first file
        first_file_info = ncinfo(sorted_files{1});
        var_names = {first_file_info.Variables.Name};
        
        % Check if the file type requires special handling for channels
        is_mcs_data = strcmp(current_type, 'ReceiverScanMCS');
        % Determine if special handling is possible for this MCS file
        has_mcs_vars = false;
        if is_mcs_data
            if ismember('Channel', var_names) && ismember('ChannelAssignment', var_names) && ismember('Data', var_names)
                has_mcs_vars = true;
            end
        end
        if has_mcs_vars
            % --- NEW DIAGNOSTIC LINE ---
            fprintf('Executing special MCS data handling block.\n');

            % --- SPECIAL HANDLING FOR MCS DATA ---
            
            % Read all data first, as it's not sequentially ordered
            temp_combined_data.data = ncread(sorted_files{1}, 'Data');
            temp_combined_data.channel = ncread(sorted_files{1}, 'Channel');
            temp_combined_data.ChannelAssignment = ncread(sorted_files{1}, 'ChannelAssignment');
            
            % --- CORRECTED: Use filename date and add hours from nc time ---
            file_date = datetime(regexp(sorted_files{1}, '(\d{8})', 'match', 'once'), 'InputFormat', 'yyyyMMdd');
            temp_time_hours = ncread(sorted_files{1}, 'time');
            temp_combined_data.time = file_date + hours(temp_time_hours);
            
            % Concatenate data from remaining files
            for k = 2:length(sorted_files)
                file_name = sorted_files{k};
                
                % --- CORRECTED CONCATENATION DIMENSION ---
                temp_combined_data.data = cat(1, temp_combined_data.data, ncread(file_name, 'Data'));
                temp_combined_data.channel = cat(1, temp_combined_data.channel, ncread(file_name, 'Channel'));
                
                % --- CORRECTED: Use filename date and add hours from nc time ---
                new_file_date = datetime(regexp(file_name, '(\d{8})', 'match', 'once'), 'InputFormat', 'yyyyMMdd');
                new_time_hours = ncread(file_name, 'time');
                new_time = new_file_date + hours(new_time_hours);
                temp_combined_data.time = cat(1, temp_combined_data.time, new_time);
            end
            
            % --- NEW DIAGNOSTIC OUTPUT TO SHOW THE MAPPING ---
            fprintf('--- Channel Assignment to Data Mapping ---\n');
            for k=1:length(temp_combined_data.ChannelAssignment)
                channel_name = strtrim(temp_combined_data.ChannelAssignment(k));
                fprintf('Channel Number %d: "%s"\n', k, channel_name);
            end
            fprintf('----------------------------------------\n');
            
            % Corrected parsing loop to add variables to the existing structure
            
            for k = 1:length(temp_combined_data.ChannelAssignment)
                current_channel_name = strrep(temp_combined_data.ChannelAssignment(k), ' ', '');
                
                % Find all indices in the combined 'channel' variable that match the current channel number (k)
                channel_indices = find(temp_combined_data.channel == k - 1);
                
                if ~isempty(channel_indices)
                    % Check if a field with this name already exists
                    data_var_name = ['data_', char(current_channel_name)];
                    time_var_name = ['time_', char(current_channel_name)];
                    
                    if isfield(temp_combined_data, data_var_name)
                         % If it exists, append a unique suffix to the name to avoid overwriting
                         data_var_name = ['data_', char(current_channel_name), num2str(k)];
                         time_var_name = ['time_', char(current_channel_name), num2str(k)];
                    end
                    
                    temp_combined_data.(data_var_name) = temp_combined_data.data(channel_indices, :, :);
                    temp_combined_data.(time_var_name) = temp_combined_data.time(channel_indices);
                end
            end
            
            % --- ADDED CODE TO OUTPUT FILE SIZES ---
            fprintf('--- File sizes for combined ReceiverScanMCS variables ---\n');
            mcs_vars_to_check = fieldnames(temp_combined_data);
            for m = 1:length(mcs_vars_to_check)
                var_name = mcs_vars_to_check{m};
                var_data = temp_combined_data.(var_name);
                var_size_bytes = numel(var_data) * whos('var_data').bytes / numel(var_data);
                fprintf('Variable: %s | Size: %.2f KB (%.2f MB)\n', var_name, var_size_bytes / 1024, var_size_bytes / (1024*1024));
            end
            fprintf('------------------------------------------------------\n');
            
        else
            % --- STANDARD HANDLING FOR OTHER DATA TYPES ---
            fprintf('Using standard data handling block.\n');
            
            % --- Corrected: Use filename date and add hours from nc time for first file ---
            file_date = datetime(regexp(sorted_files{1}, '(\d{8})', 'match', 'once'), 'InputFormat', 'yyyyMMdd');
            if ismember('time', var_names)
                temp_time_hours = ncread(sorted_files{1}, 'time');
                temp_combined_data.time = file_date + hours(temp_time_hours);
            end
            
            for k = 1:length(first_file_info.Variables)
                var_name = first_file_info.Variables(k).Name;
                if ~strcmp(var_name, 'time')
                    temp_combined_data.(var_name) = ncread(sorted_files{1}, var_name);
                end
            end
            
            % Concatenate data from remaining files
            for k = 2:length(sorted_files)
                file_name = sorted_files{k};
                current_file_info = ncinfo(file_name);
                
                % --- Corrected: Use filename date and add hours from nc time for remaining files ---
                new_file_date = datetime(regexp(file_name, '(\d{8})', 'match', 'once'), 'InputFormat', 'yyyyMMdd');
                if ismember('time', {current_file_info.Variables.Name})
                    new_time_hours = ncread(file_name, 'time');
                    new_time = new_file_date + hours(new_time_hours);
                    temp_combined_data.time = cat(1, temp_combined_data.time, new_time);
                end
                
                for l = 1:length(current_file_info.Variables)
                    var_name = current_file_info.Variables(l).Name;
                    
                    % We'll skip variables that don't have a time dimension.
                    if strcmp(var_name, 'ChannelAssignment') || strcmp(var_name, 'time')
                        continue;
                    end
                    
                    if isfield(temp_combined_data, var_name)
                        new_data = ncread(file_name, var_name);
                        
                        time_dim_index = -1;
                        if isfield(current_file_info.Variables(l), 'Dimensions')
                            for m = 1:length(current_file_info.Variables(l).Dimensions)
                                if strcmp(current_file_info.Variables(l).Dimensions(m).Name, 'time')
                                    time_dim_index = m;
                                    break;
                                end
                            end
                        end
                        
                        if time_dim_index ~= -1
                            temp_combined_data.(var_name) = cat(time_dim_index, temp_combined_data.(var_name), new_data);
                        else
                            warning('Variable %s or its "time" dimension not found. Skipping.', var_name);
                        end
                    end
                end
            end
        end
    end
    
    % Store the combined data for the current type
    all_combined_data.(current_type) = temp_combined_data;
end


%% --- NEW CODE FOR GRIDDING ---
fprintf('\nBeginning data gridding to a common time base...\n');

% 1. Create a common master time grid
% Get the min/max timestamps from all data types
field_names = fieldnames(all_combined_timestamps);
all_timestamps = cell(1, length(field_names));
for i = 1:length(field_names)
    all_timestamps{i} = all_combined_timestamps.(field_names{i});
end
min_time = min(cellfun(@(x) min(x), all_timestamps));
max_time = max(cellfun(@(x) max(x), all_timestamps));
master_time_grid = (min_time:seconds(time_interval_sec):max_time)';

% Initialize a new structure for the gridded data
gridded_data = struct();
gridded_data.master_time_grid = master_time_grid;

% 2. Loop through all combined data and interpolate
data_types_to_grid = {'ReceiverScanMCS', 'ReceiverScanWavemeter'};

for i = 1:length(data_types_to_grid)
    current_data_type = data_types_to_grid{i};
    
    gridded_data.(current_data_type) = struct();
    
    variables_to_grid = fieldnames(all_combined_data.(current_data_type));
    
    for j = 1:length(variables_to_grid)
        var_name = variables_to_grid{j};
        
        % We'll only grid variables that have a corresponding time variable.
        original_time_var_name = '';
        if contains(var_name, 'data_') && isfield(all_combined_data.(current_data_type), ['time_', var_name(6:end)])
            original_time_var_name = ['time_', var_name(6:end)];
        elseif strcmp(var_name, 'Wavelength') && isfield(all_combined_data.(current_data_type), 'time')
            original_time_var_name = 'time';
        else
            continue; % Skip variables that don't need gridding
        end
        
        original_time = all_combined_data.(current_data_type).(original_time_var_name);
        original_data = all_combined_data.(current_data_type).(var_name);
        
        % --- ADDED DIAGNOSTIC OUTPUT FOR TIME RANGES ---
        fprintf('Processing variable: %s in %s\n', var_name, current_data_type);
        fprintf('Original Time Range: %s to %s\n', datestr(min(original_time)), datestr(max(original_time)));
        fprintf('Master Time Grid Range: %s to %s\n', datestr(min(master_time_grid)), datestr(max(master_time_grid)));

        if isnumeric(original_data) || islogical(original_data)
            % Interpolate the data. Add 'extrap' to fill out-of-range values
            interpolated_data = interp1(datenum(original_time), original_data, datenum(master_time_grid), 'nearest');

            % --- NEW DIAGNOSTIC CHECK FOR ALL-ZEROS ---
            if all(interpolated_data(:) == 0)
                warning('Gridded data for "%s" is all zeros. This may be due to a time range mismatch.', var_name);
            end

            % Store the interpolated data in the new sub-structure
            gridded_data.(current_data_type).(var_name) = interpolated_data;
        else
            warning('Skipping variable "%s" in "%s" as it is not a numeric data type.', var_name, current_data_type);
        end
    end
end
% Store the final gridded structure in the master structure
all_combined_data.GriddedData = gridded_data;

fprintf('Data gridding complete. Gridded data stored in all_combined_data.GriddedData.\n');

cd(d)



%% --- CODE for Plotting select time periods of the scans ---

start_time = datetime('22 Aug 2025 18:0', 'InputFormat', 'dd MMM yyyy HH:mm')
end_time = datetime('22 Aug 2025 18:10', 'InputFormat', 'dd MMM yyyy HH:mm')

start_time = datetime('22 Aug 2025 18:17', 'InputFormat', 'dd MMM yyyy HH:mm')
end_time = datetime('22 Aug 2025 18:24', 'InputFormat', 'dd MMM yyyy HH:mm')
 
 start_time = datetime('22 Aug 2025 22:51', 'InputFormat', 'dd MMM yyyy HH:mm')
 end_time = datetime('22 Aug 2025 23:04', 'InputFormat', 'dd MMM yyyy HH:mm')
 
 start_time = datetime('25 Aug 2025 18:14', 'InputFormat', 'dd MMM yyyy HH:mm')
 end_time = datetime('25 Aug 2025 18:26', 'InputFormat', 'dd MMM yyyy HH:mm')

% start_time = datetime('26 Aug 2025 17:19', 'InputFormat', 'dd MMM yyyy HH:mm')
% end_time = datetime('26 Aug 2025 17:55', 'InputFormat', 'dd MMM yyyy HH:mm')

 start_time = datetime('25 Aug 2025 01:11:30', 'InputFormat', 'dd MMM yyyy HH:mm:ss')
 end_time = datetime('25 Aug 2025 01:14:15', 'InputFormat', 'dd MMM yyyy HH:mm:ss')
 check_time = datetime('25 Aug 2025 01:12:00', 'InputFormat', 'dd MMM yyyy HH:mm:ss')

%check_time= start_time;


figure(1)
semilogy(all_combined_data.ReceiverScanMCS.time_O2OfflineMol, sum(all_combined_data.ReceiverScanMCS.data_O2OfflineMol(:,7:end),2),'bo')
hold on
semilogy(all_combined_data.ReceiverScanMCS.time_O2OfflineComb, sum(all_combined_data.ReceiverScanMCS.data_O2OfflineComb(:,7:end),2),'ko')
semilogy(all_combined_data.ReceiverScanMCS.time_WVOnline, sum(all_combined_data.ReceiverScanMCS.data_WVOnline(:,7:end),2),'ro')
semilogy(all_combined_data.ReceiverScanMCS.time_WVOffline, sum(all_combined_data.ReceiverScanMCS.data_WVOffline(:,7:end),2),'go')
hold off
grid on
xlim([start_time, end_time])

[~, start_index] = min(abs(master_time_grid - start_time))
[~, end_index] = min(abs(master_time_grid - end_time))
[~, check_index] = min(abs(master_time_grid - check_time))


figure(2)
plot(all_combined_data.GriddedData.ReceiverScanWavemeter.Wavelength(start_index:end_index), sum(all_combined_data.GriddedData.ReceiverScanMCS.data_O2OfflineMol(start_index:end_index,7:end),2),'bo')
hold on
plot(all_combined_data.GriddedData.ReceiverScanWavemeter.Wavelength(start_index:end_index), sum(all_combined_data.GriddedData.ReceiverScanMCS.data_WVOnline(start_index:end_index,7:end),2),'ro')
plot(all_combined_data.GriddedData.ReceiverScanWavemeter.Wavelength(start_index:end_index), sum(all_combined_data.GriddedData.ReceiverScanMCS.data_O2OfflineComb(start_index:end_index,7:end),2),'ko')
plot(all_combined_data.GriddedData.ReceiverScanWavemeter.Wavelength(start_index:end_index), sum(all_combined_data.GriddedData.ReceiverScanMCS.data_WVOffline(start_index:end_index,7:end),2),'go')
hold off
grid on
xlim([770.08, 770.13])
xlim([770.105, 770.113])
%xlim([828.1, 828.4])

figure(4)
semilogy(all_combined_data.GriddedData.ReceiverScanWavemeter.Wavelength(start_index:end_index), sum(all_combined_data.GriddedData.ReceiverScanMCS.data_O2OfflineMol(start_index:end_index,7:end),2),'bo')
hold on
semilogy(all_combined_data.GriddedData.ReceiverScanWavemeter.Wavelength(start_index:end_index), sum(all_combined_data.GriddedData.ReceiverScanMCS.data_WVOnline(start_index:end_index,7:end),2),'ro')
semilogy(all_combined_data.GriddedData.ReceiverScanWavemeter.Wavelength(start_index:end_index), sum(all_combined_data.GriddedData.ReceiverScanMCS.data_O2OfflineComb(start_index:end_index,7:end),2),'ko')
semilogy(all_combined_data.GriddedData.ReceiverScanWavemeter.Wavelength(start_index:end_index), sum(all_combined_data.GriddedData.ReceiverScanMCS.data_WVOffline(start_index:end_index,7:end),2),'go')
hold off
grid on
xlim([770.08, 770.13])
%xlim([828.175, 828.32])
%xlim([828.28, 828.33])


figure(10)
%plot(all_combined_data.GriddedData.ReceiverScanMCS.data_O2OfflineMol(check_index,:),'b')
plot(all_combined_data.GriddedData.ReceiverScanMCS.data_WVOnline(check_index,:),'r')
hold on
plot(all_combined_data.GriddedData.ReceiverScanMCS.data_WVOffline(check_index,:),'g')
%plot(all_combined_data.GriddedData.ReceiverScanMCS.data_O2OfflineComb(check_index,:),'k')
hold off
grid on