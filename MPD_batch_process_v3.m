%% MPD_batch_process_v3.m: Independent Batch Processor

% Add path to utilities and project folder structure (MUST be set correctly)
addpath('/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing')

% Define the configuration function name
MPD_config_script = 'MPD_config'; 

% --- 1. Handle User Input and Load Configuration ---
if strcmp(getenv('HOSTNAME'),'eol-smaug.eol.ucar.edu') == 1 
  prompt_start = 'Enter start date (e.g., 20230411): ';
  start_date = input(prompt_start, 's');
  prompt_stop = 'Enter stop date (e.g., 20230411): ';
  stop_date = input(prompt_stop, 's');
  config = feval(MPD_config_script, start_date, stop_date);
else
  config = feval(MPD_config_script);
end

% Check for new processing flags (assuming you update MPD_config.m)
if ~isfield(config.flags, 'process_data'); config.flags.process_data = config.flags.process; end
if ~isfield(config.flags, 'run_plots'); config.flags.run_plots = config.flags.plot_multiday; end


% --- 2. Independent Data Processing Loop (Slow: Creates Data Files) ---

if config.flags.process_data == 1
    
    disp(['Starting DATA PROCESSING from ', config.dates.start_str, ' to ', config.dates.stop_str]);
    
    tStart = tic; 
    start_day = config.dates.start_day;
    stop_day = config.dates.stop_day;

    for current_day_num = start_day:stop_day
        
        file_date_str = datestr(current_day_num, 'yyyymmdd');
        disp(['Processing Day: ', file_date_str]);

        for i = 1:length(config.systems_to_process)
            job = config.systems_to_process{i};
            
            current_job_config = config;
            current_job_config.job = job;

            fprintf('    -> Running %s: %s, %s...\n', job.node, job.channels, job.correction);

            MPD_run_daily_analysis( ...
                job.channels, ...
                job.correction, ...
                job.node, ...
                file_date_str, ...
                current_job_config);
        end
    end
    
    tElapsed = toc(tStart);
    disp(['Data processing complete. Total time: ', num2str(tElapsed), ' seconds.']);

end


% --- 3. Independent Plotting Loop (Fast: Reads Data Files) ---

if config.flags.run_plots == 1
    
    disp('Starting INDEPENDENT PLOTTING (Reading saved .mat files)...');
    tStartPlot = tic; 

    % Iterate over each predefined plotting job
    for i = 1:length(config.multiday_plots)
        plot_job = config.multiday_plots{i};
        
        disp(['Plotting Job: ', plot_job.node, ...
              ' from ', plot_job.start_date, ' to ', plot_job.stop_date]);
              
        MPD_plot_utility( ...
            config.flags.save_quicklook, ... 
            config.flags.save_data, ...  
            plot_job.near, ...
            plot_job.afterpulse, ...
            plot_job.node, ...
            plot_job.start_date, ...
            plot_job.stop_date, ...
            plot_job.skip);
    end
    
    tElapsedPlot = toc(tStartPlot);
    disp(['Plotting complete. Total time: ', num2str(tElapsedPlot), ' seconds.']);

end

% --- CRITICAL: Update your MPD_config.m to enable independent control ---
% Please update MPD_config.m to use the new flags.

% In MPD_config.m, update the flags section:
% config.flags.process_data = 0;   % Set this to 0 when you ONLY want plots
% config.flags.run_plots = 1;      % Set this to 1 when you ONLY want plots