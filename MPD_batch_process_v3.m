%% MPD_batch_process_v3.m: Refactored Batch Processor

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

% Extract control variables for cleaner loop structure
start_day = config.dates.start_day;
stop_day = config.dates.stop_day;
systems = config.systems_to_process;


% --- 2. Daily Processing Loop (Modular Execution) ---

if config.flags.process == 1
    
    disp(['Starting daily processing from ', config.dates.start_str, ' to ', config.dates.stop_str]);
    
    tStart = tic; 
    
    % Iterate over each day
    for current_day_num = start_day:stop_day
        
        file_date_str = datestr(current_day_num, 'yyyymmdd');
        disp(['Processing Day: ', file_date_str]);

        % Iterate over each predefined system/channel job
        for i = 1:length(systems)
            job = systems{i};
            
            % Create a temporary config structure specific to this job run
            current_job_config = config;
            current_job_config.job = job;

            fprintf('    -> Running %s: %s, %s...\n', job.node, job.channels, job.correction);

            % *** NEW CALL: Use the refactored orchestrator ***
            MPD_run_daily_analysis( ...
                job.channels, ...
                job.correction, ...
                job.node, ...
                file_date_str, ...
                current_job_config);
            % *************************************************

        end
    end
    
    tElapsed = toc(tStart);
    disp(['Daily processing complete. Total time: ', num2str(tElapsed), ' seconds.']);

end


% --- 3. Multi-Day Plotting Loop ---

if config.flags.plot_multiday == 1
    
    % Iterate over each predefined plotting job
    for i = 1:length(config.multiday_plots)
        plot_job = config.multiday_plots{i};
        
        disp(['Starting multi-day plot for ', plot_job.node, ...
              ' from ', plot_job.start_date, ' to ', plot_job.stop_date]);
              
        % *** NEW CALL: Use the refactored plotting utility ***
        MPD_plot_utility( ...
            config.flags.save_quicklook, ... 
            config.flags.save_data, ...  
            plot_job.near, ...
            plot_job.afterpulse, ...
            plot_job.node, ...
            plot_job.start_date, ...
            plot_job.stop_date, ...
            plot_job.skip);
        % *************************************************

    end
    
end