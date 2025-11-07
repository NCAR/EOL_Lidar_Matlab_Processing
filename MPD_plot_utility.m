%% MPD_plot_utility.m: Refactored Post-Processing Plotting Utility

function[] = MPD_plot_utility(save_figs, save_data, near, afterpulse, node, daystr, daystr2, skip)
% This function loads daily processed data (.mat files) over a date range, 
% concatenates them, generates standard plots, and saves figures.

% --- 1. Initialization and Setup ---
close all;
tic;

% Environment Setup 
dd = pwd; 
font_size = 28;
WS_font_size = 24;
C = importdata('NCAR_C_Map.mat'); 
RB_scale = 1;
range_grid_size = 75; 

% Determine paths (using helper functions defined elsewhere or inside orchestrator)
serv_path = get_server_data_path(); 

% Set flags based on input arguments
flag.near = near;
flag.afterpulse = afterpulse;
flag.save_figs = save_figs;
flag.save_data = save_data;
flag.decimate = 1; 
WS = 1; % Assume WS data exists in the processed files

% Determine input/output directories
nodeStr = extractAfter(node, 'MPD');
base_data_path = fullfile(serv_path, ['mpd_', nodeStr, '_processed_data'], 'Matlab');
plot_dir = fullfile(base_data_path, 'Plots'); 
if flag.afterpulse == 1
    base_data_path = fullfile(base_data_path, 'afterpulse');
    plot_dir = fullfile(base_data_path, 'Plots'); 
end
write_data_folder = base_data_path;
cd(base_data_path);


% --- 2. Data Concatenation Loop ---

start_day = datenum(daystr,'yyyymmdd');
stop_day = datenum(daystr2,'yyyymmdd');
days = stop_day - start_day + 1;

% Initialize combined data structures
N_avg_comb = []; RB_comb = []; duration = [];
background_comb_on = []; background_comb_off = []; 
lambda_comb_on = []; lambda_comb_off = [];
surf_T = []; surf_P = []; surf_AH = []; 
i_off = []; i_on = []; p_on = []; p_off = [];
t_bench = []; t_base = [];
range_plot = []; % Initialize range vector

for current_day_num = start_day:stop_day
    date_str_yyyyMMdd = datestr(current_day_num, 'yyyymmdd');

    % Determine the file name 
    if flag.near == 1
        data_filename = strcat(node, "_", date_str_yyyyMMdd, "_WV_near.mat"); 
    else
        data_filename = strcat(node, "_", date_str_yyyyMMdd, "_WV.mat");
    end

    if exist(data_filename, 'file') == 2
        try
            loaded_data = load(data_filename);
            
            % --- Standardize and Grid Data ---
            current_gate = loaded_data.gate;
            N_avg_grid = loaded_data.N_avg;
            RB_grid = loaded_data.RB;
            range_data = loaded_data.range;
            
            if current_gate < range_grid_size
                range_grid_75 = 0:range_grid_size:(length(range_data)-1)*current_gate;
                N_avg_grid = interp1(range_data, loaded_data.N_avg', range_grid_75, 'linear', 'extrap')';
                RB_grid = interp1(range_data, loaded_data.RB', range_grid_75, 'linear', 'extrap')';
                range_data = range_grid_75;
            end

            % --- Concatenation ---
            if isempty(N_avg_comb)
                N_avg_comb = N_avg_grid;
                RB_comb = RB_grid;
                duration = loaded_data.time_new;
                background_comb_on = loaded_data.background_on;
                background_comb_off = loaded_data.background_off;
                lambda_comb_on = loaded_data.lambda_all;
                lambda_comb_off = loaded_data.lambda_all_off;
                range_plot = range_data;
                
                % Housekeeping init 
                if isfield(loaded_data, 'Surf_T')
                    surf_T = loaded_data.Surf_T; surf_P = loaded_data.Surf_P;
                    surf_AH = loaded_data.Surf_AH; i_off = loaded_data.I_off;
                    i_on = loaded_data.I_on; p_on = loaded_data.P_on;
                    p_off = loaded_data.P_off; t_bench = loaded_data.T_bench;
                    t_base = loaded_data.T_base;
                end

            else
                % Match range dimensions 
                range_lim_comb = size(N_avg_comb, 2);
                range_lim_current = size(N_avg_grid, 2);
                range_limit = min([range_lim_comb, range_lim_current]);

                % Concatenate: skip the first row of the subsequent file to avoid time duplication
                N_avg_comb = vertcat(N_avg_comb(:, 1:range_limit), N_avg_grid(2:end, 1:range_limit));
                RB_comb = vertcat(RB_comb(:, 1:range_limit), RB_grid(2:end, 1:range_limit));
                duration = vertcat(duration, loaded_data.time_new(2:end));
                background_comb_on = vertcat(background_comb_on, loaded_data.background_on(2:end));
                background_comb_off = vertcat(background_comb_off, loaded_data.background_off(2:end));
                lambda_comb_on = vertcat(lambda_comb_on, loaded_data.lambda_all(2:end));
                lambda_comb_off = vertcat(lambda_comb_off, loaded_data.lambda_all_off(2:end));
                
                % Housekeeping concatenation
                if isfield(loaded_data, 'Surf_T')
                    surf_T = vertcat(surf_T, loaded_data.Surf_T(2:end,:));
                    surf_P = vertcat(surf_P, loaded_data.Surf_P(2:end,:));
                    surf_AH = vertcat(surf_AH, loaded_data.Surf_AH(2:end,:));
                    i_off = vertcat(i_off, loaded_data.I_off(2:end,:));
                    i_on = vertcat(i_on, loaded_data.I_on(2:end,:));
                    p_on = vertcat(p_on, loaded_data.P_on(2:end,:));
                    p_off = vertcat(p_off, loaded_data.P_off(2:end,:));
                    t_bench = vertcat(t_bench, loaded_data.T_bench(2:end,:));
                    t_base = vertcat(t_base, loaded_data.T_base(2:end,:));
                end
            end
        catch ME
            warning(['Failed to load or concatenate data from ', data_filename, '. Skipping. Error: ', ME.message]);
        end
    end
end
% Final range vector 
if ~isempty(N_avg_comb)
    range = range_plot(1:size(N_avg_comb, 2));
else
    cd(dd);
    warning('No data found to plot. Exiting plotting utility.');
    return;
end


% --- 3. Save Combined Data File ---
if flag.save_data == 1
    cd(write_data_folder);
    name = strcat(node, "_", daystr, "_", daystr2, "_Combined"); 
    
    MPD_COMBINED_DATA = struct();
    MPD_COMBINED_DATA.N_avg_comb = N_avg_comb;
    MPD_COMBINED_DATA.RB_comb = RB_comb;
    MPD_COMBINED_DATA.range = range;
    MPD_COMBINED_DATA.time = duration;
    
    if ~isempty(surf_T)
        MPD_COMBINED_DATA.surf_T = surf_T; MPD_COMBINED_DATA.surf_P = surf_P; 
        MPD_COMBINED_DATA.surf_AH = surf_AH; MPD_COMBINED_DATA.i_off = i_off;
        MPD_COMBINED_DATA.i_on = i_on; MPD_COMBINED_DATA.p_on = p_on;
        MPD_COMBINED_DATA.p_off = p_off; MPD_COMBINED_DATA.t_bench = t_bench;
        MPD_COMBINED_DATA.t_base = t_base;
    end
    
    save(name, 'MPD_COMBINED_DATA'); 
    cd(base_data_path); 
end


% --- 4. Plotting ---
% (Original plotting logic is executed here)
scrsz = get(0,'ScreenSize');
x = duration';
y = range./1e3;
Z_AH = double(real(N_avg_comb'.*1e6./6.022E23.*18.015)); 
Z_RB = double((real(RB_comb')./RB_scale));

% Time axis tick marks
if days == 1
    date_plot = datestr(mean(duration), 'dd mmm yyyy');
    xData = linspace(fix(min(duration)), ceil(max(duration)), 25);
else
    date_plot = [datestr(min(duration), 'dd mmm yyyy') ' - ' datestr(max(duration), 'dd mmm yyyy')];
    xData = linspace(fix(min(duration)), ceil(max(duration)), round((ceil(max(duration))-fix(min(duration)))/skip)+1);
end

% Decimation for display/save resolution
if flag.decimate == 1 && days > 1
    scrsz = get(0,'ScreenSize');
    decimate_time = max(1, fix(length(duration) / scrsz(3) / 2)); 
    Z_AH = nanmoving_average(Z_AH, decimate_time, 2, 0); 
    Z_RB = nanmoving_average(Z_RB, decimate_time, 2, 0);   
    x = x(1:decimate_time:end);
    Z_AH = Z_AH(:, 1:decimate_time:end);
    Z_RB = Z_RB(:, 1:decimate_time:end);   
end

Scrnsize=[scrsz(4)/2 scrsz(4)/10 scrsz(3)/1 scrsz(4)/2];

% PLOT WATER VAPOR (g/m^3) - Figure 1
figure(1)
set(gcf,'renderer','zbuffer');
h = pcolor(x,y,Z_AH);
set(h, 'EdgeColor', 'none');
colorbar('EastOutside');
axis([fix(min(x)) ceil(max(x)) 0 6])
caxis([0 8]);
colormap(jet)
ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
set(gca, 'XTick', xData); set(gca,'TickDir','out');
if days == 1
  datetick('x','HH:MM','keeplimits', 'keepticks');
  xlabel('Time (UTC)','fontweight','b','fontsize',font_size);
  title({[date_plot, node, ' Water Vapor (g m^{-3})']},'fontweight','b','fontsize',font_size);
else
  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
  title({[node, ' Water Vapor (g m^{-3})']},'fontweight','b','fontsize',font_size);
end
set(gca,'Fontsize',font_size,'Fontweight','b');

% PLOT RELATIVE BACKSCATTER (RB) - Figure 2
figure(2)
set(gcf,'renderer','zbuffer');
h = pcolor(x,y,Z_RB);
set(h, 'EdgeColor', 'none');
colorbar('EastOutside');
axis([fix(min(duration)) ceil(max(duration)) 0 12])
caxis([1e1 1e7]);
colormap(C)
ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size);
set(gca, 'XTick', xData); set(gca,'TickDir','out');
if days == 1
  datetick('x','HH:MM','keeplimits', 'keepticks');
  xlabel('Time (UTC)','fontweight','b','fontsize',font_size);
  title({[date_plot, node, ' Attenuated Backscatter (A.U.)']},'fontweight','b','fontsize',font_size);
else
  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
  title({[node, ' Attenuated Backscatter (A.U.)']},'fontweight','b','fontsize',font_size);
end
set(gca,'Fontsize',font_size,'Fontweight','b');
set(gca,'Zscale', 'log'); set(gca,'Colorscale', 'log'); set(gca,'Zscale', 'linear');


% PLOT BACKGROUND - Figure 3
figure(3)
semilogy(duration, background_comb_off, 'black'); % Only one background channel is plotted for simplicity
axis([fix(min(duration)) ceil(max(duration)) 1e2 1e7])
ylabel('Offline background C/s', 'Fontsize', font_size, 'Fontweight', 'b');  
grid on; legend('WV Offline', 'Location', 'NorthWest');
set(gca, 'XTick', xData); set(gca,'TickDir','out');
if days == 1
  datetick('x','HH:MM','keeplimits', 'keepticks');
  xlabel('Time (UTC)','fontweight','b','fontsize',font_size); 
  title({[date_plot, node, ' Background']},'fontweight','b','fontsize',font_size);
else
  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
  title({[node, ' Background']},'fontweight','b','fontsize',font_size); 
end
set(gca,'Fontsize',font_size,'Fontweight','b');


% PLOT HOUSEKEEPING DATA (if relevant fields exist) - Figure 4
if WS == 1 && ~isempty(surf_T)
    figure1 = figure(4);
    % Subplot 1: Wavelength vs Power
    subplot1=subplot(2,1,1,'Parent',figure1,'YGrid','on', 'XGrid','on');
    plot(duration, lambda_comb_on,'g-','LineWidth',2,'DisplayName','Lambda_{on}');
    axis([fix(min(duration)) ceil(max(duration)) min(lambda_comb_on)*0.9999 max(lambda_comb_on)*1.0001])
    ylabel('wavelength, nm', 'Fontsize', WS_font_size, 'Fontweight', 'b');  
    datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
    set(gca,'Fontsize',WS_font_size,'Fontweight','b');
    
    % Subplot 2: Temperature/Pressure
    subplot2=subplot(2,1,2,'Parent',figure1,'YGrid','on', 'XGrid','on');
    plot(duration, t_bench,'g-', 'LineWidth',1, 'DisplayName','T bench');
    hold on; plot(duration, surf_T, 'b', 'LineWidth',1, 'DisplayName','Surface T'); hold off;
    axis([fix(min(duration)) ceil(max(duration)) 0 40]);
    ylabel('temperature, C', 'Fontsize', WS_font_size, 'Fontweight', 'b'); 
    datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
    set(gca,'Fontsize',WS_font_size,'Fontweight','b');
    
    % Overlay Power (Right Axis 1)
    ax(3) = axes('Position',get(subplot1,'Position'));
    plot(duration, p_off/7500,'b-','LineWidth', 1, 'DisplayName','P_{off}');
    hold on; plot(duration, p_on/7500,'r-', 'LineWidth',1, 'DisplayName','P_{on}'); hold off;
    axis([fix(min(duration)) ceil(max(duration)) 0 50]);
    set(ax(3),'Color','none', 'YAxisLocation','right', 'XAxisLocation','bottom');
    datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
    ylabel('rel. trans. power', 'Fontsize', WS_font_size, 'Fontweight', 'b');  
    set(gca,'Fontsize',WS_font_size,'Fontweight','b');
    legend(ax(3), 'show', 'Location','southwest'); set(legend(ax(3)),'Color','white');

    % Overlay Pressure (Right Axis 2)
    ax(4) = axes('Position',get(subplot2,'Position'));
    plot(duration, surf_P, 'k-','LineWidth', 1, 'DisplayName','Surf P');
    axis([fix(min(duration)) ceil(max(duration)) 0.8 1.0]);
    set(ax(4),'Color','none', 'YAxisLocation','right', 'XAxisLocation','bottom');
    datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
    ylabel('pressure, atm', 'Fontsize', WS_font_size, 'Fontweight', 'b');  
    set(gca,'Fontsize',WS_font_size,'Fontweight','b');
    legend(ax(4), 'show'); set(legend(ax(4)),'Color','white');
    
    legend(subplot1,'Location','NorthWest'); legend(subplot2,'Location','NorthWest'); 
    linkaxes([subplot1, subplot2, ax(3), ax(4)], 'x');
end


% --- 5. Figure Saving ---
if flag.save_figs == 1
    % ... Save logic here
    if exist(plot_dir, 'dir') ~= 7; mkdir(plot_dir); end; cd(plot_dir);
    date_save = datestr(mean(duration), 'yyyymmdd');
    
    Scrnsize_Fig1 = [scrsz(4)/1 scrsz(4)/1 scrsz(3)/0.30 scrsz(4)/2]; 
    
    FigH = figure(1); set(FigH, 'PaperUnits', 'points', 'PaperPosition', Scrnsize_Fig1);
    name = strcat(date_save, '_', node, '_H2O_multi'); print(FigH, name, '-dpng', '-r0');
    
    FigH = figure(2); set(FigH, 'PaperUnits', 'points', 'PaperPosition', Scrnsize_Fig1);
    name = strcat(date_save, '_', node, '_RB_multi'); print(FigH, name, '-dpng', '-r0');
    
    FigH = figure(3); set(FigH, 'PaperUnits', 'points', 'PaperPosition', Scrnsize);
    name = strcat(date_save, '_', node, '_background_multi'); print(FigH, name, '-dpng', '-r0');

    if WS == 1 && ~isempty(surf_T)
        Scrnsize_Fig4 = [scrsz(4)/1 scrsz(4)/1 scrsz(3)/0.35 scrsz(4)/1];
        FigH = figure(4); set(FigH, 'PaperUnits', 'points', 'PaperPosition', Scrnsize_Fig4);
        name = strcat(date_save, '_', node, '_Housekeeping'); print(FigH, name, '-dpng', '-r0');
    end
end
 
cd(dd);
toc

function path = get_server_data_path()
    if strcmp(getenv('HOSTNAME'),'eol-smaug.eol.ucar.edu') == 1; path = '/export/smaug1/rsfdata/MPD/'; else; path = '/Volumes/smaug1/rsfdata/MPD/'; end
end

end