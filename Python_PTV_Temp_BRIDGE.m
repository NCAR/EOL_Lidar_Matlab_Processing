
clear all; close all

skip = 1
node = 'MPD04';
addpath '/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing/matplotlib/';

serv_path = '/Volumes/eol/sci/mhayman';
cd(strcat(serv_path,'/DIAL/Processed_Data/BRIDGE_2025/ptv0.2/'))
d_read_data = pwd; % get the current path
plot_path = '/Users/spuler/Desktop';
cd(strcat(plot_path,'/mpd/Plots'))
d_save_data = pwd; %set the plot save path
flag.save_data = 0;  %save data at end of processing (0=off 1=on)
low_range_mask = 0;
cd(d_read_data);

[Pythonfilename, Pythondir] = uigetfile('*.*','Select the sonde file', 'MultiSelect', 'on');
jj=1;
  
 variable{1} = 'time';
 variable{2} = 'range';
 variable{3} = 'Temperature_PTV';
 variable{4} = 'Temperature_PTV_mask';
 variable{5} = 'Temperature_PTV_uncertainty';
 variable{6} = 'Absolute_Humidity_PTV'; 
 variable{7} = 'Absolute_Humidity_PTV_mask';
 variable{8} = 'Absolute_Humidity_PTV_uncertainty'; 
 variable{9} = 'Aerosol_Backscatter_Coefficient_PTV';
 variable{10} = 'Aerosol_Backscatter_Coefficient_PTV_mask';
 variable{11} = 'Aerosol_Backscatter_Coefficient_PTV_uncertainty';
 variable{12} = 'Backscatter_Photon_Counts_828'; 

  variable{16} = 'Absolute_Humidity_Standard'; 
  variable{17} = 'Absolute_Humidity_Standard_mask';
  variable{18} = 'Absolute_Humidity_Standard_uncertainty'; 

  variable{26} = 'Absolute_Humidity_MultiPulse'; 
  variable{27} = 'Absolute_Humidity_MultiPulse_mask';
  variable{28} = 'Absolute_Humidity_MultiPulse_uncertainty'; 

 variable{30} =  'Surface_Temperature';
 variable{31} =  'Surface_Pressure';
 variable{32} =  'Surface_Absolute_Humidity';
    
 variable{40} = 'Absolute_Humidity_ERA5';
 variable{41} = 'Temperature_ERA5';

for jj = 1:size(Pythonfilename,2)
  filename = Pythonfilename{jj};
  %date = filename(end-15:end-10);
  %n = datenum(date, 'yymmdd');
  date = filename(end-15:end-8);
  n = datenum(date, 'yyyymmdd');
  
  ncid = netcdf.open(filename, 'NC_NOWRITE');
    %ncdisp(filename, '/', 'min') % use this to display all variables

    time{jj} = ncread(filename,variable{1}); 
    alt{jj} = ncread(filename,variable{2});

    T{jj}  = ncread(filename,variable{3});  
    T_mask{jj} = ncread(filename,variable{4}); 
    T_var{jj} = ncread(filename,variable{5}); 
    T{jj}(T_mask{jj} == 1) = nan;
    T_var{jj}(T_mask{jj} == 1) = nan; 
    T_model{jj}  = ncread(filename,variable{41});  


    T_surf{jj} =  ncread(filename,variable{30});
    P_surf{jj} =  ncread(filename,variable{31});
    AH_surf{jj} =  ncread(filename,variable{32});
    
    AH{jj}  = ncread(filename,variable{6});  
    AH_mask{jj} = ncread(filename,variable{7}); 
    AH_var{jj} = ncread(filename,variable{8}); 
    AH{jj}(AH_mask{jj} == 1) = nan;
    AH_var{jj}(AH_mask{jj} == 1) = nan; 
    AH{jj}(AH_var{jj} > 25) = nan;

    AH_1{jj}  = ncread(filename,variable{16});  
    AH_1_mask{jj} = ncread(filename,variable{17}); 
    AH_1_var{jj} = ncread(filename,variable{18});
    AH_1{jj}(AH_1_mask{jj} == 1) = nan;
    AH_1_var{jj}(AH_1_mask{jj} == 1) = nan; 
    AH_1{jj}(AH_1_var{jj} > 3) = nan;

    AH_2{jj}  = ncread(filename,variable{26});  
    AH_2_mask{jj} = ncread(filename,variable{27}); 
    AH_2_var{jj} = ncread(filename,variable{28});
    AH_2{jj}(AH_2_mask{jj} == 1) = nan;
    AH_2_var{jj}(AH_2_mask{jj} == 1) = nan; 
    AH_2{jj}(AH_2_var{jj} > 3) = nan;

    AH_3{jj}  = ncread(filename,variable{40});  

    ABC{jj}  = ncread(filename,variable{9});   
    ABC_mask{jj} = ncread(filename,variable{10}); 
    ABC_var{jj} = ncread(filename,variable{11});
    ABC{jj}(ABC_mask{jj} == 1) = nan;
   % ABC_var{jj}(ABC_mask{jj} == 1) = nan; 
   % ABC{jj}(ABC_var{jj} > 1e-11) = nan;
   % ABC{jj}(ABC{jj} < 1e-12) = nan;
    Counts{jj} = ncread(filename,variable{12});   
    
  netcdf.close(ncid); 

  %convert from Unix time to date number (days since Jan 0 0000) 
  duration{jj} =  n+double(time{jj}/3600/24);
end

% ----------------------------------------------------------------------
% NEW PREALLOCATION SECTION
% ----------------------------------------------------------------------

% Determine the size of the final combined arrays
num_files = size(Pythonfilename, 2);
if num_files == 0
    % Handle case where no files were selected
    warning('No files selected. Exiting preallocation section.');
else
    % Total number of time steps (columns/rows)
    total_timesteps = sum(cellfun('size', duration, 1)); 

    % Number of range gates (rows for 2D data)
    num_ranges = size(AH{1}, 1); 

    % Preallocate 2D arrays (Range x Time)
    comb_AH = nan(num_ranges, total_timesteps);
    comb_AH_var = nan(num_ranges, total_timesteps);
    comb_AH_1 = nan(num_ranges, total_timesteps);
    comb_AH_1_var = nan(num_ranges, total_timesteps);
    comb_AH_2 = nan(num_ranges, total_timesteps);
    comb_AH_2_var = nan(num_ranges, total_timesteps);
    comb_AH_3 = nan(num_ranges, total_timesteps);
    comb_ABC = nan(num_ranges, total_timesteps);
    comb_ABC_var = nan(num_ranges, total_timesteps);
    comb_T = nan(num_ranges, total_timesteps);
    comb_T_var = nan(num_ranges, total_timesteps);
    comb_T_model = nan(num_ranges, total_timesteps);
    comb_Counts = nan(num_ranges, total_timesteps);

    % Preallocate 1D arrays (Time x 1)
    comb_duration = nan(total_timesteps, 1);
    comb_T_surf = nan(total_timesteps, 1); 
    comb_P_surf = nan(total_timesteps, 1); 
    comb_AH_surf = nan(total_timesteps, 1); 

    % Loop counter for filling the preallocated arrays
    start_col = 1;

% ----------------------------------------------------------------------
% REVISED COMBINATION LOOP (Lines 137-155 are replaced by this)
% ----------------------------------------------------------------------

    for jj = 1:num_files
        
        % Determine the size of the current file's time dimension
        current_timesteps = size(AH{jj}, 2); 
        end_col = start_col + current_timesteps - 1;

        % Fill the preallocated arrays by indexing
        comb_duration(start_col:end_col) = duration{jj};
        
        % 2D data (Range x Time) - columns correspond to time
        comb_AH(:, start_col:end_col) = AH{jj};
        comb_AH_var(:, start_col:end_col) = AH_var{jj};
        comb_AH_1(:, start_col:end_col) = AH_1{jj};
        comb_AH_1_var(:, start_col:end_col) = AH_1_var{jj};
        comb_AH_2(:, start_col:end_col) = AH_2{jj};
        comb_AH_2_var(:, start_col:end_col) = AH_2_var{jj};
        comb_AH_3(:, start_col:end_col) = AH_3{jj};
        comb_ABC(:, start_col:end_col) = ABC{jj};
        comb_ABC_var(:, start_col:end_col) = ABC_var{jj};
        comb_T(:, start_col:end_col) = T{jj};
        comb_T_var(:, start_col:end_col) = T_var{jj};
        comb_T_model(:, start_col:end_col) = T_model{jj};
        comb_Counts(:, start_col:end_col) = Counts{jj};
        
        % 1D surface data (Time x 1) - rows correspond to time
        comb_T_surf(start_col:end_col) = T_surf{jj};
        comb_P_surf(start_col:end_col) = P_surf{jj};
        comb_AH_surf(start_col:end_col) = AH_surf{jj};
        
        % Update the start column for the next iteration
        start_col = end_col + 1;
    end
end

%  for jj = 1:size(Pythonfilename,2)
%   if jj == 1
%       comb_duration = duration{jj};
%       comb_AH = AH{jj};
%       comb_AH_var = AH_var{jj};
%       comb_AH_1 = AH_1{jj};
%       comb_AH_1_var = AH_1_var{jj};
%       comb_AH_2 = AH_2{jj};
%       comb_AH_2_var = AH_2_var{jj};
%       comb_AH_3 = AH_3{jj};
%       comb_ABC = ABC{jj};
%       comb_ABC_var = ABC_var{jj};
%       comb_T = T{jj};
%       comb_T_var = T_var{jj};
%       comb_T_model = T_model{jj};
%       comb_T_surf = T_surf{jj}; 
%       comb_P_surf = P_surf{jj}; 
%       comb_AH_surf = AH_surf{jj}; 
%       comb_Counts = Counts{jj}; 
%   else
%       comb_duration = [comb_duration; duration{jj}];
%       % find the maximum range to accumulate (in case it changes)
%       % max_range = min(cellfun('size',alt,1))
%       comb_AH = [comb_AH AH{jj}];
%       comb_AH_var = [comb_AH_var AH_var{jj}];
%       comb_AH_1 = [comb_AH_1 AH_1{jj}];
%       comb_AH_1_var = [comb_AH_1_var AH_1_var{jj}];
%       comb_AH_2 = [comb_AH_2 AH_2{jj}];
%       comb_AH_2_var = [comb_AH_2_var AH_2_var{jj}];
%       comb_AH_3 = [comb_AH_3 AH_3{jj}];
%       comb_ABC = [comb_ABC ABC{jj}];
%       comb_ABC_var = [comb_ABC_var ABC_var{jj}];
%       comb_T = [comb_T T{jj}];
%       comb_T_var = [comb_T_var T_var{jj}];
%       comb_T_model = [comb_T_model T_model{jj}];
%       comb_T_surf = [comb_T_surf; T_surf{jj}];
%       comb_P_surf = [comb_P_surf; P_surf{jj}];
%       comb_AH_surf = [comb_AH_surf; AH_surf{jj}];
%       comb_Counts = [comb_Counts Counts{jj}];
%   end
% end

scrsz = get(0,'ScreenSize');
date=datestr(n, 'yyyy-mmm-dd');
plot_size1 = [scrsz(4)/1.5 scrsz(4)/10 scrsz(3)/1.5 scrsz(4)/3];
font_size = 16;

x = comb_duration;
y = (alt{1}./1000);
xData =  linspace( fix(min(x)),  ceil(max(x)), round((ceil(max(x))-fix(min(x)))/skip)+1 );
xData_m =  linspace( fix(min(x)),  ceil(max(x)), round((ceil(max(x))-fix(min(x)))/skip*24/2)+1 );


% % plot the T
 % Z = real(comb_T-273.15);
  Z = real(comb_T);
  figure1 = figure('Position',plot_size1);
  set(gcf,'renderer','zbuffer');
  h = pcolor(x, y, Z);
  set(h, 'EdgeColor', 'none'); 
  axis xy; colorbar('EastOutside'); 
%  caxis([-15 30]);
  clim([268 298]);
  axis([fix(min(x)) ceil(max(x)) 0 6]) 
%  shading interp
  set(gca, 'XTick',  xData)
  set(gca,'XMinorTick','on')
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  hh = title({[node, ' Temp, PTV (K)']},...
       'fontweight','b','fontsize',font_size);  
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
 % colormap(jet)
  colormap(plasma)
  set(gca,'Fontsize',font_size,'Fontweight','b');

% plot the atmospheric backscatter coefficient 
  Z = real(comb_ABC);
  figure2 = figure('Position',plot_size1);
  set(gcf,'renderer','zbuffer');
  h = pcolor(x, y, Z);
  set(h, 'EdgeColor', 'none'); 
  axis xy; 
  colorbar('EastOutside'); 
  %caxis([0 12]);
  %ylim([-.1 6]);
  axis([fix(min(x)) ceil(max(x)) 0 6])
  %  shading interp
  set(gca, 'XTick',  xData)
  set(gca,'XMinorTick','on')
  xAx = get(gca,'XAxis');
  xAx.MinorTickValues=xData_m;
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  set(gca,'Zscale', 'log')
  set(gca,'Colorscale', 'log')
  set(gca,'Zscale', 'linear')
  hh = title({[node, ' Aerosol Backscatter Coefficient, PTV (m^{-1} sr^{-1})']},...
       'fontweight','b','fontsize',font_size);     
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
 % datetick('x','HH:MM');
 % colormap(jet)
  colormap(viridis)
  set(gca,'Fontsize',font_size,'Fontweight','b');
  clim([5e-9 1e-6]);



% plot attenuated backscatter 
  Z = real(comb_Counts);
  figure3 = figure('Position',plot_size1);
  set(gcf,'renderer','zbuffer');
  h = pcolor(x, y, Z);
  set(h, 'EdgeColor', 'none'); 
  axis xy; 
  colorbar('EastOutside'); 
  %caxis([0 12]);
  %ylim([-.1 6]);
  axis([fix(min(x)) ceil(max(x)) 0 6])
  %  shading interp
  set(gca, 'XTick',  xData)
  set(gca,'XMinorTick','on')
  xAx = get(gca,'XAxis');
  xAx.MinorTickValues=xData_m;
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  set(gca,'Zscale', 'log')
  set(gca,'Colorscale', 'log')
  set(gca,'Zscale', 'linear')
  hh = title({[node, ' Attenuated Backscatter, 828 nm ']},...
       'fontweight','b','fontsize',font_size);     
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
 % datetick('x','HH:MM');
  colormap(jet)
 % colormap(viridis)
  set(gca,'Fontsize',font_size,'Fontweight','b');
  clim([1e2 1e6]);


% plot the AH_PTV 
  Z = real(comb_AH);
  figure4 = figure('Position',plot_size1);
  set(gcf,'renderer','zbuffer');
  h = pcolor(x, y, Z);
  set(h, 'EdgeColor', 'none'); 
  axis xy; 
  colorbar('EastOutside'); 
  clim([0 25]);
  axis([fix(min(x)) ceil(max(x)) 0 6]) 
%  shading interp
  set(gca, 'XTick',  xData)
  set(gca,'XMinorTick','on')
  xAx = get(gca,'XAxis');
  xAx.MinorTickValues=xData_m;
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  hh = title({[node, ' Absolute Humidity, PTV (g m^{-3})']},...
       'fontweight','b','fontsize',font_size);  
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
%  colormap(jet)
  colormap(CM_YlGnBu(64))
  set(gca,'Fontsize',font_size,'Fontweight','b');


% plot the AH_Standard 
  Z = real(comb_AH_1);
  figure5 = figure('Position',plot_size1);
  set(gcf,'renderer','zbuffer');
  h = pcolor(x, y, Z);
  set(h, 'EdgeColor', 'none'); 
  axis xy; 
  colorbar('EastOutside'); 
  clim([0 25]);
  axis([fix(min(x)) ceil(max(x)) 0 6]) 
%  shading interp
  set(gca, 'XTick',  xData)
  set(gca,'XMinorTick','on')
  xAx = get(gca,'XAxis');
  xAx.MinorTickValues=xData_m;
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  hh = title({[node, ' Absolute Humidity, Standard (g m^{-3})']},...
       'fontweight','b','fontsize',font_size);  
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
%  colormap(jet)
  colormap(CM_YlGnBu(64))
  set(gca,'Fontsize',font_size,'Fontweight','b');

% plot the AH_MultiPulse 
  Z = real(comb_AH_2);
  figure6 = figure('Position',plot_size1);
  set(gcf,'renderer','zbuffer');
  h = pcolor(x, y, Z);
  set(h, 'EdgeColor', 'none'); 
  axis xy; 
  colorbar('EastOutside'); 
  clim([0 25]);
  axis([fix(min(x)) ceil(max(x)) 0 6]) 
%  shading interp
  set(gca, 'XTick',  xData)
  set(gca,'XMinorTick','on')
  xAx = get(gca,'XAxis');
  xAx.MinorTickValues=xData_m;
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  hh = title({[node, ' Absolute Humidity, MultiPulse (g m^{-3})']},...
       'fontweight','b','fontsize',font_size);  
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
%  colormap(jet)
  colormap(CM_YlGnBu(64))
  set(gca,'Fontsize',font_size,'Fontweight','b');


% plot the AH_ERA5_Model 
  Z = real(comb_AH_3);
  figure7 = figure('Position',plot_size1);
  set(gcf,'renderer','zbuffer');
  h = pcolor(x, y, Z);
  set(h, 'EdgeColor', 'none'); 
  axis xy; 
  colorbar('EastOutside'); 
  clim([0 25]);
  axis([fix(min(x)) ceil(max(x)) 0 6]) 
%  shading interp
  set(gca, 'XTick',  xData)
  set(gca,'XMinorTick','on')
  xAx = get(gca,'XAxis');
  xAx.MinorTickValues=xData_m;
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  hh = title({[node, ' Absolute Humidity, ERA5 Model (g m^{-3})']},...
       'fontweight','b','fontsize',font_size);  
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
%  colormap(jet)
  colormap(CM_YlGnBu(64))
  set(gca,'Fontsize',font_size,'Fontweight','b');


% 1. Define the total number of colors (must be even, e.g., 64)
N = 64;
N_half = N / 2; % 32 colors in each half
% 2. Define the key colors in RGB [R G B]
blue = [0 0 1];
red = [1 0 0];
gray = [0.9 0.9 0.9]; % Your desired zero-point color
% Create the Blue half (Transition from Blue to Gray)
R1 = linspace(blue(1), gray(1), N_half)';
G1 = linspace(blue(2), gray(2), N_half)';
B1 = linspace(blue(3), gray(3), N_half)';
cmap_half1 = [R1, G1, B1];
% Create the Red half (Transition from Gray to Red)
R2 = linspace(gray(1), red(1), N_half)';
G2 = linspace(gray(2), red(2), N_half)';
B2 = linspace(gray(3), red(3), N_half)';
cmap_half2 = [R2, G2, B2];
% Combine the two halves
smooth_gray_cmap = [cmap_half1; cmap_half2];


% plot the AH PTV-ERA5 
  Z = real(comb_AH-comb_AH_3); 
  figure7 = figure('Position',plot_size1);
  set(gcf,'renderer','zbuffer');
  h = pcolor(x, y, Z);
  set(h, 'EdgeColor', 'none'); 
  axis xy; 
  colorbar('EastOutside'); 
  clim([-10 10]);
  axis([fix(min(x)) ceil(max(x)) 0 6]) 
%  shading interp
  set(gca, 'XTick',  xData)
  set(gca,'XMinorTick','on')
  xAx = get(gca,'XAxis');
  xAx.MinorTickValues=xData_m;
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  hh = title({[node, ' Absolute Humidity, PTV-ERA5 (g m^{-3})']},...
       'fontweight','b','fontsize',font_size);  
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
  %colormap(redblue)
  colormap(smooth_gray_cmap)
  set(gca,'Fontsize',font_size,'Fontweight','b');


% plot the AH Standard-ERA5 
  Z = real(comb_AH_1-comb_AH_3);
  figure7 = figure('Position',plot_size1);
  set(gcf,'renderer','zbuffer');
  h = pcolor(x, y, Z);
  set(h, 'EdgeColor', 'none'); 
  axis xy; 
  colorbar('EastOutside'); 
  clim([-10 10]);
  axis([fix(min(x)) ceil(max(x)) 0 6]) 
%  shading interp
  set(gca, 'XTick',  xData)
  set(gca,'XMinorTick','on')
  xAx = get(gca,'XAxis');
  xAx.MinorTickValues=xData_m;
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  hh = title({[node, ' Absolute Humidity, Standard-ERA5 (g m^{-3})']},...
       'fontweight','b','fontsize',font_size);  
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
  colormap(redblue)
  colormap(smooth_gray_cmap)
  set(gca,'Fontsize',font_size,'Fontweight','b');

% plot the AH MultiPulse-ERA5 
  Z = real(comb_AH_2-comb_AH_3);
  figure7 = figure('Position',plot_size1);
  set(gcf,'renderer','zbuffer');
  h = pcolor(x, y, Z);
  set(h, 'EdgeColor', 'none'); 
  axis xy; 
  colorbar('EastOutside'); 
  clim([-10 10]);
  axis([fix(min(x)) ceil(max(x)) 0 6]) 
%  shading interp
  set(gca, 'XTick',  xData)
  set(gca,'XMinorTick','on')
  xAx = get(gca,'XAxis');
  xAx.MinorTickValues=xData_m;
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  hh = title({[node, ' Absolute Humidity, MultiPulse-ERA5 (g m^{-3})']},...
       'fontweight','b','fontsize',font_size);  
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
%  colormap(redblue)
  colormap(smooth_gray_cmap)
  set(gca,'Fontsize',font_size,'Fontweight','b');

  % plot the T PTV-ERA5 
  Z = real(comb_T-comb_T_model);
  figure7 = figure('Position',plot_size1);
  set(gcf,'renderer','zbuffer');
  h = pcolor(x, y, Z);
  set(h, 'EdgeColor', 'none'); 
  axis xy; 
  colorbar('EastOutside'); 
  clim([-5 5]);
  axis([fix(min(x)) ceil(max(x)) 0 6]) 
%  shading interp
  set(gca, 'XTick',  xData)
  set(gca,'XMinorTick','on')
  xAx = get(gca,'XAxis');
  xAx.MinorTickValues=xData_m;
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  hh = title({[node, ' Temperature, PTV-ERA5 (K)']},...
       'fontweight','b','fontsize',font_size);  
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
%  colormap(redblue)
  colormap(smooth_gray_cmap)
  set(gca,'Fontsize',font_size,'Fontweight','b');


  % plot the AH uncertainty PTV 
  Z = real(comb_AH_var);
  figure7 = figure('Position',plot_size1);
  set(gcf,'renderer','zbuffer');
  h = pcolor(x, y, Z);
  set(h, 'EdgeColor', 'none'); 
  axis xy; 
  colorbar('EastOutside'); 
  clim([-5 5]);
  axis([fix(min(x)) ceil(max(x)) 0 6]) 
%  shading interp
  set(gca, 'XTick',  xData)
  set(gca,'XMinorTick','on')
  xAx = get(gca,'XAxis');
  xAx.MinorTickValues=xData_m;
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  hh = title({[node, ' Absolute Humidity Uncertainty, PTV (g m^{-3})']},...
       'fontweight','b','fontsize',font_size);  
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
  colormap(redblue)
  set(gca,'Fontsize',font_size,'Fontweight','b');

  % plot the AH uncertainty Standard 
  Z = real(comb_AH_1_var);
  figure7 = figure('Position',plot_size1);
  set(gcf,'renderer','zbuffer');
  h = pcolor(x, y, Z);
  set(h, 'EdgeColor', 'none'); 
  axis xy; 
  colorbar('EastOutside'); 
  clim([-5 5]);
  axis([fix(min(x)) ceil(max(x)) 0 6]) 
%  shading interp
  set(gca, 'XTick',  xData)
  set(gca,'XMinorTick','on')
  xAx = get(gca,'XAxis');
  xAx.MinorTickValues=xData_m;
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  hh = title({[node, ' Absolute Humidity Uncertainty, Standard (g m^{-3})']},...
       'fontweight','b','fontsize',font_size);  
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
  colormap(redblue)
  set(gca,'Fontsize',font_size,'Fontweight','b');

  % plot the AH uncertainty MultiPulse 
  Z = real(comb_AH_2_var);
  figure7 = figure('Position',plot_size1);
  set(gcf,'renderer','zbuffer');
  h = pcolor(x, y, Z);
  set(h, 'EdgeColor', 'none'); 
  axis xy; 
  colorbar('EastOutside'); 
  clim([-5 5]);
  axis([fix(min(x)) ceil(max(x)) 0 6]) 
%  shading interp
  set(gca, 'XTick',  xData)
  set(gca,'XMinorTick','on')
  xAx = get(gca,'XAxis');
  xAx.MinorTickValues=xData_m;
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  hh = title({[node, ' Absolute Humidity Uncertainty, MultiPulse (g m^{-3})']},...
       'fontweight','b','fontsize',font_size);  
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
  colormap(redblue)
  set(gca,'Fontsize',font_size,'Fontweight','b');

cd(strcat(plot_path,'/mpd/Plots'))
%plot_size = [1 1 1920 250];
plot_size = [100 100 1920 225];


   
  FigH = figure(1);
  % set(gca,'Fontsize',16,'Fontweight','b'); 
  % set(FigH, 'PaperUnits', 'points', 'PaperPosition',  plot_size);
  % name=strcat(date, node, ' T_PTV_multi'); 
  % print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
  FigH.Units = 'pixels'; % Ensure units are pixels for direct mapping to your old 'PaperPosition' width/height
  FigH.Position = plot_size; % x, y, width, height in pixels
  name=char(strcat(node, "_", date, ' T_PTV_multi')); 
  exportgraphics(FigH, [name, '.png'], 'Resolution', 150);

  FigH = figure(2);
  % set(gca,'Fontsize',16,'Fontweight','b'); 
  % set(FigH, 'PaperUnits', 'points', 'PaperPosition',  plot_size);
  % name=strcat(date, node, ' Back_Coeff_PTV_multi'); 
  % print(FigH, name, '-dpng', '-r0') % set at the screen resolution  drawnow;
  FigH.Units = 'pixels'; % Ensure units are pixels for direct mapping to your old 'PaperPosition' width/height
  FigH.Position = plot_size; % x, y, width, height in pixels
  name=char(strcat(node, "_", date, ' Back_Coeff_PTV_multi')); 
  exportgraphics(FigH, [name, '.png'], 'Resolution', 150);
  
  FigH = figure(3);
  % set(gca,'Fontsize',16,'Fontweight','b'); 
  % set(FigH, 'PaperUnits', 'points', 'PaperPosition',  plot_size);
  % name=strcat(date, node, ' 828_Counts_multi'); 
  % print(FigH, name, '-dpng', '-r0') % set at the screen resolution
  drawnow;
  FigH.Units = 'pixels'; % Ensure units are pixels for direct mapping to your old 'PaperPosition' width/height
  FigH.Position = plot_size; % x, y, width, height in pixels
  name=char(strcat(node, "_", date, ' 828_Counts_multi')); 
  exportgraphics(FigH, [name, '.png'], 'Resolution', 150);

  FigH = figure(4);
  % set(gca,'Fontsize',16,'Fontweight','b'); 
  % set(FigH, 'PaperUnits', 'points', 'PaperPosition',  plot_size);
  % name=strcat(date, node, ' WV_PTV_multi'); 
  % print(FigH, name, '-dpng', '-r0') % set at the screen resolution
  drawnow;
  FigH.Units = 'pixels'; % Ensure units are pixels for direct mapping to your old 'PaperPosition' width/height
  FigH.Position = plot_size; % x, y, width, height in pixels
  name=char(strcat(node, "_", date, ' WV_PTV_multi')); 
  exportgraphics(FigH, [name, '.png'], 'Resolution', 150);

  FigH = figure(5);
  % set(gca,'Fontsize',16,'Fontweight','b'); 
  % set(FigH, 'PaperUnits', 'points', 'PaperPosition',  plot_size);
  % name=strcat(date, node, ' WV_Standard_multi'); 
  % print(FigH, name, '-dpng', '-r0') % set at the screen resolution
  drawnow;
  FigH.Units = 'pixels'; % Ensure units are pixels for direct mapping to your old 'PaperPosition' width/height
  FigH.Position = plot_size; % x, y, width, height in pixels
  name=char(strcat(node, "_", date, ' WV_Standard_multi')); 
  exportgraphics(FigH, [name, '.png'], 'Resolution', 150);

  FigH = figure(6);
  % set(gca,'Fontsize',16,'Fontweight','b'); 
  % set(FigH, 'PaperUnits', 'points', 'PaperPosition',  plot_size);
  % name=strcat(date, node, ' WV_MultiPulse_multi'); 
  % print(FigH, name, '-dpng', '-r0') % set at the screen resolution
  drawnow;
  FigH.Units = 'pixels'; % Ensure units are pixels for direct mapping to your old 'PaperPosition' width/height
  FigH.Position = plot_size; % x, y, width, height in pixels
  name=char(strcat(node, "_", date, ' WV_MultiPulse_multi')); 
  exportgraphics(FigH, [name, '.png'], 'Resolution', 150);
  
  FigH = figure(7);
 % set(gca,'Fontsize',16,'Fontweight','b'); 
 % set(FigH, 'PaperUnits', 'points', 'PaperPosition',  plot_size);
 % name=strcat(date, node, ' WV_Model_multi'); 
 % print(FigH, name, '-dpng', '-r0') % set at the screen resolution
  drawnow;
  FigH.Units = 'pixels'; % Ensure units are pixels for direct mapping to your old 'PaperPosition' width/height
  FigH.Position = plot_size; % x, y, width, height in pixels
  name=char(strcat(node, "_", date, ' WV_Model_multi')); 
  exportgraphics(FigH, [name, '.png'], 'Resolution', 150);

  FigH = figure(8);
  %set(gca,'Fontsize',16,'Fontweight','b'); 
  %set(FigH, 'PaperUnits', 'points', 'PaperPosition',  plot_size);
  %name=strcat(date, node, ' WV_Diff_PTV_Model'); 
  %print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
  drawnow;
  FigH.Units = 'pixels'; % Ensure units are pixels for direct mapping to your old 'PaperPosition' width/height
  FigH.Position = plot_size; % x, y, width, height in pixels
  name=char(strcat(node, "_", date, ' WV_Diff_PTV_Model')); 
  exportgraphics(FigH, [name, '.png'], 'Resolution', 150);

  FigH = figure(9);
  %set(gca,'Fontsize',16,'Fontweight','b'); 
  %set(FigH, 'PaperUnits', 'points', 'PaperPosition',  plot_size);
  %name=strcat(date, node, ' WV_Diff_Standard_Model'); 
  %print(FigH, name, '-dpng', '-r0') % set at the screen resolution
  drawnow;
  FigH.Units = 'pixels'; % Ensure units are pixels for direct mapping to your old 'PaperPosition' width/height
  FigH.Position = plot_size; % x, y, width, height in pixels
  name=char(strcat(node, "_", date, ' WV_Diff_Standard_Model')); 
  exportgraphics(FigH, [name, '.png'], 'Resolution', 150);
  
  FigH = figure(10);
  %set(gca,'Fontsize',16,'Fontweight','b'); 
  %set(FigH, 'PaperUnits', 'points', 'PaperPosition',  plot_size);
  %name=strcat(date, node, ' WV_Diff_MultiPulse_Model'); 
  %print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
  drawnow;
  FigH.Units = 'pixels'; % Ensure units are pixels for direct mapping to your old 'PaperPosition' width/height
  FigH.Position = plot_size; % x, y, width, height in pixels
  name=char(strcat(node, "_", date, '_WV_Diff_MultiPulse_Model')); 
  exportgraphics(FigH, [name, '.png'], 'Resolution', 150);

  FigH = figure(11);
  %set(gca,'Fontsize',16,'Fontweight','b'); 
  %set(FigH, 'PaperUnits', 'points', 'PaperPosition',  plot_size);
  %name=strcat(date, node, ' WV_Diff_MultiPulse_Model'); 
  %print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
  drawnow;
  FigH.Units = 'pixels'; % Ensure units are pixels for direct mapping to your old 'PaperPosition' width/height
  FigH.Position = plot_size; % x, y, width, height in pixels
  name=char(strcat(node, "_", date, '_T_Diff_PTV_Model')); 
  exportgraphics(FigH, [name, '.png'], 'Resolution', 150);
  
  
  FigH = figure(12);
  % set(gca,'Fontsize',16,'Fontweight','b'); 
  % set(FigH, 'PaperUnits', 'points', 'PaperPosition',  plot_size);
  % name=strcat(date, node, ' WV_PTV_Uncertainty'); 
  % print(FigH, name, '-dpng', '-r0') % set at the screen resolution
  drawnow;
  FigH.Units = 'pixels'; % Ensure units are pixels for direct mapping to your old 'PaperPosition' width/height
  FigH.Position = plot_size; % x, y, width, height in pixels
  name=char(strcat(node, "_", date, ' WV_PTV_Uncertainty')); 
  exportgraphics(FigH, [name, '.png'], 'Resolution', 150);

  FigH = figure(13);
  % set(gca,'Fontsize',16,'Fontweight','b'); 
  % set(FigH, 'PaperUnits', 'points', 'PaperPosition',  plot_size);
  % name=strcat(date, node, ' WV_Standard_Uncertainty'); 
  % print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
  drawnow;
  FigH.Units = 'pixels'; % Ensure units are pixels for direct mapping to your old 'PaperPosition' width/height
  FigH.Position = plot_size; % x, y, width, height in pixels
  name=char(strcat(node, "_", date, ' WV_Standard_Uncertainty')); 
  exportgraphics(FigH, [name, '.png'], 'Resolution', 150);
  
  FigH = figure(14);
  % set(gca,'Fontsize',16,'Fontweight','b'); 
  % set(FigH, 'PaperUnits', 'points', 'PaperPosition',  plot_size);
  % name=strcat(date, node, ' WV_MultiPulse_Uncertainty'); 
  % print(FigH, name, '-dpng', '-r0') % set at the screen resolution
  drawnow;
  FigH.Units = 'pixels'; % Ensure units are pixels for direct mapping to your old 'PaperPosition' width/height
  FigH.Position = plot_size; % x, y, width, height in pixels
  name=char(strcat(node, "_", date, ' WV_MultiPulse_Uncertainty')); 
  exportgraphics(FigH, [name, '.png'], 'Resolution', 150);

 % 
 % if flag.save_data == 1
 % 
 %  cd(d_save_data);     
 %  range = alt{1}'; 
 %  N_avg_comb = (comb_AH./1e6.*6.022E23./18.015);
 %  % duration = duration;
 %  % test1 = datestr(comb_duration(1:2), 'yyyy-mmm-dd HH:MM');
 %  % test2 = Thermo_BLH2(1:2);
 %  test3 = [comb_duration Thermo_BLH1 Thermo_BLH2 Thermo_BLH3];  
 %  figure(1)
 %  plot(test3(:,1), test3(:,2))
 %  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
 % 
 % 
 % 
 % 
 %  %cd('/Users/spuler/Desktop/WV_DIAL_data') % point to the directory where data is stored 
 %  name=strcat(date, '_combined');
 %  if strcmp(node,'MPD01')==1
 %      MPD01.N_avg_comb = N_avg_comb;
 %     % MPD01.RB_comb = RB_comb;
 %      MPD01.range = range;
 %      MPD01.time = duration;
 %      save(name, 'MPD01')
 %  elseif strcmp(node,'MPD02')==1
 %      MPD02.N_avg_comb = N_avg_comb;
 %    %  MPD02.RB_comb = RB_comb;
 %      MPD02.range = range;
 %      MPD02.time = duration;
 %      save(name, 'MPD02')
 %  elseif strcmp(node,'MPD03')==1
 %      MPD03.N_avg_comb = N_avg_comb;
 %    %  MPD03.RB_comb = RB_comb;
 %      MPD03.range = range;
 %      MPD03.time = duration;
 %      save(name, 'MPD03')
 %  elseif strcmp(node,'MPD04')==1
 %      MPD04.N_avg_comb = N_avg_comb;
 %    %  MPD04.RB_comb = RB_comb;
 %      MPD04.range = range;
 %      MPD04.time = duration;
 %      save(name, 'MPD04')
 %   elseif strcmp(node,'MPD05')==1
 %      MPD05.N_avg_comb = N_avg_comb;
 %   %   MPD05.RB_comb = RB_comb;
 %      MPD05.range = range;
 %      MPD05.time =  duration;
 %      save(name, 'MPD05')
 %  end
 % 
 % end
 % 