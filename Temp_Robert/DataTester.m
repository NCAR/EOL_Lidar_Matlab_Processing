
clear; close all; clc

%% Default time and range
Time    = 150:300:(60*60*24-150);  % 5 minute resolution
Range   = 0:75:6e3;

%% Defining the data source
Systems     = {'05'};
Dates2Load  = datenum({'20210620';'20210623'},'yyyymmdd');
MatlabPath  = pwd;
serv_path = '/Volumes/fog1/rsfdata/MPD';
%MatlabPath =
%(strcat(serv_path,'/mpd_05_processed_data/Quickload/FullProcessingBootStrap'));
MatlabPath = (strcat(serv_path,'/mpd_05_processed_data/Quickload/ReProcessing'));
MatlabFile  = 'mpd%s.%s.Matlab.mat';    % File format


%% Loading water vapor data
Var2Load = {strcat('Retrievals.WaterVapor.', {'TimeStamp';'Range';'Smoothed2';});
            strcat('Retrievals.WaterVapor.', {'TimeStamp';'Range';'Mask'})
            strcat('Retrievals.Python.MPD.Humidity.', {'TimeStamp';'Range';'Value'});
            strcat('Retrievals.Python.MPD.Humidity.', {'TimeStamp';'Range';'Mask'})};
               
fprintf('-----------Loading WV Data-----------\n')   
WaterVapor = LoaderMPD(Systems,Dates2Load,Var2Load,MatlabPath,MatlabFile,Time,Range);
WaterVapor = WaterVapor{1};

% Applying the WV masks
WaterVapor{1}(WaterVapor{2}>0) = nan;
WaterVapor{3}(WaterVapor{4}>0) = nan;


%% Loading HSRL data
Var2Load = {strcat('Retrievals.HSRL.', {'TimeStamp';'Range';'Smoothed'});
            strcat('Retrievals.HSRL.', {'TimeStamp';'Range';'ABC'});
            strcat('Retrievals.HSRL.', {'TimeStamp';'Range';'Mask'});
            strcat('Retrievals.Python.MPD.BackRatio.', {'TimeStamp';'Range';'Value'});
            strcat('Retrievals.Python.MPD.BSCoefficient.', {'TimeStamp';'Range';'Value'});
            strcat('Retrievals.Python.MPD.BackRatio.', {'TimeStamp';'Range';'Mask'})};
 
fprintf('-----------Loading HSRL Data-----------\n')   
HSRL = LoaderMPD(Systems,Dates2Load,Var2Load,MatlabPath,MatlabFile,Time,Range);
HSRL = HSRL{1};

HSRL{1}(HSRL{3}>0) = nan;
HSRL{2}(HSRL{3}>0) = nan;
HSRL{4}(HSRL{6}>0) = nan;
HSRL{5}(HSRL{6}>0) = nan;

%% Loading temperature data
Var2Load = {strcat('Retrievals.Temperature.', {'TimeStamp';'Range';'Smoothed'});
            strcat('Retrievals.Temperature.', {'TimeStamp';'Range';'VarianceSm'});
            strcat('Retrievals.Python.NCIP.Temperature.', {'TimeStamp';'Range';'Value'})};

fprintf('-----------Loading Temp Data-----------\n')   
Temperature = LoaderMPD(Systems,Dates2Load,Var2Load,MatlabPath,MatlabFile,Time,Range);
Temperature = Temperature{1};

Temperature{1}(Temperature{2}>5^2) = nan;

%% Plotting
TimeArray = ((0.5:1:size(WaterVapor{1,1},2)).*(Time(2)-Time(1)))./60./60./24;

figure(1)
ToPlot = [1,3];
Titles = {'Matlab Absoluate Humidity','Python Absolute Humidity'};
for m=1:1:length(ToPlot)
    subplot(length(ToPlot),1,m)
    pcolor(TimeArray,Range./1e3,WaterVapor{ToPlot(m),1}); shading flat; colorbar; caxis([0,10])
  %  colormap(gca,CM_viridis(64))
    xlabel('Days after start'); ylabel('Range [km]'); title(Titles{m});
end

figure(2)
ToPlot = [1,4];
Titles = {'Matlab Backscatter Ratio','Python Backscatter Ratio'};
for m=1:1:length(ToPlot)
    subplot(length(ToPlot),1,m)
    pcolor(TimeArray,Range,HSRL{ToPlot(m),1}); shading flat; colorbar; caxis([0,10])
%    colormap(gca,flipud(CM_magma(64)))
    xlabel('Days after start'); ylabel('Range [km]'); title(Titles{m});
end

figure(3)
ToPlot = [1,3];
Titles = {'MicroPulse DIAL Temperature','NCEP Temperature'};
for m=1:1:length(ToPlot)
    subplot(length(ToPlot),1,m)
    pcolor(TimeArray,Range,Temperature{ToPlot(m),1}-273.15); shading flat; colorbar; caxis([-5,30])
%    colormap(gca,CM_cividis(64))
    xlabel('Days after start'); ylabel('Range [km]'); title(Titles{m});
end


% added by Spuler to integrate with other legacy programs

figure(4)
pcolor(TimeArray,Range,Temperature{ToPlot(1),1}); shading flat; colorbar; caxis([240,320])
colormap(jet)

font_size = 16;
scrsz = [1  1  1920 1200];
Scrnsize=[scrsz(4)/2 scrsz(4)/10 scrsz(3)/1 scrsz(4)/2];
x = Dates2Load(1)+TimeArray; 
skip=1;
xData =  linspace( fix(min(x)),  ceil(max(x)), round((ceil(max(x))-fix(min(x)))/skip)+1 );

figure('Position',Scrnsize)
set(gcf,'renderer','zbuffer');
h = pcolor(x,Range/1000,Temperature{ToPlot(1),1}-273.15); colorbar; caxis([-5,30])
%colormap(jet)
%colormap(gca,CM_viridis(64))
set(h, 'EdgeColor', 'none');
colorbar('EastOutside');
ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
set(gca, 'XTick',  xData)
set(gca,'TickDir','out');
set(gca,'TickLength',[0.005; 0.0025]);
datetick('x','dd-mmm-yy')
hh = title({['MPD',Systems{1}, ' Temperature (C)']},'fontweight','b','fontsize',font_size);
set(gca,'Fontsize',font_size,'Fontweight','b');

plot_path = '/Users/spuler/Desktop/mpd/Plots/';
cd(plot_path) % point to the directory where data is stored 
FigH = figure(5);
set(gca,'Fontsize',16,'Fontweight','b');  
set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920 250]);       %1500 x 300
name=strcat(date, Systems{1}, ' Temp_Matlab_multi_avg'); 
print(FigH, name, '-dpng', '-r0') % set at the screen resolution 


FigH = figure(3);
%set(gca,'Fontsize',16,'Fontweight','b');  
set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 500 400]);       %1500 x 300
name=strcat(date, Systems{1}, ' Temp_Matlab_NCEP_avg'); 
print(FigH, name, '-dpng', '-r0') % set at the screen resolution 

Temp_comb_avg = Temperature{ToPlot(1),1}';
Temp_comb = Temperature{ToPlot(1),1}';
Temp_comb_var = Temperature{2}';
range = Range;
duration = x';
lapse_rate = 0.0065; %standard atmosphere lapse rate
