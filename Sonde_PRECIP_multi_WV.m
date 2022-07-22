elevation= 25; %MPD04 was at 25m elevation at Hsinchu weather station
flag.plot_overlay = 0; %plot sondes on the time vs hieght AH plot
flag.data_type = 0;  % 0=matlab WV, 1=python WV, 2=raman WV, 3=PTV WV
sonde_end_int = 15; % integration time (in min) +/- around the sound launch  
WV_min = 0;
WV_max = 30;

%d=pwd;

%cd('/Volumes/documents/MPD/Sondes_CSU')
%plot_path = '/Volumes/documents/MPD/Plots/';
%cd('/Volumes/eol/fog1/rsfdata/MPD/mpd_ancillary_data/radiosondes/CSU_PrePRECIP')
cd('/Volumes/Macintosh HD/Users/spuler/Desktop/mpd/sondes');
plot_path = '/Volumes/Macintosh HD/Users/spuler/Desktop/mpd/Plots';
%cd('/Volumes/eol/sci/voemel/data/radiosondes/boulder/ncdf')
[sondefilename, sondedir] = uigetfile('*.*','Select the sonde file', 'MultiSelect', 'on');
%flag.MR = 0; % instead of absolute humidity plot the mixing ratio
jj=1;
%cd= d;

for jj = 1:size(sondefilename,2)
    cd('/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_processing/')
    %cd('/Users/lroot/Documents/GitHub/EOL_Lidar_Matlab_processing/')
   if flag.data_type == 2    % add the following three lines for Raman 
      range_grid_size = 60;  
      N_avg_comb = (comb_Raman_AH./1e6.*6.022E23./18.015);
      duration = comb_Raman_duration;
   elseif flag.data_type == 1  % add the following lines for Python processing of the MPD 
      %range_grid_size = 37.5; 
      range_grid_size = diff(alt{1}');
      range_grid_size = range_grid_size(1,1)
      N_avg_comb = (comb_AH./1e6.*6.022E23./18.015);
      duration = comb_duration;
      range_grid_in = alt{1}';
   elseif flag.data_type == 0
       range_grid_size = diff(range(1:2)) 
       range_grid_in = range;
       comb_AH_var = N_error_comb.*1e6./6.022E23.*18.015;
    elseif flag.data_type == 3
      %range_grid_size = 37.5; 
      range_grid_size = diff(alt{1}');
      range_grid_size = range_grid_size(1,1)
      N_avg_comb = (comb_AH./1e6.*6.022E23./18.015);
      comb_AH_var = N_avg_comb*0;
      duration = comb_duration;
      range_grid_in = alt{1}';
   end
   [xx(jj,:), yy(jj,:), range_grid] = Sonde_read_PRECIP_files(jj, elevation, sondedir, sondefilename,  N_avg_comb, duration, range_grid_size, range_grid_in, comb_AH_var, sonde_end_int, plot_path, flag); 
  % [xx(jj,:), yy(jj,:), range_grid] = Sonde_read_CSU_files_T_lapse(jj, elevation, sondedir, sondefilename,  N_avg_comb, duration, range_grid_size, range_grid_in, comb_AH_var, sonde_end_int, T_lapse, plot_path, flag); 
   %Sonde_DIAL_comparison_funct_v6(N_H2O, sonde_top, sonde_range, t, date, T_sonde, P_sonde, sonde_stop, shift, error_threshold, Wind_speed, save_figs, ID_sonde);
   %Sonde_DIAL_comparison_funct_Python(N_H2O, sonde_top, sonde_range, t, date, T_sonde, P_sonde, sonde_stop, shift, error_threshold, Wind_speed, save_figs)
 end
    
    
% create a line plot of the sonde vs MPD data
%scrsz = get(0,'ScreenSize');
scrsz = [1  1  1920 1200];
Scrsize=[scrsz(4)/1 scrsz(4)/1 scrsz(3)/1.5 scrsz(4)/1.5];
font_size = 14;
bins = WV_max*4; % bin size is x 0.1 x 0.1 g/m^2
%bin_min = 1;
%bin_max = 400;

% remove the surface station data 
x_no_surface = xx(:,2:end); % Sondes (assume as truth) 
y_no_surface = yy(:,2:end); % MPD  

% remove the data below 225m (at 75 m gates, 1=0, 2=75, 3=150)  
x_no_surface = xx(:,4:end); % Sondes (assume as truth) 
y_no_surface = yy(:,4:end); % MPD  

%calculate the RMS using the sonde data as truth
RMSE_numerator = nansum((yy - xx).^2);
RMSE_denominator = sum(~isnan((yy - xx).^2));
RMSE = sqrt(RMSE_numerator./RMSE_denominator);
RMSPE = nanmean((yy - xx)./xx).^2;
% calcuate the Mean Bias Error
mean_error_numerator = nansum(yy - xx);
mean_error_denominator = sum(~isnan((yy - xx)));
mean_error = (mean_error_numerator./mean_error_denominator);
% calcuate the Percent Error
percent_error_numerator = nansum(abs(real((xx - yy)./yy)));
percent_error_denominator = sum(~isnan(real((xx - yy)./yy)));
percent_error = (percent_error_numerator./percent_error_denominator)*100;

figure(20)
subplot(1,4,1)
plot(percent_error, range_grid)
  set(gca,'Fontsize',20,'Fontweight','b'); %   
  grid on
  set(gca,'XMinorGrid','on','YMinorGrid','on')
  xlim([0 25]);
  ylim([0,6]);
  ylabel('Altitude, AGL (km)'); 
  xlabel('Mean Percent Error'); 
subplot(1,4,2)
plot(mean_error, range_grid)
  set(gca,'Fontsize',20,'Fontweight','b'); %   
  grid on
  set(gca,'XMinorGrid','on','YMinorGrid','on')
  xlim([-4 4]);
  ylim([0,6]);
  ylabel('Altitude, AGL (km)'); 
  xlabel('Mean Error (g m^{-3})'); 
subplot(1,4,3)
plot(RMSE, range_grid)
  set(gca,'Fontsize',20,'Fontweight','b'); %   
  grid on
  set(gca,'XMinorGrid','on','YMinorGrid','on')
  xlim([0 4]);
  ylim([0,6]);
  ylabel('Altitude, AGL (km)'); 
  xlabel('RMSE (g m^{-3})');
subplot(1,4,4)
plot(RMSE_denominator, range_grid)
  set(gca,'Fontsize',20,'Fontweight','b'); % 
  grid on
  set(gca,'XMinorGrid','on','YMinorGrid','on')
 % xlim([0 45]);
  ylim([0,6]);
  ylabel('Altitude, AGL (km)'); 
  xlabel('Samples used');



figure(21)
subplot(1,3,1)
plot(mean_error, range_grid)
  set(gca,'Fontsize',20,'Fontweight','b'); %   
  grid on
  set(gca,'XMinorGrid','on','YMinorGrid','on')
  xlim([-4 4]);
  ylim([0,6]);
  ylabel('Altitude, AGL (km)'); 
  xlabel('Mean Error (g m^{-3})'); 
subplot(1,3,2)
plot(RMSE, range_grid)
  set(gca,'Fontsize',20,'Fontweight','b'); %   
  grid on
  set(gca,'XMinorGrid','on','YMinorGrid','on')
  xlim([0 4]);
  ylim([0,6]);
  ylabel('Altitude, AGL (km)'); 
  xlabel('RMSE (g m^{-3})');
subplot(1,3,3)
plot(RMSE_denominator, range_grid)
  set(gca,'Fontsize',20,'Fontweight','b'); % 
  grid on
  set(gca,'XMinorGrid','on','YMinorGrid','on')
 % xlim([0 45]);
  ylim([0,6]);
  ylabel('Altitude, AGL (km)'); 
  xlabel('Samples used');
  
  
% reshape into one long array
x_sonde = reshape(x_no_surface,1,[]); %sondes
y_MPD05 = reshape(y_no_surface,1,[]); %MPD
xx0 = WV_min:1:WV_max;
yy0 = WV_min:1:WV_max;



y_MPD05(y_MPD05<-WV_max*3)=NaN; % remove points 3x outside of plot range 
y_MPD05(y_MPD05>WV_max*3)=NaN; % remove points 3x outside of plot range 
x_sonde(x_sonde<-WV_max*3)=NaN; % remove points 3x outside of plot range 
x_sonde(x_sonde>WV_max*3)=NaN; % remove points 3xoutside of plot range 

y_MPD05=real(y_MPD05); % add to catch non-real values

% calculate the best fit (in a least-squares sense) and the correlation coefficient
idx = (isnan(y_MPD05)|isnan(x_sonde)); %remove the NaNs
%num_samples = size(x_sonde(~idx), 2); 
%fit = polyfit(x_sonde(~idx),y_MPD05(~idx),1);
%[Corr, P_test] = corrcoef(x_sonde(~idx),y_MPD05(~idx))
%Cov = cov(x_sonde(~idx),y_MPD05(~idx))
%StDev = sqrt((cov(x_sonde(~idx),y_MPD05(~idx))))
mdl = fitlm(x_sonde(~idx),y_MPD05(~idx));
Rsquared = mdl.Rsquared.Ordinary
Corr = sqrt(Rsquared) 
StDev = mdl.RMSE
num_samples = mdl.NumObservations
fit = flip(mdl.Coefficients.Estimate')


% Evaluate fit equation using polyval
xrange = WV_min:0.1:WV_max;
y_est = polyval(fit,xrange);

figure(7); %Frequency histogram 
ax = subplot(1,1,1);
%figure('Position',Scrsize);
% this is temporary (the python code should remove these 
binscatter(x_sonde,y_MPD05, [bins bins]);    
  %N=hh.NumBins;
  %hh.NumBins =[bins bins];
  %hh.XLimits = [WV_min WV_max];
  %hh.YLimits = [WV_min WV_max];
ax.CLim = [0 50];
ax.XGrid = 'on';
axh.YGrid = 'on';
colormap('parula')
oldcmap = colormap;
colormap( flipud(oldcmap) );
c=colorbar;
c.Label.String = 'Sample pairs per bin';
%title(['WV x=MPD05, ', ' y=Raman Lidar'])
xlim([WV_min WV_max]);
%pause(0.001) % there is something odd about this all running at once
%set(gca,'Fontsize',font_size,'Fontweight','b');
ylim([WV_min WV_max]); 
ylabel('MPD05 absolute humidity (g m^{-3})'); 
xlabel('Radiosonde absolute humidity (g m^{-3})'); 
set(gca,'Fontsize',font_size,'Fontweight','b');
% Display fit infor on graph
text(0.5, WV_max-1, ['Number of samples = ' num2str(num_samples(1))], 'FontSize', 22, 'Color', 'r')
text(0.5, WV_max-2, ['y = ' num2str(fit(1),3) '*x + ' num2str(fit(2),3)], 'FontSize', 22, 'Color', 'r')
text(0.5, WV_max-3, ['Corr = ' num2str(Corr,3)], 'FontSize', 22, 'Color', 'r')
text(0.5, WV_max-4, ['StDev = ' num2str(StDev,3) 'g m^{-3}'], 'FontSize', 22, 'Color', 'r')
%text(1, 15, ['R^2 = ' num2str(Rsquared,3)], 'FontSize', 22, 'Color', 'r')

hold on
plot(xx0,yy0, 'k-')  % plot the 1:1 line
plot(xrange,y_est,'r--','LineWidth',2)  % plot the least squared fit line
hold off

figure(8)
plot(x_sonde, y_MPD05, 'o')
xlim([0 20])
ylim([0 20])
xlim([WV_min WV_max]);
%pause(0.001) % there is something odd about this all running at once
%set(gca,'Fontsize',font_size,'Fontweight','b');
ylim([WV_min WV_max]); 
ylabel('MPD05 absolute humidity (g m^{-3})'); 
xlabel('Radiosonde absolute humidity (g m^{-3})'); 
set(gca,'Fontsize',font_size,'Fontweight','b');
% Display fit infor on graph
text(0.5, WV_max-1, ['Number of samples = ' num2str(num_samples(1))], 'FontSize', 22, 'Color', 'r')
text(0.5, WV_max-2, ['y = ' num2str(fit(1),3) '*x + ' num2str(fit(2),3)], 'FontSize', 22, 'Color', 'r')
text(0.5, WV_max-3, ['Corr = ' num2str(Corr,3)], 'FontSize', 22, 'Color', 'r')
text(0.5, WV_max-4, ['StDev = ' num2str(StDev,3) 'g m^{-3}'], 'FontSize', 22, 'Color', 'r')
% Add trend line to plot
hold on
plot(xx0,yy0, 'k-')  % plot the 1:1 line
plot(xrange,y_est,'r--','LineWidth',2) % plot the least squared fit line
hold off

%cd('/Volumes/documents/WV_DIAL_data/plots/') % point to the directory where data is stored 
%cd('/Users/lroot/Desktop/mpd/Plots/') % point to the directory where data is stor
cd(plot_path)

FigH = figure(1);
set(gca,'Fontsize',16,'Fontweight','b'); % 
set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920 250]);
name=strcat(date, 'AH_multi'); 
print(FigH, name, '-dpng', '-r0') % set at the screen resolution 

FigH = figure(7);
set(gca,'Fontsize',24,'Fontweight','b'); % 
set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 800 800]);
name=strcat(date, 'Sonde_hist_multi'); 
print(FigH, name, '-dpng', '-r0') % set at the screen resolution 

% FigH = figure(8);
% set(gca,'Fontsize',24,'Fontweight','b'); % 
% set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 800 800]);
% name=strcat(date, 'Sonde_multi'); 
% print(FigH, name, '-dpng', '-r0') % set at the screen resolution 

FigH = figure(20);
set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1200 800]);
name=strcat(date, 'P_error'); 
print(FigH, name, '-dpng', '-r0') % set at the screen resolution

FigH = figure(21);
set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 800 800]);
name=strcat(date, 'error_RMSE'); 
print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
