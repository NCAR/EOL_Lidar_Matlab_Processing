
x_range_grid_size = 60;  
x_values = comb_Raman_AH;
x_duration = comb_Raman_duration;
y_values = N_avg_comb.*1e6./6.022E23.*18.015;
y_duration = duration;
y_range_grid_size = range_grid_size; 

% grid MPD tp Raman data in time 
y_values_grid = interp1(y_duration, y_values, x_duration, 'nearest');
% grid Raman lidar to MPD in range 
x_values_grid = interp1(Raman_alt{1}, x_values', 0:.075:6, 'nearest');
y_values_grid = y_values_grid(:,1:size(x_values_grid,1))';

x = x_duration;
y = 0:.075:6;
Z = real(x_values_grid);
% plot the Raman AH
  figure1 = figure('Position',plot_size1);
  set(gcf,'renderer','zbuffer');
  h = pcolor(x, y, Z);
  set(h, 'EdgeColor', 'none'); 
  axis xy; colorbar('EastOutside'); 
  caxis([0 18]);
  ylim([0 6])
Z = real(y_values_grid);
  % plot MPD AH
  figure1 = figure('Position',plot_size1);
  set(gcf,'renderer','zbuffer');
  h = pcolor(x, y, Z);
  set(h, 'EdgeColor', 'none'); 
  axis xy; colorbar('EastOutside'); 
  caxis([0 18]);
  ylim([0 6])
  

% create a line plot of the sonde vs MPD data
scrsz = get(0,'ScreenSize');
Scrsize=[scrsz(4)/1 scrsz(4)/1 scrsz(3)/1.5 scrsz(4)/1.5];
font_size = 14;
WV_min = 0;
WV_max = 20;
bins = WV_max*4; % bin size is x 0.1 x 0.1 g/m^2
%bin_min = 1;
%bin_max = 400;

% reshape into one long array
x_Raman = reshape(x_values_grid,1,[]); %Raman Lidar
y_MPD05 = reshape(y_values_grid,1,[]); %MPD
xx0 = WV_min:1:WV_max;
yy0 = WV_min:1:WV_max;


y_MPD05(y_MPD05<-WV_max*2)=NaN; % remove points 3x outside of plot range 
y_MPD05(y_MPD05>WV_max*2)=NaN; % remove points 3x outside of plot range 
x_Raman(x_Raman<-WV_max*2)=NaN; % remove points 3x outside of plot range 
x_Raman(x_Raman>WV_max*2)=NaN; % remove points 3xoutside of plot range 

% calculate the best fit (in a least-squares sense) and the correlation coefficient
idx = (isnan(y_MPD05)|isnan(x_Raman)); %remove the NaNs
%num_samples = size(x_sonde(~idx), 2); 
%fit = polyfit(x_sonde(~idx),y_MPD05(~idx),1);
%[Corr, P_test] = corrcoef(x_sonde(~idx),y_MPD05(~idx))
%Cov = cov(x_sonde(~idx),y_MPD05(~idx))
%StDev = sqrt((cov(x_sonde(~idx),y_MPD05(~idx))))
mdl = fitlm(x_Raman(~idx),y_MPD05(~idx));
Rsquared = mdl.Rsquared.Ordinary
Corr = sqrt(Rsquared) 
StDev = mdl.RMSE
num_samples = mdl.NumObservations
fit = flip(mdl.Coefficients.Estimate')


% Evaluate fit equation using polyval
xx = WV_min:0.1:WV_max;
y_est = polyval(fit,xx);

figure(7) %Frequency histogram 
%figure('Position',Scrsize);
% this is temporary (the python code should remove these 
binscatter(x_Raman,y_MPD05, [bins bins]);    
  %N=hh.NumBins;
  %hh.NumBins =[bins bins];
  %hh.XLimits = [WV_min WV_max];
  %hh.YLimits = [WV_min WV_max];
%ax.CLim = [10 90];
%ax.XGrid = 'on';
%ax.YGrid = 'on';
colormap('parula')
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
text(1, 19, ['Number of samples = ' num2str(num_samples(1))], 'FontSize', 22, 'Color', 'r')
text(1, 18, ['y = ' num2str(fit(1),3) '*x + ' num2str(fit(2),3)], 'FontSize', 22, 'Color', 'r')
text(1, 17, ['Corr = ' num2str(Corr,3)], 'FontSize', 22, 'Color', 'r')
text(1, 16, ['StDev = ' num2str(StDev,3) 'g m^{-3}'], 'FontSize', 22, 'Color', 'r')
%text(1, 15, ['R^2 = ' num2str(Rsquared,3)], 'FontSize', 22, 'Color', 'r')

hold on
plot(xx0,yy0, 'k-')  % plot the 1:1 line
%plot(xx,y_est,'r--','LineWidth',2)  % plot the least squared fit line
hold off

figure(8)
plot(x_Raman, y_MPD05, 'o')
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
text(1, 19, ['Number of samples = ' num2str(num_samples(1))], 'FontSize', 22, 'Color', 'r')
text(1, 18, ['y = ' num2str(fit(1),3) '*x + ' num2str(fit(2),3)], 'FontSize', 22, 'Color', 'r')
text(1, 17, ['Corr = ' num2str(Corr,3)], 'FontSize', 22, 'Color', 'r')
text(1, 16, ['StDev = ' num2str(StDev,3) 'g m^{-3}'], 'FontSize', 22, 'Color', 'r')
% Add trend line to plot
hold on
plot(xx0,yy0, 'k-')  % plot the 1:1 line
plot(xx,y_est,'r--','LineWidth',2) % plot the least squared fit line
hold off

cd('/Volumes/documents/WV_DIAL_data/plots/') % point to the directory where data is stored 
FigH = figure(7);
set(gca,'Fontsize',30,'Fontweight','b'); % 
set(FigH, 'PaperUnits', 'points', 'PaperPosition', Scrsize);
name=strcat(date, 'Raman_hist_multi'); 
print(FigH, name, '-dpng', '-r0') % set at the screen resolution 


cd('/Volumes/documents/WV_DIAL_data/plots/') % point to the directory where data is stored 
FigH = figure(8);
set(gca,'Fontsize',30,'Fontweight','b'); % 
set(FigH, 'PaperUnits', 'points', 'PaperPosition', Scrsize);
name=strcat(date, 'Raman_multi'); 
print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
