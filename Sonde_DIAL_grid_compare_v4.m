clear all
close all
scrsz = get(0,'ScreenSize');
sonde_select = 1;
sonde_smooth = 150; % range smooth the sondes to match DIAL
save_figs = 0; 
gate = 75;
cut = 3; %number of std to consider as outliner


cd('/Volumes/documents/WV_DIAL_data/processed_data/') % point to the directory where data is stored 
file = uigetfile('*sonde.mat','Select the sonde error file', 'MultiSelect', 'on');
NF_top = round(5000/75);
j=1;

for j = 1:size(file,2)

 load(file{j})
 date =  strread(file{j}, '%18c', 1);

% this is a ugly way to catch sondes that go up and then back down 
if size(sonde_range,1)>2000
    sonde_range=sonde_range(1:2000,1);
    g_H2O = g_H2O(1:2000,1);
end

% use this bit of code to get the sonde data on the range spacing of the DIAL
idx = find(diff(sonde_range)<=0); % check for monotonic altitudes in the sonde
while isempty(idx) == 0
  idx = find(diff(sonde_range)<=0); % check for monotonic altitudes in the sonde
  %figure(12)
  %plot(diff(sonde_range))
  for i=1:size(idx,1)
    sonde_range(idx(i),1)= sonde_range(idx(i)+1,1)-1e-12;
  end
end

% convert sonde to native DIAL at 75 meter gates
g_H2O_grid = interp1(sonde_range, g_H2O, dial_range, 'linear'); 

% now smooth sonde data to 150 m  DIAL
%length_spatial = ones(1,2)/(2);
%g_H2O_grid_avg = filter2(length_spatial,  g_H2O_grid );
dd=pwd;
cd('/Users/spuler/Desktop/WV_DIAL/matlab');
%g_H2O_grid_avg1 = nanmoving_average(g_H2O_grid,150/gate/2,2,1); % smooth sonde to level of DIAL
%g_H2O_grid_avg2 = nanmoving_average(g_H2O_grid,300/gate/2,2,1); % smooth sonde to level of DIAL
%g_H2O_grid_avg3 = nanmoving_average(g_H2O_grid,600/gate/2,2,1); % smooth sonde to level of DIAL
%r1 = round(1500/gate); % index for smoothing range 1 (2km)
%r2 = round(2500/gate); % index for smoothing range 2 (3km)
%g_H2O_grid_avg=[g_H2O_grid_avg1(:,1:r1)';...
%  g_H2O_grid_avg2(:,r1+1:r2)';...
%  g_H2O_grid_avg3(:,r2+1:size(g_H2O_grid,2))'    ]';
%g_H2O_grid_avg = nanmoving_average(g_H2O_grid_avg,150/gate/2,2,1);
g_H2O_grid_avg = nanmoving_average(g_H2O_grid,sonde_smooth/gate/2,2,1); % smooth sonde to level of DIAL
cd(dd);

% grid everything to 75 meter gates
dial_grid_range = 0:0.075:10;
g_H2O_grid = interp1(dial_range, g_H2O_grid_avg, dial_grid_range, 'linear');
g_FF_grid = interp1(dial_range, g_FF_ave, dial_grid_range, 'linear');
g_error_FF_grid = interp1(dial_range, g_error_FF_ave, dial_grid_range, 'linear');

figure(1)
plot(g_H2O_grid_avg(1:219), dial_range(1:219), 'b-o', g_FF_ave(1:219), dial_range(1:219), 'r-o'); %,  g_NF_ave(1:219), dial_range(1:219), 'g-o')
ylim([0 6])

% use this remove empty DIAL or sonde profies (or only 1 point)
% if sonde_select == 1
   if isnan(nanmean(g_FF_ave(3:219)-g_H2O_grid_avg(3:219))) == 1 || nanstd(g_FF_ave(3:219)) == 0
  %   prompt = 'enter 1 to reject this sonde ';
  %   x = input(prompt);
  %   if x==1
  %    prompt = 'rejecting this sonde ';
  %    pause
      continue
    end
%  end

%DIAL_profile_std = nanstd(g_FF_ave(3:219))
%DIAL_profile_cov = nancov(g_FF_ave(3:219))
%pause


AH_diff = (g_FF_grid-g_H2O_grid); % DIAL-sonde
error = (g_FF_grid - g_H2O_grid)./g_H2O_grid; %DIAL-sonde/sonde
sonde_ID = str2num(ID_sonde);
all_profile_values = (error(~isnan(error)));
all_sonde_ID_values = (sonde_ID.*(ones(size(all_profile_values))));

if j == 1  % first time through do this
  error_FF = (g_FF_grid - g_H2O_grid)./g_H2O_grid; %DIAL-sonde/sonde
 % AH_error_tot = (g_FF_grid - g_H2O_grid)./g_H2O_grid; %DIAL-sonde/sonde
  AH_diff_tot = AH_diff ; % DIAL-sonde
  T_diff_tot = T_diff; % DIAL-sonde
  P_diff_tot = P_diff; % DIAL-sonde
  sonde_ID_tot = str2num(ID_sonde);
  all_profile_values_tot = (error(~isnan(error)));
  all_sonde_ID_values_tot = (sonde_ID.*(ones(size(all_profile_values))));
else
  error_FF = vertcat(error_FF, (g_FF_grid-g_H2O_grid)./g_H2O_grid);
 % AH_error_tot = vertcat(AH_error_tot, (g_FF_grid - g_H2O_grid)./g_H2O_grid); 
  AH_diff_tot = vertcat(AH_diff_tot, AH_diff); 
  T_diff_tot = vertcat(T_diff_tot, T_diff); 
  P_diff_tot = vertcat(P_diff_tot, P_diff);
  sonde_ID_tot = vertcat(sonde_ID_tot, str2num(ID_sonde));
  all_profile_values_tot = horzcat(all_profile_values, all_profile_values_tot);
  all_sonde_ID_values_tot = horzcat(all_sonde_ID_values, all_sonde_ID_values_tot);
end





figure(100) % show plot with shifted range
plot(g_FF_grid, dial_grid_range, 'r-o')%,  g_NF_grid, circshift(dial_grid_range,[1,shift/2-index]), 'g-o')
 hold on
 plot(g_H2O_grid, dial_grid_range, 'b', 'Linewidth', 2)
 hold off
ylim([0 4])
xlim([0 18])
title({[date, ' UTC']}, 'Fontsize', 30, 'Fontweight', 'b')
xlabel('Water Vapor (g/m^{3})',  'Fontsize', 30, 'Fontweight', 'b');
ylabel('Height (km, AGL)', 'Fontsize', 30, 'Fontweight', 'b');
% add error bars
  hold on
     dd = pwd;
     cd('/Users/spuler/Desktop/WV_DIAL/Matlab/')
     hb = herrorbar(g_FF_grid,  dial_grid_range, g_error_FF_grid(end,:)); 
     set(hb, 'Color', 'r');
     cd(dd)
  hold off
  %legend('Far Field', 'Near Field', 'Sonde' ,'Location', 'NorthEast');
  legend('Far Field', 'Sonde' ,'Location', 'NorthEast');
  grid on
  set(gca, 'FontSize', 22) 


  if save_figs==1
     % print the figures to file 
     dd = pwd;
     cd('/Volumes/documents/WV_DIAL_data/plots/') % point to the directory where data is stored 
     FigH = figure(100);
     set(FigH, 'PaperUnits', 'points', 'PaperPosition', [scrsz(4)/2 scrsz(4)/10 scrsz(3)/2 scrsz(4)/1.5]);
     name=strcat(date,'_WV_2CH'); 
     print(FigH, name, '-dpng', '-r300') % set the resolution as 300 dpi
     cd(dd)
  end

 
figure(200)
plot(AH_diff_tot , dial_grid_range, 'ro', 'MarkerSize', 3); %error_NF , dial_grid_range, 'gx')
xlim([-4 4])
ylim([0 6])
legend('\Delta AH (g/m^{3})')
xlabel('AH Diff (K) DIAL-sonde',  'Fontsize', 30, 'Fontweight', 'b');

figure(201)
plot(T_diff_tot , dial_range, 'ro', 'MarkerSize', 3); %error_NF , dial_grid_range, 'gx')
%xlim([-1 1])
ylim([0 6])
legend('\Delta T (K)')
xlabel('T Diff (K) DIAL-sonde',  'Fontsize', 30, 'Fontweight', 'b');

figure(202)
plot(P_diff_tot , dial_range, 'ro', 'MarkerSize', 3); %error_NF , dial_grid_range, 'gx')
%hist_range=repmat(dial_range,[size(P_diff_tot,1),1]);
%histogram2(P_diff_tot,hist_range,'DisplayStyle','tile','ShowEmptyBins','on'); 
%xlim([-1 1])
ylim([0 6])
legend('\Delta P (atm)')
xlabel('P Diff (atm) DIAL-sonde',  'Fontsize', 30, 'Fontweight', 'b');

profile_diff_mean = nanmedian(AH_diff_tot(:,2:4500/gate),2); % 4.5km and down 75m (ignore suface station) 
%profile_diff_mean = nanmedian(AH_diff_tot,2);
figure(203)
plot(sonde_ID_tot, profile_diff_mean, 'ro', 'MarkerSize', 3); 
%hist_range=repmat(dial_range,[size(P_diff_tot,1),1]);
%histogram2(P_diff_tot,hist_range,'DisplayStyle','tile','ShowEmptyBins','on'); 
%xlim([-1 1])
%ylim([0 6])
%legend('\Delta P (atm)')
xlabel('mean column AH Diff',  'Fontsize', 20, 'Fontweight', 'b');


profile_error_mean = nanmedian(error_FF(:,2:4500/gate),2); % 4.5km and down 75m (ignore suface station)
%profile_error_mean = nanmedian(AH_diff_tot,2);
figure(204)
plot(sonde_ID_tot, profile_error_mean, 'ro', 'MarkerSize', 3); 
%hist_range=repmat(dial_range,[size(P_diff_tot,1),1]);
%histogram2(P_diff_tot,hist_range,'DisplayStyle','tile','ShowEmptyBins','on'); 
%xlim([-1 1])
%ylim([0 6])
%legend('\Delta P (atm)')
xlabel('mean column relative error',  'Fontsize', 20, 'Fontweight', 'b');

 label = num2str(profile_error_mean(end));
 figure(205)
 plot(error_FF(end,1:40), dial_grid_range(end,1:40))
 text(double(error_FF(1,6)), double(dial_grid_range(1,6)), label);
 grid on
 xlabel('mean column relative error',  'Fontsize', 20, 'Fontweight', 'b');

 label = num2str(profile_diff_mean(end));
 figure(206)
 plot(AH_diff_tot(end,1:40), dial_grid_range(end,1:40))
 text(double(AH_diff_tot(1,6)), double(dial_grid_range(1,6)), label);
 grid on
 xlabel('mean column AH Diff',  'Fontsize', 20, 'Fontweight', 'b');
 
 %pause
end


%spuler = [sonde_ID_tot profile_error_mean];
%edges=  min(sonde_ID_tot):1:max(sonde_ID_tot)+1;



spuler = [all_sonde_ID_values_tot; all_profile_values_tot]';
low = (spuler(:,2)<3*nanstd(spuler(:,2)));
high = ((spuler(:,2))>-3*nanstd(spuler(:,2))); 
spuler_C = [spuler((low)&(high),1), spuler((low)&(high),2)];


edges=  min(all_sonde_ID_values_tot):1:max(all_sonde_ID_values_tot)+1;
N = histcounts(all_sonde_ID_values_tot, edges);
[N,bin] = histcounts(all_sonde_ID_values_tot, edges);
sonde_ID_counts = N(N>0);
sonde_ID_bin = bin(N>0);
%spuler_mean = sonde_ID_bin';
%spuler_std = sonde_ID_bin';
i = 1;
for i=1:size(sonde_ID_bin,2)
spuler(spuler(:,1)==sonde_ID_bin(i)) = i;
dummy = spuler(spuler(:,1)==i,2); 
dummy2 = dummy( ((dummy)<3*nanstd(dummy))  & ((dummy)>-3*nanstd(dummy)) );
spuler_mean(i) = nanmean(dummy2);
%spuler_std(i) = nanstd(spuler(spuler(:,1)==i,2));
spuler_std(i) = std(dummy2,'omitnan');
spuler_skew(i) = skewness(dummy2,0);
spuler_kurt(i) = kurtosis(dummy2);
end
spulertest= [spuler_mean; spuler_std; spuler_skew; spuler_kurt];
 
figure(504)
plot(spuler(:,1), spuler(:,2), 'ro', 'MarkerSize', 3); 
%ylabel('mean AH Diff (g/m^3) DIAL-sonde',  'Fontsize', 20, 'Fontweight', 'b');
ylabel('Relative error abs humidty [DIAL-sonde/sonde]',  'Fontsize', 20, 'Fontweight', 'b');
xlabel('sonde batch number',  'Fontsize', 20, 'Fontweight', 'b');
ylim([-1 1])
xlim([0 size(sonde_ID_bin,2)+1])
grid on

x = 1:size(sonde_ID_bin,2);
labels = cellstr(num2str([sonde_ID_bin; sonde_ID_counts]'));
figure(505)
%plot(x, spuler_mean, 'ro', 'MarkerSize', 3); 
errorbar(x, spuler_mean, spuler_std, 'ro', 'MarkerSize', 3)
text(double(x), double(spuler_mean)+0.05, labels, 'VerticalAlignment', 'top', ...
   'HorizontalAlignment', 'right', 'Rotation', 90);
%ylabel('Diff absoulte humidty (g/m^3) DIAL-sonde',  'Fontsize', 20, 'Fontweight', 'b');
ylabel('Relative error abs humidty [DIAL-sonde/sonde]',  'Fontsize', 20, 'Fontweight', 'b');
xlabel('sonde batch number',  'Fontsize', 20, 'Fontweight', 'b');
ylim([-1 1])
ylim([-0.25 0.25])
xlim([0 size(sonde_ID_bin,2)+1])
grid on

figure(507)
boxplot(spuler_C(:,2), spuler_C(:,1), 'Notch', 'on', ...
    'DataLim', [-0.5, 0.5], 'ExtremeMode', 'clip'); %, 'PlotStyle', 'compact')
grid on
ylabel('Relative error abs humidty [DIAL-sonde/sonde]',  'Fontsize', 20, 'Fontweight', 'b');
xlabel('sonde batch number',  'Fontsize', 20, 'Fontweight', 'b');

hist_range=repmat(dial_range,[size(P_diff_tot,1),1]);
hist_grid_range=repmat(dial_grid_range,[size(AH_diff_tot,1),1]);

yedges = 0:0.250:6;
figure(300)
histogram2(AH_diff_tot,hist_grid_range,'YBinEdges',yedges,'DisplayStyle','tile','ShowEmptyBins','on', 'Normalization', 'probability'); 
xlim([-4 4])
%ylim([0 6])
legend('\Delta AH (g/m^{3})')
xlabel('AH Diff (K) DIAL-sonde',  'Fontsize', 20, 'Fontweight', 'b');
ylabel('Height (km, AGL)', 'Fontsize', 20, 'Fontweight', 'b');

xedges = -15:1:15;
figure(301)
histogram2(T_diff_tot,hist_range,'XBinEdges',xedges, 'YBinEdges',yedges, 'DisplayStyle','tile','ShowEmptyBins','on'); 
%ylim([0 6])
%xlim([-15 15])
legend('\Delta T (K)')
xlabel('T Diff (K) DIAL-sonde',  'Fontsize', 20, 'Fontweight', 'b');
ylabel('Height (km, AGL)', 'Fontsize', 20, 'Fontweight', 'b');

figure(302)
histogram2(P_diff_tot,hist_range,'YBinEdges',yedges,'DisplayStyle','tile','ShowEmptyBins','on'); 
%xlim([-0.01 0.01])
%ylim([0 6])
legend('\Delta P (atm)')
xlabel('P Diff (atm) DIAL-sonde',  'Fontsize', 20, 'Fontweight', 'b');
ylabel('Height (km, AGL)', 'Fontsize', 20, 'Fontweight', 'b');

% remove 3 standard deviations as outliers 
 outlier = nanstd(error_FF).*cut;
 outlier2 = repmat(outlier,size(error_FF,1),1);
 mean = nanmean(error_FF);
 mean2 = repmat(mean,size(error_FF,1),1);
 for j =1:size(error_FF,1)
  error_FF(error_FF-mean2 < -outlier2) = NaN;
  error_FF(error_FF-mean2 > outlier2) = NaN;
 end
 outlier = nanstd(error_FF).*cut;
 outlier2 = repmat(outlier,size(error_FF,1),1);
 mean = nanmean(error_FF);
 mean2 = repmat(mean,size(error_FF,1),1);
 for j =1:size(error_FF,1)
  error_FF(error_FF-mean2 < -outlier2) = NaN;
  error_FF(error_FF-mean2 > outlier2) = NaN;
 end
 outlier = nanstd(error_FF).*cut;
 outlier2 = repmat(outlier,size(error_FF,1),1);
 mean = nanmean(error_FF);
 mean2 = repmat(mean,size(error_FF,1),1);
 for j =1:size(error_FF,1)
  error_FF(error_FF-mean2 < -outlier2) = NaN;
  error_FF(error_FF-mean2 > outlier2) = NaN;
 end
 
%for j =1:size(error_FF,1)
% error_FF(j, error_FF(j,:) < nanstd(error_FF,0)*-cut) = NaN;
% error_FF(j, error_FF(j,:) > nanstd(error_FF,0)*cut) = NaN;
%end
%for j =1:size(error_FF,1)
% error_FF(j, error_FF(j,:) < nanstd(error_FF,0)*-cut) = NaN;
% error_FF(j, error_FF(j,:) > nanstd(error_FF,0)*cut) = NaN;
%end
%for j =1:size(error_FF,1)
% error_FF(j, error_FF(j,:) < nanstd(error_FF,0)*-cut) = NaN;
% error_FF(j, error_FF(j,:) > nanstd(error_FF,0)*cut) = NaN;
%end

figure(22)
plot(error_FF , dial_grid_range, 'ro', 'MarkerSize', 3) ; %, error_NF , dial_grid_range, 'gx')
xlim([-1 1])
ylim([0 6])


summed_FF = (nansum(error_FF,1));
summed_diff_FF = (nansum(AH_diff_tot,1));
sq_summed_FF = (nansum(error_FF.^2,1));
NaN_number_FF =  ~isnan(error_FF);
number_FF = sum(NaN_number_FF,1);
NaN_diff_number_FF =  ~isnan(AH_diff_tot);
diff_number_FF = sum(NaN_diff_number_FF,1);


%summed_NF = (nansum(error_NF,1));
%sq_summed_NF = (nansum(error_NF.^2,1));
%NaN_number_NF =  ~isnan(error_NF);
%number_NF = sum(NaN_number_NF,1);

RMS_error_FF = sqrt(sq_summed_FF./number_FF)*100;
%RMS_error_NF = sqrt(sq_summed_NF./number_NF)*100;

RMS_sum_error_FF = nansum(RMS_error_FF) 

figure(3)
%plot((nanmean(error_FF,1))*100 , dial_range(1:219), 'r', nanmean(error_NF(:,1:NF_top),1)*100 , dial_range(:,1:NF_top), 'g', 'LineWidth', 2)
plot(summed_FF./number_FF*100 , dial_grid_range, 'r','LineWidth', 2)
%summed_NF./number_NF*100 , dial_grid_range, 'g', 'LineWidth', 2)
%legend('Far Field', 'Near Field', 'Location', 'NorthEast');
%legend('Far Field', 'Location', 'NorthEast');
hold on
plot(nanmean(error_FF,1)*100+nanstd(error_FF,1)*100 , dial_grid_range, 'r--'); %, nanmean(error_NF,1)*100+nanstd(error_NF,1)*100 , dial_grid_range, 'g--')
plot(nanmean(error_FF,1)*100-nanstd(error_FF,1)*100 , dial_grid_range, 'r--');%, nanmean(error_NF,1)*100-nanstd(error_NF,1)*100 , dial_grid_range, 'g--')
%plot(nanmedian(error_FF,1)+nanstd(error_FF,1) , dial_range, 'r--', nanmedian(error_NF(:,1:NF_top),1)+nanstd(error_NF(:,1:NF_top),1) , dial_range(:,1:NF_top), 'g--')
%plot(nanmedian(error_FF,1)-nanstd(error_FF,1) , dial_range, 'r--', nanmedian(error_NF(:,1:NF_top),1)-nanstd(error_NF(:,1:NF_top),1) , dial_range(:,1:NF_top), 'g--')
hold off
xlim([-50 50])
grid on
ylim([0 4.5])
%title('mean relative error', 'Fontsize', 20, 'Fontweight', 'b')
xlabel('mean relative error, %',  'Fontsize', 30, 'Fontweight', 'b');
ylabel('Height (km, AGL)', 'Fontsize', 30, 'Fontweight', 'b');
set(gca, 'FontSize', 22) 
xData=[-50; -15; -5; 0; 5; 15; 50];   
set(gca, 'XTick',  xData)

figure(4)
plot(RMS_error_FF , dial_grid_range, 'r',  'LineWidth', 2)%
%RMS_error_NF , dial_grid_range, 'g', 'LineWidth', 2)
%plot(nanmedian(error_FF,1) , dial_range, 'r', nanmedian(error_NF(:,1:NF_top),1) , dial_range(:,1:NF_top), 'g',  'LineWidth', 2)
%legend('Far Field', 'Near Field', 'Location', 'Northeast');
xlim([0 200])
grid on
ylim([0 4])
xlabel('rms relative error, %',  'Fontsize', 30, 'Fontweight', 'b');
ylabel('Height (km, AGL)', 'Fontsize', 30, 'Fontweight', 'b');
set(gca, 'FontSize', 22) 
xData=[0; 50; 100; 150; 200];
set(gca, 'XTick',  xData)

figure(31)
plot(summed_diff_FF./diff_number_FF , dial_grid_range, 'r','LineWidth', 2)
xlim([-4 4])
grid on
ylim([0 4.5])
%title('mean relative error', 'Fontsize', 20, 'Fontweight', 'b')
xlabel('WV Diff (g/m^3)',  'Fontsize', 30, 'Fontweight', 'b');
ylabel('Height (km, AGL)', 'Fontsize', 30, 'Fontweight', 'b');
set(gca, 'FontSize', 18) 
xData=[-4; -3; -2; -1; 0; 1; 2; 3; 4];   
set(gca, 'XTick',  xData)


if save_figs==1
  scrsz = get(0,'ScreenSize');
  FigH = figure(506);
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', [scrsz(4)/2 scrsz(4)/10 scrsz(3)/2 scrsz(4)/1.5]);
  name=strcat('mean_sonde_error'); 
  print(FigH, name, '-dpng', '-r300') % set the resolution as 300 dpi

  scrsz = get(0,'ScreenSize');
  FigH = figure(31);
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', [scrsz(4)/2 scrsz(4)/10 scrsz(3)/2 scrsz(4)/1.5]);
  name=strcat('DIAL_sonde_diff'); 
  print(FigH, name, '-dpng', '-r300') % set the resolution as 300 dpi

  FigH = figure(4);
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', [scrsz(4)/2 scrsz(4)/10 scrsz(3)/2 scrsz(4)/1.5]);
  name=strcat('rms_sonde_error'); 
  print(FigH, name, '-dpng', '-r300') % set the resolution as 300 dpi
end

