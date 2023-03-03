clear all; close all;

dd = pwd; % get the current path
date = '22 Aug 2022'; % Last day of a collocated test 

profile_time = '21 Aug 2022 18:00';
plot_histograms = 1; 
mask_range_level = 0;
plot_path = '/Users/spuler/Desktop/mpd/Plots/';

if strcmp(getenv('HOSTNAME'),'fog.eol.ucar.edu')
   serv_path = '/export/fog1/rsfdata/MPD/'; % when running on server
elseif strcmp(getenv('HOSTNAME'),'')
    serv_path = '/Users/spuler/Desktop/'; % when running on server   
else 
   serv_path = '/Volumes/eol/fog1/rsfdata/MPD/'; % 
end

serv_path = '/Volumes/eol/fog1/rsfdata/MPD/'; % 
   
% cd(strcat(serv_path, 'mpd_01_processed_data/Matlab')) % point to the directory where data is stored

% load(strcat('MPD01_',date,'_Matlab_combined.mat'))

cd(strcat(serv_path, 'mpd_02_processed_data/Matlab')) % point to the directory where data is stored
load(strcat('MPD02_',date,'_Matlab_combined.mat'))

cd(strcat(serv_path, 'mpd_03_processed_data/Matlab')) % point to the directory where data is stored
load(strcat('MPD03_',date,'_Matlab_combined.mat'))

cd(strcat(serv_path, 'mpd_04_processed_data/Matlab')) % point to the directory where data is stored
load(strcat('MPD04_',date,'_Matlab_combined.mat'))

% cd(strcat(serv_path, 'mpd_05_processed_data/Matlab')) % point to the directory where data is stored
% load(strcat('MPD05_',date,'_Matlab_combined.mat'))
    
% just for 05 Oct 2020 remove the lowest 600m from MPD 05 during when
% WFOV receiver was blocked
% MPD05.N_avg_comb(4800:7500,1:8) = NaN;  % matlab
% MPD05.N_avg_comb(2475:3675,1:17) = NaN;   %python 
% MPD05.N_avg_comb( p_start:p_stop,1:25) = NaN;   %python 

WV_min = 0;
WV_max = 30;
bins = WV_max*10; % bin size is x 0.1 x 0.1 g/m^2
if bins > 250
    bins = 250;
end
bin_min = 1;
bin_max = 250;

% plot Narrow water vapor in g/m^3
 figure('Position',[1,1,2048,760])
 %Z = double(real(N_avg_comb'.*1e6./6.022E23.*18.015));  %number density in mol/cm3(1e6 cm3/m3)/(N_A mol/mole)*(18g/mole)
 % Z(isnan(Z)) = -1;
 set(gcf,'renderer','zbuffer');
 x = MPD02.time;
 y = MPD02.range./1e3;
 
 subplot(3,1,1)
 Z_AH = double(real(MPD02.N_avg_comb'.*1e6./6.022E23.*18.015));
 h = pcolor(x,y,Z_AH);
 set(h, 'EdgeColor', 'none');
 colorbar('EastOutside');
 axis([fix(min(x)) ceil(max(x)) 0 6])
 caxis([0 25]);
 colormap(jet)
 ylabel('Height (km, AGL)','fontweight','b','fontsize',12);
 datetick('x','dd-mmm-yy HH:MM','keeplimits', 'keepticks');
 hh = title({['MPD02, Water Vapor (g m^{-3})']},'fontweight','b','fontsize',12);
 set(gca,'Fontsize',10,'Fontweight','b'); % 
 
 subplot(3,1,2)
 Z_AH = double(real(MPD03.N_avg_comb'.*1e6./6.022E23.*18.015));
 h = pcolor(x,y,Z_AH);
 set(h, 'EdgeColor', 'none');
 colorbar('EastOutside');
 axis([fix(min(x)) ceil(max(x)) 0 6])
 caxis([0 25]);
 colormap(jet)
 ylabel('Height (km, AGL)','fontweight','b','fontsize',12);
 datetick('x','dd-mmm-yy HH:MM','keeplimits', 'keepticks');
 hh = title({['MPD03, Water Vapor (g m^{-3})']},'fontweight','b','fontsize',12);
 set(gca,'Fontsize',10,'Fontweight','b'); % 
 
 subplot(3,1,3)
 Z_AH = double(real(MPD04.N_avg_comb'.*1e6./6.022E23.*18.015));
 h = pcolor(x,y,Z_AH);
 set(h, 'EdgeColor', 'none');
 colorbar('EastOutside');
 axis([fix(min(x)) ceil(max(x)) 0 6])
 caxis([0 25]);
 colormap(jet)
 ylabel('Height (km, AGL)','fontweight','b','fontsize',12);
 datetick('x','dd-mmm-yy HH:MM','keeplimits', 'keepticks');
 hh = title({['MPD04, Water Vapor (g m^{-3})']},'fontweight','b','fontsize',12);
 set(gca,'Fontsize',10,'Fontweight','b'); % 
 

cd('/Users/spuler/Desktop') % point to the directory where data is stor
FigH = figure(1);
%$set(gca,'Fontsize',30,'Fontweight','b'); % 
set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1280 800]);
name=strcat('Python_Gen5_intercomparison');
print(FigH, name, '-dpng', '-r300')  
 
  
% plot profiles 
 p_start = find(MPD04.time >= datenum(profile_time), 1, 'first');
 p_stop = p_start+250;
 MPD02_AH_profile = mean(MPD02.N_avg_comb(p_start:p_stop,:),1,'omitnan').*1e6./6.022E23.*18.015; 
 MPD03_AH_profile = mean(MPD03.N_avg_comb(p_start:p_stop,:),1, 'omitnan').*1e6./6.022E23.*18.015; 
 MPD04_AH_profile = mean(MPD04.N_avg_comb(p_start:p_stop,:),1, 'omitnan').*1e6./6.022E23.*18.015; 
 figure(3)
 plot(MPD02_AH_profile, MPD02.range./1e3, '-o', 'DisplayName','MPD02'); 
 hold on
 plot(MPD03_AH_profile, MPD03.range./1e3, '-o', 'DisplayName','MPD03'); 
 plot(MPD04_AH_profile, MPD04.range./1e3, '-o', 'DisplayName','MPD04'); 
 hold off
 axis([0 25 0 4])
 xlabel('Absolute Humidity (g m^{-3})','fontweight','b','fontsize',12);
 ylabel('Height (km, AGL)','fontweight','b','fontsize',12);
 legend('show', 'Location','southwest')
 title(profile_time)


% remove the lowest range bin which has the surface station data
% then use the loop to clip the ranges back in 75m increments
j=1
%for j = 1:10     %cutoff the 75m range bins from 1bin(75m) to 10bin(750m)
%jj

range_mask = repmat(MPD02.range, size(MPD02.time,1),1); 

MPD02.N_avg_comb(range_mask<=mask_range_level) = NaN; 
MPD03.N_avg_comb(range_mask<=mask_range_level) = NaN; 
MPD04.N_avg_comb(range_mask<=mask_range_level) = NaN; 



 xx{1} = real(reshape(MPD02.N_avg_comb,1,[]).*1e6./6.022E23.*18.015);
 xx{2} = real(reshape(MPD03.N_avg_comb,1,[]).*1e6./6.022E23.*18.015);
 xx{3} = real(reshape(MPD04.N_avg_comb,1,[]).*1e6./6.022E23.*18.015);

 xx0 = WV_min:1:WV_max;
 y0 = WV_min:1:WV_max;
 y90 = WV_min:1*0.9:WV_max*0.9;
 y110 = WV_min:1*1.1:WV_max*1.1;

  xfit = WV_min:0.1:WV_max;

  if plot_histograms == 1 
      figure('Position',[1 1 1200 1200]);
   end
  
    %k=3; % these are set for just looking at MPD03 and MPD04
    %m=3;

    for k=1:2 %row
      for m=1:2  %column
         if (m>=k)== 1
          
        %  coarse removal of outliers
          xx{k}(xx{k}<-WV_max | xx{k}>WV_max)= NaN; % remove points outside of plot range 
          xx{m+1}(xx{m+1}<-WV_max | xx{m+1}>WV_max)=NaN; % remove points outside of plot range 
        %
         % calculate the best fit (in a least-squares sense) and the correlation coefficient
          idx = (isnan(xx{k})|isnan(xx{m+1})); %remove the NaNs
          num_samples = size(xx{k}(~idx), 2); 
          X = [xx{k}(~idx)' xx{m+1}(~idx)'];
          fit = polyfit(X(:,1),X(:,2),1);
          [Corr, P_test] = corrcoef(X)
          Cov = cov(X)
          %[R,s1] = corrcov(Cov) %statistical toolbox method to find R and std in one step
          StDev = sqrt((cov(X)))
          %StDev = sqrt(diag(Cov))
          % Evaluate fit equation using polyval
          y_est = polyval(fit,xfit);  
           
         if plot_histograms == 1 
           ax=subplot(2,2,((k-1)*2)+m);
           binscatter(xx{k},xx{m+1}, [bins bins]);    
           %N=h.NumBins;
           %h.NumBins =[bins bins];
           %h.XLimits = [WV_min WV_max];
           %h.YLimits = [WV_min WV_max];
           %ax.Colormap = 'parula';
          % ax.ColorScale = 'log'; 
           ax.CLim = [bin_min bin_max];
          % ax.CLim = [1 5000];
           ax.XGrid = 'on';
           ax.YGrid = 'on';
           colormap('parula')
           oldcmap = colormap;
           colormap( flipud(oldcmap) );

           colorbar(gca,'off')
         % plot(x0,y90, 'r:')
         % plot(x0,y110, 'r:')
           xlim([WV_min WV_max])
           ylim([WV_min WV_max])
           title(['WV x=',num2str(k),' y=', num2str(m+1)])
         %  title(['WV x=',num2str(k),' y=', num2str(m+1), ', range ', num2str(j)])
         % xlabel(gca, 'absolute humidity')
         % ylabel(gca, 'absolute humidity')
        
         % Display fit infor on graph
           text(WV_min+0.25, WV_max-0.75, ['samples = ' num2str(num_samples(1))], 'FontSize', 8, 'Color', 'r')
           text(WV_min+0.25, WV_max-1.5, ['y = ' num2str(fit(1),3) '*x + ' num2str(fit(2),3)], 'FontSize', 8, 'Color', 'r')
           text(WV_min+0.25, WV_max-2.25, ['Corr = ' num2str(Corr(2),3)], 'FontSize', 8, 'Color', 'r')
           text(WV_min+0.25, WV_max-3.0, ['StDev = ' num2str(StDev(2),3) 'g m^{-3}'], 'FontSize', 8, 'Color', 'r')
           hold on
           plot(xx0,y0, 'k-') % plot the 1:1 line
           plot(xfit,y_est,'r--','LineWidth',2)  % plot the least squared fit line
           hold off
 
           hp4 = get(subplot(2,2,4),'Position')
           colorbar('Position', [hp4(1)+hp4(3)+0.09  hp4(2)  0.01  (hp4(2)+hp4(3))*2]) % x , y, width, height
         end
  
%        j/Corr_R_inc  
%        data{count,j/Corr_R_inc} = Corr(2);
%        count=count+1
         
         end
      end
    end
%  end
  
  
cd(plot_path) % point to the directory
FigH = figure(1);
set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 800 800]);
name=strcat(date, 'Self_comparison'); 
print(FigH, name, '-dpng', '-r300') % set at the screen resolution 

FigH = figure(2);
set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 800 800]);
name=strcat(date, 'Self_comparison_hist'); 
print(FigH, name, '-dpng', '-r300') % set at the screen resolution 

FigH = figure(3);
set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 400 400]);
name=strcat(date, 'Self_comparison_profile'); 
print(FigH, name, '-dpng', '-r300') % set at the screen resolution 


