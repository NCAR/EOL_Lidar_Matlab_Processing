clear all; close all;

dd = pwd; % get the current path
%date = '10 Apr 2019'; % Last day of a five-unit side-by-side test 
date = '2019-Apr-10'; %  Python date 
%date = '05 Oct 2020'; %   
date = '2022-Aug-18'; %  Python date 

% p_start = 2900; %1-Oct-2020 profile WFOV blocked
% p_start = 6200; %3-Oct-2020 profile WFOV open
% p_stop = p_start+250; 


plot_histograms = 0; 
mask_range_level = 0;

if strcmp(getenv('HOSTNAME'),'fog.eol.ucar.edu')
   serv_path = '/export/fog1/rsfdata/MPD/'; % when running on server
elseif strcmp(getenv('HOSTNAME'),'')
    serv_path = '/Users/spuler/Desktop/'; % when running on server   
else 
   serv_path = '/Volumes/eol/fog1/rsfdata/MPD/'; % 
end

serv_path = '/Volumes/eol/fog1/rsfdata/MPD/'; % 
   
cd(strcat(serv_path, 'mpd_01_processed_data/Matlab')) % point to the directory where data is stored
%cd(strcat(serv_path, '/mpd/intercomparison/mpd01')) % point to the directory where data is stored 
load(strcat(date,'_combined.mat'))

cd(strcat(serv_path, 'mpd_02_processed_data/Matlab')) % point to the directory where data is stored
%cd(strcat(serv_path, '/mpd/intercomparison/mpd02')) % point to the directory where data is stored 
load(strcat(date,'_combined.mat'))

cd(strcat(serv_path, 'mpd_03_processed_data/Matlab')) % point to the directory where data is stored
%cd(strcat(serv_path, '/mpd/intercomparison/mpd03')) % point to the directory where data is stored 
load(strcat(date,'_combined.mat'))

cd(strcat(serv_path, 'mpd_04_processed_data/Matlab')) % point to the directory where data is stored
%cd(strcat(serv_path, '/mpd/intercomparison/mpd04')) % point to the directory where data is stored 
load(strcat(date,'_combined.mat'))

cd(strcat(serv_path, 'mpd_05_processed_data/Matlab')) % point to the directory where data is stored
%cd(strcat(serv_path, '/mpd/intercomparison/mpd05')) % point to the directory where data is stored 
load(strcat(date,'_combined.mat'))
    
% just for 05 Oct 2020 remove the lowest 600m from MPD 05 during when
% WFOV receiver was blocked
% MPD05.N_avg_comb(4800:7500,1:8) = NaN;  % matlab
% MPD05.N_avg_comb(2475:3675,1:17) = NaN;   %python 
% MPD05.N_avg_comb( p_start:p_stop,1:25) = NaN;   %python 

WV_min = 0;
WV_max = 8;
bins = WV_max*10; % bin size is x 0.1 x 0.1 g/m^2
%bins = WV_max*40; % bin size is x 0.1 x 0.1 g/m^2
bin_min = 1;
bin_max = 4000;
bin_max = 2000;

% make sure they same set to same range 
%  range_limit = min([size(MPD03.N_avg_comb,2) size(MPD03.N_avg_comb,2)])
%  MPD03.N_avg_comb = real(MPD03.N_avg_comb(:,1:range_limit));   
%  MPD04.N_avg_comb = real(MPD04.N_avg_comb(:,1:range_limit)); 


% plot Narrow water vapor in g/m^3
 figure('Position',[1,1,2048,760])
 %Z = double(real(N_avg_comb'.*1e6./6.022E23.*18.015));  %number density in mol/cm3(1e6 cm3/m3)/(N_A mol/mole)*(18g/mole)
 % Z(isnan(Z)) = -1;
 set(gcf,'renderer','zbuffer');
 x = MPD05.time;
 y = MPD05.range./1e3;
 
 subplot(5,1,1)
 Z_AH = double(real(MPD01.N_avg_comb'.*1e6./6.022E23.*18.015));
 h = pcolor(x,y,Z_AH);
 set(h, 'EdgeColor', 'none');
 colorbar('EastOutside');
 axis([fix(min(x)) ceil(max(x)) 0 6])
 caxis([0 8]);
 colormap(jet)
 ylabel('Height (km, AGL)','fontweight','b','fontsize',12);
 datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
 hh = title({['MPD01, Water Vapor (g m^{-3})']},'fontweight','b','fontsize',12);
 set(gca,'Fontsize',10,'Fontweight','b'); % 
 
 subplot(5,1,2)
 Z_AH = double(real(MPD02.N_avg_comb'.*1e6./6.022E23.*18.015));
 h = pcolor(x,y,Z_AH);
 set(h, 'EdgeColor', 'none');
 colorbar('EastOutside');
 axis([fix(min(x)) ceil(max(x)) 0 6])
 caxis([0 8]);
 colormap(jet)
 ylabel('Height (km, AGL)','fontweight','b','fontsize',12);
 datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
 hh = title({['MPD02, Water Vapor (g m^{-3})']},'fontweight','b','fontsize',12);
 set(gca,'Fontsize',10,'Fontweight','b'); % 
 
 subplot(5,1,3)
 Z_AH = double(real(MPD03.N_avg_comb'.*1e6./6.022E23.*18.015));
 h = pcolor(x,y,Z_AH);
 set(h, 'EdgeColor', 'none');
 colorbar('EastOutside');
 axis([fix(min(x)) ceil(max(x)) 0 6])
 caxis([0 8]);
 colormap(jet)
 ylabel('Height (km, AGL)','fontweight','b','fontsize',12);
 datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
 hh = title({['MPD03, Water Vapor (g m^{-3})']},'fontweight','b','fontsize',12);
 set(gca,'Fontsize',10,'Fontweight','b'); % 
 
 subplot(5,1,4)
 Z_AH = double(real(MPD04.N_avg_comb'.*1e6./6.022E23.*18.015));
 h = pcolor(x,y,Z_AH);
 set(h, 'EdgeColor', 'none');
 colorbar('EastOutside');
 axis([fix(min(x)) ceil(max(x)) 0 6])
 caxis([0 8]);
 colormap(jet)
 ylabel('Height (km, AGL)','fontweight','b','fontsize',12);
 datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
 hh = title({['MPD04, Water Vapor (g m^{-3})']},'fontweight','b','fontsize',12);
 set(gca,'Fontsize',10,'Fontweight','b'); % 
 
 subplot(5,1,5)
 Z_AH = double(real(MPD05.N_avg_comb'.*1e6./6.022E23.*18.015));
 h = pcolor(x,y,Z_AH);
 set(h, 'EdgeColor', 'none');
 colorbar('EastOutside');
 axis([fix(min(x)) ceil(max(x)) 0 6])
 caxis([0 8]);
 colormap(jet)
 ylabel('Height (km, AGL)','fontweight','b','fontsize',12); 
 datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
 hh = title({['MPD05, Water Vapor (g m^{-3})']},'fontweight','b','fontsize',12);
 set(gca,'Fontsize',10,'Fontweight','b'); % 
 
cd('/Users/spuler/Desktop') % point to the directory where data is stor
FigH = figure(1);
%$set(gca,'Fontsize',30,'Fontweight','b'); % 
set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1280 800]);
name=strcat('Python_Gen5_intercomparison');
print(FigH, name, '-dpng', '-r300')  
 
  
% plot profiles 
MPD01_AH_profile = nanmean(MPD01.N_avg_comb(p_start:p_stop,:),1).*1e6./6.022E23.*18.015; 
MPD02_AH_profile = nanmean(MPD02.N_avg_comb(p_start:p_stop,:),1).*1e6./6.022E23.*18.015; 
MPD03_AH_profile = nanmean(MPD03.N_avg_comb(p_start:p_stop,:),1).*1e6./6.022E23.*18.015; 
MPD04_AH_profile = nanmean(MPD04.N_avg_comb(p_start:p_stop,:),1).*1e6./6.022E23.*18.015; 
MPD05_AH_profile = nanmean(MPD05.N_avg_comb(p_start:p_stop,:),1).*1e6./6.022E23.*18.015; 
figure(101)
plot(MPD01_AH_profile, MPD01.range./1e3, '-o', 'DisplayName','MPD01'); 
hold on
plot(MPD02_AH_profile, MPD02.range./1e3, '-o', 'DisplayName','MPD02'); 
plot(MPD03_AH_profile, MPD03.range./1e3, '-o', 'DisplayName','MPD03'); 
plot(MPD04_AH_profile, MPD04.range./1e3, '-o', 'DisplayName','MPD04'); 
plot(MPD05_AH_profile, MPD05.range./1e3, '-o', 'DisplayName','MPD05'); 
hold off
axis([0 8 0 3.5])
xlabel('Absolute Humidity (g m^{-3})','fontweight','b','fontsize',12);
ylabel('Height (km, AGL)','fontweight','b','fontsize',12);
legend('show', 'Location','southwest')


% remove the lowest range bin which has the surface station data
% then use the loop to clip the ranges back in 75m increments
%j=1
%for j = 1:10     %cutoff the 75m range bins from 1bin(75m) to 10bin(750m)
%jj

range_mask = repmat(MPD05.range, size(MPD05.time,1),1); 

MPD01.N_avg_comb(range_mask<=mask_range_level) = NaN; 
MPD02.N_avg_comb(range_mask<=mask_range_level) = NaN; 
MPD03.N_avg_comb(range_mask<=mask_range_level) = NaN; 
MPD04.N_avg_comb(range_mask<=mask_range_level) = NaN; 
MPD05.N_avg_comb(range_mask<=mask_range_level) = NaN; 


 xx{1} = real(reshape(MPD01.N_avg_comb,1,[]).*1e6./6.022E23.*18.015);
 xx{2} = real(reshape(MPD02.N_avg_comb,1,[]).*1e6./6.022E23.*18.015);
 xx{3} = real(reshape(MPD03.N_avg_comb,1,[]).*1e6./6.022E23.*18.015);
 xx{4} = real(reshape(MPD04.N_avg_comb,1,[]).*1e6./6.022E23.*18.015);
 xx{5} = real(reshape(MPD05.N_avg_comb,1,[]).*1e6./6.022E23.*18.015);

 
%xx{1}(xx{1}>10|xx{1}<-0.5)=NaN;
%xx{2}(xx{2}>10|xx{2}<-0.5)=NaN;
%xx{3}(xx{3}>10|xx{3}<-0.5)=NaN;
%xx{4}(xx{4}>10|xx{4}<-0.5)=NaN;
%xx{5}(xx{5}>10|xx{5}<-0.5)=NaN;


%   Corr_R_start = 150; 
%   Corr_R_thick = 150;
%   Corr_R_inc = 150; 
%   Corr_R_stop = 4000; % end at 6000 m
%   j= Corr_R_start;
%    
%  for j = Corr_R_start:Corr_R_inc:Corr_R_stop     %cutoff the 75m range bins from 1bin(75m) to 10bin(750m)
%    count=1
% 
%    
%    MPD01.N_avg_slice = MPD01.N_avg_comb;
%    MPD02.N_avg_slice = MPD02.N_avg_comb;
%    MPD03.N_avg_slice = MPD03.N_avg_comb;
%    MPD04.N_avg_slice = MPD04.N_avg_comb;
%    MPD05.N_avg_slice = MPD05.N_avg_comb;
%    
%    MPD01.N_avg_slice((range_mask <= j-Corr_R_thick/2) | (range_mask >= j+Corr_R_thick/2))= NaN; 
%    MPD02.N_avg_slice((range_mask <= j-Corr_R_thick/2) | (range_mask >= j+Corr_R_thick/2))= NaN; 
%    MPD03.N_avg_slice((range_mask <= j-Corr_R_thick/2) | (range_mask >= j+Corr_R_thick/2))= NaN; 
%    MPD04.N_avg_slice((range_mask <= j-Corr_R_thick/2) | (range_mask >= j+Corr_R_thick/2))= NaN; 
%    MPD05.N_avg_slice((range_mask <= j-Corr_R_thick/2) | (range_mask >= j+Corr_R_thick/2))= NaN; 
% 
%    
%    subplot(5,1,5)
%    Z_AH = double(real(MPD05.N_avg_slice'.*1e6./6.022E23.*18.015));
%    h = pcolor(x,y,Z_AH);
%    set(h, 'EdgeColor', 'none');
%    colorbar('EastOutside');
%    axis([fix(min(x)) ceil(max(x)) 0 6])
%    caxis([0 8]);
%    colormap(jet)
%    ylabel('Height (km, AGL)','fontweight','b','fontsize',12); 
%    datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
%    hh = title({['MPD05, Water Vapor (g m^{-3})']},'fontweight','b','fontsize',12);
%    set(gca,'Fontsize',10,'Fontweight','b'); %   
   
%    xx{1} = real(reshape(MPD01.N_avg_slice,1,[]).*1e6./6.022E23.*18.015);
%    xx{2} = real(reshape(MPD02.N_avg_slice,1,[]).*1e6./6.022E23.*18.015);
%    xx{3} = real(reshape(MPD03.N_avg_slice,1,[]).*1e6./6.022E23.*18.015);
%    xx{4} = real(reshape(MPD04.N_avg_slice,1,[]).*1e6./6.022E23.*18.015);
%    xx{5} = real(reshape(MPD05.N_avg_slice,1,[]).*1e6./6.022E23.*18.015);

   xx0 = WV_min:1:WV_max;
   y0 = WV_min:1:WV_max;
   y90 = WV_min:1*0.9:WV_max*0.9;
   y110 = WV_min:1*1.1:WV_max*1.1;

  %  scrsz = get(0,'ScreenSize');
  %  Scrsize=[scrsz(4)/1 scrsz(4)/1 scrsz(3)/1.5 scrsz(4)/1.5];
  %  Scrsize=[scrsz(4)/1 scrsz(4)/1 scrsz(3)/1 scrsz(4)/1];

    xfit = WV_min:0.1:WV_max;
    if plot_histograms == 1 
      figure('Position',[1 1 1200 1200]);
    end
    k=1;
    m=1;

    %k=3; % these are set for just looking at MPD03 and MPD04
    %m=3;

    for k=1:4 %row
      for m=1:4  %column
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
           ax=subplot(4,4,((k-1)*4)+m);
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
           text(WV_min+0.25, WV_max-0.5, ['samples = ' num2str(num_samples(1))], 'FontSize', 8, 'Color', 'r')
           text(WV_min+0.25, WV_max-1.0, ['y = ' num2str(fit(1),3) '*x + ' num2str(fit(2),3)], 'FontSize', 8, 'Color', 'r')
           text(WV_min+0.25, WV_max-1.5, ['Corr = ' num2str(Corr(2),3)], 'FontSize', 8, 'Color', 'r')
           text(WV_min+0.25, WV_max-2.0, ['StDev = ' num2str(StDev(2),3) 'g m^{-3}'], 'FontSize', 8, 'Color', 'r')
           hold on
           plot(xx0,y0, 'k-') % plot the 1:1 line
           plot(xfit,y_est,'r--','LineWidth',2)  % plot the least squared fit line
           hold off
 
           hp4 = get(subplot(4,4,16),'Position')
           colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.01  (hp4(2)+hp4(3))*3]) % x , y, width, height
         end
  
%        j/Corr_R_inc  
%        data{count,j/Corr_R_inc} = Corr(2);
%        count=count+1
         
         end
      end
    end
%  end
  
  
%   % Plot the data
%    Test_Range = Corr_R_start:Corr_R_inc:Corr_R_stop;
%    for k=1:1:size(data,2)
%      Mean_Cov{k} = nanmean([data{:,k}]);
%      Std_Cov{k} = nanstd([data{:,k}]);
%    end
% 
%   figure(202)
%   hold on
%   plot([Mean_Cov{:}], Test_Range/1000, 'r', 'LineWidth', 2, 'DisplayName', 'Gen4 2019 test')
%   hold on
%   eb(1) = errorbar([Mean_Cov{:}],Test_Range/1000,[Std_Cov{:}], 'horizontal', 'LineStyle', 'none', 'HandleVisibility','off');
%   set(eb, 'color', 'r', 'LineWidth', 2)
%   grid on
%   %ylim([0 3.5])
%   %xlim([0 1])
%   axis([0 1 0 4])
%   xlabel('Correlation Coefficient','fontweight','b','fontsize',12);
%   ylabel('Height (km, AGL)','fontweight','b','fontsize',12);
%   set(gca,'Fontsize',16,'Fontweight','b'); % 
%   legend

%cd(strcat(serv_path, 'mpd_03_processed_data/Plots')) % point to the directory where data is stored
FigH = figure(3);
%set(gca,'Fontsize',30,'Fontweight','b'); % 
set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 800 800]);
name=strcat(date, 'Self_comparison_hist_multi', num2str(j)); 
print(FigH, name, '-dpng', '-r300') % set at the screen resolution 


%cd(dd) % point back to original directory