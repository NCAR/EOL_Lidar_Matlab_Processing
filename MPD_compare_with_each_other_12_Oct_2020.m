clear all; close all;

dd = pwd; % get the current path
%date = '10 Apr 2019'; % Last day of a five-unit side-by-side test  
date = '12 Oct 2020'; %   


if strcmp(getenv('HOSTNAME'),'fog.eol.ucar.edu')
   serv_path = '/export/fog1/rsfdata/MPD/'; % when running on server
elseif strcmp(getenv('HOSTNAME'),'')
    serv_path = '/Users/spuler/Desktop/'; % when running on server   
else 
   serv_path = '/Volumes/eol/fog1/rsfdata/MPD/'; % 
end


cd(strcat(serv_path, 'mpd_01_processed_data/Matlab')) % point to the directory where data is stored 
load(strcat(date,'_combined.mat'))
cd(strcat(serv_path, 'mpd_02_processed_data/Matlab')) % point to the directory where data is stored  
load(strcat(date,'_combined.mat'))
cd(strcat(serv_path, 'mpd_03_processed_data/Matlab')) % point to the directory where data is stored 
load(strcat(date,'_combined.mat'))
cd(strcat(serv_path, 'mpd_04_processed_data/Matlab')) % point to the directory where data is stored 
load(strcat(date,'_combined.mat'))
cd(strcat(serv_path, 'mpd_05_processed_data/Matlab')) % point to the directory where data is stored 
load(strcat(date,'_combined.mat'))
    
WV_min = 0;
WV_max = 8;
bins = WV_max*10; % bin size is x 0.1 x 0.1 g/m^2
%bins = WV_max*40; % bin size is x 0.1 x 0.1 g/m^2
bin_min = 1;
bin_max = 4000;

% make sure they same set to same range 
  range_limit = min([size(MPD03.N_avg_comb,2) size(MPD03.N_avg_comb,2)])
  MPD03.N_avg_comb = real(MPD03.N_avg_comb(:,1:range_limit));   
  MPD04.N_avg_comb = real(MPD04.N_avg_comb(:,1:range_limit)); 

  
% just for 06 or 12 Oct 2020 remove the lowest 600m from MPD 05 during when
% WFOV receiver was blocked
MPD05.N_avg_comb(4800:7500,1:8) = NaN; 
  

 p_start = 6250;
 p_stop = p_start+500; 
 
%blank the data at 225m 
%MPD05.N_avg_comb(p_start:p_stop,1:1) = NaN;  
 
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

 MPD05.N_avg_comb(p_start:p_stop,1:1) = NaN; 

 %Scrnsize = get(0,'screensize')
 % plot Narrow water vapor in g/m^3
 figure('Position',[1,1,12560/2,1440/2])
 %Z = double(real(N_avg_comb'.*1e6./6.022E23.*18.015));  %number density in mol/cm3(1e6 cm3/m3)/(N_A mol/mole)*(18g/mole)
 % Z(isnan(Z)) = -1;
 set(gcf,'renderer','zbuffer');
 x = MPD05.time;
 y = MPD05.range./1e3;
 Z_AH = double(real(MPD05.N_avg_comb'.*1e6./6.022E23.*18.015));
 h = pcolor(x,y,Z_AH);
 set(h, 'EdgeColor', 'none');
 colorbar('EastOutside');
 axis([fix(min(x)) ceil(max(x)) 0 6])
 caxis([0 8]);
 colormap(jet)
  ylabel('Height (km, AGL)','fontweight','b','fontsize',12); 

  font_size = 28; % use this for 2019b version
  scrsz = [1  1  1024 768];
  Scrnsize = [scrsz(4)/1 scrsz(4)/1 scrsz(3)/0.30 scrsz(4)/2]; % use for ILRC really long plot
  FigH = figure(1);
  set(gca,'Fontsize',font_size,'Fontweight','b'); % use for Perdigao BAMS plots 
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', Scrnsize);
  name=strcat(date, 'H2O_multi'); 
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution 


% remove the lowest range bin which has the surface station data
% then use the loop to clip the ranges back in 75m increments
j=1
%for j = 1:10     %cutoff the 75m range bins from 1bin(75m) to 10bin(750m)
%jj
MPD01.N_avg_comb(:,1:j) = NaN; 
MPD02.N_avg_comb(:,1:j) = NaN; 
MPD03.N_avg_comb(:,1:j) = NaN; 
MPD04.N_avg_comb(:,1:j) = NaN;
MPD05.N_avg_comb(:,1:j) = NaN;

% just for 12 Oct 2020 remove the lowest 600m from MPD 05 during 
% MPD05.N_avg_comb(4800:7500,1:8) = NaN; 
 
xx{1} = real(reshape(MPD01.N_avg_comb,1,[]).*1e6./6.022E23.*18.015);
xx{2} = real(reshape(MPD02.N_avg_comb,1,[]).*1e6./6.022E23.*18.015);
xx{3} = real(reshape(MPD03.N_avg_comb,1,[]).*1e6./6.022E23.*18.015);
xx{4} = real(reshape(MPD04.N_avg_comb,1,[]).*1e6./6.022E23.*18.015);
xx{5} = real(reshape(MPD05.N_avg_comb,1,[]).*1e6./6.022E23.*18.015);

xx0 = WV_min:1:WV_max;
y0 = WV_min:1:WV_max;
y90 = WV_min:1*0.9:WV_max*0.9;
y110 = WV_min:1*1.1:WV_max*1.1;

scrsz = get(0,'ScreenSize');
Scrsize=[scrsz(4)/1 scrsz(4)/1 scrsz(3)/1.5 scrsz(4)/1.5];
Scrsize=[scrsz(4)/1 scrsz(4)/1 scrsz(3)/1 scrsz(4)/1];

xfit = WV_min:0.1:WV_max;
figure('Position',Scrsize);
k=1;
m=1;

%k=3; % these are set for just looking at MPD03 and MPD04
%m=3;

for k=1:4 %row
    for m=1:4  %column
      if (m>=k)== 1
        
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
          
  
        ax=subplot(4,4,((k-1)*4)+m);
        binscatter(xx{k},xx{m+1}, [bins bins]);    
        %N=h.NumBins;
        %h.NumBins =[bins bins];
        %h.XLimits = [WV_min WV_max];
        %h.YLimits = [WV_min WV_max];
        %ax.Colormap = 'parula';
       % ax.ColorScale = 'log'; 
        ax.CLim = [bin_min bin_max];
%        ax.CLim = [1 5000];
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        colormap('parula')
        oldcmap = colormap;
        colormap( flipud(oldcmap) );

        colorbar(gca,'off')
     %   plot(x0,y90, 'r:')
     %   plot(x0,y110, 'r:')
        xlim([WV_min WV_max])
        ylim([WV_min WV_max])
        title(['WV x=',num2str(k),' y=', num2str(m+1)])
     %   xlabel(gca, 'absolute humidity')
     %   ylabel(gca, 'absolute humidity')
        
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
    end
end


cd(strcat(serv_path, 'mpd_03_processed_data/Plots')) % point to the directory where data is stored 
FigH = figure(2);
%set(gca,'Fontsize',30,'Fontweight','b'); % 
set(FigH, 'PaperUnits', 'points', 'PaperPosition', Scrsize);
name=strcat(date, 'Self_comparison_hist_multi', num2str(j)); 
print(FigH, name, '-dpng', '-r300') % set at the screen resolution 

close all

%end

cd(dd) % point back to original directory