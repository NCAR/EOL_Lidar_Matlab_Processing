clear all; close all;

dd = pwd; % get the current path
date = '10 Sep 2019'; % Last day of a five-unit side-by-side test  
%cd('/Volumes/documents/WV_DIAL_data/MPD1_processed_data') % point to the directory where data is stored 
%load(strcat(date,'_combined.mat'))
cd('/Volumes/documents/WV_DIAL_data/MPD2_processed_data') % point to the directory where data is stored 
load(strcat(date,'_combined.mat'))
%cd('/Volumes/documents/WV_DIAL_data/MPD3_processed_data') % point to the directory where data is stored 
%load(strcat(date,'_combined.mat'))
cd('/Volumes/documents/WV_DIAL_data/MPD4_processed_data') % point to the directory where data is stored 
load(strcat(date,'_combined.mat'))
%cd('/Volumes/documents/WV_DIAL_data/MPD5_processed_data') % point to the directory where data is stored 
%load(strcat(date,'_combined.mat'))
    
WV_min = 0;
WV_max = 15;
bins = WV_max*10; % bin size is x 0.1 x 0.1 g/m^2
bin_min = 1;
bin_max = 2500;

xx{1} = reshape(MPD02.N_avg_comb,1,[]).*1e6./6.022E23.*18.015;
xx{2} = reshape(MPD04.N_avg_comb,1,[]).*1e6./6.022E23.*18.015;
%xx{2} = reshape(MPD02.N_avg_comb,1,[]).*1e6./6.022E23.*18.015;
%xx{3} = reshape(MPD03.N_avg_comb,1,[]).*1e6./6.022E23.*18.015;
%xx{4} = reshape(MPD04.N_avg_comb,1,[]).*1e6./6.022E23.*18.015;
%xx{5} = reshape(MPD05.N_avg_comb,1,[]).*1e6./6.022E23.*18.015;
xx0 = WV_min:1:WV_max;
y0 = WV_min:1:WV_max;
y90 = WV_min:1*0.9:WV_max*0.9;
y110 = WV_min:1*1.1:WV_max*1.1;

scrsz = get(0,'ScreenSize');
Scrsize=[scrsz(4)/1 scrsz(4)/1 scrsz(3)/1.5 scrsz(4)/1.5];

xfit = WV_min:0.1:WV_max;
figure('Position',Scrsize);
k=1;
m=1;

for k=1:4 %row
    for m=1:4 %column
      if (m>=k)== 1
        
        % calculate the best fit (in a least-squares sense) and the correlation coefficient
        idx = (isnan(xx{k})|isnan(xx{m+1})); %remove the NaNs
        num_samples = size(xx{k}(~idx), 2); 
        fit = polyfit(xx{k}(~idx),xx{m+1}(~idx),1);
        [Corr, P_test] = corrcoef(xx{k}(~idx),xx{m+1}(~idx))
        Cov = cov(xx{k}(~idx),xx{m+1}(~idx));
        StDev = sqrt((cov(xx{k}(~idx),xx{m+1}(~idx))))
        % Evaluate fit equation using polyval
        y_est = polyval(fit,xfit);  
          
          
        ax=subplot(4,4,((k-1)*4)+m);
        binscatter(xx{k},xx{m+1}, [bins bins]);    
        %N=h.NumBins;
        %h.NumBins =[bins bins];
        %h.XLimits = [WV_min WV_max];
        %h.YLimits = [WV_min WV_max];
        %ax.Colormap = 'parula';
        %ax.ColorScale = 'log'; 
        ax.CLim = [bin_min bin_max];
        ax.CLim = [10 1000];
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        colormap('parula')
        %colorbar(gca,'off')
     %   plot(x0,y90, 'r:')
     %   plot(x0,y110, 'r:')
        xlim([WV_min WV_max])
        ylim([WV_min WV_max])
        title(['WV x=',num2str(k),' y=', num2str(m+1)])
     %   xlabel(gca, 'absolute humidity')
     %   ylabel(gca, 'absolute humidity')
        
     % Display fit infor on graph
     text(0.5, 9.5, ['samples = ' num2str(num_samples(1))], 'FontSize', 10, 'Color', 'r')
     text(0.5, 8.75, ['y = ' num2str(fit(1),3) '*x + ' num2str(fit(2),3)], 'FontSize', 10, 'Color', 'r')
     text(0.5, 8.0, ['Corr = ' num2str(Corr(2),3)], 'FontSize', 10, 'Color', 'r')
     text(0.5, 7.25, ['StDev = ' num2str(StDev(2),3) 'g m^{-3}'], 'FontSize', 10, 'Color', 'r')
     hold on
     hold on 
     plot(xx0,y0, 'k-') % plot the 1:1 line
     plot(xfit,y_est,'r--','LineWidth',2)  % plot the least squared fit line
     hold off
 
      end
    end
end


cd('/Volumes/documents/WV_DIAL_data/plots/') % point to the directory where data is stored 
FigH = figure(1);
%set(gca,'Fontsize',30,'Fontweight','b'); % 
set(FigH, 'PaperUnits', 'points', 'PaperPosition', Scrsize);
name=strcat(date, 'Self_comparison_hist_multi'); 
print(FigH, name, '-dpng', '-r0') % set at the screen resolution 


cd(dd) % point back to original directory