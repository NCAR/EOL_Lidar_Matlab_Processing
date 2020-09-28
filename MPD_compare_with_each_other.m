clear all; close all;

dd = pwd; % get the current path
%date = '10 Apr 2019'; % Last day of a five-unit side-by-side test  
date = '28 Sep 2020'; %   


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
WV_max = 10;
bins = WV_max*10; % bin size is x 0.1 x 0.1 g/m^2
%bins = WV_max*40; % bin size is x 0.1 x 0.1 g/m^2
bin_min = 1;
bin_max = 4000;

% make sure they same set to same range 
  range_limit = min([size(MPD03.N_avg_comb,2) size(MPD03.N_avg_comb,2)])
  MPD03.N_avg_comb = real(MPD03.N_avg_comb(:,1:range_limit));   
  MPD04.N_avg_comb = real(MPD04.N_avg_comb(:,1:range_limit)); 

% clip the range bin
for kk = 3:10     %cutoff the 75m range bins from 3bin(225m) to 10bin(750m)
kk
MPD01.N_avg_comb(:,1:kk) = NaN; 
MPD02.N_avg_comb(:,1:kk) = NaN; 
MPD03.N_avg_comb(:,1:kk) = NaN; 
MPD04.N_avg_comb(:,1:kk) = NaN; 
MPD05.N_avg_comb(:,1:kk) = NaN; 


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
FigH = figure(1);
%set(gca,'Fontsize',30,'Fontweight','b'); % 
set(FigH, 'PaperUnits', 'points', 'PaperPosition', Scrsize);
name=strcat(date, 'Self_comparison_hist_multi', kk); 
print(FigH, name, '-dpng', '-r300') % set at the screen resolution 

close all

end

cd(dd) % point back to original directory