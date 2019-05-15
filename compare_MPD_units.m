clear all; close all;

dd = pwd; % get the current path
date = '10 Apr 2019'; % Last day of a five-unit side-by-side test  
cd('/Volumes/documents/WV_DIAL_data/MPD1_processed_data') % point to the directory where data is stored 
load(strcat(date,'_combined.mat'))
cd('/Volumes/documents/WV_DIAL_data/MPD2_processed_data') % point to the directory where data is stored 
load(strcat(date,'_combined.mat'))
cd('/Volumes/documents/WV_DIAL_data/MPD3_processed_data') % point to the directory where data is stored 
load(strcat(date,'_combined.mat'))
cd('/Volumes/documents/WV_DIAL_data/MPD4_processed_data') % point to the directory where data is stored 
load(strcat(date,'_combined.mat'))
cd('/Volumes/documents/WV_DIAL_data/MPD5_processed_data') % point to the directory where data is stored 
load(strcat(date,'_combined.mat'))

WV_min = 0;
WV_max = 6;
bins = 150;
bin_min = 1;
bin_max = 200;

x{1} = reshape(MPD01.N_avg_comb,1,[]).*1e6./6.022E23.*18.015;
x{2} = reshape(MPD02.N_avg_comb,1,[]).*1e6./6.022E23.*18.015;
x{3} = reshape(MPD03.N_avg_comb,1,[]).*1e6./6.022E23.*18.015;
x{4} = reshape(MPD04.N_avg_comb,1,[]).*1e6./6.022E23.*18.015;
x{5} = reshape(MPD05.N_avg_comb,1,[]).*1e6./6.022E23.*18.015;
x0 = WV_min:1:WV_max;
y0 = WV_min:1:WV_max;
y90 = WV_min:1*0.9:WV_max*0.9;
y110 = WV_min:1*1.1:WV_max*1.1;

scrsz = get(0,'ScreenSize');
size=[scrsz(4)/1 scrsz(4)/1 scrsz(3)/1.5 scrsz(4)/1.5];

f=figure('Position',size);
f.Color = 'w';
for k=1:4 %row
    for m=1:4 %column
      if (m>=k)== 1
        ax=subplot(4,4,((k-1)*4)+m);
        h=binscatter(x{k},x{m+1});    
        N= h.NumBins;
        h.NumBins =[bins bins];
        h.XLimits = [WV_min WV_max];
        h.YLimits = [WV_min WV_max];
        %ax.Colormap = 'parula';
        ax.ColorScale = 'log'; 
        ax.CLim = [bin_min bin_max];
        ax.XGrid = 'on';
        ax.YGrid = 'on';
        colormap(gca,'parula')
        colorbar(gca,'off')
        hold on 
        plot(x0,y0, 'r--')
     %   plot(x0,y90, 'r:')
     %   plot(x0,y110, 'r:')
        xlim(gca,[WV_min WV_max])
        ylim(gca,[WV_min WV_max])
        title(gca, ['WV x=',num2str(k),' y=', num2str(m+1)])
     %   xlabel(gca, 'absolute humidity')
     %   ylabel(gca, 'absolute humidity')
        hold off
      end
    end
end

cd(dd) % point back to original directory