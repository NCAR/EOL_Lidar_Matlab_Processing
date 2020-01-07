%clear all; close all;

dd = pwd; % get the current path
Raman_date = '29 Apr 2019'; % Last day of a five-unit side-by-side test  
cd('/Volumes/documents/WV_DIAL_data/MPD5_processed_data') % point to the directory where data is stored 
load(strcat(Raman_date,'_combined.mat'))

%   Grid Raman lidar data to the MPD05 range and time
    R.range = double(Raman_alt*1000)';
    R.time = Raman_duration;
    R.N = double(Raman_AH./1e6.*6.022E23./18.015)';
    [xq,yq] = meshgrid(MPD05.range, MPD05.time);
    vq = griddata(R.range,R.time,R.N, xq, yq, 'nearest');
    Raman.N_avg_comb=vq;
%   truncate data set to include only last day of the 5 days
    Raman.N_avg_comb = Raman.N_avg_comb(720-144:end,1:100);
    MPD05.N_avg_comb = MPD05.N_avg_comb(720-144:end,1:100);
%   plot the gridded data
    figure(201)
    xp = MPD05.time(720-144:end);
    yp = MPD05.range(1:100);
    zp = (Raman.N_avg_comb.*1e6./6.022E23.*18.015)';
    set(gcf,'renderer','zbuffer');
    h = pcolor(xp, yp, zp);
    set(h, 'EdgeColor', 'none');
    colorbar('EastOutside');
    colormap(jet)
    caxis([0 18]);
    ylim([0 6000])
    title({[Raman_date,' SGP Raman Lidar absolute humidity']},...
       'fontweight','b','fontsize',font_size)
   ylabel('range (km)','fontweight','b','fontsize',font_size); 
   datetick('x','HH:MM', 'keeplimits');%, 'keepticks');
   xlabel('Time (UTC)','fontweight','b','fontsize',font_size);

    
WV_min = 0;
WV_max = 16;
bins = WV_max*10; % bin size is x 0.1 x 0.1 g/m^2
bin_min = 1;
bin_max = 100;

xx{1} = reshape(Raman.N_avg_comb,1,[]).*1e6./6.022E23.*18.015;
xx{2} = reshape(MPD05.N_avg_comb,1,[]).*1e6./6.022E23.*18.015;
xx0 = WV_min:1:WV_max;
yy0 = WV_min:1:WV_max;

scrsz = get(0,'ScreenSize');
Scrsize=[scrsz(4)/1 scrsz(4)/1 scrsz(3)/1.5 scrsz(4)/1.5];
font_size = 14;

figure('Position',Scrsize);
binscatter(xx{2},xx{1}, [bins bins]);    
%N=hh.NumBins;
%hh.NumBins =[bins bins];
%hh.XLimits = [WV_min WV_max];
%hh.YLimits = [WV_min WV_max];
ax.CLim = [bin_min bin_max];
ax.XGrid = 'on';
ax.YGrid = 'on';
colormap('parula')
c=colorbar;
c.Label.String = 'Sample pairs per bin';
hold on 
plot(xx0,yy0, 'r--')
hold off
title(['WV x=MPD05, ', ' y=Raman Lidar'])
xlim([WV_min WV_max]);
pause(0.001) % there is something odd about this all running at once
set(gca,'Fontsize',font_size,'Fontweight','b');
ylim([WV_min WV_max]); 

cd(dd) % point back to original directory