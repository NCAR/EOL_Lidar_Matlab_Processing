clear all; close all;

%dd = pwd; % get the current path
date = '29 Apr 2019'; % Last day of combined data  
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

MPD01.location = [36.311 -97.927 386.0]; %lat, lon, elevation
MPD02.location = [36.880 -97.086 340.0]; %lat, lon, elevation
MPD03.location = [36.820 -97.820 337.0]; %lat, lon, elevation
MPD04.location = [36.373 -97.069 311.0]; %lat, lon, elevation
MPD05.location = [36.605 -97.486 317.0]; %lat, lon, elevation


MPD01.loc_array = repmat(MPD01.location, length(MPD01.range),1)';
MPD02.loc_array = repmat(MPD02.location, length(MPD02.range),1)';
MPD03.loc_array = repmat(MPD03.location, length(MPD03.range),1)';
MPD04.loc_array = repmat(MPD04.location, length(MPD04.range),1)';
MPD05.loc_array = repmat(MPD05.location, length(MPD05.range),1)';

MPD01.AH = MPD01.N_avg_comb.*1e6./6.022E23.*18.015;
MPD02.AH = MPD02.N_avg_comb.*1e6./6.022E23.*18.015;
MPD03.AH = MPD03.N_avg_comb.*1e6./6.022E23.*18.015;
MPD04.AH = MPD04.N_avg_comb.*1e6./6.022E23.*18.015;
MPD05.AH = MPD05.N_avg_comb.*1e6./6.022E23.*18.015;

figure(101)
font_size = 14;
S = 75;
ii = 489;
for ii = 1:size(MPD01.AH,1)
  date=datestr(MPD01.time(ii), 'dd mmm yy HH:MM');
  scatter3(MPD01.loc_array(2,:),MPD01.loc_array(1,:),(MPD01.loc_array(3,:)+MPD01.range)./1000, S, MPD01.AH(ii,:), 'filled') 
  xlabel('longitude')
  ylabel('latitude')
  zlabel('elevation, MSL (km)')
  title({'MPD absolute humidity (gm^{-3})  ', date})  
  set(gca,'Fontsize',font_size,'Fontweight','b')
  %set(gca, 'color', 'w')
  zlim([0 6])
  caxis([0, 16])
  view(-30, 25) %azimuth, elevation
  hold on
  scatter3(MPD02.loc_array(2,:),MPD02.loc_array(1,:),(MPD02.loc_array(3,:)+MPD02.range)./1000, S, MPD02.AH(ii,:), 'filled') 
  scatter3(MPD03.loc_array(2,:),MPD03.loc_array(1,:),(MPD03.loc_array(3,:)+MPD03.range)./1000, S, MPD03.AH(ii,:), 'filled') 
  scatter3(MPD04.loc_array(2,:),MPD04.loc_array(1,:),(MPD04.loc_array(3,:)+MPD04.range)./1000, S, MPD04.AH(ii,:), 'filled') 
  scatter3(MPD05.loc_array(2,:),MPD05.loc_array(1,:),(MPD05.loc_array(3,:)+MPD05.range)./1000, S, MPD05.AH(ii,:), 'filled') 
  hold off
  colorbar;
  colormap(jet)
  pause(0.01)
end

figure(102)
x=[MPD01.loc_array(2,:); MPD02.loc_array(2,:); MPD03.loc_array(2,:); MPD04.loc_array(2,:); MPD05.loc_array(2,:)];
y=[MPD01.loc_array(1,:); MPD02.loc_array(1,:); MPD03.loc_array(1,:); MPD04.loc_array(1,:); MPD05.loc_array(1,:)];
z=[MPD01.loc_array(3,:); MPD02.loc_array(3,:); MPD03.loc_array(3,:); MPD04.loc_array(3,:); MPD05.loc_array(3,:)];
v=double([MPD01.AH(ii,:); MPD02.AH(ii,:); MPD03.AH(ii,:); MPD04.AH(ii,:); MPD05.AH(ii,:)]);
[xi,yi,zi] = meshgrid(36.3:0.05:36.9, -98:0.1:-97, 300:50:6000);
vq = griddata(x,y,z,v,xi,yi,zi);
surf(xi,yi,zi,vq);
scatter3(xi,yi,zi,S,vq)
stop

x = [0,45.74,83.54]; % distance in km from MPD04 (at SGP39) to 05(CF) to 03 at(SGP E32 
y = MPD05.range./1e3;
scrsz = get(0,'ScreenSize');
Scrnsize=[scrsz(4)/2 scrsz(4)/10 scrsz(3)/1 scrsz(4)/2];

for i = 1:24
  time{i} = datenum(strcat('28 Apr 2019 0',num2str(5+i),':00')); % good
  time_idx{i} = find(MPD05.time == time{i});
  Z{i} = [MPD04.N_avg_comb(time_idx{i},:)', MPD05.N_avg_comb(time_idx{i},:)', MPD03.N_avg_comb(time_idx{i},:)'];
  AH{i} = double(real(Z{i}.*1e6./6.022E23.*18.015)); %convert to AH

 % plot Narrow water vapor in g/m^3
 figure('Position',Scrnsize)
 set(gcf,'renderer','zbuffer');
 h = pcolor(x,y,AH{i});
 set(h, 'EdgeColor', 'none');
 colorbar('EastOutside');
 axis([fix(min(x)) ceil(max(x)) 0 6])
 caxis([0 18]);
 colormap(jet)
 shading interp
 xlabel('Distance (km, along transect from MPD04)');
 ylabel('Height (km, AGL)');
 title(datestr(time{i}))
end



%