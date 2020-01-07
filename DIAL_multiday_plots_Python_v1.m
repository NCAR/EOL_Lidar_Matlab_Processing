clear all; close all;
tic

%DIAL=1;

%DIAL=2;

%DIAL=3;

%DIAL=4;

DIAL=5;
date = '18 Apr 2019'; % Five-unit side-by-side test 
%days = 4; skip = 1;
%days = 32; skip = 4;
days = 72; skip = 8;
dirname = '/Volumes/documents/WV_DIAL_data/MPD5_python_10min_150m_processed_data';
  
font_size = 36; % use this for 2018a version
%font_size = 14; % use this for 2015a version

% set to 1 for using the weather station data (after Jan 2016)
WS = 0;
%include sonde data 0=off 1=on
flag.sonde = 0;
%replot time vs range images at start of processing 0=off 1=on
flag.replot = 1;
%save figures at end of processing 0=off 1=on
flag.save_figs = 0;
%Wide field channel 0=off 1=on
flag.save_data = 0;
%plot NetCDF sonde data on top 0=off 1=on
flag.plot_sonde_data = 0;

near_field = 0;  % now the HSRL channel

%number of accumulations 
if datenum(date)>=datenum('16 Nov 2018')&& datenum(date)<=datenum('18 Dec 2018')
  accum = 32000; %  
  bin_duration = 500;  % ns 
  RB_scale = 1; % use to keep the arbitrary units of RB scale the same before
else
  accum = 14000; %this was changed before perdigao need to get the date
  bin_duration = 250;  % ns (this change from 50 to 500 on 14-June-2014)
  RB_scale = 1;
end
  
C = importdata('NCAR_C_Map.mat');
%cd('/Users/spuler/Desktop/FRAPPE_PECAN') % point to the directory where data is stored 
%cd('/Users/spuler/Desktop/WV_DIAL_data/') % point to the directory where data is stored 

if DIAL==1
  cd('/Volumes/documents/WV_DIAL_data/MPD1_processed_data') % point to the directory where data is stored 
  %bin duration in ns
elseif DIAL==2
  cd('/Volumes/documents/WV_DIAL_data/MPD2_processed_data') % point to the directory where data is stored 
  bin_duration = 250;  % ns (this change from 500 to 250 for DIAL#2 on 2014)
  near_field = 1;  % now the HSRL channel
elseif DIAL==3
  cd('/Volumes/documents/WV_DIAL_data/MPD3_processed_data') % point to the directory where data is stored 
elseif DIAL==4
  cd('/Volumes/documents/WV_DIAL_data/MPD4_processed_data') % point to the directory where data is stored 
elseif DIAL==5
  cd(dirname) % point to the directory where data is stored 
elseif DIAL==5

end


gate = round((bin_duration*1e-9*3e8/2)*10)/10
i=1;

for i=1:days
  if i==1  
    if exist(strcat(date, '.mat'))==2
      load(strcat(date, '.mat'))
    end
    range_limit = size(N_avg,2);
    RB = RB(1:size(N_avg,1), 1:size(N_avg,2)); %add this line to force RB to same size as N_avg
    % grid everything to a 75 m gate size 
     if round(diff(range(1:2)),1) == 37.5
       range_grid_75 = 0:75:(range_limit-1)*37.5; 
       N_avg = interp1(range, N_avg, range_grid_75, 'linear', 'extrap')'; 
       RB = interp1(range, RB, range_grid_75, 'linear', 'extrap')';
       range = range_grid_75;
       range_limit = range_limit/2;
     end
    N_avg_comb=N_avg;
    RB_comb=RB;
    duration=time_new;
    % if gate = 37.5 down sample to 75 m bins
  else
    date = datestr(addtodate(datenum(date), 1, 'day'), 'dd mmm yyyy');
    if exist(strcat(date, '.mat'))==2
      load(strcat(date, '.mat'))
    end
    range_limit_ch = min(size(N_avg,2), size(N_avg_comb,2)*2);
    RB = RB(1:size(N_avg,1), 1:size(N_avg,2)); %add this line to force RB to same size as N_avg
    % grid everything to a 75 m gate size 
      if round(diff(range(1:2)),1) == 37.5
         range_grid_75 = 0:75:(range_limit_ch-1)*37.5; 
         N_avg = interp1(range, N_avg, range_grid_75, 'linear', 'extrap')'; 
         RB = interp1(range, RB, range_grid_75, 'linear', 'extrap')';
         range = range_grid_75;
         range_limit = range_limit_ch/2;
      end
    range_limit = size(N_avg,2); % catch any changes in range
    N_avg_comb = vertcat(N_avg_comb(:,1:range_limit), N_avg(2:end,1:range_limit));
    RB_comb = vertcat(RB_comb(:,1:range_limit), RB(2:end,1:range_limit));
    duration = vertcat(duration, time_new(2:end));
  end
end

%stop
 %% save data
  
 if flag.save_data == 1
  %cd('/Users/spuler/Desktop/WV_DIAL_data') % point to the directory where data is stored 
  name=strcat(date, '_combined');
  if DIAL==1
      MPD01.N_avg_comb = N_avg_comb;
      MPD01.RB_comb = RB_comb;
      MPD01.range = range;
      MPD01.time = duration;
      save(name, 'MPD01')
  elseif DIAL==2
      MPD02.N_avg_comb = N_avg_comb;
      MPD02.RB_comb = RB_comb;
      MPD02.range = range;
      MPD02.time = duration;
      save(name, 'MPD02')
  elseif DIAL==3
      MPD03.N_avg_comb = N_avg_comb;
      MPD03.RB_comb = RB_comb;
      MPD03.range = range;
      MPD03.time = duration;
      save(name, 'MPD03')
  elseif DIAL==4
      MPD04.N_avg_comb = N_avg_comb;
      MPD04.RB_comb = RB_comb;
      MPD04.range = range;
      MPD04.time = duration;
      save(name, 'MPD04')
   elseif DIAL==5
      MPD05.N_avg_comb = N_avg_comb;
      MPD05.RB_comb = RB_comb;
      MPD05.range = range;
      MPD05.time = duration;
      save(name, 'MPD05')
  end

 end


scrsz = get(0,'ScreenSize');
Scrnsize=[scrsz(4)/2 scrsz(4)/10 scrsz(3)/1 scrsz(4)/2];

if days == 1
  xData =  linspace(fix(min(duration)),  ceil(max(duration)), 25);
else
  xData =  linspace( fix(min(duration)),  ceil(max(duration)), round((ceil(max(duration))-fix(min(duration)))/skip)+1 );
 % xData =  linspace(fix(min(duration)),  round(max(duration)), 36);
end
x = (duration)';
y = (range(1:range_limit)./1e3);

if flag.replot==1
 % plot Narrow water vapor in g/m^3
 figure('Position',Scrnsize)
 Z = double(real(N_avg_comb'.*1e6./6.022E23.*18.015));  %number density in mol/cm3(1e6 cm3/m3)/(N_A mol/mole)*(18g/mole)
 % Z(isnan(Z)) = -1;
 set(gcf,'renderer','zbuffer');
 h = pcolor(x,y,Z);
 set(h, 'EdgeColor', 'none');
 colorbar('EastOutside');
 axis([fix(min(x)) ceil(max(x)) 0 6])
 caxis([0 20]);
 colormap(jet)
 %colormap(C)
 %colormap(perula)
 %shading interp
 % P_t = get(hh, 'Position');
 % set(hh,'Position', [P_t(1) P_t(2)+0.1 P_t(3)])
 ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
 set(gca, 'XTick',  xData)
 set(gca,'TickDir','out');
 set(gca,'TickLength',[0.005; 0.0025]);
 if days == 1
   datetick('x','HH:MM','keeplimits', 'keepticks');
   xlabel('Time (UTC)','fontweight','b','fontsize',font_size);
   hh = title({[date,'DIAL Water Vapor (g m^{-3})']},'fontweight','b','fontsize',font_size);
 else
   datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
   hh = title({'DIAL Water Vapor (g m^{-3})'},'fontweight','b','fontsize',font_size);
 end
 set(gca,'Fontsize',font_size,'Fontweight','b');
 
 
 % plot Narrow RB
 figure('Position',Scrnsize)
 Z = double((real(RB_comb')./RB_scale));
% Z_mask = Z;
% Z_mask(RB_FF'<5) = NaN;
 % Z(isnan(Z)) = -1;
 set(gcf,'renderer','zbuffer');
 h = pcolor(x,y,Z);
%  h = pcolor(x,y,Z_mask);
 set(h, 'EdgeColor', 'none');
 colorbar('EastOutside');
 axis([fix(min(duration)) ceil(max(duration)) 0 12])
 caxis([1e3 1e9]);
 colormap(C)
 %colormap(perula)
 %shading interp
 % P_t = get(hh, 'Position');
 % set(hh,'Position', [P_t(1) P_t(2)+0.2 P_t(3)])
 ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size);
 set(gca, 'XTick',  xData)
 set(gca,'TickDir','out');
 set(gca,'TickLength',[0.005; 0.0025]);
 if days == 1
   datetick('x','HH:MM','keeplimits', 'keepticks');
   xlabel('Time (UTC)','fontweight','b','fontsize',font_size);
  % hh = title({[date,'DIAL Relative Backscatter (C/ns km^2)']},'fontweight','b','fontsize',font_size);
   hh = title({[date,'DIAL Attenuated Backscatter (A.U.)']},'fontweight','b','fontsize',font_size);
 else
   datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
 %  hh = title({'DIAL Relative Backscatter (C/ns km^2)'},'fontweight','b','fontsize',font_size);
    hh = title({'DIAL Attenuated Backscatter (A.U.)'},'fontweight','b','fontsize',font_size);
 end
 set(gca,'Fontsize',font_size,'Fontweight','b');
 set(gca,'Zscale', 'log')
 set(gca,'Colorscale', 'log')
 set(gca,'Zscale', 'linear')

 %end
 
%%
if flag.plot_sonde_data==1
  % the elevation needs to be read from the json file
    elevation= 310.0; %MPD05 was at 310m elevation at SGP
  %d=pwd;
  %cd('/Volumes/documents/WV_DIAL_data/SGP_sondes/') % point to the directory where data is stored
   cd('/scr/sci/tammy/mpd/sgp/soundings/')
   [sondefilename, sondedir] = uigetfile('*.*','Select the sonde file', 'MultiSelect', 'on');
   jj=1;
  %cd= d;
  for jj = 1:size(sondefilename,2)
      cd('/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_processing/')
      Sonde_read_nc_files(jj, elevation, sondedir, sondefilename);
  end
end 

%% 
if flag.save_figs==1
  %cd('/Users/spuler/Desktop/WV_DIAL_data/plots/') % point to the directory where data is stored 
  cd('/Volumes/documents/WV_DIAL_data/plots/') % point to the directory where data is stored 
  date=datestr(nanmean(time_new), 'yyyymmdd');
  
  %Scrnsize = [scrsz(4)/2 scrsz(4)/10 scrsz(3)/1 scrsz(4)/2]; % use for standard plots
  Scrnsize = [scrsz(4)/1 scrsz(4)/1 scrsz(3)/0.35 scrsz(4)/2.05]; % use for long plots 
  Scrnsize = [scrsz(4)/1 scrsz(4)/1 scrsz(3)/0.30 scrsz(4)/2]; % use for ILRC really long plots 
  %Scrnsize = [scrsz(4)/1 scrsz(4)/1 scrsz(3)/0.51 scrsz(4)/2]; % use for Perdigao BAMS plots 
  %Scrnsize = [scrsz(4)/1 scrsz(4)/1 scrsz(3)/2 scrsz(4)/2]; % use for day plots 
  %Scrnsize = [scrsz(4)/1 scrsz(4)/1 scrsz(3)/1 scrsz(4)/2.2]; % use for AMT sized 3-day plots (with large font)
  
  FigH = figure(1);
  set(gca,'Fontsize',36,'Fontweight','b'); % use for Perdigao BAMS plots 
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', Scrnsize);
  name=strcat(date, 'H2O_multi'); 
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
 
  FigH = figure(2);
  set(gca,'Fontsize',36,'Fontweight','b'); % use for Perdigao BAMS plots 
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', Scrnsize);
  name=strcat(date, 'RB_multi'); 
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
  
 % FigH = figure(3);
 % % set(gca,'Fontsize',36,'Fontweight','b');
 % set(FigH, 'PaperUnits', 'points', 'PaperPosition', Scrnsize);
 % name=strcat(date, 'background_multi'); 
 % print(FigH, name, '-dpng', '-r0') % set at the screen resolution
  
%  FigH = figure(5);
% %  set(gca,'Fontsize',36,'Fontweight','b');
%  set(FigH, 'PaperUnits', 'points', 'PaperPosition', Scrnsize);
%  name=strcat(date, 'column_OD_multi'); 
%  print(FigH, name, '-dpng', '-r0') % set at the screen resolution
  
  if WS==1
      size2 = [scrsz(4)/1 scrsz(4)/1 scrsz(3)/0.35 scrsz(4)/1]; % use for long plots 
      FigH = figure(6);
     % set(gca,'Fontsize',36,'Fontweight','b');
      set(FigH, 'PaperUnits', 'points', 'PaperPosition', size2);
      name=strcat(date, 'Housekeeping'); 
      print(FigH, name, '-dpng', '-r0') % set at the screen resolution
  end
   
end
 
  name=strcat(date, '_WV_2CH'); 
 
end


 

 
toc
