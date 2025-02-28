clear all; close all;
tic

 node = 'MPD05';            
 date = '03 Dec 2024';
 days = 7; skip = 1;
%      date = '18 Jul 2023';
%      days = 70; skip = 5;
 WV_max_scale = 6;
 flag.afterpulse = 0; % read in the afterpulse corrected data (0=off 1=on)
 
%  node = 'MPD02';            %post PRECIP intercomparions 
%  date = '3 Jun 2022';   
%  days = 2; skip = 1;
%  WV_max_scale = 25;
%  flag.afterpulse = 1; % read in the afterpulse corrected data (0=off 1=on)
 

%serv_path = '/Volumes/documents/MPD/';
%serv_path = '/Volumes/eol/fog1/rsfdata/MPD/';
serv_path = '/Volumes/smaug1/rsfdata/MPD/';
plot_path = '/Volumes/Macintosh HD/Users/spuler/Desktop/mpd/Plots/';
addpath '/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing/';
addpath '/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing/matplotlib/';

blank = 150; % has to be in increments of 75 (the default blank in the processing is 300m before WFOV)
flag.near = 0;  % read in the near range channel (0=off 1=on)

font_size = 36; % use this for 2018a version

WS = 1;  % set to 1 for using the weather station data (after Jan 2016)
flag.sonde = 0;  %include sonde data (0=off 1=on)
flag.replot = 1;  %replot time vs range images at start of processing (0=off 1=on)
flag.save_figs = 1; %save figures at end of processing (0=off 1=on)
flag.save_data = 0; %save data at end of processing (0=off 1=on)
flag.plot_sonde_data = 0; %plot NetCDF sonde data on top (0=off 1=on)
flag.decimate = 0; %decimate figures to the screen 2x size (1900 pixels x2)
if days<=2
    flag.decimate = 0;
end
range_grid_size = 75;  %set the size of the range gridding


RB_scale = 1;

  
C = importdata('NCAR_C_Map.mat');
%cd('/Users/spuler/Desktop/FRAPPE_PECAN') % point to the directory where data is stored 
%cd('/Users/spuler/Desktop/WV_DIAL_data/') % point to the directory where data is stored 
dd=pwd;

if strcmp(node,'MPD01')==1
    if flag.afterpulse == 1
      cd(strcat(serv_path,'/mpd_01_processed_data/Matlab/afterpulse'))    
    else
      cd(strcat(serv_path,'/mpd_01_processed_data/Matlab'))
    end
 elseif strcmp(node,'MPD02')==1
    if flag.afterpulse == 1
      cd(strcat(serv_path,'/mpd_02_processed_data/Matlab/afterpulse'))    
    else
      cd(strcat(serv_path,'/mpd_02_processed_data/Matlab'))
    end
 elseif strcmp(node,'MPD03')==1
    if flag.afterpulse == 1
      cd(strcat(serv_path,'/mpd_03_processed_data/Matlab/afterpulse'))    
    else
      cd(strcat(serv_path,'/mpd_03_processed_data/Matlab'))
    end
 elseif strcmp(node,'MPD04')==1
    if flag.afterpulse == 1
      cd(strcat(serv_path,'/mpd_04_processed_data/Matlab/afterpulse'))    
    else
      cd(strcat(serv_path,'/mpd_04_processed_data/Matlab'))
    end
 elseif strcmp(node,'MPD05')==1
    if flag.afterpulse == 1
      cd(strcat(serv_path,'/mpd_05_processed_data/Matlab/afterpulse'))    
    else
      cd(strcat(serv_path,'/mpd_05_processed_data/Matlab'))
    end
end

%gate = round((bin_duration*1e-9*3e8/2)*10)/10

test_gate = exist('gate') % check for early versions were the gate wasn't saved
if test_gate == 0
 %gate = 75
end

i=1;


% [Matlabfilename, Matlabdir] = uigetfile('*.*','Select the sonde file', 'MultiSelect', 'on');
% jj=1;


for i=1:days
  if i==1  
     if exist(strcat(node, '_', datestr(date, 'yyyymmdd'), '_WV.mat'))==2
        load(strcat(node, '_', datestr(date, 'yyyymmdd'), '_WV.mat'))
        range_limit = size(N_avg,2);
       % grid everything to a 75 m gate size 
        if gate < 75
           range_grid_75 = 0:range_grid_size:(range_limit-1)*gate; 
           N_avg = interp1(range, N_avg', range_grid_75, 'linear', 'extrap')'; 
           N_error = interp1(range, N_error', range_grid_75, 'linear', 'extrap')'; 
           RB = interp1(range, RB', range_grid_75, 'linear', 'extrap')';
           OD = interp1(range, OD', range_grid_75, 'linear', 'extrap')';
           range = range_grid_75;
           range_limit = range_limit/(75/gate);
        end
        N_avg_comb=N_avg;
        N_error_comb=N_error;
        RB_comb=RB;
        background_comb_on = background_on;
        background_comb_off = background_off;
        lambda_comb_on = lambda_all;
        lambda_comb_off = lambda_all_off;
        duration=time_new;
        % if gate = 37.5 down sample to 75 m bins
       if WS==1
         surf_T = Surf_T;
         surf_P = Surf_P;
         surf_AH = Surf_AH;
         i_off = I_off;
         i_on = I_on;
         p_on = P_on;
         p_off = P_off;
         t_bench = T_bench;
       end
     end  
  else
    date = datestr(addtodate(datenum(date), 1, 'day'), 'dd mmm yyyy');
    if exist(strcat(node, '_', datestr(date, 'yyyymmdd'), '_WV.mat'))==2
       load(strcat(node, '_', datestr(date, 'yyyymmdd'), '_WV.mat'))
       range_limit_ch = size(N_avg,2);
       % grid everything to a 75 m gate size 
       if gate < 75
         range_grid_75 = 0:range_grid_size:(range_limit_ch-1)*gate; 
         N_avg = interp1(range, N_avg', range_grid_75, 'linear', 'extrap')'; 
         RB = interp1(range, RB', range_grid_75, 'linear', 'extrap')';
  %      OD = interp1(range, OD', range_grid_75, 'linear', 'extrap')';
         range = range_grid_75;
         range_limit = range_limit_ch/(75/gate);
       end
      range_lim1 = size(N_avg_comb,2); % catch any changes in range
      range_lim2 = size(N_avg,2); % catch any changes in range
      range_limit = min([range_lim1 range_lim2]);
    
      N_avg_comb = vertcat(N_avg_comb(:,1:range_limit), N_avg(2:end,1:range_limit));
      N_error_comb = vertcat(N_error_comb(:,1:range_limit), N_error(2:end,1:range_limit));
      RB_comb = vertcat(RB_comb(:,1:range_limit), RB(2:end,1:range_limit));
 %    OD_comb= vertcat(OD_comb(:,1:range_limit), OD(2:end,1:range_limit));
      duration = vertcat(duration, time_new(2:end));
      background_comb_on = vertcat(background_comb_on, background_on(2:end));
      background_comb_off = vertcat(background_comb_off, background_off(2:end));
      lambda_comb_on = vertcat(lambda_comb_on, lambda_all(2:end));
      lambda_comb_off = vertcat(lambda_comb_off, lambda_all_off(2:end));
      if WS==1
        surf_T = vertcat(surf_T,Surf_T(2:end,:));
        surf_P = vertcat(surf_P,Surf_P(2:end,:));
        surf_AH = vertcat(surf_AH, Surf_AH(2:end,:));
        i_off = vertcat(i_off,I_off(2:end,:));
        i_on = vertcat(i_on, I_on(2:end,:));
        p_on = vertcat(p_on, P_on(2:end,:));  
        p_off = vertcat(p_off, P_off(2:end,:)); 
        t_bench = vertcat(t_bench, T_bench(2:end,:));  
      end
    end
  end
end

%stop
T_lapse = surf_T-0.0065*range; %create an array of the T based on 6.5�C/km lapse rate

% remove the lowest bins 
N_avg_comb(:,2:blank/75) = NaN; 

% corrected/scaled lowest bins 
% N_avg_comb(:,9) = N_avg_comb(:,9).*0.97; % 3% at 675 m 
% N_avg_comb(:,8) = N_avg_comb(:,8).*0.90; % 10% at 600 m 
% N_avg_comb(:,7) = N_avg_comb(:,7).*0.85;  % 15% at 525 m
% N_avg_comb(:,6) = N_avg_comb(:,6).*0.85;  % 15% at 450 m
% N_avg_comb(:,5) = N_avg_comb(:,5).*0.88;  % 12% at 375 m
% N_avg_comb(:,1:4) = N_avg_comb(:,1:4).*0.92;  % 8% < 300 m


 %% save data
  
 %% save data
  
 if flag.save_data == 1
  %cd('/Users/spuler/Desktop/WV_DIAL_data') % point to the directory where data is stored 
  %name=strcat(date, '_combined');
  if strcmp(node, 'MPD01') ==1
      MPD01.N_avg_comb = N_avg_comb;
      MPD01.RB_comb = RB_comb;
      MPD01.range = range;
      MPD01.time = duration;
      cd(strcat(serv_path,'/mpd_01_processed_data/Matlab'))
      name = strcat('MPD01_',date,'_Matlab_combined');
      save(name, 'MPD01')
  elseif strcmp(node, 'MPD02') ==1
      MPD02.N_avg_comb = N_avg_comb;
      MPD02.RB_comb = RB_comb;
      MPD02.range = range;
      MPD02.time = duration;
      cd(strcat(serv_path,'/mpd_02_processed_data/Matlab'))
      name = strcat('MPD02_',date,'_Matlab_combined');
      save(name, 'MPD02')
  elseif strcmp(node, 'MPD03') ==1
      MPD03.N_avg_comb = N_avg_comb;
      MPD03.RB_comb = RB_comb;
      MPD03.range = range;
      MPD03.time = duration;
      cd(strcat(serv_path,'/mpd_03_processed_data/Matlab'))
      name = strcat('MPD03_',date,'_Matlab_combined');
      save(name, 'MPD03')
  elseif strcmp(node, 'MPD04') ==1
      MPD04.N_avg_comb = N_avg_comb;
      MPD04.RB_comb = RB_comb;
      MPD04.range = range;
      MPD04.time = duration;
      cd(strcat(serv_path,'/mpd_04_processed_data/Matlab'))
      name = strcat('MPD04_',date,'_Matlab_combined');
      save(name, 'MPD04')
   elseif strcmp(node, 'MPD05') ==1
      MPD05.N_avg_comb = N_avg_comb;
      MPD05.RB_comb = RB_comb;
      MPD05.range = range;
      MPD05.time = duration;
      cd(strcat(serv_path,'/mpd_05_processed_data/Matlab'))
      name = strcat('MPD05_',date,'_Matlab_combined');
      save(name, 'MPD05')
  end
 
 end

 if days == 1
  xData =  linspace(fix(min(duration)),  ceil(max(duration)), 25);
else
  xData =  linspace( fix(min(duration)),  ceil(max(duration)), round((ceil(max(duration))-fix(min(duration)))/skip)+1 );
 % xData =  linspace(fix(min(duration)),  round(max(duration)), 36);
 end

% point back to original directory 
cd(dd);
x = (duration)';
y = (range(1:range_limit)./1e3);
Z_AH = double(real(N_avg_comb'.*1e6./6.022E23.*18.015));  %number density in mol/cm3(1e6 cm3/m3)/(N_A mol/mole)*(18g/mole)
Z_RB = double((real(RB_comb')./RB_scale));
 
%scrsz = get(0,'ScreenSize');
scrsz = [1  1  1920 1200];
% decimate the duraion down to the screen resolution (should be 1900 pixels)
flag.int = 0; % interpolate nans in nanmoving_average
if flag.decimate == 1 
    decimate_time = fix(length(duration)/scrsz(3)/2); %decimate to number of horiz pixels;
    decimate_range = 1; % keep native gate spacing 
    % average RB data before decimating
     Z_AH = nanmoving_average(Z_AH,decimate_time,2,flag.int);
     Z_RB = nanmoving_average(Z_RB,decimate_time,2,flag.int);
    % then decimate
     x = x(1:decimate_time:end);
     Z_AH =  Z_AH(1:decimate_range:end, 1:decimate_time:end);
     Z_RB =  Z_RB(1:decimate_range:end, 1:decimate_time:end);   
end


Scrnsize=[scrsz(4)/2 scrsz(4)/10 scrsz(3)/1 scrsz(4)/2];
if flag.replot==1
 % plot Narrow water vapor in g/m^3
 figure('Position',Scrnsize)
 %Z = double(real(N_avg_comb'.*1e6./6.022E23.*18.015));  %number density in mol/cm3(1e6 cm3/m3)/(N_A mol/mole)*(18g/mole)
 % Z(isnan(Z)) = -1;
 set(gcf,'renderer','zbuffer');
 h = pcolor(x,y,Z_AH);
 set(h, 'EdgeColor', 'none');
 colorbar('EastOutside');
 axis([fix(min(x)) ceil(max(x)) 0 6])
 caxis([0 WV_max_scale]);
 %caxis([0 5]);
 colormap(jet)
 colormap(CM_YlGnBu(64))
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
   hh = title({[node, ' Water Vapor (g m^{-3})', '  ', date]},'fontweight','b','fontsize',font_size);
  % hh = title({[date,'DIAL Water Vapor (g m^{-3})']},'fontweight','b','fontsize',font_size);
 else
   datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
   hh = title({[node, ' Water Vapor (g m^{-3})']},'fontweight','b','fontsize',font_size);
   %hh = title({'DIAL Water Vapor (g m^{-3})'},'fontweight','b','fontsize',font_size);
 end
 set(gca,'Fontsize',font_size,'Fontweight','b');
 grid on
 
 
 % plot Narrow RB
 figure('Position',Scrnsize)
% Z = double((real(RB_comb')./RB_scale));
% Z_mask = Z;
% Z_mask(RB_FF'<5) = NaN;
 % Z(isnan(Z)) = -1;
 set(gcf,'renderer','zbuffer');
 h = pcolor(x,y,Z_RB);
%  h = pcolor(x,y,Z_mask);
 set(h, 'EdgeColor', 'none');
 colorbar('EastOutside');
 %axis([fix(min(duration)) ceil(max(duration)) 0 18.25])
 axis([fix(min(duration)) ceil(max(duration)) 0 12])
 caxis([1e1 1e6]);
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

 
% plot the background
  figure('Position',Scrnsize)
 % if near_field==1 
 %   len1= length(duration);
 %   len2 = length(background_comb_on);
 %   len3 = length(background_comb_off); % HSRL molecular channel
 %
 %   len = vertcat(len1, len2, len3);
 %   plot_end = min(len);  
 %   semilogy(duration(1:plot_end), background_comb_on(1:plot_end), 'r', duration(1:plot_end), background_comb_off(1:plot_end), 'k') 
 % else
    semilogy(duration, background_comb_off, 'black')
 % end
  axis([fix(min(duration)) ceil(max(duration)) 1e2 1e7])
  ylabel('Offline background C/s', 'Fontsize', font_size, 'Fontweight', 'b');  
  grid on
  %hold on
  %  semilogy(duration, 5e6, 'green', 'Linewidth', 2)    % Add a horizontal line to show the 5Mc/s linearity limit
  %hold off
 % if near_field==1  
 % legend('WV Online  ', 'WV Offline  ', 'HSRL molecular  ', 'HSRL combined  ','Location', 'NorthWest');
 % else
 legend('WV Online  ', 'WV Offline  ', 'Location', 'NorthWest');
 % end
  set(gca, 'XTick',  xData)
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  if days == 1
    datetick('x','HH:MM','keeplimits', 'keepticks');
    xlabel('Time (UTC)','fontweight','b','fontsize',font_size); 
    hh = title({[date,'Background']},'fontweight','b','fontsize',font_size);
  else
    datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
    hh = title({'DIAL Background'},'fontweight','b','fontsize',font_size); 
  end
  set(gca,'Fontsize',font_size,'Fontweight','b');
 
  
  
% plot column OD for the Narrow  
  % figure('Position',Scrnsize)
  % Z = double(real(OD_comb'));
  % set(gcf,'renderer','zbuffer');
  % h = pcolor(x,y,Z);
  % set(h, 'EdgeColor', 'none');
  % colorbar('EastOutside');
  % axis([fix(min(duration)) ceil(max(duration)) 0 12])
  % caxis([-0.1 2]);
  % colormap(C)
  %title({('  Column Optical Depth, Narrow Field  ')},...
  %     'fontweight','b','fontsize',font_size)
  % xlabel('Time (UTC)','fontweight','b','fontsize',font_size); 
  % ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
  % set(gca, 'XTick',  xData)
  % set(gca,'TickDir','out');
  % set(gca,'TickLength',[0.005; 0.0025]);
  % if days == 1
  %   datetick('x','HH:MM','keeplimits', 'keepticks');
  %   xlabel('Time (UTC)','fontweight','b','fontsize',font_size); 
  % else
  %   datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
  % end
  % set(gca,'Fontsize',font_size,'Fontweight','b');
%
%  
%  plot1 = OD_comb(:,round(2500/gate));
%  plot2 = OD_comb(:,round(5000/gate));
%  %window = ones(1,1)/1;
%  %plot1s = filter2(window, plot1);
%  %plot2s = filter2(window, plot2);
%  dd = pwd; % get the current path
%  cd('/Users/spuler/Desktop/WV_DIAL/Matlab/')
%    plot1s = ndnanfilter(plot1,'rectwin', 20); 
%    plot2s = ndnanfilter(plot2,'rectwin', 20);
%  cd(dd)
  
% % plot column optical depth at 5km and 2.5km as function of time
%  figure('Position',Scrnsize)
%  plot(duration, plot1s, 'b',duration, plot2s, 'k', 'LineWidth', 2)
%  axis([fix(min(duration)) ceil(max(duration)) 0 2.5])
%  ylabel('OD', 'Fontsize', font_size, 'Fontweight', 'b');
%  grid on
%  legend('OD at 2.5km', 'OD at 5.0km','Location', 'NorthWest');
%  set(gca, 'XTick',  xData)
%  set(gca,'TickDir','out');
%  set(gca,'TickLength',[0.005; 0.0025]);
%  if days == 1
%    datetick('x','HH:MM','keeplimits', 'keepticks');
%    xlabel('Time (UTC)','fontweight','b','fontsize',font_size); 
%    hh = title({[date,'DIAL Column OD']},'fontweight','b','fontsize',font_size);
%  else
%    datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
%       hh = title({'DIAL Column OD'},'fontweight','b','fontsize',font_size); 
%  end
%  set(gca,'Fontsize',font_size,'Fontweight','b');
  
 if WS==1
     font_size = 14
 % plot housekeeping data
   figure1 = figure('Position',Scrnsize);
   subplot1=subplot(2,1,1,'Parent',figure1,'YGrid','on', 'XGrid','on');
   box(subplot1,'on');
   hold(subplot1,'all');
  %plot(duration, (i_off),'b','LineWidth',2,'DisplayName','i_{off}') % these plot diode Temps
  %plot(duration, (i_on),'r','LineWidth',2, 'DisplayName','i_{on}')
  plot(duration, (lambda_comb_on),'k','LineWidth',2,'DisplayName','Lambda_{on}') % these plot diode Temps
  %plot(duration, (t_hsrl),'g','LineWidth',2, 'DisplayName','T_{hsrl}')
   axis([fix(min(duration)) ceil(max(duration)) 828.190 828.200])
   YTick = [100 120 140 160 180];
   ylabel('wavelength, nm', 'Fontsize', font_size, 'Fontweight', 'b');  
   datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
   set(gca,'Fontsize',font_size,'Fontweight','b');
   % Plot the temperature data
   subplot2=subplot(2,1,2,'Parent',figure1,'YGrid','on', 'XGrid','on');
   box(subplot2,'on');
   hold(subplot2,'all');
   plot(duration, t_bench,'r', 'LineWidth',1, 'DisplayName','Bench T')
   plot(duration, surf_T, 'b', 'LineWidth',1, 'DisplayName','Surface T')
 %  axis([fix(min(duration)) ceil(max(duration)) -inf inf]);   % -20 40])
  axis([fix(min(duration)) ceil(max(duration)) 0 40]);   % -20 40]) %PRECIP
      YTick = [-25 0 25 50];
   ylabel('Temp, C', 'Fontsize', font_size, 'Fontweight', 'b'); 
   datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
   set(gca,'Fontsize',font_size,'Fontweight','b');
   % Create legend
   legend(subplot1,'show','Location','southwest'); 
   legend(subplot2,'show','Location','southwest');
   %link the x axis for all 3 subplots
   ax(1)=subplot(2,1,1);
   ax(2)=subplot(2,1,2);
   % plot power on right y-axis of the upper plot (% assumes 5% pickoff)
   ax(3) = axes('Position',get(ax(1),'Position'));
   plot(duration, (p_off/7500),'b-','LineWidth', 1, 'DisplayName','P_{off}') % changed from 0.05 to 0.0425
   hold on
   plot(duration, p_on/7500,'r-', 'LineWidth',1, 'DisplayName','P_{on}')
   %plot(duration, p_hsrl/2500,'g--', 'LineWidth',1, 'DisplayName','P_{hsrl}')
   axis([fix(min(duration)) ceil(max(duration)) 0 50])
   %ax(3).YTick = [20 22.5 25 27.5 30 32.5 35 37.5 40];
   set(ax(3),'Color','none')
   set(ax(3),'YAxisLocation','right')
   set(ax(3),'XAxisLocation','bottom')
   datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
   ylabel('rel. transmit power', 'Fontsize', font_size, 'Fontweight', 'b');  
   set(gca,'Fontsize',font_size,'Fontweight','b');
   %legend('show', 'Location','southwest')
   set(legend(ax(3)),'Color','white','Location','southeast')
   %change backgroud color to transparent
   
   % plot Surface pressure right y-axis of the lower plot
   ax(4) = axes('Position',get(ax(2),'Position'));
   plot(duration, surf_P, 'k-','LineWidth', 1, 'DisplayName','Surf P') 
   axis([fix(min(duration)) ceil(max(duration)) -inf inf])
 %  axis([fix(min(duration)) ceil(max(duration)) 0.975 1])  %PRECIP
   set(ax(4),'Color','none')
   set(ax(4),'YAxisLocation','right')
   set(ax(4),'XAxisLocation','bottom')
   datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
   ylabel('surface pressure, atm', 'Fontsize', font_size, 'Fontweight', 'b');  
   set(gca,'Fontsize',font_size,'Fontweight','b');
   legend('show')
   set(legend(ax(3)),'Color','white')
   %change backgroud color to transparent
      
%   legend(ax(1),'Location','NorthWest') 
%   legend(ax(2),'Location','NorthWest') 
   linkaxes(ax, 'x');
   hold off;
 end
 
%%
if flag.plot_sonde_data==1
  % the elevation needs to be read from the json file
  %  elevation= 310.0; %MPD05 was at 310m elevation at SGP
  %d=pwd;
  %cd('/Volumes/documents/WV_DIAL_data/SGP_sondes/') % point to the directory where data is stored
   cd('/scr/sci/tammy/mpd/sgp/soundings/')
   [sondefilename, sondedir] = uigetfile('*.*','Select the sonde file', 'MultiSelect', 'on');
   jj=1;
  %cd(dd);
  for jj = 1:size(sondefilename,2)
      cd('/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_processing/')
      Sonde_read_nc_files(jj, elevation, sondedir, sondefilename);
  end
end 

%% 
if flag.save_figs==1
  %cd('/Users/spuler/Desktop/WV_DIAL_data/plots/') % point to the directory where data is stored 
  %cd('/Volumes/documents/WV_DIAL_data/plots/') % point to the directory where data is stored 
  %cd('/Users/spuler/Desktop/mpd_05_processed_data/Plots/') % point to the directory where data is stored   
  %cd('/Users/lroot/Desktop/mpd/Plots/') % point to the directory where data is stored 
  cd(plot_path)
  date=datestr(mean(time_new,'omitnan'), 'yyyymmdd');
  

  FigH = figure(1);
  set(gca,'Fontsize',16,'Fontweight','b');  
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920 250]);
  %set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 750 250]);
  name=strcat(date, node, ' WV_Matlab_multi'); 
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
 
  FigH = figure(2);
  set(gca,'Fontsize',16,'Fontweight','b');  
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920 250]);
  name=strcat(date, node, ' RB_Matlab_multi'); 
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
  
%   FigH = figure(3);
%   set(gca,'Fontsize',16,'Fontweight','b');
%   set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920 250]);
%   name=strcat(date, node, 'background_multi'); 
%   print(FigH, name, '-dpng', '-r0') % set at the screen resolution
  
%  FigH = figure(5);
% %  set(gca,'Fontsize',36,'Fontweight','b');
%  set(FigH, 'PaperUnits', 'points', 'PaperPosition', Scrnsize);
%  name=strcat(date, 'column_OD_multi'); 
%  print(FigH, name, '-dpng', '-r0') % set at the screen resolution
  
%   if WS==1
%       %size2 = [scrsz(4)/1 scrsz(4)/1 scrsz(3)/0.35 scrsz(4)/1]; % use for long plots 
%       size2 = [0 0 scrsz(3) scrsz(4)/3]; % use for long plots 
%       FigH = figure(4);
%      % set(gca,'Fontsize',36,'Fontweight','b');
%       set(FigH, 'PaperUnits', 'points', 'PaperPosition', size2);
%       name=strcat(date, node, 'Housekeeping'); 
%       print(FigH, name, '-dpng', '-r0') % set at the screen resolution
%   end
   
end
 
  name=strcat(date, '_WV_2CH'); 
 
end


 
toc
