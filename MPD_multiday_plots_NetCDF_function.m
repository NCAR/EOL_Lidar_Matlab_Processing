function[] = MPD_multiday_plots_NetCDF_function(save_figs, save_data, near, afterpulse, node, daystr, daystr2, skip)
%clear all; 
%start_date = '20190405';
%stop_date = '20190410';
%save_figs = 1; save_data=0; near= 0; afterpulse=0; node='MPD3'; daystr=start_date; daystr2=stop_date;
close all;
tic
dd = pwd; % get the current path

date = datestr(datenum(daystr,'yyyymmdd'), 'dd mmm yyyy');
days = (datenum(daystr2,'yyyymmdd') - datenum(daystr,'yyyymmdd'))+1;
%skip = 2;

flag.near = near;  % read in the near range channel (0=off 1=on)
flag.afterpulse = afterpulse; % read in the afterpulse corrected data (0=off 1=on)

font_size = 28; % use this for 2019b version
WS_font_size = 24; %12; % use this for 2019b version
%font_size = 36; % use this for 2018a version
%font_size = 14; % use this for 2015a version
%font_size = 16; % use this for 2015a versionexit
%font_size = 28; % use this for 2014a version

% set to 1 for using the weather station data (after Jan 2016)
WS = 1;
%include sonde data (0=off 1=on)
flag.sonde = 0;
%replot time vs range images at start of processing (0=off 1=on)
flag.replot = 1;
%save figures at end of processing (0=off 1=on)
flag.save_figs = save_figs;
%save data at end of processing (0=off 1=on)
flag.save_data = save_data;
%plot NetCDF sonde data on top (0=off 1=on)
flag.plot_sonde_data = 0;
%decimate figures to the screen 2x size (1900 pixels x2)
flag.decimate = 1;
%set the size of the range gridding
range_grid_size = 75;


RB_scale = 1;
C = importdata('NCAR_C_Map.mat');


if strcmp(getenv('HOSTNAME'),'fog.eol.ucar.edu')
   serv_path = '/export/fog1/rsfdata/MPD/'; % when running on server
else
   serv_path = '/Volumes/eol/fog1/rsfdata/MPD/'; % 
end

if strcmp(node,'MPD01')==1
    cd(strcat(serv_path, 'mpd_01_processed_data/Matlab'))  
    plot_dir = (strcat(serv_path, 'mpd_01_processed_data/Plots'));
  elseif strcmp(node,'MPD02')==1
    cd(strcat(serv_path, 'mpd_02_processed_data/Matlab'))
    plot_dir = (strcat(serv_path, 'mpd_02_processed_data/Plots'));
  elseif strcmp(node,'MPD03')==1
    cd(strcat(serv_path, 'mpd_03_processed_data/Matlab'))
    plot_dir = (strcat(serv_path, 'mpd_03_processed_data/Plots'));
    if flag.afterpulse == 1
      cd(strcat(serv_path, 'mpd_03_processed_data/Matlab/afterpulse'))
      plot_dir = (strcat(serv_path, 'mpd_03_processed_data/Matlab/afterpulse'));
    end
  elseif strcmp(node,'MPD04')==1
    cd(strcat(serv_path, 'mpd_04_processed_data/Matlab'))
    plot_dir = (strcat(serv_path, 'mpd_04_processed_data/Plots'));
    if flag.afterpulse == 1
      cd(strcat(serv_path, 'mpd_04_processed_data/Matlab/afterpulse'))
      plot_dir = (strcat(serv_path, 'mpd_04_processed_data/Matlab/afterpulse'));
    end
  elseif strcmp(node,'MPD05')==1
    cd(strcat(serv_path, 'mpd_05_processed_data/Matlab'))
    plot_dir = (strcat(serv_path, 'mpd_05_processed_data/Plots'));
end

nodeStr = extractAfter(node, 'MPD');
write_data_folder = strcat(serv_path, 'mpd_', nodeStr, '_processed_data/Matlab') 


%gate = round((bin_duration*1e-9*3e8/2)*10)/10

test_gate = exist('gate') % check for early versions were the gate wasn't saved
if test_gate == 0
 gate = 75
end

i=1;


for i=1:days
  if i==1  
    % load the near range data or regular data depending on flag  
    if flag.near == 1
      if exist(strcat(date, '_near.mat'))==2
        load(strcat(date, '_near.mat'))
      end
    else
      if exist(strcat(date, '.mat'))==2
        load(strcat(date, '.mat'))
      end
    end 
    range_limit_N_avg = size(N_avg,2);
    % grid everything to a 75 m gate size 
      if gate < 75
         range_grid_75 = 0:range_grid_size:(range_limit_N_avg-1)*gate; 
         N_avg_grid = interp1(range, N_avg', range_grid_75, 'linear', 'extrap')'; 
         RB_grid = interp1(range, RB', range_grid_75, 'linear', 'extrap')';
       %  OD = interp1(range, OD', range_grid_75, 'linear', 'extrap')';
         range = range_grid_75;
     end
    N_avg_comb=N_avg_grid;
    RB_comb=RB_grid;
  %  OD_comb=OD;
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
      t_base = T_base;
    end
  else
    date = datestr(addtodate(datenum(date), 1, 'day'), 'dd mmm yyyy');
    if flag.near == 1
      if exist(strcat(date, '_near.mat'))==2
        load(strcat(date, '_near.mat'))
      end
    else
      if exist(strcat(date, '.mat'))==2
        load(strcat(date, '.mat'))
      end
    end  
    range_limit_N_avg = size(N_avg,2);
    % grid everything to a 75 m gate size 
      if gate < 75
         range_grid_75 = 0:range_grid_size:(range_limit_N_avg-1)*gate; 
         N_avg_grid = interp1(range, N_avg', range_grid_75, 'linear', 'extrap')'; 
         RB_grid = interp1(range, RB', range_grid_75, 'linear', 'extrap')';
  %       OD = interp1(range, OD', range_grid_75, 'linear', 'extrap')';
         range = range_grid_75;
      end
    range_lim_comb = size(N_avg_comb,2); % catch any changes in range
    if gate < 75
      range_lim = size(N_avg,2)/(75/gate); % catch any changes in range
    else
      range_lim = size(N_avg,2); % catch any changes in range
    end
      range_limit = min([range_lim_comb range_lim]);
    
    
    N_avg_comb = vertcat(N_avg_comb(:,1:range_limit), N_avg_grid(2:end,1:range_limit));
    RB_comb = vertcat(RB_comb(:,1:range_limit), RB_grid(2:end,1:range_limit));
 %   OD_comb= vertcat(OD_comb(:,1:range_limit), OD(2:end,1:range_limit));
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
      t_base = vertcat(t_base, T_base(2:end,:));  
    end
  end
end

%stop


 %% save data
  
 if flag.save_data == 1
  %cd('/Users/spuler/Desktop/WV_DIAL_data') % point to the directory where data is stored 
  cd(write_data_folder)
  name=strcat(date, '_combined');
  if strcmp(node,'MPD01')==1
      MPD01.N_avg_comb = N_avg_comb;
      MPD01.RB_comb = RB_comb;
      MPD01.range = range;
      MPD01.time = duration;
      save(name, 'MPD01')
  elseif strcmp(node,'MPD02')==1
      MPD02.N_avg_comb = N_avg_comb;
      MPD02.RB_comb = RB_comb;
      MPD02.range = range;
      MPD02.time = duration;
      save(name, 'MPD02')
  elseif strcmp(node,'MPD03')==1
      MPD03.N_avg_comb = N_avg_comb;
      MPD03.RB_comb = RB_comb;
      MPD03.range = range;
      MPD03.time = duration;
      save(name, 'MPD03')
  elseif strcmp(node,'MPD04')==1
      MPD04.N_avg_comb = N_avg_comb;
      MPD04.RB_comb = RB_comb;
      MPD04.range = range;
      MPD04.time = duration;
      save(name, 'MPD04')
   elseif strcmp(node,'MPD05')==1
      MPD05.N_avg_comb = N_avg_comb;
      MPD05.RB_comb = RB_comb;
      MPD05.range = range;
      MPD05.time = duration;
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
 
scrsz = get(0,'ScreenSize');
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
 caxis([0 8]);
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
   hh = title({[date, node, ' Water Vapor (g m^{-3})']},'fontweight','b','fontsize',font_size);
 else
   datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
   hh = title({[node, ' Water Vapor (g m^{-3})']},'fontweight','b','fontsize',font_size);
 end
 set(gca,'Fontsize',font_size,'Fontweight','b');
 
 
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
 axis([fix(min(duration)) ceil(max(duration)) 0 12])
 caxis([1e1 1e7]);
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
   hh = title({[date, node, ' Attenuated Backscatter (A.U.)']},'fontweight','b','fontsize',font_size);
 else
   datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
 %  hh = title({'DIAL Relative Backscatter (C/ns km^2)'},'fontweight','b','fontsize',font_size);
    hh = title({[node, ' Attenuated Backscatter (A.U.)']},'fontweight','b','fontsize',font_size);
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
    hh = title({[date, node, ' Background']},'fontweight','b','fontsize',font_size);
  else
    datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
    hh = title({[node, ' Background']},'fontweight','b','fontsize',font_size); 
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
 %    font_size = 14
 % plot housekeeping data
   figure1 = figure('Position',Scrnsize);
   subplot1=subplot(2,1,1,'Parent',figure1,'YGrid','on', 'XGrid','on');
   box(subplot1,'on');
   hold(subplot1,'all');
  %plot(duration, (i_off),'b','LineWidth',2,'DisplayName','i_{off}') % these plot diode Temps
  %plot(duration, (i_on),'r','LineWidth',2, 'DisplayName','i_{on}')
  plot(duration, (lambda_comb_on),'g-','LineWidth',2,'DisplayName','Lambda_{on}') % these plot diode Temps
  %plot(duration, (t_hsrl),'g','LineWidth',2, 'DisplayName','T_{hsrl}')
   axis([fix(min(duration)) ceil(max(duration)) 828.1 828.3])
   YTick = [100 120 140 160 180];
   ylabel('wavelength, nm', 'Fontsize', WS_font_size, 'Fontweight', 'b');  
   datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
   set(gca,'Fontsize',WS_font_size,'Fontweight','b');
   % Plot the temperature data
   subplot2=subplot(2,1,2,'Parent',figure1,'YGrid','on', 'XGrid','on');
   box(subplot2,'on');
   hold(subplot2,'all');
   plot(duration, t_bench,'g-', 'LineWidth',1, 'DisplayName','T bench')
   plot(duration, t_base,'g', 'LineWidth',1, 'DisplayName','T base')
   plot(duration, surf_T, 'b', 'LineWidth',1, 'DisplayName','Surface T')
   axis([fix(min(duration)) ceil(max(duration)) 0 40]);% -inf inf]);   % -20 40])
      YTick = [-25 0 25 50];
   ylabel('temperature, C', 'Fontsize', WS_font_size, 'Fontweight', 'b'); 
   datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
   set(gca,'Fontsize',WS_font_size,'Fontweight','b');
   % Create legend
   legend(subplot1,'show'); 
   legend(subplot2,'show');
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
   ylabel('rel. trans. power', 'Fontsize', WS_font_size, 'Fontweight', 'b');  
   set(gca,'Fontsize',WS_font_size,'Fontweight','b');
   legend('show', 'Location','southwest')
   set(legend(ax(3)),'Color','white')
   %change backgroud color to transparent
   
   % plot Surface pressure right y-axis of the lower plot
   ax(4) = axes('Position',get(ax(2),'Position'));
   plot(duration, surf_P, 'k-','LineWidth', 1, 'DisplayName','Surf P') 
   axis([fix(min(duration)) ceil(max(duration)) 0.8 1.0]);% -inf inf]); 
   set(ax(4),'Color','none')
   set(ax(4),'YAxisLocation','right')
   set(ax(4),'XAxisLocation','bottom')
   datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
   ylabel('pressure, atm', 'Fontsize', WS_font_size, 'Fontweight', 'b');  
   set(gca,'Fontsize',WS_font_size,'Fontweight','b');
   legend('show')
   set(legend(ax(3)),'Color','white')
   %change backgroud color to transparent
      
   legend(ax(1),'Location','NorthWest') 
   legend(ax(2),'Location','NorthWest') 
   linkaxes(ax, 'x');
   hold off;
 end
 
%%
if flag.plot_sonde_data==1
  % the elevation needs to be read from the json file
    elevation= 310.0; %MPD05 was at 310m elevation at SGP
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
  cd(plot_dir)
    
  date=datestr(nanmean(time_new), 'yyyymmdd');
  
  %Scrnsize = [scrsz(4)/2 scrsz(4)/10 scrsz(3)/1 scrsz(4)/2]; % use for standard plots
  Scrnsize = [scrsz(4)/1 scrsz(4)/1 scrsz(3)/0.35 scrsz(4)/2.05]; % use for long plots 
  Scrnsize = [scrsz(4)/1 scrsz(4)/1 scrsz(3)/0.30 scrsz(4)/2]; % use for ILRC really long plots
  %Scrnsize = [0 0 scrsz(3) scrsz(4)/2]; % x y width height
  %Scrnsize = [scrsz(4)/1 scrsz(4)/1 scrsz(3)/0.51 scrsz(4)/2]; % use for Perdigao BAMS plots 
  %Scrnsize = [scrsz(4)/1 scrsz(4)/1 scrsz(3)/2 scrsz(4)/2]; % use for day plots 
  %Scrnsize = [scrsz(4)/1 scrsz(4)/1 scrsz(3)/1 scrsz(4)/2.2]; % use for AMT sized 3-day plots (with large font)
  
  FigH = figure(1);
  set(gca,'Fontsize',font_size,'Fontweight','b'); % use for Perdigao BAMS plots 
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', Scrnsize);
  if near == 1
      name=strcat(date, 'H2O_multi_near');
  else
      name=strcat(date, 'H2O_multi');
  end
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
 
  FigH = figure(2);
  set(gca,'Fontsize',font_size,'Fontweight','b'); % use for Perdigao BAMS plots 
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', Scrnsize);
  if near == 1
      name=strcat(date,  'RB_multi_near');
  else
      name=strcat(date,  'RB_multi');
  end
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
  
  FigH = figure(3);
  set(gca,'Fontsize',font_size,'Fontweight','b');
  % set(gca,'Fontsize',36,'Fontweight','b');
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', Scrnsize);
  name=strcat(date, 'background_multi'); 
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution
  
%  FigH = figure(5);
% %  set(gca,'Fontsize',36,'Fontweight','b');
%  set(FigH, 'PaperUnits', 'points', 'PaperPosition', Scrnsize);
%  name=strcat(date, 'column_OD_multi'); 
%  print(FigH, name, '-dpng', '-r0') % set at the screen resolution
  
  if WS==1
      %size2 = [scrsz(4)/1 scrsz(4)/1 scrsz(3)/0.35 scrsz(4)/1]; % use for long plots 
      %size2 = [0 0 scrsz(3) scrsz(4)/3]; % use for long plots  
      size2 = [scrsz(4)/1 scrsz(4)/1 scrsz(3)/0.35 scrsz(4)/1]; % use for ILRC really long plots
      FigH = figure(4);
      set(gca,'Fontsize',WS_font_size,'Fontweight','b');
      set(FigH, 'PaperUnits', 'points', 'PaperPosition', size2);
      %set(gca,'Fontsize',font_size,'Fontweight','b');
      %set(FigH, 'PaperUnits', 'points', 'PaperPosition', Scrnsize);
      if near == 1
        name=strcat(date,  'Housekeeping_near');
      else
        name=strcat(date,  'Housekeeping');
      end 
      print(FigH, name, '-dpng', '-r0') % set at the screen resolution
  end
   
end
 
  name=strcat(date, '_WV_2CH'); 
 
end


cd(dd) % point back to original directory 
toc
%end