clear all; close all;
tic

 node = 'MPD05';
 date = '20210620';   
 days = 36; skip = 4;
 date = '20210803';  
 days = 4; skip = 1;
 lapse_rate = 0.0085;
     
serv_path = '/Volumes/documents/MPD';
plot_path = '/Volumes/documents/MPD/Plots/'; 
%flag.near = 0;  % read in the near range channel (0=off 1=on)
font_size = 36; % use this for 2018a version

%replot time vs range images at start of processing (0=off 1=on)
flag.replot = 1;
%save figures at end of processing (0=off 1=on)
flag.save_figs = 1;
%save data at end of processing (0=off 1=on)
flag.save_data = 0;
%plot NetCDF sonde data on top (0=off 1=on)
flag.plot_sonde_data = 0;
%decimate figures to the screen 2x size (1900 pixels x2)
flag.decimate = 1;
  
C = importdata('NCAR_C_Map.mat');
dd=pwd;

if strcmp(node,'MPD01')==1
   cd(strcat(serv_path,'/mpd_01_processed_data//Matlab_temp'))
 elseif strcmp(node,'MPD02')==1
   cd(strcat(serv_path,'/mpd_02_processed_data//Matlab_temp'))
 elseif strcmp(node,'MPD03')==1
   cd(strcat(serv_path,'/mpd_03_processed_data//Matlab_temp'))
 elseif strcmp(node,'MPD04')==1
   cd(strcat(serv_path,'/mpd_04_processed_data//Matlab_temp'))
 elseif strcmp(node,'MPD05')==1
   cd(strcat(serv_path,'/mpd_05_processed_data//Matlab_temp'))
end

i=1;


for i=1:days
  if i==1  
    if exist(strcat('mpd05.', date, '.Matlab.mat'))==2
      load(strcat('mpd05.', date, '.Matlab.mat'),'Retrievals','Data')  
      % reading in the 'Data' variable slows down the read considerably
      % save the file as save(strcat('mpd05.', date, '.Matlab.mat'), '-struct', 'Retreivals','Data');
      % then can load with the following line more quickly  
    end
    Temp_comb = Retrievals.Temperature.Value';
    Temp_comb_avg = Retrievals.Temperature.Smoothed';
    duration =  (datenum(date, 'yyyymmdd') + Retrievals.Temperature.TimeStamp/3600/24)';
    range = Retrievals.Temperature.Range;
    % following four lines are to build a lapse rate temp field based on surface station data
      T_surf = Data.TimeSeries.WeatherStation.Temperature;
      WS_time = Data.TimeSeries.WeatherStation.TimeStamp/24;
      T_surf_grid = interp1( WS_time, T_surf, Retrievals.Temperature.TimeStamp/3600/24, 'linear')';
      T_lapse =  T_surf_grid-lapse_rate*range;

%     figure(1)
%      plot(Retrievals.Temperature.TimeStamp/3600/24)
%      hold on
%      plot(Data.TimeSeries.WeatherStation.TimeStamp/24)


  else   
    date = datestr(addtodate(datenum(date, 'yyyymmdd'), 1, 'day'), 'yyyymmdd');
    if exist(strcat('mpd05.', date, '.Matlab.mat'))==2
      load(strcat('mpd05.', date, '.Matlab.mat'),'Retrievals', 'Data')
    end
    % following four lines are to build a lapse rate temp field based on surface station data
      T_surf = Data.TimeSeries.WeatherStation.Temperature;
      WS_time = Data.TimeSeries.WeatherStation.TimeStamp/24;
      T_surf_grid = interp1( WS_time, T_surf, Retrievals.Temperature.TimeStamp/3600/24, 'linear')';
      Temp_lapse =  T_surf_grid-lapse_rate*range;
      
    T_lapse = vertcat(T_lapse, Temp_lapse);
    Temp_comb = vertcat(Temp_comb, Retrievals.Temperature.Value');
    Temp_comb_avg = vertcat(Temp_comb_avg, Retrievals.Temperature.Smoothed');
    duration = vertcat(duration, (datenum(date, 'yyyymmdd') + Retrievals.Temperature.TimeStamp/3600/24)');
  end
end

%stop

 
  if days == 1
   xData =  linspace(fix(min(duration)),  ceil(max(duration)), 25);
 else
   xData =  linspace( fix(min(duration)),  ceil(max(duration)), round((ceil(max(duration))-fix(min(duration)))/skip)+1 );
  % xData =  linspace(fix(min(duration)),  round(max(duration)), 36);
  end

% point back to original directory 
cd(dd);

scrsz = [1  1  1920 1200];
Scrnsize=[scrsz(4)/2 scrsz(4)/10 scrsz(3)/1 scrsz(4)/2];


 x = duration';
 y = range./1000;
 Z = Temp_comb_avg'; 
 
 
%if flag.replot==1
 % plot Narrow water vapor in g/m^3
 figure('Position',Scrnsize)
 %Z = double(real(N_avg_comb'.*1e6./6.022E23.*18.015));  %number density in mol/cm3(1e6 cm3/m3)/(N_A mol/mole)*(18g/mole)
 % Z(isnan(Z)) = -1;
 set(gcf,'renderer','zbuffer');
 h = pcolor(x,y,Z);
 set(h, 'EdgeColor', 'none');
 colorbar('EastOutside');
 axis([fix(min(x)) ceil(max(x)) 0 6])
 caxis([240 320]);
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
   hh = title({[date,'MPD05 temperature (K)']},'fontweight','b','fontsize',font_size);
 else
   datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
   hh = title({[node, ' Temperature (K)']},'fontweight','b','fontsize',font_size);
   %hh = title({'DIAL Water Vapor (g m^{-3})'},'fontweight','b','fontsize',font_size);
 end
 set(gca,'Fontsize',font_size,'Fontweight','b');
 
 
 Z = Temp_comb'; 

 % plot averaged Temp K
 figure('Position',Scrnsize)
 %Z = double(real(N_avg_comb'.*1e6./6.022E23.*18.015));  %number density in mol/cm3(1e6 cm3/m3)/(N_A mol/mole)*(18g/mole)
 % Z(isnan(Z)) = -1;
 set(gcf,'renderer','zbuffer');
 h = pcolor(x,y,Z);
 set(h, 'EdgeColor', 'none');
 colorbar('EastOutside');
 axis([fix(min(x)) ceil(max(x)) 0 6])
 caxis([240 320]);
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
   hh = title({[date,'MPD05 temperature (K)']},'fontweight','b','fontsize',font_size);
 else
   datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
   hh = title({[node, ' Temperature (K)']},'fontweight','b','fontsize',font_size);
   %hh = title({'DIAL Water Vapor (g m^{-3})'},'fontweight','b','fontsize',font_size);
 end
 set(gca,'Fontsize',font_size,'Fontweight','b');
 
 
%  % plot Narrow RB
%  figure('Position',Scrnsize)
% % Z = double((real(RB_comb')./RB_scale));
% % Z_mask = Z;
% % Z_mask(RB_FF'<5) = NaN;
%  % Z(isnan(Z)) = -1;
%  set(gcf,'renderer','zbuffer');
%  h = pcolor(x,y,Z_RB);
% %  h = pcolor(x,y,Z_mask);
%  set(h, 'EdgeColor', 'none');
%  colorbar('EastOutside');
%  axis([fix(min(duration)) ceil(max(duration)) 0 12])
%  caxis([1e1 1e7]);
%  colormap(C)
%  %colormap(perula)
%  %shading interp
%  % P_t = get(hh, 'Position');
%  % set(hh,'Position', [P_t(1) P_t(2)+0.2 P_t(3)])
%  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size);
%  set(gca, 'XTick',  xData)
%  set(gca,'TickDir','out');
%  set(gca,'TickLength',[0.005; 0.0025]);
%  if days == 1
%    datetick('x','HH:MM','keeplimits', 'keepticks');
%    xlabel('Time (UTC)','fontweight','b','fontsize',font_size);
%   % hh = title({[date,'DIAL Relative Backscatter (C/ns km^2)']},'fontweight','b','fontsize',font_size);
%    hh = title({[date,'DIAL Attenuated Backscatter (A.U.)']},'fontweight','b','fontsize',font_size);
%  else
%    datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
%  %  hh = title({'DIAL Relative Backscatter (C/ns km^2)'},'fontweight','b','fontsize',font_size);
%     hh = title({'DIAL Attenuated Backscatter (A.U.)'},'fontweight','b','fontsize',font_size);
%  end
%  set(gca,'Fontsize',font_size,'Fontweight','b');
%  set(gca,'Zscale', 'log')
%  set(gca,'Colorscale', 'log')
%  set(gca,'Zscale', 'linear')
% 
%  
% % plot the background
%   figure('Position',Scrnsize)
%  % if near_field==1 
%  %   len1= length(duration);
%  %   len2 = length(background_comb_on);
%  %   len3 = length(background_comb_off); % HSRL molecular channel
%  %
%  %   len = vertcat(len1, len2, len3);
%  %   plot_end = min(len);  
%  %   semilogy(duration(1:plot_end), background_comb_on(1:plot_end), 'r', duration(1:plot_end), background_comb_off(1:plot_end), 'k') 
%  % else
%     semilogy(duration, background_comb_off, 'black')
%  % end
%   axis([fix(min(duration)) ceil(max(duration)) 1e2 1e7])
%   ylabel('Offline background C/s', 'Fontsize', font_size, 'Fontweight', 'b');  
%   grid on
%   %hold on
%   %  semilogy(duration, 5e6, 'green', 'Linewidth', 2)    % Add a horizontal line to show the 5Mc/s linearity limit
%   %hold off
%  % if near_field==1  
%  % legend('WV Online  ', 'WV Offline  ', 'HSRL molecular  ', 'HSRL combined  ','Location', 'NorthWest');
%  % else
%  legend('WV Online  ', 'WV Offline  ', 'Location', 'NorthWest');
%  % end
%   set(gca, 'XTick',  xData)
%   set(gca,'TickDir','out');
%   set(gca,'TickLength',[0.005; 0.0025]);
%   if days == 1
%     datetick('x','HH:MM','keeplimits', 'keepticks');
%     xlabel('Time (UTC)','fontweight','b','fontsize',font_size); 
%     hh = title({[date,'Background']},'fontweight','b','fontsize',font_size);
%   else
%     datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
%     hh = title({'DIAL Background'},'fontweight','b','fontsize',font_size); 
%   end
%   set(gca,'Fontsize',font_size,'Fontweight','b');
%   
%   
% % plot column OD for the Narrow  
%   % figure('Position',Scrnsize)
%   % Z = double(real(OD_comb'));
%   % set(gcf,'renderer','zbuffer');
%   % h = pcolor(x,y,Z);
%   % set(h, 'EdgeColor', 'none');
%   % colorbar('EastOutside');
%   % axis([fix(min(duration)) ceil(max(duration)) 0 12])
%   % caxis([-0.1 2]);
%   % colormap(C)
%   %title({('  Column Optical Depth, Narrow Field  ')},...
%   %     'fontweight','b','fontsize',font_size)
%   % xlabel('Time (UTC)','fontweight','b','fontsize',font_size); 
%   % ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
%   % set(gca, 'XTick',  xData)
%   % set(gca,'TickDir','out');
%   % set(gca,'TickLength',[0.005; 0.0025]);
%   % if days == 1
%   %   datetick('x','HH:MM','keeplimits', 'keepticks');
%   %   xlabel('Time (UTC)','fontweight','b','fontsize',font_size); 
%   % else
%   %   datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
%   % end
%   % set(gca,'Fontsize',font_size,'Fontweight','b');
% %
% %  
% %  plot1 = OD_comb(:,round(2500/gate));
% %  plot2 = OD_comb(:,round(5000/gate));
% %  %window = ones(1,1)/1;
% %  %plot1s = filter2(window, plot1);
% %  %plot2s = filter2(window, plot2);
% %  dd = pwd; % get the current path
% %  cd('/Users/spuler/Desktop/WV_DIAL/Matlab/')
% %    plot1s = ndnanfilter(plot1,'rectwin', 20); 
% %    plot2s = ndnanfilter(plot2,'rectwin', 20);
% %  cd(dd)
%   
% % % plot column optical depth at 5km and 2.5km as function of time
% %  figure('Position',Scrnsize)
% %  plot(duration, plot1s, 'b',duration, plot2s, 'k', 'LineWidth', 2)
% %  axis([fix(min(duration)) ceil(max(duration)) 0 2.5])
% %  ylabel('OD', 'Fontsize', font_size, 'Fontweight', 'b');
% %  grid on
% %  legend('OD at 2.5km', 'OD at 5.0km','Location', 'NorthWest');
% %  set(gca, 'XTick',  xData)
% %  set(gca,'TickDir','out');
% %  set(gca,'TickLength',[0.005; 0.0025]);
% %  if days == 1
% %    datetick('x','HH:MM','keeplimits', 'keepticks');
% %    xlabel('Time (UTC)','fontweight','b','fontsize',font_size); 
% %    hh = title({[date,'DIAL Column OD']},'fontweight','b','fontsize',font_size);
% %  else
% %    datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
% %       hh = title({'DIAL Column OD'},'fontweight','b','fontsize',font_size); 
% %  end
% %  set(gca,'Fontsize',font_size,'Fontweight','b');
%   
%  if WS==1
%      font_size = 14
%  % plot housekeeping data
%    figure1 = figure('Position',Scrnsize);
%    subplot1=subplot(2,1,1,'Parent',figure1,'YGrid','on', 'XGrid','on');
%    box(subplot1,'on');
%    hold(subplot1,'all');
%   %plot(duration, (i_off),'b','LineWidth',2,'DisplayName','i_{off}') % these plot diode Temps
%   %plot(duration, (i_on),'r','LineWidth',2, 'DisplayName','i_{on}')
%   plot(duration, (lambda_comb_on),'b','LineWidth',2,'DisplayName','Lambda_{on}') % these plot diode Temps
%   %plot(duration, (t_hsrl),'g','LineWidth',2, 'DisplayName','T_{hsrl}')
%    axis([fix(min(duration)) ceil(max(duration)) -inf inf])
%    YTick = [100 120 140 160 180];
%    ylabel('seed Temp, C', 'Fontsize', font_size, 'Fontweight', 'b');  
%    datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
%    set(gca,'Fontsize',font_size,'Fontweight','b');
%    % Plot the temperature data
%    subplot2=subplot(2,1,2,'Parent',figure1,'YGrid','on', 'XGrid','on');
%    box(subplot2,'on');
%    hold(subplot2,'all');
%    %plot(duration, T_bench,'r', 'LineWidth',1, 'DisplayName','T bench')
%    plot(duration, surf_T, 'b', 'LineWidth',1, 'DisplayName','Surface T')
%    axis([fix(min(duration)) ceil(max(duration)) -inf inf]);   % -20 40])
%       YTick = [-25 0 25 50];
%    ylabel('temperature, C', 'Fontsize', font_size, 'Fontweight', 'b'); 
%    datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
%    set(gca,'Fontsize',font_size,'Fontweight','b');
%    % Create legend
%    legend(subplot1,'show'); 
%    legend(subplot2,'show');
%    %link the x axis for all 3 subplots
%    ax(1)=subplot(2,1,1);
%    ax(2)=subplot(2,1,2);
%    % plot power on right y-axis of the upper plot (% assumes 5% pickoff)
%    ax(3) = axes('Position',get(ax(1),'Position'));
%    plot(duration, (p_off/7500),'b-','LineWidth', 1, 'DisplayName','P_{off}') % changed from 0.05 to 0.0425
%    hold on
%    plot(duration, p_on/7500,'r-', 'LineWidth',1, 'DisplayName','P_{on}')
%    %plot(duration, p_hsrl/2500,'g--', 'LineWidth',1, 'DisplayName','P_{hsrl}')
%    axis([fix(min(duration)) ceil(max(duration)) 0 50])
%    %ax(3).YTick = [20 22.5 25 27.5 30 32.5 35 37.5 40];
%    set(ax(3),'Color','none')
%    set(ax(3),'YAxisLocation','right')
%    set(ax(3),'XAxisLocation','bottom')
%    datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
%    ylabel('rel. transmit power', 'Fontsize', font_size, 'Fontweight', 'b');  
%    set(gca,'Fontsize',font_size,'Fontweight','b');
%    legend('show', 'Location','southwest')
%    set(legend(ax(3)),'Color','white')
%    %change backgroud color to transparent
%    
%    % plot Surface pressure right y-axis of the lower plot
%    ax(4) = axes('Position',get(ax(2),'Position'));
%    plot(duration, surf_P, 'k-','LineWidth', 1, 'DisplayName','Surf P') 
%    axis([fix(min(duration)) ceil(max(duration)) -inf inf])
%    set(ax(4),'Color','none')
%    set(ax(4),'YAxisLocation','right')
%    set(ax(4),'XAxisLocation','bottom')
%    datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
%    ylabel('surface pressure, atm', 'Fontsize', font_size, 'Fontweight', 'b');  
%    set(gca,'Fontsize',font_size,'Fontweight','b');
%    legend('show')
%    set(legend(ax(3)),'Color','white')
%    %change backgroud color to transparent
%       
%    legend(ax(1),'Location','NorthWest') 
%    legend(ax(2),'Location','NorthWest') 
%    linkaxes(ax, 'x');
%    hold off;
%  end
%  
% %%
% if flag.plot_sonde_data==1
%   % the elevation needs to be read from the json file
%     elevation= 310.0; %MPD05 was at 310m elevation at SGP
%   %d=pwd;
%   %cd('/Volumes/documents/WV_DIAL_data/SGP_sondes/') % point to the directory where data is stored
%    cd('/scr/sci/tammy/mpd/sgp/soundings/')
%    [sondefilename, sondedir] = uigetfile('*.*','Select the sonde file', 'MultiSelect', 'on');
%    jj=1;
%   %cd(dd);
%   for jj = 1:size(sondefilename,2)
%       cd('/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_processing/')
%       Sonde_read_nc_files(jj, elevation, sondedir, sondefilename);
%   end
% end 
% 
% %% 
if flag.save_figs==1
  %cd('/Users/spuler/Desktop/WV_DIAL_data/plots/') % point to the directory where data is stored 
  %cd('/Volumes/documents/WV_DIAL_data/plots/') % point to the directory where data is stored 
  %cd('/Users/spuler/Desktop/mpd_05_processed_data/Plots/') % point to the directory where data is stored   
  %cd('/Users/lroot/Desktop/mpd/Plots/') % point to the directory where data is stored 
  cd(plot_path) % point to the directory where data is stored 

 
  FigH = figure(1);
  set(gca,'Fontsize',16,'Fontweight','b');  
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920 250]);
  name=strcat(date, ' Temp_Matlab_multi'); 
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
 
%   FigH = figure(2);
%   set(gca,'Fontsize',16,'Fontweight','b');  
%   set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920 250]);
%   name=strcat(date, ' RB_Matlab_multi'); 
%   print(FigH, name, '-dpng', '-r0') % set at the screen resolution    
end



 

 
toc
