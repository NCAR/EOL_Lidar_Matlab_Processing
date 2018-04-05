function[] = DIAL_sonde_comparison_funct_Python(N_H2O, sonde_top, sonde_range, t, date, T_sonde, P_sonde, sonde_stop, shift,error_threshold, Wind_speed, save_figs)
%clear all;
%close all;


scrsz = get(0,'ScreenSize');


file = date(1:end-7)
 if strcmp(date(end-4:end),'00:00')
     'time shift may not work properly because the hour is 00:00'
 end
 
%date = '10 Jul 2014';
%time = '02:00'; 


%include sonde data 0=off 1=on
sonde = 1;
%replot time vs range images at start of processing 0=off 1=on
replot = 0;
%save figures at end of processing 0=off 1=on
%save_figs = 1; 
% range offset
range_offset = 0.0;  % this is corrected in the data now.   correct for the pulse length (see WVDIAL_gate_designation_offset.numbers)
                        % post PECAN this will be -37.5m (no more delay)
%bin duration in ns
bin_duration = 500;  % ns 
%number of accumulations 
accum = 10000;  

gate = round((bin_duration*1e-9*3e8/2)*10)/10;

C = importdata('/Users/spuler/Desktop/WV_DIAL/Matlab/NCAR_C_Map.mat');
%cd('/Users/spuler/Desktop/WV_DIAL_data/') % point to the directory where data is stored 
dd = pwd; % get the current path
cd('/Volumes/documents/WV_DIAL_data/processed_data') % point to the directory where data is stored 
%cd('/Users/spuler/Desktop/FRAPPE_PECAN') % point to the directory where data is stored 

%data = load(strcat(file, 'FF.mat'));
data = load(strcat(file, 'Python.mat'));

N_avg_FF=data.N_avg;
N_error_FF=data.N_error;
%RB_FF=data.RB;
%OD_FF=data.OD;
%background_FF = data.background_off;
%time_new_FF = data.time_new;
range = data.range;
time_new =  data.time_new;
%T = data.T;
%P = data.P;


 for i=1:length(time_new)
   if (t<=time_new(i)) == 1 
     start_time_i = i -shift; %shift time back in time number of profiles
     break
   end
 end
 
 for i=1:length(time_new)
   if (sonde_stop<=time_new(i)) == 1 
     end_time_i = i -shift; %shift time back in time number of profiles
     break
   end
 end
 
 start_time = int16(start_time_i);  
 end_time = int16(end_time_i);

  g_FF_ave = nanmean(N_avg_FF(start_time:end_time,:));
  g_error_FF_ave = nanmean(N_error_FF(start_time:end_time,:));
    
  % calcuate how many 1 min profiles are being averaged (using 30 sec profiles)
  ave_num_FF = sum(~isnan(N_avg_FF(start_time:end_time,:)),1)/2; %./60;
  % correct error for longer averaging time
  g_error_FF_ave = bsxfun(@times, 1./(sqrt(ave_num_FF)), g_error_FF_ave);

  g_H2O = N_H2O.*1e6./6.022E23.*18.015;
  dial_range = range./1e3-range_offset;
  
  g_error_FF_ave(isnan(g_FF_ave)) = nan;
  g_FF_ave(isnan(g_error_FF_ave)) = nan;
  g_FF_ave(g_error_FF_ave./g_FF_ave>error_threshold) = nan; % remove high error values
  g_FF_ave(g_error_FF_ave./g_FF_ave<-error_threshold) = nan; % remove high error values
  
    
  figure(10)
  plot(g_FF_ave, dial_range, 'r', 'LineWidth', 2)
  ylim([0 6])
  xlim([0 18])
  title({[date, ' UTC']}, 'Fontsize', 30, 'Fontweight', 'b')
  xlabel('Water Vapor (g/m^{3})',  'Fontsize', 30, 'Fontweight', 'b');
  ylabel('Height (km, AGL)', 'Fontsize', 30, 'Fontweight', 'b');
  % overlay Near Field Channel
  figure(10)
  if sonde==1
    hold on
    if sonde_top >= 2500 
      plot(g_H2O(1:2500,1),sonde_range(1:2500,1), '-b', 'MarkerSize', 2, 'LineWidth', 2)
    else
      plot(g_H2O(1:sonde_top,1),sonde_range(1:sonde_top,1), '-b', 'MarkerSize', 2, 'LineWidth', 2)
    end
    hold off   
    %legend('1 min Far Field', '25 min Far Field', '1 min Near Field', '25 min Near Field', 'Sonde','Location', 'NorthEast');
    legend('Far Field', 'Sonde','Location', 'NorthEast');
  else
   legend('Far Field', 'Location', 'NorthEast');
  end
  grid on
  % add error bars
  hold on
    cd('/Users/spuler/Desktop/WV_DIAL/Matlab/')
    hb = herrorbar(g_FF_ave, dial_range, g_error_FF_ave); 
    set(hb, 'Color', 'r');
  hold off
  legend('Far Field', 'Sonde' ,'Location', 'NorthEast');
  grid on
  set(gca, 'FontSize', 22) 
    

if save_figs==1
% print the figures to file 
  %cd('/Users/spuler/Desktop/WV_DIAL_data/plots/') % point to the directory where data is stored 
  cd('/Volumes/documents/WV_DIAL_data/plots/') % point to the directory where data is stored 
  
  date=datestr(nanmean(t), 'yyyymmdd HHMM');
  FigH = figure(10);
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', [scrsz(4)/2 scrsz(4)/10 scrsz(3)/2 scrsz(4)/1.5]);
  name=strcat(date, '_WV_2CH'); 
  print(FigH, name, '-dpng', '-r300') % set the resolution as 300 dpi

%  FigH = figure(20);
%  set(FigH, 'PaperUnits', 'points', 'PaperPosition', [scrsz(4)/2 scrsz(4)/10 scrsz(3)/2 scrsz(4)/1.5]);
%  name=strcat(date, '_RB_2CH'); 
%  print(FigH, name, '-dpng', '-r300') % set the resolution as 300 dpi
end
  

if sonde==1  
%  % use this section to help set the temperature and pressure profile
  %T0 = 273+25; % set to match the sounding
  %T = T0-0.0065.*range; % set to match the sounding
  % pressure in atmospheres
  %P0  = 0.83; % set to match the sounding
  %P = P0.*(T0./T).^-5.5;   % set this value to match sounding

  %compare estimated T with measured T profile
  %figure(12)
  %plot(T_sonde(1:sonde_top,1)+273,sonde_range(1:sonde_top,1), 'g')
  %hold on
  %plot(T,range/1000)
  %ylim([0 8])
  %title({[date, ' Temp']}, 'Fontsize', 30, 'Fontweight', 'b')
  %hold off

  %compare estimated P with measured P profile
  %figure(13)
  %plot(P_sonde(1:sonde_top,1)./1013.25,sonde_range(1:sonde_top,1), 'g')
  %hold on
  %plot(P,range/1000)
  %ylim([0 8])
  %xlim([0 1])
  %title({[date, ' Pressure']}, 'Fontsize', 30, 'Fontweight', 'b')
  %hold off
end   



 name=strcat(date, '_WV_2CH'); 
 
  
 


%cd('/Users/spuler/Desktop/WV_DIAL_data') % point to the directory where data is stored 
  cd('/Volumes/documents/WV_DIAL_data/processed_data') % point to the directory where data is stored 
  name=strcat(date, '_sonde');
  save(name, 'g_FF_ave', 'dial_range', 'g_H2O', 'sonde_range',...
       'g_error_FF_ave')
  