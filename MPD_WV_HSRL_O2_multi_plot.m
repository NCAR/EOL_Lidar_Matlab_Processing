clear all; close all;
tic

 node = 'MPD01';
 date = '08 Oct 2025';   
 days = 8; skip = 1;
 flag.afterpulse = 0; % read in the afterpulse corrected data (0=off 1=on)
 WV_max_scale = 12;

 % node = 'MPD04';
 % date = '23 Aug 2025';
 % days = 2 ; skip = 1;
 % WV_max_scale = 25;
 % 
 % node = 'MPD04';
 % date = '7 May 2025';
 % days = 90; skip = 9;
 % WV_max_scale = 25;

 % 
 
 
serv_path = '/Volumes/smaug1/rsfdata/MPD/';
plot_path = '/Volumes/Macintosh HD/Users/spuler/Desktop/mpd/Plots/';
C = importdata('NCAR_C_Map.mat');
addpath '/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing/';
addpath '/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing/matplotlib/';
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

%% read and combine the data into a single file
for i=1:days
  if i==1  
     if exist(strcat(node, '_', datestr(date, 'yyyymmdd'), '_Backscatter Coefficient.mat'))==2
        load(strcat(node, '_', datestr(date, 'yyyymmdd'), '_Backscatter Coefficient.mat'))
         BSR_comb = BSR;
         beta_bs_comb = beta_bs;
         alpha_O2_comb = alpha_O2;
         duration=x';
         end
     if exist(strcat(node, '_', datestr(date, 'yyyymmdd'), '_WV.mat'))==2
        load(strcat(node, '_', datestr(date, 'yyyymmdd'), '_WV.mat'))
        range_limit = size(N_avg,2);
         range_limit = 490
        N_avg_comb=N_avg;
        N_error_comb=N_error;
        RB_comb=RB;
        background_comb_on = background_on;
        background_comb_off = background_off;
        lambda_comb_on = lambda_all;
        lambda_comb_off = lambda_all_off;
        durationWV=time_new;
        surf_T = Surf_T;
        surf_P = Surf_P;
        surf_AH = Surf_AH;
        i_off = I_off;
        i_on = I_on;
        p_WVon = P_on;
        p_WVoff = P_off;
        t_bench = T_bench;
     end 
  else
    date = datestr(addtodate(datenum(date), 1, 'day'), 'dd mmm yyyy');
    if exist(strcat(node, '_', datestr(date, 'yyyymmdd'), '_Backscatter Coefficient.mat'))==2
       load(strcat(node, '_', datestr(date, 'yyyymmdd'), '_Backscatter Coefficient.mat'))
       BSR_comb = vertcat(BSR_comb(:,1:range_limit), BSR(2:end,1:range_limit));
       beta_bs_comb = vertcat(beta_bs_comb(:,1:range_limit), beta_bs(2:end,1:range_limit));
       alpha_O2_comb = vertcat(alpha_O2_comb(:,1:range_limit), alpha_O2(2:end,1:range_limit));
       duration = vertcat(duration, x');
    end
    if exist(strcat(node, '_', datestr(date, 'yyyymmdd'), '_WV.mat'))==2
        load(strcat(node, '_', datestr(date, 'yyyymmdd'), '_WV.mat'))
        N_avg_comb = vertcat(N_avg_comb(:,1:range_limit), N_avg(2:end,1:range_limit));
        N_error_comb = vertcat(N_error_comb(:,1:range_limit), N_error(2:end,1:range_limit));
        RB_comb = vertcat(RB_comb(:,1:range_limit), RB(2:end,1:range_limit));
        durationWV = vertcat(durationWV, time_new(2:end));
        background_comb_on = vertcat(background_comb_on, background_on(2:end));
        background_comb_off = vertcat(background_comb_off, background_off(2:end));
        lambda_comb_on = vertcat(lambda_comb_on, lambda_all(2:end));
        lambda_comb_off = vertcat(lambda_comb_off, lambda_all_off(2:end));
        surf_T = vertcat(surf_T,Surf_T(2:end,:));
        surf_P = vertcat(surf_P,Surf_P(2:end,:));
        surf_AH = vertcat(surf_AH, Surf_AH(2:end,:));
        i_off = vertcat(i_off,I_off(2:end,:));
        i_on = vertcat(i_on, I_on(2:end,:));
        p_WVon = vertcat(p_WVon, P_on(2:end,:));  
        p_WVoff = vertcat(p_WVoff, P_off(2:end,:)); 
        t_bench = vertcat(t_bench, T_bench(2:end,:));  
    end
    
  end
end

%force range to starting size
y = y(1:range_limit);

%% plot data
cd(dd);
xData =  linspace( fix(min(duration)),  ceil(max(duration)), round((ceil(max(duration))-fix(min(duration)))/skip)+1 );

  figure(1)
  x = (durationWV)';
  Z = real(double((((BSR_comb))')));
  font_size = 14;
  set(gcf,'renderer','zbuffer');
  h = pcolor(x,y,Z);
  set(h, 'EdgeColor', 'none');
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  set(gca, 'XTick',  xData) 
  colorbar('EastOutside');
  axis([fix(min(duration))  ceil(max(duration)) 0 6])
  caxis([1 100]);
  hh = title({[node, ' ', ' Backscatter Ratio']},'fontweight','b','fontsize',font_size);  
  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size);
  set(gca,'Fontsize',font_size,'Fontweight','b');
  set(gca,'Zscale', 'log')
  set(gca,'Colorscale', 'log')
  set(gca,'Zscale', 'linear')
  colormap(jet)
  %colormap(viridis)
  grid on

  figure(2)
  Z = real(double((beta_bs_comb')));
  font_size = 14;
  set(gcf,'renderer','zbuffer');
  h = pcolor(x,y,Z);
  set(h, 'EdgeColor', 'none');
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  set(gca, 'XTick',  xData)
  colorbar('EastOutside');
  axis([fix(min(duration))  ceil(max(duration)) 0 6])
  caxis([1e-8 1e-4]);
  hh = title({[node, ' ', ' Backscatter Coefficient [m^{-1}sr^{-1}]']},'fontweight','b','fontsize',font_size);  
  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size);
  set(gca,'Fontsize',font_size,'Fontweight','b');
  set(gca,'Zscale', 'log')
  set(gca,'Colorscale', 'log')
  set(gca,'Zscale', 'linear')
  colormap(jet)
  colormap(viridis)

  figure(3)
  Z = alpha_O2_comb';
  font_size = 14;
  set(gcf,'renderer','zbuffer');
  h = pcolor(x,y,Z);
  set(h, 'EdgeColor', 'none');
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  set(gca, 'XTick',  xData)
  colorbar('EastOutside');
  axis([fix(min(duration))  ceil(max(duration)) 0 6])
  caxis([1e-6 3e-4]);
  hh = title({[node, ' ', ' O2 absorption coeff [m^{-1}]']},'fontweight','b','fontsize',font_size);  
  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size);
  set(gca,'Fontsize',font_size,'Fontweight','b');
   
  x = (durationWV)';
  y = (range(1:range_limit)./1e3);
  Z_AH = double(real(N_avg_comb'.*1e6./6.022E23.*18.015));  %number density in mol/cm3(1e6 cm3/m3)/(N_A mol/mole)*(18g/mole)
  Z_RB = double((real(RB_comb')));
  
  % plot Narrow water vapor in g/m^3
  figure(4)
  set(gcf,'renderer','zbuffer');
  h = pcolor(x,y,Z_AH);
  set(h, 'EdgeColor', 'none');
  colorbar('EastOutside');
  axis([fix(min(x)) ceil(max(x)) 0 6])
  caxis([0 WV_max_scale]);
  colormap(jet)
  colormap(CM_YlGnBu(64))
   ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
  set(gca, 'XTick',  xData)
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  hh = title({[node, ' Water Vapor (g m^{-3})']},'fontweight','b','fontsize',font_size);
  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
  set(gca,'Fontsize',font_size,'Fontweight','b');
  grid on
  
 
  % plot RB
  figure(5)
  set(gcf,'renderer','zbuffer');
  h = pcolor(x,y,Z_RB);
  set(h, 'EdgeColor', 'none');
  colorbar('EastOutside');
  axis([fix(min(duration)) ceil(max(duration)) 0 12])
  caxis([1e1 1e6]);
  colormap(C)
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size);
  set(gca, 'XTick',  xData)
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  hh = title({[node, ' Attenuated Backscatter (A.U.)']},'fontweight','b','fontsize',font_size);
  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
  set(gca,'Fontsize',font_size,'Fontweight','b');
  set(gca,'Zscale', 'log')
  set(gca,'Colorscale', 'log')
  set(gca,'Zscale', 'linear')
   
   figure6 = figure('Position',[1 1 1920 250*4]);
   font_size = 16;
   subplot1=subplot(4,1,1,'Parent',figure6,'YGrid','on', 'XGrid','on');
   box(subplot1,'on');
 %  hold(subplot1,'all');
      x = (durationWV)';
%       Z = real(double((((BSR_comb))')));
      Z = real(double((beta_bs_comb')));
      set(gcf,'renderer','zbuffer');
      h = pcolor(x,y,Z);
      set(h, 'EdgeColor', 'none');
      set(gca,'TickDir','out');
      set(gca,'TickLength',[0.005; 0.0025]);
      set(gca, 'XTick',  xData) 
      colorbar('EastOutside');
      axis([fix(min(duration))  ceil(max(duration)) 0 12])
%       caxis([1 10]);
%       hh = title({[node, ' ', ' Backscatter Ratio']},'fontweight','b','fontsize',font_size);  
      caxis([1e-8 1e-3]);
      hh = title({[node, ' ', ' Backscatter Coefficient [m^{-1}sr^{-1}]']},'fontweight','b','fontsize',font_size);  
      datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
      ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size);
      set(gca,'Fontsize',font_size,'Fontweight','b');
      set(gca,'Zscale', 'log')
      set(gca,'Colorscale', 'log')
      set(gca,'Zscale', 'linear')
      colormap(subplot1, viridis)
      colormap(subplot1, jet)
      grid on
      box(subplot1,'off'); 
  subplot2=subplot(4,1,2,'Parent',figure6,'YGrid','on', 'XGrid','on');
  box(subplot2,'on');
      x = (durationWV)';
      y = (range(1:range_limit)./1e3);
      Z_AH = double(real(N_avg_comb'.*1e6./6.022E23.*18.015));  %number density in mol/cm3(1e6 cm3/m3)/(N_A mol/mole)*(18g/mole)
      set(gcf,'renderer','zbuffer');
      h = pcolor(x,y,Z_AH);
      set(h, 'EdgeColor', 'none');
      colorbar('EastOutside');
      axis([fix(min(x)) ceil(max(x)) 0 6])
      caxis([0 WV_max_scale]);
      colormap(subplot2, CM_YlGnBu(64))
      ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
      set(gca, 'XTick',  xData)
      set(gca,'TickDir','out');
      set(gca,'TickLength',[0.005; 0.0025]);
      hh = title({[node, ' Water Vapor (g m^{-3})']},'fontweight','b','fontsize',font_size);
      datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
      set(gca,'Fontsize',font_size,'Fontweight','b');
      grid on
  subplot3=subplot(4,1,3,'Parent',figure6,'YGrid','on', 'XGrid','on');
  box(subplot3,'on'); 
      x = (durationWV)';    
      Z = alpha_O2_comb';
      set(gcf,'renderer','zbuffer');
      h = pcolor(x,y,Z);
      set(h, 'EdgeColor', 'none');
      set(gca,'TickDir','out');
      set(gca,'TickLength',[0.005; 0.0025]);
      set(gca, 'XTick',  xData)
      colorbar('EastOutside');
      axis([fix(min(duration))  ceil(max(duration)) 0 6])
      caxis([1e-6 3e-4]);
      hh = title({[node, ' ', ' O2 absorption coeff [m^{-1}]']},'fontweight','b','fontsize',font_size);  
      datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
      ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size);
      set(gca,'Fontsize',font_size,'Fontweight','b');
      colormap(subplot3, parula)
      %colormap(subplot3, jet)
  subplot4=subplot(4,1,4,'Parent',figure6,'YGrid','on', 'XGrid','on');
  box(subplot4,'on'); 
      x = (durationWV)';
      set(gcf,'renderer','zbuffer');
      h = pcolor(x,y,Z_RB);
      set(h, 'EdgeColor', 'none');
      colorbar('EastOutside');
      axis([fix(min(duration)) ceil(max(duration)) 0 12])
      caxis([1e1 1e6]);
      colormap(subplot4, C)
      ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size);
      set(gca, 'XTick',  xData)
      set(gca,'TickDir','out');
      set(gca,'TickLength',[0.005; 0.0025]);
      hh = title({[node, ' Attenuated Backscatter (A.U.)']},'fontweight','b','fontsize',font_size);
      datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
      set(gca,'Fontsize',font_size,'Fontweight','b');
      set(gca,'Zscale', 'log')
      set(gca,'Colorscale', 'log')
      set(gca,'Zscale', 'linear')
      
 
      


 % plot housekeeping data
   font_size = 14
   figure1 = figure('Position',[1 1 1920 250*2]);
   subplot1=subplot(2,1,1,'Parent',figure1,'YGrid','on', 'XGrid','on');
   box(subplot1,'on');
   hold(subplot1,'all');
   plot(durationWV, (lambda_comb_on),'k','LineWidth',2,'DisplayName','Lambda_{on}') % these plot diode Temps
   axis([fix(min(durationWV)) ceil(max(durationWV)) 828.196 828.207])
   YTick = [100 120 140 160 180];
   ylabel('wavelength, nm', 'Fontsize', font_size, 'Fontweight', 'b');  
   datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
   set(gca,'Fontsize',font_size,'Fontweight','b');
   % Plot the temperature data
   subplot2=subplot(2,1,2,'Parent',figure1,'YGrid','on', 'XGrid','on');
   box(subplot2,'on');
   hold(subplot2,'all');
   plot(durationWV, t_bench,'r', 'LineWidth',1, 'DisplayName','Bench T')
   plot(durationWV, surf_T, 'b', 'LineWidth',1, 'DisplayName','Surface T')
  axis([fix(min(duration)) ceil(max(duration)) 10 40]);   % -20 40]) %PRECIP
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
   plot(durationWV, (p_WVoff/7500),'b-','LineWidth', 1, 'DisplayName','P_{off}') % changed from 0.05 to 0.0425
   hold on
   plot(durationWV, p_WVon/7500,'r-', 'LineWidth',1, 'DisplayName','P_{on}')
   %plot(duration, p_hsrl/2500,'g--', 'LineWidth',1, 'DisplayName','P_{hsrl}')
   axis([fix(min(durationWV)) ceil(max(durationWV)) 0 50])
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
   plot(durationWV, surf_P, 'k-','LineWidth', 1, 'DisplayName','Surf P') 
   axis([fix(min(durationWV)) ceil(max(durationWV)) -inf inf])
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

figure(101)
semilogy(durationWV, background_comb_off)
datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
title('WV offline Background')



      
      
 %% save figure 
  cd(plot_path)
   
   FigH = figure(1);
   drawnow;
   FigH.Units = 'pixels'; % Ensure units are pixels for direct mapping to your old 'PaperPosition' width/height
   FigH.Position = [100 100 1920 250]; % x, y, width, height in pixels
   name=char(strcat(node, "_", date, '_Backscatter_Ratio_comb')); 
   exportgraphics(FigH, [name, '.png'], 'Resolution', 150);
%    
   FigH = figure(2);
   drawnow;
   FigH.Units = 'pixels'; % Ensure units are pixels for direct mapping to your old 'PaperPosition' width/height
   FigH.Position = [100 100 1920 300]; % x, y, width, height in pixels
   name=char(strcat(node, "_", date, '_Backscatter_Coeff_comb')); 
   exportgraphics(FigH, [name, '.png'], 'Resolution', 150);
  
%   FigH = figure(3);
%   set(gca,'Fontsize',16,'Fontweight','b'); 
%   set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920 250]);       
%   name=strcat(node, "_", date, "_O2_extinction_Coeff_comb");
%   print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
%   
   FigH = figure(4);
 %  drawnow;
   FigH.Units = 'pixels'; % Ensure units are pixels for direct mapping to your old 'PaperPosition' width/height
   FigH.Position = [100 100 1920 300]; % x, y, width, height in pixels
   name=char(strcat(node, "_", date, '_WV_comb')); 
 %  exportgraphics(FigH, [name, '.png'], 'Resolution', 150);


% %  
%   FigH = figure(5);
%   set(gca,'Fontsize',16,'Fontweight','b');  
%   set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920 250]);
%   name=strcat(node, "_", date, '_RB_comb'); 
%   print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
%   
   FigH = figure(6);
   drawnow;
   FigH.Units = 'pixels'; % Ensure units are pixels for direct mapping to your old 'PaperPosition' width/height
   FigH.Position = [0 0 1920 1100]; % x, y, width, height in pixels
   name=char(strcat(node, "_", date, '_all_comb')); 
  exportgraphics(FigH, [name, '.png'], 'Resolution', 150);
   % print(FigH, name, '-dpng', '-r75') % set at the screen resolution 

  
  % FigH = figure(7);
  %  FigH.Units = 'pixels'; % Ensure units are pixels for direct mapping to your old 'PaperPosition' width/height
  %  FigH.Position = [100 100 1920 550]; % x, y, width, height in pixels
  %  name=char(strcat(node, "_", date, '_house_comb')); 
  %  exportgraphics(FigH, [name, '.png'], 'Resolution', 150);


  
 cd(dd);  % point back to original directory 

 toc