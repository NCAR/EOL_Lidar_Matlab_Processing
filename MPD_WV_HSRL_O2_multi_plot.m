clear all; close all;
tic

 node = 'MPD05';
 date = '26 Oct 2023';   
 days = 7; skip = 1;
 flag.afterpulse = 1; % read in the afterpulse corrected data (0=off 1=on)
 WV_max_scale = 10;
 
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
        p_on = P_on;
        p_off = P_off;
        t_bench = T_bench;
     end 
  else
    date = datestr(addtodate(datenum(date), 1, 'day'), 'dd mmm yyyy');
    if exist(strcat(node, '_', datestr(date, 'yyyymmdd'), '_Backscatter Coefficient.mat'))==2
       load(strcat(node, '_', datestr(date, 'yyyymmdd'), '_Backscatter Coefficient.mat'))
       BSR_comb = vertcat(BSR_comb, BSR);
       beta_bs_comb = vertcat(beta_bs_comb, beta_bs);
       alpha_O2_comb = vertcat(alpha_O2_comb, alpha_O2);
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
        p_on = vertcat(p_on, P_on(2:end,:));  
        p_off = vertcat(p_off, P_off(2:end,:)); 
        t_bench = vertcat(t_bench, T_bench(2:end,:));  
    end
    
  end
end



%% plot data
cd(dd);
xData =  linspace( fix(min(duration)),  ceil(max(duration)), round((ceil(max(duration))-fix(min(duration)))/skip)+1 );

  figure(1)
  x = (duration)';
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
  caxis([1 10]);
  hh = title({[node, ' ', ' Backscatter Ratio']},'fontweight','b','fontsize',font_size);  
  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size);
  set(gca,'Fontsize',font_size,'Fontweight','b');
  set(gca,'Zscale', 'log')
  set(gca,'Colorscale', 'log')
  set(gca,'Zscale', 'linear')
  colormap(jet)
  colormap(viridis)
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
  caxis([1e-8 1e-6]);
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
      x = (duration)';
      Z = real(double((((BSR_comb))')));
      set(gcf,'renderer','zbuffer');
      h = pcolor(x,y,Z);
      set(h, 'EdgeColor', 'none');
      set(gca,'TickDir','out');
      set(gca,'TickLength',[0.005; 0.0025]);
      set(gca, 'XTick',  xData) 
      colorbar('EastOutside');
      axis([fix(min(duration))  ceil(max(duration)) 0 6])
      caxis([1 10]);
      hh = title({[node, ' ', ' Backscatter Ratio']},'fontweight','b','fontsize',font_size);  
      datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
      ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size);
      set(gca,'Fontsize',font_size,'Fontweight','b');
      set(gca,'Zscale', 'log')
      set(gca,'Colorscale', 'log')
      set(gca,'Zscale', 'linear')
      colormap(subplot1, viridis)
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
      x = (duration)';    
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
      
  
 %% save figure 
  cd(plot_path)
   
%   FigH = figure(1);
%   set(gca,'Fontsize',16,'Fontweight','b'); 
%   set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920 250]);      
%   name=strcat(node, "_", date, "_Backscatter_Ratio_comb");
%   print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
%    
% %   FigH = figure(2);
% %   set(gca,'Fontsize',16,'Fontweight','b'); 
% %   set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920 250]);      
% %   name=strcat(node, "_", date, "_Backscatter_Coeff_comb");
% %   print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
% %   
%   FigH = figure(3);
%   set(gca,'Fontsize',16,'Fontweight','b'); 
%   set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920 250]);       
%   name=strcat(node, "_", date, "_O2_extinction_Coeff_comb");
%   print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
%   
%   FigH = figure(4);
%   set(gca,'Fontsize',16,'Fontweight','b');  
%   set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920 250]);
%   name=strcat(node, "_", date, '_WV_comb'); 
%   print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
%  
%   FigH = figure(5);
%   set(gca,'Fontsize',16,'Fontweight','b');  
%   set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920 250]);
%   name=strcat(node, "_", date, '_RB_comb'); 
%   print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
%   
  FigH = figure(6);
  %set(gca,'Fontsize',16,'Fontweight','b');  
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920 275*4]);
  name=strcat(node, "_", date, '_all_comb'); 
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
  
  
 cd(dd);  % point back to original directory 

 toc