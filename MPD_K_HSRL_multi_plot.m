clear all; close all;
tic

 node = 'MPD05';
 date = '7 Apr 2024';   
 days = 4; skip = 1;
% date = '29 Aug 2023';   
% days = 3; skip = 1;
%      date = '18 Jul 2023';
%      days = 70; skip = 5;
 flag.afterpulse = 1; % read in the afterpulse corrected data (0=off 1=on)
 
 
%serv_path = '/Volumes/fog1/rsfdata/MPD/';
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
  else
    date = datestr(addtodate(datenum(date), 1, 'day'), 'dd mmm yyyy');
    if exist(strcat(node, '_', datestr(date, 'yyyymmdd'), '_Backscatter Coefficient.mat'))==2
       load(strcat(node, '_', datestr(date, 'yyyymmdd'), '_Backscatter Coefficient.mat'))
       BSR_comb = vertcat(BSR_comb, BSR);
       beta_bs_comb = vertcat(beta_bs_comb, beta_bs);
       alpha_O2_comb = vertcat(alpha_O2_comb, alpha_O2);
       duration = vertcat(duration, x');
    end
  end
end


%% plot data

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
  %caxis([1 300]);
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

%   colormap(flipud(hot))
%   colormap(hot)
%   colormap(parula)
%   caxis([1e-0 1e1]);
%   ylim([0 3])
%   xlim([datenum('16-Jul-2021 17:00', 'dd-mmm-yyyy HH:MM') datenum('16-Jul-2021 19:00', 'dd-mmm-yyyy HH:MM')])
%   xData =  linspace(datenum('16-Jul-2021 17:00', 'dd-mmm-yyyy HH'), datenum('16-Jul-2021 19:00', 'dd-mmm-yyyy HH'), 3);
%   set(gca, 'XTick',  xData) 
%   datetick('x','dd-mmm-yy HH:MM','keeplimits', 'keepticks');
  
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
  caxis([1e-8 1e-5]);
  hh = title({[node, ' ', ' Backscatter Coefficient [m^{-1}sr^{-1}]']},'fontweight','b','fontsize',font_size);  
  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
%  xlabel('Time (UTC)','fontweight','b','fontsize',font_size);
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
%  xlabel('Time (UTC)','fontweight','b','fontsize',font_size);
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size);
  set(gca,'Fontsize',font_size,'Fontweight','b');
%   set(gca,'Zscale', 'log')
%   set(gca,'Colorscale', 'log')
%   set(gca,'Zscale', 'linear')
  %colormap(jet)
  
  
 %% save figure 
  cd(plot_path)
   
  FigH = figure(1);
  set(gca,'Fontsize',16,'Fontweight','b'); 
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920 250]);      % 1500 300
%   set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 750 250]);       % 1500 300
  name=strcat(node, "_", date, "_Backscatter_Ratio_comb");
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
   
   FigH = figure(2);
   set(gca,'Fontsize',16,'Fontweight','b'); 
   set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920 250]);      % 1500 300
   name=strcat(node, "_", date, "_Backscatter_Coeff_comb");
   print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
%   
  FigH = figure(3);
  set(gca,'Fontsize',16,'Fontweight','b'); 
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920 250]);       % 1500 300
  name=strcat(node, "_", date, "_O2_extinction_Coeff_comb");
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
  
 cd(dd);  % point back to original directory 

 toc