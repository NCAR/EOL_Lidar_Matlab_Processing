clear all; close all;
tic

 node = 'MPD05';
 date = '18 Jun 2021';   
 days = 21; skip = 3;
 date = '10 Dec 2021';   
 days = 4; skip = 1;
%  date = '06 Jul 2021';   
%  days = 2; skip = 1;
 
 
serv_path = '/Volumes/fog1/rsfdata/MPD/';
plot_path = '/Volumes/Macintosh HD/Users/spuler/Desktop/mpd/Plots/';
C = importdata('NCAR_C_Map.mat');
dd=pwd;

if strcmp(node,'MPD01')==1
    cd(strcat(serv_path,'/mpd_01_processed_data//Matlab')) 
 elseif strcmp(node,'MPD02')==1
    cd(strcat(serv_path,'/mpd_02_processed_data//Matlab')) 
 elseif strcmp(node,'MPD03')==1
    cd(strcat(serv_path,'/mpd_03_processed_data//Matlab')) 
 elseif strcmp(node,'MPD04')==1
    cd(strcat(serv_path,'/mpd_04_processed_data//Matlab'))
 elseif strcmp(node,'MPD05')==1
    cd(strcat(serv_path,'/mpd_05_processed_data//Matlab'))
end

%% read and combine the data into a single file
for i=1:days
  if i==1  
     if exist(strcat(node, '_', datestr(date, 'yyyymmdd'), '_Backscatter Coefficient.mat'))==2
        load(strcat(node, '_', datestr(date, 'yyyymmdd'), '_Backscatter Coefficient.mat'))
         BSR_comb = BSR;
         beta_bs_comb = beta_bs;
         duration=x';
     end
  else
    date = datestr(addtodate(datenum(date), 1, 'day'), 'dd mmm yyyy');
    if exist(strcat(node, '_', datestr(date, 'yyyymmdd'), '_Backscatter Coefficient.mat'))==2
       load(strcat(node, '_', datestr(date, 'yyyymmdd'), '_Backscatter Coefficient.mat'))
       BSR_comb = vertcat(BSR_comb, BSR);
       beta_bs_comb = vertcat(beta_bs_comb, beta_bs);
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
  axis([fix(min(duration))  ceil(max(duration)) 0 12])
  caxis([1e-1 1e3]);
  %caxis([0 20]);
  hh = title({[node, ' ', ' Backscatter Ratio']},'fontweight','b','fontsize',font_size);  
  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size);
  set(gca,'Fontsize',font_size,'Fontweight','b');
  set(gca,'Zscale', 'log')
  set(gca,'Colorscale', 'log')
  set(gca,'Zscale', 'linear')
  colormap(jet)
 % colormap(flipud(hot))
  
  
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
  axis([fix(min(duration))  ceil(max(duration)) 0 12])
  caxis([1e-7 1e-3]);
  hh = title({[node, ' ', ' Backscatter Coefficient [m^{-1}sr^{-1}]']},'fontweight','b','fontsize',font_size);  
  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
  xlabel('Time (UTC)','fontweight','b','fontsize',font_size);
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size);
  set(gca,'Fontsize',font_size,'Fontweight','b');
  set(gca,'Zscale', 'log')
  set(gca,'Colorscale', 'log')
  set(gca,'Zscale', 'linear')
  colormap(jet)

 %% save figure 
  cd(plot_path)
   
  FigH = figure(1);
  set(gca,'Fontsize',16,'Fontweight','b'); 
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1500 300]);
  name=strcat(node, "_", date, "_Backscatter_Ratio_comb");
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
   
  FigH = figure(2);
  set(gca,'Fontsize',16,'Fontweight','b'); 
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1500 300]);
  name=strcat(node, "_", date, "_Backscatter_Coeff_comb");
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
  
 cd(dd);  % point back to original directory 

 toc