addpath('./natsortfiles')
clear all; close all

serv_path1 = '/Volumes/Macintosh HD/Users/spuler/Desktop/mpd/';  % need to use finder and connect to server 'smb://cit.ucar.edu/eol'
%serv_path = '/Volumes/fog1/rsfdata/MPD/mpd_field_campaigns/2022_PRECIP/ProcessedData/';
%serv_path = '/Volumes/eol/sci/mhayman/DIAL/Processed_Data/PRECIP/low_alt/qc_masked';
%serv_path = '/Volumes/eol/sci/mhayman/DIAL/Processed_Data/PRECIP/const_smooth/';
%serv_path = '/Volumes/eol/sci/mhayman/DIAL/Processed_Data/PRECIP/no_smooth/';
%serv_path = '/Volumes/eol/sci/mhayman/DIAL/Processed_Data/PRECIP/10min_37m/qc_masked/';
%serv_path = '/Volumes/eol/sci/mhayman/DIAL/Processed_Data/PRECIP/final_adj_mask/qc_masked/';
%serv_path = '/Volumes/eol/sci/mhayman/DIAL/Processed_Data/PRECIP/full_altitude/qc_masked/';
serv_path = '/Volumes/eol/sci/mhayman/DIAL/Processed_Data/PRECIP/final_adj_mask_iii/qc_masked/';
serv_path = '/Volumes/eol/sci/mhayman/DIAL/Processed_Data/PRECIP/full_alt_mask_iii/qc_masked/';
serv_path = '/Volumes/eol/sci/mhayman/DIAL/Processed_Data/ECLIPSE/test_1.2/qc_masked/';
serv_path = '/Volumes/eol/sci/mhayman/DIAL/Processed_Data/ECLIPSE/test_1.2/qc_masked/';
serv_path = '/Volumes/eol/sci/mhayman/DIAL/Processed_Data/NOAA2023/1.0/';
plot_path = '/Volumes/Macintosh HD/Users/spuler/Desktop/mpd/Plots/';

node = 'mpd05';
cd(strcat(serv_path))
addpath /Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing
d_read_data = pwd; % get the current path
cd(strcat(serv_path1,'/Plots'))
d_save_data = pwd; % set the plot save path
flag.save_data = 1;  % save data at end of processing (0=off 1=on)
flag.decimate = 0;  % decimate time dimension to fit screen (needed for long timespans)
flag.grid75 = 0; % grid data to 75 m
skip = 3
mask_var = sqrt(25); %standard is 25
force_low_range_mask_off = 0;  %removes mask from low range  

cd(d_read_data);
[Pythonfilename, Pythondir] = uigetfile(strcat(node,'*.*'),'Select the sonde file', 'MultiSelect', 'on');
serv_path = Pythondir;
jj=1;
% sort file names into order (https://fr.mathworks.com/matlabcentral/fileexchange/47434-natural-order-filename-sort)
Pythonfilename = natsortfiles(Pythonfilename); 
% variable{1} = 'Aerosol_Backscatter_Coefficient_variance';
 variable{1} = 'time';
 variable{2} = 'range';
 variable{3} = 'Absolute_Humidity'; 
 variable{4} = 'Absolute_Humidity_mask';
 variable{5} = 'Absolute_Humidity_variance'; 
% variable{6} = 'Backscatter_Ratio';
% variable{7} = 'Backscatter_Ratio_mask';
% variable{8} = 'Backscatter_Ratio_variance'; 
  variable{6} = 'Aerosol_Backscatter_Coefficient';
  variable{7} = 'Aerosol_Backscatter_Coefficient_mask';
  variable{8} = 'Aerosol_Backscatter_Coefficient_variance';
% variable{6} = 'Relative_Backscatter';
 variable{10} = 'Surface_AbsHum';
 variable{11} = 'Absolute_Humidity_mask_layers';
    
for jj = 1:size(Pythonfilename,2)
  try
     filename = Pythonfilename{jj};
  catch
    filename = Pythonfilename;
  end
  file_date = filename(7:14); % changed to full year
  n = datenum(file_date, 'yyyymmdd');
  n_test{jj} = datenum(file_date, 'yyyymmdd');
  ncid = netcdf.open(filename, 'NC_NOWRITE');
    %ncdisp(filename, '/', 'min') % use this to display all variables
    %ncdisp(filename) % use this to display all variables
    time{jj} = ncread(filename,variable{1}); 
    alt{jj} = ncread(filename,variable{2});
    AH{jj}  = ncread(filename,variable{3});  
    AH_mask{jj} = ncread(filename,variable{4}); 
    AH_var{jj} = sqrt(ncread(filename,variable{5}));
    AH_sur{jj} = ncread(filename,variable{10});  
    AH_mask_layers{jj} = ncread(filename,variable{11}); %each mask layer applied to the Absolute Humidity data.  The order is reflective of when they are applied.'
    
%     if force_low_range_mask_off == 1
%       AH_mask{jj}(6:9,:)= 0; %as a test turn off the low range mask
%     end
    AH{jj}(AH_mask{jj} == 1) = nan;
    AH{jj}(AH_var{jj} >= mask_var) = nan;
    AH_var{jj}(AH_mask{jj} == 1) = nan; 
    AH{jj}(isnan(AH_var{jj})) = nan;

    %     AH_mask_layer_1 = squeeze(AH_mask_layers{jj}(1,:,:));
%     AH_mask_layer_2 = squeeze(AH_mask_layers{jj}(2,:,:));
%     AH_mask_layer_3 = squeeze(AH_mask_layers{jj}(3,:,:));
%     AH_mask_layer_4 = squeeze(AH_mask_layers{jj}(4,:,:));
%     AH_mask_layer_5 = squeeze(AH_mask_layers{jj}(5,:,:));
%     AH_mask_all = AH_mask{jj};
%    AH_mask_layers_all = zeros(size(AH_mask_all));
%    AH_mask_layers_all((AH_mask_layer_1 == 0) & (AH_mask_layer_2 == 0) & (AH_mask_layer_4 == 0)) = 1;
%    AH{jj}(AH_mask_layers_all==0) = nan;  %mask layers
%   %  mask the absolute humidity data based on its variance
%   AH{jj}(AH_var{jj} >= mask_var) = nan;  %standard is 25 (5) 
%   %  AH_var{jj}(AH_var{jj} >= mask_var) = nan; 


% figure(10)
% mesh(AH_mask_all)
% figure(11)
% mesh(AH_mask_layer_1)
% figure(12)
% mesh(AH_mask_layer_2)
% figure(13)
% mesh(AH_mask_layer_3)
% figure(14)
% mesh(AH_mask_layer_4)
% figure(15)
% mesh(AH_mask_layers_all)


    
  ABC{jj}  = ncread(filename,variable{6});    

   %     try
%       ABC_mask{jj} = ncread(filename,variable{7}); 
%       ABC_var{jj} = ncread(filename,variable{8});
%       
%       ABC{jj}(ABC{jj} < sqrt(ABC_var{jj})) = nan;   
%   %    ABC{jj}(ABC_mask{jj} == 1) = nan;
%   %    ABC_var{jj}(ABC_mask{jj} == 1) = nan; 
%     catch
% %      ABC_mask{jj} = ncread(filename,variable{4}); 
% %      ABC_var{jj} = ncread(filename,variable{5}); 
%     end
    
    
  netcdf.close(ncid); 
  %convert from Unix time to date number (days since Jan 0 0000) 
  duration{jj} =  n+double(time{jj}/3600/24);

  %range_grid = (0:75:12000);
  %AH_grid{jj} = interp1( AH{jj}(:,1) ,AH{jj}(:,2:end), range_grid, 'nearest', 'extrap'); 
 % for jj=1:12  (there are problems with the data for the first three days)

  % mask the lowest gates if desired
 % low_mask = repmat(alt{jj}, 1, size(AH{jj}, 2)); 
 % AH{jj}(low_mask < low_range_mask)= nan;
 % AH_var{jj}(low_mask < low_range_mask)= nan;
 
  if jj == 1
      comb_duration = duration{jj};
      comb_AH = AH{jj}';
      comb_AH_var = AH_var{jj}';
      comb_ABC = ABC{jj}';
      comb_AH_sur = AH_sur{jj};
%       comb_ABC_var = ABC_var{jj}';
  else
      comb_duration = [comb_duration; duration{jj}];
      % find the maximum range to accumulate (in case it changes)
      max_range = min(cellfun('size',alt,1))
      comb_AH = [comb_AH(:,1:max_range); AH{jj}(1:max_range,:)'];
      comb_AH_var = [comb_AH_var(:,1:max_range); AH_var{jj}(1:max_range,:)'];
      comb_ABC = [comb_ABC(:,1:max_range); ABC{jj}(1:max_range,:)'];
      comb_AH_sur = vertcat(comb_AH_sur, AH_sur{jj});
%       comb_ABC_var = [comb_ABC_var(:,1:max_range); ABC_var{jj}(1:max_range,:)'];
  end
  
  % remove non-physical values AH < - 0.5 (allow for some slight negative noise)
  comb_AH(comb_AH < -0.5)= nan;
  
 % end
  
end

scrsz = get(0,'ScreenSize');
date=datestr(n, 'yyyy-mmm-dd');
plot_size1 = [scrsz(4)/1.5 scrsz(4)/10 scrsz(3)/1.5 scrsz(4)/3];
font_size = 14;

if flag.grid75 == 1
% grid to 75m AGL
 r_grid = (0:75:alt{1}(end));
 range = alt{1}';
    comb_AH = interp1(range, comb_AH', r_grid , 'linear', 'extrap')';
    comb_AH_var = interp1(range, comb_AH_var', r_grid , 'linear', 'extrap')'; 
    comb_ABC = interp1(range, comb_ABC', r_grid , 'linear', 'extrap')'; 
    alt{1} = r_grid';
end

x = comb_duration;
y = alt{1}/1000;
Z_AH = real(comb_AH)';
Z_AH(2,:) = comb_AH_sur';  % add in the AbsHum weather station data (this data starts below the ground so need to use bin 2 and 3)
Z_AH(3,:) = comb_AH_sur';  % add in the AbsHum weather station data
Z_ABC = real(comb_ABC)';


flag.int = 0; % interpolate nans in nanmoving_average
if flag.decimate == 1 
  decimate_time = fix(length(comb_duration)/scrsz(3)/2); %decimate to number of horiz pixels;
  decimate_range = 1; % keep native gate spacing 
  % average data before decimating
   Z_AH = nanmoving_average(Z_AH,decimate_time,2,flag.int);
   Z_ABCH = nanmoving_average(Z_ABC,decimate_time,2,flag.int);
  % then decimate in time
   x = x(1:decimate_time:end);
   Z_AH =  Z_AH(1:decimate_range:end, 1:decimate_time:end);
   Z_ABC =  Z_ABC(1:decimate_range:end, 1:decimate_time:end);
end
    
xData =  linspace( fix(min(x)),  ceil(max(x)), round((ceil(max(x))-fix(min(x)))/skip)+1 );

% plot the AH
  figure1 = figure('Position',plot_size1);
  set(gcf,'renderer','zbuffer');
  h = pcolor(x, y, Z_AH);
  set(h, 'EdgeColor', 'none'); 
  axis xy; colorbar('EastOutside'); 
  caxis([0 15]);
  %caxis([0 6]);
  axis([fix(min(x)) ceil(max(x)) 0 6]) 
%  shading interp
  set(gca, 'XTick',  xData)
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  hh = title({[node, ' Absolute Humidity (g m^{-3})']},...
       'fontweight','b','fontsize',font_size);  
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
  colormap(jet)
  set(gca,'Fontsize',font_size,'Fontweight','b');
 
   
%  % plot the atmospheric backscatter coefficient 

  figure2 = figure('Position',plot_size1);
  set(gcf,'renderer','zbuffer');
  h = pcolor(x, y, Z_ABC);
  set(h, 'EdgeColor', 'none'); 
  axis xy; colorbar('EastOutside'); 
  caxis([0 12]);
  caxis([1e-8 1e-3]);
%  caxis([1e3 1e8]);  % for relative backscatter
%  caxis([1e-1 1e4]); % backscatter ratio
  axis([fix(min(x)) ceil(max(x)) 0 18.25])
  axis([fix(min(x)) ceil(max(x)) 0 12])
%  axis([fix(min(x)) ceil(max(x)) 0 6])
  %  shading interp
  set(gca, 'XTick',  xData)
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  set(gca,'Zscale', 'log')
  set(gca,'Colorscale', 'log')
  set(gca,'Zscale', 'linear')
  hh = title({[node, ' Attenuated Backscatter ']},...
       'fontweight','b','fontsize',font_size);    
%  hh = title({[node, ' Aerosol Backscatter Coefficient m^{-1} sr^{-1}']},...
%       'fontweight','b','fontsize',font_size);     
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
%  colormap(jet)
  set(gca,'Fontsize',font_size,'Fontweight','b');
  
% %   plot the backcatter ratio 
%   Z = real(comb_ABC)';
%   figure2 = figure('Position',plot_size1);
%   set(gcf,'renderer','zbuffer');
%   h = pcolor(x, y, Z);
%   set(h, 'EdgeColor', 'none'); 
%   axis xy; colorbar('EastOutside'); 
%   caxis([0 12]);
%   caxis([1e-1 1e3]);
%   axis([fix(min(x)) ceil(max(x)) 0 8])
%   %  shading interp
%   set(gca, 'XTick',  xData)
%   set(gca,'TickDir','out');
%   set(gca,'TickLength',[0.005; 0.0025]);
%   set(gca,'Zscale', 'log')
%   set(gca,'Colorscale', 'log')
%   set(gca,'Zscale', 'linear')
%   hh = title({[node, ' Backscatter Ratio']},...
%        'fontweight','b','fontsize',font_size);     
%   ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
%   datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
%   colormap(jet)
%   colormap(flipud(hot))
%   set(gca,'Fontsize',font_size,'Fontweight','b');
 
 
%   [minValue, closestIndex] = min(abs(comb_duration - datenum('21-Jul-2021 10:00:00')))
%   figure(20)
%   BSR_profile = Z(:,closestIndex); 
%   plot(BSR_profile, y);
%   save ('BSR_at_21Jun2021_1000','BSR_profile', 'y') 


% make AH with slider 
 skip = 1
 dx = 3;  % # days to plot at time
 xData =  linspace( fix(min(x)),  ceil(max(x)), round((ceil(max(x))-fix(min(x)))/skip)+1 ); 
 plot_size2 = [scrsz(4) scrsz(4) scrsz(3) scrsz(4)/3];
 figure3 = figure('Position',plot_size2);
 set(gcf,'renderer','zbuffer');
 h = pcolor(x, y, Z_AH);
 set(h, 'EdgeColor', 'none'); 
 axis xy; colorbar('EastOutside'); 
 caxis([0 15]);
  %caxis([0 6]);
 axis([fix(min(x)) ceil(max(x)) 0 6]) 
%  shading interp
 set(gca, 'XTick',  xData)
 set(gca,'TickDir','out');
 set(gca,'TickLength',[0.005; 0.0025]);
 hh = title({[node, ' Absolute Humidity (g m^{-3})']},...
       'fontweight','b','fontsize',font_size);  
 ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
 datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
 colormap(jet)
 set(gca,'Fontsize',font_size,'Fontweight','b');
 grid on
 grid minor
 a=gca;
 set(gcf, 'doublebuffer','on');
 set(a,'xlim',[fix(min(x)) fix(min(x))+dx]);
 %%%%%Generate constants for use in uicontrol initialization
 pos=get(a,'position');
 Newpos=[pos(1) pos(2)-0.1 pos(3) 0.03];    
 S=['set(gca,''xlim'',get(gcbo,''value'')+[0 ' num2str(dx) '])'];
 h=uicontrol('style','slider',...
      'units','normalized','position',Newpos,...
      'value', fix(min(x)),...
      'callback',S,'min',fix(min(x)),'max',ceil(max(x))-dx);  



  
 if flag.save_data == 1
  
  cd(d_save_data);     
  range = alt{1}'; 
  N_avg_comb = (comb_AH./1e6.*6.022E23./18.015);
  duration = duration;

  cd(plot_path); 
  FigH = figure(1);
  set(gca,'Fontsize',16,'Fontweight','b'); 
%  set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1500 300]);
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920 250]);
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920 185]);
  name=strcat(date, node, ' WV_Python_multi'); 
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
  
  FigH = figure(2);
  set(gca,'Fontsize',16,'Fontweight','b'); 
 % set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1500 300]);
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920 250]);  
  name=strcat(date, node, ' RB_Python_multi'); 
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution
  
  
  %cd('/Users/spuler/Desktop/WV_DIAL_data') % point to the directory where data is stored 
  serv_path = '/Volumes/eol/smaug1/rsfdata/MPD/';
  name=strcat(date, '_combined');
  if strcmp(node,'mpd01')==1
     cd(strcat(serv_path, 'mpd_01_processed_data/Matlab')) 
      MPD01.N_avg_comb = N_avg_comb;
     % MPD01.RB_comb = RB_comb;
      MPD01.range = y;
      MPD01.time = x;
      save(name, 'MPD01')
  elseif strcmp(node,'mpd02')==1
      cd(strcat(serv_path, 'mpd_02_processed_data/Matlab')) 
      MPD02.N_avg_comb = N_avg_comb;
    %  MPD02.RB_comb = RB_comb;
      MPD02.range = y;
      MPD02.time = x;
      save(name, 'MPD02')
  elseif strcmp(node,'mpd03')==1
      cd(strcat(serv_path, 'mpd_03_processed_data/Matlab')) 
      MPD03.N_avg_comb = N_avg_comb;
    %  MPD03.RB_comb = RB_comb;
      MPD03.range = y;
      MPD03.time = x;
      save(name, 'MPD03')
  elseif strcmp(node,'mpd04')==1
      cd(strcat(serv_path, 'mpd_04_processed_data/Matlab')) 
      MPD04.N_avg_comb = N_avg_comb;
    %  MPD04.RB_comb = RB_comb;
      MPD04.range = y;
      MPD04.time = x;
      save(name, 'MPD04')
   elseif strcmp(node,'mpd05')==1
       cd(strcat(serv_path, 'mpd_05_processed_data/Matlab')) 
       MPD05.N_avg_comb = N_avg_comb;
   %   MPD05.RB_comb = RB_comb;
      MPD05.range = y;
      MPD05.time =  x;
      save(name, 'MPD05')
  end

 end
  