%cd /scr/sci/spuler/mpd/sgp/raman_lidar
clear all; close all
%serv_path = '/Volumes/eol/sci/spuler';
serv_path = '/Users/spuler/Desktop';
%serv_path = '/Volumes/eol/fog1/rsfdata/MPD';
%cd(strcat(serv_path,'/mpd/mpd_05_processed_data/python'))
cd(strcat(serv_path,'/SOCRATES_data'))
%serv_path = '/Volumes/eol/sci/mhayman/DIAL/Processed_Data';
%cd(strcat(serv_path,'/mpd/Marshall/mpd03_5min'))
%cd(strcat(serv_path, '/MPDSGP/t_res_10min/'));  %This was the original SGP 10 min data
%cd(strcat(serv_path,'/MPD_Gen5_Pub/'));
d_read_data = pwd; % get the current path

%serv_path = '/Volumes/eol/fog1/rsfdata/MPD';
cd(strcat(serv_path,'/mpd/Plots'))
d_save_data = pwd; %set the plot save path
flag.save_plots = 1;  %save plots (0=off 1=on)
node = 'GV-HSRL';
low_range_mask = 0;

 plot_start_time =  '19-Jan-2018 04:20:00';
 plot_end_time = '19-Jan-2018 04:35:00';

cd(d_read_data);
[Pythonfilename, Pythondir] = uigetfile('*.*','Select the sonde file', 'MultiSelect', 'on');
jj=1;

 %variable{1} = 'Aerosol_Backscatter_Coefficient_variance';
  
 variable{1} = 'time';
 variable{2} = 'range';

 variable{3} = 'Volume_Linear_Depolarization_Ratio'; 
 variable{4} = 'Volume_Linear_Depolarization_Ratio_mask';
 variable{5} = 'Volume_Linear_Depolarization_Ratio_variance'; 
 variable{6} = 'Aerosol_Backscatter_Coefficient';
 variable{7} = 'Aerosol_Backscatter_Coefficient_mask';
 variable{8} = 'Aerosol_Backscatter_Coefficient_variance';
 
 variable{9} = 'Molecular_Backscatter_Channel'; 
 variable{10} = 'Merged_Combined_Channel';   
 variable{11} = 'Cross_Polarization_Channel';  


for jj = 1:size(Pythonfilename,2)
  filename = Pythonfilename{jj};
  %date = filename(end-15:end-10);
  %n = datenum(date, 'yymmdd');
  date = filename(end-15:end-8);
  n = datenum(date, 'yyyymmdd');
  
  ncid = netcdf.open(filename, 'NC_NOWRITE');
    %ncdisp(filename, '/', 'min') % use this to display all variables
    ncdisp(filename) % use this to display all variables
    time{jj} = ncread(filename,variable{1}); 
    alt{jj} = ncread(filename,variable{2});
    aircraft_alt{jj} = ncread(filename,'GGALT'); 
    pointing{jj} = ncread(filename,'TelescopeDirection');
 
    ABC{jj}  = ncread(filename,variable{6});   
    ABC_mask{jj} = ncread(filename,variable{7}); 
    ABC_var{jj} = ncread(filename,variable{8});
    ABC{jj}(ABC_mask{jj} == 1) = nan;
    ABC_var{jj}(ABC_mask{jj} == 1) = nan; 
    
    
    LDR{jj}  = ncread(filename,variable{3});  
    LDR_mask{jj} = ncread(filename,variable{4}); 
    LDR_var{jj} = ncread(filename,variable{5}); 
    LDR{jj}(LDR_mask{jj} == 1) = nan;
    LDR_var{jj}(LDR_mask{jj} == 1) = nan; 
 
     % add a QC step
     mol_counts{jj} = ncread(filename,variable{9}); 
     comb_counts{jj} = ncread(filename,variable{10});
     cross_counts{jj} = ncread(filename,variable{11});
    % LDR{jj}(ABC{jj} <= 1e-9) = nan;
    % ABC{jj}(ABC{jj} <= 1e-9) = nan; % remove unrealisitc backscatter values
      LDR{jj}(mol_counts{jj} < 1) = nan;
      ABC{jj}(mol_counts{jj} < 1) = nan; % remove data with < photon counted
%      LDR{jj}(comb_counts{jj} < 10) = nan;
%      ABC{jj}(comb_counts{jj} < 10) = nan; % remove data with < 1 photon counted
%      LDR{jj}(cross_counts{jj} < 10) = nan;
%      ABC{jj}(cross_counts{jj} < 10) = nan; % remove data with < 1 photon counted
     
     
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
      comb_aircraft_alt = aircraft_alt{jj};
      comb_pointing = pointing{jj};
       comb_LDR = LDR{jj}';
       comb_LDR_var = LDR_var{jj}';
      comb_ABC = ABC{jj}';
      comb_ABC_var = ABC_var{jj}';
  else
      comb_duration = [comb_duration; duration{jj}];
      comb_aircraft_alt = [comb_aircraft_alt; aircraft_alt{jj}];
      comb_pointing = [comb_pointing; pointing{jj}];
      % find the maximum range to accumulate (in case it changes)
      max_range = min(cellfun('size',alt,1))
      comb_LDR = [comb_LDR(:,1:max_range); LDR{jj}(1:max_range,:)'];
      comb_LDR_var = [comb_LDR_var(:,1:max_range); LDR_var{jj}(1:max_range,:)'];
      comb_ABC = [comb_ABC(:,1:max_range); ABC{jj}(1:max_range,:)'];
      comb_ABC_var = [comb_ABC_var(:,1:max_range); ABC_var{jj}(1:max_range,:)'];
  end
  
  % remove non-physical values AH < - 0.5 (allow for some slight negative noise)
%   comb_AH(comb_AH < -0.5)= nan;
  
 % end
  
end

scrsz = get(0,'ScreenSize');
date=datestr(n, 'yyyy-mmm-dd');
plot_size1 = [scrsz(4)/1.5 scrsz(4)/10 scrsz(3)/1.5 scrsz(4)/3];
font_size = 14;

% change from range to aircraft coordinates
lidar_range = (alt{1}(1:max_range,1));
lidar_range_all = repmat(lidar_range, 1, length(comb_aircraft_alt));
air_alt = repmat(comb_aircraft_alt, 1, length(lidar_range));

%(comb_pointing==0)
global_range = air_alt'+lidar_range_all;
for k=1:length(comb_pointing)
    if comb_pointing(k) == 0
        global_range(:,k) = air_alt(k,:)'-lidar_range_all(:,k);
    end
end


x = comb_duration;
%y = double(alt{1}(1:max_range,1)')/1000;
y = global_range./1000;
% Z = real(comb_AH)';
skip = 1
xData =  linspace( fix(min(x)),  ceil(max(x)), round((ceil(max(x))-fix(min(x)))/skip)+1 );


 
 % plot the atmospheric backscatter coefficient 
  Z = real(comb_ABC)';
  figure1 = figure('Position',plot_size1);
  set(gcf,'renderer','zbuffer');
  h = pcolor(x, y, Z);
  set(h, 'EdgeColor', 'none');  axis xy; colorbar('EastOutside'); 
  caxis([1e-8 1e-3]);
  %axis([fix(min(x)) ceil(max(x)) 0 6])
  %  shading interp
  set(gca, 'XTick',  xData); set(gca,'TickDir','out'); set(gca,'TickLength',[0.005; 0.0025]);
  set(gca,'Zscale', 'log'); set(gca,'Colorscale', 'log'); set(gca,'Zscale', 'linear')
  hh = title({[date, ' Aerosol Backscatter Coefficient m^{-1} sr^{-1}']},...
       'fontweight','b','fontsize',font_size);     
  ylabel('Height (km, MSL)','fontweight','b','fontsize',font_size); 
 % datetick('x','HH:MM','keeplimits', 'keepticks');
  datetick('x','HH:MM');
  colormap(jet)
  set(gca,'Fontsize',font_size,'Fontweight','b');
  ylim([-.1 2]);
  xlim([datenum(plot_start_time) datenum(plot_end_time)]);
 

 
% % plot the Linear Depolarization
   figure2 = figure('Position',plot_size1);
   Z = real(comb_LDR)';
   set(gcf,'renderer','zbuffer');
   h = pcolor(x, y, Z);
   set(h, 'EdgeColor', 'none');  axis xy; colorbar('EastOutside'); 
   caxis([0 1]);
%   axis([fix(min(x)) ceil(max(x)) 0 6]) 
% %  shading interp
   set(gca, 'XTick',  xData); set(gca,'TickDir','out'); set(gca,'TickLength',[0.005; 0.0025]);
%    caxis([0.1 1]);
%    set(gca,'Zscale', 'log'); set(gca,'Colorscale', 'log'); set(gca,'Zscale', 'linear')
   hh = title({[date, ' Linear Depolarization Ratio']},...
        'fontweight','b','fontsize',font_size);  
   ylabel('Height (km, MSL)','fontweight','b','fontsize',font_size); 
   datetick('x','HH:MM');
   colormap(jet)
   set(gca,'Fontsize',font_size,'Fontweight','b');
  ylim([-.1 2]);
  xlim([datenum(plot_start_time) datenum(plot_end_time)]);
 
 
 if flag.save_plots == 1
  
%   cd(d_save_data);     
%   range = alt{1}'; 
%   N_avg_comb = (comb_AH./1e6.*6.022E23./18.015);
%   duration = duration;

  FigH = figure(1);
  set(gca,'Fontsize',16,'Fontweight','b'); 
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920 400]);
  name=strcat(date, node, ' GVHSRL_backscatter_coeff'); 
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
  
  FigH = figure(2);
  set(gca,'Fontsize',16,'Fontweight','b'); 
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920 400]);
  name=strcat(date, node, ' GVHSRL_Linear_depol_ratio'); 
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution
  
  
%   %cd('/Users/spuler/Desktop/WV_DIAL_data') % point to the directory where data is stored 
%   name=strcat(date, '_combined');
%   if strcmp(node,'MPD01')==1
%       MPD01.N_avg_comb = N_avg_comb;
%      % MPD01.RB_comb = RB_comb;
%       MPD01.range = range;
%       MPD01.time = duration;
%       save(name, 'MPD01')
%   elseif strcmp(node,'MPD02')==1
%       MPD02.N_avg_comb = N_avg_comb;
%     %  MPD02.RB_comb = RB_comb;
%       MPD02.range = range;
%       MPD02.time = duration;
%       save(name, 'MPD02')
%   elseif strcmp(node,'MPD03')==1
%       MPD03.N_avg_comb = N_avg_comb;
%     %  MPD03.RB_comb = RB_comb;
%       MPD03.range = range;
%       MPD03.time = duration;
%       save(name, 'MPD03')
%   elseif strcmp(node,'MPD04')==1
%       MPD04.N_avg_comb = N_avg_comb;
%     %  MPD04.RB_comb = RB_comb;
%       MPD04.range = range;
%       MPD04.time = duration;
%       save(name, 'MPD04')
%    elseif strcmp(node,'MPD05')==1
%       MPD05.N_avg_comb = N_avg_comb;
%    %   MPD05.RB_comb = RB_comb;
%       MPD05.range = range;
%       MPD05.time =  duration;
%       save(name, 'MPD05')
%   end

 end
  