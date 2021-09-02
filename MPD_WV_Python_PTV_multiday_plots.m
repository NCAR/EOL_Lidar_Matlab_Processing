%cd /scr/sci/spuler/mpd/sgp/raman_lidar
clear all; close all
%serv_path = '/Volumes/eol/sci/spuler';
%serv_path = '/Volumes/documents/MPD';
serv_path = '/Volumes/eol/fog1/rsfdata/MPD/';  % need to use finder and connect to server 'smb://cit.ucar.edu/eol'
% cd(strcat(serv_path,'/mpd_05_processed_data/python'))
cd(strcat(serv_path,'/mpd_05_processed_data/PTV'))
d_read_data = pwd; % get the current path

serv_path = '/Volumes/Macintosh HD/Users/spuler/Desktop/mpd';
cd(strcat(serv_path,'/Plots'))
d_save_data = pwd; %set the plot save path
flag.save_data = 1;  %save data at end of processing (0=off 1=on)
node = 'MPD05';
low_range_mask = 0;
skip = 1

cd(d_read_data);
[Pythonfilename, Pythondir] = uigetfile('*.*','Select the sonde file', 'MultiSelect', 'on');
jj=1;

 %variable{1} = 'Aerosol_Backscatter_Coefficient_variance';
  
 variable{1} = 'time';
 variable{2} = 'range';
 variable{3} = 'Absolute_Humidity'; 
 variable{4} = 'Absolute_Humidity_mask';
% variable{5} = 'Absolute_Humidity_variance'; 
 variable{5} = 'Absolute_Humidity_uncertainty'; 
%  variable{6} = 'Aerosol_Backscatter_Coefficient';
%  variable{7} = 'Aerosol_Backscatter_Coefficient_mask';
%  variable{8} = 'Aerosol_Backscatter_Coefficient_variance';
 variable{9} = 'Absolute_Humidity_ic1';  % intitial conditions starting at 0 g/m^3 
 variable{10} = 'Absolute_Humidity_ic2'; % intitial conditions starting at 20 g/m^3 
    

for jj = 1:size(Pythonfilename,2)
  filename = Pythonfilename{jj};
%  date = filename(end-15:end-10);
  date = filename(7:12);
  n = datenum(date, 'yymmdd');
  ncid = netcdf.open(filename, 'NC_NOWRITE');
    %ncdisp(filename, '/', 'min') % use this to display all variables
    ncdisp(filename) % use this to display all variables
    time{jj} = ncread(filename,variable{1}); 
    alt{jj} = ncread(filename,variable{2});
  
    AH{jj}  = ncread(filename,variable{3});  
%    ic1{jj}  = ncread(filename,variable{9});  
%    ic2{jj}  = ncread(filename,variable{10});  
%    AH{jj} = (ic1{jj} + ic2{jj})/2; 
     
    AH_mask{jj} = ncread(filename,variable{4}); 
    AH_var{jj} = ncread(filename,variable{5}); 
    AH{jj}(AH_mask{jj} == 1) = nan;
%     AH_var{jj}(AH_mask{jj} == 1) = nan; 
    % mask the absolute humidity data based on its variance
     AH{jj}(AH_var{jj} >= 1.5) = nan;
%     AH_var{jj}(AH_var{jj} >= 5) = nan; 
        
%     ABC{jj}  = ncread(filename,variable{6});   
%     ABC_mask{jj} = ncread(filename,variable{7}); 
%     ABC_var{jj} = ncread(filename,variable{8});
%     ABC{jj}(ABC_mask{jj} == 1) = nan;
%     ABC_var{jj}(ABC_mask{jj} == 1) = nan; 
    
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
%       comb_AH_var = AH_var{jj}';
%       comb_ABC = ABC{jj}';
%       comb_ABC_var = ABC_var{jj}';
  else
      comb_duration = [comb_duration; duration{jj}];
      % find the maximum range to accumulate (in case it changes)
      max_range = min(cellfun('size',alt,1))
      comb_AH = [comb_AH(:,1:max_range); AH{jj}(1:max_range,:)'];
%       comb_AH_var = [comb_AH_var(:,1:max_range); AH_var{jj}(1:max_range,:)'];
%       comb_ABC = [comb_ABC(:,1:max_range); ABC{jj}(1:max_range,:)'];
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

x = comb_duration;
y = double(alt{1}(1:max_range,1)')/1000;
Z = real(comb_AH)';
xData =  linspace( fix(min(x)),  ceil(max(x)), round((ceil(max(x))-fix(min(x)))/skip)+1 );


% plot the AH
  figure1 = figure('Position',plot_size1);
  set(gcf,'renderer','zbuffer');
  h = pcolor(x, y, Z);
  set(h, 'EdgeColor', 'none'); 
  axis xy; colorbar('EastOutside'); 
  caxis([0 12]);
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
%   Z = real(comb_ABC)';
%   figure2 = figure('Position',plot_size1);
%   set(gcf,'renderer','zbuffer');
%   h = pcolor(x, y, Z);
%   set(h, 'EdgeColor', 'none'); 
%   axis xy; colorbar('EastOutside'); 
%   caxis([0 12]);
%   caxis([1e-8 1e-3]);
%   axis([fix(min(x)) ceil(max(x)) 0 6])
%   %  shading interp
%   set(gca, 'XTick',  xData)
%   set(gca,'TickDir','out');
%   set(gca,'TickLength',[0.005; 0.0025]);
%   set(gca,'Zscale', 'log')
%   set(gca,'Colorscale', 'log')
%   set(gca,'Zscale', 'linear')
%   hh = title({[node, ' Aerosol Backscatter Coefficient m^{-1} sr^{-1}']},...
%        'fontweight','b','fontsize',font_size);     
%   ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
%   datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
%   colormap(jet)
%   set(gca,'Fontsize',font_size,'Fontweight','b');
 
 
 

 if flag.save_data == 1
  
  cd(d_save_data);     
  range = alt{1}'; 
  N_avg_comb = (comb_AH./1e6.*6.022E23./18.015);
  duration = duration;

  FigH = figure(1);
  set(gca,'Fontsize',16,'Fontweight','b'); 
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920 250]);
  name=strcat(date, node, ' WV_Python_multi'); 
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
  
%   FigH = figure(2);
%   set(gca,'Fontsize',16,'Fontweight','b'); 
%   set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920 250]);
%   name=strcat(date, node, ' RB_Python_multi'); 
%   print(FigH, name, '-dpng', '-r0') % set at the screen resolution
  
  
  %cd('/Users/spuler/Desktop/WV_DIAL_data') % point to the directory where data is stored 
  name=strcat(date, '_combined');
  if strcmp(node,'MPD01')==1
      MPD01.N_avg_comb = N_avg_comb;
     % MPD01.RB_comb = RB_comb;
      MPD01.range = range;
      MPD01.time = duration;
      save(name, 'MPD01')
  elseif strcmp(node,'MPD02')==1
      MPD02.N_avg_comb = N_avg_comb;
    %  MPD02.RB_comb = RB_comb;
      MPD02.range = range;
      MPD02.time = duration;
      save(name, 'MPD02')
  elseif strcmp(node,'MPD03')==1
      MPD03.N_avg_comb = N_avg_comb;
    %  MPD03.RB_comb = RB_comb;
      MPD03.range = range;
      MPD03.time = duration;
      save(name, 'MPD03')
  elseif strcmp(node,'MPD04')==1
      MPD04.N_avg_comb = N_avg_comb;
    %  MPD04.RB_comb = RB_comb;
      MPD04.range = range;
      MPD04.time = duration;
      save(name, 'MPD04')
   elseif strcmp(node,'MPD05')==1
      MPD05.N_avg_comb = N_avg_comb;
   %   MPD05.RB_comb = RB_comb;
      MPD05.range = range;
      MPD05.time =  duration;
      save(name, 'MPD05')
  end

 end
  