%cd /scr/sci/spuler/mpd/sgp/raman_lidar
clear all; close all
%serv_path = '/Volumes/eol/sci/spuler';
%serv_path = '/Users/spuler/Desktop';
%serv_path = '/Volumes/eol/fog1/rsfdata/MPD';
%cd(strcat(serv_path,'/mpd_05_processed_data'))
serv_path = '/Volumes/eol/sci/mhayman/DIAL/Processed_Data';
%cd(strcat(serv_path,'/mpd/Marshall/mpd03_10min')) % this is the 10min averaged data
%cd(strcat(serv_path,'/mpd/Marshall/mpd03_5min'))
%cd(strcat(serv_path,'/mpd/sgp/mpd05'))
%cd(strcat(serv_path, '/MPDSGP/t_res_10min/'));  %This was the original SGP 10 min data
cd(strcat(serv_path,'/MPD_Gen5_Pub/'));
d_read_data = pwd; % get the current path

serv_path = '/Volumes/eol/fog1/rsfdata/MPD';
cd(strcat(serv_path,'/mpd_05_processed_data/Matlab'))
flag.save_data = 1;  %save data at end of processing (0=off 1=on)
MPD = 5;
low_range_mask = 0;
d_save_data = pwd; % get the current path


cd(d_read_data);
[Pythonfilename, Pythondir] = uigetfile('*.*','Select the sonde file', 'MultiSelect', 'on');
jj=1;

 variable{1} = 'Absolute_Humidity';
 %variable{2} = 'time_Absolute_Humidity';
 variable{2} = 'time';
 %variable{3} = 'range_Absolute_Humidity';
 variable{3} = 'range';
 variable{4} = 'Absolute_Humidity_mask';
 variable{5} = 'Absolute_Humidity_variance'; 
 
 
% variable{5} = 'Attenuated_Backscatter';
% variable{6} = 'time_Attenuated_Backscatter ';
% variable{7} = 'range_Attenuated_Backscatter';
     

for jj = 1:size(Pythonfilename,2)
  filename = Pythonfilename{jj};
  date = filename(end-15:end-10);
  n = datenum(date, 'yymmdd');
  ncid = netcdf.open(filename, 'NC_NOWRITE');
    %ncdisp(filename, '/', 'min') % use this to display all variables
    ncdisp(filename) % use this to display all variables
    AH{jj}  = ncread(filename,variable{1});   
    AH_time{jj} = ncread(filename,variable{2}); 
    AH_alt{jj} = ncread(filename,variable{3});
    AH_mask{jj} = ncread(filename,variable{4}); 
    AH_var{jj} = ncread(filename,variable{5}); 
    AH{jj}(AH_mask{jj} == 1) = nan;
    AH_var{jj}(AH_mask{jj} == 1) = nan; 
  netcdf.close(ncid); 
  %convert from Unix time to date number (days since Jan 0 0000) 
  AH_duration{jj} =  n+double(AH_time{jj}/3600/24);

  %range_grid = (0:75:12000);
  %AH_grid{jj} = interp1( AH{jj}(:,1) ,AH{jj}(:,2:end), range_grid, 'nearest', 'extrap'); 
 % for jj=1:12  (there are problems with the data for the first three days)

  % mask the lowest gates if desired
  low_mask = repmat(AH_alt{jj}, 1, size(AH{jj}, 2)); 
  AH{jj}(low_mask < low_range_mask)= nan;
  AH_var{jj}(low_mask < low_range_mask)= nan;
 
  if jj == 1
      comb_AH_duration = AH_duration{jj};
      comb_AH = AH{jj}';
      comb_AH_var = AH_var{jj}';
  else
      comb_AH_duration = [comb_AH_duration; AH_duration{jj}];
      % find the maximum range to accumulate (in case it changes)
      max_range = min(cellfun('size',AH_alt,1))
      comb_AH = [comb_AH(:,1:max_range); AH{jj}(1:max_range,:)'];
      comb_AH_var = [comb_AH_var(:,1:max_range); AH_var{jj}(1:max_range,:)'];
  end
  
  % remove non-physical values AH < - 0.5 (allow for some slight negative noise)
  comb_AH(comb_AH < -0.5)= nan;
  
 % end
  
end

scrsz = get(0,'ScreenSize');
date=datestr(n, 'yyyy-mmm-dd');
plot_size1 = [scrsz(4)/1.5 scrsz(4)/10 scrsz(3)/1.5 scrsz(4)/3];
font_size = 14;


x = comb_AH_duration;
y = double(AH_alt{1}(1:max_range,1)')/1000;
Z = real(comb_AH)';


% plot the AH
  figure1 = figure('Position',plot_size1);
  set(gcf,'renderer','zbuffer');
  h = pcolor(x, y, Z);
  set(h, 'EdgeColor', 'none'); 
  axis xy; colorbar('EastOutside'); 
  caxis([0 16]);
  caxis([0 6]);
  axis([fix(min(x)) ceil(max(x)) 0 6])
  title({'MPD Absolute Humidity (g m^{-3})'},...
       'fontweight','b','fontsize',font_size)
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
  %datetick('x','HH:MM:SS');
  %datetick('x','mm/dd', 'keeplimits');%, 'keepticks');
  datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
  %xlabel('Day','fontweight','b','fontsize',font_size);
  colormap(jet)
%  shading interp
 set(gca,'TickDir','out');
 set(gca,'TickLength',[0.005; 0.0025]);
 set(gca,'Fontsize',font_size,'Fontweight','b');

  %datetick('x','mm/dd', 'keeplimits');
  %datetick('x','dd-mmm-yy','keeplimits', 'keepticks');
  %Scrnsize = [scrsz(4)/2 scrsz(4)/10 scrsz(3)/1 scrsz(4)/1.5]; % use for standard plots
  %Scrnsize = [scrsz(4)/1 scrsz(4)/1 scrsz(3)/0.30 scrsz(4)/2]; % use for ILRC really long plots
 % cd('/Volumes/documents/WV_DIAL_data/plots/') % point to the directory where data is stored 
 
 
  FigH = figure(1);
  set(gca,'Fontsize',16,'Fontweight','b'); 
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1600 400]);
  name=strcat(date, 'Python_multi'); 
%  print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
 

 if flag.save_data == 1
  
  cd(d_save_data);     
  range = AH_alt{1}'; 
  N_avg_comb = (comb_AH./1e6.*6.022E23./18.015);
  duration = comb_AH_duration;

  %cd('/Users/spuler/Desktop/WV_DIAL_data') % point to the directory where data is stored 
  name=strcat(date, '_combined');
  if MPD==1
      MPD01.N_avg_comb = N_avg_comb;
     % MPD01.RB_comb = RB_comb;
      MPD01.range = range;
      MPD01.time = duration;
      save(name, 'MPD01')
  elseif MPD==2
      MPD02.N_avg_comb = N_avg_comb;
    %  MPD02.RB_comb = RB_comb;
      MPD02.range = range;
      MPD02.time = duration;
      save(name, 'MPD02')
  elseif MPD==3
      MPD03.N_avg_comb = N_avg_comb;
    %  MPD03.RB_comb = RB_comb;
      MPD03.range = range;
      MPD03.time = duration;
      save(name, 'MPD03')
  elseif MPD==4
      MPD04.N_avg_comb = N_avg_comb;
    %  MPD04.RB_comb = RB_comb;
      MPD04.range = range;
      MPD04.time = duration;
      save(name, 'MPD04')
   elseif MPD==5
      MPD05.N_avg_comb = N_avg_comb;
   %   MPD05.RB_comb = RB_comb;
      MPD05.range = range;
      MPD05.time =  duration;
      save(name, 'MPD05')
  end

 end
  