%cd /scr/sci/spuler/mpd/sgp/raman_lidar
clear all; close all

serv_path = '/Volumes/eol/sci/';
cd(strcat(serv_path,'spuler/mpd/sgp/raman_lidar/'))
[Ramanfilename, Ramandir] = uigetfile('*.*','Select the sonde file', 'MultiSelect', 'on');
jj=1;
%cd= d;

   variable{1} = 'base_time'; % time (seconds since 1970-1-1)
   variable{2} = 'time_offset'; % time (seconds since 1970-1-1)
   variable{3} = 'height'; % height above ground level (km)
   variable{4} = 'mr_merged'; % 'High and Low merged water vapor mixing ratio'
   variable{5} = 'temp_sonde'; % 'sonde temperature' 
   variable{6} = 'pres_sonde'; % 'sonde pressure' 

for jj = 1:size(Ramanfilename,2)
  filename = Ramanfilename{1};
  date = filename(end-15:end-10);
  n = datenum(date, 'yymmdd');
  ncid = netcdf.open(filename, 'NC_NOWRITE');
    %ncdisp(filename, '/', 'min') % use this to display all variables
    ncdisp(filename) % use this to display all variables
    base_time{jj}  = ncread(filename,variable{1});   
    Raman_time{jj} = ncread(filename,variable{2}); 
    Raman_alt{jj} = ncread(filename,variable{3});
    Raman_MR{jj} = ncread(filename,variable{4}); 
    temp_sonde{jj} = ncread(filename,variable{5}); 
    pres_sonde{jj} = ncread(filename,variable{6}); 
  netcdf.close(ncid); 
  %convert from Unix time to date number (days since Jan 0 0000) 
  Raman_duration{jj} = datenum(datetime(int64(Raman_time{jj}) + int64(base_time{jj}), 'convertfrom', 'posixtime'));
  % to convert from mixing ratio to absolute humidity to multiply by the
  % density of dry air.  At 20C density is 1.2041 kg.m^3
  % the mixing ratio is in kg/kg * kg/m^3 = kg/m^3
  % density p=Press/(R*Temp)
  % where R = 8.314E-2 m^3 mbar K^-1 mol^1, M.air = 28.97 g/mol
  density{jj} = pres_sonde{jj}./(8.31447E-2*temp_sonde{jj})*28.97;
  Raman_AH{jj} = density{jj}.*Raman_MR{jj}.*1e-3; %to convert from kg/m# to g/m^3

  if jj == 1
      comb_Raman_duration = Raman_duration{jj};
      comb_Raman_AH = Raman_AH{jj}';
  else
      comb_Raman_duration = [comb_Raman_duration; Raman_duration{jj}];
      comb_Raman_AH = [comb_Raman_AH; Raman_AH{jj}'];
  end
  
end


scrsz = get(0,'ScreenSize');
date=datestr(n, 'yyyy-mmm-dd');
plot_size1 = [scrsz(4)/1.5 scrsz(4)/10 scrsz(3)/1.5 scrsz(4)/3];
font_size = 14;

x = comb_Raman_duration;
y = double(Raman_alt{1}');
Z = real(comb_Raman_AH)';

x = Raman_duration{1};
y = double(Raman_alt{1}');
Z = real(Raman_AH{2});

 
% plot the Raman AH
  figure1 = figure('Position',plot_size1);
  set(gcf,'renderer','zbuffer');
  h = pcolor(x, y, Z);
  set(h, 'EdgeColor', 'none'); 
  axis xy; colorbar('EastOutside'); 
  caxis([0 18]);
  %ylim([0 6])
  %title({[date,' SGP Raman Lidar absolute humidity']},...
  %     'fontweight','b','fontsize',font_size)
  ylabel('range (km)','fontweight','b','fontsize',font_size); 
  %datetick('x','HH:MM:SS');
  datetick('x','HH:MM', 'keeplimits');%, 'keepticks');
  xlabel('Time (UTC)','fontweight','b','fontsize',font_size);
  colormap(jet)
%  shading interp
  set(gca,'Fontsize',font_size,'Fontweight','b');

  datetick('x','HH:MM', 'keeplimits');
  Scrnsize = [scrsz(4)/2 scrsz(4)/10 scrsz(3)/1 scrsz(4)/1.5]; % use for standard plots
  cd('/Volumes/documents/WV_DIAL_data/plots/') % point to the directory where data is stored 
  FigH = figure(1);
  set(gca,'Fontsize',36,'Fontweight','b'); 
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', Scrnsize);
  name=strcat(date, 'Raman_MR'); 
  print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
 
