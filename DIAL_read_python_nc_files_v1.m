clear all
close all

flag.save_fig = 1;
flag.save_data = 1;

write_data_folder = uipickfiles('num',1,'out', 'char', 'prompt', ...
    'select folder to store data',  'FilterSpec', '/Volumes/documents/WV_DIAL_data/');
% take care to pick resolution 
files = uipickfiles('prompt', 'select data files to process', 'FilterSpec', '/scr/sci/mhayman/DIAL/Processed_Data/MPDSGP/t_res_5min/');

j=1;

for j = 1:size(files,2)
    folder = (files{j});
    date = textscan(folder(end-15:end-9), '%6f'); date=date{1};  % read date of file
    n = datenum(num2str(date), 'yymmdd');
    ncid = netcdf.open(folder, 'NC_NOWRITE');
    variable{1} = 'Absolute_Humidity';
    variable{2} = 'Attenuated_Backscatter';
    for i = 1:size(variable,2)
      var_units{i} = ncreadatt(folder, variable{i},'units');
      var_units{i} = erase(var_units{i}, '$');
      var_time{i} = ncread(folder,horzcat('time_',variable{i}));   
      var_range{i} = ncread(folder,horzcat('range_',variable{i}));   
      var{i} = ncread(folder, variable{i});
      
      var_variance{i} = ncread(folder, horzcat(variable{i},'_variance')); 
      var_mask{i} = ncread(folder,horzcat(variable{i},'_mask')); 
      x{i} = n+double(var_time{i}/3600/24);
      y{i} = double(var_range{i}');
    end  
  netcdf.close(ncid);
  scrsz = get(0,'ScreenSize');
  date_plot=datestr(n, 'dd mmm yyyy');
  plot_size1 = [scrsz(4)/1.5 scrsz(4)/10 scrsz(3)/1.5 scrsz(4)/3];
  font_size = 14;
  xData =  linspace(fix(min(x{1})),  ceil(max(x{1})), 25);

  % add masking
   for i = 1:size(variable,2)
     var{i}(var_mask{i} == 1) = nan;
   end

   for i = 1:size(variable,2) 
  %   plot variable 1
    figure1 = figure('Position',plot_size1);
    set(gcf,'renderer','zbuffer');
    h = pcolor(x{i}, y{i}/1000, real(log10(var{i})));
    set(h, 'EdgeColor', 'none'); 
    axis xy; colorbar('EastOutside'); %caxis([0 10]);
    title({[date_plot ,' ',replace(variable{i}, '_', ' '), ' [',var_units{i},']']},...
       'fontweight','b','fontsize',font_size)%,'Interpreter', 'none')
    ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
    axis([n n+1 0 12])
    set(gca, 'XTick',  xData)
    datetick('x','HH','keeplimits', 'keepticks');
    colormap(jet)
%    shading interp
    set(gca,'Fontsize',font_size,'Fontweight','b');
   end
  %caxis([0 2.5]);
 
  %plot combined backscatter and humidity
  figure10 = figure('visible', 'off','Position',[scrsz(4)/1.5 scrsz(4)/10 scrsz(3)/1.5 scrsz(4)/1.5]);
  set(figure10, 'visible', 'on', 'PaperUnits', 'points', 'PaperPosition', [0 0 1280 800]);
  subplot1=subplot(2,1,2,'Parent',figure10);
  box(subplot1,'on');
  set(gcf,'renderer','zbuffer');
  h = pcolor(x{1},y{1}/1000, var{1}); %real(log10(var{1})));
  set(h, 'EdgeColor', 'none');
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  set(gca, 'XTick',  xData)
  colorbar('EastOutside');
  axis([fix(min(x{1})) fix(min(x{1}))+1 0 6])
  caxis([1 20]);
  %caxis([1 5]);
  datetick('x','HH','keeplimits', 'keepticks');
  colormap(jet)
  %shading interp 
  hh =  title({[date_plot ,' ',replace(variable{1}, '_', ' '), ' [',var_units{1},']']},...
       'fontweight','b','fontsize',font_size)%,'Interpreter', 'none')
  P_t = get(hh, 'Position');
  set(hh,'Position', [P_t(1) P_t(2)+0.2 P_t(3)])
  xlabel('Time (UTC)','fontweight','b','fontsize',font_size);
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size);
  set(gca,'Fontsize',font_size,'Fontweight','b');

  set(gca,'Zscale', 'linear')
  set(gca,'Colorscale', 'linear')
  set(gca,'Zscale', 'linear')

  subplot1=subplot(2,1,1,'Parent',figure10);
  box(subplot1,'on'); %(number density in mol/cm3)(1e6 cm3/m3)/(N_A mol/mole)*(18g/mole)
  set(gcf,'renderer','zbuffer');
  h = pcolor(x{2},y{2}/1000,var{2});
  set(h, 'EdgeColor', 'none');
  set(gca, 'XTick',  xData)
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  colorbar('EastOutside');
  axis([fix(min(x{2})) fix(min(x{2}))+1 0 12])
  caxis([1e3 1e9]);
  datetick('x','HH','keeplimits', 'keepticks');
  colormap(jet)
  %shading interp
  hh =  title({[date_plot ,' ',replace(variable{2}, '_', ' '), ' [',var_units{2},']']},...
       'fontweight','b','fontsize',font_size)%,'Interpreter', 'none')
  P_t = get(hh, 'Position');
  set(hh,'Position', [P_t(1) P_t(2)+0.2 P_t(3)])
  xlabel('Time (UTC)','fontweight','b','fontsize',font_size); 
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size); 
  set(gca,'Fontsize',font_size,'Fontweight','b');
 
  set(gca,'Zscale', 'log')
  set(gca,'Colorscale', 'log')
  set(gca,'Zscale', 'linear')
  
  
    if flag.save_fig == 1
      %cd('/Volumes/documents/WV_DIAL_data/plots/') % point to the directory where data is stored 
      cd(write_data_folder)
      Plt_size = [scrsz(4)/1.5 scrsz(4)/10 scrsz(3)/1.5 scrsz(4)/1.5];
      FigH = figure(3);
      set(FigH, 'PaperUnits', 'points', 'PaperPosition', Plt_size);
      name=strcat(date_plot, 'Python'); 
      print(FigH, name, '-dpng', '-r0') % set the resolution as 300 dpiFigH = figure(1);
    end

    if flag.save_data == 1
      cd(write_data_folder)
      name=strcat(date_plot);
       time_new = x{1};
       range = y{1};
       RB = var{2}; 
       N_avg = var{1}./1e6.*6.022E23./18.015; % convert to number density      
      save(name, 'N_avg', 'RB', 'range', 'time_new')
    end
 close all
end
