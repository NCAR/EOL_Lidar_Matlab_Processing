clear all; close all

serv_path = '/Users/spuler/Desktop';
cd(strcat(serv_path,'/mpd/V01')) %test_1
%cd(strcat(serv_path,'/mpd/V02')) 

[Etalonfilename, Etalondir] = uigetfile('*.*','Select the sonde file', 'MultiSelect', 'on');
variable{1} = 'WavemeterWavelength';
variable{2} = 'WavemeterTime'; 
variable{3} = 'DetectorSignal';
variable{4} = 'DetectorTime';
  
jj=1;
 
for jj = 1:size(Etalonfilename,2)
  filename = Etalonfilename{jj};
  date = filename(end-15:end-10);
  n = datenum(date, 'yymmdd');
  ncid = netcdf.open(filename, 'NC_NOWRITE');
    %ncdisp(filename, '/', 'min') % use this to display all variables
    ncdisp(filename) % use this to display all variables

    AOI{jj} = ncreadatt(filename,'/', 'AOI');
    EtalonType{jj} =  ncreadatt(filename, '/', 'EtalonType');
    Expansion{jj} =  ncreadatt(filename, '/', 'Expansion');
    FiberType{jj} =  ncreadatt(filename, '/', 'FiberType');
    
    WavemeterWavelength{jj}  = ncread(filename,variable{1});
    WavemeterTime{jj}  = ncread(filename,variable{2});  
    DetectorSignal{jj} = ncread(filename,variable{3}); 
    DetectorTime{jj} = ncread(filename,variable{4});
    
  netcdf.close(ncid); 
  
  % interpolate the wavelength position to the detector signal time grid
  WavelengthTime{jj} = interp1(WavemeterTime{jj}, WavemeterWavelength{jj}, DetectorTime{jj}, 'nearest');
end


figure(1)
plot(DetectorTime{1},DetectorSignal{1})

figure(2)
plot(WavemeterTime{1},WavemeterWavelength{1})

  figure(3) 
  % AO1=2.25, 4x, SM, FWHM=8.6pm
  wave_comb = vertcat(WavelengthTime{1:4});   
  signal_comb = vertcat(DetectorSignal{1:4});   
  [C,ia,idx] = unique(wave_comb );  %find the unique wavelengths and their indices 
  mean_signal = accumarray(idx,signal_comb,[],@mean); % calc mean at unique wavelengths
  max(mean_signal)/2;
  idx = find(mean_signal > max(mean_signal)/2);
  idx_diff = find(diff(idx)>2);
  FWHM = C(idx(1)+idx_diff)- C(idx(1))
  plot(C-.004,mean_signal, 'ko', ...      
     'DisplayName', ['AOI=', num2str(AOI{1}), num2str(Expansion{1}), num2str(FiberType{1})])
  legend

  % AO1=2.25, 4x, MM, FWHM=x.xpm
  hold on
  wave_comb = vertcat(WavelengthTime{5:8});   
  signal_comb = vertcat(DetectorSignal{5:8});
  [C,ia,idx] = unique(wave_comb );  %find the unique wavelengths and their indices 
  mean_signal = accumarray(idx,signal_comb,[],@mean); % calc mean at unique wavelengths 
  max(mean_signal)/2;
  idx = find(mean_signal > max(mean_signal)/2);
  idx_diff = find(diff(idx)>2);
  FWHM = C(idx(1)+idx_diff)- C(idx(1))
  plot(C,mean_signal, 'ro', ...    
     'DisplayName', ['AOI= ', num2str(AOI{5}), num2str(Expansion{5}), num2str(FiberType{5})])
  legend
  hold off

  % AO1=2.25, 2x, MM, FWHM=17.6pm
  hold on
  wave_comb = vertcat(WavelengthTime{9:12});   
  signal_comb = vertcat(DetectorSignal{9:12});
  [C,ia,idx] = unique(wave_comb );  %find the unique wavelengths and their indices 
  mean_signal = accumarray(idx,signal_comb,[],@mean); % calc mean at unique wavelengths 
  max(mean_signal)/2;
  idx = find(mean_signal > max(mean_signal)/2);
  idx_diff = find(diff(idx)>2);
  FWHM = C(idx(1)+idx_diff)- C(idx(1))
  plot(C,mean_signal, 'go', ...   
    'DisplayName', ['AOI=', num2str(AOI{9}), num2str(Expansion{9}), num2str(FiberType{9})])
  legend  
  hold off
  
  % AO1=2.25, 2x, SM, FWHM=x.xpm
  hold on
  wave_comb = vertcat(WavelengthTime{13:19});   
  signal_comb = vertcat(DetectorSignal{13:19});
  [C,ia,idx] = unique(wave_comb);  %find the unique wavelengths and their indices 
  mean_signal = accumarray(idx,signal_comb,[],@mean); % calc mean at unique wavelengths 
  max(mean_signal)/2;
  idx = find(mean_signal > max(mean_signal)/2);
  idx_diff = find(diff(idx)>2);
  FWHM = C(idx(1)+idx_diff)- C(idx(1))
  plot(C,mean_signal, 'bo', ...   
    'DisplayName', ['AOI=', num2str(AOI{13}), num2str(Expansion{13}), num2str(FiberType{13})])
  legend
  hold off
  
 grid on
title('Etalon V02')
ylabel('rel. trans. power'); 
xlabel('wavelength (nm)'); 
xlim([769.72  770.1])
ylim([1e-6 4e-4])
set(gca,'Fontsize',20,'Fontweight','b'); %

ylim([1e-5 1e-3])
ylim([0 4e-4])
xlim([769.74  769.79]) 
  
  
  
%  figure(4) 
  % AO1=0.0, 4x, MM, FWHM=x.xpm
  hold on
  wave_comb = vertcat(WavelengthTime{20:22});   
  signal_comb = vertcat(DetectorSignal{20:22});
  [C,ia,idx] = unique(wave_comb);  %find the unique wavelengths and their indices 
  mean_signal = accumarray(idx,signal_comb,[],@mean); % calc mean at unique wavelengths 
  max(mean_signal)/2;
  idx = find(mean_signal > max(mean_signal)/2);
  idx_diff = find(diff(idx)>3);
  FWHM = C(idx(1)+idx_diff)- C(idx(1))
  plot(C,mean_signal, 'b+', ...   
    'DisplayName', ['AOI=', num2str(AOI{20}), num2str(Expansion{20}), num2str(FiberType{20})])
  legend
  hold off

 % AO1=0.0, 4x, SM, FWHM=x.xpm
  hold on
  wave_comb = vertcat(WavelengthTime{23:26});   
  signal_comb = vertcat(DetectorSignal{23:26});
  [C,ia,idx] = unique(wave_comb);  %find the unique wavelengths and their indices 
  mean_signal = accumarray(idx,signal_comb,[],@mean); % calc mean at unique wavelengths 
  max(mean_signal)/2;
  idx = find(mean_signal > max(mean_signal)/2);
  idx_diff = find(diff(idx)>3);
  FWHM = C(idx(1)+idx_diff)- C(idx(1))
  plot(C,mean_signal, 'g+', ...     
    'DisplayName', ['AOI=', num2str(AOI{23}), num2str(Expansion{23}), num2str(FiberType{23})])
  legend
  hold off
  
  % AO1=0.0, 2x, SM, FWHM=8.1pm
  hold on
  wave_comb = vertcat(WavelengthTime{27:32});   
  signal_comb = vertcat(DetectorSignal{27:32});
  [C,ia,idx] = unique(wave_comb);  %find the unique wavelengths and their indices 
  mean_signal = accumarray(idx,signal_comb,[],@mean); % calc mean at unique wavelengths 
  max(mean_signal)/2;
  idx = find(mean_signal > max(mean_signal)/2);
  idx_diff = find(diff(idx)>3);
  FWHM = C(idx(1)+idx_diff)- C(idx(1))
  plot(C,mean_signal, 'r+', ... 
    'DisplayName', ['AOI=', num2str(AOI{27}), num2str(Expansion{27}), num2str(FiberType{27})])
  legend
  hold off

  % AO1=0.0, 2x, MM, FWHM=x.xpm  
  hold on
  wave_comb = vertcat(WavelengthTime{33:38});   
  signal_comb = vertcat(DetectorSignal{33:38});
  [C,ia,idx] = unique(wave_comb);  %find the unique wavelengths and their indices 
  mean_signal = accumarray(idx,signal_comb,[],@mean); % calc mean at unique wavelengths 
  max(mean_signal)/2;
  idx = find(mean_signal > max(mean_signal)/2);
  idx_diff = find(diff(idx)>3);
  FWHM = C(idx(1)+idx_diff)- C(idx(1))
  plot(C,mean_signal, 'k+', ... 
    'DisplayName', ['AOI=', num2str(AOI{33}), num2str(Expansion{33}), num2str(FiberType{33})])
  legend  
  hold off
 

grid on
title('Etalon V01')
ylabel('rel. trans. power'); 
xlabel('wavelength (nm)'); 
xlim([769.72  770.15])
ylim([1e-6 4e-4])
set(gca,'Fontsize',20,'Fontweight','b'); %

%ylim([1e-5 1e-3])
%ylim([0 4e-4])
%xlim([769.74  769.79])

  FigH = figure(3);
%   set(gca,'Fontsize',36,'Fontweight','b'); 
   set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1  1  1600 1200]);
   name=strcat('Etalon_V01_data'); 
   print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
 


  