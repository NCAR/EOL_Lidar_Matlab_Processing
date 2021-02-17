clear all; close all

serv_path = '/Users/spuler/Desktop';
cd(strcat(serv_path,'/mpd/V02')) 

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
  % AO1=0, 2x, MM, FWHM=8.2pm
  wave_comb = vertcat(WavelengthTime{1:10});   
  signal_comb = vertcat(DetectorSignal{1:10});   
  [C,ia,idx] = unique(wave_comb );  %find the unique wavelengths and their indices 
  mean_signal = accumarray(idx,signal_comb,[],@mean); % calc mean at unique wavelengths
  max(mean_signal)/2;
  idx = find(mean_signal > max(mean_signal)/2);
  idx_diff = find(diff(idx)>2);
  FWHM = C(idx(1)+idx_diff)- C(idx(1))
  semilogy(C-.004,mean_signal, 'ko', ...      
     'DisplayName', ['AOI=', num2str(AOI{1}), num2str(Expansion{1}), num2str(FiberType{1})])
  legend

  % AO1=0, 2x, SM, FWHM=8.1pm
  hold on
  wave_comb = vertcat(WavelengthTime{11:15});   
  signal_comb = vertcat(DetectorSignal{11:15});
  [C,ia,idx] = unique(wave_comb );  %find the unique wavelengths and their indices 
  mean_signal = accumarray(idx,signal_comb,[],@mean); % calc mean at unique wavelengths 
  max(mean_signal)/2;
  idx = find(mean_signal > max(mean_signal)/2);
  idx_diff = find(diff(idx)>2);
  FWHM = C(idx(1)+idx_diff)- C(idx(1))
  semilogy(C,mean_signal, 'ro', ...    
     'DisplayName', ['AOI= ', num2str(AOI{11}), num2str(Expansion{11}), num2str(FiberType{11})])
  legend
  hold off

  % AO1=0, 4x, SM, FWHM=7.4pm
  hold on
  wave_comb = vertcat(WavelengthTime{16:21});   
  signal_comb = vertcat(DetectorSignal{16:21});
  [C,ia,idx] = unique(wave_comb );  %find the unique wavelengths and their indices 
  mean_signal = accumarray(idx,signal_comb,[],@mean); % calc mean at unique wavelengths 
  max(mean_signal)/2;
  idx = find(mean_signal > max(mean_signal)/2);
  idx_diff = find(diff(idx)>2);
  FWHM = C(idx(1)+idx_diff)- C(idx(1))
  semilogy(C,mean_signal, 'go', ...   
    'DisplayName', ['AOI=', num2str(AOI{16}), num2str(Expansion{16}), num2str(FiberType{16})])
  legend  
  hold off
  
  % AO1=0, 4x, MM, FWHM=8.1pm
  hold on
  wave_comb = vertcat(WavelengthTime{22:27});   
  signal_comb = vertcat(DetectorSignal{22:27});
  [C,ia,idx] = unique(wave_comb);  %find the unique wavelengths and their indices 
  mean_signal = accumarray(idx,signal_comb,[],@mean); % calc mean at unique wavelengths 
  max(mean_signal)/2;
  idx = find(mean_signal > max(mean_signal)/2);
  idx_diff = find(diff(idx)>2);
  FWHM = C(idx(1)+idx_diff)- C(idx(1))
  semilogy(C,mean_signal, 'bo', ...   
    'DisplayName', ['AOI=', num2str(AOI{22}), num2str(Expansion{22}), num2str(FiberType{22})])
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
  % AO1=4.7, 4x, MM, FWHM=27.3pm
  hold on
  wave_comb = vertcat(WavelengthTime{28:36});   
  signal_comb = vertcat(DetectorSignal{28:36});
  [C,ia,idx] = unique(wave_comb);  %find the unique wavelengths and their indices 
  mean_signal = accumarray(idx,signal_comb,[],@mean); % calc mean at unique wavelengths 
  max(mean_signal)/2;
  idx = find(mean_signal > max(mean_signal)/2);
  idx_diff = find(diff(idx)>3);
  FWHM = C(idx(1)+idx_diff)- C(idx(1))
  semilogy(C,mean_signal, 'b+', ...   
    'DisplayName', ['AOI=', num2str(AOI{28}), num2str(Expansion{28}), num2str(FiberType{28})])
  legend
  hold off

 % AO1=4.7, 4x, SM, FWHM=8.6pm
  hold on
  wave_comb = vertcat(WavelengthTime{37:48});   
  signal_comb = vertcat(DetectorSignal{37:48});
  [C,ia,idx] = unique(wave_comb);  %find the unique wavelengths and their indices 
  mean_signal = accumarray(idx,signal_comb,[],@mean); % calc mean at unique wavelengths 
  max(mean_signal)/2;
  idx = find(mean_signal > max(mean_signal)/2);
  idx_diff = find(diff(idx)>3);
  FWHM = C(idx(1)+idx_diff)- C(idx(1))
  semilogy(C,mean_signal, 'g+', ...     
    'DisplayName', ['AOI=', num2str(AOI{37}), num2str(Expansion{37}), num2str(FiberType{37})])
  legend
  hold off
  
  % AO1=4.7, 4x, SM, FWHM=8.1pm
  hold on
  wave_comb = vertcat(WavelengthTime{49:55});   
  signal_comb = vertcat(DetectorSignal{49:55});
  [C,ia,idx] = unique(wave_comb);  %find the unique wavelengths and their indices 
  mean_signal = accumarray(idx,signal_comb,[],@mean); % calc mean at unique wavelengths 
  max(mean_signal)/2;
  idx = find(mean_signal > max(mean_signal)/2);
  idx_diff = find(diff(idx)>3);
  FWHM = C(idx(1)+idx_diff)- C(idx(1))
  semilogy(C,mean_signal, 'r+', ... 
    'DisplayName', ['AOI=', num2str(AOI{49}), num2str(Expansion{49}), num2str(FiberType{49})])
  legend
  hold off

  % AO1=4.7, 2x, MM, FWHM=29.7pm  
  hold on
  wave_comb = vertcat(WavelengthTime{56:68});   
  signal_comb = vertcat(DetectorSignal{56:68});
  [C,ia,idx] = unique(wave_comb);  %find the unique wavelengths and their indices 
  mean_signal = accumarray(idx,signal_comb,[],@mean); % calc mean at unique wavelengths 
  max(mean_signal)/2;
  idx = find(mean_signal > max(mean_signal)/2);
  idx_diff = find(diff(idx)>3);
  FWHM = C(idx(1)+idx_diff)- C(idx(1))
   semilogy(C,mean_signal, 'k+', ... 
    'DisplayName', ['AOI=', num2str(AOI{56}), num2str(Expansion{56}), num2str(FiberType{56})])
  legend  
  hold off
 

grid on
title('Etalon V02')
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
   name=strcat('Etalon_V02_data'); 
   print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
 


  