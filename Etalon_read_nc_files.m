
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



%for ii=1:10
  wave_comb = vertcat(WavelengthTime{1:10});   
  signal_comb = vertcat(DetectorSignal{1:10});   
  [C,ia,idx] = unique(wave_comb );  %find the unique wavelengths and their indices (do not sort)
  mean_signal = accumarray(idx,signal_comb,[],@mean); % calc mean at unique wavelengths
  %semilogy(WavelengthTime{ii},DetectorSignal{ii}, 'ko--', ...
  figure(3) % AO1=0, 2x, MM, FWHM=8.2pm
  %hold on
  semilogy(C-.004,mean_signal, 'ko', ...      
     'DisplayName', ['AOI=', num2str(AOI{ii}), num2str(Expansion{ii}), num2str(FiberType{ii})])
  %hold off
  max(mean_signal)/2;
  idx = find(mean_signal > max(mean_signal)/2);
  idx_diff = find(diff(idx)>1);
  FWHM = C(idx(1)+idx_diff)- C(idx(1))
%end
legend

figure(4) % AO1=0, 2x, SM, FWHM=8.1pm
hold on
for ii=11:15
  [C,ia,idx] = unique(WavelengthTime{ii}); 
  mean_signal = accumarray(idx,DetectorSignal{ii},[],@mean);  
  max(mean_signal)/2;
  idx = find(mean_signal > max(mean_signal)/2);
  idx_diff = find(diff(idx)>1);
  FWHM = C(idx(1)+idx_diff)- C(idx(1))
  %semilogy(WavelengthTime{ii},DetectorSignal{ii}, 'ro', ...
  semilogy(C,mean_signal, 'ro', ...    
     'DisplayName', ['AOI= ', num2str(AOI{ii}), num2str(Expansion{ii}), num2str(FiberType{ii})])
  hold off

  max(mean_signal)/2;
  idx = find(mean_signal > max(mean_signal)/2);
  idx_diff = find(diff(idx)>1);
  FWHM = C(idx(1)+idx_diff)- C(idx(1))
  
end
legend

figure(5)  % AO1=0, 4x, SM, FWHM=7.4pm
hold on
for ii=16:21
  [C,ia,idx] = unique(WavelengthTime{ii});
  mean_signal = accumarray(idx,DetectorSignal{ii},[],@mean); 
  max(mean_signal)/2;
  idx = find(mean_signal > max(mean_signal)/2);
  idx_diff = find(diff(idx)>1);
  FWHM = C(idx(1)+idx_diff)- C(idx(1))
  
  %semilogy(WavelengthTime{ii},DetectorSignal{ii}, 'go', ...
  semilogy(C,mean_signal, 'go', ...   
    'DisplayName', ['AOI=', num2str(AOI{ii}), num2str(Expansion{ii}), num2str(FiberType{ii})])
  hold off
  
  max(mean_signal)/2;
  idx = find(mean_signal > max(mean_signal)/2);
  idx_diff = find(diff(idx)>1);
  FWHM = C(idx(1)+idx_diff)- C(idx(1))
  
end
legend

figure(6) % AO1=0, 4x, MM, FWHM=8.3pm
hold on
for ii=22:27
  [C,ia,idx] = unique(WavelengthTime{ii});
  mean_signal = accumarray(idx,DetectorSignal{ii},[],@mean); 
  %semilogy(WavelengthTime{ii},DetectorSignal{ii}, 'bo', ...
  semilogy(C,mean_signal, 'bo', ...   
    'DisplayName', ['AOI=', num2str(AOI{ii}), num2str(Expansion{ii}), num2str(FiberType{ii})])
  hold off
  
  max(mean_signal)/2;
  idx = find(mean_signal > max(mean_signal)/2);
  idx_diff = find(diff(idx)>1);
  FWHM = C(idx(1)+idx_diff)- C(idx(1))
  
end
legend

figure(7)
hold on
for ii=28:36
  [C,ia,idx] = unique(WavelengthTime{ii},'stable');
  mean_signal = accumarray(idx,DetectorSignal{ii},[],@mean); 
  %semilogy(WavelengthTime{ii},DetectorSignal{ii}, 'b+', ...
  semilogy(C,mean_signal, 'b+', ...   
    'DisplayName', ['AOI=', num2str(AOI{ii}), num2str(Expansion{ii}), num2str(FiberType{ii})])
  hold off
end
legend

figure(8)
hold on
for ii=37:48
  [C,ia,idx] = unique(WavelengthTime{ii},'stable');
  mean_signal = accumarray(idx,DetectorSignal{ii},[],@mean);    
  %semilogy(WavelengthTime{ii},DetectorSignal{ii}, 'g+', ...
  semilogy(C,mean_signal, 'g+', ...     
    'DisplayName', ['AOI=', num2str(AOI{ii}), num2str(Expansion{ii}), num2str(FiberType{ii})])
  hold off
end
legend

figure(9) %4x SM
hold on
for ii=49:55
  [C,ia,idx] = unique(WavelengthTime{ii},'stable');
  mean_signal = accumarray(idx,DetectorSignal{ii},[],@mean); 
 % semilogy(WavelengthTime{ii},DetectorSignal{ii}, 'r+', ...
  semilogy(C,mean_signal, 'r+', ... 
    'DisplayName', ['AOI=', num2str(AOI{ii}), num2str(Expansion{ii}), num2str(FiberType{ii})])
 hold off
end
legend

figure(10)
hold on
for ii=56:68
  [C,ia,idx] = unique(WavelengthTime{ii},'stable');
   mean_signal = accumarray(idx,DetectorSignal{ii},[],@mean); 
  % semilogy(WavelengthTime{ii},DetectorSignal{ii}, 'k+', ...
   semilogy(C,mean_signal, 'k+', ... 
    'DisplayName', ['AOI=', num2str(AOI{ii}), num2str(Expansion{ii}), num2str(FiberType{ii})])
  hold off
end
legend

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

  FigH = figure(3);
%   set(gca,'Fontsize',36,'Fontweight','b'); 
   set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1  1  1600 1200]);
   name=strcat('Etalon_V02_data'); 
   print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
 


  