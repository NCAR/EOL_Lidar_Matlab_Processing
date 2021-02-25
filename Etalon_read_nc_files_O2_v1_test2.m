clear all; close all

serv_path = '/Users/spuler/Desktop';
cd(strcat(serv_path,'/mpd/V01_test_2')) 
d = pwd;

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
    EtalonTemperature{jj} =  ncreadatt(filename, '/', 'EtalonTemperature');
    
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
  % AO1=0.0, 2x, MM, T=20.0, FWHM=9.5pm, FSR=0.3122nm
  ind_start = 1;   ind_end = 1;
  step = 0.0005 %spline step size in nm
  cd('/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing')
  [C,mean_signal,FWHM, FSR, s, xx,locs, loc, val] = etlaon_fit(ind_start,ind_end, WavelengthTime,DetectorSignal, step);
  cd(d)
  hold on
  plot(C,mean_signal, 'ko', ...      
     'DisplayName',    [num2str(AOI{ind_start}), char(176), 'AOI, ', num2str(Expansion{ind_start}), ', ', num2str(FiberType{ind_start}), ', T=',  num2str(EtalonTemperature{ind_start}), char(176),'C, ', num2str(FWHM,2),'pm FWHM, ', num2str(FSR,4),'nm FSR'])
  %hold on
  %  plot(xx, s, 'r') 
  %  plot(xx(locs(loc)),val,'rv', 'MarkerFaceColor', 'r');
  FWHM
  FSR
  legend
  
  % AO1=0.0, 2x, SM, T=20.0, FWHM=8.2pm, FSR=0.3135
  ind_start = 2;   ind_end = 3;
  step = 0.00075 %spline step size in nm
  cd('/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing')
  [C,mean_signal,FWHM, FSR, s, xx,locs, loc, val] = etlaon_fit(ind_start,ind_end, WavelengthTime,DetectorSignal, step);
  cd(d)
  hold on
  plot(C,mean_signal, 'ro', ...      
     'DisplayName',    [num2str(AOI{ind_start}), char(176), 'AOI, ', num2str(Expansion{ind_start}), ', ', num2str(FiberType{ind_start}), ', T=',  num2str(EtalonTemperature{ind_start}), char(176),'C, ', num2str(FWHM,2),'pm FWHM, ', num2str(FSR,4),'nm FSR'])
  %hold on
  %  plot(xx, s, 'r') 
  %  plot(xx(locs(loc)),val,'rv', 'MarkerFaceColor', 'r');
  FWHM
  FSR
  legend

  % AO1=0.0, 4x, MM, T=20.0, FWHM=9.0pm, FSR=0.3127
  ind_start = 4;   ind_end = 6;
  step = 0.00075 %spline step size in nm
  cd('/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing')
  [C,mean_signal,FWHM, FSR, s, xx,locs, loc, val] = etlaon_fit(ind_start,ind_end, WavelengthTime,DetectorSignal, step);
  cd(d)
  hold on
  plot(C,mean_signal, 'go', ...      
     'DisplayName',    [num2str(AOI{ind_start}), char(176), 'AOI, ', num2str(Expansion{ind_start}), ', ', num2str(FiberType{ind_start}), ', T=',  num2str(EtalonTemperature{ind_start}), char(176),'C, ', num2str(FWHM,2),'pm FWHM, ', num2str(FSR,4),'nm FSR'])
  %hold on
  %  plot(xx, s, 'r') 
  %  plot(xx(locs(loc)),val,'rv', 'MarkerFaceColor', 'r');
  FWHM
  FSR
  legend
  
  % AO1=2.1, 4x, MM, T=20.0, FWHM=14.5pm, FSR=0.3120
  ind_start = 7;   ind_end = 7;
  step = 0.0005 %spline step size in nm
  cd('/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing')
  [C,mean_signal,FWHM, FSR, s, xx,locs, loc, val] = etlaon_fit(ind_start,ind_end, WavelengthTime,DetectorSignal, step);
  cd(d)
  hold on
  plot(C,mean_signal, 'bo', ...      
     'DisplayName',    [num2str(AOI{ind_start}), char(176), 'AOI, ', num2str(Expansion{ind_start}), ', ', num2str(FiberType{ind_start}), ', T=',  num2str(EtalonTemperature{ind_start}), char(176),'C, ', num2str(FWHM,2),'pm FWHM, ', num2str(FSR,4),'nm FSR'])
  %hold on
  %  plot(xx, s, 'r') 
  %  plot(xx(locs(loc)),val,'rv', 'MarkerFaceColor', 'r');
  FWHM
  FSR
  legend
  
  % AO1=2.1, 2x, MM, T=20.0, FWHM=17.3pm, FSR=0.3120
  ind_start = 8;   ind_end = 8;
  step = 0.00075 %spline step size in nm
  cd('/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing')
  [C,mean_signal,FWHM, FSR, s, xx,locs, loc, val] = etlaon_fit(ind_start,ind_end, WavelengthTime,DetectorSignal, step);
  cd(d)
  hold on
  plot(C,mean_signal, 'mo', ...      
     'DisplayName',    [num2str(AOI{ind_start}), char(176), 'AOI, ', num2str(Expansion{ind_start}), ', ', num2str(FiberType{ind_start}), ', T=',  num2str(EtalonTemperature{ind_start}), char(176),'C, ', num2str(FWHM,2),'pm FWHM, ', num2str(FSR,4),'nm FSR'])
  %hold on
  %  plot(xx, s, 'r') 
  %  plot(xx(locs(loc)),val,'rv', 'MarkerFaceColor', 'r');
  FWHM
  FSR
  legend

  % AO1=2.1, 2x, MM, T=20.0, FWHM=16.5pm, FSR=0.3120nm  ** overfilled launch
  ind_start = 9;   ind_end = 9;
  step = 0.00075 %spline step size in nm
  cd('/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing')
  [C,mean_signal,FWHM, FSR, s, xx,locs, loc, val] = etlaon_fit(ind_start,ind_end, WavelengthTime,DetectorSignal, step);
  cd(d)
  hold on
  plot(C,mean_signal, 'co', ...      
     'DisplayName',    [num2str(AOI{ind_start}), char(176), 'AOI, ', num2str(Expansion{ind_start}), ', ', num2str(FiberType{ind_start}), ', T=',  num2str(EtalonTemperature{ind_start}), char(176),'C, ', num2str(FWHM,2),'pm FWHM, ', num2str(FSR,4),'nm FSR'])
  %hold on
  %  plot(xx, s, 'r') 
  %  plot(xx(locs(loc)),val,'rv', 'MarkerFaceColor', 'r');
  FWHM
  FSR
  legend
  
  % AO1=2.1, 4x, MM, T=20.0, FWHM=13.9pm, FSR=0.3120nm ** overfilled launch
  ind_start = 10;   ind_end = 10;
  step = 0.0005 %spline step size in nm
  cd('/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing')
  [C,mean_signal,FWHM, FSR, s, xx,locs, loc, val] = etlaon_fit(ind_start,ind_end, WavelengthTime,DetectorSignal, step);
  cd(d)
  hold on
  plot(C,mean_signal, 'yo', ...      
     'DisplayName',    [num2str(AOI{ind_start}), char(176), 'AOI, ', num2str(Expansion{ind_start}), ', ', num2str(FiberType{ind_start}), ', T=',  num2str(EtalonTemperature{ind_start}), char(176),'C, ', num2str(FWHM,2),'pm FWHM, ', num2str(FSR,4),'nm FSR'])
  %hold on
  %  plot(xx, s, 'r') 
  %  plot(xx(locs(loc)),val,'rv', 'MarkerFaceColor', 'r');
  FWHM
  FSR
  legend
  
  
grid on
title('Etalon V01')
ylabel('rel. trans. power'); 
xlabel('wavelength (nm)'); 
xlim([769.72  770.2])
ylim([1e-6 4e-4])
set(gca,'Fontsize',20,'Fontweight','b'); %

 xlim([769.65  770.25])
 ylim([1e-6 3e-4])

 
  

  FigH = figure(3);
%   set(gca,'Fontsize',36,'Fontweight','b'); 
   set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1  1  1600 1200]);
   name=strcat('Etalon_V01_test2'); 
   print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
 
%   FigH = figure(4);
% %   set(gca,'Fontsize',36,'Fontweight','b'); 
%    set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1  1  1600 1200]);
%    name=strcat('Etalon_V01_combined'); 
%    print(FigH, name, '-dpng', '-r0') % set at the screen resolution

  