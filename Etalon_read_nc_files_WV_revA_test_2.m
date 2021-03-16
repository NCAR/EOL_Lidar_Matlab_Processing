clear all; close all

serv_path = '/Users/spuler/Desktop';
cd(strcat(serv_path,'/mpd/Etalon/WV_revA_test_2')) 
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
  % AO1=0.0, 4x, MM, T=27.0, FWHM=2.5pm, FSR=99.5pmm
  ind_start = 1;   ind_end = 1;
  step = 0.00025 %spline step size in nm
  cd('/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing')
  [C,mean_signal,FWHM, FSR, s, xx,locs, loc, val] = etlaon_fit(ind_start,ind_end, WavelengthTime,DetectorSignal, step);
  cd(d)
  plot(C,mean_signal, 'ko', ...      
     'DisplayName',    [num2str(AOI{ind_start}), char(176), 'AOI, ', num2str(Expansion{ind_start}), ', ', num2str(FiberType{ind_start}), ', T=',  num2str(EtalonTemperature{ind_start}), char(176),'C, ', num2str(FWHM,4),'pm FWHM, ', num2str(FSR,4),'pm FSR'])
%   hold on
%     plot(xx, s, 'r') 
%     plot(xx(locs(loc)),val,'rv', 'MarkerFaceColor', 'r');
  FWHM
  FSR
  legend
  
  % AO1=0.0, 4x, SM, T=27.0, FWHM=2.3pm, FSR=99.3pm
  ind_start = 2;   ind_end = 2;
  step = 0.00025 %spline step size in nm
  cd('/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing')
  [C,mean_signal,FWHM, FSR, s, xx,locs, loc, val] = etlaon_fit(ind_start,ind_end, WavelengthTime,DetectorSignal, step);
  cd(d)
  hold on
  plot(C,mean_signal, 'ro', ...      
     'DisplayName',    [num2str(AOI{ind_start}), char(176), 'AOI, ', num2str(Expansion{ind_start}), ', ', num2str(FiberType{ind_start}), ', T=',  num2str(EtalonTemperature{ind_start}), char(176),'C, ', num2str(FWHM,4),'pm FWHM, ', num2str(FSR,4),'pm FSR'])
%   hold on
%    plot(xx, s, 'k') 
%    plot(xx(locs(loc)),val,'rv', 'MarkerFaceColor', 'k');
  FWHM
  FSR
  legend

  % AO1=0.0, 2x, SM, T=27.0, FWHM=1.9pm, FSR=99.5pm 
  ind_start = 3;   ind_end = 3;
  step = 0.00025 %spline step size in nm
  cd('/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing')
  [C,mean_signal,FWHM, FSR, s, xx,locs, loc, val] = etlaon_fit(ind_start,ind_end, WavelengthTime,DetectorSignal, step);
  cd(d)
  hold on
  plot(C,mean_signal, 'go', ...      
     'DisplayName',    [num2str(AOI{ind_start}), char(176), 'AOI, ', num2str(Expansion{ind_start}), ', ', num2str(FiberType{ind_start}), ', T=',  num2str(EtalonTemperature{ind_start}), char(176),'C, ', num2str(FWHM,4),'pm FWHM, ', num2str(FSR,4),'pm FSR'])
  %hold on
  %  plot(xx, s, 'r') 
  %  plot(xx(locs(loc)),val,'rv', 'MarkerFaceColor', 'r');
  FWHM
  FSR
  legend
  
  % AO1=0.0, 2x, MM, T=27.0, FWHM=2.3pm, FSR=99.6pm 
  ind_start = 4;   ind_end = 4;
  step = 0.00025 %spline step size in nm
  cd('/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing')
  [C,mean_signal,FWHM, FSR, s, xx,locs, loc, val] = etlaon_fit(ind_start,ind_end, WavelengthTime,DetectorSignal, step);
  cd(d)
  hold on
  plot(C,mean_signal, 'bo', ...      
     'DisplayName',    [num2str(AOI{ind_start}), char(176), 'AOI, ', num2str(Expansion{ind_start}), ', ', num2str(FiberType{ind_start}), ', T=',  num2str(EtalonTemperature{ind_start}), char(176),'C, ', num2str(FWHM,4),'pm FWHM, ', num2str(FSR,4),'pm FSR'])
  %hold on
  %  plot(xx, s, 'r') 
  %  plot(xx(locs(loc)),val,'rv', 'MarkerFaceColor', 'r');
  FWHM
  FSR
  legend
  
  % AO1=0.0, 2x, MM, T=27.0, FWHM=2.0pm, FSR=99.5pm ** 50%Overfilled
  ind_start = 5;   ind_end = 5;
  step = 0.00025; %spline step size in nm
  cd('/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing')
  [C,mean_signal,FWHM, FSR, s, xx,locs, loc, val] = etlaon_fit(ind_start,ind_end, WavelengthTime,DetectorSignal, step);
  cd(d)
  hold on
  plot(C,mean_signal, 'mo', ...      
     'DisplayName',    [num2str(AOI{ind_start}), char(176), 'AOI, ', num2str(Expansion{ind_start}), ', ', num2str(FiberType{ind_start}), ', T=',  num2str(EtalonTemperature{ind_start}), char(176),'C, ', num2str(FWHM,4),'pm FWHM, ', num2str(FSR,4),'pm FSR'])
%   hold on
%    plot(xx, s, 'k') 
%    plot(xx(locs(loc)),val,'rv', 'MarkerFaceColor', 'k');
  FWHM
  FSR
  legend

  % AO1=0.0, 4x, MM, T=27.0, FWHM=2.4pm, FSR=99.5pm ** 50%Overfilled
  ind_start = 6;   ind_end = 6;
  step = 0.00025; %spline step size in nm
  cd('/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing')
  [C,mean_signal,FWHM, FSR, s, xx,locs, loc, val] = etlaon_fit(ind_start,ind_end, WavelengthTime,DetectorSignal, step);
  cd(d)
  hold on
  plot(C,mean_signal, 'co', ...      
     'DisplayName',    [num2str(AOI{ind_start}), char(176), 'AOI, ', num2str(Expansion{ind_start}), ', ', num2str(FiberType{ind_start}), ', T=',  num2str(EtalonTemperature{ind_start}), char(176),'C, ', num2str(FWHM,4),'pm FWHM, ', num2str(FSR,4),'pm FSR'])
  %hold on
  %  plot(xx, s, 'r') 
  %  plot(xx(locs(loc)),val,'rv', 'MarkerFaceColor', 'r');
  FWHM
  FSR
  legend
  
  % AO1=0.0, 4x, MM, T=20.0, FWHM=2.3pm, FSR=99.5pm 
  ind_start = 7;   ind_end = 7;
  step = 0.00025 %spline step size in nm
  cd('/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing')
  [C,mean_signal,FWHM, FSR, s, xx,locs, loc, val] = etlaon_fit(ind_start,ind_end, WavelengthTime,DetectorSignal, step);
  cd(d)
  hold on
  plot(C,mean_signal, 'yo', ...      
     'DisplayName',    [num2str(AOI{ind_start}), char(176), 'AOI, ', num2str(Expansion{ind_start}), ', ', num2str(FiberType{ind_start}), ', T=',  num2str(EtalonTemperature{ind_start}), char(176),'C, ', num2str(FWHM,4),'pm FWHM, ', num2str(FSR,4),'pm FSR'])
  %hold on
  %  plot(xx, s, 'r') 
  %  plot(xx(locs(loc)),val,'rv', 'MarkerFaceColor', 'r');
  FWHM
  FSR
  legend
  
  
  % AO1=0.0, 4x, MM, T=20.0, FWHM=2.0pm, FSR=99.3pm **50% Overfilled
  ind_start = 8;   ind_end = 8;
  step = 0.00025 %spline step size in nm
  cd('/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing')
  [C,mean_signal,FWHM, FSR, s, xx,locs, loc, val] = etlaon_fit(ind_start,ind_end, WavelengthTime,DetectorSignal, step);
  cd(d)
  hold on
  plot(C,mean_signal, 'k+', ...      
     'DisplayName',    [num2str(AOI{ind_start}), char(176), 'AOI, ', num2str(Expansion{ind_start}), ', ', num2str(FiberType{ind_start}), ', T=',  num2str(EtalonTemperature{ind_start}), char(176),'C, ', num2str(FWHM,4),'pm FWHM, ', num2str(FSR,4),'pm FSR'])
  %hold on
  %  plot(xx, s, 'r') 
  %  plot(xx(locs(loc)),val,'rv', 'MarkerFaceColor', 'r');
  FWHM
  FSR
  legend
  
  % AO1=0.0, 4x, SM, T=20.0, FWHM=2.0pm, FSR=99.4pm 
  ind_start = 9;   ind_end = 9;
  step = 0.00025 %spline step size in nm
  cd('/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing')
  [C,mean_signal,FWHM, FSR, s, xx,locs, loc, val] = etlaon_fit(ind_start,ind_end, WavelengthTime,DetectorSignal, step);
  cd(d)
  hold on
  plot(C,mean_signal, 'r+', ...      
     'DisplayName',    [num2str(AOI{ind_start}), char(176), 'AOI, ', num2str(Expansion{ind_start}), ', ', num2str(FiberType{ind_start}), ', T=',  num2str(EtalonTemperature{ind_start}), char(176),'C, ', num2str(FWHM,4),'pm FWHM, ', num2str(FSR,4),'pm FSR'])
  %hold on
  %  plot(xx, s, 'r') 
  %  plot(xx(locs(loc)),val,'rv', 'MarkerFaceColor', 'r');
  FWHM
  FSR
  legend
  
  % AO1=0.0, 2x, SM, T=20.0, FWHM=2.0pm, FSR=99.4pm   
  ind_start = 10;   ind_end = 10;
  step = 0.00025 %spline step size in nm
  cd('/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing')
  [C,mean_signal,FWHM, FSR, s, xx,locs, loc, val] = etlaon_fit(ind_start,ind_end, WavelengthTime,DetectorSignal, step);
  cd(d)
  hold on
  plot(C,mean_signal, 'g+', ...      
     'DisplayName',    [num2str(AOI{ind_start}), char(176), 'AOI, ', num2str(Expansion{ind_start}), ', ', num2str(FiberType{ind_start}), ', T=',  num2str(EtalonTemperature{ind_start}), char(176),'C, ', num2str(FWHM,4),'pm FWHM, ', num2str(FSR,4),'pm FSR'])
  %hold on
  %  plot(xx, s, 'r') 
  %  plot(xx(locs(loc)),val,'rv', 'MarkerFaceColor', 'r');
  FWHM
  FSR
  legend
  
  % AO1=0.0, 2x, MM, T=20.0, FWHM=2.0pm, FSR=99.2pm   
  ind_start = 11;   ind_end = 11;
  step = 0.00025 %spline step size in nm
  cd('/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing')
  [C,mean_signal,FWHM, FSR, s, xx,locs, loc, val] = etlaon_fit(ind_start,ind_end, WavelengthTime,DetectorSignal, step);
  cd(d)
  hold on
  plot(C,mean_signal, 'b+', ...      
     'DisplayName',    [num2str(AOI{ind_start}), char(176), 'AOI, ', num2str(Expansion{ind_start}), ', ', num2str(FiberType{ind_start}), ', T=',  num2str(EtalonTemperature{ind_start}), char(176),'C, ', num2str(FWHM,4),'pm FWHM, ', num2str(FSR,4),'pm FSR'])
  %hold on
  %  plot(xx, s, 'r') 
  %  plot(xx(locs(loc)),val,'rv', 'MarkerFaceColor', 'r');
  FWHM
  FSR
  legend
  
  % AO1=0.0, 2x, MM, T=20.0, FWHM=2.0pm, FSR=99.3pm  **50% overfilled 
  ind_start = 12;   ind_end = 12;
  step = 0.00025 %spline step size in nm
  cd('/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing')
  [C,mean_signal,FWHM, FSR, s, xx,locs, loc, val] = etlaon_fit(ind_start,ind_end, WavelengthTime,DetectorSignal, step);
  cd(d)
  hold on
  plot(C,mean_signal, 'm+', ...      
     'DisplayName',    [num2str(AOI{ind_start}), char(176), 'AOI, ', num2str(Expansion{ind_start}), ', ', num2str(FiberType{ind_start}), ', T=',  num2str(EtalonTemperature{ind_start}), char(176),'C, ', num2str(FWHM,4),'pm FWHM, ', num2str(FSR,4),'pm FSR'])
  %hold on
  %  plot(xx, s, 'r') 
  %  plot(xx(locs(loc)),val,'rv', 'MarkerFaceColor', 'r');
  FWHM
  FSR
  legend
  
  
grid on
title('WV Etalon SN OP-10275')
ylabel('rel. trans. power'); 
xlabel('wavelength (nm)'); 
%xlim([769.72  770.2])
ylim([1e-6 9e-4])
set(gca,'Fontsize',20,'Fontweight','b'); %

xlim([828.14  828.32])
% ylim([1e-6 3e-4])

 
  

  FigH = figure(3);
%   set(gca,'Fontsize',36,'Fontweight','b'); 
   set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1  1  1600 1200]);
   name=strcat('WV_Etalon_revA_test2'); 
   print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
 
%   FigH = figure(4);
% %   set(gca,'Fontsize',36,'Fontweight','b'); 
%    set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1  1  1600 1200]);
%    name=strcat('Etalon_V01_combined'); 
%    print(FigH, name, '-dpng', '-r0') % set at the screen resolution

  