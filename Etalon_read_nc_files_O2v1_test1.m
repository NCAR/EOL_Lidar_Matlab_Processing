clear all; close all

serv_path = '/Users/spuler/Desktop';
cd(strcat(serv_path,'/mpd/Etalon/O2_V01')) %test_1
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
  % AO1=2.25, 4x, SM, FWHM=8.7pm, FSR=0.3120nm
  ind_start = 1;   ind_end = 4;
  step = 0.0003 %spline step size in nm
  cd('/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing')
  [C,mean_signal,FWHM, FSR, s, xx,locs, loc, val] = etlaon_fit(ind_start,ind_end, WavelengthTime,DetectorSignal, step);
  cd(d)
  hold on
  plot(C,mean_signal, 'ko', ...      
      'DisplayName',    [num2str(AOI{ind_start}), char(176), 'AOI, ', num2str(Expansion{ind_start}), ', ', num2str(FiberType{ind_start}), ', T=24.6', char(176),'C, ', num2str(FWHM,2),'pm FWHM, ', num2str(FSR,4),'pm FSR'])
%   hold on
%     plot(xx, s, 'r') 
%     plot(xx(locs(loc)),val,'rv', 'MarkerFaceColor', 'r');
  FWHM
  FSR
  legend

  % AO1=2.25, 4x, MM, FWHM=14.5pm, FSR=0.3126nm
  ind_start = 5;   ind_end = 8;
  step = 0.0005 %spline step size in nm
  cd('/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing')
  [C,mean_signal,FWHM, FSR, s, xx,locs, loc, val] = etlaon_fit(ind_start,ind_end, WavelengthTime,DetectorSignal, step);
  cd(d)
  hold on
  plot(C,mean_signal, 'ro', ...      
     'DisplayName',    [num2str(AOI{ind_start}), char(176), 'AOI, ', num2str(Expansion{ind_start}), ', ', num2str(FiberType{ind_start}), ', T=24.6', char(176),'C, ', num2str(FWHM,2),'pm FWHM, ', num2str(FSR,4),'pm FSR'])
  %hold on
  %  plot(xx, s, 'r') 
  %  plot(xx(locs(loc)),val,'rv', 'MarkerFaceColor', 'r');
  FWHM
  FSR
  legend

  % AO1=2.25, 2x, MM, FWHM=18.0pm, FSR=0.3127nm 
  hold on
  ind_start = 9;   ind_end = 12;
  step = 0.00075 %spline step size in nm
  cd('/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing')
  [C,mean_signal,FWHM, FSR, s, xx,locs, loc, val] = etlaon_fit(ind_start,ind_end, WavelengthTime,DetectorSignal, step);
  cd(d)
  hold on
  plot(C,mean_signal, 'go', ...      
      'DisplayName',    [num2str(AOI{ind_start}), char(176), 'AOI, ', num2str(Expansion{ind_start}), ', ', num2str(FiberType{ind_start}), ', T=24.6', char(176),'C, ', num2str(FWHM,2),'pm FWHM, ', num2str(FSR,4),'pm FSR'])
  %hold on
  %  plot(xx, s, 'r') 
  %  plot(xx(locs(loc)),val,'rv', 'MarkerFaceColor', 'r');
  FWHM
  FSR
  legend
  
  % AO1=2.25, 2x, SM, FWHM=10.5pm, FSR=0.3124nm
  hold on
  ind_start = 13;   ind_end = 19;
  step = 0.00025 %spline step size in nm
  cd('/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing')
  [C,mean_signal,FWHM, FSR, s, xx,locs, loc, val] = etlaon_fit(ind_start,ind_end, WavelengthTime,DetectorSignal, step);
  cd(d)
  hold on
  plot(C,mean_signal, 'bo', ...      
      'DisplayName',    [num2str(AOI{ind_start}), char(176), 'AOI, ', num2str(Expansion{ind_start}), ', ', num2str(FiberType{ind_start}), ', T=24.6', char(176),'C, ', num2str(FWHM,2),'pm FWHM, ', num2str(FSR,4),'pm FSR'])
%   hold on
%     plot(xx, s, 'r') 
%     plot(xx(locs(loc)),val,'rv', 'MarkerFaceColor', 'r');
  FWHM
  FSR
  legend
  

  % AO1=0.0, 4x, MM, FWHM=9.0pm, FSR=0.3120nm
  ind_start = 20;   ind_end = 22;
  step = 0.00025 %spline step size in nm
  cd('/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing')
  [C,mean_signal,FWHM, FSR, s, xx,locs, loc, val] = etlaon_fit(ind_start,ind_end, WavelengthTime,DetectorSignal, step);
  cd(d)
  hold on
  plot(C,mean_signal, 'mo', ...      
      'DisplayName',    [num2str(AOI{ind_start}), char(176), 'AOI, ', num2str(Expansion{ind_start}), ', ', num2str(FiberType{ind_start}), ', T=24.6', char(176),'C, ', num2str(FWHM(1),2),'pm FWHM, ', num2str(FSR,4),'pm FSR'])
  %hold on
  %  plot(xx, s, 'r') 
  %  plot(xx(locs(loc)),val,'rv', 'MarkerFaceColor', 'r');
  FWHM(1)
  FSR
  legend

 % AO1=0.0, 4x, SM, FWHM=9.0pm, FSR=0.3125nm
  ind_start = 23;   ind_end = 26;
  step = 0.00025 %spline step size in nm
  cd('/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing')
  [C,mean_signal,FWHM, FSR, s, xx,locs, loc, val] = etlaon_fit(ind_start,ind_end, WavelengthTime,DetectorSignal, step);
  cd(d)
  hold on
  plot(C,mean_signal, 'co', ...      
      'DisplayName',    [num2str(AOI{ind_start}), char(176), 'AOI, ', num2str(Expansion{ind_start}), ', ', num2str(FiberType{ind_start}), ', T=24.6', char(176),'C, ', num2str(FWHM,2),'pm FWHM, ', num2str(FSR,4),'pm FSR'])
%   hold on
%    plot(xx, s, 'r') 
%    plot(xx(locs(loc)),val,'rv', 'MarkerFaceColor', 'r');
  FWHM
  FSR
  legend
  
  % AO1=0.0, 2x, SM, FWHM=10.0pm, FSR=0.3121 nm
  ind_start = 27;   ind_end = 32;
  step = 0.0005 %spline step size in nm
  cd('/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing')
  [C,mean_signal,FWHM, FSR, s, xx,locs, loc, val] = etlaon_fit(ind_start,ind_end, WavelengthTime,DetectorSignal, step);
  cd(d)
  hold on
  plot(C,mean_signal, 'yo', ...      
      'DisplayName',    [num2str(AOI{ind_start}), char(176), 'AOI, ', num2str(Expansion{ind_start}), ', ', num2str(FiberType{ind_start}), ', T=24.6', char(176),'C, ', num2str(FWHM(1),2),'pm FWHM, ', num2str(FSR,4),'pm FSR'])
%   %hold on
%     plot(xx, s, 'r') 
%     plot(xx(locs(loc)),val,'rv', 'MarkerFaceColor', 'r');
  FWHM(1)
  FSR
  legend

  % AO1=0.0, 2x, MM, FWHM=9.5pm, FSR=0.3120nm  
  ind_start = 33;   ind_end = 38;
  step = 0.0005 %spline step size in nm
  cd('/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_Processing')
  [C,mean_signal,FWHM, FSR, s, xx,locs, loc, val] = etlaon_fit(ind_start,ind_end, WavelengthTime,DetectorSignal, step);
  cd(d)
  hold on
  plot(C,mean_signal, 'k+', ...      
     'DisplayName',    [num2str(AOI{ind_start}), char(176), 'AOI, ', num2str(Expansion{ind_start}), ', ', num2str(FiberType{ind_start}), ', T=24.6', char(176),'C, ', num2str(FWHM,2),'pm FWHM, ', num2str(FSR,4),'pm FSR'])
%   hold on
%    plot(xx, s, 'b') 
%    plot(xx(locs(loc)),val,'rv', 'MarkerFaceColor', 'r');
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

%ylim([1e-5 1e-3])
%ylim([0 4e-4])
%xlim([769.74  769.79])

  FigH = figure(3);
%   set(gca,'Fontsize',36,'Fontweight','b'); 
   set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1  1  1600 1200]);
   name=strcat('O2_Etalon_V01_test1'); 
   print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
 


  