%clear all; close all;
g = 9.8076; %(m/s) 
H_v = 2501000; % (J/kg)
R_sd = 287; %(J/kg K)
R_sw = 461.5; %(J/kg K)
c_pd = 1003.5; %(J/kg K)
m_air = 28.9647; m_wv = 18; 
R = [0 1 2 3 4 5 6 7 8];  %(km) range

 model= 'US Standard Atmosphere';
r_vmr = [7.75E-3 6.07E-3 4.63E-3 3.18E-3 2.16E-3 1.4E-3 9.25E-4 5.72E-4 3.67E-4]; %water vapor volume mass ratio
%p = [1 0.887 0.785 0.692 0.609 0.533 0.466 0.406 0.352]; %(atm) pressure
T_r = [288.2 281.7 275.2 268.7 262.7 255.7 249.2 242.7 236.2]; %(K) temperature

%  model=  'Midlatidue Summer Atmosphere';
% r_vmr = [1.880000e-02 1.380000e-02 9.680000e-03 5.980000e-03 3.810000e-03 2.230000e-03 1.510000e-03 1.020000e-03 6.460000e-04]; %water vapor volume mass ratio
% %p = [1.013e+03 9.020e+02 8.020e+02 7.100e+02 6.280e+02 5.540e+02 4.870e+02 4.260e+02 3.720e+02]; %(mb) pressure
% T_r = [294.200 289.700 285.200 279.200 273.200 267.200 261.200 254.700 248.200]; %(K) temperature
% 
%  model= 'Midlatidue Winter Atmosphere';
%  r_vmr = [4.320000e-03 3.450000e-03 2.790000e-03 2.090000e-03 1.280000e-03 8.240000e-04 5.100000e-04 2.320000e-04 1.080000e-04]; %water vapor volume mass ratio
%  %p = [  1.018000e+03  8.973000e+02 7.897000e+02 6.938000e+02 6.081000e+02 5.313000e+02 4.627000e+02 4.016000e+02 3.473000e+02]; %(mb) pressure
%  T_r = [272.200 268.700 265.200 261.700 255.700 249.700 243.700 237.700 231.700]; %(K) temperature
 
 r_mr = r_vmr.*(m_air/m_wv); % convert to water vapor mass ratio
% r_mr_ave = r_mr;
% for j=1:(size(r_vmr,2)-1)  
%    r_mr_ave(j) = (2/3.*r_mr(j)+ 1/3.*r_mr(j+1))./2;
% end
% r_mr = r_mr_ave;


% read radiosondes data
step = 0.075 ; 
R = 0:step:8;
sonde_end_int = 30; % time to integrate the MPD data in comparison with sonde profile
cd('/Volumes/fog1/rsfdata/MPD/mpd_ancillary_data/radiosondes/CSU_PrePRECIP')
plot_path = '/Volumes/Macintosh HD/Users/spuler/Desktop/mpd/Plots';
[sondefilename, sondedir] = uigetfile('*.*','Select the sonde file', 'MultiSelect', 'on');
jj=1;
elevation= 1574;

for jj = 1:size(sondefilename,2)
   cd('/Users/spuler/Documents/GitHub/EOL_Lidar_Matlab_processing/')
   [sonde_MR, sonde_T, sonde_P, range_grid, sonde_date, sonde_time, duration_sonde] = Sonde_read_CSU_files_basic(jj, elevation, sondedir, sondefilename, step); 

%  model =  strcat(sonde_date, ' ', sonde_time, ' Sonde') 
  T_r = sonde_T+273.15; %convert sonde surface temperature from C to K
  r_mr = sonde_MR/1000; % convert sonde mixing ratio in g/kg to kg/gk
  L_ad = 9.8; %(K/km) dry adiabtaic lapse rate 
  T_ad = T_r(1)-L_ad*R  %T profile using adiabatic lapse rate


 % use MPD to measure 
  model =  strcat(sonde_date, ' ', sonde_time, ' MPD'); 
 % select the MPD AH, and surface T and P profiles
  [minValue, closestIndex] = min(abs(min(duration_sonde-sonde_end_int/24/60)-duration))
  [minValue, closestIndex_end] = min(abs(min(duration_sonde+sonde_end_int/24/60)-duration))
  MPD_AH = mean(N_avg_comb(closestIndex:closestIndex_end,:),1, 'omitnan').*1e6./6.022E23.*18.015;
  T_surface = mean(surf_T(closestIndex:closestIndex_end,:),1, 'omitnan')+273.15;
  P_surface = mean(surf_P(closestIndex:closestIndex_end,:),1, 'omitnan');
%  % convert AH to mixing ratio
   T_ad = T_surface-L_ad*R;  %T profile using adiabatic lapse rate   
   T_st = T_surface-6.5*R; % estimate T profile using suface station and standard atmosphereic lappse rate 
   P_est = P_surface.*(T_surface./T_st).^-5.2199;   % estimate pressure (atm) using suface station and hydrostatic 
   density = P_est./(8.31447E-2*T_st)*28.97;
   r_mr = MPD_AH(1:length(density))./density./1000;  % convert AH(g/m^3) to MR (kg/kg)
 
   figure(101)
   plot(T_st,R)
   hold on
   plot(sonde_T+273.15,R) %temp convert (C) to (K)
   hold off
   
   figure(102)
   plot(P_est,R)
   hold on
   plot(sonde_P*.000986923,R) %pressure convert (mb) to (atm)
   hold off
   
   figure(103)
   plot(r_mr,R)
   hold on
   plot(sonde_MR/1000,R) %mixing ratio convert (g/kg) to (kg/kg)
   hold off



T_0 = T_ad; %initial temperature guess


L_m = g.*((1+(H_v.*r_mr)./(R_sd.*T_0))./(c_pd+(H_v.^2.*r_mr)./(R_sw.*T_0.^2))).*1000;  %(K/km) moisture corrected lapse rate   
L_m(isnan(L_m))=L_ad; %set the nans to the standard atm lapse rate
T_1 = T_0; 
for i = 2:size(r_mr,2)
  T_1(i) = T_1(i-1)-L_m(i-1)/(1/step)
end

L_m = g.*((1+(H_v.*r_mr)./(R_sd.*T_1))./(c_pd+(H_v.^2.*r_mr)./(R_sw.*T_1.^2))).*1000; %(K/km) moisture corrected lapse rate   
L_m(isnan(L_m))=L_ad; %set the nans to the standard atm lapse rate
T_2 = T_1; 
for i = 2:size(r_mr,2)
  T_2(i) = T_2(i-1)-L_m(i-1)/(1/step)
end

L_m = g.*((1+(H_v.*r_mr)./(R_sd.*T_2))./(c_pd+(H_v.^2.*r_mr)./(R_sw.*T_2.^2))).*1000; %(K/km) moisture corrected lapse rate   
L_m(isnan(L_m))=L_ad; %set the nans to the standard atm lapse rate
T_3 = T_2; 
for i = 2:size(r_mr,2)
  T_3(i) = T_3(i-1)-L_m(i-1)/(1/step)
end

L_m = g.*((1+(H_v.*r_mr)./(R_sd.*T_3))./(c_pd+(H_v.^2.*r_mr)./(R_sw.*T_3.^2))).*1000; %(K/km) moisture corrected lapse rate   
L_m(isnan(L_m))=L_ad; %set the nans to the standard atm lapse rate
T_4 = T_3; 
for i = 2:size(r_mr,2)
  T_4(i) = T_4(i-1)-L_m(i-1)/(1/step)
end

figure(201)
plot(T_0,R,'k')
hold on
plot(T_1,R, 'r', 'LineWidth',2)
plot(T_2,R, 'g','LineWidth',2 )
%plot(T_3,R , 'b', 'LineWidth',2)
%plot(T_4,R)
plot(T_r,R, 'ko')
hold off
ylim([0 8])
xlim([180 300])
ylabel('range (R)')
xlabel('Temperature (K)')
title(model)
 legend('dry adiabatic T_0', 'moist corrected T_1', 'moist corrected T_2', 'actual model T', 'Location', 'SouthWest');

figure(202)
plot(T_r-T_2,R, 'b','LineWidth',2 )
ylim([0 8])
xlim([-20 20])
ylabel('range (R)')
xlabel('Temperature (K)')
title(model)
grid on
legend('actual - moist corrected T_2', 'Location', 'SouthWest');


cd('/Users/spuler/Desktop') % point to the directory to print
FigH = figure(201);
set(gca,'Fontsize',18,'Fontweight','b'); % 
set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 800 800]);
name=strcat(model, 'wv_coorected_T');
print(FigH, name, '-dpng', '-r300')  

cd('/Users/spuler/Desktop') % point to the directory to print
FigH = figure(202);
set(gca,'Fontsize',18,'Fontweight','b'); % 
set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 800 800]);
name=strcat(model, 'diff_wv_coorected_T');
print(FigH, name, '-dpng', '-r300')  

end