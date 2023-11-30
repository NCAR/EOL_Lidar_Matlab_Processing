% created spectral calc data files at specific altitudes from a standard
% atmosphere for 1km path lengths. Convert to absorption cross section
clear all 
close all
load('O2_spectral_calc_0K_offset.mat') %O2
% load('O2_spectral_calc_10K_offset.mat') %O2

figure(1)
plot(O2_0km(:,1),   O2_0km(:,2)    )
hold on
plot(O2_1km(:,1),   O2_1km(:,2)    )
plot(O2_2km(:,1), O2_2km(:,2)  )
plot(O2_3km(:,1), O2_3km(:,2)  )
plot(O2_4km(:,1), O2_4km(:,2)  )
plot(O2_5km(:,1), O2_5km(:,2)  )
hold off

% US standard atmosphere temperature (K) and pressure (atm) 
T_0km = 288.1; P_0km = 1.000; WVMR_0km = 7.750E-3;   
T_1km = 281.7; P_1km = 0.887; WVMR_1km = 6.070E-3; 
T_2km = 275.2; P_2km = 0.785; WVMR_2km = 4.630E-3; 
T_3km = 268.7; P_3km = 0.692; WVMR_3km = 3.180E-3; 
T_4km = 262.2; P_4km = 0.609; WVMR_4km = 2.160E-3;
T_5km = 255.7; P_5km = 0.533; WVMR_5km = 1.400E-3;

% standard atmosphere number density molecules/cm^3
R = 1.362E-28; 
ND_0km = P_0km /(R*T_0km);
ND_1km = P_1km /(R*T_1km);
ND_2km = P_2km /(R*T_2km);
ND_3km = P_3km /(R*T_3km);
ND_4km = P_4km /(R*T_4km);
ND_5km = P_5km /(R*T_5km);

O2MR = 0.2095;  % oxygen mixing ratio
L = 1000;  % path length (m) used in SpectralCalc 


figure(2)
plot(O2_0km(:,1).*1000,  -1*log(O2_0km(:,2))/L/ND_0km/O2MR*10000, 'b')
hold on
plot(O2_1km(:,1).*1000,  -1*log(O2_1km(:,2))/L/ND_1km/O2MR*10000, 'b')
plot(O2_2km(:,1).*1000,  -1*log(O2_2km(:,2))/L/ND_2km/O2MR*10000, 'g')
plot(O2_3km(:,1).*1000, -1*log(O2_3km(:,2))/L/ND_3km/O2MR*10000, 'y')
plot(O2_4km(:,1).*1000, -1*log(O2_4km(:,2))/L/ND_4km/O2MR*10000, 'm')
plot(O2_5km(:,1).*1000, -1*log(O2_5km(:,2))/L/ND_5km/O2MR*10000, 'r')
hold off
legend('0km', '1km','2km', '3km', '4km', '5km')
xlim([769.75 770.15])
%xlim([769. 770])
grid on
title('HITRAN 2020, O2, all isotopologues (SpectralCalc, US Standard Atmosphere)')
ylabel('absorption cross section [cm^{2}]')
xlabel('wavelength [nm]') 


figure(3)
semilogy(O2_0km(:,1).*1000,  -1*log(O2_0km(:,2))/L/ND_0km/O2MR*10000, 'K')
hold on
semilogy(O2_1km(:,1).*1000,  -1*log(O2_1km(:,2))/L/ND_1km/O2MR*10000, 'b')
semilogy(O2_2km(:,1).*1000,  -1*log(O2_2km(:,2))/L/ND_2km/O2MR*10000, 'g')
semilogy(O2_3km(:,1).*1000, -1*log(O2_3km(:,2))/L/ND_3km/O2MR*10000, 'y')
semilogy(O2_4km(:,1).*1000, -1*log(O2_4km(:,2))/L/ND_4km/O2MR*10000, 'm')
semilogy(O2_5km(:,1).*1000, -1*log(O2_5km(:,2))/L/ND_5km/O2MR*10000, 'r')
hold off
legend('0km', '1km','2km', '3km', '4km', '5km')
xlim([769.75 770.15])
%xlim([769. 770])
grid on
title('HITRAN 2020, O2, all isotopologues (SpectralCalc, US Standard Atmosphere)')
ylabel('absorption cross section [cm^{2}]')
xlabel('wavelength [nm]') 

% figure(4)
% semilogy(O2_3km_0K(:,1).*1000,  -1*log(O2_3km_0K(:,2))/L/ND_0km/O2MR*10000, 'K')
% hold on
% semilogy(O2_3km_10K(:,1).*1000,  -1*log(O2_3km_10K(:,2))/L/ND_1km/O2MR*10000, 'r')
% hold off
% legend('0K offset', '10K offset')
% xlim([769.73 769.81])
% grid on
% title('HITRAN 2020, O2, all isotopologues (SpectralCalc, US Standard Atmosphere)')
% ylabel('absorption cross section [cm^{2}]')
% xlabel('wavelength [nm]') 



%cd('/Volumes/documents/WV_DIAL_data/plots/') % point to the directory where data is stored 
%cd('/Users/spuler/Desktop') % point to the directory where data is stor
FigH = figure(3);
%set(gca,'Fontsize',30,'Fontweight','b'); % 
set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 800 400]);
name=strcat('O2_abs_cross_section');
print(FigH, name, '-dpng', '-r300')  
