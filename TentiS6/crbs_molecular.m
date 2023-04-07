% Written by John Smith
% October 21st, 2010
% University of Colorado at Boulder, CIRES
% John.A.Smith@Colorado.EDU
% MATLAB version 7.10.0.59 (R2010a) 64-bit
% Adapted from "Coherent Rayleigh-Brillouin Scattering"
% by Xingguo Pan
% Sets and computes all relevant parameters
% for s6 and s7 models given in Xingguo Pan's
% dissertation entitled "Coherent Rayleigh-
% Brillouin Scattering", 2003, Princeton Univ.
% 
% downloaded from www.mathworks.com/matlabcentral/fileexchange/
% adapted/modified for backscatter in air by Scott Spuler Nov 2021 

% set the temperature of the gas
tem=262.7; % US standard model at 4km
p_atm=0.609; 
tem=281.7; % US standard model at 1km
p_atm=0.887;
% tem=287; % 
% p_atm=0.82; % 

% set laser parameters
lambda=1570e-9;
angle=(180)*(pi/180);  % Set at 180 for backscatter

% create the domain of xi
N=1000;
xi_lef=-5.0d0;
xi_rgt=5.0d0;
xi=linspace(xi_lef,xi_rgt,N);
% *** 
% 'xi' is dimensionless and should be scaled by k*v0/(2*pi) [Hz],
% where 'v0' is the most probable gas velocity, or sqrt(2*kb*T/m),
% and 'k' is 4*pi/lambda*sin(scatter_angle/2)
% ***

% set fundamental constants
kb=1.3806503e-23;

% The Matlab File Exchange version had an error.  The k term needed to use angle/2 
k=sin(angle/2)*4*pi/lambda;

% Modified to update mass, viscosity, bulk viscosity, and thermal conductivity for air 
% set N2 gas quantities
%m_m=(1.66053886e-27)*28.013;
    m_m=(1.66053886e-27)*28.97; % adjusted for air
%viscosity=17.63e-6; 
   T_0 = 273; % K
   eta_0 = 1.716e-5; %Pa s
   S_eta = 111; %K
   viscosity= eta_0*(tem/T_0)^(3/2)*(T_0+S_eta)/(tem+S_eta); % Sutherland law
%bulk_vis=viscosity*0.73; 
     bulk_vis=viscosity*0.71; % used by Luke Colberg/MSU
     bulk_vis=0.86e-5+(1.29e-7*(tem-250)); % Wang et al. Molecular Physics 2021
%thermal_cond=25.2e-3;
     kappa_0 = 0.0241; % 1W/m K
     S_th = 194; %K
     thermal_cond=kappa_0*(tem/T_0)^(3/2)*(T_0+S_th)/(tem+S_th); % Sutherland law
c_int=1.0;
% compute most probable gas velocity
v0=sqrt(2*kb*tem/m_m);
% convert pressures and densities
p_pa=p_atm*1.01325e5;
p_torr=p_pa*0.00750061683;
n0=p_pa/(tem*kb);
% compute and set RBS model input parameters
c_tr=3/2;
y=n0*kb*tem/(k*v0*viscosity);
gamma_int=c_int/(c_tr+c_int);
rlx_int=1.5*bulk_vis/(viscosity*gamma_int);
eukenf=m_m*thermal_cond/(viscosity*kb*(c_tr+c_int));
% run the code
%[cohsig7,sptsig7]=crbs7(y,rlx_int,eukenf,c_int,c_tr,xi);
[cohsig6,sptsig6]=crbs6(y,rlx_int,eukenf,c_int,c_tr,xi);
% OUTPUTS-
% **cohsig7: coherent RBS spectrum using s7 model**
% **cohsig6: coherent RBS spectrum using s6 model**
% **sptsig7: spontaneous RBS spectrum using s7 model**
% **sptsig6: spontaneous RBS spectrum using s6 model**


%% 
% Plot results against Rayleigh Doppler
% ignore for the PCA 

figure(2)
GHz = xi*k*v0/(2*pi);  %scaled by k*v0/(2*pi) [Hz]
Norm_Spon_Tenti = sptsig6./trapz(sptsig6);
%Norm_Coh_Tenti = cohsig6./trapz(cohsig6);
plot(GHz, Norm_Spon_Tenti, 'r')
% hold on
% plot(GHz, Norm_Coh_Tenti, 'g')
% hold off

% calculate the RD only spectrum
const.k_B = 1.380649e-23; % (J/K, or Pa m^3/K)
const.c = 299792458; % (m/s) (exact)  
const.N_A = 6.022E23; % (/mol) Avagadros number
const.M = 28.97E-3; % (kg/mol) air molecular mass per mol 
const.m = const.M./const.N_A; % (kg) mass of a air molecule

T = tem;
% from Fiocco and DeWolf 1968
K0 = 2*pi/(lambda)/const.c;
K = K0 + K0; % in backscatter 
lam = (lambda-(0.01*1e-9)):(0.00002*1e-9):(lambda+(0.01*1e-9));  
RD = sqrt(const.m./(2*pi*K^2*const.k_B.*T))...
         .*exp((-1*const.m./(2*K^2*const.k_B.*T)).*(2*pi*((1./lam)-(1./lambda))).^2);

hold on
figure(2)
Norm_RD = RD./trapz(RD,2);
% Norm_RD = RD./max(RD);
plot(const.c.*(lam-lambda)/lambda^2, Norm_RD, 'k')
legend('Spontaneous RDB', 'RD-only')
%legend('Spontaneous RDB', 'Coherent RDB', 'RD-only')
% legend('Tenti-S6 RDB (h=1km)', 'Tenti-S6 RDB (h=4km)')
xlabel('Frequency')
hold off
grid on
xlim([-3e9 3e9]);
% xlim([-1.5e9 1.5e9]);

% save the figure as a png   
% cd('/Users/spuler/Desktop') % point to the directory 
% FigH = figure(2);
% % set(gca,'Fontsize',30,'Fontweight','b'); % 
% set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 800 600]);
% name=strcat('Tenti_S6_comparison_T', num2str(tem), 'P', num2str(p_atm));
% print(FigH, name, '-dpng', '-r300')  



