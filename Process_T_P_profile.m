function [T, P] = Process_T_P_profile(range, Surf_T, Surf_P, flag)
% Process the HSRL channels 
% Calculate the backscatter ratio (BSR) and backscatter coefficient using receiver scan data
 
 tic
 dd = pwd;


%% Calculate T & P profiles from surface station

 const.k_B = 1.380649e-23; % (J/K, or Pa m^3/K)
 const.K_B = const.k_B * 9.869e-6; % (atm m^3/K)
 const.N_A = 6.022E23; % (/mol) Avagadros number
 const.M = 28.97E-3; % (kg/mol) air molecular mass per mol 
 const.m = const.M./const.N_A; % (kg) mass of a air molecule
 const.g = 9.80665; % (m/s^2) gravitational constant
 const.c = 299792458; % (m/s) (exact)  
 const.R = const.k_B*const.N_A; % (J/K mol) universal gas constant

 % Calculate temperature and pressure profile based on surface measurement and 
 % assuming a standard lapse rate (-6.5 deg/km) for the entire troposphere
 lapse = 0.0065; % (K/m) standard atmosphere lapse rate, dry adiabatic is 9.8 K/km
 if flag.WS == 1  % use surface values if they exist
   T0 = nanmedian(Surf_T)+273.15
   P0 = nanmedian(Surf_P)
   T = (Surf_T+273.15)-lapse.*range;
   P = Surf_P.*((Surf_T+273.15)./T).^-((const.M*const.g)/(const.R*lapse));   % barometric formula
 else
   T0 = 290; % surface temperature
   P0  = 1; % surface pressure in atm 
   T = T0-lapse.*range; % temperature as function of range with lapse rate
   P = P0.*(T0./T).^-((const.M*const.g)/(const.R*lapse));   % barometric formula 
 end
 
 toc

cd(dd)

  
end






