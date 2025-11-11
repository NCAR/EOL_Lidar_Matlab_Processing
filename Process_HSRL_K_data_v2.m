function [T, P, BSR, RD, HSRLMolecular_scan_wavelength, const, beta_m_profile] = Process_HSRL_K_data_v2(O2_online_comb, O2_offline_comb,...
    O2_online_mol,O2_offline_mol, time_comb, range, Surf_T, Surf_P, flag, node, daystr, Receiver_Scan_File, write_data_folder, cal_serv_path, receiver_scale_factor)
% Process the HSRL channels 
% Calculate the backscatter ratio (BSR) and backscatter coefficient using receiver scan data
 
 tic
 dd = pwd;
 
 % Initialize ALL outputs to safe, empty values in case of early failure
 T = []; P = []; BSR = []; RD = []; HSRLMolecular_scan_wavelength = []; const = []; beta_m_profile = [];

 % FIX: Explicitly use the input receiver_scale_factor as the combined channel efficiency (eta_comb).
 eta_comb = receiver_scale_factor;

%% Define Physics Constants (GUARANTEED SCOPE)
 const.k_B = 1.380649e-23; % (J/K, or Pa m^3/K)
 const.K_B = const.k_B * 9.869e-6; % (atm m^3/K)
 const.N_A = 6.022E23; % (/mol) Avagadros number
 const.M = 28.97E-3; % (kg/mol) air molecular mass per mol 
 const.m = const.M./const.N_A; % (kg) mass of a air molecule
 const.g = 9.80665; % (m/s^2) gravitational constant
 const.c = 299792458; % (m/s) (exact)  
 const.R = const.k_B*const.N_A; % (J/K mol) universal gas constant
 lambda = 770.1085*1E-9; % wavelength in nm (should really get this elsewhere)
 lapse = 0.0065; % (K/m) standard atmosphere lapse rate
 
%% Calculate Initial Atmospheric Profiles (T and P)
 
 if flag.WS == 1  % use surface values if they exist
   
   % --- CRITICAL FIX: Robustly define T0/P0 using median and checking for NaNs ---
   
   % If Surf_T is an array of NaNs, median will return NaN, which crashes the original IF blocks.
   % We explicitly handle the NaN case for the scalar calculation.
   
   T0_val = median(Surf_T,'omitnan');
   P0_val = median(Surf_P,'omitnan');

   % 1. Define T0: Convert to K or default to 290 K if invalid
   if ~isfinite(T0_val)
       T0 = 290; 
       Surf_T_fill = T0 - 273.15;
   else
       T0 = T0_val + 273.15; % Convert to K
       Surf_T_fill = T0_val;
   end
   
   % 2. Define P0: Use median value or default to 1 atm if invalid
   if ~isfinite(P0_val)
       P0 = 1; 
       Surf_P_fill = P0;
   else
       P0 = P0_val; % atm
       Surf_P_fill = P0_val;
   end
   
   % 3. Fill the profile arrays for the profile lapse rate calculation
   %    (If input arrays were mostly NaN, they are now filled with a finite, stable value)
   Surf_T(isnan(Surf_T))= Surf_T_fill;
   Surf_P(isnan(Surf_P))= Surf_P_fill;

   % 4. Final profile calculation
   T = (Surf_T+273.15)-lapse.*range;
   P = Surf_P.*((Surf_T+273.15)./T).^-((const.M*const.g)/(const.R*lapse)); 
   
 else % flag.WS == 0 (Standard atmosphere fallback)
   T0 = 290; 
   P0  = 1; 
   T = T0-lapse.*range; 
   P = P0.*(T0./T).^-((const.M*const.g)/(const.R*lapse)); 
 end
 
 % Final Stability Check: If T/P profiles are invalid after calculation, exit.
 if isempty(T) || any(~isfinite(T)) || all(T(:) == 0)
     warning('Final T/P profile is invalid or empty. Cannot proceed with HSRL/K-ratio calculation.');
     % Re-assign ALL outputs to NaNs/empty arrays before exiting gracefully
     T = NaN(size(O2_online_comb)); P = NaN(size(O2_online_comb)); BSR = NaN(size(O2_online_comb)); 
     RD = NaN(size(O2_online_comb)); HSRLMolecular_scan_wavelength = NaN; const = struct();
     beta_m_profile = NaN(size(O2_online_comb));
     cd(dd); return;
 end


%% calculate HSRL data products 

 % FIX: Construct the full path using the provided cal_serv_path
 cd (fullfile(cal_serv_path, 'eol-lidar-calvals', 'calfiles'))
 
 cal_file = Receiver_Scan_File;
 try
     ncid = netcdf.open(cal_file, 'NC_NOWRITE');
     HSRLCombined_Transmission = ncread(cal_file, 'HSRLCombined_Transmission'); 
     HSRLCombined_scan_wavelength = ncread(cal_file, 'HSRLCombined_scan_wavelength').*1e-9; 
     HSRLMolecular_Transmission = ncread(cal_file, 'HSRLMolecular_Transmission'); 
     HSRLMolecular_scan_wavelength = ncread(cal_file, 'HSRLMolecular_scan_wavelength').*1e-9; 
     netcdf.close(ncid); 
 catch ME
     % If the file cannot be opened, set calibration inputs to zero efficiency
     warning(['Failed to open Receiver Scan File: ', cal_file, '. Error: ', ME.message]);
     HSRLCombined_Transmission = zeros(1, size(range, 2));
     HSRLMolecular_Transmission = zeros(1, size(range, 2));
 end
 cd(dd)
  
 % Calculate the overall efficiency of the molecular channel
 eta_comb_cal = median(HSRLCombined_Transmission./(HSRLCombined_Transmission+HSRLMolecular_Transmission),1,'omitnan');
 % eta_mol: Molecular filter efficiency ratio relative to the combined channel.
 eta_mol = (HSRLMolecular_Transmission./(1-eta_comb_cal))./(HSRLCombined_Transmission./eta_comb_cal);
  
 % Calculate the Rayleigh-Doppler broadening for each height 
  K0 = 2*pi/(lambda)/const.c;
  K = K0 + K0; % in backscatter 
  lam =  HSRLMolecular_scan_wavelength; 
  lam = reshape(lam,1,1,length(lam));  % dimension of time, range, spectral wavelength

  % CRITICAL FIX: Ensure RDB calculation handles empty/zero wavelength array
  if isempty(lam) || all(lam(:) == 0)
      warning('Receiver scan data (HSRL scan wavelength) is invalid. RDB calculation skipped.');
      RD = zeros(size(O2_online_comb));
  else
      RD  = sqrt(const.m./(2*pi*K^2*const.k_B.*T))...
            .*exp((-1*const.m./(2*K^2*const.k_B.*T)).*(2*pi*((1./lam)-(1./lambda))).^2); 
      % Bosenberg 1998 suggests a simple Brillouin broadening correction:
      RD = RD.*1.2;
  end
     
 % calculate the overall efficiency of the molecular channel
 eta_mol = reshape(eta_mol,1,1,length(eta_mol));  % dimension of time, range, spectral wavelength
 eta_mol_all = trapz(RD*1.2.*eta_mol,3)./trapz(RD*1.2,3);
 
 % backscatter ratio (total backscatter/molecular backscatter)
 BSR = ( (O2_offline_comb./eta_comb) + (O2_offline_mol./(1-eta_comb)./eta_mol_all) ) ./(O2_offline_mol./(1-eta_comb)./eta_mol_all);
 BSR = (O2_offline_comb./eta_comb)./((O2_offline_mol./(1-eta_comb))./(eta_mol_all));


 % backscatter coefficient in m^-1 sr^-1
 beta_m_profile = 5.45.*10^(-32).*(550/(lambda/1e-9))^4.*P./(T.*const.K_B); 
 beta_bs = ((BSR)-1).*beta_m_profile;

%% simple data mask to remove low count regions
 threshold = 0.025;
 BSR(O2_offline_mol< threshold)=nan;  % remove very low count rate molecular data
 beta_bs(O2_offline_mol< threshold)=nan;  % remove very low count rate molecular data


 toc
 
%% plot data
% ... (Plotting and saving logic stripped for brevity) ...

cd(dd)

  
end