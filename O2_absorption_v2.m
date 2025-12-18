function [] = O2_absorption_v2(const, T, P, O2_online_comb, O2_offline_comb, ...
    time_comb, range, BSR, RD, HSRLMolecular_scan_wavelength, N_WV, beta_m_profile, O2_on_wavelength, node, daystr, write_data_folder, flag)

d = pwd;

% --- ARCHITECTURAL STABILITY CHECK ---
% If T profile is empty, non-finite (NaN), or a single zero (symptom of upstream HSRL failure), exit gracefully.
% This prevents the "Invalid use of operator" crash when HSRL/K-ratio failed to produce outputs.
if isempty(T) || ~any(isfinite(T(:))) || all(T(:) == 0)
    warning('O2 Absorption v2 skipped: Invalid or empty Temperature profile input from HSRL step.');
    cd(d);
    return;
end
% -------------------------------------

% calculate the density of air
vmr_O2_d =  0.2095; % volume mixing ratio of oxygen, dry

% Read in the number density of water vapor and measured via WV DIAL 
% convert it to volume mixing ratio  
try
vmr_WV = N_WV./P.*(const.K_B.*T);
catch
 N_WV = 0; %assuming zero 
 vmr_WV = N_WV./P.*(const.K_B.*T);
end

% volume mixing ratio of oxygen (corrected for WV)
vmr_O2 = vmr_O2_d*(1-vmr_WV);
% number denisty of oxygen (corrected for WV)
N_O2 =  P./(const.k_B.*T).*vmr_O2; % (m^-3) number density  

gate = range(2);
% override the native gate spacing (needs to be changed in the future)
gates2use =  225/gate;
%if strcmp(node,'MPD04')==1
%    gates2use = 6*5;
%end
gate_int = range(2)*gates2use;


% zeroth order absorption coefficient 
% which is the measured extinction of oxgyen

 Inside = (circshift(O2_online_comb,[0, gates2use/2]).* circshift(O2_offline_comb,[0,-gates2use/2]))./...
          (circshift(O2_online_comb,[0,-gates2use/2]).* circshift(O2_offline_comb,[0, gates2use/2]));
 del_cross = single(1./(2*gate_int));
 alpha_O2_0th = del_cross.*log(Inside);  % absorption coef in m^-1

% another way to calculate the DIAL equation with log 
  pad_on = [O2_online_comb O2_online_comb(:,end)]';   
  pad_off = [O2_offline_comb O2_offline_comb(:,end)]';       
  der_on = (diff(pad_on)./(gate))';
  der_off = (diff(pad_off)./(gate))';
  Inside_v2 = ((1./O2_online_comb).*der_on - (1./O2_offline_comb).*der_off);
  alpha_O2_0th2 = double(-1./(2).*(Inside_v2));  % absorption coef in m^-1
 % alpha_O2_0th_sg = sgolayfilt(alpha_O2_0th2,3,5);
  
 
 % for Rayleigh scattering theory, the lidar ratio is equal to 8*pi/3
 % attenuation due to molecules is 
 % alpha_m_off = 8*pi/3*beta_m_profile;
 % alpha_0 =   alpha_m_off+alpha_O2;  
 % alpha_O2_off = Spectra.O2Offline.AbsorptionObserved
 
 % smooth the result over the range 
 alpha_O2_avg = nanmoving_average(alpha_O2_0th,gates2use/2,2,0);
 alpha_O2_avg2 = nanmoving_average(alpha_O2_0th2,gates2use/2,2,0);
% alpha_O2_avg2 = nanmoving_average(alpha_O2_0th_sg,gates2use/2,2,0);
  
 % ... (Plotting, file saving logic skipped for brevity) ...
  
  if flag.save_data == 1
    cd(write_data_folder)
    alpha_O2 = alpha_O2_avg2; % Using the smoothed average
    name=strcat(node, "_", daystr, "_Backscatter Coefficient");
    save(name, 'alpha_O2', '-append')  
  end
   
   
   cd(d)
   
end