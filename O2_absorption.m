function [] = O2_absorption(const, T, P, O2_online_comb, O2_offline_comb, ...
    time_comb, range, BSR, RD, HSRLMolecular_scan_wavelength, N_WV, beta_m_profile)


%  O2_absorption(T, P, O2_online_comb, O2_offline_comb, ...
%      time_comb, range, BSR, RD, HSRLMolecular_scan_wavelength)

%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%outputArg1 = inputArg1;
%outputArg2 = inputArg2;

% calculate the density of air
vmr_O2_d =  0.2095; % volume mixing ratio of oxygen, dry

% Read in the number density of water vapor and measured via WV DIAL 
% convert it to volume mixing ratio  
% note this is not the same a mixing ratio which is the mass ratio of water
% vapor which is the mass of water vapor 18 gm/mol divided by the mass of
% dry air (28.97 gm/mole) 
% N_WV = 0; %assuming zero 
vmr_WV = N_WV./P.*(const.K_B.*T);
% volume mixing ratio of oxygen (corrected for WV)
vmr_O2 = vmr_O2_d*(1-vmr_WV);
% number denisty of oxygen (corrected for WV)
N_O2 =  P./(const.k_B.*T).*vmr_O2; % (m^-3) number density  

gate = range(2);

 Inside = (O2_online_comb.*(circshift(O2_offline_comb, [0, -1])))./...
     ((circshift(O2_online_comb, [0, -1])).*O2_offline_comb);
 del_cross = single(1./(2.*gate*100));
 % zeroth order absorption coefficient 
 % which is the number density of oxygen times absorption cross section
 % currently ignoring the molecular absoprtion coefficient at the offline
 alpha_m_off = 0;
 % for Rayleigh scattering theory, the lidar ratio is equal to 8*pi/3
 alpha_m_off = 8*pi/3*beta_m_profile;
 
 %test = (del_cross.*log(Inside));
 
 alpha_0 =   alpha_m_off-(del_cross.*log(Inside));  

  figure(10)
  x = (time_comb)';
  y = (range./1e3);
  Z = real(double((((alpha_0))')));
  font_size = 14;
  xData =  linspace(fix(min(time_comb)),  ceil(max(time_comb)), 25);
  set(gcf,'renderer','zbuffer');
  h = pcolor(x,y,Z);
  set(h, 'EdgeColor', 'none');
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  set(gca, 'XTick',  xData)
  colorbar('EastOutside');
  axis([fix(min(time_comb)) fix(min(time_comb))+1 0 6])
  caxis([1e-7 1e-5]);
  datetick('x','HH','keeplimits', 'keepticks');
  xlabel('Time (UTC)','fontweight','b','fontsize',font_size);
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size);
  set(gca,'Fontsize',font_size,'Fontweight','b');
  set(gca,'Zscale', 'log')
  set(gca,'Colorscale', 'log')
  set(gca,'Zscale', 'linear')
  colormap(jet)
  
   % plot a single profile
   figure(20)
   [minValue, closestIndex] = min(abs(time_comb - datenum('21-Jul-2021 10:00:00')))
   alpha_0_profile = alpha_0(closestIndex, :); 
   %semilogx(alpha_0_profile, (range./1e3));
   plot(alpha_0_profile, (range./1e3));
   ylim([0 6])
   %xlim([1e-7 1e-5])
   xlim([-1e-6 1e-5])
  
  % the first order correction to the absorption coefficient
  Delta_alpha_1 = 0.5.*alpha_0; 
  
  % BSR is the backscatter ratio in time vs range
  % RD is the Rayleigh Doppler lineshape (wavelength vs range)
  % create a normalized backscatter lineshape (assumes laser is a delta
  % function)
  g_x = BSR + BSR.*RD;
  
  figure(2)
 % [value, index] = max(RD(:,1))
%  plot(HSRLMolecular_scan_wavelength, g_x(1,1,:)./(RD(index,1)) )
  g_x_profile =  squeeze(g_x(closestIndex,:,:)); 
  plot(HSRLMolecular_scan_wavelength, g_x_profile(1500/gate,:))     
  hold on
  plot(HSRLMolecular_scan_wavelength, g_x_profile(3000/gate,:))
  hold off
  
  
  trapz(RD(:,1)./trapz(RD(:,1)));

  
end

