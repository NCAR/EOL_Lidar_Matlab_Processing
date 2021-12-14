function [] = O2_absorption(const, T, P, O2_online_comb, O2_offline_comb, ...
    time_comb, range, BSR, RD, HSRLMolecular_scan_wavelength, N_WV, beta_m_profile, O2_on_wavelength, node, daystr)

 d = pwd;

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
% override the native gate spacing
gates2use = 10;
gate_int = range(2)*gates2use;


% Perturbative retrieval 
% zeroth order term Eq. 11 from Repasky et al. 2019 Optics Express
 Inside = (O2_online_comb.*(circshift(O2_offline_comb, [0, -1*gates2use])))./...
     ((circshift(O2_online_comb, [0, -1*gates2use])).*O2_offline_comb);
 del_cross = single(1./(2.*gate_int*100));
 % zeroth order absorption coefficient 
 % which is the number density of oxygen times absorption cross section

 % for Rayleigh scattering theory, the lidar ratio is equal to 8*pi/3
 alpha_m_off = 8*pi/3*beta_m_profile;
 
 %test = (del_cross.*log(Inside));
 
 alpha_0 =   alpha_m_off-(del_cross.*log(Inside));  
 O2_abs = del_cross.*log(Inside);  
 
 
 % smooth the result over the range 
  O2_abs_avg = O2_abs;
  O2_abs_avg = nanmoving_average(O2_abs_avg,gates2use/2,2,0);
 
 
  figure(10)
  x = (time_comb)';
  y = (range./1e3);
  %Z = real(double((((alpha_0))')));
  Z = real(double((((O2_abs))')));
  Z = real(double((((O2_abs_avg))')));
  font_size = 14;
  xData =  linspace(fix(min(time_comb)),  ceil(max(time_comb)), 25);
  set(gcf,'renderer','zbuffer');
  h = pcolor(x,y,Z);
  set(h, 'EdgeColor', 'none');
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  set(gca, 'XTick',  xData)
  colorbar('EastOutside');
  axis([fix(min(time_comb)) fix(min(time_comb))+1 0 4])
  caxis([1e-8 20e-7]);
  %caxis([1e-6 3e-6]);
  datetick('x','HH','keeplimits', 'keepticks');
  xlabel('Time (UTC)','fontweight','b','fontsize',font_size);
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size);
  set(gca,'Fontsize',font_size,'Fontweight','b');
%   set(gca,'Zscale', 'log')
%   set(gca,'Colorscale', 'log')
%   set(gca,'Zscale', 'linear')
 % colormap(jet)
  
figure(400)
plot(O2_on_wavelength(600:650),'ro')
ylim([769.79 769.81])

figure(401)
plot(O2_on_wavelength(600:650), nanmedian(Z(15:50, 600:650),1), 'bo')
% hold on
% plot(O2_on_wavelength(600:650), Z(40, 600:650), 'ro')
% plot(O2_on_wavelength(600:650), Z(50, 600:650), 'go')
% plot(O2_on_wavelength(600:650), Z(60, 600:650), 'mo')
% hold off



 
   % plot a single profile
   figure(20)
   [minValue, closestIndex] = min(abs(time_comb - datenum('13-Dec-2021 17:00:00')))
   alpha_0_profile = alpha_0(closestIndex, :); 
   O2_abs_profile = O2_abs(closestIndex, :); 
 %  plot(alpha_0_profile, (range./1e3));
   plot(O2_abs_profile, (range./1e3));
   hold on
   [minValue, closestIndex] = min(abs(time_comb - datenum('13-Dec-2021 17:55:00')))
   alpha_0_profile = alpha_0(closestIndex, :);
   O2_abs_profile = O2_abs(closestIndex, :); 
%    plot(alpha_0_profile, (range./1e3)); 
   plot(O2_abs_profile, (range./1e3));
   ylim([0 4])
   %xlim([1e-7 1e-5])
   %xlim([1e-6 5e-6])
   hold off
   
   
   
   
   
   serv_path1 = '/Volumes/Macintosh HD/Users/spuler/Desktop/mpd/';  
   cd(strcat(serv_path1,'/Plots'))
   
   FigH = figure(10);
   set(gca,'Fontsize',16,'Fontweight','b'); 
   set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920/2 250]);
   name=strcat(node, "_", daystr, "_O2_absorption");
   print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
   
   cd(d)
   
  
   
end
  
%   % BSR is the backscatter ratio in time vs range
%   % RD is the Rayleigh Doppler lineshape (wavelength vs range)
%   % create a normalized backscatter lineshape (assumes laser is a delta
%   % function)
%   % Eq 12 from Repasky et al. 2019 Optics Express
%   g_x = BSR + BSR.*RD;
%   
%   figure(2)
%  % [value, index] = max(RD(:,1))
% %  plot(HSRLMolecular_scan_wavelength, g_x(1,1,:)./(RD(index,1)) )
%   g_x_profile =  squeeze(g_x(closestIndex,:,:)); 
%   plot(HSRLMolecular_scan_wavelength, g_x_profile(1500/gate,:))     
%   hold on
%   plot(HSRLMolecular_scan_wavelength, g_x_profile(3000/gate,:))
%   hold off

%   % Eq 13 from Repasky et al. 2019 Optics Express
%   % Normalized lineshap of the transmission function of the optical filter in the MPD receiver
%   E_nu = 1;  
%   % 
%   % Molecular absorption lineshape with maxium value of 1 at line center
%   f = 
%   % this seems to be some type
%   Tm_0th = exp(int(alpha_0)*f
% 
%   Delta_W_1st = g_x*E*Tr_filter
%   where Tr_filter = 1 
%   
%   
%   % Eq 13 from Repasky et al. 2019 Optics Express
%   % the first order correction to the absorption coefficient
%   Delta_alpha_1 = 0.5.*alpha_0; 
%   
%   
%   trapz(RD(:,1)./trapz(RD(:,1)));
% 
%   
% % end

