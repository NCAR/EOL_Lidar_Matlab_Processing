function [] = O2_absorption(const, T, P, O2_online_comb, O2_offline_comb, ...
    time_comb, range, BSR, RD, HSRLMolecular_scan_wavelength, N_WV, beta_m_profile, O2_on_wavelength, node, daystr, write_data_folder, flag)

 d = pwd;

%  O2_absorption(T, P, O2_online_comb, O2_offline_comb, ...
%      time_comb, range, BSR, RD, HSRLMolecular_scan_wavelength)

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
gates2use = 4;
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
  alpha_O2_0th_sg = sgolayfilt(alpha_O2_0th2,3,5);
  
 
 % for Rayleigh scattering theory, the lidar ratio is equal to 8*pi/3
 % attenuation due to molecules is 
 % alpha_m_off = 8*pi/3*beta_m_profile;
 % alpha_0 =   alpha_m_off+alpha_O2;  
 % alpha_O2_off = Spectra.O2Offline.AbsorptionObserved
 
 % smooth the result over the range 
 alpha_O2_avg = nanmoving_average(alpha_O2_0th,gates2use/2,2,0);
% alpha_O2_avg2 = nanmoving_average(alpha_O2_0th2,gates2use/2,2,0);
 alpha_O2_avg2 = nanmoving_average(alpha_O2_0th_sg,gates2use/2,2,0);
  
  figure(11)
  x = (time_comb)';
  y = (range./1e3);
  %Z = real(double((((alpha_0))')));
  %Z = real(double((((O2_abs))')));
  Z = real(double((((alpha_O2_avg))')));
  font_size = 14;
  xData =  linspace(fix(min(time_comb)),  ceil(max(time_comb)), 25);
  set(gcf,'renderer','zbuffer');
  h = pcolor(x,y,Z);
  set(h, 'EdgeColor', 'none');
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  set(gca, 'XTick',  xData)
  colorbar('EastOutside');
  axis([fix(min(time_comb)) fix(min(time_comb))+1 0 8])
  caxis([1e-6 3e-4]);
  %caxis([1e-6 3e-6]);
  hh = title({[node, ' ', daystr, ' O2 absorption coeff [m^{-1}]']},'fontweight','b','fontsize',font_size);  
  datetick('x','HH','keeplimits', 'keepticks');
  xlabel('Time (UTC)','fontweight','b','fontsize',font_size);
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size);
  set(gca,'Fontsize',font_size,'Fontweight','b');
%   set(gca,'Zscale', 'log')
%   set(gca,'Colorscale', 'log')
%   set(gca,'Zscale', 'linear')
%   colormap(jet)
 
 
[~, startInd] = min(abs(time_comb - datenum('15-Dec-2021 03:15:00')));
[~, endInd] = min(abs(time_comb - datenum('15-Dec-2021 05:15:00')));

% 
% figure(400)
% plot(O2_on_wavelength(startInd:endInd),'ro')
% ylim([769.79 769.80])
% 
% figure(401)
% plot(O2_on_wavelength(startInd:endInd), nanmedian(Z(10:30, startInd:endInd),1), 'bo')
% hold on
% plot(O2_on_wavelength(startInd:endInd), nanmedian(Z(30:50, startInd:endInd),1), 'ro')
% plot(O2_on_wavelength(startInd:endInd), nanmedian(Z(50:70, startInd:endInd),1), 'go')
% xlim([769.79 769.80])
% ylim([0 2e-6])
% grid on
% hold off
%  

   
   
%    serv_path1 = '/Volumes/Macintosh HD/Users/spuler/Desktop/mpd/';  
%    cd(strcat(serv_path1,'/Plots'))
%    
%    FigH = figure(10);
%    set(gca,'Fontsize',16,'Fontweight','b'); 
%    set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920/2 250]);
%    name=strcat(node, "_", daystr, "_O2_absorption");
%    print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
   
  
  if flag.save_data == 1
    cd(write_data_folder)
    alpha_O2 = Z';
    name=strcat(node, "_", daystr, "_Backscatter Coefficient");
%    try
      save(name, 'alpha_O2', '-append')  
%    catch
%      save(name, 'alpha_O2')  
%    end
  end
   
   
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

