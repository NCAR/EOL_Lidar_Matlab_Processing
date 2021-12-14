function [T, P, BSR, RD, HSRLMolecular_scan_wavelength, const, beta_m_profile] = Process_HSRL_K_data(O2_online_comb, O2_offline_comb,...
    O2_online_mol,O2_offline_mol, time_comb, range, Surf_T, Surf_P, flag, node, daystr, Receiver_Scan_File, write_data_folder)
%Process_HSRL_K_data calcuated the backscatter coefficient
%   Detailed explanation goes here

 d = pwd;

  figure(102)
  semilogx(O2_online_comb(100,:), range/1000, 'k')
  hold on
  semilogx(O2_offline_comb(100,:), range/1000, 'm')
  semilogx(O2_online_mol(100,:), range/1000, 'g')
  semilogx(O2_offline_mol(100,:), range/1000, 'b')
  ylim([0 6])
  legend('online comb', 'offline comb', 'online mol', 'offline mol')
  % there is more signal sent to the molecular detector
  % to correct for split, and make signal levels match, comb*(mol/comb = 1.67) 
  O2_geometeric_correction = nanmean(O2_online_comb./O2_online_mol,1);
  %semilogx(O2_online_comb(100,:).*O2_geometeric_correction, range/1000, 'k-');
  hold off

  figure(103)
  plot(O2_online_comb(100,:)./O2_online_mol(100,:), range/1000)
  hold on
  plot(O2_geometeric_correction, range/1000)
  xlim([0 1])
  ylim([0 8])
  hold off
  
 % The plot above shows the combined efficiency, rel to molecular 
   O2_geometeric_correction(O2_geometeric_correction==0)=nan;
   O2_geometeric_correction(O2_geometeric_correction==Inf)=nan;
   O2_geometeric_correction(O2_geometeric_correction==-Inf)=nan;
   eta_comb = nanmedian(O2_geometeric_correction) 
   eta_comb = 0.5988  % use for now as the above can blow up on clouds
  
 %plot relative backscatter
  figure(104)
  x = (time_comb)';
  y = (range./1e3);
  Z = double(real(O2_offline_comb'));
  font_size = 14;
  xData =  linspace(fix(min(time_comb)),  ceil(max(time_comb)), 25);
  set(gcf,'renderer','zbuffer');
  h = pcolor(x,y,Z);
  set(h, 'EdgeColor', 'none');
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  set(gca, 'XTick',  xData)
  colorbar('EastOutside');
  axis([fix(min(time_comb)) fix(min(time_comb))+1 0 12])
  caxis([1e1 1e6]);
  datetick('x','HH','keeplimits', 'keepticks');
  xlabel('Time (UTC)','fontweight','b','fontsize',font_size);
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size);
  set(gca,'Fontsize',font_size,'Fontweight','b');
  set(gca,'Zscale', 'log')
  set(gca,'Colorscale', 'log')
  set(gca,'Zscale', 'linear')

  
   %% Calculate Backscatter from Molecules
% h = waitbar(0,'Line Fiting');

% Calculate temperature and pressure profile based on surface measurement and 
% asuming a standard lapse rate (-6.5 deg/km) for the entire troposphere
if flag.WS == 1
    T0 = nanmedian(Surf_T)+273.15
    P0 = nanmedian(Surf_P)
    %T0 = Surf_T+273.15;  % use the array so it is updated in time
    %P0 = Surf_P;
else
  T0 = 290; % surface temperature
  %P0  = 1; % surface pressure in atm 
end

 
const.k_B = 1.380649e-23; % (J/K, or Pa m^3/K)
const.K_B = const.k_B * 9.869e-6; % (atm m^3/K)
const.N_A = 6.022E23; % (/mol) Avagadros number
const.M = 28.97E-3; % (kg/mol) air molecular mass per mol 
const.m = const.M./const.N_A; % (kg) mass of a air molecule
const.g = 9.80665; % (m/s^2) graviational constant
const.c = 299792458; % (m/s) (exact)  
const.R = const.k_B*const.N_A; % (J/K mol) universal gas constant 

lapse = 0.0065; % (K/m) standard atmosphere lapse rate, dry adiabatic is 9.8 K/km
lambda = 770.1085*1E-9; %wavelength in nm

% temp and pressure profiles using the median surface or changing surface values 
%T = T0-lapse.*range; % temperature as function of range with lapse rate
T = (Surf_T+273.15)-lapse.*range;
%P = P0.*(T0./T).^-5.25588;   % hydrostatic equation and ideal gas law
%P = P0.*(T0./T).^-((const.M*const.g)/(const.R*lapse));   % barometric formula 
P = Surf_P.*((Surf_T+273.15)./T).^-((const.M*const.g)/(const.R*lapse));   % barometric formula

% figure
% plot(T, range)
% figure
% plot(P, range)
% hold on
% plot(P_new, range)


beta_m_profile = 5.45.*10^(-32).*(550/(lambda/1e-9))^4.*P./(T.*const.K_B); %backscatter coefficient in m^-1 sr^-1

%beta_m = 5.45.*10^(-32).*(550/(lambda/1e-9))^4.*1./(273.15.*const.K_B) %backscatter coefficient at STP 


figure(105)
 plot(beta_m_profile(1,:), range)
 ylim([0 8])
 
%% calcuate HSRL data products 

 % only a fraction of total molecular gets through the K filter 
 % need to know filter width, molecular scatter width RB which is range dependent
 % just assume a constant for now
 eta_mol_all = 0.167;


% below is a more complete way to check this efficency 
% open the receiver scan data


% read in receiver scans
    cd /Users/spuler/Documents/GitHub/eol-lidar-calvals/calfiles
  %  if strcmp(node,'MPD05') == 1 
  %     cal_file = 'MPD05_consolidated_20211110.nc';
  %  elseif strcmp(node,'MPD01') == 1 
  %     cal_file = 'MPD01_consolidated_20211116.nc';
  %  end
  cal_file = Receiver_Scan_File;
     ncid = netcdf.open(cal_file, 'NC_NOWRITE');
     % ncdisp(cal_file) % use this to display all variables
      HSRLCombined_Transmission = ncread(cal_file, 'HSRLCombined_Transmission'); 
      HSRLCombined_scan_wavelength = ncread(cal_file, 'HSRLCombined_scan_wavelength').*1e-9; 
      HSRLMolecular_Transmission = ncread(cal_file, 'HSRLMolecular_Transmission'); 
      HSRLMolecular_scan_wavelength = ncread(cal_file, 'HSRLMolecular_scan_wavelength').*1e-9; 
    netcdf.close(ncid); 
    cd(d)
    
   % load MPD05_receiver_scan
    
      figure(252)
      plot(HSRLMolecular_scan_wavelength, HSRLMolecular_Transmission)
      hold on
      plot(HSRLMolecular_scan_wavelength, HSRLCombined_Transmission/eta_comb)
      hold off
    
      eta_mol = HSRLMolecular_Transmission./(HSRLCombined_Transmission/eta_comb);
      
      figure(253)
      plot(HSRLMolecular_scan_wavelength, eta_mol)
 
 % Calculate the Rayliegh-Doppler broadening for each height 
  
     %const.k_B = 1.380649e-23; % (J/K, or Pa m^3/K)
     % from Fiocco and DeWolf 1968
     K0 = 2*pi/(lambda)/const.c;
     K = K0 + K0; % in backscatter 
     %lam = (lambda-(0.01*1e-9)):(0.0001*1e-9):(lambda+(0.01*1e-9));  
     lam =  HSRLMolecular_scan_wavelength;
  
     RD = ones(length(time_comb), length(range), length(lam));
     for i = 1:size(Surf_T,1)
        %T = (Surf_T(i)+273.15)-lapse.*range;
        RD(i,:,:) = sqrt(const.m./(2*pi*K^2*const.k_B.*T(i,:)'))...
         .*exp((-1*const.m./(2*K^2*const.k_B.*T(i,:)')).*(2*pi*((1./lam')-(1./lambda'))).^2);
     end  

%          RD = sqrt(const.m./(2*pi*K^2*const.k_B.*T))...
%            .*exp((-1*const.m./(2*K^2*const.k_B.*T)).*(2*pi*((1./lam)-(1./lambda))).^2);

     figure(253)
     plot(HSRLMolecular_scan_wavelength, squeeze(RD(10,1,:)))
     hold on
     plot(HSRLMolecular_scan_wavelength, squeeze(RD(10,1,:)).*eta_mol)
     hold off
     
     
     eta_mol_all = ones(length(time_comb), length(range));
% calucate the overall efficieny of the molecular channel
     for i = 1:size(Surf_T,1)       
       eta_mol_all(i,:) = trapz(squeeze(RD(i,:,:)).*eta_mol',2)./trapz(squeeze(RD(i,:,:)),2);
      % eta_mol_all(1)
     end

     
 % backscatter ratio using the median temp profile   
     BSR = (O2_offline_comb./eta_comb)./(O2_offline_mol./eta_mol_all);
 
 %beta_m_array = repmat(beta_m_profile, size(BSR,1), 1);
 %beta_bs = ((BSR)-1).*beta_m_profileay;
 beta_bs = ((BSR)-1).*beta_m_profile;

 %% mask data

 BSR(O2_offline_mol<.25)=nan;  % remove very low count rate molecular data
 beta_bs(O2_offline_mol<.25)=nan;  % remove very low count rate molecular data

 %% plot data

  figure(110)
  Z = real(double((((BSR))')));
  font_size = 14;
  xData =  linspace(fix(min(time_comb)),  ceil(max(time_comb)), 25);
  set(gcf,'renderer','zbuffer');
  h = pcolor(x,y,Z);
  set(h, 'EdgeColor', 'none');
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  set(gca, 'XTick',  xData)
  colorbar('EastOutside');
  axis([fix(min(time_comb)) fix(min(time_comb))+1 0 12])
  caxis([1e-1 1e3]);
  %caxis([0 20]);
  hh = title({[node, ' ', daystr, ' Backscatter Ratio']},'fontweight','b','fontsize',font_size);  
  datetick('x','HH','keeplimits', 'keepticks');
  xlabel('Time (UTC)','fontweight','b','fontsize',font_size);
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size);
  set(gca,'Fontsize',font_size,'Fontweight','b');
  set(gca,'Zscale', 'log')
  set(gca,'Colorscale', 'log')
  set(gca,'Zscale', 'linear')
  colormap(jet)
  
  figure(111)
  Z = real(double((beta_bs')));
  font_size = 14;
  xData =  linspace(fix(min(time_comb)),  ceil(max(time_comb)), 25);
  set(gcf,'renderer','zbuffer');
  h = pcolor(x,y,Z);
  set(h, 'EdgeColor', 'none');
  set(gca,'TickDir','out');
  set(gca,'TickLength',[0.005; 0.0025]);
  set(gca, 'XTick',  xData)
  colorbar('EastOutside');
  axis([fix(min(time_comb)) fix(min(time_comb))+1 0 12])
  caxis([1e-7 1e-3]);
  datetick('x','HH','keeplimits', 'keepticks');
  xlabel('Time (UTC)','fontweight','b','fontsize',font_size);
  ylabel('Height (km, AGL)','fontweight','b','fontsize',font_size);
  set(gca,'Fontsize',font_size,'Fontweight','b');
  set(gca,'Zscale', 'log')
  set(gca,'Colorscale', 'log')
  set(gca,'Zscale', 'linear')
  colormap(jet)
  
  
  

  
%    % test the profiles of the python to the matlab 
%    figure(20) 
%    [minValue, closestIndex] = min(abs(time_comb - datenum('21-Jul-2021 10:00:00')))
%    BSR_profile_mat = BSR(closestIndex, :); 
%    plot(BSR_profile_mat, (range./1e3));
%    ylim([0 6])
%    xlim([0 10])
   
%    load BSR_at_21Jun2021_1000
%    hold on
%    plot(BSR_profile, y);
%    hold off
%    title('BSR profile on 21-Jul-2021 10:00:00')
%    legend('Matlab BSR', 'Python BSR')
  
   serv_path1 = '/Volumes/Macintosh HD/Users/spuler/Desktop/mpd/';  
   cd(strcat(serv_path1,'/Plots'))
   
%    FigH = figure(20);
%    set(gca,'Fontsize',16,'Fontweight','b'); 
%    set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 600 600]);
%    name=strcat(' Backscatter Ratio comparison'); 
%    print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
   
   FigH = figure(110);
   set(gca,'Fontsize',16,'Fontweight','b'); 
   set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920/2 250]);
   name=strcat(node, "_", daystr, "_Backscatter Ratio");
   print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
   
   FigH = figure(111);
   set(gca,'Fontsize',16,'Fontweight','b'); 
   set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920/2 250]);
   name=strcat(node, "_", daystr, "_Backscatter Coefficient");
   print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
   
 
   
 %% save data
  
  if flag.save_data == 1
    cd(write_data_folder)
    save(name, 'x', 'y', 'beta_bs', 'BSR')
  end
   
 
 
  cd(d)



  
end






