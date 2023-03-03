function [T, P, BSR, RD, HSRLMolecular_scan_wavelength, const, beta_m_profile] = Process_HSRL_K_data(O2_online_comb, O2_offline_comb,...
    O2_online_mol,O2_offline_mol, time_comb, range, Surf_T, Surf_P, flag, node, daystr, Receiver_Scan_File, write_data_folder)
% Process the HSRL channels 
% Calculate the backscatter ratio (BSR) and backscatter coefficient using receiver scan data
 
 tic
 dd = pwd;

 % Receiver beam splitter sends more signal to the molecular detector vs combined. 
 % Determine split level using the measured O2 online counts in both channels 
 O2_geometric_correction = mean(O2_online_comb(1:end-15,:)./O2_online_mol(1:end-15,:),1,'omitnan');
 O2_geometric_correction = mean(O2_online_comb(1:end-15,:)./(O2_online_mol(1:end-15,:)+O2_online_comb(1:end-15,:)),1,'omitnan');
 O2_geometric_correction = median(O2_online_comb./(O2_online_mol+O2_online_comb),1,'omitnan');
   
 % trying to make this an automatic calculated value, but it blows up in clouds 
 O2_geometric_correction(O2_geometric_correction==0)=nan;
 O2_geometric_correction(O2_geometric_correction==Inf | O2_geometric_correction==-Inf)=nan;
 try
   % the geometric correction is ignored for now  
   eta_comb = median(O2_geometric_correction, 'omitnan')  
 catch
   warning('Problem with the geometric correction, assigning a fixed value.');
   eta_comb = 0.6  % override, and use this value for now 
 end
  
      if flag.troubleshoot == 1
          figure(102)
          sample_profile = 100;
          semilogx(O2_online_comb(sample_profile,:), range/1000, 'k')
          hold on
          semilogx(O2_offline_comb(sample_profile,:), range/1000, 'm')
          semilogx(O2_online_mol(sample_profile,:), range/1000, 'g')
          semilogx(O2_offline_mol(sample_profile,:), range/1000, 'b')
          ylim([0 6])
          legend('online comb', 'offline comb', 'online mol', 'offline mol')
          %semilogx(O2_online_comb(sample_profile,:).*O2_geometeric_correction, range/1000, 'k-');
          hold off

          %Combined efficiency, rel to molecular 
          figure(103)
          %plot(O2_online_comb(100,:)./O2_online_mol(100,:), range/1000)
          %hold on
          plot(O2_geometric_correction, range/1000)
          xlim([0 2])
          ylim([0 8])
          %hold off

        % plot relative backscatter
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
      end
  
%% Calculate Backscatter from Molecules

 const.k_B = 1.380649e-23; % (J/K, or Pa m^3/K)
 const.K_B = const.k_B * 9.869e-6; % (atm m^3/K)
 const.N_A = 6.022E23; % (/mol) Avagadros number
 const.M = 28.97E-3; % (kg/mol) air molecular mass per mol 
 const.m = const.M./const.N_A; % (kg) mass of a air molecule
 const.g = 9.80665; % (m/s^2) gravitational constant
 const.c = 299792458; % (m/s) (exact)  
 const.R = const.k_B*const.N_A; % (J/K mol) universal gas constant
 lambda = 770.1085*1E-9; % wavelength in nm (should really get this elsewhere)
 
 % Calculate temperature and pressure profile based on surface measurement and 
 % assuming a standard lapse rate (-6.5 deg/km) for the entire troposphere
 lapse = 0.0065; % (K/m) standard atmosphere lapse rate, dry adiabatic is 9.8 K/km
 if flag.WS == 1  % use surface values if they exist
   T0 = median(Surf_T,'omitnan')+273.15
   P0 = median(Surf_P,'omitnan')
   Surf_T(isnan(Surf_T))= T0; % fills in missing values with median
   Surf_P(isnan(Surf_P))= P0; % fills in missing values with median
   if isnan(T0)==1
       T0=25+273.15;
       Surf_T= ones(size(Surf_T)).*T0;
       warning('No Temperature Weather Station Data');
  end
  if isnan(P0)==1
       P0=1;
       Surf_P= ones(size(Surf_P)).*P0;
       warning('No Pressure Weather Station Data');
  end 
   T = (Surf_T+273.15)-lapse.*range;
   P = Surf_P.*((Surf_T+273.15)./T).^-((const.M*const.g)/(const.R*lapse));   % barometric formula
 else
   T0 = 290; % surface temperature
   P0  = 1; % surface pressure in atm 
   T = T0-lapse.*range; % temperature as function of range with lapse rate
   P = P0.*(T0./T).^-((const.M*const.g)/(const.R*lapse));   % barometric formula 
 end
 
%% calculate HSRL data products 

 % only a fraction of total molecular gets through the K filter 
 % need to know filter width, molecular scatter width RB which is range dependent
 % just assume a constant for now
 eta_mol_all = 0.167;

 % below is a more complete way to check this efficiency 
 % open the receiver scan data and read in receiver scans
 if strcmp(getenv('HOSTNAME'),'fog.eol.ucar.edu') == 1
   cal_path = '/home/rsfdata/Processing'; % when running on server
 else
   cal_path = '/Users/spuler/Documents/GitHub'; % when running on server s
 end
 cd ([cal_path '/eol-lidar-calvals/calfiles'])
 
 cal_file = Receiver_Scan_File;
 ncid = netcdf.open(cal_file, 'NC_NOWRITE');
    % ncdisp(cal_file) % use this to display all variables
    HSRLCombined_Transmission = ncread(cal_file, 'HSRLCombined_Transmission'); 
    HSRLCombined_scan_wavelength = ncread(cal_file, 'HSRLCombined_scan_wavelength').*1e-9; 
    HSRLMolecular_Transmission = ncread(cal_file, 'HSRLMolecular_Transmission'); 
    HSRLMolecular_scan_wavelength = ncread(cal_file, 'HSRLMolecular_scan_wavelength').*1e-9; 
 netcdf.close(ncid); 
 cd(dd)
  
 % sometimes the transmissions are shifted several bins?  Just a check on correcting  
 %HSRLMolecular_Transmission = circshift(HSRLMolecular_Transmission,-5);
 %HSRLCombined_Transmission = circshift(HSRLCombined_Transmission,5);

 % this should be the same eta_comb ratio as measured above but not
 % necessarily if the receiver has been adjusted since the calibration
 eta_comb_cal = median(HSRLCombined_Transmission./(HSRLCombined_Transmission+HSRLMolecular_Transmission),1,'omitnan');
 eta_mol = (HSRLMolecular_Transmission./(1-eta_comb_cal))./(HSRLCombined_Transmission./eta_comb_cal);
  
%     if flag.troubleshoot == 1
       figure(252)
      % plot(HSRLMolecular_scan_wavelength, HSRLMolecular_Transmission)
      % hold on
      % plot(HSRLCombined_scan_wavelength, HSRLCombined_Transmission)
       plot(HSRLCombined_scan_wavelength, HSRLCombined_Transmission/eta_comb_cal)
       hold on
       plot(HSRLMolecular_scan_wavelength, HSRLMolecular_Transmission./(1-eta_comb_cal))
       hold off
       figure(253)
       plot(HSRLMolecular_scan_wavelength, eta_mol)
%     end
      
 % Calculate the Rayleigh-Doppler broadening for each height 
 % from Fiocco and DeWolf 1968
  K0 = 2*pi/(lambda)/const.c;
  K = K0 + K0; % in backscatter 
  lam =  HSRLMolecular_scan_wavelength; 
  lam = reshape(lam,1,1,length(lam));  % dimension of time, range, spectral wavelength

  RD  = sqrt(const.m./(2*pi*K^2*const.k_B.*T))...
        .*exp((-1*const.m./(2*K^2*const.k_B.*T)).*(2*pi*((1./lam)-(1./lambda))).^2); 
 % Bosenberg 1998 suggests a simple Brillouin broadening correction: RD width x 1.2 ~ RDB width  
 % (valid in the lower 10km of a standard atmosphere)
  RD = RD.*1.2;  
  
  %   if flag.troubleshoot == 1  
         % Plot RD spectrum and receiver transmission
         figure(254)
         plot(HSRLMolecular_scan_wavelength, squeeze(RD(10,1,:)))
         hold on
         plot(HSRLMolecular_scan_wavelength, squeeze(RD(10,1,:)).*eta_mol)
         legend('RD', 'RD*receiver_eff') 
         hold off
  %   end
     
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
 threshold = 0.25;
 BSR(O2_offline_mol< threshold)=nan;  % remove very low count rate molecular data
 beta_bs(O2_offline_mol< threshold)=nan;  % remove very low count rate molecular data


 toc
 %% plot data

  figure(110)
  x = (time_comb)';
  y = (range./1e3);
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
  
    
%% save the plots 
  
%    serv_path1 = '/Volumes/Macintosh HD/Users/spuler/Desktop/mpd/';  
%    cd(strcat(serv_path1,'/Plots'))
%    
%    FigH = figure(110);
%    set(gca,'Fontsize',16,'Fontweight','b'); 
%    set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920/2 250]);
%    name=strcat(node, "_", daystr, "_Backscatter Ratio");
%    print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
%    
%    FigH = figure(111);
%    set(gca,'Fontsize',16,'Fontweight','b'); 
%    set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1920/2 250]);
%    name=strcat(node, "_", daystr, "_Backscatter Coefficient");
%    print(FigH, name, '-dpng', '-r0') % set at the screen resolution 
   
 
   
%  %% save data
  name=strcat(node, "_", daystr, "_Backscatter Coefficient");
  if flag.save_data == 1
    cd(write_data_folder)
    save(name, 'x', 'y', 'beta_bs', 'BSR')
  end
   
 
cd(dd)

  
end






