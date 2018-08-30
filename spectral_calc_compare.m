oxygen = 0;
width = 0;
cd /Users/spuler/Desktop
if oxygen == 1   
  fileID(1) = fopen('Spectral_Calc/O2_3km_500m_path_0K_offset.txt'); % normal temp file
  if width == 1
    fileID(2) = fopen('Spectral_Calc/O2_3km_500m_path_0K_offset_1.8GHz.txt'); % Gaussian "Instrument Function" 
  else
    fileID(2) = fopen('Spectral_Calc/O2_3km_500m_path_10K_offset.txt'); % high temp file
  end
  %fileID(2) = fopen('Spectral_Calc/O2_20km_500m_path_50K_offset.txt'); % high alt file
else
  fileID(1) = fopen('Spectral_Calc/H2O_3km_500m_path_0K_offset.txt'); % normal temp file
  fileID(2) = fopen('Spectral_Calc/H2O_3km_500m_path_10K_offset.txt'); % high temp file
  %fileID(2) = fopen('Spectral_Calc/H2O_20km_500m_path_50K_offset.txt'); % high alt file
end

i=1;
for i = 1:2
  C(i) = textscan(fileID(i), '%q', 88);
  if width == 1
    if i == 1
      C(i) = textscan(fileID(i), '%q', 88);
    else
      C(i) = textscan(fileID(i), '%q', 99);
    end
  end
  Header = (C{i});
  IndexC = strfind(C{i}, 'Temperature');
  Index = find(not(cellfun('isempty', IndexC)));
  T = str2num(Header{Index+3,1}) % temperaute (K)
  IndexC = strfind(C{i}, 'Pressure');
  Index = find(not(cellfun('isempty', IndexC)));
  P = str2num(Header{Index+3,1})*0.000986923 % pressure (convert from mb to atm)
  IndexC = strfind(C{i}, 'Isotopes');  
  Index = find(not(cellfun('isempty', IndexC)));
  m = str2num(Header{Index+1,1}) % mixing ratio
  IndexC = strfind(C{i}, 'Length');  
  Index = find(not(cellfun('isempty', IndexC)));
  l = str2num(Header{Index+3,1})/100 % length (convert from cm to m)
  data = fscanf(fileID(i), '%f %f', [2 Inf])';
  fclose(fileID(i));
  x{i} = data(:,1).*1000; % wavelength in nm
  Trn{i} = data(:,2); % Transmittance
  Abs{i} = 1-Trn{i}; % Absorption
  k_B = 1.3626E-28; % Boltzmann constant in atm K-1 m^3:  
  N = P/(k_B*T) % number density of all air molecules based on T and P
  sigma{i}=-1*log(Trn{i})/l/N/m*100^2;  %convert transmission to abs cross section in cm2
end

% read in text files from spectral cal 
% x1,y1, O2 absorption at 1km alt, 1km path at standard atmosphere 
% x2,y2, O2 absorption at 1km alt, 1km path at standard atmosphere + 10K
% the spectral spacing is different for each text file so interpolate to commmon grid 

x_grid=(x{1}(1):0.0001:x{1}(end))';
sigma1_grid=interp1(x{1},sigma{1},x_grid);
sigma2_grid=interp1(x{2},sigma{2},x_grid);
Trn1_grid=interp1(x{1},Trn{1},x_grid, 'spline');
Trn2_grid=interp1(x{2},Trn{2},x_grid, 'spline');
Abs1_grid=interp1(x{1},Abs{1},x_grid, 'spline');
Abs2_grid=interp1(x{2},Abs{2},x_grid, 'spline');
rel_diff_grid = ((Abs2_grid)-(Abs1_grid))./(Abs1_grid); % relative diff Abs
diff_grid=((Abs2_grid)-(Abs1_grid)); % delta Abs
rel_diff_grid = (log(Trn2_grid)-log(Trn1_grid))./log(Trn1_grid); % relative diff n

figure(100)
AH_T0 = (-log(Trn1_grid)./sigma1_grid./100./l).*1e6./6.022E23.*18.015;
AH_T = (-log(Trn2_grid)./sigma1_grid./100./l).*1e6./6.022E23.*18.015;
AH_diff = AH_T-AH_T0;
plot(x_grid,AH_T0, 'k')
hold on
plot(x_grid,AH_T, 'r')
hold off
%fontsize(12)
%figure(2)
%semilogy(x{1},sigma{1},x{2},sigma{2}) 
%figure(3)
%plot(x_grid,diff_grid)

scrsz = get(0,'ScreenSize');
figure1 = figure('visible', 'on','Position',[scrsz(4)/1.5 scrsz(4)/10 scrsz(3)/1.5 scrsz(4)/1.5]);
subplot1=subplot(3,1,1,'Parent',figure1);
%plot(x{1},Abs{1},x{2},Abs{2}) % plot absorption
%hold on
plot(x{1},-log(Trn{1}),x{2},-log(Trn{2})) % plot absorption
grid on
if oxygen ==1
  title('Oxygen (3km AGL, 500 path)')
else
  title('Water vapor (3km AGL, 500 path)')
end
ylabel('Absorption')
xlabel('wavelength (nm)')
if width == 1
  legend('narrow laser','molecular scatter')
else
  legend('{\Delta}0K','{\Delta}10K')
end
subplot1=subplot(3,1,2,'Parent',figure1);
semilogy(x{1},sigma{1},x{2},sigma{2})

grid on
ylabel('Absorption Cross Section cm{^2}')
xlabel('wavelength (nm)')
if width == 1
  legend('narrow laser','molecular scatter')
else
  legend('{\Delta}0K','{\Delta}10K')
end
subplot1=subplot(3,1,3,'Parent',figure1);
%plot(x_grid,diff_grid, 'k')
%ylabel('{\Delta} Abs')
plot(x_grid,rel_diff_grid, 'k')
ylabel('{\Delta}A/A')
%plot(x_grid,AH_diff, 'k')
%ylabel('{\Delta}AH (g m^{-3})')
%xlabel('wavelength (nm)')
grid on

%link the x axis for all 2 subplots
 ax(1)=subplot(3,1,1);
 ax(2)=subplot(3,1,2);
 ax(3)=subplot(3,1,3);
 linkaxes([ax(1),ax(2),ax(3)],'x')
if oxygen ==1
  xlim([769.76 769.84])
else
  xlim([828.15 828.23])
end

%mean(x_grid(1969:1989))
%1-mean(y2_grid(1969:1989))
%x_grid(1969)
%1-y2_grid(1969)
%x_grid(1979)
%1-y2_grid(1979)
%x_grid(1989)
%1-y2_grid(1989)

%Q=trapz(x_grid(1969:1999), %integrate from 828.197 to 828.199 