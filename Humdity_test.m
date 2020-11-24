%clear all 
%close all


% import Dranetz data as a matrix and save as 'D' 
% open D, Hum, and WS data

time_unix = D(:,1); % the first column is unit time
time = time_unix/86400+datenum(1970,1,1,6,0,0); % convert from unix time and correct for MST DST

scrsz = get(0,'ScreenSize');
Scrnsize=[scrsz(4)/2 scrsz(4)/10 scrsz(3)/1 scrsz(4)/2];

figure1 = figure('Position',Scrnsize);
   subplot1=subplot(2,1,1,'Parent',figure1,'YGrid','on', 'XGrid','on');
   box(subplot1,'on');
   hold(subplot1,'all');  
   
%figure(1)
plot(time, D(:,17), 'k', 'DisplayName', 'HVAC current') % HVAC
hold on
%plot(time, D(:,14), 'g') % instrument
plot(time, D(:,23), 'r', 'DisplayName', 'Wall Heater current')  % wall heater
%plot(time, D(:,20), 'b') % window heater
hold off

axis([datenum(2020, 9, 30, (14), 0 ,0) datenum(2020, 10, 1, (7), 0, 0) 0 60])
datetick('x','mmm dd yyyy HH','keeplimits', 'keepticks');
%xlabel('date')
ylabel('current (A)')

FigH = figure(1);
set(gca,'Fontsize',18,'Fontweight','b');


%Humidity = [Humidity_11; ...
%            Humidity_12; Humidity_13; Humidity_14; Humidity_15;...
%            Humidity_16; Humidity_17; Humidity_18; Humidity_19;...
%            Humidity_20; Humidity_21; Humidity_22; Humidity_23;...
%            Humidity_24; Humidity_25; Humidity_26; Humidity_27;... 
%            Humidity_28; Humidity_29; Humidity_30]; 


time_0 = datenum(2020,09,30, 11, 0, 0); % convert from unix time
time_h = time_0:(1/893/24):time_0+(size(Humidity,1)-1)/893/24;

%WS = [WS_11; ...
%      WS_12; WS_13; WS_14; WS_15;...
%      WS_16; WS_17; WS_18; WS_19;...
%      WS_20; WS_21; WS_22; WS_23;...
%      WS_24; WS_25; WS_26; WS_27;... 
%      WS_28; WS_29; WS_30]; 
  
time_ws = time_0:(1/715/24):time_0+(size(WS,1)-1)/715/24;

%figure(2)
subplot2=subplot(2,1,2,'Parent',figure1,'YGrid','on', 'XGrid','on');
   box(subplot2,'on');
   hold(subplot2,'all');  

plot(time_ws, WS(:,1), 'k', 'DisplayName', 'outside temp') % outside temp
hold on
%plot(time_h, Humidity(:,1), 'r') % Hum internal temp
plot(time_h, Humidity(:,2), 'r', 'DisplayName', 'internal temp') % Hum internal temp
%plot(time_h, Humidity(:,3), 'b') % RH
plot(time_h, Humidity(:,4), 'g', 'DisplayName', 'dew point') % dew point
hold off

axis([datenum(2020, 9, 30, (14), 0 ,0) datenum(2020, 10, 1, (7), 0, 0) 0 30])
datetick('x','mmm dd yyyy HH','keeplimits', 'keepticks');
%xlabel('date')
ylabel('temp (C)')


FigH = figure(1);
set(gca,'Fontsize',18,'Fontweight','b')

%   % Create legend
%   legend(subplot1,'show'); 
%   legend(subplot2,'show');


  FigH = figure(1);
  set(gca,'Fontsize',18,'Fontweight','b');
  set(FigH, 'PaperUnits', 'points', 'PaperPosition', Scrnsize);
  print(FigH, 'Humidity_test', '-dpng', '-r0') % set at the screen resolution




%
figure(10)
plot(time, D(:,17), 'k') % HVAC
hold on
plot(time_ws, WS(:,1), 'k') % outside temp
plot(time_h, Humidity(:,3), 'b') % RH
plot(time_h, Humidity(:,4), 'g', 'DisplayName', 'dew point') % dew point
plot(time, D(:,23), 'r')  % wall heater
hold off
axis([datenum(2020, 9, 30, (14), 0 ,0) datenum(2020, 10, 1, (7), 0, 0) 0 90])
datetick('x','mmm dd yyyy HH','keeplimits', 'keepticks');