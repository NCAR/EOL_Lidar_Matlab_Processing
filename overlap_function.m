clear all
close all
O_x= [50;100;200;300;400;500;750;1000;1250;1500; 2000;3000;4000;5000;6000;8000;12000]; 
O_y_near = [2.47E-2; 9.90E-2; 3.99E-1; 8.72E-1; 1.00E+0; 1.00E+0 ;1.00E+0; 1.00E+0; 1.00E+0; 1.00E+0; 1.00E+0; 1.00E+0; 1.00E+0; 1.00E+0; 1.00E+0; 1.00E+0; 1.00E+0]; %near range receiver overlap function
O_y_primary = [7e-7; 1.5e-5; 2.77e-4; 1.38e-3; 3.98e-3; 8.89e-3; 3.72e-2; 1.06e-1; 2.08e-1; 3.39e-1; 6.61e-1; 9.53e-1; 9.74e-1; 9.86e-1; 9.92e-1; 1.00E+0; 1.00E+0]; %primary receiver overlap function

figure(1)
lg = loglog(O_y_primary, O_x, 'b',  'DisplayName', 'Primary')
lg(1).LineWidth = 2;
hold on
lg = loglog(O_y_near, O_x, 'r', 'DisplayName', 'Near')
lg(1).LineWidth = 2;
%loglog((O_y_near*.1+O_y_primary*.9), O_x, 'k', 'DisplayName', 'Combined')
grid on
legend('Location','northwest')
xlabel('Overlap function')
ylabel('Range (m)')
xlim([1e-7  2])
ylim([50 12000])
hold off

% the sigal is relative to the overlap function, O, and the Area
A_1 = 935; %cm^2 primary
A_2 = 14; %cm^2 WFOV near
eta = 0.90
A_near_p = A_2/A_1


figure(2)
lg = loglog(O_y_primary, O_x, 'b',  'DisplayName', 'Primary');
lg(1).LineWidth = 2;
hold on
lg = loglog(O_y_near*A_near_p, O_x, 'r', 'DisplayName', 'Near')
lg(1).LineWidth = 2;
lg = loglog((1/A_1)*(O_y_primary*A_1*eta + O_y_near*A_2*(1-eta)), O_x, 'k', 'DisplayName', 'Combined')
lg(1).LineWidth = 2;
grid on
legend('Location','northwest')
xlabel('Relative signal')
ylabel('Range [m]')
xlim([1e-7  2])
ylim([50 12000])
hold off

  FigH = figure(1);
   set(gca,'Fontsize',24,'Fontweight','b'); 
   set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1  1  800 600]);
   name=strcat('overlap_function'); 
   print(FigH, name, '-dpng', '-r0') % set at the screen resoluti


  FigH = figure(2);
   set(gca,'Fontsize',24,'Fontweight','b'); 
   set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1  1  800 600]);
   name=strcat('overlap_relative_signal'); 
   print(FigH, name, '-dpng', '-r0') % set at the screen resoluti

