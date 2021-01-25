clear all
close all

O_x = [50;100;200;300;400;500;750;1000;1250;1500; 2000;3000;4000;5000;6000;8000;12000]; %range in meters
O_y_near = [2.47E-2; 9.90E-2; 3.99E-1; 8.72E-1; 1.00E+0; 1.00E+0 ;1.00E+0; 1.00E+0; 1.00E+0; 1.00E+0; 1.00E+0; 1.00E+0; 1.00E+0; 1.00E+0; 1.00E+0; 1.00E+0; 1.00E+0]; %near range receiver overlap function
O_y_primary = [7e-7; 1.5e-5; 2.77e-4; 1.38e-3; 3.98e-3; 8.89e-3; 3.72e-2; 1.06e-1; 2.08e-1; 3.39e-1; 6.61e-1; 9.53e-1; 9.74e-1; 9.86e-1; 9.92e-1; 1.00E+0; 1.00E+0]; %primary receiver overlap function

O_y_combined = O_y_primary*0.9 + O_y_near*.1;  

figure(1)
loglog(O_y_primary, O_x, 'b', 'DisplayName','Primary', 'LineWidth', 2)
hold on
loglog(O_y_near, O_x, 'r', 'DisplayName','Near', 'LineWidth', 2)
loglog(O_y_combined, O_x, 'k', 'DisplayName','Combined', 'LineWidth', 2)
grid on
ylabel('Range (m)'); 
xlabel('Overlap function'); 
legend('Location', 'NorthWest');
axis([1e-7 2e0 50 12000])

scrsz = [1  1  1920 1200]
Scrsize=[scrsz(4)/1 scrsz(4)/1 scrsz(3)/1.5 scrsz(4)/1.5];

%cd('/Volumes/documents/WV_DIAL_data/plots/') % point to the directory where data is stored 
cd('/Users/spuler/Desktop') % point to the directory where data is stor
FigH = figure(1);
set(gca,'Fontsize',30,'Fontweight','b'); % 
set(FigH, 'PaperUnits', 'points', 'PaperPosition', Scrsize);
name=strcat('Gen5_overlap_function');
print(FigH, name, '-dpng', '-r300') 