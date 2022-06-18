% printing figures

scrsz = [1  1  1920 1200]
Scrsize=[scrsz(4)/1 scrsz(4)/1 scrsz(3)/1.5 scrsz(4)/1.5];

%cd('/Volumes/documents/WV_DIAL_data/plots/') % point to the directory where data is stored 
cd('/Users/spuler/Desktop') % point to the directory where data is stor
FigH = figure(1);
%set(gca,'Fontsize',30,'Fontweight','b'); % 
set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 3200/1.5 2000/1.5]);
name=strcat('WV_Beta_T_comparison');
print(FigH, name, '-dpng', '-r300')  



%cd('/Volumes/documents/WV_DIAL_data/plots/') % point to the directory where data is stored 
cd('/Users/spuler/Desktop') % point to the directory where data is stor
FigH = figure(202);
%set(gca,'Fontsize',30,'Fontweight','b'); % 
set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 800 800]);
name=strcat('Corr_coeff_vs_range_Gen4_to_Gen5_matlab');
print(FigH, name, '-dpng', '-r300')  


%cd('/Volumes/documents/WV_DIAL_data/plots/') % point to the directory where data is stored 
cd('/Users/spuler/Desktop') % point to the directory where data is stor
FigH = figure(1);
%set(gca,'Fontsize',30,'Fontweight','b'); % 
set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1280 800]);
name=strcat('MPD_Sonde_SGP');
print(FigH, name, '-dpng', '-r300')  


%cd('/Volumes/documents/WV_DIAL_data/plots/') % point to the directory where data is stored 
cd('/Users/spuler/Desktop') % point to the directory where data is stor
FigH = figure(28);
set(gca,'Fontsize',24,'Fontweight','b'); % 
set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 800 800]);
name=strcat('MPD_Net_Demo_Raman_sonde');
print(FigH, name, '-dpng', '-r300')  


%cd('/Volumes/documents/WV_DIAL_data/plots/') % point to the directory where data is stored 
cd('/Users/spuler/Desktop') % point to the directory where data is stor
FigH = figure(1);
%set(gca,'Fontsize',16,'Fontweight','b'); % 
set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1280 800]);
name=strcat('5_unit_intercomparison');
print(FigH, name, '-dpng', '-r300') 

 
%cd('/Volumes/documents/WV_DIAL_data/plots/') % point to the directory where data is stored 
cd('/Users/spuler/Desktop') % point to the directory where data is stor
FigH = figure(1);
set(gca,'Fontsize',16,'Fontweight','b'); % 
set(FigH, 'PaperUnits', 'points', 'PaperPosition', [1 1 1600 250]);
name=strcat('Raman_SGP');
print(FigH, name, '-dpng', '-r300')  

