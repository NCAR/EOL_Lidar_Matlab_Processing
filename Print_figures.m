% printing figures

scrsz = [1  1  1920 1200]
Scrsize=[scrsz(4)/1 scrsz(4)/1 scrsz(3)/1.5 scrsz(4)/1.5];

%cd('/Volumes/documents/WV_DIAL_data/plots/') % point to the directory where data is stored 
cd('/Users/spuler/Desktop') % point to the directory where data is stor
FigH = figure(1);
%set(gca,'Fontsize',30,'Fontweight','b'); % 
set(FigH, 'PaperUnits', 'points', 'PaperPosition', Scrsize);
name=strcat('Figure_8');
print(FigH, name, '-dpng', '-r300')  