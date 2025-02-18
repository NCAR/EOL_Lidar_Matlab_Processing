start = 1116;
stop = 1242;

figure(2)
plot(range, sum(Online_Temp_Spatial_Avg(start:stop,:)),'r')
hold on
plot(range, sum(Offline_Temp_Spatial_Avg(start:stop,:)),'k')
xlim([0 18000])
grid on

figure(3)
plot(range, mean(Online_Temp_Spatial_Avg(start:stop,:)),'r')
hold on
plot(range, mean(Offline_Temp_Spatial_Avg(start:stop,:)),'k')
xlim([0 18000])
grid on

%QE 
figure(4)
plot(range, 1+(mean(Online_Temp_Spatial_Avg(start:stop,:))./mean(Online_Temp_Spatial_back(start:stop,:))),'r')
hold on
plot(range, 1+(mean(Offline_Temp_Spatial_Avg(start:stop,:))./mean(Offline_Temp_Spatial_back(start:stop,:))),'k')
xlim([0 18000])
grid on

figure(5)
plot(range, mean(Online_Temp_Spatial_Avg_c(start:stop,:)),'r')
hold on
plot(range, mean(Offline_Temp_Spatial_Avg_c(start:stop,:)),'k')
xlim([0 18000])
grid on