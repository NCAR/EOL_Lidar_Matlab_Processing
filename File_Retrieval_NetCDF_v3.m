function[online_merged,offline_merged,folder,flag] = File_Retrieval_NetCDF_v3(flag, bins, folder)

dd = pwd; % get the current path
%cd /scr/eldora1/wvdial_2_data/2018
%folder = '20180818';
cd(folder)
MCS.dirListing = dir(strcat('MCSsample','*')); %ncdisp('MCSsample000000.nc', '/', 'min') 
LL.dirListing = dir(strcat('LLsample','*')); %ncdisp('LLsample000000.nc', '/', 'min')
%Etalon.dirListing = dir(strcat('Etalonsample','*')); %ncdisp('Etalonsample000000.nc', '/', 'min')
HKeep.dirListing = dir(strcat('HKeepsample','*')); %ncdisp('HKeepsample000000.nc', '/', 'min')
Pow.dirListing = dir(strcat('Powsample','*')); %ncdisp('Powsample000000.nc', '/', 'min')
WS.dirListing = dir(strcat('WSsample','*')); %ncdisp('WSsample000000.nc', '/', 'min')

d=1;

for d = 1:length(MCS.dirListing)
  %read in the MCS photo count data  
  MCS.filename=MCS.dirListing(d).name;
  MCS.data = ncread(MCS.filename,'Data'); 
  MCS.time = ncread(MCS.filename,'time'); 
  MCS.channel = ncread(MCS.filename,'Channel'); 
  if d>1   % sum the days nc data files into a single array
    MCS.time1=[MCS.time1;MCS.time];
    MCS.data1=[MCS.data1;MCS.data];
    MCS.channel1=[MCS.channel1;MCS.channel];
  else
    MCS.time1=MCS.time;
    MCS.data1=MCS.data;
    MCS.channel1=MCS.channel;
  end
end

for d = 1:length(LL.dirListing)
  %read in the laser locking wavelength data 
  LL.filename=LL.dirListing(d).name;
  LL.time = ncread(LL.filename,'time');
  LL.wavelength = ncread(LL.filename,'Wavelength');
  LL.current =  ncread(LL.filename,'Current');
  LL.name = h5read(LL.filename,'/LaserName'); % note the special type
  if d>1   % sum the days nc data files into a single array
    LL.time1=[LL.time1;LL.time];
    LL.wavelength1=[LL.wavelength1;LL.wavelength];
    LL.current1=[LL.current1;LL.current];
    LL.name1 = [LL.name1;LL.name];
  else
    LL.time1 = LL.time;
    LL.wavelength1 = LL.wavelength;
    LL.current1 = LL.current;
    LL.name1 = LL.name;
  end
end


for d = 1:length(HKeep.dirListing)
  %read in the HK data
  HKeep.filename = HKeep.dirListing(d).name;
  HKeep.time = ncread(HKeep.filename,'time');
  HKeep.temp = ncread(HKeep.filename,'Temperature');
  if d>1   % sum the days nc data files into a single array
    if size(HKeep.temp,2)== size(HKeep.temp1,2) % catch addition of thermocouple error
      HKeep.time1 = [HKeep.time1; HKeep.time];
      HKeep.temp1 = [HKeep.temp1; HKeep.temp];
    end
  else
      HKeep.time1 = HKeep.time;
      HKeep.temp1 = HKeep.temp;
  end
end

for d = 1:length(Pow.dirListing)
  %read in the power montioring data
  Pow.filename = Pow.dirListing(d).name;
  Pow.time = ncread(Pow.filename,'time');
  Pow.power = ncread(Pow.filename,'Power');
  Pow.channel = h5read(Pow.filename,'/ChannelAssignment');
  if d>1   % sum the days nc data files into a single array
    Pow.time1 = [Pow.time1; Pow.time];
    Pow.power1 = [Pow.power1; Pow.power];
  else
    Pow.time1 = Pow.time;
    Pow.power1 = Pow.power;
  end
end

for d = 1:length(WS.dirListing)
  %read in the weather station data
  WS.filename = WS.dirListing(d).name;
  WS.time = ncread(WS.filename,'time');
  WS.temp = ncread(WS.filename,'Temperature');
  WS.press = ncread(WS.filename,'Pressure');
  WS.relhum = ncread(WS.filename,'RelHum');
    if d>1   % sum the days nc data files into a single array
    WS.time1 = [WS.time1; WS.time];
    WS.temp1 = [WS.temp1; WS.temp];
    WS.press1 = [WS.press1; WS.press];
    WS.relhum1 = [WS.relhum1; WS.relhum];
  else
    WS.time1 = WS.time;
    WS.temp1 = WS.temp;
    WS.press1 = WS.press;
    WS.relhum1 = WS.relhum;
  end
end

% combine time and data 
MCS.all = [MCS.time1,MCS.data1];
% parse out the channels
MCS.online = MCS.all(MCS.channel1==0,:); % 0 is the online
MCS.offline = MCS.all(MCS.channel1==8,:); % 8, first demux channel is offline
MCS.combined = MCS.all(MCS.channel1==2,:); % 2 is the combined HSRL channel
MCS.molecular = MCS.all(MCS.channel1==3,:); % 3 is the molecular HSRL channel

if isfield(LL,'time') == 1
  LL.all = [LL.time1, LL.wavelength1, LL.current1];
  LL.online = LL.all(strcmp(LL.name1,'WVOnline'),:);
  LL.offline = LL.all(strcmp(LL.name1,'WVOffline'),:);
end
%LL.hsrl =  LL.all(strcmp(LL.name1,'HSRL'),:);
if isfield(HKeep,'time') == 1
  HKeep.all = [HKeep.time1, HKeep.temp1];
end
if isfield(Pow,'time') == 1
  Pow.all = [Pow.time1, Pow.power1];
end
if isfield(WS,'time') == 1
  WS.all = [WS.time1, WS.temp1, WS.press1, WS.relhum1];
else
  flag.WS = 0
end

% grid to a fixed time base
ave_time = 2;  % time grid in seconds
time_grid = (floor(min(MCS.time1)):1/60/60*(ave_time):ceil(max(MCS.time1)))';
LL.online_grid = interp1(LL.online(:,1), LL.online(:,2:end), time_grid, 'nearest', 'extrap'); 
LL.offline_grid = interp1(LL.offline(:,1), LL.offline(:,2:end), time_grid, 'nearest', 'extrap'); 
%LL.hsrl_grid = interp1(LL.hsrl(:,1), LL.hsrl(:,2:end), time_grid, 'nearest', 'extrap');
[x, ia, ic] = unique(MCS.online(:,1),'rows');  % remove duplicate time points
uA = MCS.online(ia,:);  % apply those to the other rows
MCS.online_grid = interp1(uA(:,1), uA(:,2:end), time_grid, 'nearest', 'extrap'); 
[x, ia, ic] = unique(MCS.offline(:,1),'rows');  % remove duplicate time points
uA = MCS.offline(ia,:);  % apply those to the other rows
MCS.offline_grid = interp1(uA(:,1), uA(:,2:end), time_grid, 'nearest', 'extrap'); 
%MCS.combined_grid = interp1(MCS.combined(:,1), MCS.combined(:,2:end), time_grid, 'nearest', 'extrap'); 
%MCS.molecular_grid = interp1(MCS.molecular(:,1), MCS.molecular(:,2:end), time_grid, 'nearest', 'extrap');
if isfield(HKeep,'time') == 1
  HKeep.temp1 = interp1(HKeep.all(:,1), HKeep.all(:,2), time_grid, 'nearest', 'extrap');
  if size(HKeep.all,2)==3
    HKeep.temp2 = interp1(HKeep.all(:,1), HKeep.all(:,3), time_grid, 'nearest', 'extrap');
  else
    HKeep.temp2 = interp1(HKeep.all(:,1), HKeep.all(:,2), time_grid, 'nearest', 'extrap');
  end
end
[x, ia, ic] = unique(Pow.all(:,1),'rows');  % remove duplicate time points
uA = Pow.all(ia,:);  % apply those to the other rows
Pow.online = interp1(uA(:,1), uA(:,2), time_grid, 'nearest', 'extrap');
Pow.offline = interp1(uA(:,1), uA(:,8), time_grid, 'nearest', 'extrap');
%Pow.online = HKeep.temp1;
%Pow.offline = HKeep.temp2;
if isfield(WS,'time') == 1
  WS.online = interp1(WS.all(:,1), WS.all(:,2:end), time_grid, 'nearest', 'extrap');
end

% combine all of the data into original data format
if isfield(HKeep,'time') == 1 && isfield(WS,'time') == 1
  online_merged = [time_grid, LL.online_grid, Pow.online, HKeep.temp1, WS.online, MCS.online_grid];
  offline_merged = [time_grid, LL.offline_grid, Pow.offline,  HKeep.temp2, WS.online, MCS.offline_grid];
else
  online_merged = [time_grid, LL.online_grid, Pow.online, MCS.online_grid];
  offline_merged = [time_grid, LL.offline_grid, Pow.offline, MCS.offline_grid];
end
%hsrl_com_merged = [time_grid, LL.hsrl_grid, Pow.online, WS.online, MCS.combined_grid];
%hsrl_mol_merged = [time_grid, LL.hsrl_grid, HKeep.temp1, WS.online, MCS.molecular_grid];

cd(dd) % point back to original directory

end

