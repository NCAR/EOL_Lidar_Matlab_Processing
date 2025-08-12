function[HSRL_Mol_merged, HSRL_Combined_merged, MCS] = MPD_File_Retrieval_NetCDF_HSRL(flag, MCS, folder, read_time_in)


dd = pwd; % get the current path
%cd /scr/eldora1/wvdial_2_data/2018
%folder = '20180818';
cd(folder)
MCSsample.dirListing = dir(strcat('MCS','*')); %ncdisp('MCS_03_20200519_010000.nc', '/', 'min') 
LL.dirListing = dir(strcat('LL','*')); %ncdisp('LLsample000000.nc', '/', 'min')
%Etalon.dirListing = dir(strcat('Etalonsample','*')); %ncdisp('Etalonsample000000.nc', '/', 'min')
HKeep.dirListing = dir(strcat('HKeep','*')); %ncdisp('HKeepsample000000.nc', '/', 'min')
Pow.dirListing = dir(strcat('Pow','*')); %ncdisp('Powsample000000.nc', '/', 'min')
WS.dirListing = dir(strcat('WS','*')); %ncdisp('WSsample000000.nc', '/', 'min')
Etalon.dirListing = dir(strcat('Etalon','*')); %ncdisp('Etalon_03_20230724_080000.nc', '/', 'min')

d=1;

for d = 1:length(MCSsample.dirListing)
  %read in the MCS photo count data  
  MCSsample.filename=MCSsample.dirListing(d).name;
 % MCSsample.data = ncread(MCSsample.filename,'Data'); 
  MCSsample.data = h5read(MCSsample.filename,'/Data');  % changed to the h5read
  MCSsample.time = ncread(MCSsample.filename,'time'); 
  MCSsample.channel = ncread(MCSsample.filename,'Channel'); 
  MCSsample.ChannelAssignment = h5read(MCSsample.filename,'/ChannelAssignment');   % note the special type of read
  MCSsample.nsPerBin = ncread(MCSsample.filename,'nsPerBin'); 
  MCSsample.NBins = ncread(MCSsample.filename,'NBins'); 
  MCSsample.ProfilesPerHist = ncread(MCSsample.filename,'ProfilesPerHist'); 
  
  MCS.bins = median(MCSsample.NBins);  % overide the number from the calfiles
  MCS.bin_duration = median(MCSsample.nsPerBin);  % overide the number from the calfiles
  MCS.accum = median(MCSsample.ProfilesPerHist);   % overide the number from the calfiles
  

  if d>1   % sum the days nc data files into a single array
    if size(MCSsample.data,2)==  size(MCSsample.data1,2) % check the MSC didn't change
      MCSsample.time1=[MCSsample.time1;MCSsample.time];
      MCSsample.data1=[MCSsample.data1;MCSsample.data];
      MCSsample.channel1=[MCSsample.channel1;MCSsample.channel];
    end
  else
    MCSsample.time1=MCSsample.time;
    MCSsample.data1=MCSsample.data;
    MCSsample.channel1=MCSsample.channel;
  end
end



for d = 1:length(LL.dirListing)
  %read in the laser locking wavelength data 
  LL.filename=LL.dirListing(d).name;
  LL.time = ncread(LL.filename,'time');
  LL.wavelength = ncread(LL.filename,'Wavelength');
  LL.wavediff = ncread(LL.filename,'WaveDiff');
  LL.current =  ncread(LL.filename,'Current');
  LL.name = h5read(LL.filename,'/LaserName'); % note the special type
  if d>1   % sum the days nc data files into a single array
    LL.time1=[LL.time1;LL.time];
    LL.wavelength1=[LL.wavelength1;LL.wavelength];
    LL.wavediff1=[LL.wavediff1;LL.wavediff]; 
    LL.current1=[LL.current1;LL.current];
    LL.name1 = [LL.name1;LL.name];
  else
    LL.time1 = LL.time;
    LL.wavelength1 = LL.wavelength;
    LL.wavediff1 = LL.wavediff;    
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
    try  % catch addition of thermocouple error this will only allow three thermocouples to be read
      HKeep.time1 = [HKeep.time1; HKeep.time];
      HKeep.temp1 = [HKeep.temp1(:,1:end); HKeep.temp(:,1:end)];
    catch
      HKeep.temp1 = [[HKeep.temp1(:,1:1) HKeep.temp1(:,1:1) HKeep.temp1(:,1:1)]; [HKeep.temp(:,1:1) HKeep.temp(:,1:1) HKeep.temp(:,1:1)]];
    end
  else
      HKeep.time1 = HKeep.time;
      HKeep.temp1 = HKeep.temp;
  end
end

 for d = 1:length(Pow.dirListing)
 %for d = 1:14
  %read in the power montioring data
  Pow.filename = Pow.dirListing(d).name;
 % if exist(Pow.filename) == 2 
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
 % end
 end

for d = 1:length(Etalon.dirListing)
  %read in the power montioring data
  Etalon.filename = Etalon.dirListing(d).name;
 % if exist(Pow.filename) == 2 
    Etalon.time = ncread(Etalon.filename,'time');
    Etalon.temp = ncread(Etalon.filename,'Temperature');
    Etalon.name = h5read(Etalon.filename,'/EtalonNum');
    if d>1   % sum the days nc data files into a single array
      Etalon.time1 = [Etalon.time1; Etalon.time];
      Etalon.temp1 = [Etalon.temp1; Etalon.temp];
      Etalon.name1 = [Etalon.name1; Etalon.name];
    else
      Etalon.time1 = Etalon.time;
      Etalon.temp1 = Etalon.temp;
      Etalon.name1 = Etalon.name;   
    end
 % end
end
 

for d = 1:length(WS.dirListing)
  %read in the weather station data
  WS.filename = WS.dirListing(d).name;
  %ncdisp(WS.filename, '/', 'min') % use this to display all variables
  WS.time = ncread(WS.filename,'time');
  WS.temp = ncread(WS.filename,'Temperature');
  WS.press = ncread(WS.filename,'Pressure');
  WS.relhum = ncread(WS.filename,'RelHum');
  WS.abshum = ncread(WS.filename,'AbsHum');
    if d>1   % sum the days nc data files into a single array
    WS.time1 = [WS.time1; WS.time];
    WS.temp1 = [WS.temp1; WS.temp];
    WS.press1 = [WS.press1; WS.press];
    WS.relhum1 = [WS.relhum1; WS.relhum];
    WS.abshum1 = [WS.abshum1; WS.abshum];
  else
    WS.time1 = WS.time;
    WS.temp1 = WS.temp;
    WS.press1 = WS.press;
    WS.relhum1 = WS.relhum;
    WS.abshum1 = WS.abshum;
  end
end
%figure(10)
%plot(WS.time1, WS.abshum1)
%figure(11)
%plot(WS.time1, WS.temp1)
%figure(12)
%plot(WS.time1, WS.press1)
%figure(13)
%plot(WS.time1, WS.relhum1)
% remove bad data

try
   WS.temp1(WS.temp1<-1000)=NaN;
   WS.press1(WS.press1<-1000)=NaN;
   WS.relhum1(WS.relhum1<-1000)=NaN;  
   WS.abshum1(WS.abshum1<-1000)=NaN; 
%   LL.wavelength1(LL.wavelength1<-1000)=NaN; 
catch
end


% combine time and data 
MCSsample.all = [MCSsample.time1,MCSsample.data1];
% parse out the channels based on searches for channel mapping and populate
% if it exists
index.wvonline = find(contains(MCSsample.ChannelAssignment,'WVOnline') & not(contains(MCSsample.ChannelAssignment,'WVOnlineLow')))-1;
if isempty(index.wvonline) == 0
   MCSsample.online = MCSsample.all(MCSsample.channel1==index.wvonline,:); % 
   MCSsample.online_time = MCSsample.time1(MCSsample.channel1==index.wvonline,:); % 
end
index.wvoffline = find(contains(MCSsample.ChannelAssignment,'WVOffline') & not(contains(MCSsample.ChannelAssignment,'WVOfflineLow')))-1;
if isempty(index.wvoffline) == 0
   MCSsample.offline = MCSsample.all(MCSsample.channel1==index.wvoffline,:); % 
   MCSsample.offline_time = MCSsample.time1(MCSsample.channel1==index.wvoffline,:); % 
end
% index.wvonline_near = find(contains(MCSsample.ChannelAssignment,'WVOnlineLow'))-1;
% if isempty(index.wvonline_near) == 0
%    MCSsample.on_near = MCSsample.all(MCSsample.channel1==index.wvonline_near,:); % 
%    MCSsample.on_near_time = MCSsample.time1(MCSsample.channel1==index.wvonline_near,:); %  
% end
% index.wvoffline_near = find(contains(MCSsample.ChannelAssignment,'WVOfflineLow'))-1;
% if isempty(index.wvoffline_near) == 0
%    MCSsample.off_near = MCSsample.all(MCSsample.channel1==index.wvoffline_near,:); % 
%    MCSsample.off_near_time = MCSsample.time1(MCSsample.channel1==index.wvoffline_near,:); % 
% end

index.O2online_comb = find(contains(MCSsample.ChannelAssignment,'O2OnlineComb') & not(contains(MCSsample.ChannelAssignment,'O2OnlineCombLow')))-1;
if isempty(index.O2online_comb) == 0
   MCSsample.O2online_comb = MCSsample.all(MCSsample.channel1==index.O2online_comb,:); % 
   MCSsample.O2online_comb_time = MCSsample.time1(MCSsample.channel1==index.O2online_comb,:); % 
end

index.O2offline_comb = find(contains(MCSsample.ChannelAssignment,'O2OfflineComb') & not(contains(MCSsample.ChannelAssignment,'O2OfflineCombLow')))-1;
if isempty(index.O2offline_comb) == 0
   MCSsample.O2offline_comb = MCSsample.all(MCSsample.channel1==index.O2offline_comb,:); % 
   MCSsample.O2offline_comb_time = MCSsample.time1(MCSsample.channel1==index.O2offline_comb,:); % 
end

index.O2online_mol = find(contains(MCSsample.ChannelAssignment,'O2OnlineMol')  & not(contains(MCSsample.ChannelAssignment,'O2OnlineMolLow')))-1;
if isempty(index.O2online_mol) == 0
   MCSsample.O2online_mol = MCSsample.all(MCSsample.channel1==index.O2online_mol,:); % 
   MCSsample.O2online_mol_time = MCSsample.time1(MCSsample.channel1==index.O2online_mol,:); % 
end

index.O2offline_mol = find(contains(MCSsample.ChannelAssignment,'O2OfflineMol')  & not(contains(MCSsample.ChannelAssignment,'O2OfflineMolLow')))-1;
if isempty(index.O2offline_mol) == 0
   MCSsample.O2offline_mol = MCSsample.all(MCSsample.channel1==index.O2offline_mol,:); % 
   MCSsample.O2offline_mol_time = MCSsample.time1(MCSsample.channel1==index.O2offline_mol,:); % 
end

index.HSRL_Mol = find(contains(MCSsample.ChannelAssignment,'HSRLMol'))-1;
if isempty(index.HSRL_Mol) == 0
   MCSsample.HSRL_Mol = MCSsample.all(MCSsample.channel1==index.HSRL_Mol,:); % 
   MCSsample.HSRL_Mol_time = MCSsample.time1(MCSsample.channel1==index.HSRL_Mol,:); % 
end

index.HSRL_Combined = find(contains(MCSsample.ChannelAssignment,'HSRLCombined'))-1;
if isempty(index.HSRL_Combined) == 0
   MCSsample.HSRL_Combined = MCSsample.all(MCSsample.channel1==index.HSRL_Combined,:); % 
   MCSsample.HSRL_Combined_time = MCSsample.time1(MCSsample.channel1==index.HSRL_Combined,:); % 
end


if isempty(index.wvonline) == 0

            
            %MCS.time_step = mean(diff(MCSsample.time1))*60*60
            %figure(100)
            %plot((diff(MCSsample.online_time))*60*60)
            %median((diff(MCSsample.online_time))*60*60)
            figure(101)
            plot((diff(MCSsample.offline_time))*60*60)
            offline_timestep = median((diff(MCSsample.offline_time))*60*60)
            online_timestep = median((diff(MCSsample.online_time))*60*60)
            %read_time_in = min(offline_timestep, online_timestep)
            try
              O2offline_timestep = median((diff(MCSsample.O2offline_mol_time))*60*60)
              O2online_timestep = median((diff(MCSsample.O2online_mol_time))*60*60)
            catch
            end
            
            if isfield(LL,'time') == 1
              LL.all = [LL.time1, LL.wavelength1, LL.wavediff1, LL.current1];
              LL.online = LL.all(strcmp(LL.name1,'WVOnline'),:);
              LL.offline = LL.all(strcmp(LL.name1,'WVOffline'),:);
              LL.O2online = LL.all(strcmp(LL.name1,'O2Online'),:);
              LL.O2offline = LL.all(strcmp(LL.name1,'O2Offline'),:);
            end
            
            if isfield(Etalon,'time') == 1
              Etalon.all = [Etalon.time1, Etalon.temp1];
              Etalon.WV_temp = Etalon.all(strcmp(Etalon.name1,'WVEtalon'),:);
              Etalon.O2_temp = Etalon.all(strcmp(Etalon.name1,'O2Etalon'),:);
            end
             figure(11)
             plot(Etalon.WV_temp(:,2))
             hold on
             plot(Etalon.O2_temp(:,2))
             grid on
             
            %LL.hsrl =  LL.all(strcmp(LL.name1,'HSRL'),:);
            if isfield(HKeep,'time') == 1
              try
                HKeep.all = [HKeep.time1, HKeep.temp1];
              catch
              end
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
            ave_time = read_time_in;  % time grid in seconds
            time_grid = (floor(min(MCSsample.time1)):1/60/60*(ave_time):ceil(max(MCSsample.time1)))';
            LL.online_grid = interp1(LL.online(:,1), LL.online(:,2:end), time_grid, 'nearest', 'extrap'); 
            LL.offline_grid = interp1(LL.offline(:,1), LL.offline(:,2:end), time_grid, 'nearest', 'extrap'); 
            LL.O2online_grid = interp1(LL.O2online(:,1), LL.O2online(:,2:end), time_grid, 'nearest', 'extrap'); 
            LL.O2offline_grid = interp1(LL.O2offline(:,1), LL.O2offline(:,2:end), time_grid, 'nearest', 'extrap'); 
            %LL.hsrl_grid = interp1(LL.hsrl(:,1), LL.hsrl(:,2:end), time_grid, 'nearest', 'extrap');
            
            if isempty(index.wvonline) == 0
              [x, ia, ic] = unique(MCSsample.online(:,1),'rows');  % remove duplicate time points
              uA = MCSsample.online(ia,:);  % apply those to the other rows
              MCSsample.online_grid = interp1(uA(:,1), uA(:,2:end), time_grid, 'nearest', 'extrap'); 
            end
            
            if isempty(index.wvoffline) == 0
              [x, ia, ic] = unique(MCSsample.offline(:,1),'rows');  % remove duplicate time points
              uA = MCSsample.offline(ia,:);  % apply those to the other rows
              MCSsample.offline_grid = interp1(uA(:,1), uA(:,2:end), time_grid, 'nearest', 'extrap'); 
            end
            
            % if isempty(index.wvonline_near) == 0
            %   [x, ia, ic] = unique(MCSsample.on_near(:,1),'rows');  % remove duplicate time points
            %   uA = MCSsample.on_near(ia,:);  % apply those to the other rows
            %   MCSsample.on_near_grid = interp1(uA(:,1), uA(:,2:end), time_grid, 'nearest', 'extrap'); 
            % else
            %   MCSsample.on_near_grid = zeros(size(MCSsample.online_grid));
            % end
            % 
            % if isempty(index.wvoffline_near) == 0
            %   [x, ia, ic] = unique(MCSsample.off_near(:,1),'rows');  % remove duplicate time points
            %   uA = MCSsample.off_near(ia,:);  % apply those to the other rows
            %   MCSsample.off_near_grid = interp1(uA(:,1), uA(:,2:end), time_grid, 'nearest', 'extrap'); 
            % else
            %   MCSsample.off_near_grid = zeros(size(MCSsample.offline_grid));
            % end
            
            if isempty(index.O2online_comb) == 0
              [x, ia, ic] = unique(MCSsample.O2online_comb(:,1),'rows');  % remove duplicate time points
              uA = MCSsample.O2online_comb(ia,:);  % apply those to the other rows
              MCSsample.O2online_comb_grid = interp1(uA(:,1), uA(:,2:end), time_grid, 'nearest', 'extrap'); 
            else
              MCSsample.O2online_comb_grid = zeros(size(MCSsample.O2online_comb_grid));
            end
            
            if isempty(index.O2offline_comb) == 0
              [x, ia, ic] = unique(MCSsample.O2offline_comb(:,1),'rows');  % remove duplicate time points
              uA = MCSsample.O2offline_comb(ia,:);  % apply those to the other rows
              MCSsample.O2offline_comb_grid = interp1(uA(:,1), uA(:,2:end), time_grid, 'nearest', 'extrap'); 
            else
              MCSsample.O2offline_comb_grid = zeros(size(MCSsample.O2offline_comb_grid));
            end
            
            if isempty(index.O2online_mol) == 0
              [x, ia, ic] = unique(MCSsample.O2online_mol(:,1),'rows');  % remove duplicate time points
              uA = MCSsample.O2online_mol(ia,:);  % apply those to the other rows
              MCSsample.O2online_mol_grid = interp1(uA(:,1), uA(:,2:end), time_grid, 'nearest', 'extrap'); 
            else
              MCSsample.O2online_mol_grid = zeros(size(MCSsample.O2online_mol_grid));
            end
            
            if isempty(index.O2offline_mol) == 0
              [x, ia, ic] = unique(MCSsample.O2offline_mol(:,1),'rows');  % remove duplicate time points
              uA = MCSsample.O2offline_mol(ia,:);  % apply those to the other rows
              MCSsample.O2offline_mol_grid = interp1(uA(:,1), uA(:,2:end), time_grid, 'nearest', 'extrap'); 
            else
              MCSsample.O2offline_mol_grid = zeros(size(MCSsample.O2offline_mol_grid));
            end
            
            %MCS.combined_grid = interp1(MCS.combined(:,1), MCS.combined(:,2:end), time_grid, 'nearest', 'extrap'); 
            %MCS.molecular_grid = interp1(MCS.molecular(:,1), MCS.molecular(:,2:end), time_grid, 'nearest', 'extrap');
            if isfield(HKeep,'time') == 1
              [x, ia, ic] = unique(HKeep.all(:,1),'rows');  % remove duplicate time points
              uA = HKeep.all(ia,:);  % apply those to the other rows  
              HKeep.temp1 = interp1(uA(:,1), uA(:,2), time_grid, 'nearest');  
              %HKeep.temp1 = interp1(HKeep.all(:,1), HKeep.all(:,2), time_grid, 'nearest');  
              if size(HKeep.all,2)>=3
                HKeep.temp2 = interp1(uA(:,1), uA(:,3), time_grid, 'nearest');
              else
                HKeep.temp2 = interp1(uA(:,1), uA(:,2), time_grid, 'nearest');
              end
            end
            
            [~, ia, ~] = unique(Pow.all(:,1),'rows');  % remove duplicate time points
            uA = Pow.all(ia,:);  % apply those to the other rows
            Pow.online = interp1(uA(:,1), uA(:,2), time_grid, 'nearest', 'extrap');
            Pow.offline = interp1(uA(:,1), uA(:,8), time_grid, 'nearest', 'extrap');
            %Pow.online = HKeep.temp1;
            %Pow.offline = HKeep.temp2;
            if isfield(WS,'time') == 1
              [~, ia, ~] = unique(WS.all(:,1),'rows');  % remove duplicate time points
              uA = WS.all(ia,:);  % apply those to the other rows
              WS.online = interp1(uA(:,1), uA(:,2:end), time_grid, 'nearest', 'extrap');
            end
            
            
            % combine all of the data into original data format
            if isfield(HKeep,'time') == 1 && isfield(WS,'time') == 1
              online_merged = [time_grid, LL.online_grid, Pow.online, HKeep.temp1, WS.online, MCSsample.online_grid];
              offline_merged = [time_grid, LL.offline_grid, Pow.offline,  HKeep.temp2, WS.online, MCSsample.offline_grid];
             % on_near_merged = [time_grid, LL.online_grid, Pow.online, HKeep.temp1, WS.online, MCSsample.on_near_grid];
             % off_near_merged = [time_grid, LL.offline_grid, Pow.offline,  HKeep.temp2, WS.online, MCSsample.off_near_grid];
              O2online_comb_merged = [time_grid, LL.O2online_grid, Pow.online, HKeep.temp1, WS.online, MCSsample.O2online_comb_grid];
              O2offline_comb_merged = [time_grid, LL.O2offline_grid, Pow.offline,  HKeep.temp2, WS.online, MCSsample.O2offline_comb_grid];
              O2online_mol_merged = [time_grid, LL.O2online_grid, Pow.online, HKeep.temp1, WS.online, MCSsample.O2online_mol_grid];
              O2offline_mol_merged = [time_grid, LL.O2offline_grid, Pow.offline,  HKeep.temp2, WS.online, MCSsample.O2offline_mol_grid];
            else
              online_merged = [time_grid, LL.online_grid, Pow.online, MCSsample.online_grid];
              offline_merged = [time_grid, LL.offline_grid, Pow.offline, MCSsample.offline_grid];
              %on_near_merged = [time_grid, LL.online_grid, Pow.online, MCSsample.on_near_grid];
              %off_near_merged = [time_grid, LL.offline_grid, Pow.offline, MCSsample.off_near_grid];
              O2online_comb_merged = [time_grid, LL.online_grid, Pow.online, MCSsample.O2online_comb_grid];
              O2offline_comb_merged = [time_grid, LL.offline_grid, Pow.offline, MCSsample.O2offline_comb_grid];
              O2online_mol_merged = [time_grid, LL.online_grid, Pow.online, MCSsample.O2online_mol_grid];
              O2offline_mol_merged = [time_grid, LL.offline_grid, Pow.offline, MCSsample.O2offline_mol_grid]; 
            end
            %hsrl_com_merged = [time_grid, LL.hsrl_grid, Pow.online, WS.online, MCS.combined_grid];
            %hsrl_mol_merged = [time_grid, LL.hsrl_grid, HKeep.temp1, WS.online, MCS.molecular_grid];

elseif isempty(index.HSRL_Combined) == 0


    if isfield(LL,'time') == 1
      LL.all = [LL.time1, LL.wavelength1, LL.wavediff1, LL.current1];
      LL.hsrl = LL.all(strcmp(LL.name1,'HSRL'),:);
    end
    
    if isfield(Etalon,'time') == 1
      Etalon.all = [Etalon.time1, Etalon.temp1];
      Etalon.hsrl = Etalon.all(strcmp(Etalon.name1,'HSRLEtalon'),:);
    end
     figure(11)
     plot(Etalon.hsrl(:,2))
     grid on

    if isfield(HKeep,'time') == 1
      try
        HKeep.all = [HKeep.time1, HKeep.temp1];
      catch
      end
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
        ave_time = read_time_in;  % time grid in seconds
        time_grid = (floor(min(MCSsample.time1)):1/60/60*(ave_time):ceil(max(MCSsample.time1)))';
        LL.hsrl_grid = interp1(LL.hsrl(:,1), LL.hsrl(:,2:end), time_grid, 'nearest', 'extrap');
        
        if isempty(index.HSRL_Mol) == 0
          [x, ia, ic] = unique(MCSsample.HSRL_Mol(:,1),'rows');  % remove duplicate time points
          uA = MCSsample.HSRL_Mol(ia,:);  % apply those to the other rows
          MCSsample.HSRL_Mol_grid = interp1(uA(:,1), uA(:,2:end), time_grid, 'nearest', 'extrap'); 
        end
        
        if isempty(index.HSRL_Combined) == 0
          [x, ia, ic] = unique(MCSsample.HSRL_Combined(:,1),'rows');  % remove duplicate time points
          uA = MCSsample.HSRL_Combined(ia,:);  % apply those to the other rows
          MCSsample.HSRL_Combined_grid = interp1(uA(:,1), uA(:,2:end), time_grid, 'nearest', 'extrap'); 
        end

        if isfield(HKeep,'time') == 1
          [x, ia, ic] = unique(HKeep.all(:,1),'rows');  % remove duplicate time points
          uA = HKeep.all(ia,:);  % apply those to the other rows  
          HKeep.temp1 = interp1(uA(:,1), uA(:,2), time_grid, 'nearest');  
          %HKeep.temp1 = interp1(HKeep.all(:,1), HKeep.all(:,2), time_grid, 'nearest');  
          if size(HKeep.all,2)>=3
            HKeep.temp2 = interp1(uA(:,1), uA(:,3), time_grid, 'nearest');
          else
            HKeep.temp2 = interp1(uA(:,1), uA(:,2), time_grid, 'nearest');
          end
        end

        [~, ia, ~] = unique(Pow.all(:,1),'rows');  % remove duplicate time points
            uA = Pow.all(ia,:);  % apply those to the other rows
            Pow.online = interp1(uA(:,1), uA(:,2), time_grid, 'nearest', 'extrap');
            Pow.offline = interp1(uA(:,1), uA(:,8), time_grid, 'nearest', 'extrap');

        if isfield(WS,'time') == 1
             [~, ia, ~] = unique(WS.all(:,1),'rows');  % remove duplicate time points
             uA = WS.all(ia,:);  % apply those to the other rows
            WS.online = interp1(uA(:,1), uA(:,2:end), time_grid, 'nearest', 'extrap');
        end
            
            
        % combine all of the data into original data format
       if isfield(HKeep,'time') == 1 && isfield(WS,'time') == 1
           HSRL_Combined_merged  = [time_grid, LL.hsrl_grid, Pow.online, HKeep.temp1, WS.online, MCSsample.HSRL_Combined_grid];
           HSRL_Mol_merged = [time_grid, LL.hsrl_grid, Pow.offline,  HKeep.temp2, WS.online, MCSsample.HSRL_Mol_grid];
       else
           HSRL_Combined_merged = [time_grid, LL.hsrl_grid, Pow.online, MCSsample.HSRL_Combined_grid];
           HSRL_Mol_merged = [time_grid, LL.hsrl_grid, Pow.offline, MCSsample.HSRL_Mol_grid];
       end

end


cd(dd) % point back to original directory

end

