%% MPD_read_and_grid_data.m: Unified Data Retrieval and Gridding Service

function [data_struct, MCS] = MPD_read_and_grid_data(node, folder_in, channel, config, MCS_in)

% Output structure initialization
data_struct = struct();
MCS = MCS_in;
dd = pwd; 
cd(folder_in);

read_time_in = config.read_time_in; % seconds


try
    % --- 1. Consolidated File Listing and Reading ---
    
    % This helper function now contains the core logic from your original V6 file.
    [MCSsample, LL, HKeep, Pow, WS, Etalon] = read_and_combine_files_v6_logic(folder_in, channel);

    % Update MCS properties from data files
    if isfield(MCSsample, 'NBins') && ~isempty(MCSsample.NBins)
        MCS.bins = median(MCSsample.NBins); 
        MCS.bin_duration = median(MCSsample.nsPerBin); 
        MCS.accum = median(MCSsample.ProfilesPerHist);
    end
    
    % Handle empty data case
    if isempty(MCSsample.time1)
        error('No data found in MCS files.');
    end

    % --- 2. Channel Mapping and Separation ---
    
    % CRITICAL FIX: Robustly define MCSsample.all here.
    if isempty(MCSsample.time1) || isempty(MCSsample.data1)
        error('MCS raw data arrays (time1 or data1) are empty after file reading.');
    end
    MCSsample.all = [MCSsample.time1, MCSsample.data1]; 
    
    % Prepare containers map to align channel data with outputs
    channel_map = containers.Map;
    channel_map('WVOnline') = 'online'; 
    channel_map('WVOffline') = 'offline';
    channel_map('O2OnlineComb') = 'O2online_comb'; 
    channel_map('O2OfflineComb') = 'O2offline_comb';
    channel_map('O2OnlineMol') = 'O2online_mol'; 
    channel_map('O2OfflineMol') = 'O2offline_mol';
    channel_map('HSRLMol') = 'HSRL_Mol'; 
    channel_map('HSRLCombined') = 'HSRL_Combined';
    
    channel_names_in_data = keys(channel_map);
    
    % Extract data based on channel assignments
    for k = 1:length(channel_names_in_data)
        ch_key = channel_names_in_data{k};
        ch_var = channel_map(ch_key);
        
        index = find(strcmp(MCSsample.ChannelAssignment, ch_key)) - 1; 

        if ~isempty(index)
            raw_data = MCSsample.all(MCSsample.channel1 == index, :);
            MCSsample.(ch_var).raw = raw_data;
        else
            MCSsample.(ch_var).raw = []; 
        end
    end
    
    % ----------------------------------------------------------------------------------
    % --- CRITICAL FIX: Enforce Unique Timestamps (Data Cleaning before Interpolation) ---
    
    % 1. Clean the master MCS time array
    [MCSsample.time1, idx_mcs] = unique(MCSsample.time1, 'stable');
    MCSsample.data1 = MCSsample.data1(idx_mcs, :);
    MCSsample.channel1 = MCSsample.channel1(idx_mcs, :);
    
    % 2. Clean the auxiliary data time arrays for stable interpolation
    [LL.time1, idx_ll] = unique(LL.time1, 'stable');
    LL.wavelength1 = LL.wavelength1(idx_ll, :);
    LL.wavediff1 = LL.wavediff1(idx_ll, :);
    LL.current1 = LL.current1(idx_ll, :);
    
    [HKeep.time1, idx_hk] = unique(HKeep.time1, 'stable');
    HKeep.temp1 = HKeep.temp1(idx_hk, :);
    
    [Pow.time1, idx_pow] = unique(Pow.time1, 'stable');
    Pow.power1 = Pow.power1(idx_pow, :);
    
    [WS.time1, idx_ws] = unique(WS.time1, 'stable');
    WS.temp1 = WS.temp1(idx_ws, :);
    WS.press1 = WS.press1(idx_ws, :);
    WS.relhum1 = WS.relhum1(idx_ws, :);
    WS.abshum1 = WS.abshum1(idx_ws, :);
    
    % ----------------------------------------------------------------------------------


    % --- 3. Time Gridding & Interpolation ---

    time_grid = (floor(min(MCSsample.time1)):1/60/60*(read_time_in):ceil(max(MCSsample.time1)))';
    data_struct.time_grid = time_grid;
    N_grid = length(time_grid);
    N_time_shift = 3; % Number of LL fields (Wvl, WvDiff, Curr)
    N_temp_shift = 1; % Number of Temp fields (T1)

    % Interpolate HK and LL data (Nearest/Extrap)
    
    % --- 1. WV LL Interpolation (using LL.online) ---
    if ~isempty(LL.online) && size(LL.online,1) >= 2
        data_struct.LL_grid.online = interp1(LL.online(:,1), LL.online(:,2), time_grid, 'nearest', 'extrap'); % Wvl
        data_struct.LL_grid.online(:,2) = interp1(LL.online(:,1), LL.online(:,3), time_grid, 'nearest', 'extrap'); % WvlDiff
        data_struct.LL_grid.online(:,3) = interp1(LL.online(:,1), LL.online(:,4), time_grid, 'nearest', 'extrap'); % Curr
    else
        data_struct.LL_grid.online = zeros(N_grid, N_time_shift);
    end
    data_struct.LL_grid.offline = data_struct.LL_grid.online; % Default offline to online for WV
    
    % --- 2. O2 LL Interpolation (using LL.O2online) ---
    if ~isempty(LL.O2online) && size(LL.O2online,1) >= 2
        data_struct.LL_grid.O2online = interp1(LL.O2online(:,1), LL.O2online(:,2), time_grid, 'nearest', 'extrap'); % Wvl
        data_struct.LL_grid.O2online(:,2) = interp1(LL.O2online(:,1), LL.O2online(:,3), time_grid, 'nearest', 'extrap'); % WvlDiff
        data_struct.LL_grid.O2online(:,3) = interp1(LL.O2online(:,1), LL.O2online(:,4), time_grid, 'nearest', 'extrap'); % Curr
    else
        % Fallback for O2 if no O2 LL data is available
        data_struct.LL_grid.O2online = zeros(N_grid, N_time_shift); 
    end
    data_struct.LL_grid.O2offline = data_struct.LL_grid.O2online; % Default O2 offline to O2 online
    
    % **HKeep Interpolation (Robust Check)**
    if length(HKeep.time1) >= 2
        data_struct.HKeep_grid.temp1 = interp1(HKeep.time1, HKeep.temp1(:,1), time_grid, 'nearest', 'extrap');
        data_struct.HKeep_grid.temp2 = interp1(HKeep.time1, HKeep.temp1(:,2), time_grid, 'nearest', 'extrap');
    else
        warning('HK data incomplete (less than 2 unique points). Using placeholder arrays.');
        data_struct.HKeep_grid.temp1 = zeros(N_grid, N_temp_shift);
        data_struct.HKeep_grid.temp2 = zeros(N_grid, N_temp_shift);
    end

    % **Pow Interpolation (Robust Check)**
    if length(Pow.time1) >= 2
        data_struct.Pow_grid.online = interp1(Pow.time1, Pow.power1(:,1), time_grid, 'nearest', 'extrap');
        data_struct.Pow_grid.offline = interp1(Pow.time1, Pow.power1(:,8), time_grid, 'nearest', 'extrap');
    else
        data_struct.Pow_grid.online = zeros(N_grid, 1);
        data_struct.Pow_grid.offline = zeros(N_grid, 1);
        warning('Power data incomplete (less than 2 unique points). Using placeholder arrays.');
    end

    
    % --- CRITICAL FIX: Robust WS Data Interpolation and Stability Check ---
    
    N_ws_cols = 4; % [Temp, Press, RelHum, AbsHum]
    WS_stable_array = NaN(N_grid, N_ws_cols);
    
    % Interpolate and fill with defaults if individual sensor data is empty/missing
    if length(WS.time1) >= 2
        WS_stable_array(:,1) = interp1(WS.time1, WS.temp1, time_grid, 'nearest', 'extrap');
        WS_stable_array(:,2) = interp1(WS.time1, WS.press1, time_grid, 'nearest', 'extrap');
        WS_stable_array(:,3) = interp1(WS.time1, WS.relhum1, time_grid, 'nearest', 'extrap');
        WS_stable_array(:,4) = interp1(WS.time1, WS.abshum1, time_grid, 'nearest', 'extrap');
    else
        warning('WS data incomplete (less than 2 unique points). Using placeholder arrays.');
    end
    
    % Final stability check on WS data (Setting default values if all NaN)
    if all(isnan(WS_stable_array(:,1)))
        warning('WS Surface Temperature array is all NaN. Setting to default 290 K.');
        WS_stable_array(:,1) = 290 - 273.15; % Default (Celsius)
    end
    if all(isnan(WS_stable_array(:,2)))
        warning('WS Surface Pressure array is all NaN. Setting to default 1000 hPa.');
        WS_stable_array(:,2) = 1000;
    end

    % Final assignment for WS grid (guarantees a stable 4-column array)
    data_struct.WS_grid.online = WS_stable_array;


    % --- 4. Photon Count Gridding (Using grid_photon_counts_2D_v3) ---
    
    channel_list = keys(channel_map);
    
    for k = 1:length(channel_list)
        ch_key = channel_list{k};
        ch_var = channel_map(ch_key);
        raw_data_struct = MCSsample.(ch_var);
        
        if isfield(raw_data_struct, 'raw') && ~isempty(raw_data_struct.raw)
            raw_data = raw_data_struct.raw;
            
            % Columns 1 is time, Cols 11:end are counts (assuming V6 header structure)
            time_data_cols = [raw_data(:,1), raw_data(:,11:end)];
            
            [counts_grid, ~] = grid_photon_counts_2D_v3(time_data_cols, time_grid, read_time_in);
            data_struct.(strcat(ch_var, '_grid')) = counts_grid;

            
        else
            % Pad with zeros matching the grid size/bin count
            data_struct.(strcat(ch_var, '_grid')) = zeros(N_grid, MCS.bins); 
        end
    end

    
    % --- 5. Final Merging (Replicating the array structure for compatibility) ---
    
    % Reconstruct the assumed 10-column header for compatibility with analysis functions
    % [1:time] [2-4:LL] [5:Pow] [6:HKeep_T1] [7:WS_T] [8:WS_P] [9:WS_RH] [10:WS_AH]
    
   % header_10_cols = [data_struct.time_grid, data_struct.LL_grid.online, ...
   %                   data_struct.Pow_grid.online, data_struct.HKeep_grid.temp1, ...
   %                   data_struct.WS_grid.online]; 
   % 
   % data_struct.online_merged = [header_10_cols, data_struct.online_grid];
   % data_struct.offline_merged = [header_10_cols, data_struct.offline_grid];
   % data_struct.O2online_comb_merged = [header_10_cols, data_struct.O2online_comb_grid];
   % data_struct.O2offline_comb_merged = [header_10_cols, data_struct.O2offline_comb_grid];
   % data_struct.O2online_mol_merged = [header_10_cols, data_struct.O2online_mol_grid];
   % data_struct.O2offline_mol_merged = [header_10_cols, data_struct.O2offline_mol_grid];

   % WV Header (Uses LL_grid.online - which is the WV wavelength)
    wv_header_cols = [data_struct.time_grid, data_struct.LL_grid.online, ...
                      data_struct.Pow_grid.online, data_struct.HKeep_grid.temp1, ...
                      data_struct.WS_grid.online]; 
    
    % O2 Header (Uses LL_grid.O2online - which is the O2 wavelength)
    O2_header_cols = [data_struct.time_grid, data_struct.LL_grid.O2online, ...
                      data_struct.Pow_grid.online, data_struct.HKeep_grid.temp1, ...
                      data_struct.WS_grid.online];
                      
    % Assign WV merged arrays
    data_struct.online_merged = [wv_header_cols, data_struct.online_grid];
    data_struct.offline_merged = [wv_header_cols, data_struct.offline_grid];
    
    % Assign O2 merged arrays
    data_struct.O2online_comb_merged = [O2_header_cols, data_struct.O2online_comb_grid];
    data_struct.O2offline_comb_merged = [O2_header_cols, data_struct.O2offline_comb_grid];
    data_struct.O2online_mol_merged = [O2_header_cols, data_struct.O2online_mol_grid];
    data_struct.O2offline_mol_merged = [O2_header_cols, data_struct.O2offline_mol_grid];
    
    % The HSRL arrays may need an HSRL-specific header, but for now, use the O2 one:
    % (assuming HSRL uses the O2 laser/metadata, which is standard)
    data_struct.HSRL_Mol_merged = [O2_header_cols, data_struct.HSRL_Mol_grid];
    data_struct.HSRL_Combined_merged = [O2_header_cols, data_struct.HSRL_Combined_grid];

    
catch ME
    warning(['Error reading or processing files in ', folder_in, ': ', ME.message]);
    cd(dd);
    data_struct = []; 
    return;
end

cd(dd); 

end


% --- Internal Helper Function: CONTAINS INTEGRATED V6 FILE READING LOGIC ---
function [MCSsample, LL, HKeep, Pow, WS, Etalon] = read_and_combine_files_v6_logic(folder_in, channel)
    
    % **THIS BLOCK CONTAINS THE INTEGRATED V6 LOGIC.**
    
    dd = pwd; % Store current path before changing directory
    cd(folder_in)
    
    % Initialize structures and search file types
    MCSsample.dirListing = dir(strcat('MCS','*')); 
    LL.dirListing = dir(strcat('LL','*')); 
    HKeep.dirListing = dir(strcat('HKeep','*')); 
    Pow.dirListing = dir(strcat('Pow','*')); 
    WS.dirListing = dir(strcat('WS','*')); 
    Etalon.dirListing = dir(strcat('Etalon','*')); 
    
    % Initialize all data storage arrays (as done in V6)
    MCSsample.time1 = []; MCSsample.data1 = []; MCSsample.channel1 = []; MCSsample.ChannelAssignment = {};
    LL.time1 = []; LL.wavelength1 = []; LL.wavediff1 = []; LL.current1 = []; LL.name1 = {};
    HKeep.time1 = []; HKeep.temp1 = [];
    Pow.time1 = []; Pow.power1 = []; Pow.channel = {};
    WS.time1 = []; WS.temp1 = []; WS.press1 = []; WS.relhum1 = []; WS.abshum1 = [];
    Etalon.time1 = []; Etalon.temp1 = []; Etalon.name1 = {};

    
    % 1. Read and combine MCS data
    for d = 1:length(MCSsample.dirListing)
        MCSsample.filename=MCSsample.dirListing(d).name;
        MCSsample.data = double(h5read(MCSsample.filename,'/Data'));  % uses h5read
        MCSsample.time = ncread(MCSsample.filename,'time'); 
        MCSsample.channel = ncread(MCSsample.filename,'Channel'); 
        MCSsample.ChannelAssignment = h5read(MCSsample.filename,'/ChannelAssignment');
        MCSsample.nsPerBin = double(ncread(MCSsample.filename,'nsPerBin')); 
        MCSsample.NBins = double(ncread(MCSsample.filename,'NBins')); 
        MCSsample.ProfilesPerHist = double(ncread(MCSsample.filename,'ProfilesPerHist')); 
        
        if d>1   
            if size(MCSsample.data,2)== size(MCSsample.data1,2) 
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
    
    % 2. Read and combine LL data
    for d = 1:length(LL.dirListing)
        LL.filename=LL.dirListing(d).name;
        LL.time = ncread(LL.filename,'time');
        LL.wavelength = ncread(LL.filename,'Wavelength');
        LL.wavediff = ncread(LL.filename,'WaveDiff');
        LL.current =  ncread(LL.filename,'Current');
        LL.name = h5read(LL.filename,'/LaserName'); 
        if d>1   
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

    % 3. Read and combine HKeep data
    for d = 1:length(HKeep.dirListing)
        HKeep.filename = HKeep.dirListing(d).name;
        HKeep.time = ncread(HKeep.filename,'time');
        HKeep.temp = ncread(HKeep.filename,'Temperature');
        if d>1   
            try
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

    % 4. Read and combine Pow data
    for d = 1:length(Pow.dirListing)
        Pow.filename = Pow.dirListing(d).name;
        Pow.time = ncread(Pow.filename,'time');
        Pow.power = double(ncread(Pow.filename,'Power'));  % Changed this line for MCX int32
        Pow.channel = h5read(Pow.filename,'/ChannelAssignment');
        if d>1   
            Pow.time1 = [Pow.time1; Pow.time];
            Pow.power1 = [Pow.power1; Pow.power];
        else
            Pow.time1 = Pow.time;
            Pow.power1 = Pow.power;
        end
    end
    
    % 5. Read and combine WS data
    for d = 1:length(WS.dirListing)
        WS.filename = WS.dirListing(d).name;
        WS.time = ncread(WS.filename,'time');
        WS.temp = ncread(WS.filename,'Temperature');
        WS.press = ncread(WS.filename,'Pressure');
        WS.relhum = ncread(WS.filename,'RelHum');
        WS.abshum = ncread(WS.filename,'AbsHum');
        if d>1   
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
    
    % 6. Read and combine Etalon data 
    for d = 1:length(Etalon.dirListing)
        Etalon.filename = Etalon.dirListing(d).name;
        Etalon.time = ncread(Etalon.filename,'time');
        Etalon.temp = ncread(Etalon.filename,'Temperature');
        Etalon.name = h5read(Etalon.filename,'/EtalonNum');
        if d>1   
            Etalon.time1 = [Etalon.time1; Etalon.time];
            Etalon.temp1 = [Etalon.temp1; Etalon.temp];
            Etalon.name1 = [Etalon.name1; Etalon.name];
        else
            Etalon.time1 = Etalon.time;
            Etalon.temp1 = Etalon.temp;
            Etalon.name1 = Etalon.name;   
        end
    end
    
    % Handle bad data values (as in V6)
    try
        WS.temp1(WS.temp1<-1000)=NaN;
        WS.press1(WS.press1<-1000)=NaN;
        WS.relhum1(WS.relhum1<-1000)=NaN;  
        WS.abshum1(WS.abshum1<-1000)=NaN; 
    catch
    end
    
    % 7. Separate LL streams (as in V6)
    if isfield(LL,'name1') && ~isempty(LL.name1)
        LL.all = [LL.time1, LL.wavelength1, LL.wavediff1, LL.current1];
        
        % Check if LL.all is empty before attempting indexing
        if isempty(LL.all)
            LL.online = []; LL.offline = []; LL.O2online = []; LL.O2offline = [];
        else
            LL.online = LL.all(strcmp(LL.name1,'WVOnline'),:);
            LL.offline = LL.all(strcmp(LL.name1,'WVOffline'),:);
            LL.O2online = LL.all(strcmp(LL.name1,'O2Online'),:);
            LL.O2offline = LL.all(strcmp(LL.name1,'O2Offline'),:);
        end
        
    else
        % Provide empty placeholders if no LL data was found
        LL.online = []; LL.offline = []; LL.O2online = []; LL.O2offline = [];
    end
    
    % 8. Separate MCS channels into raw data streams (as in V6)
    
    if ~isempty(MCSsample.time1)
        MCSsample.all = [MCSsample.time1, MCSsample.data1]; % Redefined here for robustness
    else
        MCSsample.all = [];
    end
    
    % Extract channels only if MCSsample.all is not empty
    if ~isempty(MCSsample.all)
        index.wvonline = find(strcmp(MCSsample.ChannelAssignment, 'WVOnline')) - 1;
        if isempty(index.wvonline) == 0
           MCSsample.online.raw = MCSsample.all(MCSsample.channel1==index.wvonline,:);
        else; MCSsample.online.raw = []; end

        index.wvoffline = find(strcmp(MCSsample.ChannelAssignment, 'WVOffline')) - 1;
        if isempty(index.wvoffline) == 0
           MCSsample.offline.raw = MCSsample.all(MCSsample.channel1==index.wvoffline,:);
        else; MCSsample.offline.raw = []; end

        index.O2online_comb = find(strcmp(MCSsample.ChannelAssignment, 'O2OnlineComb')) - 1;
        if isempty(index.O2online_comb) == 0
           MCSsample.O2online_comb.raw = MCSsample.all(MCSsample.channel1==index.O2online_comb,:);
        else; MCSsample.O2online_comb.raw = []; end

        index.O2offline_comb = find(strcmp(MCSsample.ChannelAssignment, 'O2OfflineComb')) - 1;
        if isempty(index.O2offline_comb) == 0
           MCSsample.O2offline_comb.raw = MCSsample.all(MCSsample.channel1==index.O2offline_comb,:);
        else; MCSsample.O2offline_comb.raw = []; end

        index.O2online_mol = find(strcmp(MCSsample.ChannelAssignment, 'O2OnlineMol')) - 1;
        if isempty(index.O2online_mol) == 0
           MCSsample.O2online_mol.raw = MCSsample.all(MCSsample.channel1==index.O2online_mol,:);
        else; MCSsample.O2online_mol.raw = []; end

        index.O2offline_mol = find(strcmp(MCSsample.ChannelAssignment, 'O2OfflineMol')) - 1;
        if isempty(index.O2offline_mol) == 0
           MCSsample.O2offline_mol.raw = MCSsample.all(MCSsample.channel1==index.O2offline_mol,:);
        else; MCSsample.O2offline_mol.raw = []; end

        index.HSRL_Mol = find(strcmp(MCSsample.ChannelAssignment, 'HSRLMol')) - 1;
        if isempty(index.HSRL_Mol) == 0
           MCSsample.HSRL_Mol.raw = MCSsample.all(MCSsample.channel1==index.HSRL_Mol,:);
        else; MCSsample.HSRL_Mol.raw = []; end

        index.HSRL_Combined = find(strcmp(MCSsample.ChannelAssignment, 'HSRLCombined')) - 1;
        if isempty(index.HSRL_Combined) == 0
           MCSsample.HSRL_Combined.raw = MCSsample.all(MCSsample.channel1==index.HSRL_Combined,:);
        else; MCSsample.HSRL_Combined.raw = []; end


    else
        % If MCS data is truly empty, initialize channel raw fields to empty
        MCSsample.online.raw = []; MCSsample.offline.raw = [];
        MCSsample.O2online_comb.raw = []; MCSsample.O2offline_comb.raw = [];
        MCSsample.O2online_mol.raw = []; MCSsample.O2offline_mol.raw = [];
        MCSsample.HSRL_Mol.raw = []; MCSsample.HSRL_Combined.raw = [];
    end


    cd(dd) % Return to original directory
    
end