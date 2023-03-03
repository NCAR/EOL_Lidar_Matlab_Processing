


function [Cell] = LoaderMPD(Systems,Dates,Variables,PathGen,FileGen,DesiredTime,DesiredRange)
%
%
%
%
%
%
%% Pre-allocating data
Cell = cell(length(Systems),1);
%% Looping over systems
for m=1:1:length(Systems)
    % Pre-allocating cell elements for each variable
    Cell{m,1} = cell(length(Variables),1);
    % Looping over days
    for n=Dates(1):1:Dates(2)
        % Checking what day to load
        DateString = datestr(n,'yyyymmdd');
        fprintf(['Processing: MPD',Systems{m},' ',DateString,'\n'])
        % Determining filename
        FilePath = sprintf(PathGen,Systems{m});
        FileName = sprintf(FileGen,Systems{m},DateString);
        File     = fullfile(FilePath,FileName);
        if isfile(File)
            % Looping over variables
            for p=1:1:length(Variables)
                % Vars is cell with rows (time, range, value)
                for q=1:1:length(Variables{p})
                    Vars{q,1} = LoadVariable(File,Variables{p}{q});
                end
                % Interpolating 
                ValInterp = InterpolateGeneral(Vars{1},Vars{2},Vars{3},DesiredTime,DesiredRange);
                % Saving data
                Cell{m,1}{p,1} = [Cell{m,1}{p,1},ValInterp];
            end
        else
            for p=1:1:length(Variables)
                % Creating a nan padded array and saving
                ValInterp = meshgrid(DesiredTime,DesiredRange).*nan;
                Cell{m,1}{p,1} = [Cell{m,1}{p,1},ValInterp];
            end
        end
    end
end
end

function [Var] = LoadVariable(FileName,Varname)
%
%
%
%
%
%%
try
    if contains(FileName,'.nc')
        Var = ncread(FileName,Varname);
    elseif contains(FileName,'.mat')
        % Checking if element is part of a structure 
        if contains(Varname,'.')
            % Determining how many different levels to descend in struct
            A = split(Varname,'.');
            % Loading top level struct
            S = load(FileName,A{1});
            % Loading variable by pulling out the variable from the struct
            Var = getfield(S,A{:});
        else
            Var = load(FileName,Varname);
        end
    else
        error('Filetype not recognized')
    end
catch
    Var = nan;
end
end

function [ValInterp] = InterpolateGeneral(Time,Range,Value,DesiredTime,DesiredRange)
%
%
%
%
%
%% Checking if data is availible
if all(isnan(Time)) || all(isnan(Range)) || all(all(isnan(Value)))
   ValInterp = meshgrid(DesiredTime,DesiredRange).*nan;
   return
end

%% If availible, check the dimension and interpolate
if isempty(Range) 
    ValInterp = interp1(Time,Value,DesiredTime);
else
    [X,Y] = meshgrid(Time,Range);
    [x,y] = meshgrid(DesiredTime,DesiredRange);
    ValInterp = interp2(X,Y,double(Value),x,y);
end
end
