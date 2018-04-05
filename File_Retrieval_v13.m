function[data,data1,folder] = File_Retrieval_v13(bins, folder)

%modification....16-May 2012 by Spuler to run on Mac (forward slashes and 
% need to check on first folder because the file structure from uigetdir 
% has 1=., 2=.., and 3 = .DS_Store, and the subfolders begin in some directories)

dd = pwd; % get the current path
%cd('/scr/rsf1/vanandel/HSRL/h2o_data/2012/') % point to the directory where data is stored  
%cd('/Users/spuler/Dropbox/WV_DIAL/data/') % point to the directory where data is stored  
%cd('/Users/spuler/Desktop/WV_DIAL_data/data/2014/') % point to the directory where data is stored 
%cd('/Volumes/documents/WV_DIAL_data/data/2014/') % point to the directory where data is stored 
%cd('/scr/eldora1/wvdial_1_data/2014/') % point to the directory where data is stored 

amin=0;
data=0;
% first, you have to find the folder
%folder = uigetdir; % check the help for uigetdir to see how to specify a starting path, which makes your life easier

%remove any hidden folders 
cd(folder)
delete('._*') 

% get the names of all files. dirListing is a struct array. 
dirListing = dir(strcat(folder));
% dirListing=dir(strcat(folder,'\'));
% for jj=1:length(

%check for first folder 
test = dirListing(3).name; %check if the 3rd file (after . and ..) in the structure
if  strcmp(test,'.DS_Store')
    f=4;
else
    f=3;
end
d=f; 

cd(dd) % point back to original directory

% reallocate memory for the data based on file size?
%data = zeros(58010,1994);
%data1 = zeros(58010,1995);

% loop through the files and open. Note that dir also lists the directories, so you have to check for them.
%h=waitbar(0,'Online Retrieval');
for d = f:length(dirListing);
%    waitbar(d/(length(dirListing)),h);
% if ~dirListing(d).isdir
    dirListing_new = dir(strcat(folder,'/',dirListing(d).name,'/','Online_Raw_Data.dat'));
   
  
    if ~dirListing_new(1).isdir
      fileName = fullfile(folder,dirListing(d).name,dirListing_new(1).name); % use full path because the folder may not be the active path
      % open your file here 
      fid = fopen(fileName);
      amin=fread(fid, [(7+bins-1),inf], 'double', 'b');
      fclose(fid);   
      if d>f
        data=[data;amin'];
      else
        data=amin';
      end
    % do something
    end % if-clause
end % for-loop
%delete(h)



% loop through the files and open. Note that dir also lists the directories, so you have to check for them.
amin=0;
%h=waitbar(0,'Offline Retrieval');
for d = f:length(dirListing);
%    waitbar(d/(length(dirListing)),h);
% if ~dirListing(d).isdir
    dirListing_new = dir(strcat(folder,'/',dirListing(d).name,'/','Offline_Raw_Data.dat'));

 
    if ~dirListing_new(1).isdir
      fileName = fullfile(folder,dirListing(d).name,dirListing_new(1).name); % use full path because the folder may not be the active path
      % open your file here 
      fid = fopen(fileName);
      amin=fread(fid, [(7+bins-1),inf], 'double', 'b');
      fclose(fid);   
      if d>f
        data1=[data1;amin'];
      else
        data1=amin';
      end
    % do something
    end % if-clause
end % for-loop
%delete(h)
