function []=dta2txt(Path,Option)
%% Converts Bruker ESR files to txt files
% Scans given folder and subfolders for Bruker ESR files and converts both
% DSC and DTA files to .txt files containing data and experimantal
% conditions. Resulting files can be saved either in the same folder with
% the originals or in a new "Converted" folder.
% If no folder path is given, DTA_to_txt promts the user to enter a path
% through a GUI.
% Option specifies whether the converted files should be saved in the
% 'Original folder', 'New folder' or a 'New folder structure'.

%% Read folder

% "Path" is the folder containing the ESR Data (either directly or in
% sub-folders)

if nargin<1
    Path='/Users/SamSchott/Documents/Studium/Cambridge/PhD Project/Nagoya/ESR Data/FI_ESR_1_lowT';
end

Path = uigetdir(Path);
 
% List all Bruker ESR files in "Path" and sub-folders

sub_dir = regexp(genpath(Path),'[^:]*','match');

for i=1:length(sub_dir)
    Content=dir([sub_dir{i},'/*.DSC']);
    for j=1:length(Content)
        file_name{j,i}=Content(j).name;
    end
end

%% Convert to text files
if nargin<2
    Option = questdlg('Where would you like to save converted files?', ...
	'Location','Orginal location','New folder structure','New folder','New folder');
end

if strcmp(Option,'New folder')==1
mkdir([Path,'/Converted']);
[Nfiles,Ndir]=size(file_name);

for i=1:Ndir
    for j=1:Nfiles
        if isempty(file_name{j,i})==0
            copyfile([sub_dir{i},'/',file_name{j,i}],[Path,'/Converted/PARAM_',strrep(file_name{j,i},'DSC','txt')]);
        end
    end
end

for i=1:Ndir
    for j=1:Nfiles
        if isempty(file_name{j,i})==0
           clear Int
           clear B
           [B,Int]=BrukerRead([sub_dir{i},'/',strrep(file_name{j,i},'DSC','DTA')]);
           data = [B Int];
           save([Path,'/Converted/DATA_',strrep(file_name{j,i},'DSC','txt')],'data','-ascii');
        end
    end
end
end

% write to folder(s)

if strcmp(Option,'New folder structure')==1
    [Nfiles,Ndir]=size(file_name);
    new_path=strrep(sub_dir,Path,'');

for i=1:Ndir
    for j=1:Nfiles
        if isempty(file_name{j,i})==0
            mkdir([Path,'/Converted',new_path{i}]);
            copyfile([sub_dir{i},'/',file_name{j,i}],[Path,'/Converted',new_path{i},'/PARAM_',strrep(file_name{j,i},'DSC','txt')]);
        end
    end
end

for i=1:Ndir
    for j=1:Nfiles
        if isempty(file_name{j,i})==0
           [B,Int]=BrukerRead([sub_dir{i},'/',strrep(file_name{j,i},'DSC','DTA')]);
           data = [B Int];
           save([Path,'/Converted',new_path{i},'/DATA_',strrep(file_name{j,i},'DSC','txt')],'data','-ascii');
        end
    end
end
end

if strcmp(Option,'Original folder')==1

[Nfiles,Ndir]=size(file_name);

for i=1:Ndir
    for j=1:Nfiles
        if isempty(file_name{j,i})==0
            copyfile([sub_dir{i},'/',file_name{j,i}],[sub_dir{i},'/PARAM_',strrep(file_name{j,i},'DSC','txt')]);
        end
    end
end

for i=1:Ndir
    for j=1:Nfiles
        if isempty(file_name{j,i})==0
           [B,Int]=BrukerRead([sub_dir{i},'/',strrep(file_name{j,i},'DSC','DTA')]);
           data = [B Int];
           save([sub_dir{i},'/DATA_',strrep(file_name{j,i},'DSC','txt')],'data','-ascii');
        end
    end
end
end
