function [output_table] = BatchLoad(slice)
% BATCHLOAD Loads and returns data from multiple ESR data files
%
%   INPUTS:
%   slice   - for 2D data-sets only: return a slice of the data set, e.g.,
%             at a certain microwave power or goniometer angle
%
%   OUTPU:
%   output_table    - table with all data
%

import esr_analyses.*
import esr_analyses.utils.*

if nargin < 2
    slice = false;
end

% try automatical matching of background signals
bckgrndStrQ = input('Would you like to subtract background spectra? [y]/n ', 's');

global Path

[fileNames, Path] = uigetfile([Path, '*.DSC'], ...
    'Please select ESR data files.', 'MultiSelect', 'on');

    
if ~iscell(fileNames)
    fileNames = {fileNames};
end
    
nFiles = numel(fileNames);
output_table = table;

for i = 1:nFiles
    
    filePath = [Path fileNames{i}];    
    filePathBG = [];
    
    if ~strcmp(bckgrndStrQ, 'n')
        filePathBGcandidates = {'Vg_00.', 'Vg_0.'};
        
        for candidate=filePathBGcandidates
            filePathBGcandidate = regexprep(filePath, 'Vg_(\S+)(\.)', candidate{1});
            if exist(filePathBGcandidate, 'file') == 2
                filePathBG = filePathBGcandidate;
                break
            end
        end
  
        if filePathBG
            disp('Background file:');
            disp(filePathBG);
        else
            disp('Error: no background file found for');
            disp(filePath);
        end
    end
            
    if isempty(filePathBG)
        dset = BrukerRead(filePath);
    else
        dset = subtract_background(filePath, filePathBG);
    end
    
    pars = dset.Properties.UserData;
    
    title = matlab.lang.makeValidName(strcat(fileNames{i}));
    
    if isfield(pars, 'z_axis')
        steps = repmat(pars.z_axis, width(dset) - 1);
    end

    output_table.(strcat(title, '_x')) = dset{:,1};
    if slice && isfield(pars, 'z_axis')
        output_table.(strcat(title, '_y')) = dset{:,2:end}(:, steps == slice);
    else
        output_table.(strcat(title, '_y')) = dset{:,2:end};
    end
end

end