function [output_table] = BatchLoad(step, slice)
% BATCHLOAD Loads and returns data from multiple ESR data files
%
%   INPUTS:
%   step    - name of parameter that is stepped between files, e.g.,
%             'Temperature'
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

[fileNames, Path] = uigetfile([Path, '*.DSC'], 'Please select ESR data files.', 'MultiSelect', 'on');

nFiles = numel(fileNames);

output_table = table;

for i = 1:nFiles
    
    if iscell(fileNames)==0
        filePath = [Path fileNames];
    else
        filePath = [Path fileNames{i}];
    end
    
    filePathBG = [];
    
    if strcmp(bckgrndStrQ, 'y')
        filePathBG1 = regexprep(filePath,'Vg_(\S+)(\.)','Vg_00.');
        filePathBG2 = regexprep(filePath,'Vg_(\S+)(\.)','Vg_0.');
        if exist(filePathBG1, 'file') == 2
            filePathBG = filePathBG1;
        elseif exist(filePathBG2, 'file') == 2
            filePathBG = filePathBG2;
        else
            disp('Error: no background file found for \n');
        end
    end
            
    if isempty(filePathBG)
        [x, y, pars] = BrukerRead(filePath);
    else
        [x, y, pars] = subtract_background(filePath, filePathBG);
    end
    
    title = matlab.lang.makeValidName(strcat(step, '_', pars.(step)));
    
    output_table.(strcat(title, '_x')) = x;
    if slice && isfield(pars, 'z_axis')
        output_table.(strcat(title, '_y')) = y(:, pars.z_axis == slice);
    else
        output_table.(strcat(title, '_y')) = y;
    end
end

end