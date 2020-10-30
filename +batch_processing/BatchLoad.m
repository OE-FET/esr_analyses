function [output_table] = BatchLoad(varargin)
% BATCHLOAD Loads and returns data from multiple ESR data files
%
%   KEYWORD INPUT(S):
%   yslice    - slice to take from y-axis (integer of boolen array)
%
%   OUTPU:
%   output_table    - table with all data
%

import esr_analyses.*
import esr_analyses.utils.*

yslice = get_kwarg(varargin, 'yslice', false);

global Path

[fileNames, Path] = uigetfile([Path, '*.DSC'], ...
    'Please select ESR data files.', 'MultiSelect', 'on');

if ~iscell(fileNames)
    fileNames = {fileNames};
end

nFiles = numel(fileNames);
output_table = table;

bckgrndStrQ = input('Would you like to subtract background spectra? y/[n]  ', 's');

if strcmp(bckgrndStrQ, 'y')

    [fileNamesBG, PathBG] = uigetfile([Path, '*.DSC'], ...
        'Please select background data files.', 'MultiSelect', 'on');

    if ~iscell(fileNamesBG)
        fileNamesBG = {fileNamesBG};
    end

    nFilesBG = numel(fileNamesBG);

    if nFiles ~= nFilesBG
        error('Got %d data files but %d background files', nFiles, nFilesBG);
    end
end

for i = 1:nFiles

    filePath = [Path fileNames{i}];

    if strcmp(bckgrndStrQ, 'y')
        filePathBG = [PathBG fileNamesBG{i}];
    end

    if strcmp(bckgrndStrQ, 'y')
        dset = subtract_background(filePath, filePathBG);
    else
        dset = BrukerRead(filePath);
    end

    pars = dset.Properties.UserData;

    title = matlab.lang.makeValidName(strcat(fileNames{i}));

    output_table.(strcat(title, '_x')) = dset{:,1};
    if yslice && isfield(pars, 'y_axis')
        output_table.(strcat(title, '_y')) = dset{:,2:end}(:, yslice);
    else
        output_table.(strcat(title, '_y')) = dset{:,2:end};
    end
end

end