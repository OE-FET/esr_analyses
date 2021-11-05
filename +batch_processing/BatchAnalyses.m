function [output, output_table, errors] = BatchAnalyses(analysesFunc, col_names)
%% BATCHANALYSES Analyses ESR spectra from multiple files
%
%   INPUTS(S):
%   analysesFunc - Handle of function to use for analyses, e.g.,
%                  '@PowerSatAnalysesVoigtFit'.
%   col_names    - Cell array with names of parameters that should be
%                  prepended to the returned table.
%
%   OUTPU(S):
%   output          - cell array with all analyses results
%   output_table    - table with all fitting results
%
%   EXAMPLE USGAE:
%   [output, output_table, errors] = BatchAnalyses(@SliceAnalysesLorentzFit, {'Temperature'})
%
%

import esr_analyses.*
import esr_analyses.utils.*

global Path

[fileNames, Path] = uigetfile([Path, '*.DSC'], ...
    'Please select ESR data files.', 'MultiSelect', 'on');

if ~iscell(fileNames)
    fileNames = {fileNames};
end

nFiles = numel(fileNames);
output = cell(1, nFiles);
errors = cell(1, nFiles);

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

% get sample height and position
global_pars = get_sample_position(struct());

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

    dset.Properties.UserData.SampleL = global_pars.SampleL;
    dset.Properties.UserData.SampleH = global_pars.SampleH;

    try
        output{i} = analysesFunc(dset);
    catch ME
        errors{i} = ME;
        warning(strcat('Analyses failed for ', fileNames{i}))
        disp(ME)
        continue
    end

    out_analyses = rmfield(output{i}, {'x', 'o', 'pars'});

    try
        out_analyses = rmfield(out_analyses, {'fitres'});
    end

    if ~exist('output_table', 'var')
        output_table = struct2table(out_analyses);
    else
        output_table = [output_table; struct2table(out_analyses)];  %#ok
    end

    if i ~= nFiles
        input('Press enter to continue.');
        close all
    end
end

if nargin==2
    cols = nan([length(output), length(col_names)]);
    names = cell(size(col_names));
    for i=1:length(col_names)
        names{i} = matlab.lang.makeValidName(col_names{i});
        for j=1:nFiles
            if ~isempty(output{j}) && isfield(output{j}.pars, col_names{i})
                cols(j, i) = output{j}.pars.(col_names{i});
            end
        end
    end
    cols = rmmissing(cols);
    prepend_table = array2table(cols, 'VariableNames', names);
    output_table = [prepend_table, output_table];

end

clc; disp(output_table);

end
