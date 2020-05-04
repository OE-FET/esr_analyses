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

import esr_analyses.*
import esr_analyses.utils.*

% Try automatical matching of background signals
bckgrndStrQ = input('Would you like to subtract background spectra? [y]/n ', 's');

global Path

[fileNames, Path] = uigetfile([Path, '*.DSC'], 'Please select ESR data files.', 'MultiSelect', 'on');

if ~iscell(fileNames)
    fileNames = {fileNames};
end

nFiles = numel(fileNames);
output = cell(1, nFiles);
errors = cell(1, nFiles);

% get sample height and position
global_pars = get_sample_position(struct());

for i = 1:nFiles

    filePath = [Path fileNames{i}];
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
        dset = BrukerRead(filePath);
    else
        dset = subtract_background(filePath, filePathBG);
    end

    dset.Properties.UserData.SampleL = global_pars.SampleL;
    dset.Properties.UserData.SampleH = global_pars.SampleH;

    try
        output{i} = analysesFunc(dset);
    catch ME
        errors{i} = ME;
        warning('Analyses failed for ')
        disp(ME)
    end

    out_analyses = rmfield(output{i}, {'x', 'y', 'pars', 'fitres'});
    if ~exist('output_table', 'var')
        output_table = struct2table(out_analyses);
    else
        output_table = [output_table; struct2table(out_analyses)];  %#ok
    end
    
    if i ~= nFiles
        input('Press enter to continue.');
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
    prepend_table = array2table(cols, 'VariableNames', names);
end

output_table = [prepend_table, output_table];

clc; disp(output_table);

end