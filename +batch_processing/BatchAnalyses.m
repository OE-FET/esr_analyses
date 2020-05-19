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

[fileNames, Path] = uigetfile([Path, '*.DSC'], ...
    'Please select ESR data files.', 'MultiSelect', 'on');

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

    if ~strcmp(bckgrndStrQ, 'n')
        filePathBGcandidates = {'Vg_00.', 'Vg_0.'};
        
        for candidate=filePathBGcandidates
            filePathBGcandidate = regexprep(filePath, 'Vg_(\S+)(\.)', ...
                candidate{1});
            if exist(filePathBGcandidate, 'file') == 2
                filePathBG = filePathBGcandidate;
                break
            end
        end
  
        if filePathBG
            disp('Background file for:');
            disp(filePath);
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

    dset.Properties.UserData.SampleL = global_pars.SampleL;
    dset.Properties.UserData.SampleH = global_pars.SampleH;
    % use a default QValueErr of 100 if none is given
    dset.Properties.UserData.QValueErr = get_par(dset.Properties.UserData, ...
                                                 'QValueErr', 100);

    try
        output{i} = analysesFunc(dset);
    catch ME
        errors{i} = ME;
        warning(strcat('Analyses failed for ', fileNames{i}))
        disp(ME)
        continue
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
    cols = rmmissing(cols);
    prepend_table = array2table(cols, 'VariableNames', names);
end

output_table = [prepend_table, output_table];

clc; disp(output_table);

end