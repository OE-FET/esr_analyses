function [output, output_table] = BatchAnalyses(analysesFunc)
%% BATCH PROCESSING

% Select analyses function, e.g., PowerSatAnalysesVoigtFit
if nargin < 1
    analysesFunc = @PowerSatAnalysesNum;
end

% Try automatical matching of background signals
bckgrndStrQ = input('Would you like to subtract background spectra? [y]/n ', 's');

global Path

[fileNames, Path] = uigetfile([Path, '*.DSC'], 'Please select ESR data files.', 'MultiSelect', 'on');

nFiles = numel(fileNames);
output = cell(1, nFiles);

% get sample height and position
global_pars = get_sample_position(struct());

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
    
    pars.SampleL = global_pars.SampleL;
    pars.SampleH = global_pars.SampleH;
    
    output{i} = analysesFunc(x, y, pars);
    out_analyses = rmfield(output{i}, {'x', 'y', 'pars', 'fitres'});
    if ~exist('output_table', 'var')
        output_table = struct2table(out_analyses);
    else
        output_table = [output_table; struct2table(out_analyses)];
    end
    
    pause(2) % pause for user to view plots
end

clc; disp(output_table);

end