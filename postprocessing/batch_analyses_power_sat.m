%% BATCH PROCESSING

% Select analyses function, e.g., PowerSatAnalysesVoigtFit
analysesFunc = @PowerSatAnalysesNum;

% Try automatical matching of background signals
bckgrndStrQ = input('Would you like to automatically match background spectra? [y]/n ', 's');


global Path

[fileNames, Path] = uigetfile([Path, '*.DTA'], 'Please select ESR data files.', 'MultiSelect', 'on');

nFiles = numel(fileNames);
output = cell(1, nFiles);

for i = 1:nFiles
    
    if iscell(fileNames)==0
        filePath = [Path fileNames];
    else
        filePath = [Path fileNames{i}];
    end
    
    if ~strcmp(bckgrndStrQ,'y')
        argout = analysesFunc(filePath);
    else
        filePathBG1 = regexprep(filePath,'Vg_(\S+)(\.)','Vg_00.');
        filePathBG2 = regexprep(filePath,'Vg_(\S+)(\.)','Vg_0.');
        try
            argout = analysesFunc(filePath, filePathBG1);
        catch
            try
                argout = analysesFunc(filePath, filePathBG2);
            catch
                display(['Error: no background file found for \n', filePathBG2]);
            end
        end
    end
            
    output{i} = argout;
end
