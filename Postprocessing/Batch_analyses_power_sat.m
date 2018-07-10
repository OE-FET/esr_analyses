%% BATCH PROCESSING
global Path

% data files
[fileNames, Path] = uigetfile([Path, '*.DTA'], 'Please select ESR data files.', 'MultiSelect', 'on');

bckgrndStrQ = input('Would you like to automatically match background spectra? [y]/n: ', 's');

for i = 1:length(fileNames)
    
    if iscell(fileNames)==0
        filePath = [Path fileNames];
    else
        filePath = [Path fileNames{i}];
    end
    
    if strcmp(bckgrndStrQ,'y')==0
        argout = PowerSatAnalysesVoigtFIT(filePath);
    else
        filePathBG1 = regexprep(filePath,'Vg_(\S+)(\.)','Vg_00.');
        filePathBG2 = regexprep(filePath,'Vg_(\S+)(\.)','Vg_0.');
        try
            argout = PowerSatAnalysesVoigtFIT(filePath, filePathBG1);
        catch
            try
                argout = PowerSatAnalysesVoigtFIT(filePath, filePathBG2);
            catch
                display(['Error: no background file found for \n', filePathBG2]);
            end
        end
    end
            
    outputMatrix(i,:) = [argout.T, 100/argout.T, argout.gfactor, argout.Brms, argout.Bp2p, argout.Chi, argout.dChi, argout.NSpin, argout.dNSpin, argout.T1, argout.dT1, argout.T2, argout.dT2];
end
