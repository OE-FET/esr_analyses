% BATCH PROCESSING post aquisition

global Path

%file with QValues
[fileName, Path] = uigetfile([Path, '*.txt'],'Please select a file with Q-values.');
Path2QVal = [Path, fileName];

if fileName == 0
    return
end

[fileName, Path] = uigetfile([Path, '*.DSC'], 'Please select ESR data files.', 'MultiSelect', 'on');

if ischar(fileName)
    fileName = {fileName};
end 

for i = 1:length(fileName)
    filePath = [Path fileName{i}];
    try
        ReplaceTitleString(filePath);
        AddTemperaturePars(filePath);
        ReplaceQValue(filePath, Path2QVal);
    catch err
        disp(err);
        disp(filePath);
    end
end
