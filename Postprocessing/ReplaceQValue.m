function [] = ReplaceQValue(filePathData, filePathQValues)
%% [] = ReplaceQValue(filePathData, filePathQValues)
% Adds an independently determined Q-Value read from filePathQValues to the
% parameter list in the '.DSC' file. The Q-Value is selected by matching
% the measurement temperature and choosing the closest time stamp.
%
% Input arguments:
% filePathData - .DTA/.DSC files with measurement data (saved by Xepr)
% filePathQValues - .txt file with Q-Values in json format (saved by CustomXepr)
%
%% Input analyses 

switch nargin
    case 0
        [Name, Path] = uigetfile('*.DTA');
        filePathData = [Path, Name];
        [Name, Path] = uigetfile('*.txt');
        filePathQValues = [Path, Name];
    case 1
        [Name, Path] = uigetfile('*.txt');
        filePathQValues = [Path, Name];
end

%% Find the right Q-Value from .txt file

% Get measurement temperature from file name or .DSC entry

[folderData, fileNameData] = fileparts(filePathData);

Path2Par = [folderData '/' fileNameData '.DSC'];
Path2Data = [folderData '/' fileNameData '.DTA'];

[~, ~, Pars] = BrukerRead(Path2Data);
try
    T = str2double(strtrim(regexprep(Pars.Temperature,'K','')));
catch
    [startIndex,endIndex] = regexp(Path2Data,'[0-9]+K');
    T = str2double(Path2Data(startIndex:endIndex-1));
    Pars.Temperature = [num2str(T), ' K'];
end

% get measurement time from .DSC file
data_measurement_time = datetime([Pars.DATE ' ' Pars.TIME],'InputFormat', 'MM/dd/yy HH:mm:ss');

% read QValue for the measurement temperature from text file 

fileQValues = readtable(filePathQValues);

qValues = fileQValues.QValue;
qValueTemps = fileQValues.Temperature_K_;
qValuesTimes = fileQValues.TimeStamp;

% get all QValues for measurement temperature
qValueSlice = qValues(qValueTemps == T,:);
timeStampSlice = qValuesTimes(qValueTemps == T,:);

if isempty(qValueSlice)
    error(['No QValue data foud for T = ' T 'K'])
end

% get the time differences between measurement time and QValue time stamp
time_diff = data_measurement_time - timeStampSlice;

% select the Q-Value measurement preceeding the ESR measurement
qValueSlice = qValueSlice(time_diff > 0); % get only preceeding Q-Values
timeStampSlice = timeStampSlice(time_diff > 0, :); % get only preceeding times
time_diff = time_diff(time_diff>0); % get only positve time differences

[~, idx] = min(time_diff);

NewQValue = qValueSlice(idx);

%% Modify DSC file

% read parameter file
A = regexp( fileread(Path2Par), '\n', 'split');
[~,d] = size(A);

% find QValue entry or create new entry
if isfield(Pars, 'QValue')
    for i = 1:d
        pos = strfind(A{i}, num2str(Pars.QValue));
        if pos > 0
            I = i;
        end
    end
    % replace recorded with averaged QValue
    A{I} = sprintf('QValue             %0.1f', NewQValue);
else
    % find the right position for QValue entry
    for i = 1:d
        pos = strfind(A{i}, 'PowerAtten');
        if pos > 0
            I = i;
        end
    end
    % create new line at I+1 for QValue
    A(I+1:end+1) = A(I:end);
    A{I+1} = sprintf('QValue             %0.1f', NewQValue);
end


% write changes to DSC file
fid = fopen(Path2Par, 'w');
fprintf(fid, '%s\n', A{:});
fclose(fid);

end