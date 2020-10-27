function varargout = BrukerRead(varargin)
%BRUKERREAD Load Bruker BE3ST files (.DTA / .DSC / .YGF).
%
%   BRUKERREAD()
%   BRUKERREAD('/path/to/file')
%   dset = BRUKERREAD(...)
%   [x, o, pars] = BRUKERREAD(...)
%
%   BRUKERREAD when run without any inputs, opens a GUI so that the user can
%   open the file themselves. BRUKERREAD can also accept a path to a file as
%   an input if the path is put in 'quotes' and the extension (.DTA) is left
%   off.
%
%   SYNTAX:
%   [x, o, pars] = BrukerRead
%                  GUI load a file
%
%   dset = BrukerRead
%          GUI load a file, return as a table
%
%   [x, o, pars] = BrukerRead('/path/to/file.DTA')
%                  load x,o and info of file.DTA to the workspace
%
%   INPUT:
%   input1     - a string input to the path of a file
%
%   OUTPUT:
%   x    - x-axis
%   o    - ordinate: the measured quantity
%   pars - array of information about the loaded file
%
%   OUTPUT:
%   dset - a matlab table with the measurement data
%

%%                         Input arguments
% ========================================================================

% remember file path from previous function calls
global Path

switch nargin
    case 0
        if Path==0; Path=[]; end
        [file, directory] = uigetfile({'*.DTA;', 'Bruker File (*.DTA)'; '*.*', 'All Files (*.*)'}, 'Load Bruker file', Path);
        Path = directory;
        % if user cancels command nothing happens
        if isequal(file, 0) || isequal(directory, 0)
            Path = [];
            return
        end

        % File name/path manipulation
        address = [directory, file];
        [~, name, ~] = fileparts(address);

    case 1
        address = varargin{1};
        [directory, name, ~] = fileparts(address);

end


%%                         Parameter file
% ========================================================================

% Load .dsc file
file_dsc = [directory '/' name '.DSC'];
fid = fopen(file_dsc, 'r');

if fid < 0
    error('Both *.DSC and *.DTA files are required to open the file.')
end

% Get characters from file
string = fscanf(fid , '%c');

% Close file, free up memory
fclose(fid);

% Convert character array into useful string array (insert line breaks)
lines = strsplit(string, '\n');
parameter_list = char(lines);

% Do basic reading of parameter file
pars = param2struct(parameter_list);


%%                            Data file
% ========================================================================

% Load .dta file
fid = fopen( [directory '/' name '.DTA'], 'r', 'ieee-be.l64');

if fid < 0
    error(['File ''', name, '.DTA'' could not be opened, both *.DTA and *.DSC files are required to open the file.'])
end

[dta, ~] = fread(fid, inf, 'float64'); % TODO: add support for other formats

fclose(fid);

yNum = length(pars.IKKF);  % number of datasets per slice scan
yNumCplx = length(pars.IKKF(pars.IKKF == "CPLX"));

yNumTotal = yNum + yNumCplx;

o = [];

for i = 1:yNumTotal
    o = [o dta(i:yNumTotal:end)];
end

%%                             X-Axes
% ========================================================================

if strcmp(pars.XTYP, 'IDX')  % indexed data
    x = linspace(pars.XMIN, pars.XMIN + pars.XWID, pars.XPTS)';
elseif strcmp(pars.XTYP, 'IGD')  % data points saved in file
    % if exist, load .YGF , convert to usable matrix
    fid = fopen( [directory '/' name '.XGF'], 'r', 'ieee-be.l64');

    if fid < 0
        error('*.XGF file expected but not found.')
    end

    x = fread(fid, inf, 'float64');
else
    x = [];
end

%%                             Y-Axes
% ========================================================================

if strcmp(pars.YTYP, 'IDX')  % indexed data
    y = linspace(pars.YMIN, pars.YMIN + pars.YWID, pars.YPTS)';
elseif strcmp(pars.YTYP, 'IGD')  % data points saved in file
    % if exist, load .YGF , convert to usable matrix
    fid = fopen( [directory '/' name '.YGF'], 'r', 'ieee-be.l64');

    if fid < 0
        error('*.YGF file expected but not found.')
    end

    y = fread(fid, inf, 'float64');
else
    y = [];
end

%%                             Z-Axes
% ========================================================================

if strcmp(pars.ZTYP, 'IDX')  % indexed data
    z = linspace(pars.ZMIN, pars.ZMIN + pars.ZWID, pars.ZPTS)';
elseif strcmp(pars.ZTYP, 'IGD')  % data points saved in file
    % if exist, load .ZGF , convert to usable matrix
    fid = fopen( [directory '/' name '.ZGF'], 'r', 'ieee-be.l64');

    if fid < 0
        error('*.ZGF file expected but not found.')
    end

    z = fread(fid, inf, 'float64');
else
    z = [];
end

%%                         Reshape data
% ========================================================================

dset = table();

for i=1:size(o, 2)

    if ~isempty(z)
        dset.(join(['o', num2str(i)])) = reshape(o(:,i), length(x), length(y),  length(z));
    elseif ~isempty(y)
        dset.(join(['o', num2str(i)])) = reshape(o(:,i), length(x),  length(y));
    else
        dset.(join(['o', num2str(i)])) = o(:,i);
    end

end

%%                         Output arguments
% ========================================================================

dset = [table(x) dset];

yUnits = {};

for k = 1:length(pars.IKKF)
    if strcmp(pars.IKKF{k}, 'CPLX')
        yUnits(end+1:end+2) = {pars.IRUNI{k}; pars.IRUNI{k}};
    else
        yUnits(end+1) = {pars.IRUNI{k}};
    end
end

pars.x_axis = x;
pars.y_axis = y;
pars.z_axis = z;

dset.Properties.VariableUnits = [{pars.XUNI} yUnits];
dset.Properties.UserData = pars;

% Output results according to requests
switch nargout
    case 1
        varargout{1} = dset;
    case 3
        varargout{1} = x;
        varargout{2} = dset{:,2:end};
        varargout{3} = pars;
    otherwise
        varargout{1} = dset;
end

end

function [par_struct] = param2struct(par_list)
    %% generate info structure

    % get number of rows
    [N, ~] = size(par_list);

    % preallocate memory for analysed rows
    Keep        = zeros(N, 1);
    ParaMatrix	= cell(N, 2);

    % separate parameter names from values and save both in array
    for i = 1:N
        [parameter, value] = strtok(par_list(i, :));

        value = strtrim(value);
        ParaMatrix{i, 1} = parameter;

        % if value not numeric, paste as string, else convert to double
        if isnan(str2double(value))
            value = strip(value, 'both', "'"); % strip single quotes
            ParaMatrix{i, 2} = value;
        else
            ParaMatrix{i, 2} = str2double(value);
        end

        % Mark rows that do not contain aquisition parameters.
        % We recognize them because they start with non-letter characters
        % or are empty
        if isempty(parameter)
            Keep(i) = 0;
        else
            isvar = isletter(parameter);
            Keep(i) = isvar(1);
        end
    end

    % delete rows that do not contain aquisition parameters
    ParaMatrix = ParaMatrix(Keep==1, :);

    % save all parameters in info structure
    for i = 1:length(ParaMatrix)
        name = matlab.lang.makeValidName(ParaMatrix{i, 1});
        par_struct.(name) = ParaMatrix{i, 2};
    end

    % convert string arrays
    for name = ["IKKF", "IRFMT", "IRNAM", "IRUNI"]
        par_struct.(name) = strip(split(par_struct.(name), ','), 'both', "'");
    end

    % replace empty units by a.u.
    for k=1:length(par_struct.IRUNI)
        if isempty(par_struct.IRUNI{k})
            par_struct.IRUNI{k} = 'a.u.';
        end
    end

end
