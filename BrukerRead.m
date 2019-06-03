function varargout = BrukerRead(varargin)
%BRUKERREAD Load Bruker BE3ST files (.DTA / .DSC or .spc/.par).
%
%   BRUKERREAD()
%   BRUKERREAD('/path/to/file')
%   [x, y] = BRUKERREAD(...)
%   [x, y, Pars] = BRUKERREAD(...)
%   [x, y, z] = BRUKERREAD(...)
%   [x, y, z, Pars] = BRUKERREAD(...)
%
%   BRUKERREAD when run without any inputs, opens a GUI so that the user can
%   open the file themselves. BRUKERREAD can also accept a path to a file as
%   an input if the path is put in 'quotes' and the extension (.DTA) is left
%   off.
%
%   BRUKERREAD now also works with Bruker .spc/.par from EMX machines -
%   thanks to Muege Aksoyoglu (Uni. Freiburg, DE) who kindly donated some
%   .spc/.par files for testing
%
%   BRUKERREAD outputs a x matrix (magnetic field / time), a y matrix
%   (intensity) and an optional info field.
%
%   BRUKERREAD is fully functional and extensively tested with cw experiments
%   and 3 dimensional cw experiments where there is an additional .YGF file -
%   such as power saturation experiments. BRUKERREAD provides functionality
%   for pulsed experiments, it has been extensively tested with field swept
%   echoes, fourier induced decays and PELDOR traces. HYSCORE and ENDOR
%   functionality was added in version 13.01 but has had little testing.
%   Other complex pulsed experiments are untested and may have varying
%   results but should yield complete arrays. Please email me if you come
%   across any errors.
%
%   Inputs:
%   input1     - a string input to the path of a file
%
%   Outputs:
%   output0    - plot
%                A new figure showing the data
%   output1    - x axis
%                Magnetic field / time
%   output2    - y axis
%                Intensity
%   output3    - information
%                Array of information about the loaded file
%   output4    - z axis
%                For 3 dimensional pulsed experiments - such as HYSCORE,
%                then the intensity is in the 3rd "z" dimension
%
%   Example: 
%   [x, y] = BrukerRead
%            GUI load a file
%
%   [x, y, Pars] = BrukerRead('/path/to/file.DTA')
%                  load x,y and info of file.DTA with to the workspace and 
%                  plot x,y as a new figure
%
%   MODIFICATIONS
%   This version of BrukerRead automatically detects which parameters have
%   been saved and reads them out. This is done by the function param2struct
%   at the end of the file. Modification by Sam Schott, Oct 2016.
%

%%                         Input arguments
% ========================================================================

% remember file path from previous function calls
global Path

switch nargin
    case 0
        if Path==0; Path=[]; end
        [file, directory] = uigetfile({'*.DTA;*.spc', 'Bruker File (*.DTA,*.spc)'; '*.*', 'All Files (*.*)'}, 'Load Bruker file', Path);
        Path = directory;
        % if user cancels command nothing happens
        if isequal(file, 0) || isequal(directory, 0)
            Path = [];
            return
        end
                
        % File name/path manipulation
        address = [directory, file];
        [~, name, extension] = fileparts(address);
                      
    case 1
        address = varargin{1};
        [directory, name, extension] = fileparts(address);

end


%%                         Parameter files
% ========================================================================

% Load .dsc/.par file
switch extension
    case '.spc'
        file_par = [directory '/' name '.par'];
        fid = fopen(file_par, 'r');
        
        if fid < 0
            error('Both *.spc and *.par files are required to open the file.')            
        end
                
    case {'.DTA', '.DSC'}
        file_dsc = [directory '/' name '.DSC'];
        fid = fopen(file_dsc, 'r');
        
        if fid < 0
            error('Both *.DSC and *.DTA files are required to open the file.')
        end
                
end

% Get characters from file
string = fscanf(fid , '%c');

% Close file, free up memory
fclose(fid);

% Convert character array into useful string array (insert line breaks)
lines = strsplit(string, '\n');
parameter_list = char(lines);

% Do basic reading of parameter file
par_struct = param2struct(parameter_list);


%%                            Data files
% ========================================================================

% Load .dta/.spc file
switch extension
    case '.spc'
        fid   = fopen( [directory '/' name '.spc'], 'r');
        
        if fid < 0
            error(['File ''',name,'.spc'' could not be opened, both *.spc and *.par files are required to open the file.'])
        end
        
        [y, ~] = fread(fid, inf, 'float');
        
    case {'.DTA', '.DSC'}
        fid = fopen( [directory '/' name '.DTA'], 'r', 'ieee-be.l64');
        
        if fid < 0
            error(['File ''', name, '.DTA'' could not be opened, both *.DTA and *.DSC files are required to open the file.'])
        end
        
        [y, ~] = fread(fid, inf, 'float64');
        
        if strcmp(par_struct.IKKF, 'CPLX')
            y = complex(y(1:2:end), y(2:2:end));
        end
end

fclose(fid);


%%                      Magnetic field / X - axis
% ========================================================================

switch extension
    case '.spc'
        
        MagField.min      = par_struct.HCF - (par_struct.HSW / 2);
        MagField.max      = par_struct.HCF + (par_struct.HSW / 2);
        MagField.sampling = par_struct.HSW / par_struct.ANZ;
        
        x = (MagField.min:MagField.sampling:MagField.max)';
        
    case {'.DTA', '.DSC'}
        
        par_struct.XSTEP	= par_struct.XWID / par_struct.XPTS;
        par_struct.XMAX	= par_struct.XMIN + par_struct.XWID - par_struct.XSTEP;
        
        x = (par_struct.XMIN:par_struct.XSTEP:par_struct.XMAX)';
        
end


%%                          Data / Y - axis
% ========================================================================

% Some work required for Pulsed experiments, cw experiments fine

if strcmp(par_struct.EXPT, 'PLS')
    
    % PELDOR, FSE and FID require splitting into real and imaginary
    % channels, these should have no Y axis data.
    
    if strcmp(par_struct.YTYP, 'NODATA')
    
        z = reshape(y, 2, []);
        clear y;
        y.real = z(1, :)';
        y.imag = z(2, :)';
    
    % HYSCORE obviously have Y axis data collection and require splitting
    
    elseif strcmp(par_struct.YTYP, 'IDX')
        
        z = y;
        clear y
        
        % Create Y-axis
        % =============
        
        % Create other Y points
        par_struct.YSTEP	= par_struct.YWID / par_struct.YPTS;
        par_struct.YMAX	= par_struct.YMIN + par_struct.YWID - par_struct.YSTEP;
        
        y = (par_struct.YMIN:par_struct.YSTEP:par_struct.YMAX)';
        
        % Format Z-axis
        % =============
        
        z = reshape(z, size(y, 1), []);
        
    end
end


%%                         YGF files / Z - axis
% ========================================================================

% search directory for .YGF file
if exist([directory '/' name '.YGF'], 'file')
    
    % if exist, load .YGF , convert to usable matrix
    fid = fopen( [directory '/' name '.YGF'], 'r', 'ieee-be.l64');
    
    if fid < 0
        error('BrukerRead: a *.YGF file was found in the folder but could not be opened. BrukerRead will now abort. Please remove the file from the folder or check its permissions.')
    end
    
    [par_struct.z_axis, par_struct.z_axis_points] = fread(fid, inf, 'float64');
    
    % reshape the y-axis into columns using number of data points
    y = reshape(y, par_struct.XPTS , []);
    
end


%%                         Output arguments
% ========================================================================

% Output results according to requests
switch nargout
    case 1
        varargout{1} = y;
    case 2
        varargout{1} = x;
        varargout{2} = y;
    case 3
        if strcmp(par_struct.EXPT, 'PLS') && strcmp(par_struct.YTYP, 'IDX')
            varargout{1} = x;
            varargout{2} = y;
            varargout{3} = z;
        else
            varargout{1} = x;
            varargout{2} = y;
            varargout{3} = par_struct;
        end
        
    case 4
        varargout{1} = x;
        varargout{2} = y;
        varargout{3} = z;
        varargout{4} = par_struct;
        
    otherwise
        varargout{1} = x;
        varargout{2} = y;
end

end

function [info] = param2struct(parameter_list)
    %% generate info structure
    
    % get number of remaining rows
    [N, ~] = size(parameter_list);
    
    % preallocate memory for analysed rows, ignore last (empty) row
    Keep        = zeros(N, 1);
    ParaMatrix	= cell(N, 2);
    
    % separate parameter names from values and save both in array
    for i = 1:N
        [parameter, value] = strtok(parameter_list(i, :));
        
        value = strtrim(value);
        ParaMatrix{i, 1} = parameter;
        
        % if value not numeric, paste as string, else convert to double
        if isnan(str2double(value))
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
        info.(ParaMatrix{i, 1}) = ParaMatrix{i, 2};
    end
end
