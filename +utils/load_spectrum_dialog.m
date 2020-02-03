function dset = load_spectrum_dialog(varargin)
%LOAD_SPECTRUM_DIALOG loads a Bruker ESR spectrum and prepares it for
%further analyuses.
%
%   Uses BrukerRead to import the selected spectrum and allows the user to
%   subtract a background spectrum, if desired. Normalises the dataset and
%   determines other experimental parameters not imported by BrukerRead.
%
%   INPUT(S):
%   LOAD_SPECTRUM_DIALOG()         - opens gui to select data file and
%                                    background file, if desired
%   ...('signal_path')             - path to signal data with the option to
%                                    select background data
%   ...('signal_path', 'bg_path')  - path to signal, path to background
%   ...(dset)                      - uses given data set directly
%   ...(x, y, pars)                - uses x and y data with parameter
%                                    structure
%
%   OUTPUT(S):
%	dset
%
%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2019/06/17 10:43 $    $Revision: 1.1 $
%%

import esr_analyses.*
import esr_analyses.utils.*


err_msg = "Input must be either one or two file paths for signal " + ...
          "and background measuements," + newline + "a result table from " + ...
          "BrukerRead or (x, y, pars) data.";

switch nargin
    case 0
        str = input('Would you like to subtract a background signal? y/[n]: ', 's');
        if strcmpi(str,'y')
            dset = subtract_background;
        else
            dset = BrukerRead;
        end
    case 1
        if istable(varargin{1})
            dset = varargin{1};
        elseif ischar(varargin{1})
            str = input('Would you like to subtract a background signal? y/[n]: ', 's');
            if strcmpi(str, 'y')
                dset = subtract_background(varargin{1});
            else
                dset = BrukerRead(varargin{1});
            end
        else
            error(err_msg)
        end
    case 2
        if ischar(varargin{1}) && ischar(varargin{2})
            dset = subtract_background(varargin{1},  varargin{2});
        else
            error(err_msg)
        end
    case 3
        if isvector(varargin{1}) && isvector(varargin{2}) && isstruct(varargin{3})
            x    = varargin{1};
            y    = varargin{2};
            pars = varargin{3};

            dset = table(x, y);
            dset.Properties.UserData = pars;
        else
            error(err_msg)
        end
end

pars = dset.Properties.UserData;

% confirm or ask for missing parameters
pars = get_sample_position(pars);
pars = get_par(pars, 'Temperature', 298);
pars = get_par(pars, 'QValue');
pars = get_par(pars, 'QValueErr');

% convert temperature from string to double if necessary
if ischar(pars.Temperature)
    pars.Temperature = str2double(strtrim(regexprep(pars.Temperature,'K','')));
end

dset.Properties.UserData = pars;

dset = normalise_spectrum(dset);

end