function dset = load_spectrum_dialog(argcell)
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
%
%   OUTPUT(S):
%	dset
%
%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2019/06/17 10:43 $    $Revision: 1.1 $
%%

import esr_analyses.*
import esr_analyses.utils.*

argnum = length(argcell);

switch argnum
    case 0
        str = input('Would you like to subtract a background signal? y/[n]: ', 's');
        if strcmpi(str,'y')
            dset = subtract_background;
        else
            dset = BrukerRead;
        end
    case 1
        
        if istable(argcell{1})
            dset = argcell{1};
        elseif ischar(argcell{1})
            str = input('Would you like to subtract a background signal? y/[n]: ', 's');
            if strcmpi(str, 'y')
                dset = subtract_background(argcell{1});
            else
                dset = BrukerRead(argcell{1});
            end
        else
            error('You must give either a data table or a file path.');
        end
    case 2
        dset = subtract_background(argcell{1},  argcell{2});
end

pars = dset.Proprties.UserData;

% confirm or ask for missing parameters
pars = get_sample_position(pars);
pars = get_par(pars, 'Temperature', 298);
pars = get_par(pars, 'QValue');
pars = get_par(pars, 'QValueErr');

% convert temperature from string to double if necessary
if ischar(pars.Temperature)
    pars.Temperature = str2double(strtrim(regexprep(pars.Temperature,'K','')));
end

dset.Proprties.UserData = pars;

dset = normalise_spectrum(dset);

end