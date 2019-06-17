function [x, y, pars] = load_spectrum_dialog(argcell)
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
%   ...(x, y, pars)                - uses given data directly
%
%   OUTPUT(S):
%	x    - magnetic field values
%   y    - normalised intensity values
%   pars - experimental parameters
%
%
%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2019/06/17 10:43 $    $Revision: 1.1 $
%%
argnum = length(argcell);

switch argnum
    case 0
        str = input('Would you like to subtract a background signal? y/[n]: ', 's');
        if strcmpi(str,'y')
            [x, y, pars] = subtract_background;
        else
            [x, y, pars] = BrukerRead;
        end
    case 1
        str = input('Would you like to subtract a background signal? y/[n]: ', 's');
        if strcmpi(str, 'y')
            [x, y, pars] = subtract_background(argcell{1});
        else
            [x, y, pars] = BrukerRead(argcell{1});
        end
    case 2
        [x, y, pars] = subtract_background(argcell{1},  argcell{2});
    case 3
        x = argcell{1}; y = argcell{2}; pars = argcell{3};
end

% confirm or ask for missing parameters
pars = get_sample_position(pars);
pars = get_par(pars, 'Temperature', 298);
pars = get_par(pars, 'QValue');
pars = get_par(pars, 'QValueErr');

% convert temperature from string to double if necessary
if ischar(pars.Temperature)
    pars.Temperature = str2double(strtrim(regexprep(pars.Temperature,'K','')));
end

[x, y, pars] = normalise_spectrum(x, y, pars);

end