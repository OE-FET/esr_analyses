function [x, y, pars] = load_spectrum_dialog(argcell)

argnum = length(argcell);

switch argnum
    case 0
        str = input('Would you like to subtract a background signal? y/[n]', 's');
        if strcmp(str,'y')
            [x, y , pars] = subtract_background;
        else
            [x, y, pars]  = BrukerRead;
        end
    case 1
        str = input('Would you like to subtract a background signal? y/[n]', 's');
        if strcmp(str, 'y')
            [x,y,pars]   = subtract_background(argcell{1});
        else
            [x, y, pars] = BrukerRead(argcell{1});
        end
    case 2
        [x, y, pars] = subtract_background(argcell{1},  argcell{2});
    case 3
        x = argcell{1}; y = argcell{2}; pars = argcell{3};
end

% get sample position in cavity
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