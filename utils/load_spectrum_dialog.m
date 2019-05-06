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
            [x, y, pars] = BrukerRead;
        end
    case 2
        [x, y, pars] = subtract_background(argcell{1},  argcell{2});
    case 3
        x = argcell{1}; y = argcell{2}; pars = argcell{3};
end

% get sample position in cavity
pars = get_sample_position(pars);

% get measurement temperature
if ~isfield(pars, 'Temperature')
    T = input('Please give the temperature in K [default = 298 K]: ');
    if isempty(T); T = 298; end
    pars.Temperature = sprintf('%.1f K', T);
end

% get Q value
if ~isfield(pars, 'QValue')
    qValue = input('Please give the cavity QValue: ');
    pars.QValue = qValue;
end

[x, y, pars] = normalise_spectrum(x, y, pars);

end