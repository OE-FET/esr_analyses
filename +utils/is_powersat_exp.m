function [res] = is_powersat_exp(input)
%Returns True if the given parameters belong to a powersat experiment,
%False otherwise.

if istable(input)
    pars = input.Properties.UserData;
elseif isstruct(input)
    pars = input;
else
    error('Either dataset table or parameter structure required.')
end

if isfield(pars,'YNAM') && strcmp(pars.YNAM, 'Microwave Power')
    res = true;
else
    res = false;
end

end

