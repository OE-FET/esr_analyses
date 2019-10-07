function [res] = is_powersat_exp(pars)
%Returns True if the given parameters belong to a powersat experiment,
%False otherwise.

if isfield(pars,'YNAM') && strcmp(pars.YNAM, 'Microwave Power')
    res = true;
else
    res = false;
end

end

