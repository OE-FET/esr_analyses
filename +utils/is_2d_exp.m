function [res] = is_2d_exp(input)
%IS_2D_EXP Returns true for 2D datasets, false otherwise.

if istable(input)
    pars = input.Properties.UserData;
elseif isstruct(input)
    pars = input;
else
    error('Either dataset table or parameter structure required.')
end

if isfield(pars,'YNAM') && isfield(pars, 'y_axis')
    res = true;
else
    res = false;
end

end