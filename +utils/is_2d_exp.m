function [res] = is_2d_exp(input)
%Returns True if the specified Xepr parameters correspond to a
%2D measurement, False otherwise.

if istable(input)
    pars = intput.Properties.UserData;
elseif isstruct(input)
    pars = intput;
else
    error('Either dataset table or parameter structure required.')
end

if isfield(pars,'YNAM') && isfield(pars, 'z_axis')
    res = true;
else
    res = false;
end

end