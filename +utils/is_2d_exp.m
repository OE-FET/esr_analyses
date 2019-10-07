function [res] = is_2d_exp(pars)
%Returns True if the specified Xepr parameters correspond to a
%2D measurement, False otherwise.

if isfield(pars,'YNAM') && isfield(pars, 'z_axis')
    res = true;
else
    res = false;
end

end