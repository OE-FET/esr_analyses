function [] = assert_2d_exp(pars)
%Raises an error if the specified Xepr parameters do not correspond to a
%2D measurement.

import esr_analyses.utils.*

if ~is_2d_exp(pars)
    error('The given data is not from a 2D measurement.');
end

end

