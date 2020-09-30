function [] = assert_2d_exp(input)
%ASSERT_2D_EXP Raises an error if the specified Xepr parameters do not correspond to a
%	2D measurement.

import esr_analyses.utils.*

if ~is_2d_exp(input)
    error('The given data is not from a 2D measurement.');
end

end

