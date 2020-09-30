function [] = assert_powersat_exp(pars)
%ASSERT_POWERSAT_EXP Raises an error if the specified Xepr parameters do
%	not correspond to a power saturation measurement.

import esr_analyses.utils.*

if ~is_powersat_exp(pars)
    error('The given data is not from a power saturation measurement.');
end

end

