function Bmw = get_mw_fields(pars, varargin)
%GET_MW_FIELDS Gets microwave fields from measurement parameter structure.
%
%   Converts MW powers given in pars structure to MW magnetic field
%   amplitudes, given the measurement conditions from pars. In case of
%   power saturation measurements, a vector is returned. Otherwise,
%   GET_MW_FIELDS returns a scalar.
%
% 	SYNTAX:
% 	Bmw = GET_MW_FIELDS(pars)
% 	Bmw = GET_MW_FIELDS(pars, 'cal', 2.2)
%
%   INPUT(S):
%   pars  - structure containing measurement conditions, sample dimensions, etc.
%   'cal' - calibration factor for cavity, defaults to 2.86 to Bruker SHQE cavity
%           Oxford Instruments ESR900 cryostat under vacuum
%
%   OUTPUT(S):
%   Bmw  - Microwave magnetic field amplitude(s) in Tesla
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$

import esr_analyses.*
import esr_analyses.utils.*

%% Get microwave powers based on whether this is a PowerSat measurement
if is_powersat_exp(pars)
    Pmw = pars.y_axis*1e-3; % MW Power in W
else
    Pmw = pars.MWPW; % MW Power in W
end

cal  = get_kwarg(varargin, 'cal', 2.2*1.3);
Qref    = 8355;

% convert MWPW to magnetic field strength in Tesla

pars    = get_sample_position(pars);
f_mean  = mw_mean(pars);
Bmw     = f_mean * cal * sqrt(Pmw) * sqrt(pars.QValue/Qref) * 1E-4; % in Tesla

end