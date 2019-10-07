function [Chi, dChi] = susceptibility_calc(doubleIntArea, pars, varargin)
%SUSCEPTEBILITY_CALC Susceptebility calculation from the integrated EPR
%intensity.
%
%   Determines the magnetic susceptibility times unit volume from the
%   double-integrated intensity of a 1st harmonic cw-EPR spectrum.
%
%   The second return argument is the standard error of the susceptiblity,
%   calculated from error propagration with the error in 'doubleIntArea'
%   and the QValue error.
%
%   SYNTAX:
%   Chi = SUSCEPTEBILITY_CALC(doubleIntArea, pars)
%   [Chi, dChi] = SUSCEPTEBILITY_CALC(doubleIntArea, pars, 'dA', areaError)
%
%   INPUT(S):
%   doubleIntArea - Integrated area of absorption EPR spectrum or double
%                   integral from 1st harmonic spectrum
%   pars - structure containing measurement parameters from EPR experiment
%   dA   - standard error for doubleIntArea, used for error propagation
%
%   OUTPUT(S):
%   Chi  - magnetic susceptebility in m^3
%   dChi - standard error of magnetic susceptebility in m^3
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2019/05/06 12:58 $    $Revision: 1.1 $

import esr_analyses.*
import esr_analyses.utils.*

%% INPUT PROCESSING
pars = get_par(pars, 'ConvFact', 9.2710e-09);  % get calibration factor
pars = get_par(pars, 'QValue'); % get or ask for QValue
if nargout > 1; pars = get_par(pars, 'QValueErr'); end
if nargout > 1; dA = get_kwarg(varargin, 'dA'); end

%% CALCULATE SUSCEPTIBILITY
try
    position_correction = mw_mean(pars);
catch
    fprintf(['MW field distribution in the cavity could not be read ' ...
        'from the DSC file.\nProceeding without correcting for the '...
        'position of the sample in the cavity.\n \n']);
    position_correction = 1;
end

% get MW power(s)
if is_powersat_exp(pars)
    Pmw = pars.z_axis*1e-3; % MW Power in W
else
    Pmw = pars.MWPW;
end

% get g-factor from pars
sample_g = pars.GFactor;

% cavity and MW bridge calibration factors
k = 200/(pars.BridgeCalib * pars.ConvFact);

% normalization factor for measurement conditions
norm = (pars.QValue * sqrt(Pmw) * pars.B0MA * position_correction);

% -------------------------------------------------------------------------
Chi = mu0* k * doubleIntArea .* sample_g.^2 .* bmagn^2 ./ (3 * ...
    planck * pars.MWFQ  * norm);
% -------------------------------------------------------------------------

if nargout > 1
    dChi = dA .* Chi./doubleIntArea + pars.QValueErr .* Chi./pars.QValue;
end

end
