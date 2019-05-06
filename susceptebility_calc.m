function [Chi] = susceptebility_calc(doubleIntArea, pars)
%SUSCEPTEBILITY_CALC Susceptebility calculation from the integrated EPR
%intensity.
%
%   Determines the magnetic susceptibility times unit volume from the
%   double-integrated intensity of a 1st harmonic cw-EPR spectrum.
%
%   SYNTAX:
%   Chi = SUSCEPTEBILITY_CALC(doubleIntArea, pars)
%
%   INPUT(S):
%   doubleIntArea - Integrated area of absorption EPR spectrum or double 
%                   integral from 1st harmonic spectrum
%   pars - structure containing measurement parameters from EPR experiment
%
%   OUTPUT(S):
%   Chi - magnetic susceptebility in m^3
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2019/05/06 12:58 $    $Revision: 1.1 $

%% INPUT PROCESSING
% get cavity calibration factor from pars or prompt for input
if ~isfield(pars, 'ConvFact')
    pars.ConvFact = input(['Please input the cavity specific ' ...
        'conversion factor (from spin counting calibration):\n' ...
        '[default = 9.2710e-09]\n']);
    if isempty(pars.ConvFact)
        pars.ConvFact = 9.2710e-09;
    end
end

% get QValue from pars or promt user for input
if ~isfield(pars, 'QValue')
    pars.QValue = input('Please give the cavity Q-value: ');
end

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
if strcmp(pars.YTYP, 'IGD')
    Pmw = pars.z_axis/1000; % MW Power in W
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

end
