function [NSpin] = spincounting(doubleIntArea, pars, S)
%SPINCOUNTING Spin-counting from the integrated EPR intensity.
%   
%   Determines the absolute number of spins in a sample from the
%   double-integrated intensity of a 1st harmonic cw-EPR spectrum assuming 
%   a Curie susceptibility. This assumption is valid for independent spins
%   when T >> ?E/2. X-Band EPR typically operates around B = 350 mT while
%   a field of up to 0.5 T corresponds to ?E/2 = 28.9 ueV. This still 
%   remains far lower than the corresponding thermal energy of 450 µeV at
%   5 K. 
%
%   If the input argument S is not given, S = 1/2 is used.
%
%   SYNTAX:
%   [NSpin] = spincounting(doubleIntArea, pars)
%   [NSpin] = spincounting(doubleIntArea, pars, S)
%
%   OUTPUT(S):
%   NSpin - number of spins
%
%   DEPENDENCIES:
%   mw_mean.m
%   Natural Constants
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2018/07/05 12:58 $    $Revision: 0.1 $

%% INPUT PROCESSING
if nargin < 3; S = 1/2; end

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


% get measurement temperature
if isfield(pars, 'Temperature')
    T = str2double(strtrim(regexprep(pars.Temperature,'K','')));
else
    T = input('Please give measurement temperature in K [default = 298 K]: ');
    if isempty(T); T = 298; end
end

%% CALCULATE NSPIN
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

% cavity and MW bridge calibration factors
k = 200/(pars.BridgeCalib * pars.ConvFact);

% boltzmann factor
boltzman_factor = planck*pars.MWFQ / (boltzm*T);

% -------------------------------------------------------------------------
NSpin = k * doubleIntArea ./ (pars.QValue * sqrt(Pmw) * pars.B0MA * ...
     S*(S+1) * boltzman_factor * position_correction);
% -------------------------------------------------------------------------

end