function [NSpin, dNSpin] = spincounting(doubleIntArea, pars, varargin)
%SPINCOUNTING Spin-counting from the integrated EPR intensity.
%   
%   Determines the absolute number of spins in a sample from the
%   double-integrated intensity of a 1st harmonic cw-EPR spectrum assuming 
%   a Curie susceptibility. This assumption is valid for independent spins
%   when T >> E_Zeeman/2. X-Band EPR typically operates around B = 350 mT
%   while a field of up to 0.5 T corresponds to E_Zeeman/2 = 28.9 µeV. This  
%   still remains far lower than the corresponding thermal energy of
%   450 µeV at 5 K. 
%
%   If the input argument S is not given, S = 1/2 is used.
%
%   SYNTAX:
%   NSpin = SPINCOUNTING(doubleIntArea, pars)
%   [NSpin, dNSpin] = SPINCOUNTING(doubleIntArea, pars, 'S', S, 'dA', dA)
%
%   OUTPUT(S):
%   NSpin  - number of spins
%   dNSpin - standard error

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2019/05/06 12:58 $    $Revision: 1.1 $

%% INPUT PROCESSING
S = get_kwarg(varargin, 'S', 1/2);
if nargout > 1; dA = get_kwarg(varargin, 'dA'); end

pars = get_par(pars, 'ConvFact', 9.2710e-09);  % get calibration factor
pars = get_par(pars, 'QValue'); % get or ask for QValue
pars = get_par(pars, 'Temperature', 298); % get or ask for temperature
if nargout > 1; pars = get_par(pars, 'QValueErr'); end
if ischar(pars.Temperature)
    pars.Temperature = str2double(strtrim(regexprep(pars.Temperature,'K','')));
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
if strcmp(pars.YNAM, 'Microwave Power')
    Pmw = pars.z_axis/1000; % MW Power in W
else
    Pmw = pars.MWPW;
end

% cavity and MW bridge calibration factors
k = 200/(pars.BridgeCalib * pars.ConvFact);

% boltzmann factor
boltzman_factor = planck*pars.MWFQ / (boltzm*pars.Temperature);

% -------------------------------------------------------------------------
NSpin = k * doubleIntArea ./ (pars.QValue * sqrt(Pmw) * pars.B0MA * ...
     S*(S+1) * boltzman_factor * position_correction);
% -------------------------------------------------------------------------

if nargout > 1
    dNSpin = pars.QValueErr .* NSpin./pars.QValue + dA .* NSpin ./ doubleIntArea;
end

end
