function [x_norm, y_norm, pars] = normalise_spectrum(x, y, pars)
%NORMALISESPECTRUM Normlizes a Bruker EPR spectrum for acquisition parameters.
%
% 	ESR data is normalised for receiver gain, number of scans (Ns = 1),
% 	time constant (Tc = 1 ms). The normalised conditions correspond to the
% 	Xepr "normalised aquisition" option. The resulting spectrum is then
% 	flagged as normlised by setting Pars.Norm = 'True'.
%
%   Note: The actual normalisation performed by Xepr and the signal channel
%         depends on the signal channel / signal processing unit version.
%         This function DOES NOT check for the signal channel version but 
%         assumes a SPU which ALWAYS scales for the time constant.
%
%   INPUT(S):
%   x    - magnetic field axis
%   y    - non-normalised signal intensity
%   pars - experimental parameters
%
% 	OUTPUT(S):
% 	x_norm                        - B field [gauss]
% 	y_norm                        - normalised ESR signal intensity
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$

import esr_analyses.*
import esr_analyses.utils.*

y_norm = y;
if strcmp(pars.SctNorm, 'False') % check if y-axis already has been normalised
% -------------------------------------------------------------------------
    y_norm = 44.5423*y/(20*10^(pars.RCAG/20)*pars.AVGS*1000);
% -------------------------------------------------------------------------
end

x_norm=x;

% Flag spectrum as normalised to prevent second normalisation
pars.SctNorm = 'True';

end