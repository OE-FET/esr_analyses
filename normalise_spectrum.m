function [dset] = normalise_spectrum(dset)
%NORMALISESPECTRUM Normlizes a Bruker EPR spectrum for acquisition parameters.
%
% 	ESR data is normalised for receiver gain, number of scans (Ns = 1),
% 	time constant (Tc = 1 ms). The normalised conditions correspond to the
% 	Xepr "normalised acquisition" option. The resulting spectrum is then
% 	flagged as normalised by setting Pars.Norm = 'True'.
%
%   Note: The actual normalisation performed by Xepr and the signal channel
%         depends on the signal channel / signal processing unit version.
%         This function DOES NOT check for the signal channel version but
%         assumes a SPU which ALWAYS scales for the time constant.
%
%   INPUT:
%   dset - dataset from BrukerRead
%
% 	OUTPUT:
% 	dset - normalised dset
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$

import esr_analyses.*
import esr_analyses.utils.*

pars = dset.Properties.UserData;
y = dset{:, 2:end};

if strcmp(pars.SctNorm, 'False') % check if y-axis already has been normalised
    % ---------------------------------------------------------------------------
    dset{:,2:end} = 44.5423 * y ./ (20*10^(pars.RCAG/20) * pars.AVGS * 1000);
    % ---------------------------------------------------------------------------
end

% Flag spectrum as normalised to prevent second normalisation
dset.Properties.UserData.SctNorm = 'True';

end