function [x_norm, y_norm, pars] = normalise_spectrum(x, y, pars)
%NORMALISESPECTRUM Normlizes a Bruker EPR spectrum for aquisiation parameters.
%
% 	ESR data is normalised for receiver gain, number of scans (Ns = 1),
% 	time constant (Tc = 1 ms). The normalised conditions correspond to the
% 	Xepr "normalised aquisition" option. The resulting spectrum is then
% 	flagged as normlised by setting Pars.Norm = 'True'.
%
% 	OUTPUT(S):
% 	x_norm                        - B field [gauss]
% 	y_norm                        - normalised ESR signal intensity
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2019/05/06 12:58 $    $Revision: 1.1 $

y_norm = y;
if strcmp(pars.SctNorm, 'False') == 1 % check if y-axis already has been normalised
    %--------------------------------------------------------------------
    y_norm = 4.0134*y/(20*10^(pars.RCAG/20)*pars.AVGS*pars.SPTP*1000);
    %--------------------------------------------------------------------
end

x_norm=x;

% Flag spectrum as normalised to prevent second normalization
pars.SctNorm = 'True';

end