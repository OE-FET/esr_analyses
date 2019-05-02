function [x_norm, y_norm, pars] = normalise_spectrum(x, y, pars)
%NORMALISESPECTRUM Normlizes a Bruker EPR spectrum for aquisiation
%parameters
% 	ESR data is normalised for receiver gain, number of scans (Ns = 1),
% 	time constant (Tc = 1 ms). The normalised conditions correspond to the
% 	Xepr "normalised aquisition" option. The resulting spectrum is then
% 	flagged as normlised by setting Pars.Norm = 'True'.
%
% 	INPUT(S):
% 	NORMALISESPECTRUM(x, y, Pars)   - uses given data in (x, y) and 
%                                   experimental conditions from Pars
%
% 	OUTPUT(S): 
% 	x_norm                        - B field [gauss] 
% 	y_norm                        - normalised ESR signal intensity
%
% 	DEPENDENCIES:
% 	BrukerRead.m
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2018/07/05 12:58 $    $Revision: 0.1 $

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