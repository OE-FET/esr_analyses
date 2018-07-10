function [y, Area] = PseudoVoigtDeriv(x, a, w, x0, s)
%PSEUDOVOIGTDERIV Pseudo-Voigt derivative curve
%
% 	INPUTS:
% 	a - peak to peak amplitude
% 	w - peak to peak line width 
% 	x0 - center
% 	s - line shape ("gauss-ness")
%
%	OUTPUTS:
%	y - derivative curve data data
%	Area - area under Pseudo-Voigt peak
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2018/07/05 12:58 $    $Revision: 0.1 $

xgauss = (x - x0)/w;
agauss = a*exp(1/2);
gauss = -agauss*xgauss.*exp(-2 * xgauss.^2);

acorr = (3/4)^(3/2);
wcorr = 2*sqrt(1/3);
alorentz = a/acorr;
xlorentz = wcorr*(x-x0)/w;

lorentz = -alorentz * xlorentz ./ (1+xlorentz.^2).^2;

y = s*gauss + (1-s)*lorentz;

AreaLorentz = a*pi*w^2 / sqrt(3); 
AreaGauss =  1/4*a*sqrt((exp(1)*pi)/2)*w^2;

Area = s*AreaGauss + (1-s)*AreaLorentz;

end