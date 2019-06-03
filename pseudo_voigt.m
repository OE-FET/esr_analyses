function [y] = pseudo_voigt(x, x0, FWHM_gauss, FWHM_lorentz, varargin)
%PSEUDO_VOIGT calculates a Pseudo-Voight profile or its n-th derivative.
%
%   Returns a Pseudo-Voigt peak with an area normalized to 1 or its n-th
%   derivative.
%
%   SYNTAX:
%   [y] = PSEUDO_VOIGT(x, x0, FWHM_gauss, FWHM_lorentz)
%   [y] = PSEUDO_VOIGT(x, x0, FWHM_gauss, FWHM_lorentz, 'deriv', 1)
%
%   INPUT(S):
%   x - as x-axis values
%   x0 - position of the line center
%   FWHM_gauss - FWHM of the Gaussian component
%   FWHM_lorentz - FWHM of the Lorentzian component
%
% 	OUTPUT(S):
%   y - pseudo-voigt profile
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2019/05/06 12:58 $    $Revision: 1.1 $

%%

try n = get_kwarg(varargin, 'deriv'); catch; n = 0; end

f_G = FWHM_gauss;
f_L = FWHM_lorentz;

f = (f_G^5 + 2.69269*f_G^4*f_L + 2.42843*f_G^3*f_L^2 + 4.47163*f_G^2*f_L^3 + 0.07842*f_G*f_L^4 + f_L^5)^(1/5);

eta = 1.36603*(f_L/f) - 0.47719*(f_L/f)^2 + 0.11116*(f_L/f)^3;

y = eta * lorentzian(x, x0, f, 'deriv', n) + (1-eta) * gaussian(x, x0, f, 'deriv', n);

end