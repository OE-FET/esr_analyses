function [y] = lorentzian(x, x0, FWHM, varargin)
%LORENTZIAN calculates a 1D Gaussian profile or its nth derivative. The 
% area is normalized to 1.
%
%   This function uses an explicit experession for the 1st derivate and
%   analytical differentiation for higher derivatives. It therefore becomes
%   much slower for n > 1.
%
%   SYNTAX:
%   [y] = lorentzian(x, x0, FWHM)
%   [y] = lorentzian(x, x0, FWHM, 'deriv', n)
%
%   INPUT(S):
%   x - as x-axis values
%   x0 - position of the line center
%   FWHM - full-width-at-half-maximum 
%
% 	OUTPUT(S):
%   y - Lorentzian profile
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2018/07/05 12:58 $    $Revision: 0.1 $

%%

n = get_varargin(varargin, 'deriv', 0);

hwhm = FWHM/2;  % HWHM
x = x - x0;

if n == 0
    y = hwhm./(pi*(x.^2 + hwhm^2));
elseif n == 1
    y = -2*x*hwhm ./(pi*(x.^2 + hwhm^2).^2);
else
    syms L(x_data);
    L(x_data) = hwhm./(pi*(x_data.^2 + hwhm^2));
    dL(x_data) = diff(L(x_data), x_data, n);
    l = @(x) double(dL(x));
    y = l(x);
end
