function [y] = lorentzian(x, x0, FWHM, varargin)
%LORENTZIAN calculates a 1D Gaussian profile or its n-th derivative.
%
%   Returns a 1D Lorentzian with an area of 1 or its n-th derivative:
%
%   HWHM = FWHM/2
%   y = HWHM/(pi*(x^2 + HWHM^2))
%
%   This function uses an explicit expression for the 1st derivate and
%   analytical differentiation for higher derivatives. It therefore becomes
%   much slower for n > 1.
%
%   SYNTAX:
%   [y] = LORENTZIAN(x, x0, FWHM)
%   [y] = LORENTZIAN(x, x0, FWHM, 'deriv', n)
%
%   INPUT:
%   x       - as x-axis values
%   x0      - position of the line center
%   FWHM    - full-width-at-half-maximum
%
%   KEYWORD INPUT:
%   'deriv' - If given and > 0, this function returns the n-th
%             derivative of a Voigtian.
%
% 	OUTPUT:
%   y - Lorentzian profile
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$

import esr_analyses.*
import esr_analyses.utils.*

%%

import esr_analyses.*
import esr_analyses.utils.*

n = utils.get_kwarg(varargin, 'deriv', 0);

HWHM = FWHM/2;
x = x - x0;

if n == 0
    y = HWHM./(pi*(x.^2 + HWHM^2));
elseif n == 1
    y = -2*x*HWHM ./(pi*(x.^2 + HWHM^2).^2);
elseif n > 1
    syms L(x_data);
    L(x_data) = HWHM./(pi*(x_data.^2 + HWHM^2));
    dL(x_data) = diff(L(x_data), x_data, n);
    l = @(x) double(dL(x));
    y = l(x);
end

end
