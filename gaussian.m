function [y] = gaussian(x, x0, FWHM, varargin)
%GAUSSIAN 1D Gaussian profile or its n-th derivative.
%
%   Returns a 1D Gaussian with an area 1 or its n-th derivative:
%
%   sigma = FWHM/sqrt(8*log(2))
%   y = exp(-(x-x0)^2/(2*sigma^2)) ./ (sigma*sqrt(2*pi))
%
%   This function uses an explicit experession for the 1st derivate and
%   hermite polynomial expressions for higher orders. Evaluation therefore
%   becomes much slower for n > 1.
%
%   SYNTAX:
%   [y] = GAUSSIAN(x, x0, FWHM)
%   [y] = GAUSSIAN(x, x0, FWHM, 'deriv', n)
%
%   INPUT:
%   x - as x-axis values
%   x0 - position of the line center
%   FWHM - full-width-at-half-maximum
%
%   KEYWORD INPUT:
%   'deriv' - If given and > 0, this function returns the n-th
%             derivative of a Voigtian.
%
% 	OUTPUT:
%   y - Gaussian profile
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$

import esr_analyses.*
import esr_analyses.utils.*

%%

n = get_kwarg(varargin, 'deriv', 0);

sigma = FWHM/sqrt(8*log(2));
x = x - x0;

if n == 0
    y = exp(-x.^2/(2*sigma^2)) ./ (sigma*sqrt(2*pi));
elseif n == 1
    y = -x./sigma^2 .* exp(-x.^2/(2*sigma^2)) ./ (sigma*sqrt(2*pi));
else
    deriv_factor = hermiteH(n, x/sigma*sqrt(2)) * (-1/sigma*sqrt(2))^n;
    y = deriv_factor .* exp(-x.^2/(2*sigma^2)) ./ (sigma*sqrt(2*pi));
end

end
