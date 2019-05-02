function [y] = gaussian(x, x0, FWHM, varargin)
%GAUSSIAN calculates a 1D Gaussian profile or its nth derivative. The area
% is normalized to 1.
%
%   This function uses an explicit experession for the 1st derivate and
%   hermite polynomial expressions for higher orders. It therefore becomes
%   much slower for n > 1.
%
%   SYNTAX:
%   [y] = gaussian(x, x0, FWHM)
%   [y] = gaussian(x, x0, FWHM, 'deriv', n)
%
%   INPUT(S):
%   x - as x-axis values
%   x0 - position of the line center
%   FWHM - full-width-at-half-maximum 
%
% 	OUTPUT(S):
%   y - Gaussian profile
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2018/07/05 12:58 $    $Revision: 0.1 $

%%

n = get_varargin(varargin, 'deriv', 0);

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
