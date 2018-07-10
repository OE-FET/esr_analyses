function [y, Area] = GaussDeriv(x, a, w, x0)
%GAUSSDERIV outputs a Gaussian derivative curve
%   [y, Area] = GAUSSDERIV(X, A, W, X0) outputs the derivative of a Gaussian
%   with a peak-to-peak amplitude A, a peak-to-peak linewidth of W and a
%   center at X0.
%
%   INPUT:
%   X    - x-axis data
%   A    - peak-to-peak amplitude
%   W    - peak-to-peak linewidth
%   X0   - peak center
%
%   OUTPUT:
%   y    - Gaussian derivative curve
%   Area - Area of Gaussian
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2018/07/05 12:58 $    $Revision: 1.1 $

xnew = (x - x0)/w;
anew = a*exp(1/2);

y = -anew*xnew.*exp(-2 * xnew.^2);

Area = 1/4*a*sqrt((exp(1)*pi)/2)*w^2;

end