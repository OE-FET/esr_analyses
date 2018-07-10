function [y, Area] = LorentzDeriv(x, a, w, x0)
%LORENTZDERIV outputs a Lorentzian derivative curve
%   [y, Area] = LORENTZDERIV(X, A, W, X0) outputs the derivative of a
%   Lorentzian with a peak-to-peak amplitude A, a peak-to-peak linewidth of
%   W and a center at A.
%
%   INPUT:
%   X    - x-axis data
%   A    - peak-to-peak amplitude
%   W    - peak-to-peak linewidth
%   X0   - peak center
%
%   OUTPUT:
%   y    - Lorentzian derivative curve
%   area - Area of Lorentzian
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2018/07/05 12:58 $    $Revision: 0.1 $

acorr = (3/4)^(3/2);
wcorr = 2*sqrt(1/3);
anew = a/acorr;
xnew = wcorr*(x-x0)/w;

y = -anew * xnew ./ (1+xnew.^2).^2;

Area = a*pi*w^2 / sqrt(3); 

end