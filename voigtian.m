function [y] = voigtian(x, x0, FWHM_gauss, FWHM_lorentz, varargin)
%VOIGTiAN calculates a Voight profile or its 1st derivative. The area is
% normalized to 1.
%
%   SYNTAX:
%   [y] = voigtian(x, x0, FWHM_gauss, FWHM_lorentz)
%   [y] = voigtian(x, x0, FWHM_gauss, FWHM_lorentz, 'deriv', 1)
%
%   INPUT(S):
%   x - as x-axis values
%   x0 - position of the line center
%   FWHM_gauss - FWHM of the Gaussian component
%   FWHM_lorentz - FWHM of the Lorentzian component
%
% 	OUTPUT(S):
%   y - voigt profile
%
% 27-December-2013 N. Cherkasov
% Comments and questions to: n.b.cherkasov@gmail.com

%%

deriv = get_varargin(varargin, 'deriv', 0);

% padding with additional values for derivative
N = length(x); dx=mean(diff(x));
if deriv==0
    xNew = x;
elseif deriv==1
    xNew(2:N+1) = x;
    xNew(1) = x(1)-dx;
    xNew(N+2) = x(N)+dx;
end 

% converting to dimensionless coordinates

HWHM_lorentz = FWHM_lorentz/2;
HWHM_gauss = FWHM_gauss/2;

Real = sqrt(log(2)).*(xNew-x0)./(HWHM_gauss);
Imag = sqrt(log(2)).*(HWHM_lorentz/HWHM_gauss);

w = fadf(Real+1j*Imag);
y = sqrt(log(2)/pi)/HWHM_gauss.*real(w);

% calculate first derivative
if deriv==1
    yDeriv = (y(3:N+2)-y(1:N))/(2*dx);
    y = yDeriv;
end

y = reshape(y, size(x));

end