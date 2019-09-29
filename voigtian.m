function [y] = voigtian(x, x0, FWHM_gauss, FWHM_lorentz, varargin)
%VOIGTIAN calculates a Voight profile or its n-th derivative.
%
%   The area of the Voigtian is normalized to 1.
%
%   SYNTAX:
%   [y] = VOIGTIAN(x, x0, FWHM_gauss, FWHM_lorentz)
%   [y] = VOIGTIAN(x, x0, FWHM_gauss, FWHM_lorentz, 'deriv', 1)
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

deriv = get_kwarg(varargin, 'deriv', 0);

% padding with additional values for derivative
N = length(x); dx=mean(diff(x));
if deriv == 0
    xNew = x;
else
    xNew(1+deriv:N+deriv) = x;
    xNew(1:deriv) = x(1)-dx*(deriv:-1:1);
    xNew(N+deriv+1:N+2*deriv) = x(N)+dx*(1:deriv);
end

% convert to dimensionless coordinates and calculate voigtian from fadf
HWHM_lorentz = FWHM_lorentz/2;
HWHM_gauss = FWHM_gauss/2;

Real = sqrt(log(2)).*(xNew-x0)./(HWHM_gauss);
Imag = sqrt(log(2)).*(HWHM_lorentz/HWHM_gauss);

w = fadf(Real+1j*Imag);
y = sqrt(log(2)/pi)/HWHM_gauss.*real(w);

% calculate n-th derivative numerically
while deriv > 0
    yDeriv = (y(3:N+2)-y(1:N))/(2*dx);
    y = yDeriv;
    deriv = deriv - 1;
end

y = reshape(y, size(x));

end