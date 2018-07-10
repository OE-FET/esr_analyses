function [y] = Voigt(x, x0, FWHMGauss, FWHMLorentz, deriv)

% voigt  Calculation of the VOIGT profile 
%
%   [y] = voigt(x, x0, FWHMGauss, FWHMLorentz, deriv)
%   The function calculates a Voight profile or its derivative.
%
%
%   INPUT ARGUMENTS
%       x - as x-axis values
%       x0 - position of the line center
%       FWHMGauss - parameter of the width of the Gaussian component (full-width at half maximum)
%       FWHMLorentz - parameter of the width of the Lorentzian component (full-width at half maximum)
%       deriv - 0 for voigt profile, 1 for first derivative
%
% 	OUTPUT
%       y - array of intensities
%
%
% 27-December-2013 N. Cherkasov
% Comments and questions to: n.b.cherkasov@gmail.com

% padding with additional values for derivative
N=length(x);dx=mean(diff(x));
if deriv==0
    xNew = x;
elseif deriv==1
    xNew(2:N+1) = x;
    xNew(1) = x(1)-dx;
    xNew(N+2) = x(N)+dx;
end 

% converting to dimensionless coordinates

widthLorentz=FWHMLorentz/2;
widthGauss=FWHMGauss/2;

Real=sqrt(log(2)).*(xNew-x0)./(widthGauss);
Imag=sqrt(log(2)).*(widthLorentz/widthGauss);

w = fadf(Real+1j*Imag);
y = sqrt(log(2)/pi)/widthGauss.*real(w);

% calculate first derivative
if deriv==1
    yDeriv=(y(3:N+2)-y(1:N))/(2*dx);
    y=yDeriv;
end

end