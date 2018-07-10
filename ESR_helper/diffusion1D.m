function [y] = diffusion1D(x,x0,width,drvtv)
% voigt  Calculation of the VOIGT profile 
%
%   [y] = voigt(x,x0,widthGauss,widthLorentz,deriv)
%   The function calculates a Voight profile or its derivative.
%
%
%   INPUT ARGUMENTS
%       x - as x-axis values
%       x0 - position of the line center
%       widthGauss - parameter of the width of the Gaussian component (Half-width at half maximum)
%       FWHMLorentz - parameter of the width of the Lorentzian component (full-width at half maximum)
%       deriv - 0 for voigt profile, 1 for first derivative, 2 for second
%       derivative
%
% 	OUTPUT
%       y - array of intensities
%
%
% 27-December-2013 N. Cherkasov
% Comments and questions to: n.b.cherkasov@gmail.com

% padding with additional values for derivative
N=length(x);dx=mean(diff(x));
if drvtv==0
    xNew = x;
elseif drvtv==1
    xNew(2:N+1) = x;
    xNew(1) = x(1)-dx;
    xNew(N+2) = x(N)+dx;
end 

nNew=length(xNew);Interval=(max(xNew)-min(xNew));

% computing the function
a=10*width/Interval;
center = -(mean(xNew)-x0);
t=0:a:(nNew-1)*a;

g=exp(-a*abs(t).^(3/2)).*exp(-2*pi*1i*center*t);

y=fft(g); y=abs(fftshift(y));

% calculate first derivative
if drvtv==1
    yDeriv=(y(3:N+2)-y(1:N))/(2*dx);
    y=yDeriv;
end

end