function [y]=HFISimulation(Sys,x,Deriv)

Gamma=Sys.FWHM/2;

if nargin==2
    Deriv=1;
end

[locations, amplitudes]=HFIPeaks(Sys.n,Sys.I,Sys.A);

nPeaks=length(locations);

lorentzians=zeros(length(x),nPeaks);

if Deriv==0
    for i=1:nPeaks
        lorentzians(:,i)=amplitudes(i)*(Gamma/pi) ./ (Gamma^2+(x - locations(i)).^2);
    end
end

if Deriv==1
    for i=1:nPeaks
        lorentzians(:,i)=-amplitudes(i)*2*Gamma*(x-locations(i))./(Gamma^2+(x-locations(i)).^2).^2/pi;
    end
end

y=sum(lorentzians,2);
end