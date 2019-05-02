function yn = ESRVoigtSimulation(x, B0, T1, T2, Bmw, Brms, n, ModAmp)
%ESRVOIGTIMULATION simulates the n-th harmonic detection of a Voigtian
%resonance with x-axis modulation.  
%   Simulates the n-th harmonic detection of a Voigtian ESR resonance line.
%
%   INPUT:
%   x - vector with external magnetic fields [Gauss]
%   B0 - resonance center [Gauss]
%   T1 - spin lattice relaxation time [sec]
%   T2 - spin coherence time [sec]
%   Bmw - vector with microwave magnetic field amplitudes [Tesla]
%   ModAmp - field modulation amplitude [Gauss]
%   n - n-th harmonic detection
%
%   OUTPUT:
%   yn - simulated n-th harmonic spectrum
%
%   DEPENDENCIES:
%   field_mod_sim.m
%   voigt.m
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2018/07/05 12:58 $    $Revision: 1.1 $

% convert Gauss to Tesla
x = x*1e-4;
B0 = B0*1e-4;
Brms = Brms*1e-4;
ModAmp = ModAmp*1e-4;

% extend x-range
xstart = round(-0.2*length(x));
xstop = round(1.2*length(x));
if size(x,2)==1
    xpd = interp1(x', xstart:xstop, 'linear', 'extrap');
else
    xpd = interp1(x', xstart:xstop, 'linear', 'extrap')';
end


%% Calculate ESR voigt signal (centered around zero)

ypd_0 = zeros(size(xpd))';
ypd_n = zeros(size(xpd))';
xpd1 = xpd(1,:);
Bmw1 = Bmw(:,1);

for i=1:size(ypd_0,2)

    FWHMLorentz = 2/(gmratio*T2) * sqrt(1 + gmratio^2*Bmw1(i).^2*T1*T2);

    FWHMGauss = 2*sqrt(2*log(2))*Brms;
    FWHMGauss = max(FWHMGauss, 1e-8);  % need finite value for voigt function

    ypd_0(:,i) = Bmw1(i)/FWHMLorentz * voigtian(xpd1, B0, FWHMGauss, FWHMLorentz);
    ypd_n(:,i) = field_mod_sim(xpd1, ypd_0(:,i), ModAmp, n);
end

% truncate edges
yn = reshape(ypd_n(ismember(xpd, x)'), fliplr(size(x)));

end