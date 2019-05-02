function [yn] = ESRLorentzSimulation(x, B0, T1, T2, Bmw, n, ModAmp)
%ESRLORENTZSIMULATION simulates the n-th harmonic detection of a Lorentzian
%resonance with x-axis modulation.   
%   Simulates the n-th harmonic detection of a Lonretzian ESR resonance line:
%
%   L = Bmw / ( 1 + Bmw^2*T1*T2*gmratio^2 + (xpd-B0)^2*T2^2*gmratio^2 )
%
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
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2018/07/05 12:58 $    $Revision: 1.1 $

% convert Gauss to Tesla
x = x*1e-4;
B0 = B0*1e-4;

% extend x-range
x1 = x(1,:);

len2 = round(length(x1)/2);
xpd = padarray(x, [0, len2]);

dx = median(diff(x1));
xpd(:,1:len2) = ones(size(x,1), 1) * linspace(min(x1)-len2*dx, min(x1)-dx, len2);
xpd(:,len2+length(x1)+1:end) = ones(size(x,1),1) * linspace(max(x1)+dx, max(x1)+dx*len2, len2);

[~, Bmw] = meshgrid(xpd(1,:), Bmw(:,1));

% Calculate ESR absorption signal
xpd1 = xpd(1,:);

ypd_0 = Bmw./(1 + Bmw.^2*T1*T2*gmratio^2 + (xpd-B0).^2*T2^2*gmratio^2);
ypd_n = zeros(size(xpd));


% Calculate field modulated signal
for i=1:size(ypd_0, 1)
    ypd_n(i,:) = field_mod_sim(xpd1, ypd_0(i,:), ModAmp, n);
end

% truncate edges
yn = ypd_n(:, len2+1:len2+length(x1))';
end