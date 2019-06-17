function [yn] = ESRLorentzSimulation(x, B0, T1, T2, Bmw, modAmp, n)
%ESRLORENTZSIMULATION simulates the n-th harmonic detection of a Lorentzian
%resonance with x-axis modulation.   
%
%   Simulates the n-th harmonic detection of a Lorentzian ESR resonance
%   line. If no modulation amplitude and harmonic are given, the ESR 
%   absorption signal without field modulation is returned.
%
%   L = Bmw / ( 1 + Bmw^2*T1*T2*gmratio^2 + (xpd-B0)^2*T2^2*gmratio^2 )
%
%
%   INPUT(S):
%   x - vector with external magnetic fields [Gauss]
%   B0 - resonance center [Gauss]
%   T1 - spin lattice relaxation time [sec]
%   T2 - spin coherence time [sec]
%   Bmw - vector with microwave magnetic field amplitudes [Tesla]
%   modAmp - field modulation amplitude [Gauss]
%   n - n-th harmonic detection
%
%   OUTPUT(S):
%   yn - simulated n-th harmonic spectrum
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2018/07/05 12:58 $    $Revision: 1.1 $

% internally, we use SI units ONLY

% convert Gauss to Tesla
x = x*1e-4;
B0 = B0*1e-4;
if nargin > 5; modAmp = modAmp*1e-4; end

%% Calculate ESR voigt signal (centered around zero)
x1 = x(1,:);
y0 = zeros(size(x))';
if nargin > 5; yn = zeros(size(x))'; end
Bmw1 = Bmw(:,1);

for i=1:size(y0, 2)

    saturation_factor = sqrt(1 + gmratio^2*Bmw1(i).^2*T1*T2);
    FWHMLorentz = 2/(gmratio*T2) * saturation_factor;
    area_scaling = Bmw1(i) / saturation_factor;

    y0(:,i) = area_scaling * lorentzian(x1, B0, FWHMLorentz);

    if nargin > 5
        yn(:,i) = field_mod_sim(x1, y0(:,i), modAmp, n);
    else
        yn = y0;
    end
end

end