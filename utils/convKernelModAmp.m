function [y] = convKernelModAmp(x, Bm, n)
%% Modulation Kernel
% Kernel for convolution with ESR absorption spectrum to simulate lock-in
% detection with external field modulaiton. The kernel is the Fourier
% transform of a Bessel function of the first kind and order n. It
% simulates a n-th harmonic detection.
%
% INPUT(S):
% x - external field position
% Bm - modulation amplitude, same units as x
% n - order of harmonic deterction, n = 0,1,2
%
%
% References:
% 1. Hyde, J. S., et al. Appl. Magn. Reson. 1, 483?496 (1990).
% 2. Stoll, S. "Spectral simulations in solid-state electron paramagnetic
%    resonance" (2003)
%

%% Kernel functions
if n==0
    y = 2 ./ sqrt(Bm^2-x.^2);
elseif n==1
    y = - 2* (x/Bm) ./ sqrt(Bm^2-x.^2);
elseif n==2
    y = (4* (x/Bm).^2 -2) ./ sqrt(Bm^2-x.^2);
else
    error('Invalid order of harmonic. Please specify n=0,1,2.');
end

end