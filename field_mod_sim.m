function yModInPhase = field_mod_sim(x, y, modAmp, harmonic)
%FIELD_MOD_SIM simulates the distortions of lock-in detection on a signal.
%
%   yModInPhase = FIELD_MOD_SIM(x, y, modAmp);
%   yModInPhase = FIELD_MOD_SIM(x, y, modAmp, harmonic);
%
%   Computes the effect of x-axis modulation with amplitude ModAmp and
%   n-th harmonic detection on a signal (x,y). This code implements the
%   psuedo-modulation approach from:
%
%   Hyde, J. S., et al. Applied Magnetic Resonance 1, 483-496 (1990)
%
%   INPUT(S):
%   x      - x-axis data vector, must be uniformliy spaced
%   y      - signal vector
%   modAmp - peak-to-peak modulation amplitude [same units as x]
%   harmonic - harmonic (0, 1, 2, ...); default is 1
%
%   OUTPUT(S):
%   - yModInPhase: in-phase component of pseudo-modulated spectrum
%
%   If no output variable is given, field_mod_sim plots the
%   original and the modulated spectrum.
%
%   Example:
%
%     x = linspace(300,400,1001);
%     y = lorentzian(x,342,4);
%     field_mod_sim(x,y,20);
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$

import esr_analyses.*
import esr_analyses.utils.*

%% input checks
if (nargin == 0)
    help(mfilename)
    return
end

Display = (nargout == 0);

if (nargin < 3) || (nargin > 4), error('Wrong number of input arguments!'); end
if (nargout < 0), error('Not enough output arguments.'); end
if (nargout > 1), error('Too many output arguments.'); end


if (nargin < 4)
    harmonic = 1;
end
if numel(harmonic) ~= 1 || (harmonic < 0) || ~isreal(harmonic) || mod(harmonic, 1)
    error('Harmonic must be a positive integer (1, 2, 3, etc)!');
end

if (modAmp <= 0)
    error('Modulation amplitude (3rd argument) must be positive.');
end

n = length(x);
if length(y) ~= n
    error('x and y must have the same length!');
end

sizey = size(y);
if all(sizey ~= 1)
    error('y must be a row or column vector.');
end

isRowVector = (sizey(1) == 1);
y = y(:);

%% calculate convolution
dx = x(2) - x(1); % get spacing of x-axis, assuming uniform spacing
Ampl = modAmp / 2 / dx; % modulation amplitude in multiples of x-axis steps

NN = n; % zero padding length for fft
ffty = fft(y, NN); % fourier transform signal
ffty(ceil(NN / 2) + 1:end) = 0; % truncate upper half with zeros

% multiply with modulation kernel and fourier transform back
S = (0:NN - 1)' / NN;
yMod = ifft(ffty .* besselj(harmonic, 2 * pi * Ampl * S));
yMod = yMod(1:n); % truncate data (reverse zero padding)

yMod = (1i) ^ harmonic * yMod; % swap real and imaginary parts depending on harmonic

if isRowVector
    yMod = yMod';
end

yModInPhase = real(yMod); % get in-phase data from real part

% plot if output variable is specified
if (Display)
    subplot(3, 1, 1);
    plot(x, y);
    title('Original spectrum');
    subplot(3, 1, [ 2, 3 ]);
    plot(x, yModInPhase);
    xlabel('Magnetic field [G]');
    title(sprintf('Modulated spectrum, harmonic %d, modulation amplitude %g G', harmonic, modAmp));
end

end