function f = gb2f(g, B)
% Calculates microwave frequency in Hz from g-factor and resonance center
% in Tesla.
%
% Dependencies:
% "Natureal constants" folder
%

f = g*bmagn*B / planck;

end