function f = gb2f(g, B)
% Calculates microwave frequency in Hz from g-factor and resonance center
% in Tesla.
%
% Dependencies:
% "Natureal constants" folder
%

import esr_analyses.*
import esr_analyses.utils.*

f = g*bmagn*B / planck;

end