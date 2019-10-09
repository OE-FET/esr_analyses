function g = b2g(B, MWFQ)
%B2G Converts B value in Tesla to g-value for given MWFQ (microwave frequency).
%
%   [g] = B2G(B) converts a field value B in Tesla to a g-value using a
%   default MWFQ of 9.4e9 Hz.
%
%   [g] = B2G(B, MWFQ) uses the given MWFQ in Hz for the conversion.
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$

import esr_analyses.*
import esr_analyses.utils.*

if nargin < 2
    MWFQ = 9.4e9;
end

g = planck * MWFQ./(B * bmagn);

end