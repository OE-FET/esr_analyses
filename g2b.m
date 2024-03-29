function B = g2b(g, MWFQ)
%G2B Converts a g-factor to a magnetic field value in Tesla.
%
%   [B] = G2B(B, MWFQ) converts a g-factor to a field value B in Tesla for
%   a given MWFQ. Uses the microwave frquency in Hz given in MWFQ.
%
%   [B] = G2B(B) Same as above but uses a MWFQ of 9.4e9 Hz.
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$

import esr_analyses.*
import esr_analyses.utils.*

if nargin < 2
    MWFQ = 9.4*1e9;
end

B = planck * MWFQ ./(g * bmagn); % returns value in Tesla

end