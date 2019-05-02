function g = b2g(B, MWFQ)
%B2G Converts B value in Tesla to g-value for given MWFQ
%   [g] = B2G(B, MWFQ) converts a field value B in Tesla to a g-value for a
%   given MWFQ. Uses the microwave frquency in Hz given in MWFQ.
%
%   [g] = B2G(B) Same as above but uses a default MWFQ of 9.4e9 Hz.
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2018/07/05 12:58 $    $Revision: 0.1 $

if nargin < 2
    MWFQ = 9.4e9;
end

g = planck * MWFQ./(B * bmagn);

end