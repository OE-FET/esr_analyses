function g = b2g(B, MWFQ)
%B2G Converts B value in Tesla to g-value for given MWFQ.
%
%   [g] = B2G(B) converts a field value B in Tesla to a g-value using a
%   default MWFQ of 9.4*1e9 Hz.
%
%   [g] = B2G(B, MWFQ) uses the given MWFQ in Hz for the conversion.
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2019/05/06 12:58 $    $Revision: 1.1 $

if nargin < 2
    MWFQ = 9.4e9;
end

g = planck * MWFQ./(B * bmagn);

end