function [x, y, pars] = normalise_mwfq(x, y, pars, varargin)
%NORMALISE_MWFQ Normalises magnetic field axis to given microwave frequency
%
%	SYNTAX:
% 	[x, y, pars] = NORMALISE_MWFQ(x, y, pars, 'MWFQ', freq) normalises the
% 	data given by (x,y,pars) to a microwave frequency given by freq in Hz.
%	The new data and parameters are output again as tuple (x y, pars).
%
%	[x, y, pars] = NORMALISE_MWFQ(x, y, pars) normalises the data given by
% 	(x, y, pars) to a microwave frequency of 9.4 GHz.
%
%   KEYWORD INPUT:
%	'MWFQ'  - The microwave frequency to use. Defaults to 9.4e9 Hz.
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$

import esr_analyses.*
import esr_analyses.utils.*

freq = get_kwarg(varargin, 'MWFQ', 9.4*1e9);

B_offset = 10000 * ( g2b(gfree, freq) - g2b(gfree, pars.MWFQ) );

pars.MWFQ = freq;

x = x + B_offset;

end