function [x, y, pars] = normalize_mwfq(x, y, pars, varargin)
%NORMALIZE_MWFQ Normalises magnetic field axis to given microwave frequency
% 	[x, y, pars] = normalize_mwfq(x, y, pars, 'MWFQ', freq) normalises the
% 	data given by (x,y,pars) to a microwave frequency given by freq in Hz.  
%	The new data and parameters are output again as tuple (x y, pars).
%
%	[x, y, pars] = normalize_mwfq(x, y, pars) normalises the data given by 
% 	(x, y, pars) to a microwave frequency of 9.4 GHz.
%
% 	Dependencies:
% 	g2b.m
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2018/07/05 12:58 $    $Revision: 0.1 $

freq = get_varargin(varargin, 'MWFQ', 9.4*1e9);

B_offset = 10000* ( g2b(gfree, freq) - g2b(gfree, pars.MWFQ) );

pars.MWFQ = freq;

x = x + B_offset;

end