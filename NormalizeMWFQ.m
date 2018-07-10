function [x,y,Pars] = NormalizeMWFQ(x,y,Pars,varargin)
%NORMALIZEMWFQ Normalises magnetic field axis to given microwave frequency
% 	[x,y,Pars] = NORMALIZEMWFQ(x, y, Pars,'MWFQ', freq) normalises the data
% 	given by (x,y,Pars) to a microwave frequency given by freq. The new data
%	and parameters are output again as (x,y,Pars).
%
%	[x,y,Pars] = NORMALIZEMWFQ(x, y, Pars) normalises the data given by 
% 	(x,y,Pars) to a microwave frequency of 9.4 GHz.
%
% 	Dependencies:
% 	g2b.m
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2018/07/05 12:58 $    $Revision: 0.1 $


freq = 9.4e9; nInput=length(varargin);
for i=1:nInput
    if strcmp(varargin{i},'MWFQ')
        freq=varargin{i+1};
        nInput=nargin-2;
    end
end

B_offset = 10000* ( g2b(gfree, freq) - g2b(gfree, Pars.MWFQ) );

Pars.MWFQ = freq;

x = x + B_offset;

end