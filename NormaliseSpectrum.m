function [xNorm,yNorm,Pars] = NormaliseSpectrum(varargin)
%NORMALISESPECTRUM Normlizes a Bruker EPR spectrum for aquisiation
%parameters
% 	ESR data is normalised for receiver gain, number of scans (Ns = 1),
% 	time constant (Tc = 1 ms). The normalised conditions correspond to the
% 	Xepr "normalised aquisition" option. The resulting spectrum is then
% 	flagged as normlised by setting Pars.Norm = 'True'.
%
% 	INPUT(S):
% 	NORMALISESPECTRUM()           - prompts user for path input via GUI
% 	NORMALISESPECTRUM(Path)       - loads data from specified path
% 	NORMALISESPECTRUM(x,y,Pars)   - uses given data in (x,y) and experimental
%                                 conditions from Pars
%
% 	OUTPUT(S): 
% 	x_norm                        - B field [gauss] 
% 	y_norm                        - normalised ESR signal intensity
%
% 	DEPENDENCIES:
% 	BrukerRead.m
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2018/07/05 12:58 $    $Revision: 0.1 $

%% Input analyses
% load file
narginchk(0,3);
switch nargin
    case 0
        [x,y,Pars] = BrukerRead;
    case 1
        [x,y,Pars] = BrukerRead(varargin{1});
    case 3
        x=varargin{1};y=varargin{2};Pars=varargin{3};
    otherwise
        error('Inputting only two arguments is not accepted.');
end

%% normalise y-axis
yNorm = y;
if strcmp(Pars.SctNorm,'False') == 1 % check if y-axis already has been normalised
    %--------------------------------------------------------------------
    yNorm = 4.0134*y/(20*10^(Pars.RCAG/20)*Pars.AVGS*Pars.SPTP*1000);
    %--------------------------------------------------------------------
end

xNorm=x;
% Flag spectrum as normalised to prevent second normalization
Pars.SctNorm = 'True';
end