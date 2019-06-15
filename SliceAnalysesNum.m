function [argout] = SliceAnalysesNum(varargin)
%SLICEANALYSESNUM Numercial analyses of cw-EPR spectrum.
%
%   Numerically integrates a cw-EPR spectrum to determine the magnetic
%   susceptibility and the number of spins.
%
% 	WARNING: Using numerical integration may underestimate the tails of EPR
% 	spectra. This can lead to significant errors for long-tailed resonance
% 	shapes, such as Lorentzians, if the SNR ratio is small or the
% 	measurement range is less than 5 times the peak-to-peak linewidth.
%
%	INPUTS:
%	SLICEANALYSESNUM()         - opens GUI for file selection
%	SLICEANALYSESNUM(x,y,pars) - uses data given by (x,y,pars)
%	...('sigPath')             - reads data from file
%	...('sigPath', 'bgPath')   - reads data and background from file
%
%	OUTPUT:
%	argout  - structure containing the measurement data and fit results 
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2019/05/06 12:58 $    $Revision: 1.1 $

close all

[x, y, pars] = load_spectrum_dialog(varargin);

%%                         Perform numercial analyses
%%=========================================================================

pars.GFactor = gfactor_determination(x, y, pars, 'plot', 'y');

doubleIntArea = double_int_num(x, y, 'baseline', 'y');

Chi   = susceptebility_calc(doubleIntArea, pars);
NSpin = spincounting(doubleIntArea, pars);

%%                                Output
%%=========================================================================

argout.x        = x;
argout.y        = y;
argout.pars     = pars;

argout.Chi      = Chi;
argout.NSpin    = NSpin;

end