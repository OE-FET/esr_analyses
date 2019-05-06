function [argout] = SliceAnalysesNum(varargin)
%ESRANALYESNUM performs normalization and spin-counting of an ESR signal by
%numercial double integration.
%
%   The ESR signal is normalised according to measurement conditions. If
%   required, a background signal can be subtracted before performing the
%   analyses.
%
%   The total ESR intensity and spin susceptibility are determined by
%   numerical double-integration. 
%
%   INPUT(S):
%   ESRAnalysesNUM()            - prompts user for spectrum file
%   ...NUM('Path')              - path to file with ESR data
%   ...NUM('PathSIG','PathBG')  - path to signal, path to background
%   ...NUM(x,y,Pars)            - field, signal, and spectral params
%
%   OUTPUT(S):
%   argout                      - output structure containing (x, y, pars)
%                                 and the calculated number of spins and
%                                 susceptibility
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2018/07/05 12:58 $    $Revision: 1.1 $

close all

[x, y, pars] = load_spectrum_dialog(varargin);

%%                         Perform numercial analyses
%%=========================================================================

pars.GFactor = gfactor_determination(x, y, pars, 'plot', 'y');

doubleIntArea = double_int_num(x, y, 'baseline', 'y');

Chi = susceptebility_calc(doubleIntArea, pars);
NSpin = spincounting(doubleIntArea, pars);

%%                                Output
%%=========================================================================

argout.x        = x;
argout.y        = y;
argout.pars     = pars;

argout.Chi      = Chi;
argout.NSpin    = NSpin;

end