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
%   INPUT SYNTAX:
%	SLICEANALYSESNUM()       - opens GUI for file selection
%	SLICEANALYSESNUM(dset)   - uses data given by dset
%	...(x,o,pars)            - uses data given by [x,o,pars]
%	...('sigPath')           - reads data from file
%	...('sigPath', 'bgPath') - reads data and background from file
%
%	OUTPUT:
%	argout  - structure containing the measurement data and fit results
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$

import esr_analyses.*
import esr_analyses.utils.*

dset = load_spectrum_dialog(varargin{:});
[x,y,pars] = slice_experiment(dset);

yes = input('Would you like to perform a baseline correction? y/[n] ', 's');
if strcmp(yes, 'y')
    y = baseline_corr(x, y);
end

%%                         Perform numercial analyses
%%=========================================================================

pars.GFactor = gfactor_determination(x, y, pars, 'plot', 'y');

intArea = double_int_num(x, y);
dA = intArea*pars.QValueErr/pars.QValue;

[Chi, dChi]     = susceptibility_calc(intArea, pars, 'dA', dA);
[NSpin, dNSpin] = spincounting(intArea, pars, 'dA', dA);

%%                                Output
%%=========================================================================

argout.x        = x;
argout.y        = y;
argout.pars     = pars;

argout.Chi      = Chi;
argout.dChi     = dChi;
argout.NSpin    = NSpin;
argout.dNSpin   = dNSpin;

end