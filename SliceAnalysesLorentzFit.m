function [out_struct, out_table] = SliceAnalysesLorentzFit(varargin)
%SLICEANALYSESLORENTZFIT performs normalization and spin-counting of an
%ESR signal by fitting it to a Lorentz function.
%
%   Performs a 1D fit of a cw-EPR spectrum to determine spin cohrence time
%   T2, the magnetic susceptibility and the number of spins in the sample.
%   This program assumes a Lorentzian resonance line from a single
%   spin-ensemble.
%
%   Field modulation and the 1st harmonic detection of the cw-EPR signal,
%   together with possible distortions from overmodulation, are explicitly
%   accounted for. The MW field distribution in the cavity is taken from
%   the Bruker DSC file and is used when calculating the MW field amplitude
%   over the sample volume.
%
%   INPUT SYNTAX:
%	SLICEANALYSESLORENTZFIT()      - opens GUI for file selection
%	SLICEANALYSESLORENTZFIT(dset)  - uses data given by dset
%	...(x,o,pars)                  - uses data given by [x,o,pars]
%	...('sigPath')                 - reads data from file
%	...('sigPath', 'bgPath')       - reads data and background from file
%
%   KEYWORD INPUT(S):
%   plot        - if true, plot data and best fit at each iteration
%
%   OUTPUT(S):
%	out_struct  - structure containing the measurement data and fit results
%   out_table   - fit results in table format
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$

import esr_analyses.*
import esr_analyses.utils.*

[plotting, varargin] = get_kwarg(varargin, 'plot', true);

dset = load_spectrum_dialog(varargin{:});
[x,o,pars] = slice_experiment(dset);

yes = input('Would you like to perform a baseline correction? y/[n] ','s');
if strcmp(yes, 'y')
    o = baseline_corr(x, o);
end

%%                   Assume that we are not in saturation
%%=========================================================================
Bmw = get_mw_fields(pars);

%%                      Get starting points for fit
%%=========================================================================

fit = lorentz_fit(x, o, 'deriv', 1);

FWHM  = fit.FWHM;       % in Gauss

A0    = fit.a/(pars.B0MA*1e4 * 1e4/8 * Bmw);
B0    = fit.x0;                  % in Gauss
T1    = 1e-20;                   % in sec, assume not-saturated
T2    = 1/(gmratio * FWHM*1E-4); % in sec

var0  = [A0 B0 T2];              % starting points

%%                          Perform lorentz fit
%%=========================================================================

% create fit function and options
fitfunc = @(var, x) abs(var(1))*esr_lorentz_simulation(x, abs(var(2)), T1, ...
    abs(var(3)), Bmw, pars.B0MA*1e4, 1);
opt = optimset('TolFun', 1e-9, 'TolX', 1e-9, 'MaxFunEvals', 1e10, 'MaxIter', 1e10);

% fit model to data with Nelder Mead algorithm
fitres   = nelder_mead_fit(fitfunc, x, o, var0, opt, 'plot', plotting);
conf_int = standarderror(fitres); % estimate confidence intervals

A     = abs(fitres.coef(1));
B0    = abs(fitres.coef(2));
T2    = abs(fitres.coef(3));

dA    = abs(conf_int(1));
dB0   = abs(conf_int(2));
dT2   = abs(conf_int(3));

%%                           Plot results
%%=========================================================================

h = plot(fitres);

xlabel(h{1}(1).Parent, 'Magnetic field [G]')
ylabel(h{1}(1).Parent, 'ESR signal [a.u.]')
title(h{1}(1).Parent, pars.TITL, 'interpreter', 'none')

%%                      Susceptibility Calculation
%%=========================================================================
pars.GFactor    = b2g(B0*1e-4, pars.MWFQ);
pars.GFactorErr = dB0 * b2g(B0*1e-4, pars.MWFQ) / B0;
modScaling      = pars.B0MA*1e4 * 1e4/8; % scaling for pseudo-modulation

areaDI          = modScaling * Bmw .* A;
areaDIerror     = modScaling * Bmw .* dA;

[Chi, dChi]     = susceptibility_calc(areaDI, pars, 'dA', areaDIerror);
[NSpin, dNSpin] = spincounting(areaDI, pars, 'dA', areaDIerror);

%%                                Output
%%=========================================================================

% create output structure

out_struct = struct(...
    'x', x, 'o', o, 'pars', pars, 'fitres', fitres, ...
    'B0', B0, 'dB0', dB0, 'g', pars.GFactor, 'dg', pars.GFactorErr, ...
    'T2', T2, 'dT2', dT2, 'Chi', Chi(1), 'dChi', dChi(1), ...
    'NSpin', NSpin(1), 'dNSpin', dNSpin(1));


out_table = struct2table(rmfield(out_struct, {'x', 'y', 'pars', 'fitres'}));

disp(out_table);

end