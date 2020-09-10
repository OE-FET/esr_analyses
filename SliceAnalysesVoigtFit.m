function [out_struct, out_table] = SliceAnalysesVoigtFit(varargin)
%SLICEANALYSESVOIGTFIT performs normalization and spin-counting of an ESR
%signal by fitting it to a voigt function.
%
%   Performs a 1D fit of a cw-EPR spectrum to determine spin cohrence time
%   T2, the magnetic susceptibility and the number of spins in the sample.
%   This program assumes a Voigtian resonance line which is the convolution
%   of multiple Lorentzian spin-ensembles with identifical lifetimes T1 and
%   T2 with a Gaussian distribution of resonance fields.
%
%   Field modulation and the 1st harmonic detection of the cw-EPR signal,
%   together with possible distortions from overmodulation, are explicitly
%   accounted for. The MW field distribution in the cavity is taken from
%   the Bruker DSC file and is used when calculating the MW field amplitude
%   over the sample volume.
%
%   INPUT(S):
%	SLICEANALYSESVOIGTFIT()      - opens GUI for file selection
%	SLICEANALYSESVOIGTFIT(dset)  - uses data given by (x,y,pars)
%	...('sigPath')               - reads data from file
%	...('sigPath', 'bgPath')     - reads data and background from file
%
%   OUTPUT(S):
%	out_struct  - structure containing the measurement data and fit results
%   out_table   - fit results in table format
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$

import esr_analyses.*
import esr_analyses.utils.*

dset = load_spectrum_dialog(varargin{:});
[x,y,pars] = dset_to_tuple(dset);

yes = input('Would you like to perform a baseline correction? y/[n] ','s');
if strcmp(yes, 'y')
    y = baseline_corr(x, y);
end

%%                   Assume that we are not in saturation
%%=========================================================================
Bmw = 1e-6;

%%                      Get starting points for fit
%%=========================================================================

fit = pseudo_voigt_fit(x, y, 'deriv', 1);

FWHM_lorentz  = fit.FWHM_lorentz;       % in Gauss
FWHM_gauss    = fit.FWHM_gauss;         % in Gauss

A0   = fit.a/(pars.B0MA*1e4 * 1e4/8 * Bmw);
B0   = fit.x0;                          % in Gauss
T1   = 1e-20;                           % in sec, assume not-saturated
T2   = 1/(gmratio * FWHM_lorentz*1E-4); % in sec

var0 = [A0 B0 T2 FWHM_gauss];           % starting points

%%                          Perform voigt fit
%%=========================================================================

% create fit function and options
fitfunc = @(var, x) abs(var(1))*esr_voigt_simulation(x, abs(var(2)), T1, ...
    abs(var(3)), abs(var(4)), Bmw, pars.B0MA*1e4, 1);
opt = optimset('TolFun', 1e-9, 'TolX', 1e-9, 'PlotFcns', ...
    @optimplotfval, 'MaxFunEvals', 1e10, 'MaxIter', 1e10);

% fit model to data with Nelder Mead algorithm
fitres   = nelder_mead_fit(fitfunc, x, y, var0, opt);
conf_int = standarderror(fitres); % estimate confidence intervals

A     = abs(fitres.coef(1));
B0    = abs(fitres.coef(2));
T2    = abs(fitres.coef(3));
Brms  = abs(fitres.coef(4));

dA    = abs(conf_int(1));
dB0   = abs(conf_int(2));
dT2   = abs(conf_int(3));
dBrms = abs(conf_int(4));

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
    'x', x, 'y', y, 'pars', pars, 'fitres', fitres, ...
    'B0', B0, 'dB0', dB0, 'g', pars.GFactor, 'dg', pars.GFactorErr, ...
    'T2', T2, 'dT2', dT2, 'Brms', Brms, 'dBrms', dBrms, ...
    'Chi', Chi(1), 'dChi', dChi(1), 'NSpin', NSpin(1), 'dNSpin', dNSpin(1));


out_table = struct2table(rmfield(out_struct, {'x', 'y', 'pars', 'fitres'}));

disp(out_table);

end