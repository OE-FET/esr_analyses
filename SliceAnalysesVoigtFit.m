function [argout] = SliceAnalysesVoigtFit(varargin)
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
%   SLICEANALYSESVOIGTFIT()     - prompts user for spectrum file
%   ...FIT('Path')              - path to file with ESR data
%   ...FIT('PathSIG', 'PathBG') - path to signal, path to background
%   ...FIT(x, y, Pars)          - field, signal, and spectral params
%
%   OUTPUT(S):
%	argout  - structure containing the measurement data and fit results 
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2019/05/06 12:58 $    $Revision: 1.1 $

close all

[x, y, pars] = load_spectrum_dialog(varargin);

%%                         Calculate MW field
%%=========================================================================

Bmw = get_mw_fields(pars);

%%                      Get starting points for fit
%%=========================================================================

fit = pseudo_voigt_fit(x, y, 'deriv', 1);

FWHM_lorentz  = fit.FWHM_lorentz;       % in Gauss
FWHM_gauss    = fit.FWHM_gauss;         % in Gauss

A0   = fit.a/(pars.B0MA*1e4 * 1e4/8 * Bmw(mid));
B0   = fit.x0;                          % in Gauss
T1   = 1e-20;                           % in sec, assume not-saturated
T2   = 1/(gmratio * FWHM_lorentz*1E-4); % in sec

var0 = [A0 B0 T2 FWHM_gauss];           % starting points

%%                          Perform voigt fit
%%=========================================================================

% create fit function and options
fitfunc = @(var, x) abs(var(1))*ESRVoigtSimulation(x, T1, abs(var(2)), ...
    abs(var(3)), abs(var(4)), Bmw, pars.B0MA*1e4, 1);
opt = optimset('TolFun', 1e-9, 'TolX', 1e-9, 'PlotFcns', ...
    @optimplotfval, 'MaxFunEvals', 1e10, 'MaxIter', 1e10);

% fit model to data with Nelder Mead algorithm
fitres   = nelder_mead_fit(fitfunc, x, y, var0, opt);
conf_int = confint(fitres); % estimate confidence intervals

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

%%                      Susceptibility Calculation
%%=========================================================================
pars.GFactor   = b2g(B0*1e-4, pars.MWFQ);
modScaling     = pars.B0MA*1e4 * 1e4/8; % scaling for pseudo-modulation
doubleIntAreas = modScaling * Bmw .* A;

Chi   = susceptebility_calc(doubleIntAreas, pars);
NSpin = spincounting(doubleIntAreas, pars);

%%                                Output
%%=========================================================================

% create output structure
argout.x        = x;
argout.y        = y;
argout.pars     = pars;

argout.fitres   = fitres;

argout.A        = A;
argout.B0       = B0;
argout.T2       = T2;
argout.Brms     = Brms;

argout.dA       = dA;
argout.dB0      = dB0;
argout.dT2      = dT2;
argout.dBrms    = dBrms;

argout.Chi      = Chi(1);
argout.NSpin    = NSpin(1);

end