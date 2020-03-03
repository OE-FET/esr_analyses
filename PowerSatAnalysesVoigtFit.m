function [out_struct, out_table] = PowerSatAnalysesVoigtFit(varargin)
%POWERSATANALYSESVOIGTFIT Analyses of CW-ESR power saturation measurements
%
%   Performs a 2D fit of a cw-EPR power saturation spectrum to determine
%   spin lifetimes, the magnetic susceptibility and the number of spins in
%   the sample. This program assumes a Voigtian resonance line which is the
%   convolution of multiple Lorentzian spin-ensembles with identifical
%   lifetimes T1 and T2 with a Gaussian distribution of resonance fields.
%
%   The effects of field modulation and the 1st harmonic detection of the
%   cw-EPR signal, together with possible distortions from overmodulation
%   are explicitly included. The MW field distribution in the cavity is
%   taken from the Bruker DSC file and is used when calculating the MW
%   field amplitude over the sample volume.
%
%   INPUT SYNTAX:
%	POWERSATANALYSESVOIGTFIT()      - opens GUI for file selection
%	POWERSATANALYSESVOIGTFIT(dset)  - uses data given by (x,y,pars)
%	...('sigPath')                  - reads data from file
%	...('sigPath', 'bgPath')        - reads data and background from file
%
%   OUTPUT(S):
%	out_struct  - structure containing the measurement data and fit results
%   out_table   - fit results in table format
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$

import esr_analyses.*
import esr_analyses.utils.*

close all

dset = load_spectrum_dialog(varargin{:});
assert_powersat_exp(dset);
[x,y,pars] = dset_to_tuple(dset);

yes = input('Would you like to perform a baseline correction? y/[n] ', 's');
if strcmp(yes, 'y')
    y = baseline_corr(x, y);
end

%%                         Calculate MW fields
%%=========================================================================
Bmw = get_mw_fields(pars);

%%                      Get starting points for fit
%%=========================================================================

% perform slice fit of center spectrum
mid  = round(length(Bmw)/2);
slice_fit  = pseudo_voigt_fit(x, y(:,mid), 'deriv', 1);

% Numerically double-integrate and then fit the 2D spectrum to estimate T1*T2
DI = double_int_num(x, y, 'baseline', false);
norm = 1e9;
ft = fittype('A * 1e9 * x /sqrt(1+1e9*gmSquaredT1T2*x^2)');
pwrst_fit = fit(Bmw, DI, ft, 'StartPoint', [slice_fit.a, 1], 'Lower', [0, 0]);

FWHM_lorentz  = slice_fit.FWHM_lorentz;                % in Gauss
FWHM_gauss    = slice_fit.FWHM_gauss;                  % in Gauss

A0   = slice_fit.a/(pars.B0MA*1e4 * 1e4/8 * Bmw(mid)); % see 'modScaling'
B0   = slice_fit.x0;                                   % in Gauss
T1T2 = norm*pwrst_fit.gmSquaredT1T2 / gmratio^2;       % in sec^2
T2   = 2/(gmratio * FWHM_lorentz*1E-4);                % in sec
T1   = T1T2/T2;                                        % in sec
% RLC 21/09/19 Note: Bad fits are often attributable to inaccurate
% estimates of T1. If your fits are off then change the order-of-magnitude
% for this estimation parameter.
var0 = [A0 B0 T1 T2 FWHM_gauss];                       % starting points

%%                          Perform Voigt fit
%%=========================================================================

% grid data for fitting algorithm
[X, Y]  = meshgrid(x, Bmw);
Z       = y;

% create fit function and options
fitfunc = @(var, x) abs(var(1))*esr_voigt_simulation(x{1}, abs(var(2)), ...
    abs(var(3)), abs(var(4)), abs(var(5)), x{2}, pars.B0MA*1e4, 1);
opt = optimset('TolFun', 1e-9, 'TolX', 1e-9, 'PlotFcns', ...
    @optimplotfval, 'MaxFunEvals', 1e10, 'MaxIter', 1e10);

% fit model to data with Nelder Mead algorithm
fitres   = nelder_mead_fit(fitfunc, {X, Y}, Z, var0, opt);
conf_int = standarderror(fitres, 'accurate'); % estimate confidence intervals

A     = abs(fitres.coef(1));
B0    = abs(fitres.coef(2));
T1    = abs(fitres.coef(3));
T2    = abs(fitres.coef(4));
Brms  = abs(fitres.coef(5));

dA    = abs(conf_int(1));
dB0   = abs(conf_int(2));
dT1   = abs(conf_int(3));
dT2   = abs(conf_int(4));
dBrms = abs(conf_int(5));

%%                           Plot results
%%=========================================================================

h = plot(fitres);

xlabel(h{1}(1).Parent, 'Magnetic field [G]')
ylabel(h{1}(1).Parent, 'Microwave field [T]')
zlabel(h{1}(1).Parent, 'ESR signal [a.u.]')

xlabel(h{3}(1).Parent, 'Magnetic field [G]')
ylabel(h{3}(1).Parent, 'ESR signal [a.u.]')

%%                      Susceptibility Calculation
%%=========================================================================
pars.GFactor   = b2g(B0*1e-4, pars.MWFQ);
modScaling     = pars.B0MA*1e4 * 1e4/8; % scaling for pseudo-modulation

areaDI = modScaling * Bmw .* A;
areaDIerror = modScaling * Bmw .* dA;

[Chi, dChi]     = susceptibility_calc(areaDI, pars, 'dA', areaDIerror);
[NSpin, dNSpin] = spincounting(areaDI, pars, 'dA', areaDIerror);

%%                                Output
%%=========================================================================

% create output structure

out_struct = struct(...
    'x', x, 'y', y, 'pars', pars, 'fitres', fitres, ...
    'B0', B0, 'dB0', dB0, 'T1', T1, 'dT1', dT1, 'T2', T2, 'dT2', dT2, ...
    'Brms', Brms, 'dBrms', dBrms, 'Chi', Chi(1), 'dChi', dChi(1), ...
    'NSpin', NSpin(1), 'dNSpin', dNSpin(1));

out_table = struct2table(rmfield(out_struct, {'x', 'y', 'pars', 'fitres'}));

disp(out_table);

end
