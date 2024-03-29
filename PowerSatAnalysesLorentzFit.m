function [out_struct, out_table] = PowerSatAnalysesLorentzFit(varargin)
%POWERSATANALYSESLORENTZFIT Analyses of CW-ESR power saturation measurements
%
%   Performs a 2D fit of a cw-EPR power saturation spectrum to determine
%   spin lifetimes, the magnetic susceptibility and the number of spins in
%   the sample. This program assumes a Lorentzian resonance line, i.e. a
%   single spin-ensemble with lifetimes T1 and T2.
%
%   Field modulation and the 1st harmonic detection of the cw-EPR signal,
%   together with possible distortions from over-modulation, are explicitly
%   accounted for. The MW field distribution in the cavity is taken from
%   the Bruker DSC file and is used when calculating the MW field amplitude
%   over the sample volume.
%
%
%   SYNTAX:
%	POWERSATANALYSESLORENTZFIT()      - opens GUI for file selection
%	POWERSATANALYSESLORENTZFIT(dset)  - uses data given by dset
%	...(x,o,pars)                     - uses data given by [x,o,pars]
%	...('sigPath')                    - reads data from file
%	...('sigPath', 'bgPath')          - reads data and background from file
%
%   KEYWORD INPUT:
%   plot        - if true, plot data and best fit at each iteration
%
%   OUTPUT:
%	out_struct  - structure containing the measurement data and fit results
%   out_table   - fit results in table format
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$

import esr_analyses.*
import esr_analyses.utils.*

[plotting, varargin] = get_kwarg(varargin, 'plot', true);

dset = load_spectrum_dialog(varargin{:});
assert_powersat_exp(dset);
[x,o,pars] = dset_to_tuple(dset);

yes = input('Would you like to perform a baseline correction? y/[n] ', 's');
if strcmp(yes, 'y')
    o = baseline_corr(x, o);
end

%%                         Calculate MW fields
%%=========================================================================
Bmw = get_mw_fields(pars);

%%                      Get starting points for fit
%%=========================================================================

% perform numerical double-integration to estimate T1*T2
DI = double_int_num(x, o, 'baseline', false);
scaling = 1e4;
ft = fittype(sprintf('A * 1e9 * x /sqrt(1 + %e * gmSquaredT1T2 * x^2)', scaling));
pwrst_fit = fit(Bmw, DI, ft, 'StartPoint', [1, 1], 'Lower', [0, 0]);

% perform slice fit of non-saturated spectrum
index = sum(Bmw.^2 * scaling*pwrst_fit.gmSquaredT1T2 < 0.5);
index = max(1, index);
slice_fit  = lorentz_fit(x, o(:,index), 'deriv', 1);

FWHM = slice_fit.FWHM;                                   % in Gauss
A0   = slice_fit.a/(pars.B0MA*1e4 * 1e4/8 * Bmw(index)); % see 'modScaling'
B0   = slice_fit.x0;                                     % in Gauss
T1T2 = scaling*pwrst_fit.gmSquaredT1T2 / gmratio^2;      % in sec^2
T2   = 2/(gmratio * FWHM*1E-4);                          % in sec
T1   = T1T2/T2;                                          % in sec

var0 = [A0 B0 T1 T2];                                    % starting points

%%                          Perform Lorentz fit
%%=========================================================================

% grid data for fitting algorithm
[X, Y]  = meshgrid(x, Bmw);
O       = o;

% create fit function and options
fitfunc = @(var, x) abs(var(1))*esr_lorentz_simulation(x{1}, abs(var(2)), ...
    abs(var(3)), abs(var(4)), x{2}, pars.B0MA*1e4, 1);
opt = optimset('TolFun', 1e-5, 'TolX', 1e-5, 'MaxFunEvals', 1e10, 'MaxIter', 1e10);

% fit model to data with Nelder Mead algorithm
fitres   = nelder_mead_fit(fitfunc, {X, Y}, O, var0, opt, 'plot', plotting);
conf_int = standarderror(fitres, 'quick'); % rough estimate of confidence intervals

A     = abs(fitres.coef(1));
B0    = abs(fitres.coef(2));
T1    = abs(fitres.coef(3));
T2    = abs(fitres.coef(4));

dA    = abs(conf_int(1));
dB0   = abs(conf_int(2));
dT1   = abs(conf_int(3));
dT2   = abs(conf_int(4));

%%                           Plot results
%%=========================================================================

h = plot(fitres);

xlabel(h{1}(1).Parent, 'Magnetic field [G]')
ylabel(h{1}(1).Parent, 'Microwave field [T]')
zlabel(h{1}(1).Parent, 'ESR signal [a.u.]')

xlabel(h{3}(1).Parent, 'Magnetic field [G]')
ylabel(h{3}(1).Parent, 'ESR signal [a.u.]')

title(h{1}(1).Parent, pars.TITL, 'interpreter', 'none')
title(h{3}(1).Parent, pars.TITL, 'interpreter', 'none')

%%                      Susceptibility Calculation
%%=========================================================================
pars.GFactor    = b2g(B0*1e-4, pars.MWFQ);
pars.GFactorErr = dB0 * b2g(B0*1e-4, pars.MWFQ) / B0;
modScaling      = pars.B0MA*1e4 * 1e4/8; % scaling for pseudo-modulation

areaDI = modScaling * Bmw .* A;
areaDIerror = modScaling * Bmw .* dA;

[Chi, dChi]     = susceptibility_calc(areaDI, pars, 'dA', areaDIerror);
[NSpin, dNSpin] = spincounting(areaDI, pars, 'dA', areaDIerror);

%%                                Output
%%=========================================================================

% create output structure

out_struct = struct(...
    'x', x, 'o', o, 'pars', pars, 'fitres', fitres, ...
    'B0', B0, 'dB0', dB0, 'g', pars.GFactor, 'dg', pars.GFactorErr, ...
    'T1', T1, 'dT1', dT1, 'T2', T2, 'dT2', dT2, ...
    'Chi', Chi(1), 'dChi', dChi(1), 'NSpin', NSpin(1), 'dNSpin', dNSpin(1));

out_table = struct2table(rmfield(out_struct, {'x', 'o', 'pars', 'fitres'}));

disp(out_table);

end
