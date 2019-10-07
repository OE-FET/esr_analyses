function [out_struct, out_table] = PowerSatAnalysesLorentzFit(varargin)
%POWERSATANALYSESLORENTZFIT Analyses of CW-ESR power saturation measurements
%
%   Performs a 2D fit of a cw-EPR power saturation spectrum to determine
%   spin lifetimes, the magnetic susceptibility and the number of spins in
%   the sample. This program assumes a Lorentzian resonance line, i.e. a
%   single spin-ensemble with lifetimes T1 and T2.
%
%   Field modulation and the 1st harmonic detection of the cw-EPR signal,
%   together with possible distortions from overmodulation, are explicitly
%   accounted for. The MW field distribution in the cavity is taken from
%   the Bruker DSC file and is used when calculating the MW field amplitude
%   over the sample volume.
%
%   INPUT(S):
%   POWERSATANALYSESLORENTZFIT()      - opens gui to select data
%   ...('/path/to/file')              - path to signal data
%   ...('signal_path', 'bg_path')     - path to signal data, path to
%                                       background data
%   ...(x, y, pars)                   - uses given data directly
%
%   OUTPUT(S):
%	out_struct  - structure containing the measurement data and fit results
%   out_table   - fit results in table format
%

import esr_analyses.*
import esr_analyses.utils.*

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2019/05/06 12:58 $    $Revision: 1.1 $

close all

[x, y, pars] = load_spectrum_dialog(varargin);

assert_powersat_exp(pars)

%%                         Calculate MW fields
%%=========================================================================
Bmw = get_mw_fields(pars);

%%                      Get starting points for fit
%%=========================================================================

% perform slice fit of center spectrum
mid  = round(length(Bmw)/2);
slice_fit  = lorentz_fit(x, y(:,mid), 'deriv', 1);

% perform numerical double-integration to estimate T1*T2
DI = double_int_num(x, y, 'baseline', 'n');
ft = fittype('A * x /sqrt(1+1e7*gmSquaredT1T2*x^2)');
pwrst_fit = fit(Bmw, DI, ft, 'StartPoint', [slice_fit.a, 1], 'Lower', [0, 0]);

FWHM = slice_fit.FWHM;                                 % in Gauss
A0   = slice_fit.a/(pars.B0MA*1e4 * 1e4/8 * Bmw(mid)); % see 'modScaling'
B0   = slice_fit.x0;                                   % in Gauss
T1T2 = 1e7*pwrst_fit.gmSquaredT1T2 / gmratio^2;            % in sec^2
T2   = 2/(gmratio * FWHM*1E-4);                        % in sec
T1   = T1T2/T2;                                        % in sec

var0 = [A0 B0 T1 T2];                                  % starting points

%%                          Perform Lorentz fit
%%=========================================================================

% grid data for fitting algorithm
[X, Y]  = meshgrid(x, Bmw);
Z       = y;

% create fit function and options
fitfunc = @(var, x) abs(var(1))*esr_lorentz_simulation(x{1}, abs(var(2)), ...
    abs(var(3)), abs(var(4)), x{2}, pars.B0MA*1e4, 1);
opt = optimset('TolFun', 1e-5, 'TolX', 1e-5, 'PlotFcns', ...
    @optimplotfval, 'MaxFunEvals', 1e10, 'MaxIter', 1e10);

% fit model to data with Nelder Mead algorithm
fitres   = nelder_mead_fit(fitfunc, {X, Y}, Z, var0, opt);
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

%%                      Susceptibility Calculation
%%=========================================================================
pars.GFactor   = b2g(B0*1e-4, pars.MWFQ);
modScaling     = pars.B0MA*1e4 * 1e4/8; % scaling for pseudo-modulation

areaDI = modScaling * Bmw .* A;
areaDIerror = modScaling * Bmw .* dA;

[Chi, dChi]     = susceptebility_calc(areaDI, pars, 'dA', areaDIerror);
[NSpin, dNSpin] = spincounting(areaDI, pars, 'dA', areaDIerror);

%%                                Output
%%=========================================================================

% create output structure

out_struct = struct(...
    'x', x, 'y', y, 'pars', pars, 'fitres', fitres, ...
    'B0', B0, 'dB0', dB0, 'T1', T1, 'dT1', dT1, 'T2', T2, 'dT2', dT2, ...
    'Chi', Chi(1), 'dChi', dChi(1), 'NSpin', NSpin(1), 'dNSpin', dNSpin(1));

out_table = struct2table(rmfield(out_struct, {'x', 'y', 'pars', 'fitres'}));

disp(out_table);

end
