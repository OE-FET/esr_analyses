function [argout] = PowerSatAnalysesLorentzFit(varargin)
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
%	argout  - structure containing the measurement data and fit results
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2019/05/06 12:58 $    $Revision: 1.1 $

close all

[x, y, pars] = load_spectrum_dialog(varargin);

if ~strcmp(pars.YTYP, 'IGD')
    error('The specified file is not a 2D data file.');
end

%%                         Calculate MW fields
%%=========================================================================
Bmw = get_mw_fields(pars);

%%                      Get starting points for fit
%%=========================================================================

% get initial parameters from slice fit
mid  = round(length(Bmw)/2);
fit = lorentz_fit(x, y(:,mid), 'deriv', 1);

FWHM = fit.FWHM;               % in Gauss

A0   = fit.a/(pars.B0MA*1e4 * 1e4/8 * Bmw(mid));
B0   = fit.x0;                  % in Gauss
T2   = 1/(gmratio * FWHM*1E-4); % in sec
T1   = 10*T2;                   % in sec

var0 = [A0 B0 T1 T2];           % starting points

%%                          Perform Lorentz fit
%%=========================================================================

% grid data for fitting algorithm
[X, Y]  = meshgrid(x, Bmw);
Z       = y;

% create fit function and options
fitfunc = @(var, x) abs(var(1))*ESRLorentzSimulation(x{1}, abs(var(2)), ...
    abs(var(3)), abs(var(4)), x{2}, pars.B0MA*1e4, 1);
opt = optimset('TolFun', 1e-9, 'TolX', 1e-9, 'PlotFcns', ...
    @optimplotfval, 'MaxFunEvals', 1e10, 'MaxIter', 1e10);

% fit model to data with Nelder Mead algorithm
fitres   = nelder_mead_fit(fitfunc, {X, Y}, Z, var0, opt);
conf_int = confint(fitres); % estimate confidence intervals

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
argout.T1       = T1;
argout.T2       = T2;

argout.dA       = dA;
argout.dB0      = dB0;
argout.dT1      = dT1;
argout.dT2      = dT2;

argout.Chi      = Chi(1);
argout.NSpin    = NSpin(1);

end

