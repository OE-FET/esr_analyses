function [argout] = PowerSatAnalysesNum(varargin)
%POWERSATANALYSESNUM Numercial analyses of ESR power saturation curves.
%
%   Numerically integrates a series of MW power dependent cw-EPR spectra to
%   determine the power saturation beheviour and the product of spin
%   lifetimes T1*T2. Allows for baseline correction and smoothing before
%   the integration.
%
% 	Spin lifetimes are calculated by fitting the integrated areas to the
% 	following equation:
%
%   A = A0 * B_mw / sqrt(1 + gmratio^2 * B_mw^2 * T1 * T2)
%
% 	WARNING: Using numerical integration may underestimate the tails of EPR
% 	spectra. This can lead to significant errors for long-tailed resonance
% 	shapes, such as Lorentzians, if the SNR ratio is small or the
% 	measurement range is less than 5 times the peak-to-peak linewidth.
%
%	INPUTS:
%	POWERSATANALYSESNUM()         - opens GUI for file selection
%	POWERSATANALYSESNUM(x,y,pars) - uses data given by (x,y,pars)
%	...('sigPath')                - reads data from file
%	...('sigPath', 'bgPath')      - reads data and background from file
%
%	OUTPUT:
%	argout  - structure containing the measurement data and fit results
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2019/05/06 12:58 $    $Revision: 1.1 $

close all

[x, y, pars] = load_spectrum_dialog(varargin);

if ~strcmp(pars.YTYP, 'IGD')
    error('The specified file is not a 2D data file.');
end

g_factors = gfactor_determination(x, y, pars, 'plot', 'y');

%%                         Calculate MW fields
%%=========================================================================
Bmw = get_mw_fields(pars);

%%                      Double integration of spectra
%%=========================================================================

baseline = input('Perform base-line correction individually or as batch? i/[b]?', 's');
if baseline == 'i'
    doubleIntAreas = zeros(size(y, 2), 1);
    for i = 1:size(y, 2)
        doubleIntAreas(i) = double_int_num(x, y(:,i), 'y');
    end
else
    doubleIntAreas = double_int_num(x, y, 'y');
end

%%                      Fit power saturation curve
%%=========================================================================

% determine starting values
pars.GFactor = g_factors(ceil(end/2));
gm = pars.GFactor * bmagn / hbar;

A0 = mean(1.1*doubleIntAreas./Bmw);
T0 = mean( (A0^2*Bmw(end).^2 - doubleIntAreas(end).^2*gm^2)./ ...
    (Bmw(end).^2.*doubleIntAreas(end).^2*gm^4) );

var0 = [A0 T0];

% set up fit functin and options
fitfunc = @(var, x) abs(var(1)) * x ./sqrt(1+gm^2*x.^2*abs(var(2)));
opt = optimset('TolFun', 1e-9,'TolX', 1e-9,'PlotFcns', @optimplotfval, ...
    'MaxFunEvals', 1e9, 'MaxIter', 1e9);

% fit model to data with Nelder Mead algorithm
fitres   = nelder_mead_fit(fitfunc, Bmw, doubleIntAreas, var0, opt);
conf_int = confint(fitres); % estimate confidence intervals

A     = abs(fitres.coef(1));
T1T2  = abs(fitres.coef(2));

dA    = abs(conf_int(1));
dT1T2 = abs(conf_int(2));

%%                           Plot results
%%=========================================================================

h = plot(fitres);
set(h{1}, 'Marker', 'o');

xlabel(h{1}(1).Parent, 'Microwave field [T]')
ylabel(h{1}(1).Parent, 'ESR signal area [a.u.]')


%%                          Spin counting
%%=========================================================================
doubleIntAreasCalc = Bmw .* A;
Chi   = susceptebility_calc(doubleIntAreasCalc, pars);
NSpin = spincounting(doubleIntAreasCalc, pars);

%%                              Output
%%=========================================================================

argout.x        = x;
argout.y        = y;
argout.pars     = pars;

argout.fitres   = fitres;

argout.A        = A;
argout.T1T2     = T1T2;

argout.dA       = dA;
argout.dT1T2    = dT1T2;

argout.Chi      = Chi(1);
argout.NSpin    = NSpin(1);


end