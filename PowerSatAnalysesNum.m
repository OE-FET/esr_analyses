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
%   The magnetic susceptibilty is then determined from A0. Its accuracy
%   will therefore depend on the quality of the fit.
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

% make sure QValue and QValueErr are given
pars = get_par(pars, 'QValue');
pars = get_par(pars, 'QValueErr');

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

try
    g_factors = gfactor_determination(x, y, pars, 'plot', 'y');
catch
    fprintf('Could not determine g-factor. Using free electron value.\n')
    g_factors = ones(size(Bmw))*gfree;
end

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
fitres   = nelder_mead_fit(fitfunc, Bmw(1:end), doubleIntAreas(1:end), var0, opt);
se = standarderror(fitres); % estimate confidence intervals

A     = abs(fitres.coef(1));
T1T2  = abs(fitres.coef(2));

dT1T2 = abs(se(2));

% refine and get fit errors from standard matlab fit
if gm^2*T1T2 < 1
    fitfunc = @(A, x) A * x;

    ft = fittype(fitfunc, 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions('Method', 'NonlinearLeastSquares', ...
                      'Algorithm', 'Levenberg-Marquardt',...
                      'Display', 'Off', ...
                      'Robust', 'LAR', ...
                      'StartPoint', A);

    fitres = fit(Bmw(1:end), doubleIntAreas(1:end), ft, opts );

    A     = fitres.A;
    dA    = diff(confint(fitres))/2;
else
    fitfunc = @(A, gm2T1T2, x) A * x ./sqrt(1+x.^2*gm2T1T2);

    ft = fittype(fitfunc, 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions('Method', 'NonlinearLeastSquares', ...
                      'Algorithm', 'Levenberg-Marquardt',...
                      'Display', 'Off', ...
                      'Robust', 'LAR', ...
                      'StartPoint', [A, gm^2*T1T2]);

    fitres = fit(Bmw(1:end), doubleIntAreas(1:end), ft, opts );

    A     = fitres.A;
    T1T2  = fitres.gm2T1T2/gm^2;
    
    ci    = confint(fitres);
    dA    = diff(ci(:,1))/2;
    dT1T2    = diff(ci(:,2))/2;
    
end
%%                           Plot results
%%=========================================================================

figure();
h = plot(fitres);

xlabel(h.Parent, 'Microwave field [T]')
ylabel(h.Parent, 'ESR signal area [a.u.]')

hold on; plot(Bmw, doubleIntAreas, 'ko', 'DisplayName', 'data');


%%                          Spin counting
%%=========================================================================
areaDI = Bmw.*A;
areaDIerror = Bmw.*dA;

[Chi, dChi]     = susceptebility_calc(areaDI, pars, 'dA', areaDIerror);
[NSpin, dNSpin] = spincounting(areaDI, pars, 'dA', areaDIerror);

%%                              Output
%%=========================================================================

T = str2double(strtrim(regexprep(pars.Temperature,'K','')));

argout = struct('x', x, 'y', y, 'pars', pars, 'fitres', fitres, 'T', T, ...
                'T1T2', T1T2, 'dT1T2', dT1T2, ...
                'Chi', Chi(1), 'dChi', dChi(1), ...
                'NSpin', NSpin(1), 'dNSpin', dNSpin(1));

end
