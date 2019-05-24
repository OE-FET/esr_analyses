function [argout] = PowerSatAnalysesMultiVoigtFit(varargin)
%POWERSATANALYSESMUTLIVOIGTFIT Analyses of CW-ESR power saturation measurements
%
%   This function is equivelent to PowerSatAnalysesVoigtFit but fits
%   multiple Voigtians instead of a single one. The number of Voigtians can
%   be given as a keyword argument 'N' and defaults to two. Starting points
%   may be given as a keyword argument 'var0'. They must be a Nx5 matrix
%   with each row containing the starting points for one Voigtian in the
%   order: [A0 B0 T1 T2 FWHM_gauss].
%
%   Convergance may be bad. Take care to choose good starting points for
%   the fit!
%
%   INPUT(S):
%   POWERSATANALYSESVOIGTFIT(..., 'N', 3)       - fits three Voigtians
%   POWERSATANALYSESVOIGTFIT(..., 'var0', var0) - gives starting points
%
%   OUTPUT(S):
%	argout  - structure containing the measurement data and fit results
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2019/05/06 12:58 $    $Revision: 1.1 $

close all

if nargin > 0
    [N, varargin] = get_varargin(varargin, 'N', 2);
    [var0, varargin] = get_varargin(varargin, 'var0', []);
else
    N = 2;
    var0 = [];
end

[x, y, pars] = load_spectrum_dialog(varargin);

if ~strcmp(pars.YTYP, 'IGD')
    error('The specified file is not a 2D data file.');
end

%%                         Calculate MW fields
%%=========================================================================
Bmw = get_mw_fields(pars);

%%                      Get starting points for fit
%%=========================================================================

if isempty(var0)
    % perform slice fit of center spectrum
    mid  = round(length(Bmw)/2);
    slice_fit  = pseudo_voigt_fit(x, y(:,mid), 'deriv', 1);

    % perform numerical double integrtion to estimate T1*T2
    DI = double_int_num(x, y, 'baseline', 'n');
    ft = fittype('A * x /sqrt(1+gmSquaredT1T2*x^2)');
    pwrst_fit = fit(Bmw, DI, ft, 'StartPoint', [slice_fit.a, 1e7], 'Lower', [0, 0]);

    FWHM_lorentz  = slice_fit.FWHM_lorentz;                % in Gauss
    FWHM_gauss    = slice_fit.FWHM_gauss;                  % in Gauss

    A0   = slice_fit.a/(pars.B0MA*1e4 * 1e4/8 * Bmw(mid)); % see 'modScaling'
    B0   = slice_fit.x0;                                   % in Gauss
    T1T2 = pwrst_fit.gmSquaredT1T2 / gmratio^2;            % in sec^2
    T2   = 2/(gmratio * FWHM_lorentz*1E-4);                % in sec
    T1   = T1T2/T2;                                        % in sec

    var0 = ones(N,1)*[A0/2 B0 T1 T2 FWHM_gauss];             % starting points
end

%%                          Perform Voigt fit
%%=========================================================================

% grid data for fitting algorithm
[X, Y]  = meshgrid(x, Bmw);
Z       = y;

% create single fit function
func_single = @(v, x) abs(v(1))*ESRVoigtSimulation(x{1}, abs(v(2)), ...
    abs(v(3)), abs(v(4)), abs(v(5)), x{2}, pars.B0MA*1e4, 1);

% expand to multiple peaks
multi_fit_func = @(v, x) to_multi(func_single, v, x);

% set fit options
opt = optimset('TolFun', 1e-9, 'TolX', 1e-9, 'PlotFcns', ...
    @optimplotfval, 'MaxFunEvals', 1e10, 'MaxIter', 1e10);

% fit model to data with Nelder Mead algorithm
fitres   = nelder_mead_fit(multi_fit_func, {X, Y}, Z, var0, opt);
conf_int = confint(fitres, 'quick')'; % estimate confidence intervals

A     = abs(fitres.coef(:,1));
B0    = abs(fitres.coef(:,2));
T1    = abs(fitres.coef(:,3));
T2    = abs(fitres.coef(:,4));
Brms  = abs(fitres.coef(:,5));

dA    = full(abs(conf_int(1:N)));
dB0   = full(abs(conf_int(N+1:2*N)));
dT1   = full(abs(conf_int(2*N+1:3*N)));
dT2   = full(abs(conf_int(3*N+1:4*N)));
dBrms = full(abs(conf_int(4*N+1:5*N)));

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
pars.GFactor   = mean(b2g(B0*1e-4, pars.MWFQ));
modScaling     = pars.B0MA*1e4 * 1e4/8; % scaling for pseudo-modulation

Chi = zeros(size(A));
NSpin = zeros(size(A));

for i=1:length(A)
    doubleIntAreas = modScaling * Bmw .* A(i); % use sum of all peak areas!
    
    % get 'maximum' value, even though all values are equal...
    Chi(i)   = max(susceptebility_calc(doubleIntAreas, pars));
    NSpin(i) = max(spincounting(doubleIntAreas, pars));
end

%%                                Output
%%=========================================================================

% create output structure

argout = struct(...
    'x', x, 'y', y, 'pars', pars, 'fitres', fitres, ...
    'A', A, 'B0', B0, 'T1', T1, 'T2', T2, 'Brms', Brms, ...
    'dA', dA, 'dB0', dB0, 'dT1', dT1, 'dT2', dT2, 'dBrms', dBrms, ...
    'Chi', Chi, 'NSpin', NSpin);

end

function y = to_multi(func_single, variables, x)

[N, ~] = size(variables);

result1 = func_single(variables(1,:), x);
shape = num2cell(size(result1));

results = zeros(N, shape{:});
results(1,:,:) = result1;

for i=2:N
    results(i,:,:) = func_single(variables(i,:), x);
end

y = sum(results, 1);
y = squeeze(y);

end