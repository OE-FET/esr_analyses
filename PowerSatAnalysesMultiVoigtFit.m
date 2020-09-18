function [out_struct, out_table] = PowerSatAnalysesMultiVoigtFit(varargin)
%POWERSATANALYSESMUTLIVOIGTFIT Analyses of CW-ESR power saturation measurements
%
%   This function is equivelent to PowerSatAnalysesVoigtFit but fits
%   multiple Voigtians instead of a single one. The number of Voigtians can
%   be given as a keyword argument 'N' and defaults to one.
%
%   Starting points may be given as a keyword argument 'var0'. var0 must be
%   a Nx5 matrix with each row containing the starting points for one
%   Voigtian in the order: [A0 B0 T1 T2 FWHM_gauss]. If NaN is given for
%   any starting point, it will be substituted with a reasonable best-guess
%   value.
%
%   Convergance may be bad. Take care to choose good starting points for
%   the fit!
%
%   INPUT SYNTAX:
%	POWERSATANALYSESMUTLIVOIGTFIT()          - opens GUI for file selection
%	POWERSATANALYSESMUTLIVOIGTFIT(dset)      - uses data given by dset
%	POWERSATANALYSESMUTLIVOIGTFIT(x,o,pars)  - uses data given by [x,o,pars]
%	...('sigPath')                           - reads data from file
%	...('sigPath', 'bgPath')                 - reads data and background from file
%
%   KEYWORD INPUT(S):
%   N       - number of voigtians to fit
%   var0    - starting points
%   plot    - if true, plot data and best fit at each iteration
%   LB      - lower bound vector or array, must be the same size as x0
%
%             If no lower bounds exist for one of the variables, then
%             supply -inf for that variable.
%
%             If no lower bounds at all, then LB may be left empty.
%
%             Variables may be fixed in value by setting the corresponding
%             lower and upper bounds to exactly the same value.
%
%   UB      - upper bound vector or array, must be the same size as x0
%
%             If no upper bounds exist for one of the variables, then
%             supply +inf for that variable.
%
%             If no upper bounds at all, then LB may be left empty.
%
%             Variables may be fixed in value by setting the corresponding
%             lower and upper bounds to exactly the same value.
%
%   OUTPUT(S):
%	out_struct  - structure containing the measurement data and fit results
%   out_table   - fit results in table format
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2019/05/27 12:58 $    $Revision: 1.2 $

import esr_analyses.*
import esr_analyses.utils.*

[N, varargin] = get_kwarg(varargin, 'N', 2);
[var0, varargin] = get_kwarg(varargin, 'var0', nan(N, 5));
[plotting, varargin] = get_kwarg(varargin, 'plot', true);
[LB, varargin] = get_kwarg(varargin, 'LB', []);
[UB, varargin] = get_kwarg(varargin, 'UB', []);

if N ~= size(var0, 1)
    error('The number of starting points must match the number of Voigtians to fit.');
end

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

if any(isnan(var0), 'all')

    % perform numerical double integrtion to estimate T1*T2
    DI = double_int_num(x, y, 'baseline', false);
    scaling = 1e4;
    ft = fittype(sprintf('A * 1e9 * x /sqrt(1 + %e * gmSquaredT1T2 * x^2)', scaling));
    pwrst_fit = fit(Bmw, DI, ft, 'StartPoint', [1, 1], 'Lower', [0, 0]);
    
    % perform slice fit of non-saturated spectrum
    index = sum(Bmw.^2 * scaling*pwrst_fit.gmSquaredT1T2 < 0.5);
    index = max(1, index);
    slice_fit  = pseudo_voigt_fit(x, y(:,index), 'deriv', 1);

    FWHM_lorentz  = slice_fit.FWHM_lorentz;                  % in Gauss
    FWHM_gauss    = slice_fit.FWHM_gauss;                    % in Gauss

    A0   = slice_fit.a/(pars.B0MA*1e4 * 1e4/8 * Bmw(index)); % see 'modScaling'
    B0   = slice_fit.x0;                                     % in Gauss
    T1T2 = scaling*pwrst_fit.gmSquaredT1T2 / gmratio^2;      % in sec^2
    T2   = 2/(gmratio * FWHM_lorentz*1E-4);                  % in sec
    T1   = T1T2/T2;                                          % in sec

    auto_var0 = ones(N,1) * [A0/N B0 T1 T2 FWHM_gauss];      % starting points

    % wiggle all starting points apart from B0
    wiggle = 1 + (rand(size(auto_var0))-0.5)/100;
    auto_var0 = wiggle .* auto_var0;
    auto_var0(:,2) = [B0 B0];

    % replace NaN values in var0 with our best-guess starting points
    % keep all starting points provided by the user
    var0(isnan(var0)) = auto_var0(isnan(var0));
end

%%                          Perform Voigt fit
%%=========================================================================

% grid data for fitting algorithm
[X, Y]  = meshgrid(x, Bmw);
Z       = y;

% create single fit function
func_single = @(v, x) v(1)*esr_voigt_simulation(x{1}, v(2), v(3), v(4), v(5), x{2}, pars.B0MA*1e4, 1);

% expand to multiple peaks
multi_fit_func = @(v, x) to_multi(func_single, N, v, x);

% set fit options
opt = optimset('TolFun', 1e-9, 'TolX', 1e-9, 'MaxFunEvals', 1e10, 'MaxIter', 1e10);

% fit model to data with Nelder Mead algorithm
fitres = nelder_mead_fit(multi_fit_func, ...
    {X, Y}, Z, var0, opt,...
    'plot', plotting, 'LB', LB, 'UB', UB);
conf_int = standarderror(fitres, 'quick'); % estimate confidence intervals

A     = fitres.coef(:,1);
B0    = fitres.coef(:,2);
T1    = fitres.coef(:,3);
T2    = fitres.coef(:,4);
Brms  = fitres.coef(:,5);

dA    = full(abs(conf_int(1:N)));
dB0   = full(abs(conf_int(N+1:2*N)));
dT1   = full(abs(conf_int(2*N+1:3*N)));
dT2   = full(abs(conf_int(3*N+1:4*N)));
dBrms = full(abs(conf_int(4*N+1:5*N)));

%%                           Plot results
%%=========================================================================

% use custom stackplot showing individual peaks instead of plot(fitres)
figure();hold on;
p = {};
[p{1}, yoffsets] = stackplot(x, y, 'style', 'k.');
p{2} = stackplot(x, multi_fit_func(fitres.coef, {X, Y}), 'style', 'r', 'yoffsets', yoffsets);

legend_texts = {'Data', 'Fit'};
plot_handles = [p{1}(1), p{2}(1)];
linespecs = {'-.g', ':b', '--c', '-.m', ':y'};

for i=1:N
    c = fitres.coef(i,:);
    p{i+2} = stackplot(x, func_single(c, {X, Y}), 'yoffsets', yoffsets, 'style', linespecs{i});
    legend_texts{i+2} = ['Peak ' num2str(i)];
    plot_handles(i+2) = p{i+2}(1);
end

legend(plot_handles, legend_texts);
xlabel('Magnetic field [G]')
ylabel('ESR signal [a.u.]')
title(pars.TITL, 'interpreter', 'none')

axis tight;

%%                      Susceptibility Calculation
%%=========================================================================
g_factors     = b2g(B0*1e-4, pars.MWFQ);
g_factor_errs = dB0 .* g_factors ./ B0;
modScaling    = pars.B0MA*1e4 * 1e4/8; % scaling for pseudo-modulation

Chi = zeros(size(A)); dChi = zeros(size(A));
NSpin = zeros(size(A)); dNSpin = zeros(size(A));

for i=1:length(A) % calculate for each peak
    areaDI = modScaling * Bmw .* A(i);
    areaDIerror = modScaling * Bmw .* dA(i);

    % get 'maximum' value, even though all values are equal...
    pars.GFactor = g_factors(i); pars.GFactorErr = g_factor_errs(i);
    [Chi(i), dChi(i)]   = max(susceptibility_calc(areaDI, pars, 'dA', areaDIerror));
    [NSpin(i), dNSpin(i)] = max(spincounting(areaDI, pars, 'dA', areaDIerror));
end

%%                                Output
%%=========================================================================

% create output structure

out_struct = struct(...
    'x', x, 'y', y, 'pars', pars, 'fitres', fitres, ...
    'B0', B0.', 'dB0', dB0.', 'g', g_factors.', 'dg', g_factor_errs.', ...
    'T1', T1.', 'dT1', dT1.', 'T2', T2.', 'dT2', dT2.', ...
    'Brms', Brms.', 'dBrms', dBrms.', 'Chi', Chi.',...
    'dChi', dChi.','NSpin', NSpin.', 'dNSpin', dNSpin.');

out_table = struct2table(out_struct,'AsArray',true);

end

function y = to_multi(func_single, N, variables, x)

variables = reshape(variables, [N numel(variables)/N]);

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
