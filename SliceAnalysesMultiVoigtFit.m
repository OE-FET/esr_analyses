function [out_struct, out_table] = SliceAnalysesMultiVoigtFit(varargin)
%SLICEANALYSESMUTLIVOIGTFIT Analyses of CW-ESR measurements
%
%   This function is equivelent to SliceAnalysesVoigtFit but fits multiple
%   Voigtians instead of a single one. The number of Voigtians can be given
%   as a keyword argument 'N' and defaults to two.
%
%   Starting points may be given as a keyword argument 'var0'. var0 must be
%   a Nx5 matrix with each row containing the starting points for one
%   Voigtian in the order: [A0 B0 T2 FWHM_gauss]. If NaN is given for any
%   starting point, it will be substituted with a reasonable best-guess
%   value.
%
%   Convergance may be bad. Take care to choose good starting points for
%   the fit!
%
%   INPUT(S):
%   SLICEANALYSESMUTLIVOIGTFIT(..., 'N', 3)       - fits three Voigtians
%   SLICEANALYSESMUTLIVOIGTFIT(..., 'var0', var0) - gives starting points
%
%   OUTPUT(S):
%	out_struct  - structure containing the measurement data and fit results
%   out_table   - fit results in table format
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2019/05/27 12:58 $    $Revision: 1.2 $

import esr_analyses.*
import esr_analyses.utils.*

close all

[N, varargin] = get_kwarg(varargin, 'N', 2);
[var0, varargin] = get_kwarg(varargin, 'var0', nan(N, 4));
if N ~= size(var0, 1)
    error('The number of starting points must match the number of Voigtians to fit.');
end

dset = load_spectrum_dialog(varargin{:});
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
    % perform slice fit of center spectrum
    fit = pseudo_voigt_fit(x, y, 'deriv', 1);

    FWHM_lorentz  = fit.FWHM_lorentz;             % in Gauss
    FWHM_gauss    = fit.FWHM_gauss;               % in Gauss

    A0   = fit.a/(pars.B0MA*1e4 * 1e4/8 * Bmw);   % in a.u.
    B0   = fit.x0;                                % in Gauss
    T2   = 1/(gmratio * FWHM_lorentz*1E-4);       % in sec

    auto_var0 = ones(N,1)*[A0 B0 T2 FWHM_gauss];  % starting points

    % replace NaN values in var0 with our best-guess starting points
    % keep all starting points provided by the user
    var0(isnan(var0)) = auto_var0(isnan(var0));
end

T1 = 1e-20; % in sec, assume not-saturated

%%                          Perform Voigt fit
%%=========================================================================

% create single fit function
func_single = @(v, x) abs(v(1))*esr_voigt_simulation(x, abs(v(2)), ...
    T1, abs(v(3)), abs(v(4)), Bmw, pars.B0MA*1e4, 1);

% expand to multiple peaks
multi_fit_func = @(v, x) to_multi(func_single, v, x);

% set fit options
opt = optimset('TolFun', 1e-9, 'TolX', 1e-9, 'PlotFcns', ...
    @optimplotfval, 'MaxFunEvals', 1e10, 'MaxIter', 1e10);

% fit model to data with Nelder Mead algorithm
fitres   = nelder_mead_fit(multi_fit_func, x, y', var0, opt);
conf_int = standarderror(fitres, 'quick')'; % estimate confidence intervals

A     = abs(fitres.coef(:,1));
B0    = abs(fitres.coef(:,2));
T2    = abs(fitres.coef(:,3));
Brms  = abs(fitres.coef(:,4));

dA    = full(abs(conf_int(1:N)));
dB0   = full(abs(conf_int(N+1:2*N)));
dT2   = full(abs(conf_int(2*N+1:3*N)));
dBrms = full(abs(conf_int(3*N+1:4*N)));

%%                           Plot results
%%=========================================================================

% use custom stackplot showing individual peaks instead of plot(fitres)
figure();hold on;
p = {};
[p{1}, yoffsets] = stackplot(x, y, 'style', 'k.');
p{2} = stackplot(x, multi_fit_func(fitres.coef, x), 'style', 'r', 'yoffsets', yoffsets);

legend_texts = {'Data', 'Fit'};
plot_handles = [p{1}(1), p{2}(1)];
linespecs = {'-.g', ':b', '--c', '-.m', ':y'};

for i=1:N
    c = fitres.coef(i,:);
    p{i+2} = stackplot(x, func_single(c, x), 'yoffsets', yoffsets, 'style', linespecs{i});
    legend_texts{i+2} = ['Peak ' num2str(i)];
    plot_handles(i+2) = p{i+2}(1);
end

legend(plot_handles, legend_texts);
xlabel('Magnetic field [G]')
ylabel('ESR signal [a.u.]')

axis tight;

%%                      Susceptibility Calculation
%%=========================================================================
pars.GFactor   = mean(b2g(B0*1e-4, pars.MWFQ));
modScaling     = pars.B0MA*1e4 * 1e4/8; % scaling for pseudo-modulation

Chi = zeros(size(A)); dChi = zeros(size(A));
NSpin = zeros(size(A)); dNSpin = zeros(size(A));

for i=1:length(A) % calculate for each peak
    areaDI = modScaling * Bmw .* A(i);
    areaDIerror = modScaling * Bmw .* dA(i);

    % get 'maximum' value, even though all values are equal...
    [Chi(i), dChi(i)]   = susceptibility_calc(areaDI, pars, 'dA', areaDIerror);
    [NSpin(i), dNSpin(i)] = spincounting(areaDI, pars, 'dA', areaDIerror);
end

%%                                Output
%%=========================================================================

% create output structure

out_struct = struct(...
    'x', x, 'y', y, 'pars', pars, 'fitres', fitres, ...
    'B0', B0.', 'dB0', dB0.', 'T2', T2.', 'dT2', dT2.',...
    'Brms', Brms.', 'dBrms', dBrms.', 'Chi', Chi.',...
    'dChi', dChi.','NSpin', NSpin.', 'dNSpin', dNSpin.');

out_table = struct2table(out_struct,'AsArray',true);

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
