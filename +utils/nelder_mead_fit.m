function [fitres, output] = nelder_mead_fit(fitfunc, x, y, coef0, opt, varargin)
%NELDER_MEAD_FIT performs derivative-free non-linear least squares fitting.
%
%   SYNTAX:
%   [fitobject, gof, output] = nelder_mead_fit(fitfunc, x, y, coef0, opt)
%
%   INPUT(S):
%   fitfunc - A function handle to a fit function. The first argument must
%             contain the fit coefficients, the second argument must
%             contain a single independent variable.
%   x       - Data to fit.
%   y       - Data to fit.
%   coef0   - Vector with starting points for fit coefficients.
%   opt     - Fit options generated by optimset.
%
%   OUTPUTS(S):
%   fitres  - Object holding the fit function, the best-fit coefficients,
%             the starting points and the fit data. Call
%             standarderror(fitres) to get confidence intervals for fit
%             coefficients.
%   output  - Structure containing information associated with the
%             fitting algorithm.
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2018/07/05 12:58 $    $Revision: 0.1 $

import esr_analyses.*
import esr_analyses.utils.*

plotting = get_kwarg(varargin, 'plot', true);

% check correct dimensions of fit function
y0 = fitfunc(coef0, x);
if size(y0) ~= size(y)
    error('Dimensions of fit function do not agree with dimensions of data to fit.');
end

% define sum of squares to minimize
sumofsquares = @(coef) sum(sum( abs(fitfunc(coef, x) - y).^2  ));

% define our own plot function to call at each iteration

function stop = plot_func(coef, optimValues, ~)
    stop = false;
    x_slice = x{1}(1,:)';
    [~, yoffsets] = stackplot(x_slice, y, 'style', 'k.');
    hold on;
    stackplot(gca, x_slice, fitfunc(coef, x), 'yoffsets', yoffsets, 'style', 'r');
    set(get(gca,'Title'), 'String', sprintf('Current fit error: %g', optimValues.fval));
    hold off;
end

if plotting
    opt_plot_func = optimset('PlotFcns', @plot_func);
    opt = optimset(opt, opt_plot_func);
else
    opt_plot_func = optimset('Display', 'iter');
    opt = optimset(opt, opt_plot_func);
end

% run Nelder-Mead fit
[best_coef, sumofsquares_error, exitflag, output] = fminsearch(sumofsquares, coef0, opt);

% create fitobject structure
fitres = nelder_mead_fitobject(fitfunc, x, y, coef0, best_coef, sumofsquares_error);

% create output structure
if nargout > 1
    output.exitflag = exitflag;
end

end
