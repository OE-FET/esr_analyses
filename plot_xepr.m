function argout = plot_xepr(varargin)
%PLOT_XEPR Plots Xepr data.
%
%   In addition to PLOT2D, PLOT2D_XEPR will create axis labels and an
%   appropriate legend automatically from pars.
%
% 	PLOT_XEPR(x, y, pars)
% 	PLOT_XEPR(x, y, pars, 'style', style, ...)
%   PLOT_XEPR(ax, ...)
%
%   INPUT(S):
%   ax         - Axis handle for plot. If not given, the data is plotted in
%                the current axis, as returned by gca.
%   x, y, pars - Xepr data set
%
% 	OUTPUT(S):
%   phandle - image handle
%

import esr_analyses.*
import esr_analyses.utils.*

%% Input Analyses

if ishghandle(varargin{1}, 'axes')
    ax = varargin{1};
    x  = varargin{2};
    y  = varargin{3};
    pars  = varargin{4};
else
    ax = gca;
    x  = varargin{1};
    y  = varargin{2};
    pars  = varargin{3};
end

%% Data preparation
fig_titel = pars.TITL;
x_label = sprintf('%s [%s]', pars.XNAM, pars.XUNI);
y_label = 'ESR signal [a.u.]';

%% Plot
h = plot(ax, x, y);
xlabel(ax, x_label, 'Interpreter', 'none');
ylabel(ax, y_label, 'Interpreter', 'none');
title(ax, fig_titel, 'Interpreter', 'none');

%% Argout
if nargout > 0
    argout = h;
end

end