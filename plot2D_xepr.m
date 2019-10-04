function argout = plot2D_xepr(varargin)
%PLOT2D_XEPR Plots Xepr data as 2D contour plot.
%
%   In addition to PLOT2D, PLOT2D_XEPR will create axis labels and an
%   appropriate legend automatically from pars.
%
% 	PLOT2D_XEPR(x, y, pars)
% 	PLOT2D_XEPR(x, y, pars, 'style', style, ...)
%   PLOT2D_XEPR(ax, ...)
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

if ~isfield(pars,'YNAM')
    error('The given data is not a 2D data set.');
end

%% Data preparation
fig_titel = pars.TITL;
x_label = sprintf('%s [%s]', pars.XNAM, pars.XUNI);
y_label = sprintf('%s [%s]', pars.YNAM, pars.YUNI);
color_label = 'ESR signal [a.u.]';

%% Plot
h = plot2D(ax, x, pars.z_axis, y);
xlabel(ax, x_label, 'Interpreter', 'none');
ylabel(ax, y_label, 'Interpreter', 'none');
title(ax, fig_titel, 'Interpreter', 'none');
cbar = colorbar;
cbar.Title.String = color_label;

%% Argout
if nargout > 0
    argout = h;
end

end