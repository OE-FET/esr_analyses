function argout = plot3Dxepr(varargin)
%PLOT2DXEPR Plots 2D Xepr data as 3D line plot.
%
%   PLOT3DXEPR will create axis labels and an appropriate legend
%   automatically from pars.
%
% 	PLOT3DXEPR(x, y, pars)
%   PLOT3DXEPR(ax, ...)
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

X = ones(size(y)).*x;
Y = ones(size(y)).*pars.z_axis';
Z = y;

fig_title = pars.TITL;
x_label = sprintf('%s [%s]', pars.XNAM, pars.XUNI);
y_label = sprintf('%s [%s]', pars.YNAM, pars.YUNI);
z_label = 'ESR signal [a.u.]';

%% Plot
h = plot3(X, Y, Z);
grid on;

xlabel(ax, x_label, 'Interpreter', 'none');
ylabel(ax, y_label, 'Interpreter', 'none');
zlabel(ax, z_label, 'Interpreter', 'none');
title(ax, fig_title, 'Interpreter', 'none');

%% Argout
if nargout > 0
    argout = h;
end

end