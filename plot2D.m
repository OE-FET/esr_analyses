function argout = plot2D(varargin)
%PLOT2D Plots data as 2D contour plot.
%
% 	PLOT2D(x, y, o)
% 	PLOT2D(x, y, o, 'style', style, ...)
%   PLOT2D(ax, ...)
%
%   INPUT(S):
%   ax         - Axis handle for plot. If not given, the data is plotted in
%                the current axis, as returned by gca.
%   x, y       - x- and y-axis data.
%   o          - Ordinate data.
%
% 	OUTPUT(S):
%   phandle - plot handle
%

import esr_analyses.*
import esr_analyses.utils.*

%% Input Analyses

if ishghandle(varargin{1}, 'axes')
    ax = varargin{1};
    x  = varargin{2};
    y  = varargin{3};
    o  = varargin{4};
else
    ax = gca;
    x  = varargin{1};
    y  = varargin{2};
    o  = varargin{3};
end


%% Data preparation
if ~all(size(o) == [length(y) length(x)])
    o = reshape(o, [length(x) length(y)])';
end

if ~all(size(o) == [length(y) length(x)])
    error('Dimensions of axes and ordinate data do not match.')
end

%% Plot
h = imagesc(ax, x, y, o);
ax.YDir = 'normal';

colorbar;

%% Argout
if nargout > 0
    argout = h;
end

end