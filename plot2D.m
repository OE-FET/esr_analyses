function argout = plot2D(varargin)
%PLOT2D Plots data as 2D contour plot.
%
% 	PLOT2D(x, y)
% 	PLOT2D(x, y, 'style', style, ...)
%   STACK_PLOT(ax, ...)
%
%   INPUT(S):
%   ax         - Axis handle for plot. If not given, the data is plotted in 
%                the current axis, as returned by gca.
%   x, y       - x- and y-axis data.
%   o          - Ordinate data.
%   'style'    - Line style and color, specified as a character vector or
%                string scalar containing line style symbols, color
%                options, or both. 
%
% 	OUTPUT(S):
%   phandle - plot handle
%
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

style    = get_kwarg(varargin, 'style', '');

%% Data preparation
[X, Y] = meshgrid(x, y);

if ~all(size(o) == [length(y) length(x)])
    o = reshape(o, [length(x) length(y)])';
end

if ~all(size(o) == [length(y) length(x)])
    error('Dimensions of axes and ordinate data do not match.')
end

%% Plot

h = contourf(ax, X, Y, o, style);

colorbar;

if nargout > 1
    argout = h;
end

end