function varargout = stackplot(varargin)
%STACKPLOT Plots data as in stacked plots.
%
% 	SYNTAX:
% 	[phandle, yoffsets] = STACKPLOT(x, y)
% 	[phandle, yoffsets] = STACKPLOT(x, y, 'OptionName', OptionValue, ...)
%   [phandle, yoffsets] = STACKPLOT(ax, ...)
%
%
%   INPUT(S):
%   ax         - Axis handle for plot. If not given, the data is plotted in
%                the current axis, as returned by gca.
%   x, y       - Data to plot.
%   'xoffset'  - Horizontal offset between curves, creates a staggered
%                plot. Defaults to 0.
%   'yoffsets' - List of vertical offsets between curves. Automatically
%                determined if not given.
%   'style'    - Line style and color, specified as a character vector or
%                string scalar containing line style symbols, color
%                options, or both.
%   'rescale'  - Rescale all y-data such that max(abs(y)) = 1. Turned off
%                by default.
%
% 	OUTPUT(S):
%   phandle    - plot handles
%   yoffsets   - y-offsets used for plotting
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2019/05/06 12:58 $    $Revision: 1.1 $

import esr_analyses.*
import esr_analyses.utils.*

%% Input Analyses

if ishghandle(varargin{1}, 'axes')
    ax = varargin{1};
    x  = varargin{2};
    y  = varargin{3};
else
    ax = gca;
    x  = varargin{1};
    y  = varargin{2};
end

rescale  = get_kwarg(varargin, 'rescale', 0);

if rescale
    y = y./max(abs(y));
end

% get default offset
ydiff    = y(:, 1:end-1) + y(:, 2:end);
yoffsets = [0 max(abs(ydiff))*1.3];

xoffset  = get_kwarg(varargin, 'xoffset', 0);
yoffsets = get_kwarg(varargin, 'yoffsets', yoffsets);
style    = get_kwarg(varargin, 'style', '');

%% Prepare data with offsets
dimx = size(x);
dimy = size(y);

% if every y data has its own x-data, keep. Otherwise, expand x to matrix.
if dimx(2)  ==  1
	x = x(:, ones(1, dimy(2)));
end

XOFFSETS = (0:xoffset:xoffset*(dimy(2)-1));
if isempty(XOFFSETS) == 0
    xNew = x + XOFFSETS(ones(1, dimx(1)), :);
else
    xNew = x;
end

YOFFSETS = cumsum(yoffsets);
if ~isempty(YOFFSETS)
    yNew = y + YOFFSETS(ones(1, dimy(1)), :);
else
    yNew = y;
end

%% Plot
phandle = plot(ax, xNew, yNew, style);

set(ax.Parent, 'color', 'white');
grid on;

ymin = min(min(yNew)) - abs(0.5*max(max(y)));
ymax = max(max(yNew)) + abs(0.5*max(max(y)));
try
    axis([min(min(xNew)) max(max(xNew)) ymin ymax]);
catch
    axis tight
end

%% Output

switch nargout
    case 1
        varargout{1} = phandle;
    case 2
        varargout{1} = phandle;
        varargout{2} = yoffsets;
end

end
