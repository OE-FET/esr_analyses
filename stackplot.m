function [phandle, yoffsets] = stackplot(x, y, varargin)
%STACKPLOT Plots data as in stacked plots.
%
% 	SYNTAX:
% 	STACKPLOT(x, y)
% 	STACKPLOT(x, y, 'xoffset', xoffset, 'yoffsets', yoffsets, 'style', stylestring)
%
% 	If specified, 'yoffsets' is list of offsets between plots. It must have
%   the same length as the number of curves. If not given, yoffset will be
% 	determined automatically. If only a single value is given, it will be
% 	used for all offsets.
%
% 	xoffset allows you to create a staggered plot, offseting each plot by 
% 	the stated number of Gauss. Default is 0 Gauss.
%
% 	When a 'style' is specified, the line style form the variable 
% 	stylestring is used.
%
% 	OUTPUT(S):
%   phandle - plot handles
%   yoffsets - y-offsets used for plotting
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
%   $Date: 2019/05/06 12:58 $    $Revision: 1.1 $

%% get default offset
ydiff = y(:, 1:end-1) + y(:, 2:end);
yoffsets = [0 max(ydiff)*1.3];

%% Input Analyses

xoffset = get_varargin(varargin, 'xoffset', 0);
yoffsets = get_varargin(varargin, 'yoffsets', yoffsets);
style = get_varargin(varargin, 'style', '');

%% The actual programm 
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
if isempty(YOFFSETS) == 0
    yNew = y + YOFFSETS(ones(1, dimy(1)), :);
else
    yNew = y;
end

phandle = plot(xNew, yNew, style);

set(gcf, 'color', 'white');
grid on;

ymin = min(min(yNew)) - abs(0.5*max(max(y)));
ymax = max(max(yNew)) + abs(0.5*max(max(y)));
try
    axis([min(min(xNew)) max(max(xNew)) ymin ymax]);
catch
    axis tight
end
 
end
