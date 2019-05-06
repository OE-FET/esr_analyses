function [phandle] = stackplot(x, y, varargin)
%STACKPLOT Plots data as in stacked plots.
%
% 	SYNTAX:
% 	stackplot(x, y)
% 	stackplot(x, y, 'xoffset', xoffset, 'yoffset', yoffset, 'style', stylestring)
%
% 	Where yoffset is the degree to which plots are separated if more than
% 	one plot is plotted onto one graph. Default is 0.5 times the spectrum hieght.
%
% 	xoffset allows you to create a staggered plot, offseting each plot by the
% 	stated number of Gauss. Default is 0 Gauss.
%
% 	When a 'style' is specified, the line style form the variable stylestring
% 	is used.
%
% 	OUTPUT(S):
%   phandle - plot handles
%
% 	EXAMPLE: 
%   stackplot(x, y, 'yoffset', 1.2, 'xoffset', 5)
%   Plots y against x with each spectrum (y) being moved up 1.2
%   and across by 5 Gauss
%

% $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$
% $Date: 2018/07/05 12:58 $ $Revision: 0.1 $

%% Input Analyses

xoffset = get_varargin(varargin, 'xoffset', 0);
yoffset = get_varargin(varargin, 'yoffset', max(max(y))*0.5);
style = get_varargin(varargin, 'style', '');

%% The actual programm 
dimx = size(x);
dimy = size(y);

% if every y data has its own x-data, keep. Otherwise, expand x to matrix. 
if dimx(2)  ==  1
	x = x(:, ones(1, dimy(2)));
end

XOFFSET = (0:xoffset:xoffset*(dimy(2)-1));
if isempty(XOFFSET) == 0
    xNew = x + XOFFSET(ones(1, dimx(1)), :);
else
    xNew = x;
end

YOFFSET = (0:yoffset:yoffset*(dimy(2)-1));
if isempty(YOFFSET) == 0
    yNew = y + YOFFSET(ones(1, dimy(1)), :);
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
