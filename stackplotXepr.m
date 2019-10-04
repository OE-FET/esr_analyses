function varargout = stackplotXepr(varargin)
%STACKPLOTXEPR Plots Xepr data as in stacked plots.
%
%   In addition to STACKPLOT, STACKPLOTXEPR will create axis labels and an
%   appropriate legend automatically from pars.
%
% 	SYNTAX:
% 	[phandle, yoffsets] = STACKPLOTXEPR(x, y, pars)
% 	[phandle, yoffsets] = STACKPLOTXEPR(x, y, pars, 'OptionName', OptionValue, ...)
%   [phandle, yoffsets] = STACKPLOT(ax, ...)
%
%
%   INPUT(S):
%   ax         - Axis handle for plot. If not given, the data is plotted in
%                the current axis, as returned by gca.
%   x, y, pars - Xrpr data to plot.
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

if ishghandle(varargin{1}, 'axes')
    pars  = varargin{4};
else
    pars  = varargin{3};
end

%% Plot

[phandle, yoffsets] = stackplot(varargin{:});

ax = gca;

x_label = sprintf('%s [%s]', pars.XNAM, pars.XUNI);

for i = 1:length(pars.z_axis)
    labels{i} = sprintf('%f %s', pars.z_axis(i), pars.YUNI);
end

legend(phandle, labels)

xlabel(ax, x_label, 'Interpreter', 'none');
ylabel(ax, 'ESR signal [a.u.]');
title(ax, pars.TITL, 'Interpreter', 'none');

%% Output

switch nargout
    case 1
        varargout{1} = phandle;
    case 2
        varargout{1} = phandle;
        varargout{2} = yoffsets;
end

end
