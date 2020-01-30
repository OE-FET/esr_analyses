function varargout = stackplot_xepr(varargin)
%STACKPLOT_XEPR Plots Xepr data as in stacked plots.
%
%   In addition to STACKPLOT, STACKPLOTXEPR will create axis labels and an
%   appropriate legend automatically from pars.
%
% 	SYNTAX:
% 	phandle = STACKPLOT_XEPR(dset)
% 	phandle = STACKPLOT_XEPR(dset, 'OptionName', OptionValue, ...)
%   phandle = STACKPLOT_XEPR(ax, ...)
%
%
%   INPUT(S):
%   ax         - Axis handle for plot. If not given, the data is plotted in
%                the current axis, as returned by gca.
%   dset       - Xrpr data to plot.
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
%

%   $Author: Sam Schott, University of Cambridge <ss2151@cam.ac.uk>$

import esr_analyses.*
import esr_analyses.utils.*

if ishghandle(varargin{1}, 'axes')
    ax = varargin{1};
    dset = varargin{2};
    if nargin > 2
        opts = varargin(3:end);
    else
        opts = {};
    end
else
    ax = gca;
    dset = varargin{1};
    if nargin > 1
        opts = varargin(2:end);
    else
        opts = {};
    end
end

%% Plot

x = dset{:,1};
pars = dset.Properties.UserData;

N = width(dset) - 1;

x_label = sprintf('%s [%s]', pars.XNAM, pars.XUNI);
y_label = sprintf('%s [%s]', pars.IRNAM{1}, pars.IRUNI{1});

phandle = [];

for k = 1:N
    y = dset{:,k+1};
    
    subplot(1,N,k)
    ph_k = stackplot(x, y, opts{:});

    for i = 1:length(pars.z_axis)
        labels{i} = sprintf('%f %s', pars.z_axis(i), pars.YUNI);
    end

    legend(ph_k, labels)

    xlabel(x_label, 'Interpreter', 'none');
    ylabel(y_label);
    title(dset.Properties.VariableNames{k+1}, 'Interpreter', 'none');
    
    phandle = [phandle, ph_k];
end

sgtitle(pars.TITL, 'Interpreter', 'none');

%% Output

if nargout > 0
    varargout{1} = phandle;
end

end
