function varargout = stackplot_xepr(varargin)
%STACKPLOT_XEPR Plots Xepr data as in stacked plots.
%
%   In addition to STACKPLOT, STACKPLOT_XEPR will create axis labels and an
%   appropriate legend automatically from pars.
%
% 	SYNTAX:
% 	phandle = STACKPLOT_XEPR(dset)
% 	phandle = STACKPLOT_XEPR(dset, 'OptionName', OptionValue, ...)
%   phandle = STACKPLOT_XEPR(ax, ...)
%
%
%   INPUT:
%   ax         - Axis handle for plot. If not given, the data is plotted in
%                the current axis, as returned by gca.
%   dset       - Xepr dataset to plot.
%
%   KEYWORD INPUT:
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
% 	OUTPUT:
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
    ax = [];
    dset = varargin{1};
    if nargin > 1
        opts = varargin(2:end);
    else
        opts = {};
    end
end

rescale  = get_kwarg(varargin, 'rescale', 0);

%% Plot

x = dset{:,1};
pars = dset.Properties.UserData;

N = width(dset) - 1;

x_label = sprintf('%s [%s]', pars.XNAM, pars.XUNI);

if rescale
    y_label = sprintf('%s [%s]', pars.IRNAM{1}, pars.IRUNI{1});
else
    y_label = sprintf('%s [scaled]', pars.IRNAM{1});
end

phandle = [];

for k = 1:N
    y = dset{:,k+1};
       
    if isempty(ax)
        ax_target = subplot(1,N,k);
    else
        ax_target = ax;
    end
    ph_k = stackplot(ax_target, x, y, opts{:});

    xlabel(x_label);
    ylabel(y_label);
    title(dset.Properties.VariableNames{k+1}, 'Interpreter', 'none');

    phandle = [phandle, ph_k];
end

if is_2d_exp(dset)
    for i = 1:length(pars.y_axis)
        labels{i} = sprintf('%f %s', pars.y_axis(i), pars.YUNI);
    end
    legend(ph_k, labels, 'Location', 'best','fontsize',10)
end

sgtitle(pars.TITL, 'Interpreter', 'none');

%% Output

if nargout > 0
    varargout{1} = phandle;
end

end
