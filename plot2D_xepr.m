function phandle = plot2D_xepr(varargin)
%PLOT2D_XEPR Plots Xepr data as 2D contour plot.
%
%   In addition to PLOT2D, PLOT2D_XEPR will create axis labels and an
%   appropriate legend automatically from pars.
%
%   SYNTAX:
% 	PLOT2D_XEPR(dset)
%   PLOT2D_XEPR(ax, ...)
%
%   INPUT:
%   ax   - Axis handle for plot. If not given, the data is plotted in
%          the current axis, as returned by gca.
%   dset - Xepr data set
%
% 	OUTPUT:
%   phandle - plot handle
%

import esr_analyses.*
import esr_analyses.utils.*

%% Input Analyses

if ishghandle(varargin{1}, 'axes')
    ax   = varargin{1};
    dset = varargin{2};
else
    ax   = gca;
    dset = varargin{1};
end

x = dset{:,1};
pars = dset.Properties.UserData;
assert_2d_exp(dset)

N = width(dset) - 1;

x_label = sprintf('%s [%s]', pars.XNAM, pars.XUNI);
y_label = sprintf('%s [%s]', pars.YNAM, pars.YUNI);
color_label = sprintf('%s [%s]', pars.IRNAM{1}, pars.IRUNI{1});

for k=1:N
    y = dset{:,k+1};

    %% Plot
    subplot(1,N,k)
    plot2D(x, pars.y_axis, y);
    xlabel(x_label, 'Interpreter', 'none');
    ylabel(y_label, 'Interpreter', 'none');
    cbar = colorbar;
    cbar.Title.String = color_label;
    title(dset.Properties.VariableNames{k+1})
    axis square

end

sgtitle(pars.TITL, 'Interpreter', 'none');

%% Argout
if nargout > 0
    phandle = ax;
end

end