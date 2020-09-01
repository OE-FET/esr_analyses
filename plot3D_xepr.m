function argout = plot3D_xepr(varargin)
%PLOT2D_XEPR Plots 2D Xepr data as 3D line plot.
%
%   PLOT2D_XEPR will create axis labels and an appropriate legend
%   automatically from pars.
%
% 	PLOT2D_XEPR(x, y, pars)
%   PLOT2D_XEPR(ax, ...)
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
    ax   = varargin{1};
    dset = varargin{2};
else
    ax   = gca;
    dset = varargin{1};
end

%% Data preparation

x = dset{:,1};
pars = dset.Properties.UserData;
assert_2d_exp(dset)

N = width(dset) - 1;

x_label = sprintf('%s [%s]', pars.XNAM, pars.XUNI);
y_label = sprintf('%s [%s]', pars.YNAM, pars.YUNI);

for k=1:N
    y = dset{:,k+1};

    X = ones(size(y)).*x;
    Y = ones(size(y)).*pars.y_axis';
    Z = y;

    z_label = sprintf('%s [%s]', pars.IRNAM{k}, pars.IRUNI{k});

    %% Plot
    subplot(1,N,k)
    plot3(X, Y, Z);
    xlabel(x_label, 'Interpreter', 'none');
    ylabel(y_label, 'Interpreter', 'none');
    zlabel(z_label, 'Interpreter', 'none');
    title(dset.Properties.VariableNames{k+1})
    axis square
    grid on;

end

sgtitle(pars.TITL, 'Interpreter', 'none');

%% Argout
if nargout > 0
    argout = ax;
end

end