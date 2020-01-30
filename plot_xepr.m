function argout = plot_xepr(varargin)
%PLOT_XEPR Plots Xepr data.
%
%   In addition to PLOT2D, PLOT2D_XEPR will create axis labels and an
%   appropriate legend automatically from pars.
%
% 	PLOT_XEPR(dset)
%   PLOT_XEPR(ax, ...)
%
%   INPUT(S):
%   ax    - Axis handle for plot. If not given, the data is plotted in
%           the current axis, as returned by gca.
%   dsets - Xepr data set
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

x = dset{:,1};
pars = dset.Properties.UserData;

N = width(dset) - 1;

x_label = sprintf('%s [%s]', pars.XNAM, pars.XUNI);
y_label = sprintf('%s [%s]', pars.IRNAM{1}, pars.IRUNI{1});

for k=1:N
    y = dset{:,k+1};

    %% Plot
    subplot(1,N,k)
    plot(x, y);
    xlabel(x_label, 'Interpreter', 'none');
    ylabel(y_label, 'Interpreter', 'none');
    title(dset.Properties.VariableNames{k+1})
    axis square

end

sgtitle(pars.TITL, 'Interpreter', 'none');

%% Argout
if nargout > 0
    argout = ax;
end

end