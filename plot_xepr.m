function argout = plot_xepr(dset)
%PLOT_XEPR Plots Xepr data.
%
%   In addition to PLOT2D, PLOT2D_XEPR will create axis labels and an
%   appropriate legend automatically from pars.
%
% 	PLOT_XEPR(dset)
%   PLOT_XEPR(ax, ...)
%
%   INPUT(S):
%   dsets - Xepr data set
%
% 	OUTPUT(S):
%   phandle - image handle
%

import esr_analyses.*
import esr_analyses.utils.*


%% Input Analyses

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
    axis square;
    grid on;
end

set(gca, 'XLimSpec', 'Tight');
sgtitle(pars.TITL, 'Interpreter', 'none');

%% Argout
if nargout > 0
    argout = gca;
end

end