function [x,y,pars] = dset_to_tuple(dset, slice)

import esr_analyses.*

% Number of components
n = width(dset) - 1;

if nargin < 2
    % Only a dataset is given as a function input
    slice = [];
    if n > 1
        % Display data
        stackplot_xepr(dset);
        % Ask user for components to analyse
        slice = input(sprintf('Data set contains %i components. Please select which one you would like to use.\nYou can select multiple components as a vector.\nDefaults to [1]: ', n));
    end
    if isempty(slice)
        slice = 1;
    end
end

x = dset{:,1};

% Perform phase cycle
if length(slice) == 2
    sig_x = dset{:, slice(1)+1};
    sig_y = dset{:, slice(2)+1};
    y = phase_cycle(sig_x, sig_y, 'plot', true);
else
    y = dset{:,slice+1};
end

pars = dset.Properties.UserData;

end
