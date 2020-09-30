function [x,y,pars] = dset_to_tuple(dset, channel)
%DSET_TO_TUPLE Converts a dataset table to tuple of [x, o, pars]
%
% INPUT:
% dset    - Xepr dataset table.
% channel - If the are multiple channels of data (e.g., quadrature detection),
%           returns the channel(s) by the given index (indices).
%

import esr_analyses.*

% Number of components
n = width(dset) - 1;

if nargin < 2
    % Only a dataset is given as a function input
    channel = [];
    if n > 1
        % Display data
        stackplot_xepr(dset);
        % Ask user for components to analyse
        channel = input(sprintf('Data set contains %i components. Please select which one you would like to use.\nYou can select multiple components as a vector.\nDefaults to [1,2]: ', n));
    end
    if isempty(channel) && n > 1
        channel = [1,2];
    elseif isempty(channel)
        channel = 1;
    end
end

x = dset{:,1};

% Perform phase cycle
if length(channel) == 2
    sig_x = dset{:, channel(1)+1};
    sig_y = dset{:, channel(2)+1};
    y = phase_cycle(sig_x, sig_y, 'plot', true);
else
    y = dset{:,channel+1};
end

pars = dset.Properties.UserData;

end
