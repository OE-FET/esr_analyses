function [x,y,pars] = dset_to_tuple(dset, slice)

import esr_analyses.*

n = width(dset) - 1;

if nargin < 2
    slice = [];
    if n > 1
        slice = input(sprintf('Data set contains %i components. Please select which one you would like to use.\nIf you select two components, it is assumed that they correspond to 0deg\nand 90deg components.\nDefaults to [1]: ', n));
    end
    if isempty(slice)
        slice = 1;
    end
end

x = dset{:,1};

if isvector(slice)
    sig_x = dset{:, slice(1)+1};
    sig_y = dset{:, slice(2)+1};
    y = phase_cycle(sig_x, sig_y, 'plot', true);
else
    y = dset{:,slice+1};
end

pars = dset.Properties.UserData;

end
