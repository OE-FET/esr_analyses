function [x,o,pars] = slice_experiment(dset)
%SLICE_EXPERIMENT extracts a slice from a 2D experiment
%
% This will remove all y-axis related parameters from a powersat experiment
% so that it apprears like a 1D experiment.
%
% INPUTS:
% dset   - dataset table
%
% OUTPUTS:
% x    - x-axis data
% o    - ordinate data
% pars - measurement parameters
%


import esr_analyses.*
import esr_analyses.utils.*

[x,o,pars] = dset_to_tuple(dset);

if is_2d_exp(pars)
    stackplot_xepr(dset, 'rescale', true);
    nslice = input('Dataset is 2D. Please select a slice: [1] ');
    if isempty(nslice)
        nslice = 1;
    end
    o = o(:,nslice);
    
    if is_powersat_exp(pars)
        pars.MWPW = pars.y_axis(nslice) * 1e-3; % MW Power in W;
        for field={'y_axis', 'YPTS', 'YMIN', 'YWID', 'YTYP', 'YNAM', 'YUNI'}
            try
                pars = rmfield(pars, field);
            catch ME
                disp(ME)
            end
        end
    end
end

end